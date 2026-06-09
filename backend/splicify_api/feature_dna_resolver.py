"""
Unified DNA lookup for feature KB records.

Problem: the GenoLIB rebuild moved CDS features to protein-only storage
(feature_protein_kb.json + feature_protein.faa), leaving cloning workflows
without DNA for eGFP, Cas9, selection markers, and ~700 other CDSes. This
module is the single place downstream code calls to get DNA for an sseqid,
regardless of whether the KB stores it directly, whether we have a curated
DNA override, or whether we have to back-translate from protein.

Lookup priority for `get_feature_dna(sseqid, ...)`:
  1. `direct`         — feature_reference.fna has a nucleotide record for sseqid.
  2. `json`           — JSON entry's representative_sequence is already DNA.
  3. `curated`        — one of the curated DNA overrides in feature_cds_reference.fna.
  4. `backtranslated` — dnachisel reverse_translate from protein (cached in
                        feature_cds_cache.fna). Tagged with the organism.
  5. `missing`        — nothing available; returns sequence=None with a reason.

Cache: back-translated sequences are appended to `feature_cds_cache.fna` with
headers of the form `>{sseqid}|organism={org}|version={ver}|sha={aa_sha8}`.
Re-runs for the same (sseqid, organism) hit the FASTA cache at startup — no
network, no recompute.

Provenance is always returned so handlers can surface it in user-facing
replies ("back-translated for E. coli — not identical to the Addgene
reference"). That's a hard requirement for any clinical / paper-facing use.
"""
from __future__ import annotations

import hashlib
import logging
import os
import threading
from dataclasses import dataclass
from pathlib import Path
from . import _data
from typing import Dict, Optional, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

_BASE = _data.data_path("feature_db_data")
_DIRECT_FASTA = _BASE / "feature_reference.fna"
_MOTIF_FASTA = _BASE / "feature_motifs.fna"
_CURATED_FASTA = _BASE / "feature_cds_reference.fna"   # opt-in curated overrides
_CACHE_FASTA = _BASE / "feature_cds_cache.fna"         # back-translation cache
_PROTEIN_FASTA = _BASE / "feature_protein.faa"

# Bump when the back-translation rules change (different codon table, new
# constraints). Old cache lines with a stale version are ignored.
_BACKTRANSLATE_VERSION = 1

# Default organism when the caller doesn't specify. E. coli K12 is the right
# default for plasmid propagation context; override to "h_sapiens" for
# mammalian expression CDSes when the host is known.
_DEFAULT_ORGANISM = "e_coli"

_DNA_ALPHABET = set("ACGTN")


# ---------------------------------------------------------------------------
# Lightweight FASTA parsing (no BioPython dep) — same shape as
# plannotate_router._parse_fasta_to_dict but keyed by the first whitespace
# token of the header so cache entries `>sseqid|organism=...` are indexable
# by sseqid alone.
# ---------------------------------------------------------------------------

def _parse_fasta_to_map(path: Path) -> Dict[str, Dict[str, str]]:
    """
    Returns {sseqid: {"sequence": str, "header": str}} from a FASTA file.
    Only the first whitespace-delimited token of the header is used as key.
    Duplicate sseqids keep the first entry (caller can purge stale lines by
    version metadata in the header when relevant).
    """
    out: Dict[str, Dict[str, str]] = {}
    if not path.exists():
        return out
    cur_key: Optional[str] = None
    cur_header: str = ""
    chunks: list = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_key and cur_key not in out:
                    out[cur_key] = {"sequence": "".join(chunks), "header": cur_header}
                header_body = line[1:].strip()
                cur_header = header_body
                cur_key = header_body.split()[0].split("|")[0] if header_body else None
                chunks = []
            elif cur_key is not None:
                chunks.append(line.strip())
    if cur_key and cur_key not in out:
        out[cur_key] = {"sequence": "".join(chunks), "header": cur_header}
    return out


def _parse_cache_fasta(path: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    """
    Cache FASTA headers encode `>sseqid|organism={org}|version={ver}|sha={hex}`.
    Returns keys (sseqid, organism) → {"sequence", "version", "sha"}.
    Stale versions are dropped at load time.
    """
    out: Dict[Tuple[str, str], Dict[str, str]] = {}
    if not path.exists():
        return out
    cur: Optional[Dict[str, str]] = None
    chunks: list = []

    def _flush():
        if cur is None:
            return
        seq = "".join(chunks)
        try:
            ver = int(cur.get("version", "-1"))
        except ValueError:
            ver = -1
        if ver != _BACKTRANSLATE_VERSION:
            return
        key = (cur["sseqid"], cur.get("organism", _DEFAULT_ORGANISM))
        out[key] = {"sequence": seq, "version": ver, "sha": cur.get("sha", "")}

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                _flush()
                parts = line[1:].strip().split("|")
                meta = {"sseqid": parts[0]}
                for p in parts[1:]:
                    if "=" in p:
                        k, v = p.split("=", 1)
                        meta[k.strip()] = v.strip()
                cur = meta
                chunks = []
            elif cur is not None:
                chunks.append(line.strip())
    _flush()
    return out


# ---------------------------------------------------------------------------
# Cached FASTA state — loaded on first call, refreshed if any file's mtime
# changes (cheap check; sidesteps stale caches in long-running workers).
# ---------------------------------------------------------------------------

_lock = threading.Lock()
_cache_state: Dict[str, object] = {
    "mtimes": {},
    "direct": {},       # sseqid → {sequence, header}
    "motifs": {},       # sseqid → {sequence, header}
    "curated": {},      # sseqid → {sequence, header}
    "protein": {},      # sseqid → {sequence, header}
    "backtranslated": {},  # (sseqid, organism) → {sequence, version, sha}
}


def _mtime(path: Path) -> float:
    try:
        return path.stat().st_mtime
    except FileNotFoundError:
        return 0.0


def _refresh_if_stale() -> None:
    paths = {
        "direct": _DIRECT_FASTA,
        "motifs": _MOTIF_FASTA,
        "curated": _CURATED_FASTA,
        "protein": _PROTEIN_FASTA,
        "backtranslated": _CACHE_FASTA,
    }
    with _lock:
        changed = False
        for key, p in paths.items():
            if _cache_state["mtimes"].get(key) != _mtime(p):
                changed = True
                break
        if not changed and _cache_state.get("loaded"):
            return
        _cache_state["direct"] = _parse_fasta_to_map(_DIRECT_FASTA)
        _cache_state["motifs"] = _parse_fasta_to_map(_MOTIF_FASTA)
        _cache_state["curated"] = _parse_fasta_to_map(_CURATED_FASTA)
        _cache_state["protein"] = _parse_fasta_to_map(_PROTEIN_FASTA)
        _cache_state["backtranslated"] = _parse_cache_fasta(_CACHE_FASTA)
        _cache_state["mtimes"] = {k: _mtime(p) for k, p in paths.items()}
        _cache_state["loaded"] = True


# ---------------------------------------------------------------------------
# dnachisel integration
# ---------------------------------------------------------------------------

def _normalize_organism(organism: Optional[str]) -> str:
    if not organism:
        return _DEFAULT_ORGANISM
    o = organism.lower().replace(" ", "_").replace("-", "_")
    aliases = {
        "ecoli": "e_coli",
        "e.coli": "e_coli",
        "e_coli_k12": "e_coli",
        "human": "h_sapiens",
        "homo_sapiens": "h_sapiens",
        "mouse": "m_musculus",
        "mus_musculus": "m_musculus",
        "yeast": "s_cerevisiae",
        "saccharomyces_cerevisiae": "s_cerevisiae",
        "cho": "c_griseus",
        "cricetulus_griseus": "c_griseus",
    }
    return aliases.get(o, o)


def _is_dna(seq: str) -> bool:
    if not seq:
        return False
    # Allow ambiguity codes used by some KBs; require at least 90% {A,C,G,T,N}
    # on the first 100 chars so we don't misread a protein starting with
    # e.g. "Cys-Ala-..." (rare but possible in curated strings).
    head = seq[:100].upper()
    canonical = sum(1 for c in head if c in _DNA_ALPHABET)
    return canonical / len(head) >= 0.9


def _aa_sha8(aa: str) -> str:
    return hashlib.sha1(aa.encode("utf-8")).hexdigest()[:8]


def _back_translate(aa_seq: str, organism: str) -> str:
    """
    Back-translate an amino acid string to DNA using python_codon_tables.
    For each residue we pick the highest-frequency codon in the target
    organism, giving deterministic (cache-stable) output. dnachisel's
    `reverse_translate` is genetic-code-only (first-codon picks, no
    organism bias), so we build the mapping ourselves here.
    """
    import python_codon_tables as pct

    try:
        table = pct.get_codons_table(organism)
    except Exception:
        table = pct.get_codons_table(_DEFAULT_ORGANISM)

    # Reduce to AA -> best codon.
    best: Dict[str, str] = {}
    for aa, codons in table.items():
        if not codons:
            continue
        best[aa] = max(codons.items(), key=lambda kv: kv[1])[0]

    # Ambiguity / rarely-seen residues — fall back to a conservative pick.
    ambiguity_fallback = {
        "B": best.get("D", "GAT"),   # B = D/N → use D
        "Z": best.get("E", "GAA"),   # Z = E/Q → use E
        "X": "NNN",                   # any
        "U": "TGA",                   # selenocysteine — opal stop (caller may remap)
        "O": "TAG",                   # pyrrolysine — amber stop
    }

    out_codons: list = []
    for aa in aa_seq.upper():
        if aa in best:
            out_codons.append(best[aa])
        elif aa in ambiguity_fallback:
            out_codons.append(ambiguity_fallback[aa])
        else:
            # Residue we don't recognize (whitespace, digit, etc.) — skip.
            continue

    dna = "".join(out_codons)
    if not dna.endswith(("TAA", "TAG", "TGA")):
        # Organism-appropriate stop codon (most-frequent *).
        stop = best.get("*", "TAA") or "TAA"
        dna = dna + stop
    return dna


def _append_cache(sseqid: str, organism: str, dna: str, aa_sha: str) -> None:
    header = f">{sseqid}|organism={organism}|version={_BACKTRANSLATE_VERSION}|sha={aa_sha}"
    line = f"{header}\n{dna}\n"
    _BASE.mkdir(parents=True, exist_ok=True)
    with open(_CACHE_FASTA, "a") as f:
        f.write(line)
    # Update in-memory map so the same resolver call in this process returns
    # the cached entry without reparsing the file.
    with _lock:
        _cache_state["backtranslated"][(sseqid, organism)] = {
            "sequence": dna, "version": _BACKTRANSLATE_VERSION, "sha": aa_sha,
        }
        _cache_state["mtimes"]["backtranslated"] = _mtime(_CACHE_FASTA)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

@dataclass
class ResolvedFeature:
    sseqid: str
    sequence: Optional[str]
    provenance: str            # "direct" | "curated" | "json" | "backtranslated_{org}" | "missing"
    organism: Optional[str]    # populated when back-translated
    length: int
    notes: str = ""


def get_feature_dna(
    sseqid: str,
    *,
    protein_sequence: Optional[str] = None,
    json_representative_sequence: Optional[str] = None,
    organism: Optional[str] = None,
    allow_backtranslation: bool = True,
) -> ResolvedFeature:
    """
    Resolve DNA for a feature sseqid across all sources.

    Arguments:
      sseqid                         — canonical feature id (e.g. "egfp_rs_U57607_612")
      protein_sequence               — if the caller already has AA from the KB, pass it
                                       here to skip the .faa lookup.
      json_representative_sequence   — if the KB JSON's representative_sequence is
                                       already DNA-looking, pass it as a shortcut.
      organism                       — codon-table organism for back-translation.
                                       Aliases: ecoli, human, mouse, yeast, cho.
      allow_backtranslation          — when False, never synthesize DNA; return
                                       provenance="missing" if no DNA is available.
                                       Use False for workflows that require natural
                                       reference sequences (SDM against a cited variant).

    Returns a ResolvedFeature with provenance set so handlers can decide
    whether to warn the user about synthetic DNA.
    """
    _refresh_if_stale()
    organism = _normalize_organism(organism)

    # Tier 1 — direct .fna hit
    hit = _cache_state["direct"].get(sseqid) or _cache_state["motifs"].get(sseqid)  # type: ignore[attr-defined]
    if hit and _is_dna(hit["sequence"]):
        seq = hit["sequence"].upper()
        return ResolvedFeature(sseqid=sseqid, sequence=seq, provenance="direct",
                               organism=None, length=len(seq))

    # Tier 2 — curated overrides
    cur = _cache_state["curated"].get(sseqid)  # type: ignore[attr-defined]
    if cur and _is_dna(cur["sequence"]):
        seq = cur["sequence"].upper()
        return ResolvedFeature(sseqid=sseqid, sequence=seq, provenance="curated",
                               organism=None, length=len(seq))

    # Tier 3 — JSON already had DNA (caller passed it through)
    if json_representative_sequence and _is_dna(json_representative_sequence):
        seq = json_representative_sequence.upper()
        return ResolvedFeature(sseqid=sseqid, sequence=seq, provenance="json",
                               organism=None, length=len(seq))

    # Resolve the protein source — caller-provided wins, else .faa lookup.
    aa = (protein_sequence or "").strip()
    if not aa:
        p = _cache_state["protein"].get(sseqid)  # type: ignore[attr-defined]
        if p:
            aa = p["sequence"].strip()

    if not aa and not _is_dna(json_representative_sequence or ""):
        # If json_representative_sequence looks like protein, use it as AA.
        if json_representative_sequence and all(
            c.isalpha() for c in json_representative_sequence[:50]
        ):
            aa = json_representative_sequence.strip()

    # Tier 4 — back-translate from protein
    if allow_backtranslation and aa:
        aa_sha = _aa_sha8(aa)
        cached = _cache_state["backtranslated"].get((sseqid, organism))  # type: ignore[attr-defined]
        if cached and cached.get("sha") == aa_sha:
            seq = cached["sequence"].upper()
            return ResolvedFeature(
                sseqid=sseqid, sequence=seq,
                provenance=f"backtranslated_{organism}", organism=organism,
                length=len(seq),
                notes="Synthetic DNA — back-translated from protein KB. "
                      "Not guaranteed to match any physical reference strain.",
            )
        try:
            dna = _back_translate(aa, organism)
            _append_cache(sseqid, organism, dna, aa_sha)
            return ResolvedFeature(
                sseqid=sseqid, sequence=dna.upper(),
                provenance=f"backtranslated_{organism}", organism=organism,
                length=len(dna),
                notes="Synthetic DNA — back-translated from protein KB. "
                      "Not guaranteed to match any physical reference strain.",
            )
        except Exception as exc:
            logger.warning("back-translation failed for %s (%s): %s",
                           sseqid, organism, exc)

    # Tier 5 — nothing available
    reason = (
        "No DNA in feature_reference.fna, no curated override, no JSON DNA, "
        "and " + (
            "back-translation disabled" if not allow_backtranslation
            else "no protein sequence in feature_protein.faa either."
        )
    )
    return ResolvedFeature(sseqid=sseqid, sequence=None, provenance="missing",
                           organism=None, length=0, notes=reason)


def backfill_record(record: Dict[str, object], *, organism: Optional[str] = None,
                    allow_backtranslation: bool = True) -> Dict[str, object]:
    """
    Backfill a pLannotate-style KB record dict in place so downstream readers
    (who look at intrinsic_properties.sequence_derived.representative_sequence)
    see DNA.

    Returns the same dict with an added `dna_resolution` block containing
    sseqid, provenance, organism, notes. Non-destructive if the record already
    has DNA.
    """
    sseqid = str(record.get("sseqid") or "")
    if not sseqid:
        return record
    ip = record.setdefault("intrinsic_properties", {})
    sd = ip.setdefault("sequence_derived", {})
    existing = sd.get("representative_sequence") or ""
    protein = existing if existing and not _is_dna(existing) else None

    resolved = get_feature_dna(
        sseqid,
        protein_sequence=protein,
        json_representative_sequence=existing if _is_dna(existing) else None,
        organism=organism,
        allow_backtranslation=allow_backtranslation,
    )
    if resolved.sequence:
        sd["representative_sequence"] = resolved.sequence
    if protein:
        # Preserve the protein sequence in its own field so we don't lose it.
        sd.setdefault("protein_sequence", protein)
    record["dna_resolution"] = {
        "sseqid": resolved.sseqid,
        "provenance": resolved.provenance,
        "organism": resolved.organism,
        "length": resolved.length,
        "notes": resolved.notes,
    }
    return record
