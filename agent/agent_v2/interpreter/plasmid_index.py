"""Deterministic index over a single plasmid's annotation envelope.

A PlasmidIndex wraps the response from /plannotate/annotate_sequence_llm
(or the same payload as cached by annotate_llm_cached) and exposes
fast, no-LLM lookup methods that the interpreter agent / tools call.

All methods return JSON-serialisable dicts so results can be embedded
in Anthropic tool results without further marshalling. No raw DNA
sequences are returned except when the caller asks specifically (e.g.
the codon of a particular AA position) — bulk sequence stays inside
the index.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Optional


_AA_FULL_NAME = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic acid",
    "C": "Cysteine", "E": "Glutamic acid", "Q": "Glutamine", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
    "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
    "*": "Stop", "X": "Unknown",
}

# Heuristic mapping from module-type signatures to plasmid application
# labels. Each entry is (required_module_types: tuple, label, notes).
# Required types are matched as substrings against the module_type field.
_APPLICATION_PATTERNS = [
    (
        ("guide_expression_cassette", "mammalian_pol2_expression_cassette", "lentiviral_payload"),
        "CRISPR-Cas9 lentiviral knockout / knock-in vector",
        "U6 + sgRNA cassette inside a lentiviral payload with a Pol II–driven Cas9.",
    ),
    (
        ("guide_expression_cassette", "mammalian_pol2_expression_cassette"),
        "Mammalian CRISPR-Cas9 expression vector",
        "U6 + sgRNA cassette paired with a Pol II–driven Cas9 / payload CDS.",
    ),
    (
        ("guide_expression_cassette",),
        "sgRNA / guide-RNA expression vector",
        "Pol III sgRNA cassette without a paired Cas9 CDS — co-transfection vector.",
    ),
    (
        ("lac_alpha_disrupted_module",),
        "Cloned blue/white screening vector (insert disrupts lacZα)",
        "pBSK-derived vector with an insert replacing the α-fragment — successful clone.",
    ),
    (
        ("lac_alpha_blue_white_module",),
        "Blue/white screening cloning vector (intact lacZα)",
        "Empty pBSK-class vector ready to receive an insert.",
    ),
    (
        ("mammalian_pol2_expression_cassette", "lentiviral_payload"),
        "Lentiviral mammalian expression vector",
        "Pol II–driven payload inside a lentiviral payload module.",
    ),
    (
        ("aav_itr", "mammalian_pol2_expression_cassette"),
        "AAV mammalian expression vector",
        "AAV ITRs flanking a Pol II expression cassette.",
    ),
    (
        ("mammalian_pol2_expression_cassette",),
        "Mammalian expression vector",
        "Pol II promoter → CDS → polyA in mammalian cells.",
    ),
    (
        ("gateway_destination_cassette",),
        "Gateway destination vector",
        "attR1 / attR2 cassette ready for LR recombination.",
    ),
    (
        ("transposon",),
        "Transposon delivery vector",
        "Sleeping Beauty / piggyBac–style ITR-flanked payload.",
    ),
]


@dataclass
class PlasmidIndex:
    """Wraps one /plannotate/annotate_sequence_llm response."""

    plasmid_id: str
    name: Optional[str]
    envelope: dict[str, Any]
    sequence: str = field(default="", repr=False)

    # Lazily computed
    _modules_by_type: Optional[dict[str, list[dict[str, Any]]]] = field(default=None, repr=False)

    # ──────────────────────────────────────────────────────────────────
    # Construction
    # ──────────────────────────────────────────────────────────────────
    @classmethod
    def from_envelope(cls, plasmid_id: str, envelope: dict[str, Any], name: Optional[str] = None) -> "PlasmidIndex":
        return cls(
            plasmid_id=plasmid_id,
            name=name,
            envelope=envelope,
            sequence=envelope.get("sequence", "") or "",
        )

    # ──────────────────────────────────────────────────────────────────
    # Convenience accessors
    # ──────────────────────────────────────────────────────────────────
    def annotations(self) -> list[dict[str, Any]]:
        return self.envelope.get("annotations") or self.envelope.get("plannotate_annotations") or []

    def modules(self) -> list[dict[str, Any]]:
        return self.envelope.get("modules") or []

    def hierarchical(self) -> list[dict[str, Any]]:
        return self.envelope.get("hierarchical_annotations") or []

    def interactions(self) -> list[dict[str, Any]]:
        return self.envelope.get("interactions") or []

    def cloning_features(self) -> list[dict[str, Any]]:
        cf = self.envelope.get("cloning_features")
        if isinstance(cf, dict):
            return cf.get("features") or []
        if isinstance(cf, list):
            return cf
        return []

    def translation_annotations(self) -> list[dict[str, Any]]:
        return [h for h in self.hierarchical() if h.get("module_type") == "translation"]

    # ──────────────────────────────────────────────────────────────────
    # Public lookup methods
    # ──────────────────────────────────────────────────────────────────
    def summary(self) -> dict[str, Any]:
        """High-level overview used as a quick orientation pass for the
        interpreter agent before it digs into specific lookups."""
        mods = self.modules()
        anns = self.annotations()
        types = sorted({m.get("module_type") for m in mods if m.get("module_type")})
        rules = sorted({m.get("rule_id") for m in mods if m.get("rule_id")})
        return {
            "plasmid_id": self.plasmid_id,
            "name": self.name,
            "length_bp": len(self.sequence),
            "n_annotations": len(anns),
            "n_modules": len(mods),
            "n_interactions": len(self.interactions()),
            "n_translation_annotations": len(self.translation_annotations()),
            "n_cloning_features": len(self.cloning_features()),
            "module_types": types,
            "rule_ids": rules,
            "feature_names": sorted({a.get("name") for a in anns if a.get("name")}),
        }

    def find_modules(self, query: str) -> list[dict[str, Any]]:
        """Find modules whose module_type, rule_id, or name matches the
        query (case-insensitive substring). Returns coords + submodule
        list + metadata.

        Example queries:
            'guide_expression_cassette', 'sgrna', 'golden_gate',
            'lac_alpha_disrupted', 'POL3-GG-01', 'lentiviral_payload'
        """
        q = (query or "").lower()
        out = []
        for m in self.modules():
            mt = (m.get("module_type") or "").lower()
            rid = (m.get("rule_id") or "").lower()
            nm = (m.get("name") or "").lower()
            if q in mt or q in rid or q in nm:
                out.append({
                    "plasmid_id": self.plasmid_id,
                    "module_type": m.get("module_type"),
                    "rule_id": m.get("rule_id"),
                    "name": m.get("name"),
                    "start": m.get("start"),
                    "end": m.get("end"),
                    "strand": m.get("strand", m.get("direction")),
                    "submodules": m.get("submodules", []),
                    "metadata": m.get("metadata", {}),
                    "golden_gate": m.get("golden_gate"),
                    "notes": m.get("notes"),
                })
        return out

    def find_features(self, query: str) -> list[dict[str, Any]]:
        """Find pLannotate features matching query against name, sseqid,
        gene_name, protein_name, or aliases. Returns full annotation
        records (no raw sequence)."""
        q = (query or "").lower().strip()
        if not q:
            return []
        out = []
        for a in self.annotations():
            kb = a.get("kb_data") or {}
            haystack = [
                a.get("name") or "",
                a.get("sseqid") or "",
                kb.get("gene_name") or "",
                kb.get("protein_name") or "",
                kb.get("entry_name") or "",
                kb.get("feature_name") or "",
            ]
            haystack = " | ".join(s.lower() for s in haystack if s)
            if q in haystack:
                out.append({
                    "plasmid_id": self.plasmid_id,
                    "name": a.get("name"),
                    "type": a.get("type") or a.get("feature_type"),
                    "start": a.get("start"),
                    "end": a.get("end"),
                    "direction": a.get("direction") or a.get("strand"),
                    "description": a.get("description"),
                    "kb_data": kb,
                })
        return out

    def find_cloning_features(self, query: str) -> list[dict[str, Any]]:
        """Cloning features: restriction sites, Gateway att sites, PCR
        warnings. Query matches name / subtype / feature_family."""
        q = (query or "").lower().strip()
        if not q:
            return []
        out = []
        for c in self.cloning_features():
            haystack = " ".join(str(c.get(k, "")).lower()
                                for k in ("name", "subtype", "feature_family"))
            if q in haystack:
                out.append({
                    "plasmid_id": self.plasmid_id,
                    "name": c.get("name"),
                    "feature_family": c.get("feature_family"),
                    "subtype": c.get("subtype"),
                    "start": c.get("start"),
                    "end": c.get("end"),
                    "cut_profile": c.get("cut_profile"),
                })
        return out

    def lookup_amino_acid(self, feature_name: str, aa_index: int) -> Optional[dict[str, Any]]:
        """Resolve a 1-based AA position within a feature (or within an
        ORF directly). Walks translation annotations to find a region
        whose name matches feature_name, then returns the residue at the
        given position with the codon, ORF position, and feature-local
        position.

        Returns None if the feature isn't found in any ORF or the
        index is out of range.
        """
        if aa_index < 1:
            return None
        feat_q = (feature_name or "").lower().strip()
        if not feat_q:
            return None

        for t in self.translation_annotations():
            meta = t.get("metadata") or {}
            aa_seq = meta.get("aa_sequence") or ""
            regions = meta.get("feature_regions") or []
            orf_strand = t.get("direction") or t.get("strand") or 1
            orf_start = t.get("start")
            orf_end = t.get("end")
            orf_aa_len = meta.get("aa_length") or len(aa_seq)

            # Try every region; substring match is fine because feature
            # names tend to be unique enough at this level.
            for r in regions:
                if feat_q not in (r.get("name") or "").lower():
                    continue
                r_start = int(r.get("aa_start", 0))
                r_end = int(r.get("aa_end", 0))
                r_len = r_end - r_start + 1
                if aa_index > r_len:
                    return {
                        "ok": False,
                        "plasmid_id": self.plasmid_id,
                        "reason": (
                            f"{r.get('name')} is {r_len} aa long; aa_index {aa_index} "
                            "is past the end of the feature."
                        ),
                        "feature_aa_length": r_len,
                    }
                aa_idx_in_orf = r_start + aa_index - 1  # 1-based
                if aa_idx_in_orf < 1 or aa_idx_in_orf > len(aa_seq):
                    return None
                letter = aa_seq[aa_idx_in_orf - 1]

                # nt position of the codon for this AA. Forward-strand:
                # codon starts at orf_start + 3*(idx-1). Reverse strand:
                # codon starts at orf_end - 3*idx (last 3 nt of the
                # reverse-translated chunk).
                if orf_strand == -1:
                    codon_nt_start = orf_end - 3 * aa_idx_in_orf
                else:
                    codon_nt_start = orf_start + 3 * (aa_idx_in_orf - 1)
                codon = ""
                if self.sequence and codon_nt_start is not None:
                    raw = self.sequence[codon_nt_start:codon_nt_start + 3]
                    if orf_strand == -1:
                        # Reverse complement
                        comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
                        codon = "".join(comp.get(b, "N") for b in reversed(raw.upper()))
                    else:
                        codon = raw.upper()

                return {
                    "ok": True,
                    "plasmid_id": self.plasmid_id,
                    "feature_name": r.get("name"),
                    "letter": letter,
                    "amino_acid": _AA_FULL_NAME.get(letter, "Unknown"),
                    "codon": codon,
                    "aa_position_in_feature": aa_index,
                    "feature_aa_length": r_len,
                    "aa_position_in_orf": aa_idx_in_orf,
                    "orf_aa_length": orf_aa_len,
                    "orf_name": t.get("name"),
                    "orf_start": orf_start,
                    "orf_end": orf_end,
                    "codon_nt_start": codon_nt_start,
                    "strand": orf_strand,
                }

        # No region matched — see if the query names the ORF itself.
        for t in self.translation_annotations():
            if feat_q in (t.get("name") or "").lower():
                meta = t.get("metadata") or {}
                aa_seq = meta.get("aa_sequence") or ""
                if aa_index > len(aa_seq):
                    return {
                        "ok": False,
                        "plasmid_id": self.plasmid_id,
                        "reason": f"ORF is {len(aa_seq)} aa long; aa_index {aa_index} out of range.",
                        "feature_aa_length": len(aa_seq),
                    }
                letter = aa_seq[aa_index - 1]
                return {
                    "ok": True,
                    "plasmid_id": self.plasmid_id,
                    "feature_name": t.get("name"),
                    "letter": letter,
                    "amino_acid": _AA_FULL_NAME.get(letter, "Unknown"),
                    "aa_position_in_feature": aa_index,
                    "feature_aa_length": len(aa_seq),
                    "aa_position_in_orf": aa_index,
                    "orf_aa_length": len(aa_seq),
                    "orf_name": t.get("name"),
                }

        return None

    def expression_cassette_for(self, cds_name: str) -> Optional[dict[str, Any]]:
        """Find the expression cassette interaction that contains a CDS
        matching cds_name. Returns the promoter (upstream regulatory)
        and polyA (downstream regulatory) names + coords.

        Walks `interactions` looking for one whose components reference
        a CDS module whose name overlaps cds_name. Falls back to
        scanning `modules` for a `mammalian_pol2_expression_cassette`
        whose body covers an annotation matching cds_name.
        """
        q = (cds_name or "").lower().strip()
        if not q:
            return None

        # First check interactions for a cassette pairing.
        for ix in self.interactions():
            ix_type = (ix.get("type") or ix.get("interaction_type") or "").lower()
            if "expression" not in ix_type and "cassette" not in ix_type:
                continue
            comps = ix.get("components") or ix.get("modules") or []
            comp_names = " | ".join((c.get("name") or "").lower() for c in comps if isinstance(c, dict))
            if q not in comp_names:
                continue
            ur = next((c for c in comps if isinstance(c, dict)
                       and "upstream" in (c.get("role") or c.get("module_type") or "").lower()), None)
            dr = next((c for c in comps if isinstance(c, dict)
                       and "downstream" in (c.get("role") or c.get("module_type") or "").lower()), None)
            return {
                "plasmid_id": self.plasmid_id,
                "found_via": "interaction",
                "interaction_type": ix.get("type") or ix.get("interaction_type"),
                "cds_name": cds_name,
                "promoter": _component_summary(ur, role="promoter"),
                "polyA": _component_summary(dr, role="polya"),
                "components": [{"name": c.get("name"), "module_type": c.get("module_type"),
                                "start": c.get("start"), "end": c.get("end")} for c in comps if isinstance(c, dict)],
            }

        # Fallback: scan Pol II cassette modules; find one whose range
        # overlaps an annotation matching the CDS name.
        cds_annotation = next((a for a in self.annotations()
                                if q in (a.get("name") or "").lower()
                                and (a.get("type") or "").upper() == "CDS"), None)
        if not cds_annotation:
            return None
        a_start, a_end = int(cds_annotation.get("start", 0)), int(cds_annotation.get("end", 0))
        for m in self.modules():
            if (m.get("module_type") or "") != "mammalian_pol2_expression_cassette":
                continue
            if not (m.get("start", 0) <= a_start and m.get("end", 0) >= a_end):
                continue
            subs = m.get("submodules") or []
            ur = next((s for s in subs if "upstream" in (s.get("module_type") or "").lower()), None)
            dr = next((s for s in subs if "downstream" in (s.get("module_type") or "").lower()), None)
            return {
                "plasmid_id": self.plasmid_id,
                "found_via": "module",
                "cassette_module_id": m.get("module_id") or m.get("name"),
                "cds_name": cds_annotation.get("name"),
                "promoter": _component_summary(ur, role="promoter"),
                "polyA": _component_summary(dr, role="polya"),
                "cassette_start": m.get("start"),
                "cassette_end": m.get("end"),
            }
        return None

    def infer_application(self) -> dict[str, Any]:
        """Pattern-match the module composition to a likely plasmid
        application. Returns the best-matching label plus a confidence
        score (# of required types present / # required) and a list of
        rejected patterns so the caller can see the evidence."""
        present = {m.get("module_type") or "" for m in self.modules()}
        present_str = " | ".join(present).lower()

        rejected = []
        for required, label, note in _APPLICATION_PATTERNS:
            hits = [t for t in required if any(t in p for p in present_str.split(" | "))]
            score = len(hits) / max(1, len(required))
            if score >= 1.0:
                return {
                    "plasmid_id": self.plasmid_id,
                    "application": label,
                    "confidence": "high" if len(required) >= 2 else "medium",
                    "evidence": list(hits),
                    "notes": note,
                }
            rejected.append({"label": label, "hits": hits, "missing": [t for t in required if t not in hits]})

        return {
            "plasmid_id": self.plasmid_id,
            "application": "Unknown / unclassified",
            "confidence": "low",
            "evidence": list(present),
            "notes": (
                "No application pattern matched the module set. Module types present: "
                + (", ".join(sorted(present)) if present else "(none)")
            ),
            "rejected_patterns": rejected[:5],
        }


def _component_summary(comp: Optional[dict[str, Any]], role: str) -> Optional[dict[str, Any]]:
    """Tighten a regulatory submodule into the fields the agent needs."""
    if not comp:
        return None
    meta = comp.get("metadata") or {}
    out = {
        "name": comp.get("name"),
        "module_type": comp.get("module_type"),
        "start": comp.get("start"),
        "end": comp.get("end"),
        "strand": comp.get("strand") or comp.get("direction"),
    }
    if role == "promoter":
        out["promoter_name"] = (meta.get("promoter_name") or meta.get("name")
                                 or comp.get("name"))
        if "components" in meta:
            out["components"] = meta["components"]
    elif role == "polya":
        out["polya_signal"] = meta.get("polya_signal") or meta.get("polyA") or comp.get("name")
        if "components" in meta:
            out["components"] = meta["components"]
    return out
