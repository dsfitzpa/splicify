"""resolve_feature_position — explicit feature-relative coordinate resolution.

When a user references a feature in an uploaded plasmid (e.g. "D10A in Cas9",
"-35 box of CMV promoter", "residue 65 of GFP"), the agent should NOT do
coordinate math itself. It should call this tool to get the deterministic
plasmid position + codon + amino acid.

Returns a small digest — no raw DNA leaves the tool.
"""
from __future__ import annotations

from typing import Any, Optional

import agent_v2  # noqa: F401 — triggers path shim
from agent_v2 import attachment_kinds


_CODON = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _slice_circular(seq: str, start: int, length: int) -> str:
    """Slice `length` bp starting at `start`, wrapping if circular plasmid runs out."""
    n = len(seq)
    if start < 0:
        start = (start % n + n) % n
    if start + length <= n:
        return seq[start:start + length]
    return seq[start:] + seq[:(start + length) % n]


def _find_feature(features: list[dict[str, Any]], needle: str) -> tuple[Optional[dict], list[str]]:
    """Return (best_match, candidate_names). Exact lower-case match wins; then substring."""
    needle_l = (needle or "").strip().lower()
    if not needle_l:
        return None, []
    exact: list[dict] = []
    partial: list[dict] = []
    names: list[str] = []
    for f in features or []:
        name = (f.get("Feature") or f.get("name") or "").strip()
        names.append(name)
        nl = name.lower()
        if nl == needle_l:
            exact.append(f)
        elif needle_l in nl or nl in needle_l:
            partial.append(f)
    best = (exact + partial)[:1]
    return (best[0] if best else None), names[:30]


def _walk_intervals_to_codon(
    intervals: list[tuple[int, int]],
    strand: int,
    sequence: str,
    cumulative_bp_offset: int,
) -> tuple[Optional[int], Optional[str], bool, Optional[int], Optional[int]]:
    """Map a 0-indexed CDS-relative bp offset to a genomic coordinate + codon.

    intervals: ordered list of (start, end) half-open genomic coords for the
    CDS's exon segments. strand: +1 or -1.

    Returns (genomic_position, codon_3bp, spans_intron, lowest_plus_strand_pos,
    highest_plus_strand_pos). `genomic_position` is the 0-indexed sense-strand
    5'-most base (highest + strand coord for - strand genes). The last two
    fields are the lowest / highest + strand 0-indexed coordinates of the
    3 codon bases - construct + strand edit windows from these. For - strand
    the returned codon is the SENSE codon (reverse-complement of the genomic
    bases). spans_intron=True iff the 3 codon bases straddle an exon boundary
    inside the joined CDS.
    """
    # Sort intervals into spliced-CDS reading order (5' -> 3').
    parts = sorted(intervals, key=lambda p: p[0])
    if strand < 0:
        parts = list(reversed(parts))

    # Accumulate spliced bp positions and pick out the 3 bases at the offset.
    bp_collected: list[tuple[int, int]] = []  # (genomic_pos_0idx, base_index_in_codon)
    running = 0
    for seg_start, seg_end in parts:
        seg_len = seg_end - seg_start
        if cumulative_bp_offset >= running + seg_len:
            running += seg_len
            continue
        # Pull bases out of this segment.
        for k in range(seg_len):
            if running + k < cumulative_bp_offset:
                continue
            # local bp index within segment
            if strand >= 0:
                genomic_pos = seg_start + k
            else:
                genomic_pos = seg_end - 1 - k
            bp_collected.append((genomic_pos, len(bp_collected)))
            if len(bp_collected) == 3:
                break
        if len(bp_collected) >= 3:
            break
        running += seg_len

    if len(bp_collected) < 3:
        return (None, None, False, None, None)

    # Genomic position of the codon's first sense base.
    first_pos = bp_collected[0][0]
    lowest_plus = min(p for p, _ in bp_collected)
    highest_plus = max(p for p, _ in bp_collected)

    # Detect intron crossing: any two consecutive sense bases that are not
    # adjacent on the genome (allowing for + or - strand step of 1) means
    # an exon boundary lies between them.
    spans_intron = False
    for i in range(1, 3):
        prev = bp_collected[i - 1][0]
        cur = bp_collected[i][0]
        expected_delta = 1 if strand >= 0 else -1
        if cur - prev != expected_delta:
            spans_intron = True
            break

    # Pull the 3 bases and produce the sense codon.
    bases = "".join(sequence[p].upper() for p, _ in bp_collected)
    if strand < 0:
        # The bases above were already collected in sense-reading order
        # (seg_end-1, seg_end-2, ...) so we need to revcomp each base
        # individually rather than reverse the whole string.
        comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        bases = "".join(comp.get(b, b) for b in bases)

    return (first_pos, bases, spans_intron, lowest_plus, highest_plus)


async def _resolve_genomic(
    args: dict[str, Any],
    registry: Any,
) -> dict[str, Any]:
    """Genomic-kind feature resolution: walks exon-joined CDS, returns
    spliced-CDS-aware codon + amino acid for aa_residue queries.
    """
    aid = args.get("attachment_id")
    feature_name = (args.get("feature_name") or "").strip()
    kind = args.get("kind", "aa_residue")
    offset = int(args.get("offset") or 0)
    att = registry.get(aid)

    ann = attachment_kinds.get_genomic_annotation(aid)
    if ann is None:
        return {"ok": False, "error": f"no genomic annotation cached for {aid!r}"}

    # Match feature_name against gene, transcript_id, protein_id.
    needle = feature_name.lower()
    cds_matches = []
    for ft in ann.features:
        if ft.type != "CDS":
            continue
        names = [ft.gene, ft.transcript_id, ft.protein_id]
        if any(n and n.lower() == needle for n in names):
            cds_matches.append(ft)
    if not cds_matches:
        for ft in ann.features:
            if ft.type != "CDS":
                continue
            names = [ft.gene, ft.transcript_id, ft.protein_id]
            if any(n and (needle in n.lower() or n.lower() in needle) for n in names):
                cds_matches.append(ft)
    if not cds_matches:
        avail = sorted({(ft.gene or ft.transcript_id or "") for ft in ann.features
                          if ft.type == "CDS" and (ft.gene or ft.transcript_id)})
        return {
            "ok": False,
            "error": f"no CDS matching {feature_name!r}",
            "available_features": avail[:30],
        }
    def _len_no_stop_for_sort(t):
        t = t or ""
        return len(t[:-1] if t.endswith("*") else t)
    # If the user passed a specific protein_id / transcript_id as
    # feature_name, that exact match takes priority. Otherwise (the common
    # gene-symbol case where 'CGAS' matches multiple isoforms) prefer the
    # LONGEST translation — the canonical isoform — so residues near the
    # protein's C-terminus don't silently fall out of range on a
    # short alternative-splice isoform.
    def _is_exact_id_match(c):
        return any(needle == (n or "").lower()
                    for n in (c.transcript_id, c.protein_id))
    exact_id = next((c for c in cds_matches if _is_exact_id_match(c)), None)
    if exact_id is not None:
        cds = exact_id
        isoform_pick_reason = "exact protein/transcript ID match"
    else:
        cds_matches_sorted = sorted(
            cds_matches,
            key=lambda c: (-_len_no_stop_for_sort(c.translation), c.start),
        )
        cds = cds_matches_sorted[0]
        isoform_pick_reason = (
            f"longest of {len(cds_matches)} CDS isoforms matching "
            f"{feature_name!r} by gene name"
            if len(cds_matches) > 1 else "single CDS match"
        )

    # Strip a trailing stop codon ('*') from the translation so cds_length_aa
    # and the in-range error message report the actual protein length the
    # user thinks in terms of (CGAS canonical is 522 aa, NOT 523-with-stop).
    # Residue offset=N still maps to translation[N-1] correctly for N in
    # 1..len(translation), so the indexing path below is unchanged — only
    # the LENGTH that the bounds check uses is trimmed.
    raw_translation = cds.translation or ""
    translation = raw_translation[:-1] if raw_translation.endswith("*") else raw_translation
    spliced_len = sum(e - s for s, e in cds.intervals)
    sense = "+" if cds.strand >= 0 else "-"

    # Expose every isoform alternative so the LLM (or downstream tools)
    # can disambiguate by passing feature_name=<protein_id> on a follow-up
    # call. Sorted longest-first for legibility.
    def _len_no_stop(t):
        t = t or ""
        return len(t[:-1] if t.endswith("*") else t)
    alternative_isoforms = sorted(
        [{
            "transcript_id": c.transcript_id,
            "protein_id": c.protein_id,
            "length_aa": _len_no_stop(c.translation),
            "n_exons": len(c.intervals or []),
            "feature_start": c.start,
            "feature_end": c.end,
        } for c in cds_matches],
        key=lambda r: -(r.get("length_aa") or 0),
    )

    base = {
        "feature_name": cds.gene or cds.transcript_id or cds.protein_id,
        "transcript_id": cds.transcript_id,
        "protein_id": cds.protein_id,
        "feature_strand": sense,
        "feature_start": cds.start,
        "feature_end": cds.end,
        "feature_length_bp": cds.end - cds.start,
        "spliced_cds_length_bp": spliced_len,
        "n_exons": len(cds.intervals),
        "cds_length_aa": len(translation),
        "isoform_pick_reason": isoform_pick_reason,
        "alternative_isoforms": alternative_isoforms,
        "candidates_considered": [c.gene or c.transcript_id or c.protein_id
                                    for c in cds_matches][:10],
    }

    if kind == "feature_start":
        base["plasmid_position"] = cds.start
        base["ok"] = True
        base["note"] = "Genomic coord — 0-indexed half-open lower bound of the CDS span."
        return base
    if kind == "feature_end":
        base["plasmid_position"] = cds.end
        base["ok"] = True
        base["note"] = "Genomic coord — 0-indexed half-open upper bound of the CDS span."
        return base
    if kind == "bp_offset":
        base["plasmid_position"] = (cds.start + offset) if sense == "+" else (cds.end - offset - 1)
        base["ok"] = True
        base["note"] = "Genomic-span-relative bp offset (NOT spliced-CDS relative)."
        return base
    if kind != "aa_residue":
        return {"ok": False, "error": f"unknown kind: {kind!r}"}

    # aa_residue path
    if offset < 1:
        return {"ok": False, "error": "aa_residue offset must be >= 1"}
    if offset > len(translation):
        n_alt = len(alternative_isoforms) - 1
        alt_msg = ""
        if n_alt > 0:
            longest_alt = max((iso for iso in alternative_isoforms
                                  if iso["protein_id"] != cds.protein_id),
                                key=lambda r: r["length_aa"], default=None)
            if longest_alt:
                alt_msg = (
                    f". {n_alt} alternative isoform{'s' if n_alt > 1 else ''} also exist; "
                    f"the longest alternative is {longest_alt["protein_id"]} "
                    f"({longest_alt["length_aa"]} aa) — but it does not extend "
                    f"residue {offset} either if length < {offset}"
                )
        return {
            "ok": False,
            "error": (
                f"residue {offset} out of range — the picked isoform "
                f"{cds.protein_id or cds.transcript_id} is {len(translation)} aa "
                f"(valid residues 1..{len(translation)}){alt_msg}. "
                "Check the residue number against the protein, or pass "
                "feature_name=<protein_id> to target a specific isoform."
            ),
            **base,
        }

    cumulative = (offset - 1) * 3
    genomic_pos, codon_bases, spans_intron, lowest_plus, highest_plus = (
        _walk_intervals_to_codon(
            cds.intervals, cds.strand, att.sequence, cumulative,
        )
    )
    amino_acid = translation[offset - 1] if translation else None

    edit_start_1based = (lowest_plus + 1) if lowest_plus is not None else None
    edit_end_1based = (highest_plus + 1) if highest_plus is not None else None
    plus_strand_ref = None
    if (not spans_intron) and lowest_plus is not None and highest_plus is not None:
        plus_strand_ref = att.sequence[lowest_plus:highest_plus + 1].upper()

    if spans_intron:
        note = ("Codon spans an exon-exon boundary in the spliced CDS - "
                "design_pegrnas cannot edit across an intron. Refuse and "
                "explain to the user.")
    else:
        if sense == "-":
            alt_hint = ("For design_pegrnas / design_guides: pass edit_start_1based "
                        "and edit_end_1based directly. This is a - strand gene, so "
                        "`alt` must be the + strand sequence = revcomp(desired_sense_codon). "
                        f"E.g. to install Cys here: sense codon TGT -> alt='ACA'; sense TGC -> alt='GCA'.")
        else:
            alt_hint = ("For design_pegrnas / design_guides: pass edit_start_1based "
                        "and edit_end_1based directly. This is a + strand gene, so "
                        "`alt` is the desired codon directly (no revcomp).")
        note = (f"+ strand codon at 1-based [{edit_start_1based}..{edit_end_1based}] = "
                f"{plus_strand_ref}; sense codon = {codon_bases}. " + alt_hint)

    base.update({
        "ok": True,
        "plasmid_position": lowest_plus if lowest_plus is not None else genomic_pos,
        "edit_start_1based": edit_start_1based,
        "edit_end_1based": edit_end_1based,
        "plus_strand_ref": plus_strand_ref,
        "codon": codon_bases,
        "amino_acid": amino_acid,
        "spans_intron": spans_intron,
        "note": note,
    })
    return base


async def resolve_feature_position(
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,  # unused; kept for dispatch_with_emitters signature
) -> dict[str, Any]:
    aid = args.get("attachment_id")
    feature_name = (args.get("feature_name") or "").strip()
    kind = args.get("kind", "feature_start")
    offset = int(args.get("offset") or 0)

    att = registry.get(aid) if aid else None
    if att is None:
        return {"ok": False, "error": f"unknown attachment_id: {aid!r}"}
    if not feature_name:
        return {"ok": False, "error": "feature_name required"}

    # Branch on attachment kind. Default = plasmid (KB-driven).
    file_kind = attachment_kinds.get_kind(aid)
    if file_kind is not None and file_kind.kind == "genomic":
        return await _resolve_genomic(args, registry)

    # Annotate (cache-friendly).
    try:
        from splicify_api.annotation_cache import annotate_cached
        ann = await annotate_cached(att.sequence, circular=att.circular, depth="full")
    except Exception as e:
        return {"ok": False, "error": f"annotation failed: {type(e).__name__}: {e}"}

    features = ann.get("annotations") or ann.get("features") or []
    feat, candidates = _find_feature(features, feature_name)
    if feat is None:
        return {
            "ok": False,
            "error": f"no feature matching {feature_name!r}",
            "available_features": candidates,
        }

    # Position fields can come from either {qstart, qend} or {start, end}.
    start = feat.get("qstart") if feat.get("qstart") is not None else feat.get("start", 0)
    end = feat.get("qend") if feat.get("qend") is not None else feat.get("end", 0)
    strand = feat.get("strand", 1)
    sense = "+" if strand in (1, "+", "+1", "1") else "-"

    plasmid_position: Optional[int] = None
    codon: Optional[str] = None
    amino_acid: Optional[str] = None
    note = ""

    if kind == "feature_start":
        plasmid_position = start
    elif kind == "feature_end":
        plasmid_position = end
    elif kind == "bp_offset":
        plasmid_position = (start + offset) if sense == "+" else (end - offset - 1)
    elif kind == "aa_residue":
        if offset < 1:
            return {"ok": False, "error": "aa_residue offset must be >= 1"}
        codon_offset = (offset - 1) * 3
        seq = att.sequence
        if sense == "+":
            codon_start = start + codon_offset
            raw = _slice_circular(seq, codon_start, 3).upper()
            codon = raw
            plasmid_position = codon_start
        else:
            codon_start = end - codon_offset - 3
            raw = _slice_circular(seq, codon_start, 3).upper()
            codon = _rc(raw)
            plasmid_position = codon_start
            note = "Plasmid coordinate is the 5'-most base of the codon on the + strand; "
            note += "the codon string above is the sense (CDS-reading) codon."
        amino_acid = _CODON.get(codon.upper(), "?")
    else:
        return {"ok": False, "error": f"unknown kind: {kind!r}"}

    edit_start_1based = None
    edit_end_1based = None
    plus_strand_ref = None
    if kind == "aa_residue" and plasmid_position is not None:
        edit_start_1based = plasmid_position + 1
        edit_end_1based = plasmid_position + 3
        plus_strand_ref = _slice_circular(att.sequence, plasmid_position, 3).upper()
        if sense == "-":
            alt_hint = ("For design_pegrnas / design_guides: pass edit_start_1based "
                        "and edit_end_1based directly. This is a - strand feature, so "
                        "`alt` must be the + strand sequence = revcomp(desired_sense_codon).")
        else:
            alt_hint = ("For design_pegrnas / design_guides: pass edit_start_1based "
                        "and edit_end_1based directly. This is a + strand feature, so "
                        "`alt` is the desired codon directly.")
        note = (f"+ strand codon at 1-based [{edit_start_1based}..{edit_end_1based}] = "
                f"{plus_strand_ref}; sense codon = {codon}. " + alt_hint)

    return {
        "ok": True,
        "feature_name": feat.get("Feature") or feat.get("name"),
        "feature_start": start,
        "feature_end": end,
        "feature_strand": sense,
        "feature_length_bp": int(end) - int(start),
        "plasmid_position": plasmid_position,
        "edit_start_1based": edit_start_1based,
        "edit_end_1based": edit_end_1based,
        "plus_strand_ref": plus_strand_ref,
        "codon": codon,
        "amino_acid": amino_acid,
        "note": note,
        "candidates_considered": candidates[:10],
    }
