from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from . import _data
from typing import Any, Dict, List, Optional, Sequence, Tuple

_MODULE_LIBRARY_ROOT = _data.data_path("Module_Library_gb")
_PRIMER_CLAMP = "GGTGGT"
_PREFERRED_ENZYMES = [
    "EcoRI",
    "HindIII",
    "BamHI",
    "XhoI",
    "NheI",
    "KpnI",
    "SacI",
    "PstI",
    "SalI",
    "AgeI",
    "MluI",
    "BglII",
    "XbaI",
    "SpeI",
    "NotI",
    "AscI",
]
_MCS_TERMS = ("multiple cloning site", "mcs")
_SCREEN_TERMS = {
    "lacz": "Blue-white screening is available because the insert disrupts lacZ-alpha; recombinant colonies should be white on IPTG/X-gal plates.",
    "ccdb": "Negative selection is available because the insert replaces or disrupts ccdB; non-recombinant plasmids remain toxic in standard cloning strains.",
    "sacb": "Counter-selection is available because the insert disrupts sacB; sucrose selection can enrich recombinants if that cassette is active.",
    "gale": "Counter-selection is available because the insert disrupts galE, enabling screening in the appropriate host background.",
    "gfp": "A fluorescent screen is possible because the insert disrupts a GFP-family reporter in the cloning window.",
    "egfp": "A fluorescent screen is possible because the insert disrupts a GFP-family reporter in the cloning window.",
    "zsgreen": "A fluorescent screen is possible because the insert disrupts a GFP-family reporter in the cloning window.",
    "mcherry": "A fluorescent screen is possible because the insert disrupts a fluorescent reporter in the cloning window.",
    "mrfp": "A fluorescent screen is possible because the insert disrupts a fluorescent reporter in the cloning window.",
    "rfp": "A fluorescent screen is possible because the insert disrupts a fluorescent reporter in the cloning window.",
}
_INSERT_ROLES = {"transgene", "reporter", "nuclease"}
_NON_INSERT_ROLES = {
    "backbone",
    "promoter",
    "polya",
    "selection_marker",
    "origin",
    "guide_cassette",
}


@dataclass(frozen=True)
class RestrictionEnzyme:
    name: str
    recognition_seq: str
    overhang_seq: str
    overhang_type: str
    left_of_cut: int
    right_of_cut: int
    buffer: str = "CutSmart"


RE_DATABASE: Dict[str, RestrictionEnzyme] = {
    "EcoRI": RestrictionEnzyme("EcoRI", "GAATTC", "AATT", "5prime", 1, 5),
    "BamHI": RestrictionEnzyme("BamHI", "GGATCC", "GATC", "5prime", 1, 5),
    "HindIII": RestrictionEnzyme("HindIII", "AAGCTT", "AGCT", "5prime", 1, 5),
    "NotI": RestrictionEnzyme("NotI", "GCGGCCGC", "GGCC", "5prime", 2, 6),
    "XhoI": RestrictionEnzyme("XhoI", "CTCGAG", "TCGA", "5prime", 1, 5),
    "NheI": RestrictionEnzyme("NheI", "GCTAGC", "CTAG", "5prime", 1, 5),
    "XbaI": RestrictionEnzyme("XbaI", "TCTAGA", "CTAG", "5prime", 1, 5),
    "SpeI": RestrictionEnzyme("SpeI", "ACTAGT", "CTAG", "5prime", 1, 5),
    "SalI": RestrictionEnzyme("SalI", "GTCGAC", "TCGA", "5prime", 1, 5),
    "KpnI": RestrictionEnzyme("KpnI", "GGTACC", "GTAC", "3prime", 5, 1),
    "BglII": RestrictionEnzyme("BglII", "AGATCT", "GATC", "5prime", 1, 5),
    "AgeI": RestrictionEnzyme("AgeI", "ACCGGT", "CCGG", "5prime", 1, 5),
    "PstI": RestrictionEnzyme("PstI", "CTGCAG", "TGCA", "3prime", 5, 1),
    "MluI": RestrictionEnzyme("MluI", "ACGCGT", "CGCG", "5prime", 1, 5),
    "ClaI": RestrictionEnzyme("ClaI", "ATCGAT", "CGAT", "5prime", 2, 4),
    "AscI": RestrictionEnzyme("AscI", "GGCGCGCC", "CGCG", "5prime", 2, 6),
    "SacI": RestrictionEnzyme("SacI", "GAGCTC", "AGCT", "5prime", 1, 5),
}


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return sequence.translate(table)[::-1]


def maybe_build_restriction_cloning_design(
    design_spec: Dict[str, Any],
    resolved_modules: List[Dict[str, Any]],
    source_plasmid: Optional[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    if (design_spec.get("assembly_strategy") or "").lower() != "restriction_cloning":
        return None
    if not source_plasmid:
        return None

    record = _load_source_record(source_plasmid)
    if not record:
        return None

    insert_mod = _pick_insert_module(resolved_modules)
    if not insert_mod:
        return None

    insert_seq = re.sub(r"[^ACGT]", "", (insert_mod.get("sequence") or "").upper())
    if len(insert_seq) < 30:
        return None

    vector_seq = record["sequence"]
    vector_features = record["features"]
    window = _pick_cloning_window(vector_features)
    if not window:
        return None

    fusion_target = _find_fusion_target(vector_features, window)
    frame_sensitive = bool(fusion_target and len(insert_seq) % 3 == 0)
    enzyme_pair = _pick_vector_enzyme_pair(vector_seq, insert_seq, window, frame_sensitive=frame_sensitive)
    if not enzyme_pair:
        return None
    left_enzyme, right_enzyme, left_start, right_start = enzyme_pair
    replacement_bp = right_start + len(right_enzyme.recognition_seq) - left_start

    oriented_insert_seq = reverse_complement(insert_seq) if frame_sensitive else insert_seq

    assembled_sequence, insert_start = _assemble_sequence(
        vector_seq, oriented_insert_seq, left_enzyme, right_enzyme, left_start, right_start
    )
    delta = len(assembled_sequence) - len(vector_seq)
    insert_end = insert_start + len(left_enzyme.recognition_seq) + len(oriented_insert_seq) + len(right_enzyme.recognition_seq)

    annotations = _remap_features(
        vector_features,
        replaced_start_1b=left_start + 1,
        replaced_end_1b=right_start + len(right_enzyme.recognition_seq),
        delta=delta,
    )
    annotations.extend([
        {
            "name": f"{insert_mod.get('description') or 'Insert'} (restriction fragment)",
            "start": insert_start,
            "end": insert_end,
            "direction": 1,
            "color": "#77c3a2",
        },
        {
            "name": insert_mod.get("description") or "Insert CDS",
            "start": insert_start + len(left_enzyme.recognition_seq),
            "end": insert_start + len(left_enzyme.recognition_seq) + len(oriented_insert_seq),
            "direction": -1 if frame_sensitive else 1,
            "color": "#8de0b7",
        },
        {
            "name": f"{left_enzyme.name} sticky end ligated",
            "start": insert_start + left_enzyme.left_of_cut,
            "end": insert_start + left_enzyme.right_of_cut,
            "direction": 1,
            "color": "#f7d38a",
        },
        {
            "name": f"{right_enzyme.name} sticky end ligated",
            "start": insert_end - len(right_enzyme.recognition_seq) + right_enzyme.left_of_cut,
            "end": insert_end - len(right_enzyme.recognition_seq) + right_enzyme.right_of_cut,
            "direction": 1,
            "color": "#f6b39c",
        },
    ])

    disrupted = _find_disrupted_features(
        vector_features,
        replaced_start_1b=left_start + 1,
        replaced_end_1b=right_start + len(right_enzyme.recognition_seq),
    )
    screening_note = _screening_note(disrupted)
    forward_anneal, forward_tm = _pick_best_anneal(insert_seq, reverse=False)
    reverse_anneal, reverse_tm = _pick_best_anneal(insert_seq, reverse=True)
    if frame_sensitive:
        forward_primer = f"{_PRIMER_CLAMP}{left_enzyme.recognition_seq}{reverse_anneal}"
        reverse_primer = f"{_PRIMER_CLAMP}{right_enzyme.recognition_seq}{forward_anneal}"
    else:
        forward_primer = f"{_PRIMER_CLAMP}{left_enzyme.recognition_seq}{forward_anneal}"
        reverse_primer = f"{_PRIMER_CLAMP}{right_enzyme.recognition_seq}{reverse_anneal}"

    source_name = source_plasmid.get("filename") or source_plasmid.get("source_group") or "source backbone"
    insert_name = insert_mod.get("description") or "insert"

    reply_lines = [
        f"Restriction-cloning workflow: insert {insert_name} into {source_name} by opening the backbone across the cloning window and ligating in a directionally digested insert.",
        "",
        (
            f"Selected enzymes: {left_enzyme.name} on the 5' side and {right_enzyme.name} on the 3' side. "
            f"Both are unique in the backbone cloning window, absent from the insert sequence, and preserve directional cloning with sticky ends."
        ),
        (
            f"Backbone cut plan: digest the vector at {left_enzyme.name} ({left_start + 1}..{left_start + len(left_enzyme.recognition_seq)}) "
            f"and {right_enzyme.name} ({right_start + 1}..{right_start + len(right_enzyme.recognition_seq)}), removing {replacement_bp} bp from the native cloning window before ligation."
        ),
        (
            f"Primer design: amplify the insert with 5' tails `{_PRIMER_CLAMP}{left_enzyme.recognition_seq}` and "
            f"`{_PRIMER_CLAMP}{right_enzyme.recognition_seq}`. Final primers are forward {forward_primer} "
            f"(anneal {len(forward_anneal)} nt, Tm {forward_tm:.1f} C) and reverse {reverse_primer} "
            f"(anneal {len(reverse_anneal)} nt, Tm {reverse_tm:.1f} C)."
        ),
        (
            f"Final plasmid logic: the assembled construct keeps the backbone context, replaces the MCS payload with the insert, "
            f"and leaves the insert bracketed by regenerated {left_enzyme.name} and {right_enzyme.name} recognition sites for downstream verification."
        ),
    ]
    if frame_sensitive and fusion_target:
        reply_lines.append(
            f"{fusion_target} is translated from the reverse strand through the MCS, so the insert was reverse-complemented and aligned at a multiple-of-three offset from the backbone start codon."
        )
    if screening_note:
        reply_lines.append(screening_note)
    if disrupted:
        reply_lines.append(
            "Disrupted modules in the cloning window: " + ", ".join(sorted(disrupted)) + "."
        )

    return {
        "reply": "\n".join(reply_lines),
        "viz": {
            "type": "design",
            "title": f"Restriction Cloning: {source_name} + {insert_name}",
            "sequence": assembled_sequence,
            "topology": "circular",
            "total_length": len(assembled_sequence),
            "annotations": _dedupe_annotations(annotations),
            "restriction_sites": [
                {
                    "name": left_enzyme.name,
                    "start": insert_start,
                    "end": insert_start + len(left_enzyme.recognition_seq),
                    "direction": 1,
                    "re_type": "type2_re",
                    "recognition_seq": left_enzyme.recognition_seq,
                    "color": "#f0b35f",
                },
                {
                    "name": right_enzyme.name,
                    "start": insert_end - len(right_enzyme.recognition_seq),
                    "end": insert_end,
                    "direction": 1,
                    "re_type": "type2_re",
                    "recognition_seq": right_enzyme.recognition_seq,
                    "color": "#ef8d6d",
                },
            ],
        },
        "metadata": {
            "source_backbone": source_name,
            "insert_name": insert_name,
            "left_enzyme": left_enzyme.name,
            "right_enzyme": right_enzyme.name,
            "disrupted_features": sorted(disrupted),
        },
    }


def _pick_insert_module(resolved_modules: Sequence[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    candidates: List[Dict[str, Any]] = []
    for mod in resolved_modules:
        seq = re.sub(r"[^ACGT]", "", (mod.get("sequence") or "").upper())
        if not seq:
            continue
        if mod.get("source_type") == "base_plasmid":
            continue
        role = mod.get("role") or "other"
        if role in _NON_INSERT_ROLES:
            continue
        score = len(seq)
        if role in _INSERT_ROLES:
            score += 5000
        candidates.append({"score": score, "module": mod})
    if not candidates:
        return None
    candidates.sort(key=lambda item: item["score"], reverse=True)
    return candidates[0]["module"]


def _pick_cloning_window(features: Sequence[Dict[str, Any]]) -> Optional[Tuple[int, int]]:
    windows: List[Tuple[int, int]] = []
    for feature in features:
        label = _feature_label(feature).lower()
        if not any(term in label for term in _MCS_TERMS):
            continue
        spans = _feature_spans(feature)
        if not spans:
            continue
        start = min(a for a, _ in spans) - 1
        end = max(b for _, b in spans)
        windows.append((start, end))
    if not windows:
        return None
    windows.sort(key=lambda item: (item[0], item[1]))
    return windows[0]


def _pick_vector_enzyme_pair(
    vector_seq: str,
    insert_seq: str,
    window: Tuple[int, int],
    frame_sensitive: bool = False,
) -> Optional[Tuple[RestrictionEnzyme, RestrictionEnzyme, int, int]]:
    window_start, window_end = window
    candidates: List[Tuple[int, RestrictionEnzyme, RestrictionEnzyme, int, int]] = []
    for left_idx, left_name in enumerate(_PREFERRED_ENZYMES):
        left_enz = RE_DATABASE.get(left_name)
        if not left_enz:
            continue
        if _count_sites(insert_seq, left_enz) > 0:
            continue
        left_positions = _site_positions(vector_seq, left_enz.recognition_seq)
        if len(left_positions) != 1:
            continue
        left_start = left_positions[0]
        if not (window_start <= left_start < window_end):
            continue

        for right_idx, right_name in enumerate(_PREFERRED_ENZYMES):
            if right_name == left_name:
                continue
            right_enz = RE_DATABASE.get(right_name)
            if not right_enz:
                continue
            if _count_sites(insert_seq, right_enz) > 0:
                continue
            right_positions = _site_positions(vector_seq, right_enz.recognition_seq)
            if len(right_positions) != 1:
                continue
            right_start = right_positions[0]
            if not (window_start <= right_start < window_end):
                continue
            if left_start >= right_start:
                continue
            if left_enz.overhang_type == "blunt" or right_enz.overhang_type == "blunt":
                sticky_penalty = 10
            else:
                sticky_penalty = 0
            if left_enz.overhang_seq == right_enz.overhang_seq:
                direction_penalty = 8
            else:
                direction_penalty = 0
            replacement_bp = right_start + len(right_enz.recognition_seq) - left_start
            frame_penalty = 0
            if frame_sensitive and replacement_bp % 3 != 0:
                frame_penalty = 25
            score = (
                abs(left_idx - right_idx)
                + sticky_penalty
                + direction_penalty
                + (0 if left_enz.buffer == right_enz.buffer else 4)
                + frame_penalty
            )
            candidates.append((score, left_enz, right_enz, left_start, right_start))

    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0])
    _, left_enz, right_enz, left_start, right_start = candidates[0]
    return left_enz, right_enz, left_start, right_start


def _find_fusion_target(
    features: Sequence[Dict[str, Any]],
    window: Tuple[int, int],
) -> Optional[str]:
    window_start, window_end = window
    for feature in features:
        label = _feature_label(feature).lower()
        if "lacz-alpha" not in label:
            continue
        for start, end in _feature_spans(feature):
            if start - 1 <= window_end and end >= window_start + 1:
                return "lacZ-alpha"
    return None


def _assemble_sequence(
    vector_seq: str,
    insert_seq: str,
    left_enzyme: RestrictionEnzyme,
    right_enzyme: RestrictionEnzyme,
    left_start: int,
    right_start: int,
) -> Tuple[str, int]:
    left_site = left_enzyme.recognition_seq
    right_site = right_enzyme.recognition_seq
    assembled = (
        f"{vector_seq[:left_start]}{left_site}{insert_seq}{right_site}"
        f"{vector_seq[right_start + len(right_site):]}"
    )
    return assembled, left_start


def _load_source_record(source_plasmid: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    filename = source_plasmid.get("filename")
    if not filename:
        return None
    matches = list(_MODULE_LIBRARY_ROOT.rglob(filename))
    if not matches:
        return None
    return _load_genbank_record(matches[0])


def _parse_genbank_text(text: str) -> Dict[str, Any]:
    """Parse a GenBank text block into the designer's internal dict shape.
    Returns {"sequence": str, "features": List[feature_dict]} with features
    in the shape expected by _pick_cloning_window, _pick_vector_enzyme_pair, etc.
    """
    sequence = ""
    features: List[Dict[str, Any]] = []
    in_origin = False
    in_features = False
    current: Optional[Dict[str, Any]] = None
    current_key: Optional[str] = None

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if line.startswith("ORIGIN"):
            in_origin = True
            in_features = False
            if current:
                features.append(current)
                current = None
            continue
        if line.startswith("//"):
            if current:
                features.append(current)
            break

        if in_origin:
            match = re.match(r"^\s*\d+\s+(.+)$", line)
            if match:
                sequence += re.sub(r"\s+", "", match.group(1)).upper()
            continue

        if not in_features:
            continue

        feature_match = re.match(r"^     (\S+)\s+(.+)$", line)
        if feature_match and not line.startswith("                     /"):
            if current:
                features.append(current)
            current = {
                "type": feature_match.group(1),
                "location": feature_match.group(2).strip(),
                "qualifiers": {},
            }
            current_key = None
            continue

        qualifier_match = re.match(r'^                     /([^=]+)=?"?(.*?)(?:"?)?$', line)
        if qualifier_match and current:
            key = qualifier_match.group(1).strip()
            value = qualifier_match.group(2).strip().strip('"')
            current["qualifiers"][key] = value
            current_key = key
            continue

        continuation_match = re.match(r"^ {21}(.+)$", line)
        if continuation_match and current and current_key:
            prev_val = current["qualifiers"].get(current_key, "")
            cont_val = continuation_match.group(1).strip().strip('"')
            current["qualifiers"][current_key] = f"{prev_val} {cont_val}".strip()

    return {"sequence": sequence, "features": features}


def design_restriction_cloning(
    vector_gb_text: str,
    insert_seq: str,
    insert_name: str,
    vector_name: str = "vector",
    enzyme_override: Optional[Tuple[str, str]] = None,
) -> Dict[str, Any]:
    """Chat-facing entry point: design a restriction-cloning strategy for the
    given vector + insert. Returns {reply, viz, metadata} or
    {error: str} on failure reasons.
    """
    record = _parse_genbank_text(vector_gb_text)
    vector_seq = record["sequence"]
    vector_features = record["features"]
    if not vector_seq:
        return {"error": "Could not parse vector sequence from GenBank file."}

    insert_seq = re.sub(r"[^ACGT]", "", (insert_seq or "").upper())
    if len(insert_seq) < 30:
        return {"error": f"Insert sequence is too short ({len(insert_seq)} bp); need >=30 bp."}

    window = _pick_cloning_window(vector_features)
    if not window:
        return {"error": "No multiple cloning site (MCS) feature found in the vector annotations."}

    fusion_target = _find_fusion_target(vector_features, window)
    frame_sensitive = bool(fusion_target and len(insert_seq) % 3 == 0)

    enzyme_pair = None
    if enzyme_override:
        left_name, right_name = enzyme_override
        left_enz = RE_DATABASE.get(left_name)
        right_enz = RE_DATABASE.get(right_name)
        if not left_enz or not right_enz:
            known = ", ".join(sorted(RE_DATABASE))
            return {"error": f"Unknown enzyme in override: {enzyme_override}. Known: {known}"}
        left_positions = _site_positions(vector_seq, left_enz.recognition_seq)
        right_positions = _site_positions(vector_seq, right_enz.recognition_seq)
        if len(left_positions) != 1 or len(right_positions) != 1:
            return {"error": (
                f"Enzyme override requires single cut sites in the vector. "
                f"{left_name}: {len(left_positions)} site(s); {right_name}: {len(right_positions)} site(s)."
            )}
        if _count_sites(insert_seq, left_enz) > 0 or _count_sites(insert_seq, right_enz) > 0:
            return {"error": f"Override enzymes cut the insert; pick enzymes absent from the insert sequence."}
        left_start = left_positions[0]
        right_start = right_positions[0]
        if left_start >= right_start:
            left_start, right_start = right_start, left_start
            left_enz, right_enz = right_enz, left_enz
        enzyme_pair = (left_enz, right_enz, left_start, right_start)
    else:
        enzyme_pair = _pick_vector_enzyme_pair(
            vector_seq, insert_seq, window, frame_sensitive=frame_sensitive
        )

    if not enzyme_pair:
        return {"error": (
            "No compatible enzyme pair found. Pick two enzymes from the MCS that are "
            "unique in the vector and absent from the insert."
        )}
    left_enzyme, right_enzyme, left_start, right_start = enzyme_pair
    replacement_bp = right_start + len(right_enzyme.recognition_seq) - left_start

    oriented_insert_seq = reverse_complement(insert_seq) if frame_sensitive else insert_seq
    assembled_sequence, insert_start = _assemble_sequence(
        vector_seq, oriented_insert_seq, left_enzyme, right_enzyme, left_start, right_start
    )
    delta = len(assembled_sequence) - len(vector_seq)
    insert_end = insert_start + len(left_enzyme.recognition_seq) + len(oriented_insert_seq) + len(right_enzyme.recognition_seq)

    annotations = _remap_features(
        vector_features,
        replaced_start_1b=left_start + 1,
        replaced_end_1b=right_start + len(right_enzyme.recognition_seq),
        delta=delta,
    )
    annotations.extend([
        {
            "name": f"{insert_name} (restriction fragment)",
            "start": insert_start,
            "end": insert_end,
            "direction": 1,
            "color": "#77c3a2",
        },
        {
            "name": insert_name,
            "start": insert_start + len(left_enzyme.recognition_seq),
            "end": insert_start + len(left_enzyme.recognition_seq) + len(oriented_insert_seq),
            "direction": -1 if frame_sensitive else 1,
            "color": "#8de0b7",
        },
        {
            "name": f"{left_enzyme.name} sticky end ligated",
            "start": insert_start + left_enzyme.left_of_cut,
            "end": insert_start + left_enzyme.right_of_cut,
            "direction": 1,
            "color": "#f7d38a",
        },
        {
            "name": f"{right_enzyme.name} sticky end ligated",
            "start": insert_end - len(right_enzyme.recognition_seq) + right_enzyme.left_of_cut,
            "end": insert_end - len(right_enzyme.recognition_seq) + right_enzyme.right_of_cut,
            "direction": 1,
            "color": "#f6b39c",
        },
    ])

    disrupted = _find_disrupted_features(
        vector_features,
        replaced_start_1b=left_start + 1,
        replaced_end_1b=right_start + len(right_enzyme.recognition_seq),
    )
    screening_note = _screening_note(disrupted)

    forward_anneal, forward_tm = _pick_best_anneal(insert_seq, reverse=False)
    reverse_anneal, reverse_tm = _pick_best_anneal(insert_seq, reverse=True)
    if frame_sensitive:
        forward_primer = f"{_PRIMER_CLAMP}{left_enzyme.recognition_seq}{reverse_anneal}"
        reverse_primer = f"{_PRIMER_CLAMP}{right_enzyme.recognition_seq}{forward_anneal}"
    else:
        forward_primer = f"{_PRIMER_CLAMP}{left_enzyme.recognition_seq}{forward_anneal}"
        reverse_primer = f"{_PRIMER_CLAMP}{right_enzyme.recognition_seq}{reverse_anneal}"

    reply_lines = [
        f"Restriction cloning: insert {insert_name} into {vector_name} across the MCS.",
        "",
        (
            f"Enzymes: {left_enzyme.name} (5' side) and {right_enzyme.name} (3' side). "
            f"Both unique in the backbone MCS and absent from the insert sequence."
        ),
        (
            f"Backbone cut plan: digest at {left_enzyme.name} ({left_start + 1}..{left_start + len(left_enzyme.recognition_seq)}) "
            f"and {right_enzyme.name} ({right_start + 1}..{right_start + len(right_enzyme.recognition_seq)}), removing {replacement_bp} bp."
        ),
        (
            f"Insert PCR primers:"
            f"\n  Forward: {forward_primer} (anneal {len(forward_anneal)} nt, Tm {forward_tm:.1f} C)"
            f"\n  Reverse: {reverse_primer} (anneal {len(reverse_anneal)} nt, Tm {reverse_tm:.1f} C)"
        ),
        (
            f"Final plasmid: {len(assembled_sequence)} bp, insert bracketed by regenerated "
            f"{left_enzyme.name}/{right_enzyme.name} sites for diagnostic digestion."
        ),
    ]
    if frame_sensitive and fusion_target:
        reply_lines.append(
            f"{fusion_target} is on the reverse strand through the MCS; insert was reverse-complemented "
            f"and aligned at a mod-3 offset from the backbone start codon."
        )
    if screening_note:
        reply_lines.append(screening_note)
    if disrupted:
        reply_lines.append("Disrupted features in the cloning window: " + ", ".join(sorted(disrupted)) + ".")

    return {
        "reply": "\n".join(reply_lines),
        "viz": {
            "type": "design",
            "title": f"Restriction Cloning: {vector_name} + {insert_name}",
            "sequence": assembled_sequence,
            "topology": "circular",
            "total_length": len(assembled_sequence),
            "annotations": _dedupe_annotations(annotations),
            "restriction_sites": [
                {
                    "name": left_enzyme.name,
                    "start": insert_start,
                    "end": insert_start + len(left_enzyme.recognition_seq),
                    "direction": 1,
                    "re_type": "type2_re",
                    "recognition_seq": left_enzyme.recognition_seq,
                    "color": "#f0b35f",
                },
                {
                    "name": right_enzyme.name,
                    "start": insert_end - len(right_enzyme.recognition_seq),
                    "end": insert_end,
                    "direction": 1,
                    "re_type": "type2_re",
                    "recognition_seq": right_enzyme.recognition_seq,
                    "color": "#ef8d6d",
                },
            ],
        },
        "metadata": {
            "vector_name": vector_name,
            "insert_name": insert_name,
            "left_enzyme": left_enzyme.name,
            "right_enzyme": right_enzyme.name,
            "forward_primer": forward_primer,
            "reverse_primer": reverse_primer,
            "forward_anneal_tm": forward_tm,
            "reverse_anneal_tm": reverse_tm,
            "assembled_sequence": assembled_sequence,
            "insert_sequence": oriented_insert_seq,
            "insert_start": insert_start,
            "insert_end": insert_end,
            "disrupted_features": sorted(disrupted),
            "frame_sensitive": frame_sensitive,
        },
    }


def _load_genbank_record(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    sequence = ""
    features: List[Dict[str, Any]] = []
    in_origin = False
    in_features = False
    current: Optional[Dict[str, Any]] = None
    current_key: Optional[str] = None

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if line.startswith("ORIGIN"):
            in_origin = True
            in_features = False
            if current:
                features.append(current)
                current = None
            continue
        if line.startswith("//"):
            if current:
                features.append(current)
            break

        if in_origin:
            match = re.match(r"^\s*\d+\s+(.+)$", line)
            if match:
                sequence += re.sub(r"\s+", "", match.group(1)).upper()
            continue

        if not in_features:
            continue

        feature_match = re.match(r"^     (\S+)\s+(.+)$", line)
        if feature_match and not line.startswith("                     /"):
            if current:
                features.append(current)
            current = {
                "type": feature_match.group(1),
                "location": feature_match.group(2).strip(),
                "qualifiers": {},
            }
            current_key = None
            continue

        qualifier_match = re.match(r'^                     /([^=]+)=?"?(.*?)(?:"?)?$', line)
        if qualifier_match and current:
            key = qualifier_match.group(1).strip()
            value = qualifier_match.group(2).strip().strip('"')
            current["qualifiers"][key] = value
            current_key = key
            continue

        continuation_match = re.match(r"^ {21}(.+)$", line)
        if continuation_match and current and current_key:
            prev_val = current["qualifiers"].get(current_key, "")
            cont_val = continuation_match.group(1).strip().strip('"')
            current["qualifiers"][current_key] = f"{prev_val} {cont_val}".strip()

    return {"sequence": sequence, "features": features}


def _feature_spans(feature: Dict[str, Any]) -> List[Tuple[int, int]]:
    return [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", feature.get("location", ""))]


def _feature_label(feature: Dict[str, Any]) -> str:
    qualifiers = feature.get("qualifiers", {})
    return (
        qualifiers.get("label")
        or qualifiers.get("gene")
        or qualifiers.get("product")
        or qualifiers.get("note")
        or feature.get("type")
        or "feature"
    )


def _remap_features(
    features: Sequence[Dict[str, Any]],
    replaced_start_1b: int,
    replaced_end_1b: int,
    delta: int,
) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for feature in features:
        if feature.get("type") == "source":
            continue
        label = _feature_label(feature)
        spans = _feature_spans(feature)
        if not spans:
            continue
        strand = -1 if "complement" in feature.get("location", "") else 1
        for start, end in spans:
            if end < replaced_start_1b:
                mapped = (start, end)
            elif start > replaced_end_1b:
                mapped = (start + delta, end + delta)
            else:
                if start < replaced_start_1b:
                    out.append({
                        "name": str(label),
                        "start": start - 1,
                        "end": replaced_start_1b - 1,
                        "direction": strand,
                        "color": "#86a8d9",
                        "origin": "genbank",
                    })
                if end > replaced_end_1b:
                    out.append({
                        "name": str(label),
                        "start": replaced_end_1b + delta,
                        "end": end + delta,
                        "direction": strand,
                        "color": "#86a8d9",
                        "origin": "genbank",
                    })
                continue
            out.append({
                "name": str(label),
                "start": mapped[0] - 1,
                "end": mapped[1],
                "direction": strand,
                "color": "#86a8d9",
                "origin": "genbank",
            })
    return out


def _find_disrupted_features(
    features: Sequence[Dict[str, Any]],
    replaced_start_1b: int,
    replaced_end_1b: int,
) -> List[str]:
    disrupted: List[str] = []
    for feature in features:
        if feature.get("type") == "source":
            continue
        label = _feature_label(feature)
        for start, end in _feature_spans(feature):
            if start <= replaced_end_1b and end >= replaced_start_1b:
                disrupted.append(label)
                break
    return disrupted


def _screening_note(disrupted_features: Sequence[str]) -> Optional[str]:
    for feature in disrupted_features:
        lowered = feature.lower()
        for key, note in _SCREEN_TERMS.items():
            if key in lowered:
                return note
    return None


def _site_positions(sequence: str, recognition_seq: str) -> List[int]:
    return [match.start() for match in re.finditer(re.escape(recognition_seq), sequence)]


def _count_sites(sequence: str, enzyme: RestrictionEnzyme) -> int:
    if not sequence:
        return 0
    seq = sequence.upper()
    pattern = enzyme.recognition_seq.upper()
    rc = reverse_complement(pattern)
    hits = len(_site_positions(seq, pattern))
    if rc != pattern:
        hits += len(_site_positions(seq, rc))
    return hits


def _pick_best_anneal(template: str, reverse: bool, target_tm: float = 60.0) -> Tuple[str, float]:
    try:
        from .gibson_primers import ThermodynamicCalculator

        thermo = ThermodynamicCalculator()
        seq = reverse_complement(template[-30:]) if reverse else template[:30]
        best_seq = seq[:18]
        best_tm = thermo.calculate_tm(best_seq)
        best_delta = abs(best_tm - target_tm)

        for length in range(18, min(30, len(seq)) + 1):
            candidate = seq[:length]
            tm = thermo.calculate_tm(candidate)
            delta = abs(tm - target_tm)
            if delta < best_delta:
                best_seq = candidate
                best_tm = tm
                best_delta = delta
        return best_seq, best_tm
    except Exception:
        seq = reverse_complement(template[-24:]) if reverse else template[:24]
        candidate = seq[:24]
        gc = candidate.count("G") + candidate.count("C")
        tm = 2 * (len(candidate) - gc) + 4 * gc
        return candidate, float(tm)


def _dedupe_annotations(annotations: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    deduped: List[Dict[str, Any]] = []
    seen = set()
    for ann in annotations:
        key = (
            int(ann.get("start", -1)),
            int(ann.get("end", -1)),
            str(ann.get("name") or "").strip().lower(),
        )
        if key in seen:
            continue
        seen.add(key)
        deduped.append(dict(ann))
    return deduped
