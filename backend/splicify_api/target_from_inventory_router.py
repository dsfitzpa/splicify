"""
Target-from-Inventory Router
============================

Picks the best cloning workflow to construct a target plasmid from an inventory
of source plasmids, using full annotation (features + modules + cloning_features)
to compute per-workflow FeasibilityReports.

Spec: TARGET_FROM_INVENTORY_ROUTING.md

Public surface:
  - FeasibilityReport
  - TargetContext / InventoryContext
  - annotate_one(name, gb_text, sequence) -> Context
  - route(target_ctx, inventory_ctxs) -> (chosen, all_reports)
  - route_from_uploads(target_upload, inventory_uploads) -> dict
  - build_audit_markdown(reports, chosen) -> str

Used by chat.py when intent == "target_from_inventory".
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from fastapi import UploadFile

from .cloning.gateway_sites import scan_att_sites
from .cloning_feature_annotator import scan_cloning_features
from .inv_gib import (
    MAX_KMER_POS_PER_SEED,
    SEED_K,
    _best_exact_hit_for_inventory,
    _build_kmer_index,
    _read_seq_from_upload,
    _revcomp,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class FeasibilityReport:
    workflow: str
    feasible: bool
    score: float = 0.0
    confidence: float = 0.5
    work_estimate: int = 99
    success_estimate: float = 0.0
    rationale: str = ""
    handler_args: Dict[str, Any] = field(default_factory=dict)
    validation_mode: str = "loose"
    warnings: List[str] = field(default_factory=list)


@dataclass
class PlasmidContext:
    name: str
    sequence: str
    gb_text: Optional[str]
    annotations: List[Dict[str, Any]]
    modules: List[Dict[str, Any]]
    cloning: Any  # ScanResult from cloning_feature_annotator


TargetContext = PlasmidContext
InventoryContext = PlasmidContext


# ---------------------------------------------------------------------------
# Annotation helpers
# ---------------------------------------------------------------------------

def _seq_to_minimal_gb(seq: str, name: str) -> str:
    """Wrap a bare sequence in a minimal GenBank record so annotate_genbank can ingest it."""
    name = (name or "seq")[:16].replace(" ", "_") or "seq"
    body_lines = []
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60].lower()
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        body_lines.append(f"{i + 1:>9} {groups}")
    body = "\n".join(body_lines)
    return (
        f"LOCUS       {name:<16} {len(seq)} bp ds-DNA     circular     01-JAN-1970\n"
        f"FEATURES             Location/Qualifiers\n"
        f"ORIGIN\n{body}\n//\n"
    )


async def annotate_one(name: str, gb_text: Optional[str], sequence: str) -> PlasmidContext:
    """Run the standard annotation pipeline + cloning-feature scanner on one plasmid."""
    from .plannotate_router import AnnotateRequest, annotate_genbank

    text = gb_text or _seq_to_minimal_gb(sequence, name)
    req = AnnotateRequest(
        gb_text=text,
        session_id="router",
        options={"linear": False, "detailed": True, "file_name": name},
    )
    try:
        resp = await annotate_genbank(req)
        anns = [a.dict() if hasattr(a, "dict") else dict(a) for a in (resp.annotations or [])]
        modules = list(resp.modules or [])
    except Exception as e:
        logger.warning(f"annotate_one failed for {name}: {e}")
        anns, modules = [], []

    try:
        cloning = scan_cloning_features(sequence)
    except Exception as e:
        logger.warning(f"scan_cloning_features failed for {name}: {e}")
        cloning = None

    return PlasmidContext(
        name=name,
        sequence=sequence.upper(),
        gb_text=gb_text,
        annotations=anns,
        modules=modules,
        cloning=cloning,
    )


# ---------------------------------------------------------------------------
# Homology helpers
# ---------------------------------------------------------------------------

def _best_hit_to_target(target: str, inv_seq: str, inv_name: str = "inv") -> Optional[Any]:
    """Best exact contiguous match of inv_seq onto circular target. Returns Hit or None."""
    if not target or not inv_seq:
        return None
    target2 = target + target
    kidx = _build_kmer_index(target2, SEED_K, MAX_KMER_POS_PER_SEED)
    return _best_exact_hit_for_inventory(
        target=target,
        inv_seq=inv_seq,
        inv_name=inv_name,
        kmer_index=kidx,
        k=SEED_K,
        min_match_bp=50,  # allow short hits for the assessors
        max_total_seed_hits=4000,
    )


def _coverage_of_region(region_seq: str, inventory: List[InventoryContext]) -> Tuple[float, Optional[InventoryContext], int]:
    """For a target region, find inventory plasmid with best coverage. Returns (coverage_frac, best_ctx, match_bp)."""
    best_frac = 0.0
    best_ctx: Optional[InventoryContext] = None
    best_bp = 0
    region_len = max(1, len(region_seq))
    for inv in inventory:
        h = _best_hit_to_target(region_seq, inv.sequence, inv.name)
        if not h:
            continue
        # Clamp — circular doubled-target scan can report match_len > region_len
        # when a long inventory plasmid matches both copies of the target.
        frac = min(1.0, h.match_len / region_len)
        if frac > best_frac:
            best_frac = frac
            best_ctx = inv
            best_bp = h.match_len
    return best_frac, best_ctx, best_bp


# ---------------------------------------------------------------------------
# Per-workflow assessors
# ---------------------------------------------------------------------------

GATEWAY_PAIR_RULES = {
    # (subtype_a, subtype_b) -> (reaction_type, donor_subtypes_needed)
    ("attB1", "attB2"): ("BP", ["attP1", "attP2"]),
    ("attL1", "attL2"): ("LR", ["attR1", "attR2"]),
    ("attR1", "attR2"): ("LR", ["attL1", "attL2"]),
    ("attP1", "attP2"): ("BP", ["attB1", "attB2"]),
}


def _att_pair(target_ctx: TargetContext) -> Optional[Tuple[Any, Any]]:
    if not target_ctx.cloning:
        return None
    atts = [f for f in target_ctx.cloning.features if f.feature_family == "gateway_att"]
    if len(atts) < 2:
        return None
    # Look for any compatible pair by subtype
    for i, a in enumerate(atts):
        for b in atts[i + 1:]:
            key = (a.subtype, b.subtype) if (a.subtype, b.subtype) in GATEWAY_PAIR_RULES else (b.subtype, a.subtype)
            if key in GATEWAY_PAIR_RULES:
                return (a, b) if a.start <= b.start else (b, a)
    return None


def _att_class(subtype: str) -> str:
    """Return 'B', 'L', 'P', 'R', or '?' for an att-site subtype like 'attB1'."""
    s = (subtype or "").lower().replace("att", "")
    if not s:
        return "?"
    return s[0].upper() if s[0].upper() in {"B", "L", "P", "R"} else "?"


def _split_target_by_att(seq: str, a, b) -> Tuple[str, str]:
    """Split circular target by att pair into (cargo, backbone).

    Cargo = side bracketed by attB*/attL* cores (insert / entry side).
    Backbone = side bracketed by attP*/attR* cores (donor / destination side).
    Falls back to (between, around) when both cores are the same class.
    """
    between = seq[a.end:b.start]
    around = seq[b.end:] + seq[:a.start]
    cls_a = _att_class(getattr(a, "subtype", ""))
    cls_b = _att_class(getattr(b, "subtype", ""))
    insert_classes = {"B", "L"}
    vector_classes = {"P", "R"}
    if cls_a in insert_classes and cls_b in insert_classes:
        return between, around  # cargo between, backbone around
    if cls_a in vector_classes and cls_b in vector_classes:
        return around, between  # cargo around, backbone between
    return between, around


def assess_gateway_feasibility(target: TargetContext, inventory: List[InventoryContext]) -> FeasibilityReport:
    """Unified Gateway assessor — covers single BP/LR (att pair) and MultiSite
    (Hellsgate / pDONR P1-Pn family). Returns whichever interpretation scores
    higher. Merged 2026-04-27 from the previous separate gateway + multisite
    assessors so the router has a single Gateway entry."""
    single = _assess_single_site_gateway(target, inventory)
    multi = _assess_multisite_gateway(target, inventory)
    if multi.feasible and (not single.feasible or multi.score > single.score):
        ha = dict(multi.handler_args or {})
        ha["gateway_variant"] = "multisite"
        return FeasibilityReport(
            workflow="gateway_cloning",
            feasible=multi.feasible,
            score=multi.score,
            confidence=multi.confidence,
            work_estimate=multi.work_estimate,
            success_estimate=multi.success_estimate,
            rationale="MultiSite — " + multi.rationale,
            handler_args=ha,
            validation_mode=multi.validation_mode,
        )
    ha = dict(single.handler_args or {})
    ha["gateway_variant"] = "single"
    return FeasibilityReport(
        workflow="gateway_cloning",
        feasible=single.feasible,
        score=single.score,
        confidence=single.confidence,
        work_estimate=single.work_estimate,
        success_estimate=single.success_estimate,
        rationale=single.rationale,
        handler_args=ha,
        validation_mode=single.validation_mode,
    )


def _assess_single_site_gateway(target: TargetContext, inventory: List[InventoryContext]) -> FeasibilityReport:
    pair = _att_pair(target)
    if not pair:
        return FeasibilityReport(
            workflow="gateway_cloning", feasible=False,
            rationale="Target has no compatible gateway_att pair (>=2 att sites required).",
            validation_mode="strict",
        )
    a, b = pair
    seq = target.sequence
    payload, backbone = _split_target_by_att(seq, a, b)

    cargo_frac, cargo_inv, _ = _coverage_of_region(payload, inventory)
    backbone_frac, backbone_inv, _ = _coverage_of_region(backbone, inventory)

    pair_key = (a.subtype, b.subtype) if (a.subtype, b.subtype) in GATEWAY_PAIR_RULES else (b.subtype, a.subtype)
    reaction, needed_donor_sites = GATEWAY_PAIR_RULES[pair_key]

    native_donor: Optional[InventoryContext] = None
    for inv in inventory:
        if not inv.cloning:
            continue
        subtypes = {f.subtype for f in inv.cloning.features if f.feature_family == "gateway_att"}
        if all(s in subtypes for s in needed_donor_sites):
            native_donor = inv
            break

    feasible = (cargo_frac >= 0.5) and (native_donor is not None or backbone_frac >= 0.5)
    success = cargo_frac * (1.0 if native_donor else backbone_frac) * 0.95
    work = 1 + (0 if native_donor else 2)
    score = success * (1.0 / (1.0 + work))

    return FeasibilityReport(
        workflow="gateway_cloning",
        feasible=feasible,
        score=score,
        confidence=0.85,
        work_estimate=work,
        success_estimate=success,
        rationale=(
            f"att pair {a.subtype}+{b.subtype} -> {reaction} reaction; "
            f"cargo coverage={cargo_frac:.2f} from {cargo_inv.name if cargo_inv else 'none'}; "
            f"backbone coverage={backbone_frac:.2f} from {backbone_inv.name if backbone_inv else 'none'}; "
            f"native donor={'yes (' + native_donor.name + ')' if native_donor else 'no -- primers needed'}"
        ),
        handler_args={
            "intent": "gateway_cloning",
            "reaction_type": reaction,
            "donor_name": native_donor.name if native_donor else None,
            "cargo_source": cargo_inv.name if cargo_inv else None,
            "primers_needed": native_donor is None,
        },
        validation_mode="strict",
    )



def assess_inv_gib_feasibility(target: TargetContext, inventory: List[InventoryContext]) -> FeasibilityReport:
    """Wrap inv_gib's exact-match set-cover; gate junctions on primer_design_warning regions."""
    target2 = target.sequence + target.sequence
    kidx = _build_kmer_index(target2, SEED_K, MAX_KMER_POS_PER_SEED)
    hits = []
    for inv in inventory:
        h = _best_hit_to_target(target.sequence, inv.sequence, inv.name)
        if h:
            hits.append((inv, h))
    hits.sort(key=lambda x: x[1].match_len, reverse=True)

    L = len(target.sequence)
    covered = bytearray(L)
    chosen = []
    for inv, h in hits:
        new_cov = sum(1 for i in range(h.tgt_start, min(h.tgt_end, L)) if not covered[i])
        if new_cov >= 50:
            for i in range(h.tgt_start, min(h.tgt_end, L)):
                covered[i] = 1
            chosen.append((inv, h))

    coverage_frac = sum(covered) / max(1, L)
    n_frags = len(chosen)
    n_synth = max(0, _gap_count(covered, min_gap=20))

    # Junction designability: penalise if any junction lands inside primer_design_warning
    warning_intervals = []
    if target.cloning:
        for f in target.cloning.features:
            if f.feature_family == "primer_design_warning":
                warning_intervals.append((f.start, f.end, f.subtype))
    bad_junctions = 0
    junctions = []
    if n_frags >= 2:
        sorted_chosen = sorted(chosen, key=lambda x: x[1].tgt_start)
        for (_, h1), (_, h2) in zip(sorted_chosen, sorted_chosen[1:]):
            j = h1.tgt_end if h1.tgt_end <= L else h1.tgt_end - L
            junctions.append(j)
            for s, e, sub in warning_intervals:
                if s - 25 <= j <= e + 25:
                    bad_junctions += 1
                    break

    junction_quality = 1.0 - (bad_junctions / max(1, len(junctions))) if junctions else 1.0
    success = (0.95 ** n_frags) * coverage_frac * junction_quality
    work = max(1, n_frags * 2)
    feasible = (coverage_frac >= 0.6) and (n_frags >= 1)
    score = success * (1.0 / (1.0 + work))

    return FeasibilityReport(
        workflow="inv_gib",
        feasible=feasible,
        score=score,
        confidence=0.9,
        work_estimate=work,
        success_estimate=success,
        rationale=(
            f"{n_frags} inventory fragments cover {coverage_frac*100:.1f}%; "
            f"{n_synth} synth gaps; {bad_junctions}/{len(junctions)} junctions land in PCR warning regions"
        ),
        handler_args={"intent": "inv_gib"},
        validation_mode="loose",
        warnings=[f"junction near {sub}" for s, e, sub in warning_intervals if any(s - 25 <= j <= e + 25 for j in junctions)],
    )


def _gap_count(covered: bytearray, min_gap: int) -> int:
    n, gaps, run = 0, 0, 0
    for b in covered:
        if not b:
            run += 1
        else:
            if run >= min_gap:
                gaps += 1
            run = 0
    if run >= min_gap:
        gaps += 1
    return gaps


def assess_golden_gate_feasibility(target: TargetContext, inventory: List[InventoryContext]) -> FeasibilityReport:
    """Type IIs non-cutter intersection across the chosen fragments + their source plasmids."""
    if not target.cloning:
        return FeasibilityReport(workflow="golden_gate_primer_design", feasible=False,
                                 rationale="No cloning_features payload on target.")

    # Use inv_gib's set-cover as the fragment partition
    hits = []
    for inv in inventory:
        h = _best_hit_to_target(target.sequence, inv.sequence, inv.name)
        if h:
            hits.append((inv, h))
    hits.sort(key=lambda x: x[1].match_len, reverse=True)

    L = len(target.sequence)
    covered = bytearray(L)
    chosen: List[Tuple[InventoryContext, Any]] = []
    for inv, h in hits:
        new_cov = sum(1 for i in range(h.tgt_start, min(h.tgt_end, L)) if not covered[i])
        if new_cov >= 50:
            for i in range(h.tgt_start, min(h.tgt_end, L)):
                covered[i] = 1
            chosen.append((inv, h))

    if not chosen:
        return FeasibilityReport(workflow="golden_gate_primer_design", feasible=False,
                                 rationale="No inventory fragments match target.")

    # Non-cutter intersection: enzymes that cut neither target nor any chosen source plasmid
    target_non_cutters = set(target.cloning.non_cutters or [])
    iis_default = {"BsaI", "BsmBI", "BbsI", "SapI", "AarI", "Esp3I"}
    candidate = target_non_cutters & iis_default
    for inv, _ in chosen:
        if inv.cloning:
            candidate &= set(inv.cloning.non_cutters or [])
    candidate = candidate & iis_default

    # Tiebreak preference order
    pref = ["BsaI", "BsmBI", "BbsI", "Esp3I", "AarI", "SapI"]
    enzyme = next((e for e in pref if e in candidate), None)

    n_frags = len(chosen)
    feasible = enzyme is not None and n_frags >= 2
    success = (0.97 ** n_frags) if feasible else 0.0
    work = n_frags + 1
    score = success * (1.0 / (1.0 + work))

    return FeasibilityReport(
        workflow="golden_gate_primer_design",
        feasible=feasible,
        score=score,
        confidence=0.75,
        work_estimate=work,
        success_estimate=success,
        rationale=(
            f"{n_frags} fragments; non-cutter Type IIs intersection={sorted(candidate) or 'empty'}; "
            f"chosen enzyme={enzyme or 'none'}"
        ),
        handler_args={
            "intent": "golden_gate_primer_design",
            "enzyme": enzyme,
            "fragments": [inv.name for inv, _ in chosen],
        },
        validation_mode="loose",
    )


def assess_restriction_feasibility(target: TargetContext, inventory: List[InventoryContext]) -> FeasibilityReport:
    """
    Inventory-driven: enumerate Type II enzymes that partition target into 2-4 fragments
    at module-aware boundaries, and check if each fragment can be sourced from inventory.
    """
    if not target.cloning:
        return FeasibilityReport(workflow="restriction_cloning", feasible=False,
                                 rationale="No cloning_features payload on target.")

    L = len(target.sequence)
    cut_count = target.cloning.cut_count_per_enzyme or {}
    enzymes_2_4 = [e for e, n in cut_count.items() if 2 <= n <= 4]
    if not enzymes_2_4:
        # Try single enzyme that cuts twice via 2x scan
        enzymes_2_4 = [e for e, n in cut_count.items() if n == 2]
    if not enzymes_2_4:
        return FeasibilityReport(workflow="restriction_cloning", feasible=False,
                                 rationale="No Type II enzyme cuts target 2-4 times.")

    best_combo = None
    best_score = 0.0
    best_rationale = ""
    best_handler = {}

    # Pick the enzyme producing the highest-coverage partition
    for enz in enzymes_2_4[:8]:
        cuts = sorted([f.cut_profile.cut_top for f in target.cloning.features
                       if f.feature_family == "restriction_site_II" and f.name == enz
                       and f.cut_profile])
        if len(cuts) < 2:
            continue
        # Define fragments as (cut_i, cut_{i+1}) intervals on circular target
        fragments = []
        for i in range(len(cuts)):
            a = cuts[i]
            b = cuts[(i + 1) % len(cuts)]
            if b > a:
                frag = target.sequence[a:b]
            else:
                frag = target.sequence[a:] + target.sequence[:b]
            fragments.append((a, b, frag))

        coverages = []
        sources = []
        for a, b, frag in fragments:
            cov, src, _ = _coverage_of_region(frag, inventory)
            coverages.append(cov)
            sources.append(src.name if src else None)

        avg_cov = sum(coverages) / len(coverages)
        all_sourced = all(c >= 0.85 for c in coverages)
        success = avg_cov * (1.0 if all_sourced else 0.6)
        work = len(fragments) + 1
        s = success * (1.0 / (1.0 + work))
        if s > best_score:
            best_score = s
            best_combo = enz
            best_rationale = f"enzyme {enz} partitions target into {len(fragments)} fragments; sources={sources}; avg cov={avg_cov:.2f}"
            best_handler = {
                "intent": "restriction_cloning",
                "enzymes": [enz],
                "fragments": [{"start": a, "end": b, "source": src} for (a, b, _), src in zip(fragments, sources)],
            }

    feasible = best_combo is not None and best_score > 0.05
    return FeasibilityReport(
        workflow="restriction_cloning",
        feasible=feasible,
        score=best_score,
        confidence=0.7,
        work_estimate=3,
        success_estimate=best_score * 4,  # rough back-conversion
        rationale=best_rationale or "no viable enzyme partition found",
        handler_args=best_handler,
        validation_mode="loose",  # promote to strict if all_sourced
    )


def assess_sdm_feasibility(target: TargetContext, inventory: List[InventoryContext]) -> FeasibilityReport:
    """Find an inventory plasmid identical to target except for one deletion or single ≤40 bp change."""
    L = len(target.sequence)
    best = None
    best_diff = None
    for inv in inventory:
        h = _best_hit_to_target(target.sequence, inv.sequence, inv.name)
        if not h:
            continue
        # One-best-hit must cover almost the whole target
        cov = min(1.0, h.match_len / L)
        if cov < 0.90:
            continue
        # Diff region = target portion not covered by the hit
        diff_start = h.tgt_end if h.tgt_end <= L else h.tgt_end - L
        diff_end = h.tgt_start
        if diff_end < diff_start:
            diff_len = (L - diff_start) + diff_end
        else:
            diff_len = diff_end - diff_start

        # Compare lengths to classify
        len_inv = len(inv.sequence)
        len_diff = abs(len_inv - L)

        if len_inv > L and (len_inv - L) > 0:
            mutation_type = "deletion"  # delete extra bp from inventory to get target
            change_bp = len_inv - L
        elif L > len_inv and (L - len_inv) <= 40:
            mutation_type = "insertion"
            change_bp = L - len_inv
        elif L == len_inv and diff_len <= 40:
            mutation_type = "substitution"
            change_bp = diff_len
        elif diff_len <= 40 and len_diff <= 40:
            mutation_type = "mixed_small"
            change_bp = max(diff_len, len_diff)
        else:
            continue

        # SDM is meaningful only when there IS a change; reject identical pairs.
        if change_bp == 0:
            continue

        if best is None or change_bp < best_diff:
            best = (inv, mutation_type, change_bp, diff_start, diff_end)
            best_diff = change_bp

    if not best:
        return FeasibilityReport(workflow="sdm_design", feasible=False,
                                 rationale="No inventory plasmid is near-identical to target with ≤40 bp change.",
                                 validation_mode="strict")

    inv, mtype, bp, ds, de = best
    success = 0.98 if mtype == "deletion" else max(0.5, 0.98 - (bp / 1000))
    score = success * (1.0 / (1.0 + 1))

    return FeasibilityReport(
        workflow="sdm_design",
        feasible=True,
        score=score,
        confidence=0.85,
        work_estimate=1,
        success_estimate=success,
        rationale=f"{inv.name} matches target with single {mtype} of {bp} bp at ~{ds}-{de}",
        handler_args={
            "intent": "sdm_design",
            "template_source": inv.name,
            "mutation_type": mtype,
            "target_position_start": ds,
            "target_position_end": de,
        },
        validation_mode="strict",
    )


def assess_sgrna_golden_gate_feasibility(target: TargetContext, inventory: List[InventoryContext]) -> FeasibilityReport:
    """Pol III guide cassette in inventory + Type IIs cuts-twice-only enzyme; slot length 17-250 bp."""
    pol3_promoters = {"u6", "h1", "7sk"}
    scaffold_names = {"tracrrna", "grna scaffold", "sgrna scaffold", "pegrna scaffold", "cas9 sgrna scaffold"}

    for inv in inventory:
        # Find pol3 promoter + scaffold within 500 bp of each other
        promoters = [a for a in inv.annotations
                     if any(p in (a.get("name") or "").lower() for p in pol3_promoters)]
        scaffolds = [a for a in inv.annotations
                     if any(s in (a.get("name") or "").lower() for s in scaffold_names)]
        if not promoters or not scaffolds:
            continue

        # Find a paired promoter+scaffold within 500 bp
        cassette = None
        for p in promoters:
            p_loc = _loc_to_range(p.get("location", ""))
            if not p_loc:
                continue
            for sc in scaffolds:
                s_loc = _loc_to_range(sc.get("location", ""))
                if not s_loc:
                    continue
                gap = s_loc[0] - p_loc[1]
                if 0 <= gap <= 500:
                    cassette = (p_loc[0], s_loc[1])
                    break
            if cassette:
                break
        if not cassette:
            continue

        # Find Type IIs enzyme that cuts INSIDE cassette exactly twice and OUTSIDE zero times
        if not inv.cloning:
            continue
        per_enz_cuts: Dict[str, Tuple[int, int]] = {}
        for f in inv.cloning.features:
            if f.feature_family != "restriction_site_IIs":
                continue
            inside = cassette[0] <= f.start <= cassette[1]
            ic, oc = per_enz_cuts.get(f.name, (0, 0))
            per_enz_cuts[f.name] = (ic + (1 if inside else 0), oc + (0 if inside else 1))

        chosen_enz = None
        for enz, (ic, oc) in per_enz_cuts.items():
            if ic == 2 and oc == 0:
                chosen_enz = enz
                break
        if not chosen_enz:
            continue

        # Estimate slot length from inventory cassette flanks aligned to target
        # Simple proxy: target slot = corresponding region in target between same flanks
        cas_seq = inv.sequence[cassette[0]:cassette[1]]
        h = _best_hit_to_target(target.sequence, cas_seq, "cassette")
        if not h or (h.match_len / max(1, len(cas_seq))) < 0.7:
            continue

        # Slot length proxy: difference between target-cassette length and inventory-cassette length,
        # plus the 17-250 bp window we expect.
        slot_len_estimate = abs(len(cas_seq) - h.match_len) + 20
        if slot_len_estimate < 1 or slot_len_estimate > 300:
            continue

        guide_kind = "sgrna" if slot_len_estimate <= 30 else ("trna_array" if slot_len_estimate <= 80 else "pegrna")
        success = 0.99 if guide_kind == "sgrna" else 0.93
        work = 1 if guide_kind == "sgrna" else 2
        score = success * (1.0 / (1.0 + work))

        return FeasibilityReport(
            workflow="sgrna_golden_gate",
            feasible=True,
            score=score,
            confidence=0.8,
            work_estimate=work,
            success_estimate=success,
            rationale=(
                f"{inv.name} carries pol3 cassette {cassette[0]}-{cassette[1]} with {chosen_enz} cuts-twice-only; "
                f"slot~{slot_len_estimate} bp → {guide_kind}"
            ),
            handler_args={
                "intent": "sgrna_golden_gate",
                "vector_source": inv.name,
                "enzyme": chosen_enz,
                "guide_kind": guide_kind,
                "estimated_slot_bp": slot_len_estimate,
            },
            validation_mode="strict",
        )

    return FeasibilityReport(
        workflow="sgrna_golden_gate", feasible=False,
        rationale="No inventory plasmid carries a pol3 cassette with a cuts-twice-only Type IIs enzyme.",
        validation_mode="strict",
    )


def _loc_to_range(loc: str) -> Optional[Tuple[int, int]]:
    """Parse a GenBank location string like '100..200' or 'complement(100..200)' into (start, end)."""
    import re
    m = re.search(r"(\d+)\.\.(\d+)", loc or "")
    if not m:
        return None
    return int(m.group(1)) - 1, int(m.group(2))




# ---------------------------------------------------------------------------
# MultiSite Gateway (Hellsgate / pDONR P1-Pn family)
# ---------------------------------------------------------------------------

# Compatible MultiSite chains: ordered list of att-subtype labels expected on target,
# mapped to (reaction_type, donor signatures, destination signature)
MULTISITE_CHAINS = [
    {
        "name": "3-fragment LR Plus (P1-P2-P3-P4)",
        "target_subtypes": ["attB1", "attB2", "attB3", "attB4"],
        "alt_target_subtypes": ["attL1", "attL2", "attL3", "attL4"],
        "donor_pairs": [("attP1", "attP2"), ("attP2r", "attP3"), ("attP3", "attP4")],
        "destination_pair": ("attR4", "attR1"),
    },
    {
        "name": "3-fragment Hellsgate (P4-P1r + P1-P2 + P2r-P3)",
        "target_subtypes": ["attB4", "attB1", "attB2", "attB3"],
        "alt_target_subtypes": ["attL4", "attL1", "attL2", "attL3"],
        "donor_pairs": [("attP4", "attP1r"), ("attP1", "attP2"), ("attP2r", "attP3")],
        "destination_pair": ("attR4", "attR3"),
    },
    {
        "name": "2-fragment LR Plus (P1-P5 + P5r-P2)",
        "target_subtypes": ["attB1", "attB5", "attB2"],
        "alt_target_subtypes": ["attL1", "attL5", "attL2"],
        "donor_pairs": [("attP1", "attP5"), ("attP5r", "attP2")],
        "destination_pair": ("attR1", "attR2"),
    },
]


def _inventory_att_signatures(inv: InventoryContext) -> set:
    if not inv.cloning:
        return set()
    return {f.subtype for f in inv.cloning.features if f.feature_family == "gateway_att"}


def _assess_multisite_gateway(target: TargetContext, inventory: List[InventoryContext]) -> FeasibilityReport:
    if not target.cloning:
        return FeasibilityReport(workflow="multisite_gateway", feasible=False,
                                 rationale="No cloning_features payload on target.",
                                 validation_mode="strict")
    target_atts = [f for f in target.cloning.features if f.feature_family == "gateway_att"]
    if len(target_atts) < 3:
        return FeasibilityReport(workflow="multisite_gateway", feasible=False,
                                 rationale="MultiSite Gateway requires ≥3 att cores in target.",
                                 validation_mode="strict")
    target_subs = [f.subtype for f in sorted(target_atts, key=lambda x: x.start)]

    # Find a chain that matches target's att-subtype sequence (B-class or L-class)
    chosen_chain = None
    for chain in MULTISITE_CHAINS:
        for variant in (chain["target_subtypes"], chain["alt_target_subtypes"]):
            if all(s in target_subs for s in variant):
                chosen_chain = chain
                break
        if chosen_chain:
            break
    if not chosen_chain:
        return FeasibilityReport(
            workflow="multisite_gateway", feasible=False,
            rationale=f"Target att chain {target_subs} doesn't match any known MultiSite topology.",
            validation_mode="strict",
        )

    # Match each donor pair to an inventory plasmid
    inv_sigs = [(inv, _inventory_att_signatures(inv)) for inv in inventory]
    donor_matches = []
    for needed_pair in chosen_chain["donor_pairs"]:
        match = next((inv for inv, sig in inv_sigs if all(s in sig for s in needed_pair)), None)
        donor_matches.append((needed_pair, match))

    # Match destination
    dest_match = next((inv for inv, sig in inv_sigs
                       if all(s in sig for s in chosen_chain["destination_pair"])), None)

    n_entries = len(chosen_chain["donor_pairs"])
    matched_donors = sum(1 for _, m in donor_matches if m)
    feasible = (matched_donors == n_entries) and (dest_match is not None)
    success = ((matched_donors / n_entries) ** 1) * (1.0 if dest_match else 0.0) * 0.85
    work = n_entries + 1  # one BP per missing entry + one LR Plus
    score = success * (1.0 / (1.0 + work))

    return FeasibilityReport(
        workflow="multisite_gateway",
        feasible=feasible,
        score=score,
        confidence=0.7,
        work_estimate=work,
        success_estimate=success,
        rationale=(
            f"chain={chosen_chain['name']}; matched {matched_donors}/{n_entries} donors; "
            f"destination={dest_match.name if dest_match else 'MISSING'}"
        ),
        handler_args={
            "intent": "gateway_cloning",
            "reaction_type": "LR_plus",
            "chain": chosen_chain["name"],
            "donors": [m.name if m else None for _, m in donor_matches],
            "destination": dest_match.name if dest_match else None,
            "entries_needed": n_entries,
            "primers_needed": any(m is None for _, m in donor_matches),
        },
        validation_mode="strict",
    )



# ---------------------------------------------------------------------------
# PCR-extension Gibson hybrid (covers short coverage gaps via primer tails)
# ---------------------------------------------------------------------------

def assess_pcr_extension_gibson_feasibility(target: TargetContext, inventory: List[InventoryContext]) -> FeasibilityReport:
    """Hybrid path -- inv_gib's set-cover, but every short coverage gap (< 60 bp)
    is added as a primer-tail extension on the adjacent fragment instead of
    becoming a gBlock. Cleaner than synthesis_fallback when only a handful of
    junctions are missing homology. Wired into route() 2026-04-27."""
    L = len(target.sequence)
    if L == 0:
        return FeasibilityReport(
            workflow="pcr_extension_gibson", feasible=False,
            rationale="Empty target sequence.",
            validation_mode="loose",
        )
    covered = bytearray(L)
    fragments: List[Dict[str, Any]] = []
    for inv in inventory:
        h = _best_hit_to_target(target.sequence, inv.sequence, inv.name)
        if not h:
            continue
        s = max(0, h.tgt_start)
        e = min(L, h.tgt_start + h.match_len)
        if e <= s:
            continue
        for i in range(s, e):
            covered[i] = 1
        fragments.append({"source": inv.name, "start": s, "end": e, "length": e - s})

    if not fragments:
        return FeasibilityReport(
            workflow="pcr_extension_gibson", feasible=False,
            rationale="No inventory coverage of target -- nothing to PCR-amplify; use synthesis_fallback.",
            validation_mode="loose",
        )

    SHORT_GAP_MAX = 60
    in_gap = False
    gap_start = 0
    short_gaps: List[Tuple[int, int]] = []
    long_gaps: List[Tuple[int, int]] = []
    for i in range(L):
        if not covered[i] and not in_gap:
            in_gap = True
            gap_start = i
        elif covered[i] and in_gap:
            in_gap = False
            length = i - gap_start
            (short_gaps if length <= SHORT_GAP_MAX else long_gaps).append((gap_start, i))
    if in_gap:
        length = L - gap_start
        (short_gaps if length <= SHORT_GAP_MAX else long_gaps).append((gap_start, L))

    feasible = (not long_gaps) and bool(fragments)
    n_frags = len(fragments)
    n_extensions = len(short_gaps)
    coverage_frac = sum(covered) / L

    success = (0.95 ** n_frags) * coverage_frac
    work = (n_frags * 2) + n_extensions
    score = success * (1.0 / (1.0 + work))

    return FeasibilityReport(
        workflow="pcr_extension_gibson",
        feasible=feasible,
        score=score,
        confidence=0.7,
        work_estimate=work,
        success_estimate=success,
        rationale=(
            f"{n_frags} inventory fragment(s) cover {coverage_frac:.2f} of target; "
            f"{n_extensions} short gap(s) <= {SHORT_GAP_MAX} bp absorbed as primer extensions; "
            f"{len(long_gaps)} long gap(s) reject (would need synthesis_fallback)."
        ),
        handler_args={
            "intent": "gibson_design",
            "fragments": fragments,
            "primer_tail_extensions": [
                {"start": s, "end": e, "length": e - s} for s, e in short_gaps
            ],
        },
        validation_mode="loose",
    )


# ---------------------------------------------------------------------------
# Synthesis fallback (always feasible)
# ---------------------------------------------------------------------------

GBLOCK_MAX_BP = 1500
GBLOCK_OVERLAP_BP = 25
COST_GBLOCK_PER_BP = 0.10
COST_LONG_PER_BP = 0.30


def assess_synthesis_fallback_feasibility(target: TargetContext, inventory: List[InventoryContext]) -> FeasibilityReport:
    """Always feasible. Uncovered regions become gBlocks; covered regions become PCR amplicons; assemble via Gibson."""
    L = len(target.sequence)
    covered = bytearray(L)
    anchors: List[Dict[str, Any]] = []
    for inv in inventory:
        h = _best_hit_to_target(target.sequence, inv.sequence, inv.name)
        if not h:
            continue
        new_cov = sum(1 for i in range(h.tgt_start, min(h.tgt_end, L)) if not covered[i])
        if new_cov >= 50:
            for i in range(h.tgt_start, min(h.tgt_end, L)):
                covered[i] = 1
            anchors.append({
                "source": inv.name,
                "target_start": h.tgt_start,
                "target_end": min(h.tgt_end, L),
                "match_bp": h.match_len,
            })

    # Compute uncovered runs as synth blocks (split long ones at GBLOCK_MAX_BP)
    synth_blocks: List[Dict[str, Any]] = []
    run_start = None
    for i in range(L):
        if not covered[i]:
            if run_start is None:
                run_start = i
        else:
            if run_start is not None:
                _split_into_gblocks(synth_blocks, target.sequence, run_start, i)
                run_start = None
    if run_start is not None:
        _split_into_gblocks(synth_blocks, target.sequence, run_start, L)

    total_synth_bp = sum(b["length_bp"] for b in synth_blocks)
    n_blocks = len(synth_blocks)
    n_anchors = len(anchors)
    coverage_frac = sum(covered) / max(1, L)

    # Cost: gBlocks ≤500 bp at $0.10/bp, longer fragments at $0.30/bp
    est_cost = sum(
        b["length_bp"] * (COST_GBLOCK_PER_BP if b["length_bp"] <= 500 else COST_LONG_PER_BP)
        for b in synth_blocks
    )

    success = (0.90 ** n_blocks) * max(0.3, coverage_frac)  # never 0; fully synth still works
    work = n_blocks + n_anchors + 1
    score = success * (1.0 / (1.0 + work))

    return FeasibilityReport(
        workflow="synthesis_fallback",
        feasible=True,  # always
        score=score,
        confidence=0.85,
        work_estimate=work,
        success_estimate=success,
        rationale=(
            f"{n_anchors} inventory anchor(s) cover {coverage_frac*100:.1f}%; "
            f"{n_blocks} synth blocks ({total_synth_bp} bp, ~${est_cost:.0f}); Gibson assembly"
        ),
        handler_args={
            "intent": "plasmid_design",
            "strategy": "synth_plus_gibson",
            "synth_blocks": synth_blocks,
            "inventory_anchors": anchors,
            "total_synth_bp": total_synth_bp,
            "est_cost_usd": round(est_cost, 2),
        },
        validation_mode="loose",
    )


def _split_into_gblocks(out: List[Dict[str, Any]], seq: str, start: int, end: int) -> None:
    span_len = end - start
    if span_len <= GBLOCK_MAX_BP:
        out.append({
            "name": f"gblock_{start}_{end}",
            "target_start": start,
            "target_end": end,
            "length_bp": span_len,
            "sequence": seq[start:end],
        })
        return
    # Split into ≤GBLOCK_MAX_BP chunks with overlap
    pos = start
    idx = 0
    while pos < end:
        chunk_end = min(pos + GBLOCK_MAX_BP, end)
        out.append({
            "name": f"gblock_{start}_{end}_part{idx}",
            "target_start": pos,
            "target_end": chunk_end,
            "length_bp": chunk_end - pos,
            "sequence": seq[pos:chunk_end],
        })
        idx += 1
        if chunk_end >= end:
            break
        pos = chunk_end - GBLOCK_OVERLAP_BP


# ---------------------------------------------------------------------------
# Validation — compare an in-silico assembled product against the target
# ---------------------------------------------------------------------------

def _rotate_canonical(seq: str) -> str:
    """Rotate a circular sequence to start at its lexicographically smallest
    rotation. Deterministic canonical form — if two circular sequences encode
    the same circle, their canonical rotations are identical."""
    if not seq:
        return seq
    doubled = seq + seq
    n = len(seq)
    # Pick the smallest window of length n across all rotations.
    best = 0
    for i in range(1, n):
        if doubled[i:i + n] < doubled[best:best + n]:
            best = i
    return doubled[best:best + n]


def _strict_match(product: str, target: str) -> Dict[str, Any]:
    """Strict validation: canonical-rotate both (and try reverse complement),
    require exact sequence equality."""
    pu = (product or "").upper()
    tu = (target or "").upper()
    if len(pu) != len(tu):
        return {
            "mode": "strict", "passed": False,
            "reason": f"length mismatch: product={len(pu)} bp, target={len(tu)} bp",
        }
    p_canon = _rotate_canonical(pu)
    t_canon = _rotate_canonical(tu)
    if p_canon == t_canon:
        return {"mode": "strict", "passed": True, "reason": "sequence-identical (after circular rotation)"}
    # Try reverse-complement of product (in case assembly emitted the other strand)
    rc = _revcomp(pu)
    if _rotate_canonical(rc) == t_canon:
        return {"mode": "strict", "passed": True, "reason": "sequence-identical (reverse-complement match)"}
    # Report a simple mismatch count using aligned canonical forms.
    mismatches = sum(1 for a, b in zip(p_canon, t_canon) if a != b)
    return {
        "mode": "strict", "passed": False,
        "reason": f"{mismatches}/{len(t_canon)} positions mismatch after canonical rotation",
    }


def _module_type_set(modules: List[Dict[str, Any]]) -> set:
    return {m.get("module_type", "") for m in modules if m.get("module_type")}


def _feature_name_set(annotations: List[Dict[str, Any]], min_pi_permatch: float = 90.0) -> set:
    """Names of annotations whose pi_permatch ≥ threshold. Falls back to all names
    if pi_permatch isn't present on the payload (e.g. module-only annotations)."""
    out = set()
    for a in annotations:
        pi = a.get("pi_permatch") if isinstance(a, dict) else getattr(a, "pi_permatch", None)
        name = (a.get("name") if isinstance(a, dict) else getattr(a, "name", None)) or ""
        if not name:
            continue
        if pi is None or (isinstance(pi, (int, float)) and pi >= min_pi_permatch):
            out.add(name)
    return out


def _interaction_set(interactions: List[Dict[str, Any]]) -> set:
    """Canonical interaction signature: (rule_id, interaction_type)."""
    out = set()
    for ix in interactions or []:
        rid = ix.get("rule_id") or ""
        itype = ix.get("interaction_type") or ""
        out.add((rid, itype))
    return out



# ---------------------------------------------------------------------------
# CDS functional-context check (loose-match augmentation, v5)
#
# Verification rule:
#   An inserted CDS is "expressed" in the assembled product iff
#     1. it sits inside a module whose module_type is one of the rule-schema
#        expression-cassette types (those that wrap a CDS module), AND
#     2. the interaction-builder emits at least one expression-driving
#        interaction in which the CDS appears as a `template` participant
#        with a stimulator (upstream regulatory) and/or modifier (downstream
#        regulatory) participant.
# This avoids re-running Step 1 of the annotation pipeline; we read the
# already-emitted modules and the already-emitted interaction graph.
# ---------------------------------------------------------------------------

import re as _re

EXPRESSION_CASSETTE_TYPES = {
    "mammalian_pol2_expression_cassette",
    "mammalian_lentiviral_expression_cassette",
    "lentiviral_payload",
    "bacterial_expression_cassette",
    "bacterial_marker_cassette",
    "pol2_expression_cassette",
    "expression_cassette",
}

_REGULATORY_ROLES = {"stimulator", "modifier"}
_TEMPLATE_ROLES = {"template"}

_CDS_FEATURE_TYPES = {"CDS", "cds", "gene"}
_CDS_KB_CLASSES = {"coding_sequence", "cds", "protein_coding"}
_LOC_RE = _re.compile(r"(complement\()?\s*(\d+)\s*\.\.\s*(\d+)")


def _is_cds_class(ann: Dict[str, Any]) -> bool:
    """True if a feature is CDS-class by SO/SBOL/KB type, class, or subclass."""
    ftype = (ann.get("type") or "").strip()
    if ftype in _CDS_FEATURE_TYPES:
        return True
    kb_class = (ann.get("kb_class") or ann.get("class") or "").strip().lower()
    kb_subclass = (ann.get("kb_subclass") or ann.get("subclass") or "").strip().lower()
    if kb_class in _CDS_KB_CLASSES or kb_subclass in _CDS_KB_CLASSES:
        return True
    so_role = (ann.get("so_role") or "").lower()
    if "cds" in so_role or "coding" in so_role:
        return True
    return False


_PROMOTER_FEATURE_TYPES = {"promoter", "Promoter"}
_ENHANCER_FEATURE_TYPES = {"enhancer", "Enhancer"}
_POLYA_FEATURE_TYPES = {"polyA_signal", "polya_signal", "polya", "polyA"}
_TERMINATOR_FEATURE_TYPES = {"terminator", "Terminator"}
_RBS_FEATURE_TYPES = {"RBS", "rbs", "ribosome_binding_site", "Shine-Dalgarno"}

_PROMOTER_KB_CLASSES = {"promoter", "pol2_promoter", "pol3_promoter"}
_ENHANCER_KB_CLASSES = {"enhancer"}
_POLYA_KB_CLASSES = {"polya", "poly_a", "polya_signal"}
_TERMINATOR_KB_CLASSES = {"terminator"}
_RBS_KB_CLASSES = {"rbs", "ribosome_binding_site", "shine_dalgarno"}

# Each expression-feature kind maps to the participant role(s) the feature
# is expected to play in an expression interaction.
_INTRON_FEATURE_TYPES = {"intron", "Intron"}
_INTRON_KB_CLASSES = {"intron"}

# Distance windows for orientation-only support features.
_SUPPORT_FEATURE_WINDOW = 5000


def _is_intron_class(ann: Dict[str, Any]) -> bool:
    if _matches_kind(ann, _INTRON_FEATURE_TYPES, _INTRON_KB_CLASSES):
        return True
    name = (ann.get("name") or "").lower()
    return name.startswith("intron") or " intron" in name


_RECOMBINATION_NAME_PATTERNS = (
    "loxp", "lox66", "lox71", "lox511", "lox2272", "loxn",
    "frt", "f3 site",
    "attl", "attr", "attb", "attp",
)
_RECOMBINATION_KB_CLASSES = {
    "lox_site", "loxp_site", "frt_site",
    "att_site", "attl_site", "attr_site", "attb_site", "attp_site",
    "recombination_site",
}

_POLYPROTEIN_2A_PATTERNS = ("p2a", "t2a", "e2a", "f2a", "2a peptide", "2a self-cleaving")
_POLYPROTEIN_2A_KB_CLASSES = {"p2a", "t2a", "e2a", "f2a", "2a_peptide", "ribosomal_skip"}

# Common fusion-tag identifiers. Order matters: more-specific first so that
# e.g. "Strep-tag II" is caught before falling through to "tag".
_TAG_NAME_PATTERNS = (
    ("His-tag",  ("his-tag", "6xhis", "his6", "8xhis", "his8", "histag", "10xhis")),
    ("FLAG",     ("3xflag", "1xflag", "flag-tag", "flag epitope")),
    ("HA",       ("ha tag", "ha-tag", "ha epitope", "3xha", "ypydvpdya")),
    ("V5",       ("v5 tag", "v5-tag", "v5 epitope")),
    ("Myc",      ("myc tag", "myc-tag", "myc epitope", "eqklisee")),
    ("Strep",    ("strep-tag", "strep tag", "twin-strep")),
    ("T7",       ("t7 tag", "t7-tag", "t7 epitope")),
    ("AviTag",   ("avitag", "avi tag", "biotin acceptor")),
    ("GST",      ("gst tag", "gst-tag", "glutathione s-transferase tag")),
    ("MBP",      ("mbp tag", "mbp-tag", "maltose binding protein tag")),
    ("SUMO",     ("sumo tag", "sumo-tag", "smt3")),
    ("NLS",      ("nls", "nuclear localization", "sv40 nls", "bipartite nls")),
    ("NES",      ("nes", "nuclear export", "leptomycin")),
)
_TAG_KB_CLASSES = {"tag", "fusion_tag", "epitope_tag", "localization_signal", "nls", "nes"}


def _is_recombination_site(ann: Dict[str, Any]) -> bool:
    name = (ann.get("name") or "").lower()
    if any(p in name for p in _RECOMBINATION_NAME_PATTERNS):
        return True
    kb_class = (ann.get("kb_class") or ann.get("class") or "").strip().lower()
    kb_subclass = (ann.get("kb_subclass") or ann.get("subclass") or "").strip().lower()
    if kb_class in _RECOMBINATION_KB_CLASSES or kb_subclass in _RECOMBINATION_KB_CLASSES:
        return True
    return False


def _is_polyprotein_2a(ann: Dict[str, Any]) -> bool:
    name = (ann.get("name") or "").lower()
    if any(p in name for p in _POLYPROTEIN_2A_PATTERNS):
        return True
    kb_class = (ann.get("kb_class") or ann.get("class") or "").strip().lower()
    kb_subclass = (ann.get("kb_subclass") or ann.get("subclass") or "").strip().lower()
    if kb_class in _POLYPROTEIN_2A_KB_CLASSES or kb_subclass in _POLYPROTEIN_2A_KB_CLASSES:
        return True
    return False


def _classify_tag(ann: Dict[str, Any]) -> Optional[str]:
    """Return the canonical tag label (e.g. 'His-tag', 'FLAG'), or None."""
    name = (ann.get("name") or "").lower()
    for canonical, patterns in _TAG_NAME_PATTERNS:
        for p in patterns:
            if p in name:
                return canonical
    kb_class = (ann.get("kb_class") or ann.get("class") or "").strip().lower()
    if kb_class in _TAG_KB_CLASSES:
        return name or "tag"
    return None


# Module-level rule validation: each module_type listed here must (or may)
# participate in the named interaction rule_ids in the named role. A required
# entry that doesn't fire emits a `missing_module_interaction` warning.
MODULE_VALIDATION_RULES: Dict[str, List[Dict[str, Any]]] = {
    "lentiviral_payload": [
        {"rule_id": "INT-LENTI-UR-01", "role": "template", "required": True,
         "what": "Lentiviral upstream regulatory drives payload transcription"},
        {"rule_id": "INT-LENTI-DR-01", "role": "template", "required": True,
         "what": "Lentiviral downstream regulatory modifies payload"},
    ],
    "lentiviral_upstream_regulatory": [
        {"rule_id": "INT-LENTI-UR-01", "role": "stimulator", "required": True,
         "what": "Upstream regulatory must stimulate the payload"},
    ],
    "lentiviral_downstream_regulatory": [
        {"rule_id": "INT-LENTI-DR-01", "role": "modifier", "required": True,
         "what": "Downstream regulatory must modify the payload"},
    ],
    "insulated_expression_block": [
        {"rule_id": "INT-INS-01", "role": "modifier", "required": True,
         "what": "Insulator pair brackets the expression boundary"},
    ],
    "lac_alpha_blue_white_module": [
        {"rule_id": "INT-LACBW-01", "role": None, "required": True,
         "what": "lac promoter drives lacZα transcription"},
        # Either INT-LACBW-04 (intact MCS) or INT-LACBW-04b (post-cloning insert
        # disrupts lacZα) satisfies the disruption contract — see group_required.
        {"rule_id": "INT-LACBW-04", "role": None, "required": False,
         "what": "MCS insertion disrupts lacZα reading frame (pre-cloning state)"},
        {"rule_id": "INT-LACBW-04b", "role": None, "required": False,
         "what": "Cloned insert disrupts lacZα reading frame (post-cloning state)"},
        {"group_required": ["INT-LACBW-04", "INT-LACBW-04b"],
         "what": "Either the empty MCS or an inserted CDS must disrupt lacZα"},
        {"rule_id": "INT-LACBW-02", "role": None, "required": False,
         "what": "lac operator inhibits lac promoter (optional)"},
        {"rule_id": "INT-LACBW-03", "role": None, "required": False,
         "what": "LacI binds lac operator (optional, requires lacI gene)"},
    ],
    "mobilizable_replicon": [
        {"rule_id": "INT-MOB-01", "role": "reactant", "required": True,
         "what": "bom oriT must be paired with a replication origin"},
        {"rule_id": "INT-MOB-02", "role": "stimulator", "required": False,
         "what": "Relaxase (mob/traI) acts on bom — optional; trans-encoded in pUC"},
    ],
    "tn3_transposon": [
        {"rule_id": "INT-TRANSP-TN3-01", "role": "reactant", "required": True,
         "what": "Tn3 inverted-repeat pair must be present (both ends)"},
        {"rule_id": "INT-TRANSP-TN3-02", "role": "template", "required": False,
         "what": "Cargo CDS (bla / tnpA) enclosed by the IR pair"},
    ],
    "lac_promoter_regulatory_unit": [
        {"rule_id": "INT-LAC-OP-01", "role": "modifier", "required": True,
         "what": "lac operator must repress the lac promoter it gates"},
        {"rule_id": "INT-CAP-01", "role": "modifier", "required": False,
         "what": "CAP site provides catabolite activation when present"},
    ],
    "tet_inducible_expression_cassette": [
        {"rule_id": "INT-TET-IND-02", "role": "template", "required": True,
         "what": "Tet-responsive promoter drives a payload CDS"},
        {"rule_id": "INT-TET-IND-03", "role": "modifier", "required": True,
         "what": "tetO modulates the responsive promoter"},
        {"rule_id": "INT-TET-IND-01", "role": "template", "required": False,
         "what": "Transactivator (tTA/rtTA) binds tetO — required only if tet_regulator_cassette also present"},
    ],
    "yeast_selection_cassette": [
        {"rule_id": "INT-YEAST-SEL-01", "role": "template", "required": True,
         "what": "TEF/heterologous promoter drives the MX-cassette resistance gene"},
    ],
    "baculovirus_recombination_cassette": [
        {"rule_id": "INT-BAC-RECOMB-01", "role": "reactant", "required": True,
         "what": "Both homology arms (ORF1629 + lef2/ORF603) must be present"},
        {"rule_id": "INT-BAC-RECOMB-02", "role": "template", "required": False,
         "what": "polh/p10 expression block enclosed by homology arms"},
    ],
    # 2026-05-05: orientation-aware Gateway substrates
    "gateway_excision_module": [
        {"rule_id": "INT-REC-GW-EXC-01", "role": "reactant", "required": True,
         "what": "Same-strand compatible att pair → BP/LR clonase intramolecular deletion"},
    ],
    "gateway_inversion_module": [
        {"rule_id": "INT-REC-GW-INV-01", "role": "reactant", "required": True,
         "what": "Opposite-strand outward-pointing att pair → BP/LR clonase intramolecular inversion"},
    ],
    "gateway_intermolecular_module": [
        {"rule_id": "INT-REC-GW-INTER-01", "role": "reactant", "required": True,
         "what": "Opposite-strand inward-pointing compatible att pair → intermolecular BP/LR recombination substrate"},
    ],
    # 2026-05-06: phage-RNAP expression cassette (T7 / T3 / SP6)
    "phage_rnap_expression_cassette": [
        {"rule_id": "INT-T7-EXPR-01", "role": "template", "required": False,
         "what": "T7 phage promoter drives same-strand CDS within 200 bp downstream"},
        {"rule_id": "INT-T3-EXPR-01", "role": "template", "required": False,
         "what": "T3 phage promoter drives same-strand CDS within 200 bp downstream"},
        {"rule_id": "INT-SP6-EXPR-01", "role": "template", "required": False,
         "what": "SP6 phage promoter drives same-strand CDS within 200 bp downstream"},
        {"group_required": ["INT-T7-EXPR-01", "INT-T3-EXPR-01", "INT-SP6-EXPR-01"],
         "what": "At least one phage-RNAP transcription interaction must fire for the cassette"},
    ],
    # 2026-05-12: phage_rnap_orientation_mismatch removed. The general
    # orientation-mismatch detector in `_check_cds_functional_context`
    # (via `_detect_orientation_mismatch`) covers the same case for ANY
    # EXPRESSION_FEATURE_KINDS pair, not just T7+CDS.
}


def _module_overlaps_participant(mod: Dict[str, Any], p: Dict[str, Any]) -> bool:
    """True if the participant interval overlaps the module by >= 50% of the
    module span (or by >= 50% of the participant span when the participant
    is smaller — handles point-anchored participants)."""
    ms = int(mod.get("start", 0) or 0)
    me = int(mod.get("end", 0) or 0)
    ps = p.get("start"); pe = p.get("end")
    if ps is None or pe is None:
        return False
    try:
        ps, pe = int(ps), int(pe)
    except (TypeError, ValueError):
        return False
    if pe <= ms or me <= ps:
        return False
    overlap = min(pe, me) - max(ps, ms)
    return overlap / max(1, min(me - ms, pe - ps)) >= 0.5


def _check_module_validation(
    modules: List[Dict[str, Any]],
    interactions: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Walk every module against MODULE_VALIDATION_RULES and emit one
    warning per missing required (rule_id, role) pair."""
    warnings: List[Dict[str, Any]] = []
    for m in modules or []:
        mt = m.get("module_type") or ""
        rules = MODULE_VALIDATION_RULES.get(mt)
        if not rules:
            continue
        ms = int(m.get("start", 0) or 0)
        me = int(m.get("end", 0) or 0)
        mname = m.get("name") or mt

        # Track which rule_ids actually fired for this module — used by
        # the group_required predicate evaluated after the per-rule loop.
        fired_rule_ids: set = set()

        for rule in rules:
            # group_required: at least one of the listed rule_ids must have
            # fired. Evaluated after the main loop so we can read fired_rule_ids.
            if "group_required" in rule:
                continue
            rule_id = rule["rule_id"]
            expected_role = rule.get("role")
            required = bool(rule.get("required"))

            matched = []
            for ix in interactions or []:
                if ix.get("rule_id") != rule_id:
                    continue
                if expected_role is None:
                    matched.append(ix)
                    continue
                for p in (ix.get("participants") or []):
                    if (p.get("role") or "").lower() != expected_role.lower():
                        continue
                    if _module_overlaps_participant(m, p):
                        matched.append(ix)
                        break

            if matched:
                fired_rule_ids.add(rule_id)
            if matched or not required:
                continue

            warnings.append({
                "module_type": mt,
                "module_name": mname,
                "module_span": [ms, me],
                "rule_id": rule_id,
                "expected_role": expected_role,
                "issue": "missing_module_interaction",
                "gaps": [
                    f"required interaction {rule_id}" +
                    (f" with this module as {expected_role}" if expected_role else "") +
                    " did not fire"
                ],
                "detail": rule.get("what") or "",
                "remediation": (
                    "Verify the rule-schema cassette assembler emits this "
                    "interaction. The module is detected but the corresponding "
                    "interaction-builder entry didn't produce the expected "
                    "participant binding."
                ),
            })

        # group_required check: at least one of the listed rule_ids must have
        # fired against this module (any role). Used by lac_alpha_blue_white
        # to accept either INT-LACBW-04 (intact MCS) or INT-LACBW-04b (post-cloning).
        for rule in rules:
            group = rule.get("group_required")
            if not group:
                continue
            satisfied = False
            for rid in group:
                for ix in interactions or []:
                    if ix.get("rule_id") != rid:
                        continue
                    for p in (ix.get("participants") or []):
                        if _module_overlaps_participant(m, p):
                            satisfied = True
                            break
                    if satisfied:
                        break
                if satisfied:
                    break
            if satisfied:
                continue
            warnings.append({
                "module_type": mt,
                "module_name": mname,
                "module_span": [ms, me],
                "rule_id": "+".join(group),
                "expected_role": None,
                "issue": "missing_group_interaction",
                "gaps": [f"none of {group} fired against this module"],
                "detail": rule.get("what") or "",
                "remediation": (
                    "At least one of the listed interaction rule_ids must fire. "
                    "Check the interaction-builder for this module type."
                ),
            })

    return warnings


EXPRESSION_FEATURE_KINDS = (
    ("cds",                ("template",),),
    ("promoter",           ("stimulator",),),
    ("polya",              ("modifier",),),
    ("terminator",         ("modifier",),),
    ("rbs",                ("stimulator", "modifier"),),
    # 2026-05-03: promoted from () — enhancer/intron now require participant edges.
    ("enhancer",           ("stimulator",),),
    ("intron",             ("modifier",),),
    # v10 additions:
    ("recombination_site", ("reactant",),),
    ("polyprotein_2a",     ("stimulator",),),
    # 2026-05-03 audit-driven additions:
    ("tet_operator",       ("modifier",),),
    ("inducible_promoter", ("stimulator",),),
    ("lac_operator",       ("inhibitor",),),
    ("cap_binding_site",   ("stimulator",),),
    ("rep_origin",         ("modifier",),),
    ("oriT",               ("reactant",),),
    ("inverted_repeat",    ("reactant",),),
    ("ncrna",              ("reactant",),),
    ("rre",                ("reactant",),),
    ("misc_recomb",        ("reactant",),),
    ("homology_region",    ("reactant",),),
    # tag is a containment check, not a participant check.
    ("tag",                (),),
)


def _matches_kind(ann: Dict[str, Any], type_set: set, kb_class_set: set) -> bool:
    ftype = (ann.get("type") or "").strip()
    if ftype in type_set:
        return True
    name = (ann.get("name") or "").lower()
    for ts in type_set:
        if ts.lower() in name and len(ts) >= 3:
            # Avoid e.g. "rbs" matching "rbs1" mid-word for promoter names.
            pass
    kb_class = (ann.get("kb_class") or ann.get("class") or "").strip().lower()
    kb_subclass = (ann.get("kb_subclass") or ann.get("subclass") or "").strip().lower()
    if kb_class in kb_class_set or kb_subclass in kb_class_set:
        return True
    so_role = (ann.get("so_role") or "").lower()
    for ts in kb_class_set:
        if ts and ts in so_role:
            return True
    return False


def _is_promoter_class(ann: Dict[str, Any]) -> bool:
    if _matches_kind(ann, _PROMOTER_FEATURE_TYPES, _PROMOTER_KB_CLASSES):
        return True
    name = (ann.get("name") or "").lower()
    return "promoter" in name


def _is_enhancer_class(ann: Dict[str, Any]) -> bool:
    if _matches_kind(ann, _ENHANCER_FEATURE_TYPES, _ENHANCER_KB_CLASSES):
        return True
    name = (ann.get("name") or "").lower()
    return "enhancer" in name


def _is_polya_class(ann: Dict[str, Any]) -> bool:
    if _matches_kind(ann, _POLYA_FEATURE_TYPES, _POLYA_KB_CLASSES):
        return True
    name = (ann.get("name") or "").lower()
    return "polya" in name or "poly(a)" in name or "poly a" in name


def _is_terminator_class(ann: Dict[str, Any]) -> bool:
    if _matches_kind(ann, _TERMINATOR_FEATURE_TYPES, _TERMINATOR_KB_CLASSES):
        return True
    name = (ann.get("name") or "").lower()
    return "terminator" in name


def _is_rbs_class(ann: Dict[str, Any]) -> bool:
    if _matches_kind(ann, _RBS_FEATURE_TYPES, _RBS_KB_CLASSES):
        return True
    name = (ann.get("name") or "").lower()
    return ("rbs" in name and "rbs1" not in name) or "shine-dalgarno" in name or "shine_dalgarno" in name


def _classify_expression_feature(ann: Dict[str, Any]) -> Optional[str]:
    """Return the first matching expression-feature kind, or None."""
    # CDS first (so that 2A submodules whose annotator type=CDS still
    # get the polyprotein_2a classification when their NAME matches).
    if _is_polyprotein_2a(ann):
        return "polyprotein_2a"
    if _is_recombination_site(ann):
        return "recombination_site"
    if _classify_tag(ann) is not None:
        return "tag"
    if _is_cds_class(ann):
        return "cds"
    if _is_promoter_class(ann):
        return "promoter"
    if _is_enhancer_class(ann):
        return "enhancer"
    if _is_polya_class(ann):
        return "polya"
    if _is_terminator_class(ann):
        return "terminator"
    if _is_rbs_class(ann):
        return "rbs"
    if _is_intron_class(ann):
        return "intron"
    return None


def _ann_span(ann: Dict[str, Any]) -> Tuple[int, int, int]:
    """(start, end, strand) — parses GenBank `location` strings when start/end
    are not on the dict directly. Accepts either `strand` (GenBank-style)
    or `direction` (pLannotate-style) as the strand source."""
    start = ann.get("start")
    end = ann.get("end")
    strand_raw = ann.get("strand")
    if strand_raw is None:
        strand_raw = ann.get("direction")
    if start is None or end is None:
        loc = ann.get("location") or ""
        m = _LOC_RE.search(loc)
        if m:
            comp = m.group(1)
            start = int(m.group(2)) - 1  # GenBank 1-based -> 0-based
            end = int(m.group(3))
            if comp:
                strand_raw = -1
    try:
        start = int(start or 0)
        end = int(end or 0)
    except (TypeError, ValueError):
        start, end = 0, 0
    if strand_raw is None:
        strand_raw = 1
    try:
        strand = int(strand_raw)
    except (TypeError, ValueError):
        strand = 1 if str(strand_raw) in ("+", "1", "fwd", "forward") else -1
    return start, end, (1 if strand >= 0 else -1)


def _participant_matches_cds(
    participant: Dict[str, Any],
    cds_start: int,
    cds_end: int,
    name_hint: Optional[str] = None,
) -> bool:
    """True if an interaction participant references the same CDS feature.

    We consider a match when either (a) the participant's interval overlaps
    the CDS interval by >= 50%, or (b) the participant's name equals the CDS
    name (case-insensitive). Either signal is enough to wire the participant
    to the CDS being verified.
    """
    pname = (participant.get("name") or "").strip().lower()
    if name_hint and pname and pname == name_hint.strip().lower():
        return True
    ps = participant.get("start")
    pe = participant.get("end")
    if ps is None or pe is None:
        return False
    try:
        ps, pe = int(ps), int(pe)
    except (TypeError, ValueError):
        return False
    if pe <= cds_start or cds_end <= ps:
        return False
    overlap = min(pe, cds_end) - max(ps, cds_start)
    span = max(1, cds_end - cds_start)
    return overlap / span >= 0.5


def _participant_matches_feature(
    participant: Dict[str, Any],
    feat_start: int,
    feat_end: int,
    name_hint: Optional[str] = None,
) -> bool:
    """True iff a participant references the same feature span (>=50% interval
    overlap) or carries the same name (case-insensitive)."""
    pname = (participant.get("name") or "").strip().lower()
    if name_hint and pname and pname == name_hint.strip().lower():
        return True
    ps = participant.get("start")
    pe = participant.get("end")
    if ps is None or pe is None:
        return False
    try:
        ps, pe = int(ps), int(pe)
    except (TypeError, ValueError):
        return False
    if pe <= feat_start or feat_end <= ps:
        return False
    overlap = min(pe, feat_end) - max(ps, feat_start)
    span = max(1, feat_end - feat_start)
    return overlap / span >= 0.5


def _expression_interactions_for_feature(
    ann: Dict[str, Any],
    interactions: List[Dict[str, Any]],
    expected_roles: Tuple[str, ...],
) -> Dict[str, List[Dict[str, Any]]]:
    """Return interactions where this feature participates in any of the
    expected role(s). Result keyed by the role under which it was matched.
    """
    fs, fe, _ = _ann_span(ann)
    name_hint = ann.get("name") or ann.get("feature_name")
    by_role: Dict[str, List[Dict[str, Any]]] = {role: [] for role in expected_roles}

    for ix in interactions or []:
        for p in (ix.get("participants") or []):
            role = (p.get("role") or "").lower()
            if role not in by_role:
                continue
            if _participant_matches_feature(p, fs, fe, name_hint):
                by_role[role].append(ix)
                break  # one match per interaction is enough

    return by_role


def _expression_module_for_cds(
    cds_ann: Dict[str, Any],
    modules: List[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """Return the first expression-cassette module whose interval contains the
    CDS, or None. The module's module_type must be in EXPRESSION_CASSETTE_TYPES.
    """
    cs, ce, _ = _ann_span(cds_ann)
    for m in modules or []:
        if (m.get("module_type") or "") not in EXPRESSION_CASSETTE_TYPES:
            continue
        ms = int(m.get("start", 0) or 0)
        me = int(m.get("end", 0) or 0)
        if ms <= cs and ce <= me:
            return m
    return None


def _expression_interactions_for_cds(
    cds_ann: Dict[str, Any],
    interactions: List[Dict[str, Any]],
) -> Dict[str, List[Dict[str, Any]]]:
    """Return interactions where the CDS is a `template` participant alongside
    a regulatory partner.

    Result has two lists keyed by direction:
      - "upstream":   interactions with a stimulator partner (UR drives CDS)
      - "downstream": interactions with a modifier partner   (DR modifies CDS)
    """
    cs, ce, _ = _ann_span(cds_ann)
    name_hint = cds_ann.get("name") or cds_ann.get("feature_name")
    upstream: List[Dict[str, Any]] = []
    downstream: List[Dict[str, Any]] = []

    for ix in interactions or []:
        participants = ix.get("participants") or []
        cds_in_ix = False
        regulatory_partners: List[Dict[str, Any]] = []
        for p in participants:
            role = (p.get("role") or "").lower()
            if role in _TEMPLATE_ROLES and _participant_matches_cds(p, cs, ce, name_hint):
                cds_in_ix = True
            elif role in _REGULATORY_ROLES:
                regulatory_partners.append(p)
        if not cds_in_ix or not regulatory_partners:
            continue
        for p in regulatory_partners:
            if (p.get("role") or "").lower() == "stimulator":
                upstream.append(ix)
                break
        for p in regulatory_partners:
            if (p.get("role") or "").lower() == "modifier":
                downstream.append(ix)
                break

    return {"upstream": upstream, "downstream": downstream}


def _collect_active_participants(
    interactions: List[Dict[str, Any]],
) -> Dict[str, List[Dict[str, Any]]]:
    """Return active expression participants by role across all interactions.
    Each entry is the participant dict (carries name + start + end + strand).
    """
    out: Dict[str, List[Dict[str, Any]]] = {
        "stimulator": [], "template": [], "modifier": [], "promoter_stim": [],
    }
    for ix in interactions or []:
        for p in (ix.get("participants") or []):
            role = (p.get("role") or "").lower()
            if role in out:
                out[role].append(p)
            if role == "stimulator":
                pname = (p.get("name") or "").lower()
                if "promoter" in pname or "regulatory" in pname:
                    out["promoter_stim"].append(p)
                else:
                    # Default: any stimulator that is at the upstream-regulatory
                    # position is also accepted as a promoter for orientation
                    # comparison.
                    out["promoter_stim"].append(p)
    return out


def _participant_strand(p: Dict[str, Any]) -> int:
    s = p.get("strand", 1)
    try:
        s = int(s)
    except (TypeError, ValueError):
        s = 1
    return 1 if s >= 0 else -1


def _nearest_active_partner(
    feat_start: int,
    feat_end: int,
    candidates: List[Dict[str, Any]],
    window: int,
) -> Optional[Dict[str, Any]]:
    """Return the active participant nearest to the feature whose interval is
    within `window` bp of the feature, or None."""
    best = None
    best_gap = None
    for p in candidates:
        ps = p.get("start"); pe = p.get("end")
        if ps is None or pe is None:
            continue
        try:
            ps, pe = int(ps), int(pe)
        except (TypeError, ValueError):
            continue
        # Distance is 0 if intervals overlap.
        if pe <= feat_start:
            gap = feat_start - pe
        elif ps >= feat_end:
            gap = ps - feat_end
        else:
            gap = 0
        if gap > window:
            continue
        if best_gap is None or gap < best_gap:
            best_gap = gap
            best = p
    return best


def _paired_promoter_for_nearest_cds(
    feat_start: int,
    feat_end: int,
    interactions: List[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """Find the active template CDS nearest the feature, then return the
    stimulator participant from that CDS's expression interaction.

    The stimulator participant is the promoter inside the upstream regulatory
    module of the same cassette — exactly what enhancer / intron orientation
    should be compared against.
    """
    best_cds_ix = None
    best_gap = None
    for ix in interactions or []:
        for p in (ix.get("participants") or []):
            if (p.get("role") or "").lower() != "template":
                continue
            ps = p.get("start"); pe = p.get("end")
            if ps is None or pe is None:
                continue
            try:
                ps, pe = int(ps), int(pe)
            except (TypeError, ValueError):
                continue
            if pe <= feat_start:
                gap = feat_start - pe
            elif ps >= feat_end:
                gap = ps - feat_end
            else:
                gap = 0
            if best_gap is None or gap < best_gap:
                best_gap = gap
                best_cds_ix = ix
    if best_cds_ix is None:
        return None
    for p in (best_cds_ix.get("participants") or []):
        if (p.get("role") or "").lower() == "stimulator":
            return p
    return None


# 2026-05-12: general orientation-mismatch detection (replaces the
# phage-RNAP-specific detector). For each EXPRESSION_FEATURE_KIND feature
# whose expected interaction role isn't satisfied, search for a partner
# feature of the natural-partner kind within proximity. If a partner is
# found on the OPPOSITE strand AND no same-strand partner exists in the
# same window, the feature is mis-oriented — flipping it would form the
# expression cassette.
#
# Partnerships are directional. Each tuple is (partner_kind, side, max_gap)
# where side is "upstream" (lower coords for + strand, higher coords for -
# strand) or "downstream", relative to the feature's transcription
# direction.
_ORIENTATION_PARTNERSHIPS = {
    "cds":        [("promoter",   "upstream",   500),
                   ("polya",      "downstream", 500),
                   ("terminator", "downstream", 500)],
    "promoter":   [("cds",        "downstream", 500)],
    "polya":      [("cds",        "upstream",   500)],
    "terminator": [("cds",        "upstream",   500)],
    "rbs":        [("cds",        "downstream", 60),
                   ("promoter",   "upstream",   200)],
    "enhancer":   [("promoter",   "upstream",   2000),
                   ("promoter",   "downstream", 2000)],
    "intron":     [("cds",        "upstream",   3000),  # CDS containing the intron
                   ("cds",        "downstream", 3000)],
}


def _classify_partner(ann):
    """Lightweight kind classifier for partner search — operates directly
    on the feature dict (no need to consult interactions)."""
    t = (ann.get("type") or "").lower()
    nm = (ann.get("name") or "").lower()
    if t in ("cds", "gene", "protein", "marker", "protein_generator"):
        # Reject obvious non-CDS annotations whose type happens to be 'cds'
        # via misclassification (selection markers count as cds for our
        # orientation purposes).
        return "cds"
    if t == "promoter" or "promoter" in nm:
        return "promoter"
    if t == "polya_signal" or "polya" in t or "polya" in nm or "poly(a)" in nm:
        return "polya"
    if t == "terminator" or "terminator" in nm:
        return "terminator"
    if t == "rbs" or "rbs" in nm or "shine" in nm:
        return "rbs"
    if t == "enhancer" or "enhancer" in nm:
        return "enhancer"
    if t == "intron" or "intron" in nm:
        return "intron"
    return None


def _has_partner_on_strand(feature_ann, partner_kind, side, max_gap, annotations,
                             want_strand):
    """Return True if any feature of `partner_kind` sits on the
    `want_strand` side of `feature_ann` within `max_gap` bp."""
    fs, fe, f_strand = _ann_span(feature_ann)
    for ann in annotations:
        if ann is feature_ann:
            continue
        k = _classify_partner(ann)
        if k != partner_kind:
            continue
        ps, pe, p_strand = _ann_span(ann)
        if p_strand != want_strand:
            continue
        # Side is relative to feature's transcription direction.
        # f_strand == +1: upstream = lower coords; downstream = higher
        # f_strand == -1: upstream = higher coords; downstream = lower
        if side == "upstream":
            if f_strand >= 0:
                gap = fs - pe
            else:
                gap = ps - fe
        else:  # downstream
            if f_strand >= 0:
                gap = ps - fe
            else:
                gap = fs - pe
        if -10 <= gap <= max_gap:
            return ann
    return None


def _detect_orientation_mismatch(feature_ann, kind, annotations):
    """If the feature is missing its same-strand partner BUT an opposite-
    strand partner of the same kind exists within the window, return the
    opposite-strand partner (the orientation-mismatch evidence). Otherwise
    return None."""
    fs, fe, f_strand = _ann_span(feature_ann)
    partnerships = _ORIENTATION_PARTNERSHIPS.get(kind, [])
    for partner_kind, side, max_gap in partnerships:
        # Skip if same-strand partner exists — feature is fine
        if _has_partner_on_strand(feature_ann, partner_kind, side, max_gap,
                                    annotations, want_strand=f_strand):
            continue
        # Look for opposite-strand partner in the same window
        opp_partner = _has_partner_on_strand(
            feature_ann, partner_kind, side, max_gap, annotations,
            want_strand=-f_strand if f_strand else 1,
        )
        if opp_partner is not None:
            return {
                "partner": opp_partner,
                "partner_kind": partner_kind,
                "side": side,
                "max_gap": max_gap,
            }
    return None


def _check_cds_functional_context(
    product_ctx: "PlasmidContext",
    product_interactions: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """For every expression-feature in the product, verify it is wired into
    the rule-schema interaction graph in the appropriate role:
      - cds        -> template participant
      - promoter   -> stimulator participant
      - polya      -> modifier participant AND same strand as nearest active
                      template CDS (orientation check)
      - terminator -> modifier participant
      - rbs        -> stimulator or modifier participant
      - enhancer   -> orientation-only: same strand as nearest active
                      stimulator promoter within `_SUPPORT_FEATURE_WINDOW`
      - intron     -> orientation-only: same rule as enhancer
    The check runs strictly on the already-emitted interaction graph; it does
    not re-run Step 1 of the annotation pipeline.
    """
    warnings: List[Dict[str, Any]] = []
    annotations = product_ctx.annotations or []
    modules = product_ctx.modules or []
    interactions = product_interactions or []
    expected_role_map = dict(EXPRESSION_FEATURE_KINDS)
    active = _collect_active_participants(interactions)

    for ann in annotations:
        kind = _classify_expression_feature(ann)
        if kind is None:
            continue
        expected_roles = expected_role_map.get(kind, ())
        fs, fe, strand = _ann_span(ann)
        name = ann.get("name") or ann.get("feature_name") or f"unnamed_{kind}"

        gaps: List[str] = []
        ix_summary: Dict[str, List[str]] = {}

        if expected_roles:
            by_role = _expression_interactions_for_feature(ann, interactions, expected_roles)
            ix_summary = {
                role: [ix.get("rule_id") for ix in by_role.get(role, [])]
                for role in expected_roles
            }
            any_hit = any(by_role.get(r) for r in expected_roles)
            if not any_hit:
                gaps.append(
                    "no interaction with this feature as " +
                    " or ".join(expected_roles)
                )

        # v10: containment check for tag features (annotation-only).
        if kind == "tag":
            tag_label = _classify_tag(ann) or name
            host_cds = None
            for a2 in annotations:
                if not _is_cds_class(a2):
                    continue
                cs2, ce2, str2 = _ann_span(a2)
                if cs2 <= fs and fe <= ce2:
                    host_cds = (a2, cs2, ce2, str2)
                    break
            if host_cds is None:
                gaps.append(
                    f"{tag_label} tag is not contained inside any CDS-class "
                    f"feature — orphan tag, will not be expressed as a fusion"
                )
            else:
                _, _, _, host_strand = host_cds
                if host_strand != strand:
                    gaps.append(
                        f"{tag_label} tag is on strand "
                        f"{('+' if strand == 1 else '-')} but its parent CDS "
                        f"'{host_cds[0].get('name')}' is on strand "
                        f"{('+' if host_strand == 1 else '-')} (orientation mismatch)"
                    )
        # v10 (cont): recombination sites must appear as `reactant` in an
        # INT-REC-* interaction. Filter the role check above to that family.
        elif kind == "recombination_site":
            rec_ix = [ix for ix in interactions
                      if str(ix.get("rule_id") or "").startswith("INT-REC-")]
            by_role = _expression_interactions_for_feature(ann, rec_ix, ("reactant",))
            ix_summary = {"reactant": [ix.get("rule_id") for ix in by_role.get("reactant", [])]}
            if not by_role.get("reactant"):
                # Replace the generic gap from the participant check above so
                # the message names the recombination family.
                gaps = [
                    "no INT-REC-* interaction with this site as reactant — "
                    "unpaired recombination site (expected loxP×2 / FRT×2 / "
                    "attL×attR / attB×attP pair)"
                ]
        # v10 (cont): 2A ribosomal-skip peptides should appear as `stimulator`
        # in INT-CDS-2A-01.
        elif kind == "polyprotein_2a":
            tA_ix = [ix for ix in interactions
                     if str(ix.get("rule_id") or "").startswith("INT-CDS-2A-")]
            by_role = _expression_interactions_for_feature(ann, tA_ix, ("stimulator",))
            ix_summary = {"stimulator": [ix.get("rule_id") for ix in by_role.get("stimulator", [])]}
            if not by_role.get("stimulator"):
                gaps = [
                    "no INT-CDS-2A-* interaction with this 2A peptide as "
                    "stimulator — fused CDSes will not be co-translationally "
                    "separated"
                ]

        # v9 orientation rules:
        # - polyA orientation is enforced at interaction-calling time
        #   (_find_compatible_polyas only pairs same-strand polyAs with the
        #   cassette). The modifier-participant check above already encodes
        #   "presence + correct orientation".
        # - enhancer / intron: orientation is checked against the stimulator
        #   promoter paired with the active template CDS NEAREST to the
        #   feature — i.e. the promoter inside the upstream regulatory
        #   module of the same cassette.
        partner = None
        if kind in ("enhancer", "intron"):
            partner = _paired_promoter_for_nearest_cds(fs, fe, interactions)
            if partner is None:
                gaps.append(
                    "no active template CDS to identify the upstream "
                    "regulatory module promoter for orientation comparison"
                )
            elif _participant_strand(partner) != strand:
                gaps.append(
                    f"{kind} on strand {('+' if strand == 1 else '-')} but "
                    f"the upstream-regulatory promoter '{partner.get('name')}' "
                    f"of the nearest active CDS is on strand "
                    f"{('+' if _participant_strand(partner) == 1 else '-')} "
                    f"(orientation mismatch)"
                )

        if not gaps:
            continue

        # Diagnostic: which (if any) module wraps this feature.
        containing = None
        for m in modules:
            ms = int(m.get("start", 0) or 0)
            me = int(m.get("end", 0) or 0)
            if ms <= fs and fe <= me:
                containing = m.get("module_type")
                break

        # 2026-05-12: orientation-mismatch upgrade. If the missing-interaction
        # warning is explainable by an opposite-strand partner within
        # proximity, the natural fix is to reverse-complement the feature.
        # Promote the issue to 'orientation_mismatch' so the auto-correct
        # uses action='reverse_complement_feature' instead of
        # 'reverse_or_recheck' (both are in _AUTOCORRECT_ACTIONS, but
        # 'reverse_complement_feature' targets the right span more directly).
        _orient_evidence = _detect_orientation_mismatch(ann, kind, annotations)
        if _orient_evidence is not None:
            _p = _orient_evidence["partner"]
            _ps, _pe, _pstr = _ann_span(_p)
            gaps = list(gaps) + [
                f"opposite-strand {_orient_evidence['partner_kind']} "
                f"'{_p.get('name')}' within {_orient_evidence['max_gap']} bp on "
                f"strand {('+' if _pstr == 1 else '-')} — orientation mismatch "
                f"(reverse-complement this feature to form the expected cassette)"
            ]
            _issue = "orientation_mismatch"
            _orient_partner_summary = {
                "name": _p.get("name"),
                "kind": _orient_evidence["partner_kind"],
                "strand": _pstr,
                "start": _ps,
                "end": _pe,
                "side": _orient_evidence["side"],
            }
        else:
            _issue = "no_expression_interaction"
            _orient_partner_summary = None

        warnings.append({
            "feature": name,
            "feature_kind": kind,
            "expected_role": list(expected_roles),
            "feature_span": [fs, fe],
            "feature_strand": strand,
            "containing_module": containing,
            "expression_interactions": ix_summary,
            "nearest_active_partner": (
                {
                    "name": partner.get("name"),
                    "strand": _participant_strand(partner),
                    "start": partner.get("start"),
                    "end": partner.get("end"),
                } if partner else None
            ),
            "orientation_partner": _orient_partner_summary,
            "issue": _issue,
            "gaps": gaps,
            "remediation": (
                "For CDS/promoter/polyA/terminator/RBS: ensure the "
                "rule-schema cassette assembler emits an expression "
                "interaction with this feature in the expected role. "
                "For enhancer/intron: orient on the same strand as the "
                "active upstream-regulatory promoter."
            ),
            "detail": (
                f"{kind} '{name}' ({fs}..{fe}, strand "
                f"{'+' if strand == 1 else '-'}): " + "; ".join(gaps)
            ),
        })

    # v11: also emit module-level rule-validation warnings.
    warnings.extend(_check_module_validation(modules, interactions))

    return {
        "passed": len(warnings) == 0,
        "warnings": warnings,
    }


def _loose_match(product_ctx: PlasmidContext, target_ctx: TargetContext,
                 product_interactions: Optional[List[Dict[str, Any]]] = None,
                 target_interactions: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
    """Loose (functional) validation: re-annotate the assembled product and
    compare against target with per-class containment thresholds.

    Thresholds from TARGET_FROM_INVENTORY_ROUTING.md §Validation:
      - module_set containment ≥ 1.0 (every target module_type must exist in product)
      - feature_set containment ≥ 0.90 (named features with pi_permatch ≥ 90%)
      - interaction_set containment ≥ 1.0 (every target interaction must hold on product)
    """
    target_modules = _module_type_set(target_ctx.modules)
    product_modules = _module_type_set(product_ctx.modules)
    missing_modules = target_modules - product_modules
    module_pass = len(missing_modules) == 0

    target_features = _feature_name_set(target_ctx.annotations)
    product_features = _feature_name_set(product_ctx.annotations)
    if target_features:
        feature_containment = len(target_features & product_features) / len(target_features)
    else:
        feature_containment = 1.0
    feature_pass = feature_containment >= 0.90

    target_ixs = _interaction_set(target_interactions or [])
    product_ixs = _interaction_set(product_interactions or [])
    missing_ixs = target_ixs - product_ixs
    interaction_pass = len(missing_ixs) == 0

    functional_report = _check_cds_functional_context(product_ctx, product_interactions)
    functional_pass = functional_report["passed"]
    functional_warnings = functional_report["warnings"]

    passed = module_pass and feature_pass and interaction_pass and functional_pass

    return {
        "mode": "loose",
        "passed": passed,
        "module_containment": 1.0 if module_pass else (
            len(target_modules & product_modules) / max(1, len(target_modules))
        ),
        "feature_containment": round(feature_containment, 3),
        "interaction_containment": 1.0 if interaction_pass else (
            len(target_ixs & product_ixs) / max(1, len(target_ixs))
        ),
        "missing_modules": sorted(missing_modules),
        "missing_interactions": sorted("{}|{}".format(r, t) for r, t in missing_ixs),
        "functional_pass": functional_pass,
        "functional_warnings": functional_warnings,
        "reason": (
            "modules, features, interactions, and CDS functional context all check out"
            if passed else
            "missing: " + ", ".join(filter(None, [
                f"{len(missing_modules)} module_type(s)" if missing_modules else "",
                f"feature_containment {feature_containment:.2f} < 0.90" if not feature_pass else "",
                f"{len(missing_ixs)} interaction(s)" if missing_ixs else "",
                (f"{len(functional_warnings)} CDS functional issue(s): "
                 + "; ".join(f'{w["feature"]}->{w["issue"]}' for w in functional_warnings))
                if functional_warnings else "",
            ]))
        ),
    }


def _suggestion_for_warning(w: Dict[str, Any]) -> Dict[str, str]:
    """Generate a concrete remediation suggestion for a verifier warning."""
    issue = (w.get("issue") or "").lower()
    kind = (w.get("feature_kind") or "").lower()
    feat = w.get("feature") or w.get("module_name") or "feature"

    if issue == "orientation_mismatch":
        partner = w.get("orientation_partner") or {}
        return {
            "kind": kind or "feature",
            "subject": str(feat),
            "action": "reverse_complement_feature",
            "suggestion": (
                f"{kind} '{feat}' has its expected {partner.get('kind','partner')} "
                f"'{partner.get('name','?')}' on the OPPOSITE strand within "
                f"~{w.get('orientation_partner', {}).get('side','')}. "
                f"Reverse-complementing the {kind} forms the expected expression cassette."
            ),
        }

    if issue == "missing_module_interaction":
        return {
            "kind": "module",
            "subject": str(w.get("module_type") or feat),
            "action": "fix_rule_emission",
            "suggestion": (
                f"Module '{w.get('module_type')}' was detected but the rule "
                f"schema did not emit '{w.get('rule_id')}'. Verify the "
                f"sub-feature(s) the rule expects are present and on the "
                f"same strand; otherwise the cassette assembler will not "
                f"produce the expected interaction."
            ),
        }

    if kind == "cds":
        return {
            "kind": "cds",
            "subject": str(feat),
            "action": "reverse_or_recheck",
            "suggestion": (
                f"CDS '{feat}' has no expression interaction. Most likely "
                f"causes: insert is in reverse orientation (try "
                f"reverse-complementing the insert), or the upstream "
                f"regulatory module / Kozak is missing. Reassemble with "
                f"the CDS on the same strand as a Pol II / Pol III / T7 "
                f"promoter and a polyA / terminator downstream."
            ),
        }
    if kind == "promoter":
        return {
            "kind": "promoter",
            "subject": str(feat),
            "action": "pair_with_cds",
            "suggestion": (
                f"Promoter '{feat}' drives no detected CDS. Add a CDS "
                f"downstream on the same strand within ~3 kb, or remove "
                f"the orphan promoter."
            ),
        }
    if kind == "polya":
        return {
            "kind": "polya",
            "subject": str(feat),
            "action": "reorient_or_pair",
            "suggestion": (
                f"PolyA '{feat}' is not paired with a template CDS. "
                f"Confirm it is on the same strand as the upstream "
                f"expression cassette; if so, ensure the upstream CDS is "
                f"detected (else add a CDS module 5' of this polyA on the "
                f"same strand)."
            ),
        }
    if kind == "terminator":
        return {
            "kind": "terminator",
            "subject": str(feat),
            "action": "reorient_or_pair",
            "suggestion": (
                f"Terminator '{feat}' is not modifying any expression "
                f"cassette. Likely orphan or wrong-strand — flip to match "
                f"the upstream CDS strand or remove."
            ),
        }
    if kind in ("enhancer", "intron"):
        gap = " ".join(w.get("gaps") or [])
        if "orientation mismatch" in gap:
            return {
                "kind": kind,
                "subject": str(feat),
                "action": "reverse_complement_feature",
                "suggestion": (
                    f"Reverse-complement {kind} '{feat}' to match the "
                    f"strand of the paired upstream-regulatory promoter "
                    f"of the nearest CDS."
                ),
            }
        return {
            "kind": kind,
            "subject": str(feat),
            "action": "place_within_cassette",
            "suggestion": (
                f"Place {kind} '{feat}' within ~3 kb of an active "
                f"upstream-regulatory promoter on the same strand."
            ),
        }
    if kind == "rbs":
        return {
            "kind": "rbs",
            "subject": str(feat),
            "action": "add_or_pair_cds",
            "suggestion": (
                f"RBS '{feat}' is not driving translation of any CDS. "
                f"Place a CDS immediately downstream (within ~12 bp) on "
                f"the same strand."
            ),
        }
    if kind == "recombination_site":
        return {
            "kind": "recombination_site",
            "subject": str(feat),
            "action": "add_partner_site",
            "suggestion": (
                f"Recombination site '{feat}' is unpaired. Add the "
                f"complementary partner site (loxP×2, FRT×2, attL×attR, "
                f"or attB×attP) flanking the cassette to enable a "
                f"productive recombination reaction."
            ),
        }
    if kind == "polyprotein_2a":
        return {
            "kind": "polyprotein_2a",
            "subject": str(feat),
            "action": "place_between_cds",
            "suggestion": (
                f"2A peptide '{feat}' is not between two CDS submodules. "
                f"Place it in-frame between the two CDSes within a single "
                f"parent ORF; both flanking CDSes must be on the same "
                f"strand."
            ),
        }
    if kind == "tag":
        gap = " ".join(w.get("gaps") or [])
        if "orientation mismatch" in gap:
            return {
                "kind": "tag",
                "subject": str(feat),
                "action": "reorient_tag",
                "suggestion": (
                    f"Tag '{feat}' is on the opposite strand of its "
                    f"parent CDS. Reverse-complement the tag insert."
                ),
            }
        return {
            "kind": "tag",
            "subject": str(feat),
            "action": "place_inside_cds",
            "suggestion": (
                f"Tag '{feat}' is not contained in any CDS. Move it into "
                f"the target CDS (N-terminal or C-terminal, in-frame)."
            ),
        }

    # Fallback
    return {
        "kind": kind or "unknown",
        "subject": str(feat),
        "action": "review",
        "suggestion": w.get("remediation") or "Review feature placement and orientation.",
    }


def _revcomp(s: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
            "a": "t", "t": "a", "g": "c", "c": "g", "n": "n"}
    return "".join(comp.get(b, "N") for b in reversed(s))


_AUTOCORRECT_ACTIONS = {
    "reverse_or_recheck",
    "reverse_complement_insert",
    "reverse_complement_feature",
    "reorient_tag",
}


def _name_matches(feature_name: str, kb_names: List[str]) -> bool:
    fn = (feature_name or "").lower()
    if not fn:
        return False
    for n in kb_names or []:
        if not n:
            continue
        nl = n.lower()
        if nl == fn or nl in fn or fn in nl:
            return True
    return False


async def auto_correct_kb_part_orientation(
    target_sequence: str,
    target_annotations: List[Dict[str, Any]],
    target_modules: List[Dict[str, Any]],
    target_interactions: List[Dict[str, Any]],
    kb_resolved_part_names: List[str],
    reannotate_async,
) -> Dict[str, Any]:
    """Attempt orientation corrections for KB-referenced parts that fail
    post-assembly verification with a reverse-complement remediation.

    Args:
        target_sequence:           assembled target DNA (5' -> 3').
        target_annotations:        verifier-shape annotations on the target.
        target_modules:            rule-based modules on the target.
        target_interactions:       interaction graph for the target.
        kb_resolved_part_names:    canonical names of KB-referenced parts the
                                   user explicitly requested. ONLY warnings
                                   on features whose name matches one of these
                                   are eligible for auto-correction.
        reannotate_async:          async callable
                                   (sequence) -> (annotations, modules, interactions)
                                   used to rebuild the verifier inputs after
                                   a candidate correction.

    Returns:
        {
          "corrected_sequence":  Optional[str]   # final corrected target, or
                                                  None if nothing was applied
          "corrections":         [ ... ]         # per-attempt log
          "verification_before": Dict            # initial verifier verdict
          "verification_after":  Dict            # verdict on the (possibly
                                                  corrected) target
        }
    """
    initial = verify_target_design(
        target_sequence=target_sequence,
        target_annotations=target_annotations,
        target_modules=target_modules,
        target_interactions=target_interactions,
        target_name="target",
    )
    corrections: List[Dict[str, Any]] = []

    if initial.get("passed"):
        return {
            "corrected_sequence": None,
            "corrections": corrections,
            "verification_before": initial,
            "verification_after": initial,
        }

    current_seq = target_sequence
    current_ann = target_annotations
    current_mods = target_modules
    current_ix = target_interactions
    current_verification = initial
    applied_any = False

    for warning, suggestion in zip(
        list(current_verification.get("warnings") or []),
        list(current_verification.get("suggestions") or []),
    ):
        feat_name = warning.get("feature") or ""
        if not _name_matches(feat_name, kb_resolved_part_names):
            continue
        action = (suggestion.get("action") or "").lower()
        if action not in _AUTOCORRECT_ACTIONS:
            continue
        span = warning.get("feature_span") or warning.get("cds_span")
        if not span or span[0] is None or span[1] is None:
            continue
        s, e = int(span[0]), int(span[1])
        if e <= s or e > len(current_seq):
            continue

        candidate_seq = current_seq[:s] + _revcomp(current_seq[s:e]) + current_seq[e:]
        try:
            new_ann, new_mods, new_ix = await reannotate_async(candidate_seq)
        except Exception as exc:
            corrections.append({
                "feature": feat_name,
                "action": "reverse_complement",
                "span": [s, e],
                "result": "reverted",
                "reason": f"reannotation failed: {exc}",
            })
            continue

        new_verification = verify_target_design(
            target_sequence=candidate_seq,
            target_annotations=new_ann,
            target_modules=new_mods,
            target_interactions=new_ix,
            target_name="target",
        )
        # Did the targeted feature stop being flagged?
        still_flagged = any(
            (_name_matches(w.get("feature") or "", [feat_name]))
            and (w.get("issue") == warning.get("issue"))
            for w in (new_verification.get("warnings") or [])
        )
        # 2026-05-12: also accept when reversing the feature reduces the
        # total orientation_mismatch warning count across the plasmid —
        # the feature itself may still have a (different-issue) warning,
        # but if other features that pointed at it now form proper
        # same-strand cassettes the design is materially better.
        orient_before = sum(
            1 for w in (current_verification.get("warnings") or [])
            if (w.get("issue") or "").lower() == "orientation_mismatch"
        )
        orient_after = sum(
            1 for w in (new_verification.get("warnings") or [])
            if (w.get("issue") or "").lower() == "orientation_mismatch"
        )
        if still_flagged and orient_after >= orient_before:
            corrections.append({
                "feature": feat_name,
                "action": "reverse_complement",
                "span": [s, e],
                "result": "reverted",
                "reason": "feature still flagged after reverse-complement "
                          f"and orientation_mismatch count did not drop "
                          f"({orient_before} -> {orient_after})",
            })
            continue

        corrections.append({
            "feature": feat_name,
            "action": "reverse_complement",
            "span": [s, e],
            "result": "applied",
            "verdict_before": "FAIL",
            "verdict_after": (
                "PASS" if new_verification.get("passed") else "FAIL_other"
            ),
            "remaining_warnings_after": (
                new_verification.get("summary") or {}
            ).get("total_warnings", 0),
        })
        current_seq = candidate_seq
        current_ann = new_ann
        current_mods = new_mods
        current_ix = new_ix
        current_verification = new_verification
        applied_any = True

    return {
        "corrected_sequence": current_seq if applied_any else None,
        "corrections": corrections,
        "verification_before": initial,
        "verification_after": current_verification,
    }


def analyze_design_intent(
    intent_result: Dict[str, Any],
    user_message: str = "",
) -> Dict[str, Any]:
    """Classify the KB-resolved parts in a user's design request into the
    same expression-feature kinds the verifier uses, and report what is
    missing for a complete expression interaction network.

    Intended to run AFTER `parse_intent` and BEFORE predesign, so the
    dispatcher can warn the user up-front when the request lacks elements
    needed for expression (e.g. CDS without promoter, promoter without
    CDS, missing polyA, orphan tag).

    Returns:
        {
          "completeness":       "complete" | "incomplete" | "unknown",
          "kinds_present":      {kind: [feature_names]},
          "missing_kinds":      [kinds],
          "design_warnings":    [{kind, feature?, gap, suggestion}],
          "summary":            text summary line,
        }
    """
    kb_resolved = (intent_result or {}).get("kb_resolved") or {}
    identified = kb_resolved.get("identified") or []

    # Each identified item may carry feature name / type / kb_class.
    kinds_present: Dict[str, List[str]] = {}
    for item in identified:
        ann = {
            "name": item.get("feature_name") or item.get("name") or "",
            "type": item.get("type") or "",
            "kb_class": item.get("kb_class") or item.get("class") or "",
            "kb_subclass": item.get("kb_subclass") or item.get("subclass") or "",
        }
        kind = _classify_expression_feature(ann)
        if not kind:
            continue
        kinds_present.setdefault(kind, []).append(ann["name"] or "(unnamed)")

    # Also scan the raw message for keywords the KB lookup may have missed.
    low = (user_message or "").lower()
    if "cds" not in kinds_present:
        for token in ("gfp", "rfp", "yfp", "luciferase", "cas9", "cas12", "halotag",
                      "puror", "neor", "hygror", "blaster"):
            if token in low:
                kinds_present.setdefault("cds", []).append(token)
                break
    if "promoter" not in kinds_present and ("promoter" in low or "cmv" in low or "ef1" in low or "u6" in low or "t7" in low):
        kinds_present.setdefault("promoter", []).append("(implied)")
    if "polya" not in kinds_present and ("polya" in low or "poly(a)" in low or "polya signal" in low or "bgh" in low or "sv40" in low):
        kinds_present.setdefault("polya", []).append("(implied)")

    design_warnings: List[Dict[str, str]] = []

    # Required-pair contracts (role-based — same shape the verifier uses).
    has_cds = bool(kinds_present.get("cds"))
    has_promoter = bool(kinds_present.get("promoter"))
    has_polya = bool(kinds_present.get("polya")) or bool(kinds_present.get("terminator"))
    has_rbs = bool(kinds_present.get("rbs"))

    if has_cds and not has_promoter:
        design_warnings.append({
            "kind": "promoter",
            "feature": ", ".join(kinds_present["cds"]),
            "gap": "CDS named in request has no promoter — expression cassette incomplete",
            "suggestion": (
                "Specify a Pol II promoter (CMV / EF1α / chicken β-actin) for "
                "mammalian expression, a Pol III promoter (U6 / H1) for guide "
                "RNAs, or a T7 promoter for prokaryotic expression."
            ),
        })
    if has_cds and not has_polya:
        design_warnings.append({
            "kind": "polya",
            "feature": ", ".join(kinds_present["cds"]),
            "gap": "CDS named in request has no polyA / terminator",
            "suggestion": (
                "Specify a polyA signal (bGH / SV40 / hGH) for Pol II, or a "
                "terminator (rrnB T1T2 / T7 terminator) for prokaryotic "
                "expression."
            ),
        })
    if has_promoter and not has_cds and not has_rbs:
        design_warnings.append({
            "kind": "cds",
            "feature": ", ".join(kinds_present.get("promoter") or []),
            "gap": "Promoter named in request has no CDS or guide-RNA payload",
            "suggestion": (
                "Specify the protein/CDS to be expressed (e.g. GFP, Cas9) "
                "or the guide scaffold for Pol III cassettes."
            ),
        })
    if kinds_present.get("recombination_site") and len(kinds_present["recombination_site"]) < 2:
        design_warnings.append({
            "kind": "recombination_site",
            "feature": kinds_present["recombination_site"][0],
            "gap": "Single recombination site requested — needs a partner",
            "suggestion": (
                "Recombination requires a pair (loxP×2 / FRT×2 / attL×attR / "
                "attB×attP). Add the complementary site."
            ),
        })
    if kinds_present.get("tag") and not has_cds:
        design_warnings.append({
            "kind": "tag",
            "feature": kinds_present["tag"][0],
            "gap": "Fusion tag requested but no parent CDS",
            "suggestion": "Specify the CDS the tag should fuse to (N- or C-terminal).",
        })

    if not kinds_present:
        completeness = "unknown"
    elif design_warnings:
        completeness = "incomplete"
    else:
        completeness = "complete"

    summary = (
        f"design intent: {completeness} — "
        f"kinds present: {{{', '.join(sorted(kinds_present))}}}"
        + (f"; {len(design_warnings)} gap(s)" if design_warnings else "")
    )

    return {
        "completeness": completeness,
        "kinds_present": kinds_present,
        "missing_kinds": [w["kind"] for w in design_warnings],
        "design_warnings": design_warnings,
        "summary": summary,
    }


def verify_target_design(
    target_sequence: str,
    target_annotations: List[Dict[str, Any]],
    target_modules: List[Dict[str, Any]],
    target_interactions: List[Dict[str, Any]],
    target_name: str = "target",
) -> Dict[str, Any]:
    """Run the design verifier against an assembled target plasmid and
    decorate the warnings with concrete remediation suggestions.

    Returns:
        {
          "passed":    bool,
          "warnings":  List[Dict] (raw verifier warnings),
          "suggestions": List[Dict] (per-warning remediation suggestions),
          "summary":   {"total_warnings": int, "by_kind": {kind: count}},
        }
    """
    ctx = PlasmidContext(
        name=target_name,
        sequence=(target_sequence or "").upper(),
        gb_text=None,
        annotations=target_annotations or [],
        modules=target_modules or [],
        cloning=None,
    )
    report = _check_cds_functional_context(ctx, target_interactions or [])
    warnings = report.get("warnings") or []
    suggestions = [_suggestion_for_warning(w) for w in warnings]

    by_kind: Dict[str, int] = {}
    for w in warnings:
        k = w.get("feature_kind") or w.get("module_type") or "unknown"
        by_kind[k] = by_kind.get(k, 0) + 1

    return {
        "passed": report.get("passed", False),
        "warnings": warnings,
        "suggestions": suggestions,
        "summary": {"total_warnings": len(warnings), "by_kind": by_kind},
    }


async def validate(
    chosen: FeasibilityReport,
    target_ctx: TargetContext,
    product_sequence: Optional[str],
    product_gb_text: Optional[str] = None,
    product_interactions: Optional[List[Dict[str, Any]]] = None,
    target_interactions: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Validate an in-silico assembled product against the target.

    The validation mode is taken from the chosen FeasibilityReport:
      - strict workflows (gateway_cloning, sdm_design, sgrna_golden_gate, pcr_design,
        restriction_cloning when all fragments are sequence-identical sources):
        canonical-rotate + sequence-equal.
      - loose workflows (inv_gib, golden_gate_primer_design, synthesis_fallback):
        re-annotate product, compare module/feature/interaction sets.
    """
    mode = chosen.validation_mode if chosen else "loose"

    if not product_sequence:
        return {"mode": mode, "passed": False, "reason": "no product sequence provided"}

    if mode == "strict":
        return _strict_match(product_sequence, target_ctx.sequence)

    # Loose path: annotate the product and compare to target annotation sets.
    product_ctx = await annotate_one("product", product_gb_text, product_sequence)
    return _loose_match(product_ctx, target_ctx,
                        product_interactions=product_interactions,
                        target_interactions=target_interactions)


# ---------------------------------------------------------------------------
# Routing + audit reply
# ---------------------------------------------------------------------------

def route(target: TargetContext, inventory: List[InventoryContext]) -> Tuple[Optional[FeasibilityReport], List[FeasibilityReport]]:
    reports = [
        assess_sdm_feasibility(target, inventory),
        assess_sgrna_golden_gate_feasibility(target, inventory),
        assess_gateway_feasibility(target, inventory),  # unified: single + multisite
        assess_restriction_feasibility(target, inventory),
        assess_golden_gate_feasibility(target, inventory),
        assess_inv_gib_feasibility(target, inventory),
        assess_pcr_extension_gibson_feasibility(target, inventory),
        assess_synthesis_fallback_feasibility(target, inventory),  # always feasible
    ]
    feasible = [r for r in reports if r.feasible]
    feasible.sort(key=lambda r: (r.score, -r.work_estimate, r.success_estimate), reverse=True)
    chosen = feasible[0] if feasible else None
    return chosen, reports


def build_audit_markdown(reports: List[FeasibilityReport], chosen: Optional[FeasibilityReport]) -> str:
    lines = ["## Target-from-Inventory Routing Audit", ""]
    if chosen:
        lines.append(f"**Chosen workflow:** `{chosen.workflow}` (score={chosen.score:.3f}, validation={chosen.validation_mode})")
        lines.append(f"**Why:** {chosen.rationale}")
        if chosen.workflow == "synthesis_fallback":
            ha = chosen.handler_args or {}
            lines.append(
                f"**Synth manifest:** {len(ha.get('synth_blocks', []))} block(s), "
                f"{ha.get('total_synth_bp', 0)} bp total, est. ${ha.get('est_cost_usd', 0):.0f}"
            )
    else:
        lines.append("**No feasible workflow found.**")
    lines.append("")
    lines.append("| Workflow | Feasible | Score | Work | Success | Rationale |")
    lines.append("|---|---|---|---|---|---|")
    for r in sorted(reports, key=lambda x: (-int(x.feasible), -x.score)):
        flag = "✅" if r.feasible else "❌"
        lines.append(f"| {r.workflow} | {flag} | {r.score:.3f} | {r.work_estimate} | {r.success_estimate:.2f} | {r.rationale[:120]} |")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# UploadFile entry point used by chat.py
# ---------------------------------------------------------------------------

async def route_from_uploads(
    target_upload: UploadFile,
    inventory_uploads: Sequence[UploadFile],
) -> Dict[str, Any]:
    """Read uploads, annotate, route. Returns a dict with chosen + audit + contexts."""
    target_raw = await target_upload.read()
    target_name, target_seq = _read_seq_from_upload(target_raw, target_upload.filename)
    try:
        target_gb_text = target_raw.decode("utf-8")
    except Exception:
        target_gb_text = None

    # Reset pointer so downstream re-reads of the same UploadFile work
    try:
        await target_upload.seek(0)
    except Exception:
        pass

    target_ctx = await annotate_one(target_name, target_gb_text, target_seq)

    inventory_ctxs: List[InventoryContext] = []
    for inv in inventory_uploads:
        inv_raw = await inv.read()
        inv_name, inv_seq = _read_seq_from_upload(inv_raw, inv.filename)
        try:
            inv_gb_text = inv_raw.decode("utf-8")
        except Exception:
            inv_gb_text = None
        try:
            await inv.seek(0)
        except Exception:
            pass
        inventory_ctxs.append(await annotate_one(inv_name, inv_gb_text, inv_seq))

    chosen, all_reports = route(target_ctx, inventory_ctxs)
    audit_md = build_audit_markdown(all_reports, chosen)

    return {
        "ok": True,
        "intent": "target_from_inventory",
        "chosen_intent": chosen.workflow if chosen else None,
        "chosen_handler_args": chosen.handler_args if chosen else {},
        "chosen_validation_mode": chosen.validation_mode if chosen else None,
        "audit_markdown": audit_md,
        "reports": [asdict(r) for r in all_reports],
        "target_context": {"name": target_ctx.name, "length": len(target_ctx.sequence)},
        "inventory_contexts": [{"name": c.name, "length": len(c.sequence)} for c in inventory_ctxs],
    }
