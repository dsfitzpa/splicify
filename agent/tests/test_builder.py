"""Builder + verifier integration tests using synthetic fixtures.

For each cloning method (Gibson / Golden Gate / restriction), build
a PartSet whose junctions encode that method's compatibility, then
also build broken variants (misoriented, wrong order, missing role)
to confirm the verifier surfaces the right diagnostic codes."""
from __future__ import annotations

import asyncio

import agent_v2  # noqa: F401 — path shim
from agent_v2.builder.intent_spec import IntentSpec
from agent_v2.builder.part_set import Part, PartSet
from agent_v2.builder.virtual_construct import (
    Slot, VirtualConstruct, assess_methods,
    GIBSON_MIN_OVERLAP, TYPE_IIS_SITES, COMMON_TYPE_II,
)
from agent_v2.builder.verifier import (
    DIAG_NO_METHOD, DIAG_MODULE_MISSING, DIAG_ORIENTATION_WRONG,
    DIAG_ORDER_WRONG, DIAG_ROLE_MISSING,
)


# ─────────────────────────────────────────────────────────────────
# Fixture helpers
# ─────────────────────────────────────────────────────────────────
def mk_part(name: str, role: str, body: str, *, upstream: str = "A"*50,
            downstream: str = "T"*50, strand: int = 1, pid: str = "p1",
            start: int = 0) -> Part:
    return Part(
        name=name, role=role, source_plasmid_id=pid,
        source_start=start, source_end=start + len(body),
        source_strand=strand,
        body_sequence=body.upper(),
        upstream_junction=upstream.upper(),
        downstream_junction=downstream.upper(),
    )


def gibson_set(body_a="ATG"*20, body_b="GCT"*20, body_c="TCA"*20) -> tuple[Part, Part, Part]:
    """Three parts whose junctions overlap by 20 bp each (gibson-ready).
    Layout: [..AA..][overlap1][body_a][overlap2][body_b][overlap3][body_c][overlap1..circular]"""
    o1 = "ACGTACGTACGTACGTACGT"     # 20 bp, between C and A (circular)
    o2 = "GATCGATCGATCGATCGATC"     # between A and B
    o3 = "CCCAAATTTGGGAAACCCAA"     # between B and C
    a = mk_part("partA", "promoter", body_a,
                 upstream=o1 + "N"*30, downstream=o2 + "N"*30, pid="src1", start=100)
    b = mk_part("partB", "cds", body_b,
                 upstream=o2 + "N"*30, downstream=o3 + "N"*30, pid="src1", start=200)
    c = mk_part("partC", "polya", body_c,
                 upstream=o3 + "N"*30, downstream=o1 + "N"*30, pid="src1", start=300)
    return a, b, c


def gg_set() -> tuple[Part, Part, Part]:
    """Three-part Golden Gate (BsaI flanks every junction)."""
    bsa = "GGTCTC"
    a = mk_part("u6", "promoter", "ATG"*15,
                 upstream="N"*44 + bsa, downstream=bsa + "N"*44)
    b = mk_part("stuffer", "stuffer", "GCT"*15,
                 upstream=bsa + "N"*44, downstream=bsa + "N"*44)
    c = mk_part("scaffold", "scaffold", "TCA"*15,
                 upstream=bsa + "N"*44, downstream="N"*44 + bsa)
    return a, b, c


# ─────────────────────────────────────────────────────────────────
# assess_methods tests (deterministic, no annotation call)
# ─────────────────────────────────────────────────────────────────
def test_assess_picks_gibson_when_junctions_overlap():
    a, b, c = gibson_set()
    vc = VirtualConstruct(slots=[Slot(p) for p in [a, b, c]], topology="circular")
    m = assess_methods(vc)
    assert "gibson" in m.feasible
    assert m.pick == "gibson"


def test_assess_picks_sgrna_gg_when_pol3_stuffer_scaffold():
    a, b, c = gg_set()
    # Mark role as pol3_promoter so the sgrna_gg special case triggers
    a.role = "pol3_promoter"
    vc = VirtualConstruct(slots=[Slot(p) for p in [a, b, c]], topology="circular")
    m = assess_methods(vc)
    assert "golden_gate" in m.feasible
    assert m.pick == "sgrna_gg"


def test_assess_rejects_when_no_method_feasible():
    """Junctions with NO homology + NO Type IIs + NO restriction site."""
    # Use generic A/T-only junctions that share <15 bp homology
    a = mk_part("a", "promoter", "ATGATGATG", upstream="A"*50, downstream="C"*50)
    b = mk_part("b", "cds", "GCTGCTGCT", upstream="G"*50, downstream="T"*50)
    vc = VirtualConstruct(slots=[Slot(a), Slot(b)], topology="circular")
    m = assess_methods(vc)
    assert m.pick is None
    assert "gibson" in m.rejected


# ─────────────────────────────────────────────────────────────────
# Verifier tests with mocked annotate_fn
# ─────────────────────────────────────────────────────────────────
def _make_ann_fn(modules: list[str], interactions: list[str] = None):
    """Return an async annotate_fn that yields the requested modules /
    interactions every call."""
    interactions = interactions or []
    async def _fn(seq, *, circular=True):
        return {
            "sequence": seq, "circular": circular,
            "annotations": [{"name": "x", "start": 0, "end": 10, "direction": 1}],
            "modules": [{"module_type": m, "start": 0, "end": 10, "strand": 1, "name": m}
                        for m in modules],
            "interactions": [{"type": it} for it in interactions],
        }
    return _fn


def test_verifier_passes_when_expression_cassette_present():
    from agent_v2.builder.verifier import verify
    # Use a minimal intent so the synthetic fixture's promoter+cds+polya
    # parts are enough — no need for selection_marker / origin.
    intent = IntentSpec(
        function="expression",
        required_modules=["mammalian_pol2_expression_cassette"],
        required_interactions=[],
        required_roles=["promoter", "cds", "polya"],
    )
    a, b, c = gibson_set()
    a.role, b.role, c.role = "promoter", "cds", "polya"
    vc = VirtualConstruct(slots=[Slot(p, orientation=1) for p in [a, b, c]],
                            topology="circular")
    ann = _make_ann_fn(["mammalian_pol2_expression_cassette"])
    res = asyncio.run(verify(vc, intent, annotate_fn=ann))
    assert res.passed, res.diagnostics


def test_verifier_flags_missing_module():
    from agent_v2.builder.verifier import verify
    intent = IntentSpec.for_expression_cassette()
    a, b, c = gibson_set()
    a.role, b.role, c.role = "promoter", "cds", "polya"
    vc = VirtualConstruct(slots=[Slot(p, orientation=1) for p in [a, b, c]],
                            topology="circular")
    ann = _make_ann_fn([])     # no modules detected
    res = asyncio.run(verify(vc, intent, annotate_fn=ann))
    assert not res.passed
    assert any(d.code == DIAG_MODULE_MISSING for d in res.diagnostics)


def test_verifier_flags_orientation_when_promoter_opposite_to_cds():
    from agent_v2.builder.verifier import verify
    intent = IntentSpec.for_expression_cassette()
    a, b, c = gibson_set()
    a.role, b.role, c.role = "promoter", "cds", "polya"
    vc = VirtualConstruct(slots=[Slot(a, orientation=-1), Slot(b, orientation=1),
                                  Slot(c, orientation=1)], topology="circular")
    ann = _make_ann_fn(["mammalian_pol2_expression_cassette"])
    res = asyncio.run(verify(vc, intent, annotate_fn=ann))
    assert any(d.code == DIAG_ORIENTATION_WRONG for d in res.diagnostics)


def test_verifier_flags_role_missing():
    from agent_v2.builder.verifier import verify
    intent = IntentSpec.for_expression_cassette()
    # Drop polya
    a = mk_part("partA", "promoter", "ATG"*15)
    b = mk_part("partB", "cds", "GCT"*15)
    vc = VirtualConstruct(slots=[Slot(a), Slot(b)], topology="circular")
    ann = _make_ann_fn(["mammalian_pol2_expression_cassette"])
    res = asyncio.run(verify(vc, intent, annotate_fn=ann))
    assert any(d.code == DIAG_ROLE_MISSING and "polya" in (d.detail or "")
                for d in res.diagnostics)


# ─────────────────────────────────────────────────────────────────
# Builder loop tests
# ─────────────────────────────────────────────────────────────────
def _min_intent():
    return IntentSpec(
        function="expression",
        required_modules=["mammalian_pol2_expression_cassette"],
        required_interactions=[],
        required_roles=["promoter", "cds", "polya"],
    )


def test_builder_passes_on_happy_path():
    from agent_v2.builder.builder import build
    intent = _min_intent()
    a, b, c = gibson_set()
    a.role, b.role, c.role = "promoter", "cds", "polya"
    ps = PartSet(parts=[a, b, c])
    ann = _make_ann_fn(["mammalian_pol2_expression_cassette"])
    res = asyncio.run(build(ps, intent, annotate_fn=ann))
    assert res.success
    assert res.method_pick == "gibson"


def test_builder_fixes_orientation_in_loop():
    from agent_v2.builder.builder import build
    intent = _min_intent()
    a, b, c = gibson_set()
    a.role, b.role, c.role = "promoter", "cds", "polya"
    a.source_strand = -1     # backwards promoter — builder should flip
    ps = PartSet(parts=[a, b, c])
    ann = _make_ann_fn(["mammalian_pol2_expression_cassette"])
    res = asyncio.run(build(ps, intent, annotate_fn=ann, max_iters=5))
    assert res.success, res.unresolved_diagnostics


def test_builder_fails_when_no_assembly_method():
    from agent_v2.builder.builder import build
    intent = _min_intent()
    a = mk_part("a", "promoter", "ATG"*15, upstream="A"*50, downstream="C"*50)
    b = mk_part("b", "cds", "GCT"*15, upstream="G"*50, downstream="T"*50)
    c = mk_part("c", "polya", "TCA"*15, upstream="N"*50, downstream="N"*50)
    ps = PartSet(parts=[a, b, c])
    ann = _make_ann_fn(["mammalian_pol2_expression_cassette"])
    res = asyncio.run(build(ps, intent, annotate_fn=ann, max_iters=3))
    # Even if modules pass, junction-method check should fail.
    assert not res.success
    assert any(d.code == DIAG_NO_METHOD for d in res.unresolved_diagnostics)
