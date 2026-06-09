"""PE3b ngRNAs bypass the 32-100 bp distance band.

Anzalone 2019 + Li 2021 (easy_prime): a PE3b ngRNA's spacer overlaps the
edited bases on the opposite strand, so once the pegRNA installs the
edit the ngRNA target carries a mismatch and Cas9 stops cutting. That
self-shutoff is what avoids the simultaneous-double-strand-break risk
the 32-100 bp band exists to mitigate. PE3b ngRNAs therefore work at
ANY distance — they're typically <32 bp from the pegRNA nick because
the spacer-edit overlap forces them close. The selector must let them
through regardless of distance.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

import pytest
from splicify_api.pegrna_designer import _select_ngrna, _Sgrna


def _sg(strand: str, start: int, end: int, cut_fwd: int, spacer: str = None) -> _Sgrna:
    return _Sgrna(
        spacer=(spacer or ("A" * 20)),
        pam="AGG",
        start=start, end=end,
        strand=strand, cut_fwd=cut_fwd,
        cas9_score=70.0,
    )


def _params():
    return {"min_ngRNA_distance": 32, "max_ngRNA_distance": 100,
             "max_max_ngRNA_distance": 200}


# ---------------------------------------------------------------------------
# PE3b distance exemption — the core fix
# ---------------------------------------------------------------------------
def test_pe3b_ngRNA_below_min_distance_is_accepted():
    """An opposite-strand ngRNA at 10 bp from the pegRNA nick, with its
    spacer overlapping the substitution, must be selected even though
    10 < min_ngRNA_distance=32."""
    peg = _sg(strand="+", start=900, end=920, cut_fwd=1000)
    edit_pos0 = 1005   # 1-bp substitution at plasmid position 1006
    # ngRNA on - strand, cut at +1010 (distance = 10 bp), spacer spans
    # 1003..1023 — overlaps the edit at 1005.
    ng_pe3b = _sg(strand="-", start=1003, end=1023, cut_fwd=1010)
    out = _select_ngrna(peg, [ng_pe3b], edit_pos0, "A", "C", _params())
    assert len(out) == 1
    sg, dist, is_pe3b, edited = out[0]
    assert is_pe3b == 1
    assert abs(dist) == 10
    assert edited is not None and edited != ng_pe3b.spacer
    # Edited spacer reflects the post-edit ngRNA target (mismatch).
    assert "A" in edited or "C" in edited


def test_pe3b_ngRNA_above_max_distance_is_accepted():
    """A PE3b candidate at >100 bp also passes. (Rare, but the rule is
    'overlap dictates inclusion, not distance'.)"""
    peg = _sg(strand="+", start=900, end=920, cut_fwd=1000)
    edit_pos0 = 1200
    ng_pe3b = _sg(strand="-", start=1195, end=1215, cut_fwd=1205)
    out = _select_ngrna(peg, [ng_pe3b], edit_pos0, "A", "C", _params())
    assert len(out) == 1
    assert out[0][2] == 1   # is_pe3b


def test_non_pe3b_below_min_distance_is_rejected():
    """A non-PE3b ngRNA (spacer doesn't overlap the edit) at 10 bp must
    still be rejected by the distance filter."""
    peg = _sg(strand="+", start=900, end=920, cut_fwd=1000)
    edit_pos0 = 1500   # edit is far from the ngRNA's spacer
    ng = _sg(strand="-", start=1003, end=1023, cut_fwd=1010)
    out = _select_ngrna(peg, [ng], edit_pos0, "A", "C", _params())
    assert out == []   # 10 bp away AND no overlap -> reject


def test_non_pe3b_within_band_is_accepted():
    """Sanity: 50 bp non-PE3b candidate is fine."""
    peg = _sg(strand="+", start=900, end=920, cut_fwd=1000)
    edit_pos0 = 1500
    ng = _sg(strand="-", start=1043, end=1063, cut_fwd=1050)
    out = _select_ngrna(peg, [ng], edit_pos0, "A", "C", _params())
    assert len(out) == 1
    assert out[0][2] == 0   # non-PE3b


def test_pe3b_preferred_alongside_non_pe3b_in_same_pass():
    """Both PE3 + PE3b candidates within their respective bands are
    returned in the same iteration; downstream re-rank picks PE3b first."""
    peg = _sg(strand="+", start=900, end=920, cut_fwd=1000)
    edit_pos0 = 1005   # ngRNA at 1003..1023 will overlap
    ng_pe3b = _sg(strand="-", start=1003, end=1023, cut_fwd=1010)   # 10 bp, PE3b
    ng_pe3 = _sg(strand="-", start=1043, end=1063, cut_fwd=1050)     # 50 bp, regular
    out = _select_ngrna(peg, [ng_pe3b, ng_pe3], edit_pos0, "A", "C", _params())
    assert len(out) == 2
    pe3b_flags = sorted(c[2] for c in out)
    assert pe3b_flags == [0, 1]


def test_pe3b_only_applies_to_substitutions():
    """Insertion / deletion edits have len(ref) != len(alt); the edited-spacer
    path requires len(ref)==len(alt)>0. A would-be PE3b candidate on an
    insertion thus falls through to the distance filter and is rejected
    when too close."""
    peg = _sg(strand="+", start=900, end=920, cut_fwd=1000)
    edit_pos0 = 1005
    ng = _sg(strand="-", start=1003, end=1023, cut_fwd=1010)
    # Insertion: ref="" (len 0), alt="A" (len 1) — different lengths.
    out = _select_ngrna(peg, [ng], edit_pos0, "", "A", _params())
    assert out == []   # not PE3b-eligible, fails 32 bp min


def test_pe3b_same_strand_still_excluded():
    """Strand rule still wins: PE3b requires opposite-strand spacer."""
    peg = _sg(strand="+", start=900, end=920, cut_fwd=1000)
    edit_pos0 = 1005
    same_strand_ng = _sg(strand="+", start=1003, end=1023, cut_fwd=1010)
    out = _select_ngrna(peg, [same_strand_ng], edit_pos0, "A", "C", _params())
    assert out == []
