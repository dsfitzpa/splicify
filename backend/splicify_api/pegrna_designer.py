"""
pegRNA designer -- prime-editing guide-RNA design for the plasmid viewer.

Pure-Python port of the easy_prime pipeline (Li et al. 2021,
https://github.com/YichaoOU/easy_prime). The original pipeline ships two
pickled XGBoost regressors -- PE2 (no nicking gRNA) and PE3 (with ngRNA).
We bundle the PE3 model and reproduce the feature set faithfully; the
DeepSpCas9 score (originally a flaky web call to deepcrispr.info) is
replaced by Doench 2014 Rule Set 1 on the same 23-bp context. That swap
preserves the model's input contract -- `cas9_score` remains a 0-100
on-target ranking -- without depending on an external service.

Design contract -- given a plasmid sequence + an edit (substitution /
insertion / deletion) at a 1-based position range:

1. Scan +/- gRNA_search_space (default 200 bp) around the edit for NGG PAMs.
   Valid sgRNAs are those whose nick site is upstream of the edit by
   1..30 nt (auto-grown from 1..10).
2. For each sgRNA, sweep PBS length 10..15 and RTT length 10..20
   (grown to 50 if no candidates). The RTT contains the edited
   sequence with >= 5 nt downstream homology and must not start with C.
3. For PE3, scan opposite-strand sgRNAs. PE3b candidates (ngRNA spacer
   overlapping the edited bases on the opposite strand; pure-substitution
   edits only) are accepted at ANY distance — they self-shutoff via edit-
   induced mismatch and don't need the canonical 32-100 nt band. Non-PE3b
   PE3 ngRNAs are required to sit 32-100 nt from the pegRNA nick (auto-
   grown up to 200 nt if no canonical-band candidates pass).
4. Build the 23-feature row for each (sgRNA, PBS_len, RTT_len, ngRNA)
   combination -- including 10 RNAplfold base-pair-probability features
   on (scaffold + RTT + 14 nt PBS) -- and score with the PE3 XGBoost.
5. Re-rank using easy_prime's force_recommend_dPAM_PE3b rule: dPAM rows
   win if present; PE3b rows stay top if within 10% of best.

Returns the top 3 pegRNAs with full component sequences (spacer,
scaffold, RTT, PBS, ngRNA) and the assembled full pegRNA in the
5prime->3prime order spacer-scaffold-RTT-PBS, plus the XGBoost score
breakdown.
"""
from __future__ import annotations

import os
import pickle
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .doench_score import doench_score_30mer

SCAFFOLD = (
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTT"
    "GAAAAAGTGGCACCGAGTCGGTGC"
)
SCAFFOLD_LEN = len(SCAFFOLD)

DEFAULTS = {
    "min_PBS_length": 10,
    "max_PBS_length": 15,
    "min_RTT_length": 10,
    "max_RTT_length": 20,
    "max_max_RTT_length": 50,
    "min_distance_RTT5": 5,
    "gRNA_search_space": 200,
    "sgRNA_length": 20,
    "offset": -3,
    "PAM": "NGG",
    "max_target_to_sgRNA": 10,
    "max_max_target_to_sgRNA": 30,
    "min_ngRNA_distance": 32,
    "max_ngRNA_distance": 100,
    "max_max_ngRNA_distance": 200,
}

ATTACHED_PBS_LEN = 14

_PE3_FEATURES = [
    "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
    "cas9_score", "nick_to_pegRNA", "dPAM", "PE3b",
    "RTT_GC", "RTT_length", "PBS_GC", "PBS_length",
    "N_subsitution", "N_deletion", "N_insertions",
    "Target_pos", "Target_end_flank",
]

_MODEL_PATH = os.path.join(
    os.path.dirname(__file__), "pegrna_models", "PE3_model_final.pkl"
)
_booster_cache = None


def _load_booster():
    global _booster_cache
    if _booster_cache is None:
        import pickle, xgboost as xgb
        with open(_MODEL_PATH, "rb") as f:
            _m = pickle.load(f)
        _booster_cache = _m.get_booster()
    return _booster_cache


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    s = seq.upper()
    return sum(1 for c in s if c in "GC") / len(s)


def _clean(seq: str) -> str:
    return re.sub(r"[^ACGTNacgtn]", "", seq or "").upper()


def _rnaplfold_features(scaffold: str, rtt: str, attached_pbs: str) -> List[float]:
    try:
        import RNA  # type: ignore
    except ImportError:
        return [0.0] * 10

    seq = (scaffold + rtt + attached_pbs).upper().replace("T", "U")
    md = RNA.md()
    md.max_bp_span = 70
    md.window_size = 70
    fc = RNA.fold_compound(seq, md, RNA.OPTION_WINDOW)

    rtt_pos_start = SCAFFOLD_LEN + 1
    rtt_pos_end = SCAFFOLD_LEN + 10
    scaffold_range = range(1, SCAFFOLD_LEN + 1)
    best = [0.0] * 10

    def _cb(v, vsize, i, maxsize, what, data):
        if not (what & RNA.PROBS_WINDOW_BPP):
            return
        for j in range(1, vsize):
            try:
                p = v[j]
            except Exception:
                continue
            if p is None or p < 0.01:
                continue
            for a, b in ((i, j), (j, i)):
                if rtt_pos_start <= a <= rtt_pos_end and b in scaffold_range:
                    idx = a - rtt_pos_start
                    if p > best[idx]:
                        best[idx] = float(p)

    try:
        fc.probs_window(31, RNA.PROBS_WINDOW_BPP | RNA.PROBS_WINDOW_UP, _cb, None)
    except Exception:
        return [0.0] * 10
    return best


@dataclass
class _Sgrna:
    spacer: str
    pam: str
    start: int
    end: int
    strand: str
    cut_fwd: int
    cas9_score: float


def _doench_for_proto(sequence: str, proto_start: int, proto_end: int,
                      strand: str) -> float:
    L = len(sequence)
    if strand == "+":
        cs, ce = proto_start - 4, proto_end + 6
        if cs < 0 or ce > L:
            return 0.0
        ctx = sequence[cs:ce]
    else:
        cs, ce = proto_start - 6, proto_end + 4
        if cs < 0 or ce > L:
            return 0.0
        ctx = revcomp(sequence[cs:ce])
    if len(ctx) != 30:
        return 0.0
    try:
        d, _ = doench_score_30mer(ctx)
        return float(d)
    except Exception:
        return 0.0


def _scan_sgrnas(sequence: str, search_start: int, search_end: int) -> List[_Sgrna]:
    out: List[_Sgrna] = []
    L = len(sequence)
    fwd = sequence
    rc = revcomp(sequence)

    for i in range(20, L - 2):
        proto_start = i - 20
        proto_end = i
        if proto_end < search_start or proto_start > search_end:
            continue
        pam = fwd[i:i + 3]
        if len(pam) < 3 or pam[1] != "G" or pam[2] != "G":
            continue
        spacer = fwd[proto_start:proto_end]
        if "N" in spacer:
            continue
        cas9 = _doench_for_proto(fwd, proto_start, proto_end, "+")
        cut = proto_end - 3
        out.append(_Sgrna(spacer, pam, proto_start, proto_end, "+", cut, cas9))

    for i in range(20, L - 2):
        proto_start_rc = i - 20
        proto_end_rc = i
        pam = rc[i:i + 3]
        if len(pam) < 3 or pam[1] != "G" or pam[2] != "G":
            continue
        spacer = rc[proto_start_rc:proto_end_rc]
        if "N" in spacer:
            continue
        fwd_end = L - proto_start_rc
        fwd_start = L - proto_end_rc
        if fwd_end < search_start or fwd_start > search_end:
            continue
        cas9 = _doench_for_proto(fwd, fwd_start, fwd_end, "-")
        cut = fwd_start + 3
        out.append(_Sgrna(spacer, pam, fwd_start, fwd_end, "-", cut, cas9))

    return out


def _apply_edit(seq: str, pos0: int, ref: str, alt: str) -> str:
    return seq[:pos0] + alt + seq[pos0 + len(ref):]


@dataclass
class _Candidate:
    sgrna: _Sgrna
    pbs_seq: str
    pbs_len: int
    rtt_seq: str
    rtt_len: int
    target_pos_in_rtt: int
    target_end_flank: int
    is_dpam: int
    rnaplfold: List[float]
    ngrna: Optional[_Sgrna] = None
    nick_to_peg: Optional[int] = None
    is_pe3b: int = 0
    ngrna_spacer_edited: Optional[str] = None


def _enumerate_pbs(sgrna: _Sgrna, sequence: str, params: Dict) -> List[Tuple[int, str]]:
    out: List[Tuple[int, str]] = []
    L = len(sequence)
    for plen in range(params["min_PBS_length"], params["max_PBS_length"] + 1):
        if sgrna.strand == "+":
            s, e = sgrna.cut_fwd - plen, sgrna.cut_fwd
            if s < 0 or e > L:
                continue
            pbs = revcomp(sequence[s:e])
        else:
            s, e = sgrna.cut_fwd, sgrna.cut_fwd + plen
            if s < 0 or e > L:
                continue
            pbs = sequence[s:e]
        out.append((plen, pbs))
    return out


def _enumerate_rtt(sgrna, sequence, edit_pos0, ref, alt, params):
    L = len(sequence)
    del_len = max(0, len(ref) - len(alt))
    min_d5 = params["min_distance_RTT5"]
    rtt_cap = params["max_RTT_length"] + max(0, del_len)
    max_cap = params["max_max_RTT_length"] + max(0, del_len)
    out: List[Tuple[int, str, int, int]] = []

    while rtt_cap <= max_cap:
        candidates: List[Tuple[int, str, int, int]] = []
        for rlen in range(params["min_RTT_length"], rtt_cap + 1):
            if sgrna.strand == "+":
                gs, ge = sgrna.cut_fwd, sgrna.cut_fwd + rlen
                if ge > L:
                    continue
                upper = sgrna.cut_fwd + rlen - min_d5 - max(1, len(ref))
                if not (sgrna.cut_fwd <= edit_pos0 <= upper):
                    continue
                template = sequence[gs:ge]
                rel = edit_pos0 - gs
                edited = _apply_edit(template, rel, ref, alt)
                rtt = revcomp(edited)
                edit_end_in_edited = rel + max(1, len(alt))
                tp = len(edited) - edit_end_in_edited + 1
                ef = rel
            else:
                gs, ge = sgrna.cut_fwd - rlen, sgrna.cut_fwd
                if gs < 0:
                    continue
                lower = sgrna.cut_fwd - rlen + min_d5
                upper = sgrna.cut_fwd - max(1, len(ref))
                if not (lower <= edit_pos0 <= upper):
                    continue
                template = sequence[gs:ge]
                rel = edit_pos0 - gs
                edited = _apply_edit(template, rel, ref, alt)
                rtt = edited
                tp = rel + 1
                ef = (sgrna.cut_fwd - 1) - (edit_pos0 + max(0, len(ref) - 1))
            if not rtt or rtt[0] == "C":
                continue
            candidates.append((rlen, rtt, tp, ef))
        if candidates:
            out = candidates
            break
        rtt_cap += 5
    return out


def _is_dpam(rtt: str) -> int:
    rev = revcomp(rtt)
    if len(rev) < 6:
        return 0
    candidate = rev[3:6]
    return 0 if (candidate[1:] == "GG") else 1


def _attached_pbs(sgrna: _Sgrna, sequence: str) -> str:
    L = len(sequence)
    n = ATTACHED_PBS_LEN
    if sgrna.strand == "+":
        s, e = max(0, sgrna.cut_fwd - n), sgrna.cut_fwd
        return revcomp(sequence[s:e])
    else:
        s, e = sgrna.cut_fwd, min(L, sgrna.cut_fwd + n)
        return sequence[s:e]


def _select_ngrna(peg, all_sgrnas, edit_pos0, ref, alt, params):
    opposite = "-" if peg.strand == "+" else "+"
    max_dist = params["max_ngRNA_distance"]
    max_max = params["max_max_ngRNA_distance"]

    while max_dist <= max_max:
        candidates: List[Tuple[_Sgrna, int, int, Optional[str]]] = []
        for s in all_sgrnas:
            if s.strand != opposite:
                continue
            dist = s.cut_fwd - peg.cut_fwd
            if peg.strand == "-":
                dist = -dist
            # Detect PE3b FIRST (Anzalone 2019, easy_prime
            # force_recommend_dPAM_PE3b re-rank): an ngRNA whose spacer
            # overlaps the edited bases on the opposite strand stops
            # cutting once the edit is installed (the spacer no longer
            # matches), avoiding the simultaneous-double-nick risk that
            # the 32-100 bp distance band exists to mitigate. PE3b
            # candidates therefore BYPASS the distance filter — they
            # typically sit <32 bp from the pegRNA nick.
            is_pe3b = 0
            edited_spacer: Optional[str] = None
            if len(ref) == len(alt) and len(ref) > 0:
                edit_window = range(edit_pos0, edit_pos0 + len(ref))
                if any(s.start <= p < s.end for p in edit_window):
                    is_pe3b = 1
                    if s.strand == "+":
                        rel = edit_pos0 - s.start
                        if 0 <= rel and rel + len(ref) <= len(s.spacer):
                            edited_spacer = (
                                s.spacer[:rel] + alt + s.spacer[rel + len(ref):]
                            )
                    else:
                        rev_window = revcomp(alt)
                        rel = (s.end - (edit_pos0 + len(ref)))
                        if 0 <= rel and rel + len(ref) <= len(s.spacer):
                            edited_spacer = (
                                s.spacer[:rel] + rev_window + s.spacer[rel + len(ref):]
                            )
            # Distance filter applies ONLY to non-PE3b candidates. PE3
            # ngRNAs depend on heteroduplex resolution before the second
            # nick lands, which the 32-100 bp band optimises for; PE3b
            # ngRNAs depend on edit-induced mismatch self-shutoff and
            # are insensitive to distance.
            if not is_pe3b:
                if abs(dist) > max_dist or abs(dist) < params["min_ngRNA_distance"]:
                    continue
            candidates.append((s, dist, is_pe3b, edited_spacer))
        if candidates:
            return candidates
        max_dist += 20
    return []


def design_pegrnas(
    sequence: str,
    edit_start_1based: int,
    edit_end_1based: int,
    alt: str,
    edit_type: str = "substitution",
    n_results: int = 3,
    use_pe3: bool = True,
    params: Optional[Dict] = None,
) -> Dict:
    p = dict(DEFAULTS)
    if params:
        p.update(params)

    seq = _clean(sequence)
    L = len(seq)
    if L == 0:
        return {"ok": False, "error": "Empty plasmid sequence", "pegrnas": []}

    s1 = max(1, int(edit_start_1based))
    e1 = max(s1, int(edit_end_1based))
    if e1 > L:
        return {"ok": False, "error": f"Edit end {e1} > plasmid length {L}", "pegrnas": []}
    edit_pos0 = s1 - 1
    edit_len = e1 - s1 + 1
    ref = seq[edit_pos0:edit_pos0 + edit_len]
    alt_clean = _clean(alt) if alt else ""

    if edit_type == "deletion":
        alt_clean = ""
    elif edit_type == "insertion":
        ref = ""
        edit_len = 0
        if not alt_clean:
            return {"ok": False, "error": "Insertion requires non-empty alt", "pegrnas": []}
    elif edit_type == "substitution":
        if len(alt_clean) != len(ref):
            return {"ok": False, "error": "Substitution requires len(alt) == len(ref)", "pegrnas": []}

    span = p["gRNA_search_space"]
    ss = max(0, edit_pos0 - span)
    se = min(L, edit_pos0 + max(edit_len, 1) + span)
    all_sgrnas = _scan_sgrnas(seq, ss, se)

    valid_pegs: List[Tuple[_Sgrna, int]] = []
    max_d = p["max_target_to_sgRNA"]
    max_max_d = p["max_max_target_to_sgRNA"]
    while max_d <= max_max_d:
        for s in all_sgrnas:
            if s.strand == "+":
                td = edit_pos0 - s.cut_fwd
            else:
                td = (s.cut_fwd - 1) - (edit_pos0 + max(edit_len, 1) - 1)
            if 1 <= td <= max_d:
                valid_pegs.append((s, td))
        if valid_pegs:
            break
        max_d += 5

    if not valid_pegs:
        return {"ok": False, "error": "No valid pegRNA spacer found within search window",
                "pegrnas": [], "summary": {"n_sgrnas_scanned": len(all_sgrnas), "n_valid_pegRNAs": 0}}

    rows: List[Dict] = []
    candidates: List[_Candidate] = []
    n_subst = (sum(1 for a, b in zip(ref, alt_clean) if a != b)
               if edit_type == "substitution" else 0)
    n_ins = len(alt_clean) if edit_type == "insertion" else 0
    n_del = len(ref) if edit_type == "deletion" else 0

    for peg, _td in valid_pegs:
        pbs_options = _enumerate_pbs(peg, seq, p)
        rtt_options = _enumerate_rtt(peg, seq, edit_pos0, ref, alt_clean, p)
        if not pbs_options or not rtt_options:
            continue
        attached = _attached_pbs(peg, seq)
        if use_pe3:
            ng_opts = _select_ngrna(peg, all_sgrnas, edit_pos0, ref, alt_clean, p)
            if not ng_opts:
                continue
        else:
            ng_opts = [(None, 0, 0, None)]

        for plen, pbs in pbs_options:
            for rlen, rtt, tp, ef in rtt_options:
                dpam = _is_dpam(rtt)
                fold10 = _rnaplfold_features(SCAFFOLD, rtt, attached)
                for ng, dist, is_pe3b, edited_ngsp in ng_opts:
                    cand = _Candidate(
                        sgrna=peg, pbs_seq=pbs, pbs_len=plen,
                        rtt_seq=rtt, rtt_len=rlen,
                        target_pos_in_rtt=tp, target_end_flank=ef,
                        is_dpam=dpam, rnaplfold=fold10,
                        ngrna=ng, nick_to_peg=dist, is_pe3b=is_pe3b,
                        ngrna_spacer_edited=edited_ngsp,
                    )
                    candidates.append(cand)
                    feats = {f"{i}": fold10[i] for i in range(10)}
                    feats.update({
                        "cas9_score": peg.cas9_score,
                        "nick_to_pegRNA": float(dist if dist is not None else 0),
                        "dPAM": float(dpam),
                        "PE3b": float(is_pe3b),
                        "RTT_GC": gc_fraction(rtt),
                        "RTT_length": float(rlen),
                        "PBS_GC": gc_fraction(pbs),
                        "PBS_length": float(plen),
                        "N_subsitution": float(n_subst),
                        "N_deletion": float(n_del),
                        "N_insertions": float(n_ins),
                        "Target_pos": float(tp),
                        "Target_end_flank": float(ef),
                    })
                    rows.append(feats)

    if not rows:
        return {"ok": False, "error": "No valid pegRNA candidates after PBS/RTT/ngRNA sweep",
                "pegrnas": [], "summary": {
                    "n_sgrnas_scanned": len(all_sgrnas),
                    "n_valid_pegRNAs": len(valid_pegs),
                    "n_candidates": 0,
                }}

    import xgboost as _xgb
    X = pd.DataFrame(rows)[_PE3_FEATURES]
    _booster = _load_booster()
    preds = _booster.predict(_xgb.DMatrix(X.values))

    has_dpam = any(c.is_dpam for c in candidates)
    keep_idx = [i for i, c in enumerate(candidates) if (c.is_dpam or not has_dpam)]
    pool = [(i, candidates[i], float(preds[i])) for i in keep_idx]
    max_eff = max((eff for _, _, eff in pool), default=0.0)

    def _rank_key(item):
        _i, c, eff = item
        rank = 0
        if c.is_dpam:
            rank += 1
        if c.is_pe3b and max_eff > 0 and (max_eff - eff) <= 0.1 * max_eff:
            rank += 1
        return (-rank, -eff)

    pool.sort(key=_rank_key)
    top = pool[:max(1, int(n_results))]

    out_list: List[Dict] = []
    for rank_i, (_i, c, eff) in enumerate(top, start=1):
        full_pegrna = c.sgrna.spacer + SCAFFOLD + c.rtt_seq + c.pbs_seq
        ng_dict = None
        if c.ngrna is not None:
            ng_dict = {
                "spacer": c.ngrna_spacer_edited or c.ngrna.spacer,
                "original_spacer": c.ngrna.spacer,
                "pam": c.ngrna.pam,
                "start": c.ngrna.start,
                "end": c.ngrna.end,
                "strand": c.ngrna.strand,
                "cut_fwd": c.ngrna.cut_fwd,
                "cas9_score": round(c.ngrna.cas9_score, 2),
                "nick_to_pegRNA": c.nick_to_peg,
                "is_pe3b": bool(c.is_pe3b),
            }
        out_list.append({
            "rank": rank_i,
            "name": f"pegRNA_{rank_i}_{c.sgrna.spacer[:6]}",
            "predicted_efficiency": round(eff, 3),
            "spacer": c.sgrna.spacer,
            "pam": c.sgrna.pam,
            "spacer_start": c.sgrna.start,
            "spacer_end": c.sgrna.end,
            "direction": 1 if c.sgrna.strand == "+" else -1,
            "strand": c.sgrna.strand,
            "cut_fwd": c.sgrna.cut_fwd,
            "cas9_score": round(c.sgrna.cas9_score, 2),
            "rtt": c.rtt_seq,
            "rtt_length": c.rtt_len,
            "rtt_gc": round(gc_fraction(c.rtt_seq), 3),
            "pbs": c.pbs_seq,
            "pbs_length": c.pbs_len,
            "pbs_gc": round(gc_fraction(c.pbs_seq), 3),
            "scaffold": SCAFFOLD,
            "full_pegrna": full_pegrna,
            "full_pegrna_length": len(full_pegrna),
            "is_dpam": bool(c.is_dpam),
            "is_pe3b": bool(c.is_pe3b),
            "ngrna": ng_dict,
            "edit_type": edit_type,
            "edit_ref": ref,
            "edit_alt": alt_clean,
            "edit_start_1based": s1,
            "edit_end_1based": e1,
            "score_components": {
                "cas9_score": round(c.sgrna.cas9_score, 2),
                "rtt_gc": round(gc_fraction(c.rtt_seq), 3),
                "rtt_length": c.rtt_len,
                "pbs_gc": round(gc_fraction(c.pbs_seq), 3),
                "pbs_length": c.pbs_len,
                "dpam": bool(c.is_dpam),
                "pe3b": bool(c.is_pe3b),
                "nick_to_pegRNA": c.nick_to_peg,
                "rnaplfold_max": round(max(c.rnaplfold), 4),
                "rnaplfold_mean": round(sum(c.rnaplfold) / 10.0, 4),
                "target_pos_in_rtt": c.target_pos_in_rtt,
                "target_end_flank": c.target_end_flank,
            },
        })

    return {
        "ok": True,
        "pegrnas": out_list,
        "summary": {
            "n_sgrnas_scanned": len(all_sgrnas),
            "n_valid_pegRNAs": len(valid_pegs),
            "n_candidates": len(candidates),
            "n_returned": len(out_list),
            "use_pe3": use_pe3,
            "edit_type": edit_type,
            "edit_ref": ref,
            "edit_alt": alt_clean,
            "edit_start_1based": s1,
            "edit_end_1based": e1,
            "model": "easy_prime PE3 XGBoost (Li et al. 2021) + Doench 2014 cas9_score",
            "scaffold_length": SCAFFOLD_LEN,
        },
    }
