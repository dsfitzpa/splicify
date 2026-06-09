"""emit_guides_gb — 6th output emitter, for the CRISPR pipeline.

Takes the source GenBank for a target attachment, appends SeqFeature
entries for every designed sgRNA / pegRNA / ngRNA / primer, and emits a
new GenBank file. The frontend's CircularPlasmidViewer can render this
directly via /agent_v2/annotate-on-upload + the standard viz pipeline.

If no source GenBank text was stashed at upload time (rare — the v2
router stashes it by default since iter 38), fall back to constructing a
minimal GenBank record from the registered sequence.
"""
from __future__ import annotations

import base64
import io
import pathlib
from typing import Any, Optional

import agent_v2  # noqa: F401 - triggers path shim
from agent_v2.outputs import prefixed_filename, derive_descriptor


_ADDED_BY_TAG = "crispr_v2"


def _make_feature(ftype: str, start: int, end: int, strand: int,
                   qualifiers: dict[str, list[str]]):
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    # Tag every feature this emitter adds so the frontend can color them
    # distinctly from pre-existing annotations on the source .gb.
    qualifiers = dict(qualifiers)
    qualifiers.setdefault("added_by", [_ADDED_BY_TAG])
    return SeqFeature(
        location=SimpleLocation(int(start), int(end), strand=int(strand) or 1),
        type=ftype,
        qualifiers=qualifiers,
    )


async def emit_guides_gb(
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,
) -> dict[str, Any]:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from agent_v2 import attachment_kinds

    import time as __t, sys as __sys
    __t0 = __t.monotonic()
    print(f"[emit_guides_gb] start target_attachment_id={args.get('target_attachment_id')} pegrnas={len(args.get('pegrnas') or [])} guides={len(args.get('guides') or [])} primers={len(args.get('primers') or [])}", file=__sys.stderr, flush=True)
    target_aid = args.get("target_attachment_id")
    att = registry.get(target_aid) if target_aid else None
    if att is None:
        return {"ok": False, "error": f"unknown target_attachment_id: {target_aid!r}"}

    # Prefer the source GenBank text the router stashed at upload time.
    cached = None
    fk = attachment_kinds.get_kind(target_aid)
    if fk is not None:
        cached_obj = None
        with attachment_kinds._LOCK:  # private but stable internal API
            cached_obj = attachment_kinds._CACHE.get(target_aid)
        if cached_obj is not None and cached_obj.gb_text:
            cached = cached_obj.gb_text

    if cached:
        try:
            record = next(SeqIO.parse(io.StringIO(cached), "genbank"))
        except Exception:
            cached = None

    if not cached:
        # Build a minimal record from the bare sequence in the registry.
        topology = "circular" if att.circular else "linear"
        record = SeqRecord(
            Seq(att.sequence),
            id=att.name or "attachment",
            name=(att.name or "attachment").replace(" ", "_")[:16],
            description=att.name or "attachment",
            annotations={"molecule_type": "DNA", "topology": topology},
        )

    # Append features.
    added = 0
    for g in args.get("guides") or []:
        qual = {
            "label": [g.get("name", "sgRNA")],
            "note": [
                f"sgRNA spacer={g.get('spacer','')} pam={g.get('pam','')} "
                f"score={g.get('score','')} method={g.get('score_method','')} "
                f"gc={g.get('gc_fraction','')} n_off={g.get('n_offtargets','')}"
            ],
        }
        record.features.append(_make_feature(
            "misc_RNA", g.get("start", 0), g.get("end", 0),
            g.get("direction", 1), qual,
        ))
        added += 1

    for p in args.get("pegrnas") or []:
        qual = {
            "label": [p.get("name", "pegRNA")],
            "note": [
                f"pegRNA spacer={p.get('spacer','')} pam={p.get('pam','')} "
                f"rank={p.get('rank','')} eff={p.get('predicted_efficiency','')} "
                f"rtt={p.get('rtt','')} pbs={p.get('pbs','')} "
                f"is_dpam={p.get('is_dpam','')} is_pe3b={p.get('is_pe3b','')} "
                f"edit={p.get('edit_type','')} {p.get('edit_ref','')}>{p.get('edit_alt','')}"
            ],
        }
        record.features.append(_make_feature(
            "misc_RNA", p.get("spacer_start", 0), p.get("spacer_end", 0),
            p.get("direction", 1), qual,
        ))
        added += 1
        ng = p.get("ngrna") or None
        if ng and ng.get("start") is not None and ng.get("end") is not None:
            qng = {
                "label": [f"ngRNA_for_{p.get('name','peg')}"],
                "note": [
                    f"ngRNA spacer={ng.get('spacer','')} pam={ng.get('pam','')} "
                    f"score={ng.get('cas9_score','')} "
                    f"nick_to_pegRNA={ng.get('nick_to_pegRNA','')} "
                    f"is_pe3b={ng.get('is_pe3b','')}"
                ],
            }
            record.features.append(_make_feature(
                "misc_RNA", ng.get("start"), ng.get("end"),
                -1 if ng.get("strand") == "-" else 1, qng,
            ))
            added += 1

    for pr in args.get("primers") or []:
        application = pr.get("application", "")
        if pr.get("left_pos_plasmid") is not None and pr.get("left_annealing"):
            anneal_len = len(pr["left_annealing"])
            qf = {
                "label": [pr.get("name_fwd") or f"{pr.get('pair_label','primer')}_F"],
                "note": [
                    f"primer({application}) tm={pr.get('left_tm','')} "
                    f"product_size={pr.get('product_size','')} "
                    f"adapter={pr.get('left_adapter','')}"
                ],
            }
            lp = int(pr["left_pos_plasmid"])
            record.features.append(_make_feature("primer_bind", lp, lp + anneal_len, +1, qf))
            added += 1
        if pr.get("right_pos_plasmid") is not None and pr.get("right_annealing"):
            anneal_len = len(pr["right_annealing"])
            qr = {
                "label": [pr.get("name_rev") or f"{pr.get('pair_label','primer')}_R"],
                "note": [
                    f"primer({application}) tm={pr.get('right_tm','')} "
                    f"product_size={pr.get('product_size','')} "
                    f"adapter={pr.get('right_adapter','')}"
                ],
            }
            rp = int(pr["right_pos_plasmid"])
            # right_pos is the 3' end of the reverse primer on the + strand;
            # the annealing region runs upstream on the - strand from rp.
            record.features.append(_make_feature(
                "primer_bind", max(0, rp - anneal_len), rp, -1, qr,
            ))
            added += 1

    # Ensure record.annotations has molecule_type set (Biopython refuses to
    # write GenBank without it).
    if "molecule_type" not in record.annotations:
        record.annotations["molecule_type"] = "DNA"

    out_buf = io.StringIO()
    SeqIO.write([record], out_buf, "genbank")
    gb_text = out_buf.getvalue()

    _descriptor = derive_descriptor(args)
    file_envelope = {
        "fileName": prefixed_filename("guides.gb", _descriptor),
        "dataBase64": base64.b64encode(gb_text.encode("utf-8")).decode("ascii"),
    }
    written_path: Optional[str] = None
    if output_dir is not None:
        out = pathlib.Path(output_dir) / prefixed_filename("guides.gb", _descriptor)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(gb_text)
        written_path = str(out)

    print(f"[emit_guides_gb] done in {int((__t.monotonic() - __t0)*1000)} ms; n_features_added={added}", file=__sys.stderr, flush=True)
    return {
        "ok": True,
        "file": file_envelope,
        "n_features_added": added,
        "length_bp": len(record.seq),
        "topology": record.annotations.get("topology", "linear"),
        "written_path": written_path,
    }
