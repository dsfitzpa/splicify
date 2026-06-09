"""emit_guides_csv — 5th output emitter, for the CRISPR pipeline.

One CSV that captures every designed guide (sgRNA / pegRNA / ngRNA),
every editing-readout primer (Sanger + Illumina), and (optionally) the
cloning oligo pair the agent computed for each sgRNA. Schema is wide-
because-flat: every row has the same column set, with type-specific
columns populated only on the relevant rows.

The agent calls this AFTER design_guides / design_pegrnas / design_primers
have run. Pass it the digested entries directly; the emitter does not
re-run any design. Cloning oligos are passed in by the agent (computed
to match the user's destination vector — pX330/BbsI, lentiCRISPR/BsmBI,
etc.); the emitter doesn't synthesise them.
"""
from __future__ import annotations

import base64
import csv
import io
import pathlib
from typing import Any, Optional

import agent_v2  # noqa: F401 - triggers path shim
from agent_v2.outputs import prefixed_filename, derive_descriptor


# Column order — wide-but-flat. type-specific fields are blank on irrelevant rows.
_COLUMNS = [
    "type",                  # sgRNA | pegRNA | ngRNA | primer_sanger_fwd | primer_sanger_rev | primer_illumina_fwd | primer_illumina_rev | cloning_oligo
    "name",
    "target_attachment_id",
    "spacer",
    "pam",
    "start",
    "end",
    "direction",
    "score",
    "score_method",
    "gc_fraction",
    "n_offtargets",
    "rtt",
    "pbs",
    "scaffold",
    "full_pegrna",
    "predicted_efficiency",
    "is_dpam",
    "is_pe3b",
    "edit_type",
    "edit_ref",
    "edit_alt",
    "primer_sequence",
    "annealing_sequence",
    "adapter",
    "tm",
    "product_size",
    "application",
    "region_start_plasmid",
    "region_end_plasmid",
    "cloning_oligo_top",
    "cloning_oligo_bottom",
    "notes",
]


def _row(**kwargs) -> dict[str, Any]:
    """Build a row with every column key present (blanks where unset)."""
    return {c: kwargs.get(c, "") for c in _COLUMNS}


async def emit_guides_csv(
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,
) -> dict[str, Any]:
    guides = args.get("guides") or []
    pegrnas = args.get("pegrnas") or []
    primers = args.get("primers") or []
    cloning_oligos = args.get("cloning_oligos") or []

    rows: list[dict[str, Any]] = []

    for g in guides:
        rows.append(_row(
            type="sgRNA",
            name=g.get("name", ""),
            target_attachment_id=g.get("target_attachment_id", ""),
            spacer=g.get("spacer", ""),
            pam=g.get("pam", ""),
            start=g.get("start", ""),
            end=g.get("end", ""),
            direction=g.get("direction", ""),
            score=g.get("score", ""),
            score_method=g.get("score_method", ""),
            gc_fraction=g.get("gc_fraction", ""),
            n_offtargets=g.get("n_offtargets", ""),
            notes=g.get("notes", ""),
        ))

    for p in pegrnas:
        rows.append(_row(
            type="pegRNA",
            name=p.get("name", ""),
            target_attachment_id=p.get("target_attachment_id", ""),
            spacer=p.get("spacer", ""),
            pam=p.get("pam", ""),
            start=p.get("spacer_start", ""),
            end=p.get("spacer_end", ""),
            direction=p.get("direction", ""),
            score=p.get("cas9_score", ""),
            score_method="doench2014",
            rtt=p.get("rtt", ""),
            pbs=p.get("pbs", ""),
            scaffold=p.get("scaffold", ""),
            full_pegrna=p.get("full_pegrna", ""),
            predicted_efficiency=p.get("predicted_efficiency", ""),
            is_dpam=p.get("is_dpam", ""),
            is_pe3b=p.get("is_pe3b", ""),
            edit_type=p.get("edit_type", ""),
            edit_ref=p.get("edit_ref", ""),
            edit_alt=p.get("edit_alt", ""),
            notes=p.get("notes", ""),
        ))
        ngrna = p.get("ngrna")
        if ngrna:
            rows.append(_row(
                type="ngRNA",
                name=f"ngRNA_for_{p.get('name', '')}",
                target_attachment_id=p.get("target_attachment_id", ""),
                spacer=ngrna.get("spacer", ""),
                pam=ngrna.get("pam", ""),
                start=ngrna.get("start", ""),
                end=ngrna.get("end", ""),
                direction="-1" if ngrna.get("strand") == "-" else "1",
                score=ngrna.get("cas9_score", ""),
                score_method="doench2014",
                is_pe3b=ngrna.get("is_pe3b", ""),
                notes=f"paired with {p.get('name', '')}; nick_to_pegRNA={ngrna.get('nick_to_pegRNA', '')}",
            ))

    for pr in primers:
        application = pr.get("application", "")
        # Forward
        rows.append(_row(
            type=f"primer_{application}_fwd" if application else "primer_fwd",
            name=pr.get("name_fwd") or f"{pr.get('pair_label', 'primer')}_F",
            target_attachment_id=pr.get("target_attachment_id", ""),
            primer_sequence=pr.get("left_primer", ""),
            annealing_sequence=pr.get("left_annealing", pr.get("left_primer", "")),
            adapter=pr.get("left_adapter", ""),
            tm=pr.get("left_tm", ""),
            product_size=pr.get("product_size", ""),
            application=application,
            region_start_plasmid=pr.get("region_start_plasmid", ""),
            region_end_plasmid=pr.get("region_end_plasmid", ""),
            start=pr.get("left_pos_plasmid", ""),
            notes=pr.get("notes", ""),
        ))
        # Reverse
        rows.append(_row(
            type=f"primer_{application}_rev" if application else "primer_rev",
            name=pr.get("name_rev") or f"{pr.get('pair_label', 'primer')}_R",
            target_attachment_id=pr.get("target_attachment_id", ""),
            primer_sequence=pr.get("right_primer", ""),
            annealing_sequence=pr.get("right_annealing", pr.get("right_primer", "")),
            adapter=pr.get("right_adapter", ""),
            tm=pr.get("right_tm", ""),
            product_size=pr.get("product_size", ""),
            application=application,
            region_start_plasmid=pr.get("region_start_plasmid", ""),
            region_end_plasmid=pr.get("region_end_plasmid", ""),
            end=pr.get("right_pos_plasmid", ""),
            notes=pr.get("notes", ""),
        ))

    for oligo in cloning_oligos:
        rows.append(_row(
            type="cloning_oligo",
            name=oligo.get("name", ""),
            spacer=oligo.get("spacer", ""),
            cloning_oligo_top=oligo.get("oligo_top", ""),
            cloning_oligo_bottom=oligo.get("oligo_bottom", ""),
            notes=oligo.get("notes", "sgRNA Golden Gate cloning oligo pair"),
        ))

    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=_COLUMNS, quoting=csv.QUOTE_MINIMAL)
    writer.writeheader()
    for r in rows:
        writer.writerow(r)
    csv_text = buf.getvalue()

    _descriptor = derive_descriptor(args)
    file_envelope = {
        "fileName": prefixed_filename("guides.csv", _descriptor),
        "dataBase64": base64.b64encode(csv_text.encode("utf-8")).decode("ascii"),
    }
    written_path: Optional[str] = None
    if output_dir is not None:
        out = pathlib.Path(output_dir) / prefixed_filename("guides.csv", _descriptor)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(csv_text)
        written_path = str(out)

    return {
        "ok": True,
        "file": file_envelope,
        "n_rows": len(rows),
        "n_sgRNAs": len(guides),
        "n_pegRNAs": len(pegrnas),
        "n_primer_pairs": len(primers),
        "n_cloning_oligos": len(cloning_oligos),
        "written_path": written_path,
    }
