"""
File generation for AI Plasmid Design results.
Produces CSV, GenBank (.gb), and plain-text protocol files
encoded as base64 for the frontend download buttons.
"""
from __future__ import annotations

import base64
import csv
import io
from typing import Any, Dict, List, Tuple, Optional


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _b64(data: bytes) -> str:
    return base64.b64encode(data).decode("ascii")


def _file(name: str, mime: str, data: bytes) -> Dict[str, str]:
    return {"fileName": name, "mimeType": mime, "dataBase64": _b64(data)}


def _csv_bytes(rows: List[List[Any]], header: List[str]) -> bytes:
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(header)
    w.writerows(rows)
    return buf.getvalue().encode("utf-8")


def _fmt(val: Any, decimals: int = 1) -> str:
    """Format a number or return '' for None."""
    if val is None:
        return ""
    try:
        return f"{float(val):.{decimals}f}"
    except (TypeError, ValueError):
        return str(val)


def _make_genbank(
    seq: str,
    name: str,
    description: str,
    annotations: List[Dict[str, Any]],
    topology: str = "circular",
) -> bytes:
    """
    Build a minimal GenBank-format file using BioPython.
    Falls back to a simple FASTA if BioPython is unavailable.
    """
    try:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio import SeqIO

        record = SeqRecord(
            Seq(seq),
            id=name[:16] or "sequence",
            name=(name[:16] or "sequence"),
            description=description[:80],
        )
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = topology

        for ann in annotations:
            start = int(ann.get("start", 0))
            end = int(ann.get("end", start + 1))
            # Clamp to sequence length
            start = max(0, min(start, len(seq)))
            end = max(start + 1, min(end, len(seq)))
            direction = ann.get("direction", 1)
            strand = -1 if direction == -1 else 1
            feat_type = ann.get("feat_type", "misc_feature")
            qualifiers: Dict[str, Any] = {"label": [ann.get("name", "feature")]}
            if ann.get("sequence"):
                qualifiers["note"] = [f"sequence={ann['sequence']}"]
            if ann.get("tm") is not None:
                qualifiers["note"] = qualifiers.get("note", []) + [f"Tm={ann['tm']:.1f}C"]
            feature = SeqFeature(
                FeatureLocation(start, end, strand=strand),
                type=feat_type,
                qualifiers=qualifiers,
            )
            record.features.append(feature)

        buf = io.StringIO()
        SeqIO.write(record, buf, "genbank")
        return buf.getvalue().encode("utf-8")

    except Exception:
        # Fallback: simple FASTA
        lines = [f">{name} {description}"]
        for i in range(0, len(seq), 60):
            lines.append(seq[i:i + 60])
        return "\n".join(lines).encode("utf-8")


# ---------------------------------------------------------------------------
# Gibson Assembly files
# ---------------------------------------------------------------------------

def build_gibson_files(result: Dict[str, Any]) -> List[Dict[str, str]]:
    """
    Generate downloadable files for a gibson_design result:
    - primers_by_fragment.csv
    - junctions.csv
    - assembled_construct.gb
    - protocol.txt
    """
    files: List[Dict[str, str]] = []

    primers = result.get("primers_by_fragment", [])
    junctions = result.get("junctions", [])
    viz = result.get("viz", {}) or {}
    assembled_seq = viz.get("sequence", "")
    assembly = result.get("assembly", "circular")

    # ── 1. primers_by_fragment.csv ─────────────────────────────────────────
    header = [
        "Fragment", "Needs Primers",
        "Forward Primer (5'→3')", "Reverse Primer (5'→3')",
        "Fwd Anneal Tm (°C)", "Rev Anneal Tm (°C)",
        "Fwd Anneal Seq", "Rev Anneal Seq",
        "Fwd Extension Seq", "Rev Extension Seq",
        "Fwd Total Score", "Rev Total Score",
        "Fwd Uniqueness Score", "Rev Uniqueness Score",
        "Fwd Hairpin dG (kcal/mol)", "Rev Hairpin dG (kcal/mol)",
        "Fwd Warnings", "Rev Warnings",
    ]
    rows = []
    for p in primers:
        frag = p.get("fragment", "")
        needs = p.get("needs_primers", False)
        if needs:
            rows.append([
                frag, "Yes",
                p.get("forward_primer", ""),
                p.get("reverse_primer", ""),
                _fmt(p.get("forward_anneal_tm")),
                _fmt(p.get("reverse_anneal_tm")),
                p.get("forward_anneal_seq", ""),
                p.get("reverse_anneal_seq", ""),
                p.get("forward_extension_seq", ""),
                p.get("reverse_extension_seq", ""),
                _fmt(p.get("forward_extension_total_score"), 1),
                _fmt(p.get("reverse_extension_total_score"), 1),
                _fmt(p.get("forward_extension_uniqueness_score"), 1),
                _fmt(p.get("reverse_extension_uniqueness_score"), 1),
                _fmt(p.get("forward_extension_hairpin_dg"), 2),
                _fmt(p.get("reverse_extension_hairpin_dg"), 2),
                p.get("forward_extension_warnings", ""),
                p.get("reverse_extension_warnings", ""),
            ])
        else:
            rows.append([frag, f"No ({p.get('reason', '')})", *[""] * 16])

    files.append(_file(
        "primers_by_fragment.csv", "text/csv",
        _csv_bytes(rows, header),
    ))

    # ── 2. junctions.csv ──────────────────────────────────────────────────
    j_header = [
        "From Fragment", "To Fragment", "Source",
        "Overlap Sequence", "Length (bp)", "Tm (°C)", "Score (0-100)",
        "Left bp", "Right bp",
        "GC Content (%)", "Uniqueness Score",
        "Hairpin dG (kcal/mol)", "Self-Dimer dG (kcal/mol)",
        "Off-Target Count", "Warnings",
    ]
    j_rows = []
    for j in junctions:
        ov_seq = j.get("overlap_sequence", "")
        gc = (
            round((ov_seq.upper().count("G") + ov_seq.upper().count("C")) / len(ov_seq) * 100)
            if ov_seq else ""
        )
        j_rows.append([
            j.get("from", ""), j.get("to", ""), j.get("source", ""),
            ov_seq,
            j.get("overlap_length", ""),
            _fmt(j.get("overlap_tm")),
            _fmt(j.get("overlap_score")),
            j.get("left_bp", ""), j.get("right_bp", ""),
            gc,
            _fmt(j.get("overlap_uniqueness_score")),
            _fmt(j.get("overlap_hairpin_dg"), 2),
            _fmt(j.get("overlap_self_dimer_dg"), 2),
            j.get("overlap_off_target_count", ""),
            j.get("overlap_warnings", ""),
        ])

    files.append(_file(
        "junctions.csv", "text/csv",
        _csv_bytes(j_rows, j_header),
    ))

    # ── 3. assembled_construct.gb ─────────────────────────────────────────
    if assembled_seq:
        ann_for_gb = []
        for ann in viz.get("annotations", []):
            feat_type = (
                "primer_bind" if "primer" in ann.get("name", "").lower()
                else "misc_feature" if ann.get("type") == "overlap"
                else "misc_feature"
            )
            ann_for_gb.append({**ann, "feat_type": feat_type})

        gb_bytes = _make_genbank(
            seq=assembled_seq,
            name="Assembled_Construct",
            description=f"Gibson {assembly} assembly — {len(primers)} fragment(s)",
            annotations=ann_for_gb,
            topology=assembly,
        )
        files.append(_file("assembled_construct.gb", "application/octet-stream", gb_bytes))

    # ── 4. protocol.txt ────────────────────────────────────────────────────
    protocol_lines = [
        "Gibson Assembly Protocol",
        "=" * 40,
        f"Assembly type: {assembly.capitalize()}",
        f"Number of fragments: {len(primers)}",
        f"Number of junctions: {len(junctions)}",
        f"Assembled construct length: {len(assembled_seq)} bp",
        "",
        "PRIMERS",
        "-" * 40,
    ]
    for p in primers:
        frag = p.get("fragment", "?")
        if p.get("needs_primers"):
            fwd = p.get("forward_primer", "")
            rev = p.get("reverse_primer", "")
            fwd_tm = _fmt(p.get("forward_anneal_tm"))
            rev_tm = _fmt(p.get("reverse_anneal_tm"))
            protocol_lines += [
                f"Fragment: {frag}",
                f"  Forward primer ({len(fwd)} bp, Tm {fwd_tm}°C): {fwd}",
                f"  Reverse primer ({len(rev)} bp, Tm {rev_tm}°C): {rev}",
                "",
            ]
        else:
            protocol_lines += [f"Fragment: {frag}  — no primers needed ({p.get('reason', '')})", ""]

    protocol_lines += [
        "JUNCTIONS (OVERLAPS)",
        "-" * 40,
    ]
    for j in junctions:
        protocol_lines += [
            f"{j.get('from')} → {j.get('to')}",
            f"  Source: {j.get('source', '?')}",
            f"  Overlap: {j.get('overlap_sequence', '')} ({j.get('overlap_length', '?')} bp, Tm {_fmt(j.get('overlap_tm'))}°C, score {_fmt(j.get('overlap_score'))})",
            "",
        ]

    protocol_lines += [
        "RECOMMENDED PROTOCOL",
        "-" * 40,
        "1. Order primers as listed above (standard desalted, 25 nmol).",
        "2. Amplify each fragment by PCR using its designated primers.",
        "3. Verify PCR products by gel electrophoresis.",
        "4. Assemble using Gibson Assembly Master Mix (NEB E2611) or equivalent:",
        "   - Combine 50 ng of each fragment in equimolar ratios",
        "   - Total volume: 5 µL DNA + 15 µL 2× master mix = 20 µL",
        "   - Incubate 50°C, 60 minutes",
        "5. Transform 2 µL into competent cells.",
        "6. Verify correct assembly by colony PCR and Sanger sequencing.",
    ]

    files.append(_file(
        "protocol.txt", "text/plain",
        "\n".join(protocol_lines).encode("utf-8"),
    ))

    # ── parts_order.csv — bulk-aware cost breakdown ───────────────────────
    try:
        from .cloning.lab_profile import DEFAULT_LAB_PROFILE
        # Mine primers + synthesis fragments out of the gibson result.
        primer_rows: List[Dict[str, Any]] = []
        for pr in primers or []:
            if not pr.get("needs_primers"):
                continue
            frag = pr.get("fragment", "frag")
            fp = pr.get("forward_primer") or ""
            rp = pr.get("reverse_primer") or ""
            if fp:
                primer_rows.append({"name": f"{frag}_F", "sequence": fp})
            if rp:
                primer_rows.append({"name": f"{frag}_R", "sequence": rp})
        synth_rows: List[Dict[str, Any]] = []
        for sf in result.get("synthesis_fragments", []) or []:
            synth_rows.append({
                "name": sf.get("name", "synthesis"),
                "length_bp": sf.get("length") or sf.get("length_bp") or 0,
            })
        n_pcr = sum(1 for pr in (primers or []) if pr.get("needs_primers"))
        files.append(build_parts_order_csv(
            workflow="gibson",
            lab=DEFAULT_LAB_PROFILE,
            primers=primer_rows,
            synthesis_fragments=synth_rows,
            pcr_count=n_pcr,
            gel_count=n_pcr,
            sequencing_reads=DEFAULT_LAB_PROFILE.sequencing_reads_per_construct,
        ))
    except Exception as _e:
        # Never block file delivery on a costing error.
        pass

    return files



# ---------------------------------------------------------------------------
# parts_order.csv — unified per-construct cost breakdown
# ---------------------------------------------------------------------------

def _parts_order_rows(
    *,
    workflow: str,
    lab: Any,
    primers: List[Dict[str, Any]],
    synthesis_fragments: List[Dict[str, Any]],
    pcr_count: int = 0,
    gel_count: int = 0,
    sequencing_reads: int = 1,
) -> Tuple[List[List[Any]], List[str]]:
    """Build (rows, header) for the parts_order.csv breakdown.

    Each `primers` entry should carry: name, sequence (str). Length is
    derived from len(sequence). Each `synthesis_fragments` entry should
    carry: name, length_bp.
    """
    header = [
        "Item",
        "Type",
        "Length (bp)",
        "Quantity",
        "Bulk Price (USD)",
        "Bulk Pack (rxns)",
        "Per-Rxn / Per-Item Cost (USD)",
        "Total (USD)",
        "Notes",
    ]
    rows: List[List[Any]] = []
    grand_total = 0.0

    # Primers — length-scaled at $0.24/bp.
    for p in primers:
        seq = (p.get("sequence") or "").strip()
        L = len(seq)
        if L == 0:
            continue
        per_item = lab.primer_cost(L)
        rows.append([
            p.get("name", "primer"),
            "Primer",
            L,
            1,
            "",
            "",
            f"{per_item:.4f}",
            f"{per_item:.2f}",
            f"$0.24/bp · 25 nmol desalted",
        ])
        grand_total += per_item

    # Synthesis fragments — tiered $/bp.
    for s in synthesis_fragments:
        L = int(s.get("length_bp") or s.get("length") or 0)
        if L <= 0:
            continue
        rate = lab.synthesis_cost_per_bp(L)
        total = lab.synthesis_cost(L)
        rows.append([
            s.get("name", "synthesis_fragment"),
            "Synthesis",
            L,
            1,
            "",
            "",
            f"{rate:.4f}",
            f"{total:.2f}",
            f"${rate}/bp tier (0.5-1.8: $0.07 / 1.8-3.2: $0.08 / 3.2-5.0: $0.09)",
        ])
        grand_total += total

    # PCR (one Q5 rxn per fragment).
    if pcr_count > 0:
        q5 = lab.catalog["neb_q5_rxn"]
        per = q5.cost_per_rxn
        rows.append([
            q5.name,
            "Reagent (PCR)",
            "",
            pcr_count,
            f"{q5.bulk_price_usd:.2f}",
            f"{q5.bulk_rxns}",
            f"{per:.4f}",
            f"{per * pcr_count:.2f}",
            q5.usage_note,
        ])
        grand_total += per * pcr_count

    # Workflow reagent bundle (cells + plate + assembly + gel consumables).
    for r, qty in lab.reagent_lines(workflow):
        per = r.cost_per_rxn
        line_total = per * qty
        # Display "qty" as quantity per construct — fractional fine.
        qty_display = f"{qty:g}"
        rows.append([
            r.name,
            "Reagent",
            "",
            qty_display,
            f"{r.bulk_price_usd:.2f}",
            f"{r.bulk_rxns}",
            f"{per:.4f}",
            f"{line_total:.2f}",
            r.usage_note,
        ])
        grand_total += line_total

    # Sequencing — ONT plasmid read.
    if sequencing_reads > 0:
        sr = lab.catalog["ont_seq_read"]
        per = sr.cost_per_rxn
        rows.append([
            sr.name,
            "Sequencing",
            "",
            sequencing_reads,
            f"{sr.bulk_price_usd:.2f}",
            f"{sr.bulk_rxns}",
            f"{per:.4f}",
            f"{per * sequencing_reads:.2f}",
            sr.usage_note,
        ])
        grand_total += per * sequencing_reads

    # Grand total row.
    rows.append(["", "", "", "", "", "", "GRAND TOTAL", f"{grand_total:.2f}", ""])
    return rows, header


def build_parts_order_csv(
    *,
    workflow: str,
    lab: Any,
    primers: List[Dict[str, Any]],
    synthesis_fragments: List[Dict[str, Any]] = None,
    pcr_count: int = 0,
    gel_count: int = 0,
    sequencing_reads: int = 1,
) -> Dict[str, str]:
    """Return a {fileName, mimeType, dataBase64} payload for parts_order.csv."""
    synthesis_fragments = synthesis_fragments or []
    rows, header = _parts_order_rows(
        workflow=workflow,
        lab=lab,
        primers=primers,
        synthesis_fragments=synthesis_fragments,
        pcr_count=pcr_count,
        gel_count=gel_count,
        sequencing_reads=sequencing_reads,
    )
    return _file("parts_order.csv", "text/csv", _csv_bytes(rows, header))



def build_restriction_files(design: Dict[str, Any]) -> List[Dict[str, str]]:
    """Return downloadable files for a restriction-cloning design.

    Minimum viable shape so callers in chat.py can assign the result to a
    `files` list and emit primers.csv + assembled_construct.gb when the
    design carries them. Designs without primers / viz still return an
    empty list rather than crashing.
    """
    files: List[Dict[str, str]] = []
    if not isinstance(design, dict):
        return files

    primers = design.get("primers") or []
    if primers:
        header = [
            "Name", "Sequence", "Length (bp)", "Tm (°C)",
            "GC %", "Role", "Notes",
        ]
        rows = []
        for p in primers:
            seq = (p.get("sequence") or "").upper()
            rows.append([
                p.get("name", ""),
                seq,
                len(seq),
                _fmt(p.get("tm")),
                _fmt(p.get("gc_content")),
                p.get("role", ""),
                p.get("notes", ""),
            ])
        files.append(_file("primers.csv", "text/csv", _csv_bytes(rows, header)))

    viz = design.get("viz") or {}
    assembled_seq = viz.get("sequence", "")
    if assembled_seq:
        ann_for_gb = []
        for ann in viz.get("annotations", []) or []:
            ann_for_gb.append({
                "name": ann.get("name", "feature"),
                "start": ann.get("start", 0),
                "end": ann.get("end", 0),
                "direction": ann.get("direction", 1),
                "feat_type": "misc_feature",
            })
        topology = viz.get("topology") or "circular"
        gb_bytes = _make_genbank(
            assembled_seq,
            name=design.get("name", "Restriction_Product"),
            description=design.get("description", "Restriction cloning product"),
            annotations=ann_for_gb,
            topology=topology,
        )
        files.append(_file(
            "assembled_construct.gb",
            "application/octet-stream",
            gb_bytes,
        ))

    return files


# ---------------------------------------------------------------------------
# PCR files
# ---------------------------------------------------------------------------

def build_pcr_files(result: Dict[str, Any], template_name: str = "Template") -> List[Dict[str, str]]:
    """
    Generate downloadable files for a pcr_design result:
    - primers.csv
    - amplicon.gb
    """
    files: List[Dict[str, str]] = []

    template_seq = result.get("fragments_in", "")
    left_primer = result.get("left_primer", "")
    right_primer = result.get("right_primer", "")
    left_tm = result.get("left_tm")
    right_tm = result.get("right_tm")
    product_size = result.get("product_size", "")
    left_pos = result.get("left_pos") or {}
    right_pos = result.get("right_pos") or {}
    excluded = result.get("excluded_region") or {}

    # ── 1. primers.csv ─────────────────────────────────────────────────────
    header = [
        "Template", "Left Primer (5'→3')", "Right Primer (5'→3')",
        "Left Tm (°C)", "Right Tm (°C)", "Product Size (bp)",
        "Left Length (bp)", "Right Length (bp)",
        "Left Start", "Right End",
        "Excluded Start", "Excluded Length",
        "Left Hairpin Tm (°C)", "Right Hairpin Tm (°C)",
        "Left Self-Dimer Tm (°C)", "Right Self-Dimer Tm (°C)",
    ]
    l_scores = result.get("left_scores") or {}
    r_scores = result.get("right_scores") or {}
    rows = [[
        template_name,
        left_primer, right_primer,
        _fmt(left_tm), _fmt(right_tm),
        product_size,
        left_pos.get("len", ""), right_pos.get("len", ""),
        left_pos.get("start", ""),
        right_pos.get("start_3prime", ""),
        excluded.get("start", ""), excluded.get("length", ""),
        _fmt(l_scores.get("hairpin_th")), _fmt(r_scores.get("hairpin_th")),
        _fmt(l_scores.get("any_th")), _fmt(r_scores.get("any_th")),
    ]]
    files.append(_file("primers.csv", "text/csv", _csv_bytes(rows, header)))

    # ── 2. amplicon.gb ─────────────────────────────────────────────────────
    if template_seq:
        ann_for_gb: List[Dict[str, Any]] = []
        if left_pos:
            ann_for_gb.append({
                "name": f"Left_Primer ({left_pos.get('len', '')} bp, Tm {_fmt(left_tm)}°C)",
                "start": left_pos.get("start", 0),
                "end": left_pos.get("start", 0) + left_pos.get("len", 0),
                "direction": 1,
                "sequence": left_primer,
                "feat_type": "primer_bind",
            })
        if right_pos:
            r3 = right_pos.get("start_3prime", 0)
            r_len = right_pos.get("len", 0)
            ann_for_gb.append({
                "name": f"Right_Primer ({r_len} bp, Tm {_fmt(right_tm)}°C)",
                "start": max(0, r3 - r_len + 1),
                "end": r3 + 1,
                "direction": -1,
                "sequence": right_primer,
                "feat_type": "primer_bind",
            })
        if excluded.get("length"):
            ann_for_gb.append({
                "name": f"Excluded_Region ({excluded['length']} bp)",
                "start": excluded.get("start", 0),
                "end": excluded.get("start", 0) + excluded["length"],
                "direction": 0,
                "feat_type": "misc_feature",
            })

        gb_bytes = _make_genbank(
            seq=template_seq,
            name=template_name[:16],
            description=f"PCR — product {product_size} bp",
            annotations=ann_for_gb,
            topology="linear",
        )
        files.append(_file("amplicon.gb", "application/octet-stream", gb_bytes))

    return files


def build_batch_pcr_files(result: Dict[str, Any]) -> List[Dict[str, str]]:
    """
    Generate downloadable files for a multi_pcr_design result:
    - batch_primers.csv  (all templates in one file)
    - One amplicon.gb per template
    """
    files: List[Dict[str, str]] = []
    results = result.get("results", [])
    if not results:
        return files

    # ── batch_primers.csv ──────────────────────────────────────────────────
    header = [
        "Template", "Left Primer (5'→3')", "Right Primer (5'→3')",
        "Left Tm (°C)", "Right Tm (°C)", "Product Size (bp)",
        "Left Length (bp)", "Right Length (bp)",
        "Left Start", "Right End",
    ]
    rows = []
    for r in results:
        lp = r.get("left_pos") or {}
        rp = r.get("right_pos") or {}
        rows.append([
            r.get("template_name", f"Template_{r.get('template_index', '')}"),
            r.get("left_primer", ""), r.get("right_primer", ""),
            _fmt(r.get("left_tm")), _fmt(r.get("right_tm")),
            r.get("product_size", ""),
            lp.get("len", ""), rp.get("len", ""),
            lp.get("start", ""), rp.get("start_3prime", ""),
        ])
    files.append(_file("batch_primers.csv", "text/csv", _csv_bytes(rows, header)))

    # ── one GB per template ────────────────────────────────────────────────
    for r in results:
        tname = r.get("template_name", f"Template_{r.get('template_index', 0) + 1}")
        sub_files = build_pcr_files(r, template_name=tname)
        # Rename the GB file to include template name
        for f in sub_files:
            if f["fileName"] == "amplicon.gb":
                f["fileName"] = f"{tname}_amplicon.gb"
            elif f["fileName"] == "primers.csv":
                continue  # already in batch CSV
            files.append(f)

    return files


# ---------------------------------------------------------------------------
# Inventory Gibson files
# ---------------------------------------------------------------------------

def build_inv_gib_files(result: Dict[str, Any]) -> List[Dict[str, str]]:
    """
    Generate downloadable files for an inv_gib result:
    - inv_gib_results.csv
    """
    files: List[Dict[str, str]] = []

    fragments = result.get("fragments_in", [])
    summary = result.get("inv_gib_summary", {})

    # ── inv_gib_results.csv ────────────────────────────────────────────────
    header = [
        "Fragment", "Type", "Start (bp)", "End (bp)", "Direction",
        "Length (bp)", "Source Inventory",
    ]
    rows = []
    for frag in fragments:
        if not isinstance(frag, dict):
            continue
        name = frag.get("name", "")
        frag_type = "Synthesis Gap" if "Synthesis" in name else "Inventory Fragment"
        rows.append([
            name, frag_type,
            frag.get("target_start", frag.get("start", "")),
            frag.get("target_end", frag.get("end", "")),
            frag.get("source_orientation", "+"),
            frag.get("length_bp", ""),
            frag.get("source_inventory", ""),
        ])

    # Summary row
    if summary:
        rows.append([])
        rows.append(["SUMMARY", "", "", "", "", "", ""])
        rows.append(["Target length (bp)", "", "", "", "", summary.get("target_len", ""), ""])
        rows.append(["Inventory fragments", "", "", "", "", summary.get("emitted_inventory_fragments", ""), ""])
        rows.append(["Synthesis gaps", "", "", "", "", summary.get("synth_gap_count", ""), ""])
        covered = summary.get("covered_bp", 0)
        total = summary.get("target_len", 1) or 1
        rows.append(["Coverage (%)", "", "", "", "", f"{covered / total * 100:.1f}%", ""])

    files.append(_file(
        "inv_gib_results.csv", "text/csv",
        _csv_bytes(rows, header),
    ))

    return files


# ---------------------------------------------------------------------------
# Plasmid design files
# ---------------------------------------------------------------------------

def build_plasmid_design_files(viz: Dict[str, Any], title: str = "Designed_Plasmid") -> List[Dict[str, str]]:
    """
    Generate downloadable files for plasmid_design results:
    - designed_plasmid.gb (final assembled sequence with merged annotations)
    """
    files: List[Dict[str, str]] = []
    seq = (viz or {}).get("sequence", "") or ""
    if not seq:
        return files

    topology = (viz or {}).get("topology", "circular") or "circular"
    ann_for_gb: List[Dict[str, Any]] = []
    for ann in (viz or {}).get("annotations", []) or []:
        ann_for_gb.append({
            **ann,
            "feat_type": ann.get("feat_type", "misc_feature"),
        })

    gb_bytes = _make_genbank(
        seq=seq,
        name=title,
        description=f"AI plasmid design ({topology})",
        annotations=ann_for_gb,
        topology=topology,
    )
    files.append(_file("designed_plasmid.gb", "application/octet-stream", gb_bytes))
    return files

def build_sdm_files(result: Dict[str, Any]) -> List[Dict[str, str]]:
    """Build downloadable files for SDM design."""
    files = []
    
    primer_design = result.get("primer_design")
    if not primer_design:
        return files
    
    import csv
    from io import StringIO
    
    # Primers CSV
    csv_buffer = StringIO()
    writer = csv.writer(csv_buffer)
    writer.writerow(["Primer Name", "Sequence (5' to 3')", "Length (bp)", "Tm (C)", "Notes"])
    
    fwd = primer_design.get("forward_primer", "")
    rev = primer_design.get("reverse_primer", "")
    fwd_tm = primer_design.get("forward_tm", 0)
    rev_tm = primer_design.get("reverse_tm", 0)
    strategy = result.get("primer_strategy", "back_to_back")
    
    fwd_name = "SDM_Forward" if strategy != "single_primer" else "SDM_Mutagenic"
    writer.writerow([fwd_name, fwd, len(fwd), f"{fwd_tm:.1f}", "Order with standard desalting"])
    writer.writerow(["SDM_Reverse", rev, len(rev), f"{rev_tm:.1f}", "Order with standard desalting"])
    
    files.append({
        "name": "sdm_primers.csv",
        "content": csv_buffer.getvalue(),
        "type": "text/csv",
    })
    
    # Protocol summary
    mutation_type = result.get("mutation_type", "unknown")
    old_seq = primer_design.get("old_sequence", "")
    new_seq = primer_design.get("new_sequence", "")
    edit_start = primer_design.get("edit_start", 0)
    
    steps_text = ""
    for step in result.get("steps", []):
        steps_text += f"{step.get('step_number', 0)}. {step.get('description', '')}\n"
    
    protocol = f"""# Q5 Site-Directed Mutagenesis Protocol

## Mutation Details
- Type: {mutation_type}
- Position: {edit_start}
- Original sequence: {old_seq if old_seq else "(none - insertion)"}
- New sequence: {new_seq if new_seq else "(none - deletion)"}
- Strategy: {strategy}

## Primers
- Forward: 5'-{fwd}-3' (Tm: {fwd_tm:.1f}C)
- Reverse: 5'-{rev}-3' (Tm: {rev_tm:.1f}C)

## PCR Reaction Setup (25 uL)
| Component | Volume |
|-----------|--------|
| Q5 High-Fidelity 2X Master Mix | 12.5 uL |
| Forward Primer (10 uM) | 1.25 uL |
| Reverse Primer (10 uM) | 1.25 uL |
| Template DNA (1-25 ng) | 1 uL |
| Nuclease-free H2O | 9 uL |

## PCR Cycling
| Step | Temperature | Time | Cycles |
|------|-------------|------|--------|
| Initial Denaturation | 98C | 30 sec | 1 |
| Denaturation | 98C | 10 sec | 25 |
| Annealing | 68C | 30 sec | 25 |
| Extension | 72C | 30 sec/kb | 25 |
| Final Extension | 72C | 2 min | 1 |
| Hold | 4C | - | - |

## KLD Treatment
1. Mix 1 uL PCR product + 5 uL 2X KLD Reaction Buffer + 1 uL 10X KLD Enzyme Mix + 3 uL H2O
2. Incubate at room temperature for 5 minutes
3. Transform 5 uL into competent cells

## Protocol Steps
{steps_text}

## Estimated Cost: ${result.get('metrics', {}).get('total_cost_usd', 0):.2f}
## Estimated Time: {result.get('metrics', {}).get('total_calendar_days', 0):.1f} days
"""
    
    files.append({
        "name": "sdm_protocol.md",
        "content": protocol,
        "type": "text/markdown",
    })
    
    return files
