"""emit_parts_order — second of four output emitters.

Writes parts_order.csv: an 8-column inventory-and-ordering manifest. Each
part is classified as one of:
  - inventory: already on hand (resolved from an uploaded .gb)
  - order_addgene: Addgene plasmid (has addgene_id)
  - order_synthesis: >=40 bp synthesis fragment (Twist gBlock by default)
  - order_oligo: <40 bp oligo (IDT by default)

Vendor / cost / lead-time defaults are rough but overridable per part
(vendor_hint, cost_usd, lead_time_days, notes).
"""
from __future__ import annotations

import csv

def _lookup_sequence_for_part(part: dict, registry) -> str:
    """When the part references a registered attachment (by
    attachment_id, by name, or by source_plasmid) AND we have the
    sequence in the registry, return it. Otherwise return "" so the
    parts_order.csv row gets a sequence column populated whenever we
    actually know it, instead of always being blank."""
    if registry is None:
        return ""
    aid = part.get("attachment_id") or part.get("source_attachment_id")
    if aid and hasattr(registry, "get"):
        att = registry.get(aid)
        if att and getattr(att, "sequence", None):
            return att.sequence
    pname = (part.get("name") or "").strip().lower()
    if pname and hasattr(registry, "items"):
        for a in registry.items.values():
            if (getattr(a, "name", "") or "").strip().lower() == pname:
                return a.sequence
    return ""



import io
import pathlib
from typing import Any, Optional
from agent_v2.outputs import prefixed_filename, derive_descriptor


_CSV_HEADER = [
    "part_id", "name", "source", "length_bp", "vendor",
    "est_cost_usd", "est_lead_time_days", "notes",
]


def _classify_source(part: dict[str, Any]) -> str:
    if part.get("addgene_id"):
        return "order_addgene"
    origin = (part.get("origin") or "").lower()
    if origin in {"inventory", "user_supplied", "input"}:
        return "inventory"
    length = int(part.get("length_bp") or 0)
    if origin == "designed_oligo" or length < 40:
        return "order_oligo"
    return "order_synthesis"


def _default_vendor_cost_lead(source: str, length: int, part: dict[str, Any]) -> tuple[str, float, int]:
    """Return (vendor, est_cost_usd, est_lead_time_days)."""
    vendor = part.get("vendor_hint") or ""
    cost = part.get("cost_usd")
    lead = part.get("lead_time_days")
    if source == "inventory":
        return (vendor or "(in inventory)", float(cost if cost is not None else 0.0),
                int(lead if lead is not None else 0))
    if source == "order_addgene":
        return (vendor or f"Addgene #{part.get('addgene_id', '?')}",
                float(cost if cost is not None else 85.0),
                int(lead if lead is not None else 14))
    if source == "order_oligo":
        # ~$0.18/bp typical for desalted oligos at 25 nmol scale
        return (vendor or "IDT oligo",
                float(cost if cost is not None else round(0.18 * length, 2)),
                int(lead if lead is not None else 2))
    # order_synthesis
    return (vendor or "Twist gBlock",
            float(cost if cost is not None else round(0.07 * length, 2)),
            int(lead if lead is not None else 5))


def _build_csv(parts: list[dict[str, Any]]) -> str:
    buf = io.StringIO()
    w = csv.writer(buf, quoting=csv.QUOTE_MINIMAL)
    w.writerow(_CSV_HEADER)
    for p in parts:
        length = int(p.get("length_bp") or 0)
        source = _classify_source(p)
        vendor, cost, lead = _default_vendor_cost_lead(source, length, p)
        w.writerow([
            p.get("part_id", ""),
            p.get("name", ""),
            source,
            length,
            vendor,
            f"{cost:.2f}",
            lead,
            p.get("notes", ""),
        ])
    return buf.getvalue()




def _auto_derive_parts(args: dict[str, Any]) -> list[dict[str, Any]]:
    """Build a parts list from cached emit_guides_csv data so the LLM
    doesn't have to compose it on every run. Each cloning_oligo entry
    yields top + bottom IDT oligo parts; each primer yields a forward +
    reverse IDT oligo part. Empty lists are tolerated. Used as a fallback
    when emit_parts_order is called without an explicit parts arg."""
    parts: list[dict[str, Any]] = []

    for c in (args.get("cloning_oligos") or []):
        if not isinstance(c, dict):
            continue
        name = c.get("name") or "spacer"
        notes_base = c.get("notes") or ""
        top = c.get("oligo_top") or ""
        bot = c.get("oligo_bottom") or ""
        if top:
            parts.append({
                "part_id": f"{name}_top",
                "name": f"{name}_top",
                "length_bp": len(top),
                "origin": "designed_oligo",
                "sequence": top,
                "notes": notes_base or "spacer cloning top oligo",
            })
        if bot:
            parts.append({
                "part_id": f"{name}_bottom",
                "name": f"{name}_bottom",
                "length_bp": len(bot),
                "origin": "designed_oligo",
                "sequence": bot,
                "notes": notes_base or "spacer cloning bottom oligo",
            })

    for p in (args.get("primers") or []):
        if not isinstance(p, dict):
            continue
        label = p.get("pair_label") or p.get("name_fwd") or "primer"
        app = (p.get("application") or "primer").lower()
        fwd_seq = p.get("left_primer") or p.get("primer_sequence")
        rev_seq = p.get("right_primer")
        left_tm = p.get("left_tm") or p.get("tm")
        right_tm = p.get("right_tm") or p.get("tm")
        product = p.get("product_size")
        if fwd_seq:
            parts.append({
                "part_id": f"{label}_{app}_F",
                "name": p.get("name_fwd") or f"{label}_F",
                "length_bp": len(fwd_seq),
                "origin": "designed_oligo",
                "sequence": fwd_seq,
                "notes": f"{app} fwd; Tm={left_tm or '?'}; amplicon={product or '?'} bp",
            })
        if rev_seq:
            parts.append({
                "part_id": f"{label}_{app}_R",
                "name": p.get("name_rev") or f"{label}_R",
                "length_bp": len(rev_seq),
                "origin": "designed_oligo",
                "sequence": rev_seq,
                "notes": f"{app} rev; Tm={right_tm or '?'}; amplicon={product or '?'} bp",
            })

    return parts


async def emit_parts_order(
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,
) -> dict[str, Any]:
    parts = args.get("parts") or []
    if not isinstance(parts, list):
        return {"ok": False, "error": "parts must be a list"}
    auto_derived = False
    if not parts:
        # Fall back to the cached emit_guides_csv payload so the LLM
        # doesn't have to compose the parts array on every run.
        parts = _auto_derive_parts(args)
        auto_derived = True
    csv_text = _build_csv(parts)

    import base64
    _descriptor = derive_descriptor(args)
    file_envelope = {
        "fileName": prefixed_filename("parts_order.csv", _descriptor),
        "dataBase64": base64.b64encode(csv_text.encode("utf-8")).decode("ascii"),
    }

    written_path: Optional[str] = None
    if output_dir is not None:
        out = pathlib.Path(output_dir) / prefixed_filename("parts_order.csv", _descriptor)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(csv_text)
        written_path = str(out)

    counts: dict[str, int] = {}
    for p in parts:
        s = _classify_source(p)
        counts[s] = counts.get(s, 0) + 1

    return {
        "ok": True,
        "file": file_envelope,
        "n_parts": len(parts),
        "auto_derived_parts": auto_derived,
        "source_counts": counts,
        "written_path": written_path,
    }
