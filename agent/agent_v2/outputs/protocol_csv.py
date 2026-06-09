"""emit_protocol — third of four output emitters.

Writes protocol.csv: a per-step wet-lab playbook keyed by assembly method.
Six method templates: gibson, gateway, restriction, sdm, sgrna_gg, golden_gate.
Unknown methods fall back to a minimal generic template.

The agent can override any step (or append custom steps) by passing
`custom_steps`: a list of partial dicts indexed by `step_num`.
"""
from __future__ import annotations

import base64
import csv
import io
import pathlib
from typing import Any, Optional
from agent_v2.outputs import prefixed_filename, derive_descriptor


_CSV_HEADER = [
    "step_num", "category", "inputs", "output", "instrument",
    "time_min", "temp_C", "reagents", "notes",
]


def _step(category: str, inputs: str, output: str, instrument: str,
          time_min: int, temp_C: str, reagents: str, notes: str) -> dict[str, Any]:
    return {
        "category": category, "inputs": inputs, "output": output,
        "instrument": instrument, "time_min": time_min, "temp_C": temp_C,
        "reagents": reagents, "notes": notes,
    }


# All templates share the final transform/incubate/miniprep/verify/sanger tail
# except SDM (no ligation) and sgRNA-GG (single Sanger).
_TRANSFORM_BLOCK = [
    _step("transform", "2 uL reaction into 50 uL competent cells",
          "Transformed cells", "heat block, ice, 37 C shaker",
          90, "ice -> 42 C 30s -> ice -> 37 C SOC",
          "NEB 5-alpha or DH5alpha competent cells, SOC media",
          "30 min ice, 30 s heat shock at 42 C, 5 min ice, 60 min SOC recovery at 37 C"),
    _step("transform", "Recovered cells (100-200 uL)",
          "Plates with colonies", "incubator",
          960, "37", "LB agar + appropriate antibiotic",
          "Overnight ~16 h"),
    _step("prep", "Pick 4-8 single colonies",
          "Miniprep DNA", "centrifuge, vortex",
          240, "37 (overnight) / RT (prep)",
          "LB + antibiotic, miniprep kit (Qiagen / Zymo)",
          "5 mL overnight cultures -> miniprep, elute in 30 uL"),
    _step("verify", "Miniprep DNA + diagnostic enzymes",
          "Digest banding pattern", "thermocycler / heat block, gel rig",
          60, "37", "Restriction enzymes + CutSmart buffer",
          "Run on 1% agarose; compare to in-silico digest"),
    _step("sanger", "Validated minipreps + sequencing primers",
          "Sequence verification", "(external)",
          1440, "n/a", "Sequencing primers (3-5 covering junctions + insert)",
          "GeneWiz / Plasmidsaurus / Eurofins; submit at end of day"),
]


PROTOCOLS: dict[str, list[dict[str, Any]]] = {
    "gibson": [
        _step("PCR", "Each fragment template + overlap primers",
              "PCR product per fragment", "thermocycler",
              90, "98/56-65/72 cycling", "Q5 or Phusion polymerase, primers, dNTPs",
              "Repeat per fragment; design primers with 18-25 bp overlap"),
        _step("prep", "PCR products", "Cleaned amplicons", "spin column",
              15, "RT", "Zymo Clean & Concentrator-5",
              "Or DpnI digest if template was a plasmid"),
        _step("ligation", "Vector + cleaned inserts (equimolar, 0.02 pmol each)",
              "Gibson reaction", "thermocycler",
              60, "50", "NEBuilder HiFi 2x master mix",
              "1:2 vector:insert for 2-fragment; equimolar for 3+"),
    ] + _TRANSFORM_BLOCK,

    "gateway": [
        _step("PCR", "Template + attB1/attB2 (BP) or use entry clone (LR)",
              "attB-flanked PCR product or entry clone", "thermocycler",
              90, "98/56/72 cycling", "Q5 polymerase, attB primers (BP only)",
              "Skip if doing an LR with an existing entry clone"),
        _step("prep", "PCR product", "Cleaned amplicon", "spin column",
              15, "RT", "Zymo Clean & Concentrator-5", "BP only"),
        _step("ligation", "BP/LR clonase II + entry/dest vector",
              "Recombinant clone", "heat block / RT bench",
              60, "25", "Gateway BP or LR Clonase II Plus",
              "1 uL clonase, 5-10 fmol substrate"),
        _step("prep", "Recombinant reaction", "Inactivated reaction",
              "heat block", 10, "37", "Proteinase K (2 mg/mL)",
              "Stops the recombinase before transformation"),
    ] + _TRANSFORM_BLOCK,

    "restriction": [
        _step("digest", "Vector + insert + restriction enzymes",
              "Linearised vector + insert fragment", "thermocycler / heat block",
              60, "37", "NEB enzymes + CutSmart buffer",
              "Single or double digest depending on vector design"),
        _step("prep", "Digested vector + insert", "Gel-purified fragments",
              "gel rig + spin column",
              60, "RT", "1% agarose, QIAquick Gel Extraction",
              "Optional: CIP/SAP-treat vector to dephosphorylate ends"),
        _step("ligation", "Vector + insert (1:3 molar ratio)",
              "Ligation reaction", "heat block",
              60, "16 (overnight) or 25 (1 h)",
              "T4 DNA ligase + ligase buffer",
              "Overnight at 16 C improves yield for sticky ends"),
    ] + _TRANSFORM_BLOCK,

    "sdm": [
        _step("PCR", "Whole-plasmid template + mutagenic primer pair",
              "Linear mutated PCR product", "thermocycler",
              180, "98/55-65/72 cycling, 18 cycles",
              "Q5 or Phusion polymerase, primers, dNTPs",
              "Mutation in middle of primer; 10-15 bp matched flanks each side"),
        _step("digest", "PCR product + DpnI", "Methylated template removed",
              "heat block", 60, "37", "DpnI + CutSmart",
              "Removes original methylated plasmid template"),
    ] + _TRANSFORM_BLOCK,

    "sgrna_gg": [
        _step("prep", "Pair of complementary oligos with BsmBI/BsaI overhangs",
              "Annealed duplex", "thermocycler",
              10, "95 -> 25 ramp",
              "T4 ligase buffer (or NEB2 + 1 mM ATP)",
              "Anneal: 95 C for 5 min, ramp to 25 C at 0.1 C/s"),
        _step("prep", "Annealed duplex + T4 PNK", "Phosphorylated duplex",
              "heat block", 30, "37",
              "T4 PNK + ATP + ligase buffer",
              "Skip if oligos were ordered with 5' phosphate"),
        _step("ligation", "Vector + duplex + Type-IIs enzyme + T4 ligase",
              "Recombinant vector", "thermocycler",
              60, "37/16 cycling x 25, then 50/80",
              "BsmBI/Esp3I or BsaI + T4 ligase + buffer + ATP",
              "Cycle: 37 C 5 min / 16 C 5 min x 25; final 50 C 5 min, 80 C 5 min"),
    ] + _TRANSFORM_BLOCK[:-1] + [
        _step("sanger", "Validated minipreps + protospacer-flanking primer",
              "Guide sequence verification", "(external)",
              1440, "n/a", "U6 forward primer (or scaffold reverse)",
              "Single Sanger read across the protospacer is enough"),
    ],

    "golden_gate": [
        _step("PCR", "Each fragment template + Type-IIs-flanked primers",
              "PCR product per fragment", "thermocycler",
              90, "98/56-65/72 cycling", "Q5 / Phusion polymerase, primers",
              "Internal Type-IIs sites must be removed/avoided"),
        _step("prep", "PCR products", "Cleaned amplicons", "spin column",
              15, "RT", "Zymo Clean & Concentrator-5", ""),
        _step("ligation", "Vector + fragments + Type-IIs enzyme + T4 ligase",
              "One-pot Golden Gate reaction", "thermocycler",
              120, "37/16 cycling x 30, then 55/80",
              "BsaI/BsmBI/Esp3I + T4 ligase + buffer + ATP",
              "Cycle: 37 C 1.5 min / 16 C 3 min x 30; final 55 C 5 min, 80 C 5 min"),
    ] + _TRANSFORM_BLOCK,

    "none": [],
}


def _generic_template() -> list[dict[str, Any]]:
    return [
        _step("prep", "(method not recognised)", "(see workflow_trace.txt)",
              "n/a", 0, "n/a", "n/a",
              "Falls back to a generic transform + verify tail."),
    ] + _TRANSFORM_BLOCK


def _apply_custom_steps(steps: list[dict[str, Any]],
                         custom: Optional[list[dict[str, Any]]]) -> list[dict[str, Any]]:
    if not custom:
        return steps
    out = [dict(s) for s in steps]
    for c in custom:
        idx = c.get("step_num")
        if idx is None or idx < 1:
            out.append({k: v for k, v in c.items() if k != "step_num"})
            continue
        if 1 <= idx <= len(out):
            for k, v in c.items():
                if k == "step_num":
                    continue
                out[idx - 1][k] = v
        else:
            out.append({k: v for k, v in c.items() if k != "step_num"})
    return out


def _build_csv(steps: list[dict[str, Any]]) -> str:
    buf = io.StringIO()
    w = csv.writer(buf, quoting=csv.QUOTE_MINIMAL)
    w.writerow(_CSV_HEADER)
    for i, s in enumerate(steps, start=1):
        w.writerow([
            i,
            s.get("category", ""),
            s.get("inputs", ""),
            s.get("output", ""),
            s.get("instrument", ""),
            s.get("time_min", ""),
            s.get("temp_C", ""),
            s.get("reagents", ""),
            s.get("notes", ""),
        ])
    return buf.getvalue()


async def emit_protocol(
    args: dict[str, Any],
    registry: Any,
    *,
    output_dir: Optional[str] = None,
) -> dict[str, Any]:
    method = (args.get("assembly_method") or "none").lower()
    template = PROTOCOLS.get(method)
    if template is None:
        template = _generic_template()
    steps = _apply_custom_steps(template, args.get("custom_steps"))
    csv_text = _build_csv(steps)

    _descriptor = derive_descriptor(args)
    file_envelope = {
        "fileName": prefixed_filename("protocol.csv", _descriptor),
        "dataBase64": base64.b64encode(csv_text.encode("utf-8")).decode("ascii"),
    }

    written_path: Optional[str] = None
    if output_dir is not None:
        out = pathlib.Path(output_dir) / prefixed_filename("protocol.csv", _descriptor)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(csv_text)
        written_path = str(out)

    total_min = sum(int(s.get("time_min") or 0) for s in steps)
    categories = sorted({s.get("category", "") for s in steps if s.get("category")})

    return {
        "ok": True,
        "file": file_envelope,
        "method": method,
        "n_steps": len(steps),
        "total_time_min": total_min,
        "categories": categories,
        "written_path": written_path,
    }
