"""
Adapters that turn raw chat.py inputs (uploads + parsed message) into a
WorkflowInput. One per dispatch shape; chat.py picks the right one based on
intent + which uploads are present.

Step 3 of the cloning-workflow input normalization (workflow_input.py).
Adapters land progressively as each handler migrates; the legacy path stays
behind `if workflow_input is None: ...` until all eight handlers are done.

Each adapter is responsible only for *reshaping* — never for invoking the
annotation pipeline, the part resolver, or the cloning router. Those are the
handler's job. The adapter promises:

    raw bytes/dicts → a populated PlasmidView + WorkflowInput

…and nothing more.
"""
from __future__ import annotations

from io import StringIO
from typing import Any, Dict, List, Optional

from fastapi import UploadFile

import re

from .canonical_request import AssemblyPlan, CloningRequest
from .workflow_input import (
    PlasmidView,
    WorkflowInput,
    build_from_single_upload,
)


# ---------------------------------------------------------------------------
# annotate_gb — single uploaded GenBank, no parts, no inventory
# ---------------------------------------------------------------------------

async def build_for_annotate_gb(
    *,
    file: UploadFile,
    canonical_request: Optional[CloningRequest],
    session_id: str,
) -> WorkflowInput:
    """Build a WorkflowInput for the `annotate_gb` intent.

    Reads the uploaded GenBank, parses it for sequence + topology + name, and
    wraps the result as a PlasmidView with role="target". Parse errors are
    carried on `target.metadata["parse_error"]` so the handler can format the
    user-facing message instead of the adapter raising. Defaults match the
    legacy chat.py block exactly:

        - plasmid_name default: \"Annotated Plasmid\"  (NOT the filename)
        - sequence on parse failure: \"\"
        - is_circular default: False

    The annotation pipeline call (`/plannotate/annotate_sequence_llm`) stays
    in the handler — this adapter only ships the inputs.
    """
    raw = await file.read()
    try:
        gb_text = raw.decode("utf-8")
    except UnicodeDecodeError:
        gb_text = raw.decode("latin-1", errors="replace")

    sequence = ""
    is_circular = False
    plasmid_name = "Annotated Plasmid"
    parse_error: Optional[str] = None
    try:
        from Bio import SeqIO
        record = SeqIO.read(StringIO(gb_text), "genbank")
        sequence = str(record.seq)
        is_circular = record.annotations.get("topology", "linear") == "circular"
        plasmid_name = record.name or "Annotated Plasmid"
    except Exception as exc:
        parse_error = str(exc)

    metadata = {}
    if parse_error:
        metadata["parse_error"] = parse_error

    target_view = PlasmidView(
        name=plasmid_name,
        sequence=sequence,
        role="target",
        source="uploaded",
        source_file=getattr(file, "filename", None),
        topology="circular" if is_circular else "linear",
        gb_text=gb_text,
        metadata=metadata,
    )
    return build_from_single_upload(
        canonical_request=canonical_request,
        intent="annotate_gb",
        session_id=session_id,
        target_view=target_view,
    )


# ---------------------------------------------------------------------------
# gibson_design — extracted DNA fragments + intent-LLM gibson params
# ---------------------------------------------------------------------------

def build_for_gibson_design(
    *,
    seq_data: Dict[str, Any],
    intent_result: Dict[str, Any],
    canonical_request: Optional[CloningRequest],
    session_id: str,
) -> WorkflowInput:
    """Build a WorkflowInput from the legacy Gibson dispatch path.

    `seq_data` is the output of `extractors.extract_sequences(message)`:
        - `sequences`: ordered list of unique uppercase DNA strings
        - `fragments`: dict of {name: seq} when extract found named labels
                       (Frag1: ATGC..., Template_2: ..., etc.)
        - `count`: total

    Named fragments win over positional ones — that mirrors the legacy
    `build_fragment_objects(seq_data)` helper.

    Tuning knobs come from `intent_result["gibson_design"]["primer_params"]`
    + `canonical_request.gibson_params` if present (canonical wins on
    collision; canonical_request is the normalized + reconciled view).
    """
    sequences = list(seq_data.get("sequences") or [])
    fragments_dict = dict(seq_data.get("fragments") or {})

    parts: List[PlasmidView] = []
    if fragments_dict:
        for i, (name, seq) in enumerate(fragments_dict.items()):
            parts.append(PlasmidView(
                name=name,
                sequence=seq,
                role="fragment",
                source="prompt_sequence",
                topology="linear",
                metadata={"fragment_index": i},
            ))
    else:
        for i, seq in enumerate(sequences):
            parts.append(PlasmidView(
                name=f"Fragment_{i+1}",
                sequence=seq,
                role="fragment",
                source="prompt_sequence",
                topology="linear",
                metadata={"fragment_index": i},
            ))

    gib_params = intent_result.get("gibson_design") or {}
    primer_params: Dict[str, Any] = dict(gib_params.get("primer_params") or {})
    topology = gib_params.get("assembly") or "circular"
    if topology not in ("linear", "circular"):
        topology = "circular"

    if canonical_request and canonical_request.gibson_params:
        cr_params = canonical_request.gibson_params or {}
        cr_pp = cr_params.get("primer_params")
        if isinstance(cr_pp, dict):
            primer_params.update(cr_pp)
        if cr_params.get("assembly") in ("linear", "circular"):
            topology = cr_params["assembly"]

    workflow_args: Dict[str, Any] = {
        "primer_params": primer_params,
        "assembly": topology,
    }

    constraints = list((canonical_request.constraints if canonical_request else []) or [])
    notes = list((canonical_request.normalizer_notes if canonical_request else []) or [])

    return WorkflowInput(
        intent="gibson_design",
        workflow="gibson_design",
        session_id=session_id,
        target=None,
        parts=parts,
        vector=None,
        constraints=constraints,
        assembly=AssemblyPlan(method="gibson", topology=topology),
        workflow_args=workflow_args,
        provenance={"normalizer_notes": notes},
    )


# ---------------------------------------------------------------------------
# golden_gate_primer_design — extracted DNA fragments OR KB-resolved parts
# ---------------------------------------------------------------------------

def build_for_golden_gate_primer_design(
    *,
    message: str,
    seq_data: Dict[str, Any],
    intent_result: Dict[str, Any],
    canonical_request: Optional[CloningRequest],
    session_id: str,
) -> WorkflowInput:
    """Build a WorkflowInput for the golden_gate_primer_design dispatch.

    Fragment resolution mirrors the legacy chat.py block:
      1. If `seq_data["count"] >= 2`, use the extracted sequences as positional
         fragments (Fragment_1, Fragment_2, ...).
      2. Otherwise, extract part candidates from the prompt and resolve them
         against the pLannotate feature KB. Carry the audit trail
         (`identification_report`) on `provenance` so the handler can render
         the "Identified from pLannotate feature KB:" prefix.

    Workflow knobs (workflow_type, enzyme, target_tm) come from
    `intent_result["golden_gate"]`, with `canonical_request.golden_gate_params`
    winning on collision.
    """
    fragments: List[Dict[str, Any]] = []
    fragments_source = "prompt_sequence"
    identification_report: Optional[Dict[str, Any]] = None

    if (seq_data.get("count") or 0) >= 2:
        for i, seq in enumerate(seq_data.get("sequences") or []):
            fragments.append({"name": f"Fragment_{i+1}", "sequence": seq})
    else:
        # Lazy import to avoid module-load circular: chat.py imports adapters,
        # adapters import KB helpers from chat.py at call time. Python resolves
        # this because import happens after both modules are loaded.
        from .chat import extract_part_candidates, identify_features_from_kb
        candidates = extract_part_candidates(message)
        identified = identify_features_from_kb(candidates) if candidates else []
        fragments = [
            {"name": f["name"], "sequence": f["sequence"]}
            for f in identified if f.get("sequence")
        ]
        identified_queries = {f["query"] for f in identified}
        unresolved_candidates = [c for c in candidates if c["name"] not in identified_queries]
        identification_report = {
            "candidates": candidates,
            "identified": identified,
            "unresolved": unresolved_candidates,
        }
        fragments_source = "kb_lookup"

    parts: List[PlasmidView] = []
    for i, frag in enumerate(fragments):
        parts.append(PlasmidView(
            name=frag.get("name") or f"Fragment_{i+1}",
            sequence=frag.get("sequence") or "",
            role="fragment",
            source=fragments_source,
            topology="linear",
            metadata={"fragment_index": i},
        ))

    gg_params = intent_result.get("golden_gate") or {}
    workflow_type = gg_params.get("workflow_type") or "multi_fragment"
    enzyme = gg_params.get("enzyme") or "BsaI"
    target_tm = gg_params.get("target_tm")

    if canonical_request and canonical_request.golden_gate_params:
        cr_params = canonical_request.golden_gate_params or {}
        if cr_params.get("workflow_type"):
            workflow_type = cr_params["workflow_type"]
        if cr_params.get("enzyme"):
            enzyme = cr_params["enzyme"]
        if cr_params.get("target_tm") is not None:
            target_tm = cr_params["target_tm"]

    workflow_args: Dict[str, Any] = {
        "workflow_type": workflow_type,
        "enzyme": enzyme,
    }
    if target_tm is not None:
        workflow_args["target_tm"] = float(target_tm)

    constraints = list((canonical_request.constraints if canonical_request else []) or [])
    notes = list((canonical_request.normalizer_notes if canonical_request else []) or [])
    provenance: Dict[str, Any] = {"normalizer_notes": notes}
    if identification_report is not None:
        provenance["identification_report"] = identification_report

    return WorkflowInput(
        intent="golden_gate_primer_design",
        workflow="golden_gate_primer_design",
        session_id=session_id,
        target=None,
        parts=parts,
        vector=None,
        constraints=constraints,
        assembly=AssemblyPlan(method="golden_gate", topology="circular"),
        workflow_args=workflow_args,
        provenance=provenance,
    )


# ---------------------------------------------------------------------------
# restriction_cloning — vector file + KB-lookup insert + enzyme override
# ---------------------------------------------------------------------------

async def build_for_restriction_cloning(
    *,
    file: Optional[UploadFile],
    inventory_files: Optional[List[UploadFile]],
    has_target: bool,
    has_inventory: bool,
    seq_data: Dict[str, Any],
    intent_result: Dict[str, Any],
    canonical_request: Optional[CloningRequest],
    session_id: str,
) -> WorkflowInput:
    """Build a WorkflowInput for restriction_cloning.

    Vector resolution: `file` (target) preferred, first `inventory_files` as
    fallback. Reads bytes, parses the GenBank for sequence + topology + locus
    name. Vector name precedence matches the legacy block: prompt-provided →
    parsed locus → filename → "vector".

    Insert resolution: pull (insert_name, insert_seq, vector_name,
    enzyme_override) via the existing `_restriction_inputs_from_canonical`
    helper (lazy-imported from chat.py — same module-level helper the legacy
    block used). If only a name is given, KB-look it up via
    `identify_features_from_kb`. Resolved insert lands in `parts` with
    role="insert" and source="kb_lookup" or "prompt_sequence".

    Enzyme override is carried on workflow_args["enzyme_override"] as the
    designer expects: a `(left, right)` tuple or None.

    KB record (if any) is stashed on `provenance["identified_record"]` so the
    handler can render the audit prefix and trigger the synthesis-first path
    for back-translated KB hits.
    """
    # Lazy imports to avoid the chat.py / adapters circular at module load.
    from .chat import _restriction_inputs_from_canonical, identify_features_from_kb

    insert_name, insert_seq_param, vector_name, enzyme_override = (
        _restriction_inputs_from_canonical(canonical_request, intent_result, seq_data)
    )

    # ----- Vector file → PlasmidView with gb_text -----
    vector_view: Optional[PlasmidView] = None
    if has_target and file is not None:
        vector_file_obj = file
    elif has_inventory and inventory_files:
        vector_file_obj = inventory_files[0]
    else:
        vector_file_obj = None

    if vector_file_obj is not None:
        raw = await vector_file_obj.read()
        try:
            vector_gb_text = raw.decode("utf-8")
        except UnicodeDecodeError:
            vector_gb_text = raw.decode("latin-1", errors="replace")

        v_seq = ""
        v_topology = "circular"
        v_locus_name: Optional[str] = None
        try:
            from Bio import SeqIO
            from io import StringIO
            rec = SeqIO.read(StringIO(vector_gb_text), "genbank")
            v_seq = str(rec.seq)
            v_topology = (
                "circular"
                if rec.annotations.get("topology", "linear") == "circular"
                else "linear"
            )
            v_locus_name = rec.name or None
        except Exception:
            pass

        if not vector_name:
            vector_name = (
                v_locus_name
                or getattr(vector_file_obj, "filename", None)
                or "vector"
            )

        vector_view = PlasmidView(
            name=vector_name,
            sequence=v_seq,
            role="vector",
            source="uploaded",
            source_file=getattr(vector_file_obj, "filename", None),
            topology=v_topology,
            gb_text=vector_gb_text,
        )

    # ----- Insert resolution -----
    parts: List[PlasmidView] = []
    identified_record: Optional[Dict[str, Any]] = None
    insert_seq: Optional[str] = None
    resolved_label: Optional[str] = None

    if insert_seq_param and all(ch in "ACGTN" for ch in insert_seq_param.upper()):
        insert_seq = re.sub(r"[^ACGT]", "", insert_seq_param.upper())
        resolved_label = insert_name or "insert"
    elif insert_name:
        candidates = [{"name": insert_name, "feature_type": None}]
        identified = identify_features_from_kb(candidates)
        if identified and identified[0].get("sequence"):
            identified_record = identified[0]
            insert_seq = identified_record["sequence"]
            resolved_label = identified_record["name"]

    if insert_seq and resolved_label:
        meta: Dict[str, Any] = {}
        if identified_record:
            meta = {k: v for k, v in identified_record.items() if k != "sequence"}
        parts.append(PlasmidView(
            name=resolved_label,
            sequence=insert_seq,
            role="insert",
            source="kb_lookup" if identified_record else "prompt_sequence",
            topology="linear",
            metadata=meta,
        ))

    workflow_args: Dict[str, Any] = {
        "vector_name": vector_name or "",
        "insert_name": insert_name,
    }
    if enzyme_override is not None:
        workflow_args["enzyme_override"] = enzyme_override

    provenance: Dict[str, Any] = {
        "normalizer_notes": list((canonical_request.normalizer_notes if canonical_request else []) or []),
    }
    if identified_record is not None:
        provenance["identified_record"] = identified_record

    constraints = list((canonical_request.constraints if canonical_request else []) or [])

    return WorkflowInput(
        intent="restriction_cloning",
        workflow="restriction_cloning",
        session_id=session_id,
        target=None,
        parts=parts,
        vector=vector_view,
        constraints=constraints,
        assembly=AssemblyPlan(method="restriction", topology="circular"),
        workflow_args=workflow_args,
        provenance=provenance,
    )


# ---------------------------------------------------------------------------
# sgrna_golden_gate — gRNA + vector (uploaded or default lentiCRISPR v2)
# ---------------------------------------------------------------------------

async def build_for_sgrna_golden_gate(
    *,
    message: str,
    file: Optional[UploadFile],
    inventory_files: Optional[List[UploadFile]],
    has_target: bool,
    has_inventory: bool,
    seq_data: Dict[str, Any],
    intent_result: Dict[str, Any],
    canonical_request: Optional[CloningRequest],
    session_id: str,
) -> WorkflowInput:
    r"""Build a WorkflowInput for sgrna_golden_gate.

    Resolves (in priority order, mirroring the legacy block):
      - gRNA sequence: seq_data first → regex on message → LLM-parsed
        sgrna_params.grna_sequence. Adapter does NOT validate length /
        alphabet — handler owns the user-facing error messages.
      - gRNA name: regex `\(([A-Za-z0-9_-]+)\)` on message; default "sgRNA".
      - Vector: target file → first inventory_files → default lentiCRISPR v2
        loaded via cloning.sgrna_oligo_designer.load_lenticrispr_v2().
      - Enzyme: intent.sgrna.enzyme → canonical_request.sgrna_params.enzyme
        → "BsmBI".

    Outputs:
      - workflow_input.vector → PlasmidView with sequence, gb_text, and
        metadata["features"] (the parse_genbank_features list the assembler
        consumes via its `original_features` arg).
      - workflow_input.workflow_args → grna_sequence, grna_name, enzyme.
    """
    from .cloning.sgrna_oligo_designer import (
        load_lenticrispr_v2,
        parse_genbank_features,
        parse_genbank_sequence,
    )

    sgrna_params = intent_result.get("sgrna") or {}

    # ----- gRNA sequence (priority chain) -----
    grna_seq: Optional[str] = None
    if (seq_data.get("count") or 0) >= 1:
        first_seq = ((seq_data.get("sequences") or [""])[0] or "").upper()
        if 17 <= len(first_seq) <= 30 and set(first_seq) <= set("ACGTN"):
            grna_seq = first_seq
    if not grna_seq:
        m = re.search(r"\b([ACGT]{17,30})\b", message.upper())
        if m:
            grna_seq = m.group(1)
    if not grna_seq:
        llm_grna = sgrna_params.get("grna_sequence")
        if isinstance(llm_grna, str):
            llm_grna = llm_grna.strip().upper()
            if 17 <= len(llm_grna) <= 30 and set(llm_grna) <= set("ACGTN"):
                grna_seq = llm_grna

    # ----- gRNA name -----
    grna_name = "sgRNA"
    name_match = re.search(r"\(([A-Za-z0-9_-]+)\)", message)
    if name_match:
        grna_name = name_match.group(1)

    # ----- Vector (uploaded → inventory → default) -----
    plasmid_file = None
    if has_target and file is not None:
        plasmid_file = file
    elif has_inventory and inventory_files:
        plasmid_file = inventory_files[0]

    vector_seq: Optional[str] = None
    vector_features = None
    vector_name = "lentiCRISPR v2"
    vector_source = "uploaded"
    vector_filename: Optional[str] = None
    vector_gb_text: Optional[str] = None

    if plasmid_file is not None:
        try:
            await plasmid_file.seek(0)
        except Exception:
            pass
        try:
            raw = await plasmid_file.read()
            if raw:
                try:
                    vector_gb_text = raw.decode("utf-8")
                except UnicodeDecodeError:
                    vector_gb_text = raw.decode("latin-1", errors="replace")
                vector_seq = parse_genbank_sequence(vector_gb_text)
                vector_features = parse_genbank_features(vector_gb_text)
                vector_name = plasmid_file.filename or "uploaded vector"
                vector_filename = plasmid_file.filename
        except Exception:
            vector_seq = None  # fall through to default

    if not vector_seq:
        try:
            vector_gb_text = load_lenticrispr_v2()
            vector_seq = parse_genbank_sequence(vector_gb_text)
            vector_features = parse_genbank_features(vector_gb_text)
            vector_name = "lentiCRISPR v2"
            vector_source = "default"
            vector_filename = None
        except Exception:
            vector_seq = None
            vector_gb_text = None
            vector_features = None

    vector_view: Optional[PlasmidView] = None
    if vector_seq:
        vector_view = PlasmidView(
            name=vector_name,
            sequence=vector_seq,
            role="vector",
            source=vector_source,
            source_file=vector_filename,
            topology="circular",
            gb_text=vector_gb_text,
            metadata={"features": vector_features or []},
        )

    # ----- Enzyme -----
    # Only resolve to a specific enzyme if the user (or canonical request)
    # explicitly named one. Otherwise leave as None so predesign tries all
    # supported Type IIs enzymes (BsmBI / BbsI / BsaI) — vectors like
    # px330 / px459 / pU6-BbsI cassettes would otherwise be rejected by a
    # silent BsmBI default.
    enzyme: Optional[str] = sgrna_params.get("enzyme") or None
    if canonical_request and canonical_request.sgrna_params:
        cr = canonical_request.sgrna_params or {}
        if cr.get("enzyme"):
            enzyme = cr["enzyme"]

    workflow_args: Dict[str, Any] = {
        "grna_sequence": grna_seq,
        "grna_name": grna_name,
        "enzyme": enzyme,
    }

    constraints = list((canonical_request.constraints if canonical_request else []) or [])
    notes = list((canonical_request.normalizer_notes if canonical_request else []) or [])

    return WorkflowInput(
        intent="sgrna_golden_gate",
        workflow="sgrna_golden_gate",
        session_id=session_id,
        target=None,
        parts=[],
        vector=vector_view,
        constraints=constraints,
        assembly=AssemblyPlan(method="golden_gate", topology="circular"),
        workflow_args=workflow_args,
        provenance={"normalizer_notes": notes},
    )


# ---------------------------------------------------------------------------
# plasmid_design — natural-language description → semantic-retrieval handler
# ---------------------------------------------------------------------------

def build_for_plasmid_design(
    *,
    message: str,
    intent_result: Dict[str, Any],
    canonical_request: Optional[CloningRequest],
    session_id: str,
) -> WorkflowInput:
    """Build a WorkflowInput for plasmid_design (describe-plasmid handler).

    The describe-plasmid handler runs entirely on a `PlasmidSpec` — there are
    no parts, no vector, no inventory. The adapter builds the spec via
    `plasmid_spec.build_plasmid_spec(message=, intent_result=)` and surfaces
    it three ways so downstream consumers can pick whichever shape they need:

      - `workflow_input.plasmid_spec`          → the spec's `.to_dict()` form
                                                 (canonical for serialization)
      - `provenance["plasmid_spec_object"]`    → the live PlasmidSpec object
                                                 (handler reads attributes)
      - `provenance["message"]`                → the raw user prompt
      - `provenance["intent_result"]`          → the intent-LLM parse

    Parts / vector / target stay None — there are no uploads in this flow.
    """
    from .plasmid_spec import build_plasmid_spec

    spec = build_plasmid_spec(message=message, intent_result=intent_result)

    constraints = list((canonical_request.constraints if canonical_request else []) or [])
    notes = list((canonical_request.normalizer_notes if canonical_request else []) or [])

    return WorkflowInput(
        intent="plasmid_design",
        workflow="plasmid_design",
        session_id=session_id,
        target=None,
        parts=[],
        vector=None,
        plasmid_spec=spec.to_dict(),
        spec_diff=None,
        constraints=constraints,
        assembly=AssemblyPlan(method="plasmid_design", topology="circular"),
        workflow_args={},
        provenance={
            "plasmid_spec_object": spec,
            "message": message,
            "intent_result": intent_result,
            "normalizer_notes": notes,
        },
    )


# ---------------------------------------------------------------------------
# sdm_design — single uploaded plasmid + feature/codon/seq mutation request
# ---------------------------------------------------------------------------

async def build_for_sdm_design(
    *,
    file: Optional[UploadFile],
    inventory_files: Optional[List[UploadFile]],
    has_target: bool,
    has_inventory: bool,
    intent_result: Dict[str, Any],
    canonical_request: Optional[CloningRequest],
    session_id: str,
) -> WorkflowInput:
    """Build a WorkflowInput for sdm_design.

    Resolves the plasmid to mutate: target file preferred, first inventory
    file as fallback. Parses the GenBank for sequence + features. The handler
    runs SDMMutationParser, conditional annotation augmentation, and the
    SDMOperator on the resolved target.

    Outputs:
      - workflow_input.target → PlasmidView with sequence, gb_text,
        annotations=features (the per-feature dict shape the SDM parser
        consumes via plasmid_features arg).
      - workflow_input.target.metadata: optional {"parse_error": str},
        {"record_name": str} when parse succeeded.
      - workflow_input.target = None when no file was uploaded.
      - workflow_input.workflow_args["sdm_params"] = intent_result["sdm"]
        + canonical_request.sdm_params (canonical wins on collision).
    """
    sdm_file = None
    if has_target and file is not None:
        sdm_file = file
    elif has_inventory and inventory_files:
        sdm_file = inventory_files[0]

    target_view: Optional[PlasmidView] = None
    if sdm_file is not None:
        try:
            await sdm_file.seek(0)
        except Exception:
            pass
        raw = await sdm_file.read()
        try:
            file_text = raw.decode("utf-8")
        except UnicodeDecodeError:
            file_text = raw.decode("latin-1", errors="replace")

        plasmid_seq = ""
        features: List[Dict[str, Any]] = []
        record_name: Optional[str] = None
        parse_error: Optional[str] = None
        try:
            from Bio import SeqIO
            from io import StringIO
            record = SeqIO.read(StringIO(file_text), "genbank")
            plasmid_seq = str(record.seq).upper()
            record_name = record.name or None
            for feat in record.features:
                feat_dict = {
                    "type": feat.type,
                    "start": int(feat.location.start),
                    "end": int(feat.location.end),
                    "strand": feat.location.strand,
                    "name": (
                        feat.qualifiers.get("label",
                            feat.qualifiers.get("gene", [""]))[0]
                        if feat.qualifiers.get("label") or feat.qualifiers.get("gene")
                        else ""
                    ),
                }
                features.append(feat_dict)
        except Exception as exc:
            parse_error = str(exc)

        meta: Dict[str, Any] = {}
        if parse_error:
            meta["parse_error"] = parse_error
        if record_name:
            meta["record_name"] = record_name

        target_view = PlasmidView(
            name=record_name or "plasmid",
            sequence=plasmid_seq,
            role="target",
            source="uploaded",
            source_file=getattr(sdm_file, "filename", None),
            topology="circular",
            gb_text=file_text,
            annotations=features,
            metadata=meta,
        )

    sdm_params = dict(intent_result.get("sdm") or {})
    if canonical_request and canonical_request.sdm_params:
        cr = canonical_request.sdm_params or {}
        sdm_params.update(cr)
    workflow_args: Dict[str, Any] = {"sdm_params": sdm_params}

    constraints = list((canonical_request.constraints if canonical_request else []) or [])
    notes = list((canonical_request.normalizer_notes if canonical_request else []) or [])

    return WorkflowInput(
        intent="sdm_design",
        workflow="sdm_design",
        session_id=session_id,
        target=target_view,
        parts=[],
        vector=None,
        constraints=constraints,
        assembly=AssemblyPlan(method="sdm", topology="circular"),
        workflow_args=workflow_args,
        provenance={"normalizer_notes": notes},
    )


# ---------------------------------------------------------------------------
# gateway_cloning — donor + insert classification via Step-2.75 cloning features
# ---------------------------------------------------------------------------

def _extract_gateway_op_features(record) -> List[Dict[str, Any]]:
    """Extract BioPython features in the GatewayOperator-shape dict.
    Mirrors the legacy chat.py `extract_features_from_record` closure."""
    features: List[Dict[str, Any]] = []
    keep_types = ("CDS", "promoter", "terminator", "misc_feature",
                  "rep_origin", "protein_bind")
    for feat in record.features:
        if feat.type not in keep_types:
            continue
        name = feat.qualifiers.get(
            "label",
            feat.qualifiers.get("gene", feat.qualifiers.get("product", [""])),
        )
        if isinstance(name, list):
            name = name[0] if name else feat.type
        color = feat.qualifiers.get("ApEinfo_fwdcolor", ["#999999"])
        if isinstance(color, list):
            color = color[0]
        strand = feat.location.strand if hasattr(feat.location, "strand") else None
        direction = 1 if strand == 1 else (-1 if strand == -1 else 1)
        features.append({
            "name": name,
            "type": feat.type,
            "start": int(feat.location.start),
            "end": int(feat.location.end),
            "direction": direction,
            "color": color,
        })
    return features


def _att_subtypes_from_cloning_features(cloning: Any) -> List[str]:
    """Pull gateway_att subtype labels from the cloning_features payload
    emitted by Step 2.75 of the annotation pipeline. The pipeline already
    runs scan_att_sites internally so the gateway handler can read these
    instead of re-scanning each plasmid."""
    if not cloning:
        return []
    # cloning may be a dict (annotate_cached return) or a ScanResult-like
    # object — both expose `features` with `feature_family` + `subtype`/`name`.
    if isinstance(cloning, dict):
        features_iter = cloning.get("features") or []
    else:
        features_iter = getattr(cloning, "features", None) or []

    subtypes: List[str] = []
    for f in features_iter:
        family = (
            f.get("feature_family") if isinstance(f, dict)
            else getattr(f, "feature_family", None)
        )
        if family != "gateway_att":
            continue
        subtype = (
            f.get("subtype") if isinstance(f, dict)
            else getattr(f, "subtype", None)
        )
        if not subtype:
            subtype = (
                f.get("name") if isinstance(f, dict)
                else getattr(f, "name", None)
            )
        if subtype:
            subtypes.append(subtype)
    return subtypes


_GATEWAY_INSERT_PATTERNS = (
    r"insertion of ([A-Za-z0-9_\-]+)",
    r"clone ([A-Za-z0-9_\-]+) into",
    r"insert ([A-Za-z0-9_\-]+)",
    r"design primers for ([A-Za-z0-9_\-]+)",
)
_GATEWAY_INSERT_STOPWORDS = frozenset({
    "this", "that", "the", "into", "via", "using", "with", "vector",
    "plasmid", "primers", "insertion", "insert", "cloning", "clone",
    "design", "primer", "gene", "dna", "sequence",
})


async def build_for_gateway_cloning(
    *,
    message: str,
    file: Optional[UploadFile],
    inventory_files: Optional[List[UploadFile]],
    intent_result: Dict[str, Any],
    canonical_request: Optional[CloningRequest],
    session_id: str,
) -> WorkflowInput:
    """Build a WorkflowInput for gateway_cloning.

    Per-plasmid annotation runs once via annotate_cached(depth="full"). Step
    2.75 (cloning_feature_annotator.scan_cloning_features) is part of that
    pipeline and already produces the gateway_att hits. The adapter projects
    those subtypes onto each part's metadata["att_subtypes"] so the handler
    never has to call scan_att_sites again.

    Role classification matches the legacy block:
      - attP present OR record name contains "pdonr"/"donor" → role="vector"
      - else → role="insert"

    When no insert is found in any uploaded file, the adapter scans the
    user message for "insertion of X" / "clone X into" / etc., resolves
    candidates against the pLannotate KB, and adds the resolved hits to
    parts and to workflow_args["modules"]. The KB audit lands on
    provenance["gateway_identified_inserts"].
    """
    from .annotation_cache import annotate_llm_cached
    from .chat import identify_features_from_kb

    parts: List[PlasmidView] = []
    modules: List[Dict[str, Any]] = []
    file_count = 0
    vector_found = False
    insert_found = False

    async def _project_upload(upload: UploadFile, idx_label: str) -> None:
        nonlocal vector_found, insert_found
        try:
            await upload.seek(0)
        except Exception:
            pass
        raw = await upload.read()
        try:
            file_text = raw.decode("utf-8")
        except UnicodeDecodeError:
            file_text = raw.decode("latin-1", errors="replace")

        from Bio import SeqIO
        from io import StringIO
        record = SeqIO.read(StringIO(file_text), "genbank")
        seq = str(record.seq)

        # Single annotation pass — Step 2.75 already produced the att hits.
        # Route through annotate_llm_cached (full LLM pipeline including
        # Step 2.75) so cloning_features + gateway_att hits are present.
        # annotate_cached(depth="full") routes to /annotate_sequence_with_hierarchy
        # which OMITS cloning_features.
        try:
            ann = await annotate_llm_cached(seq, circular=True)
        except Exception:
            ann = {}
        cloning_features = ann.get("cloning_features") or {}
        att_subtypes = _att_subtypes_from_cloning_features(cloning_features)
        site_type_prefixes = {s[:4] for s in att_subtypes}

        rec_name_lower = (record.name or "").lower()
        if "attP" in site_type_prefixes or "pdonr" in rec_name_lower or "donor" in rec_name_lower:
            role = "vector"
            vector_found = True
        else:
            role = "insert"
            insert_found = True

        op_features = _extract_gateway_op_features(record)
        canonical_id = record.name or f"plasmid_{idx_label}"

        modules.append({
            "canonical_id": canonical_id,
            "sequence": seq,
            "role": role,
            "description": f"{record.name} ({len(seq)} bp)",
            "features": op_features,
        })
        parts.append(PlasmidView(
            name=canonical_id,
            sequence=seq,
            role=role,
            source="uploaded",
            source_file=getattr(upload, "filename", None),
            topology="circular",
            length=len(seq),
            gb_text=file_text,
            annotations=op_features,
            cloning_features=cloning_features,
            metadata={"att_subtypes": att_subtypes},
        ))

    if file is not None:
        file_count += 1
        await _project_upload(file, str(file_count))
    if inventory_files:
        for inv in inventory_files:
            file_count += 1
            await _project_upload(inv, str(file_count))

    # KB-resolved inserts from message — only when no insert came from files.
    gateway_identified_inserts: List[Dict[str, Any]] = []
    if not insert_found:
        potential: List[str] = []
        for pat in _GATEWAY_INSERT_PATTERNS:
            potential.extend(re.findall(pat, message, re.IGNORECASE))
        potential = [f for f in potential if f.lower() not in _GATEWAY_INSERT_STOPWORDS]
        if potential:
            kb_candidates = [{"name": f, "feature_type": None} for f in potential]
            identified = identify_features_from_kb(kb_candidates)
            for feature in identified:
                if not feature.get("sequence"):
                    continue
                insert_found = True
                gateway_identified_inserts.append(feature)
                modules.append({
                    "canonical_id": feature["name"],
                    "sequence": feature["sequence"],
                    "role": "insert",
                    "description": (
                        f"{feature['name']} ({feature['length']} bp) "
                        f"[{feature['feature_id']}, {feature['feature_type']}]"
                    ),
                    "features": [],
                })
                parts.append(PlasmidView(
                    name=feature["name"],
                    sequence=feature["sequence"],
                    role="insert",
                    source="kb_lookup",
                    topology="linear",
                    metadata={
                        "att_subtypes": [],
                        **{k: v for k, v in feature.items() if k != "sequence"},
                    },
                ))

    workflow_args: Dict[str, Any] = {
        "modules": modules,
        "vector_found": vector_found,
        "insert_found": insert_found,
    }
    constraints = list((canonical_request.constraints if canonical_request else []) or [])
    notes = list((canonical_request.normalizer_notes if canonical_request else []) or [])

    provenance: Dict[str, Any] = {"normalizer_notes": notes}
    if gateway_identified_inserts:
        provenance["gateway_identified_inserts"] = gateway_identified_inserts

    return WorkflowInput(
        intent="gateway_cloning",
        workflow="gateway_cloning",
        session_id=session_id,
        target=None,
        parts=parts,
        vector=None,
        constraints=constraints,
        assembly=AssemblyPlan(method="gateway", topology="circular"),
        workflow_args=workflow_args,
        provenance=provenance,
    )


__all__ = [
    "build_for_annotate_gb",
    "build_for_gibson_design",
    "build_for_golden_gate_primer_design",
    "build_for_restriction_cloning",
    "build_for_sgrna_golden_gate",
    "build_for_plasmid_design",
    "build_for_sdm_design",
    "build_for_gateway_cloning",
]
