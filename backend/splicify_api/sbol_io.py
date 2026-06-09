"""
SBOL3 export / import.

Converts the output of `annotate_sequence_llm` (hierarchical_annotations +
interactions + sequence) into an SBOL3 Document, and vice versa. Uses the
`sbol3` Python library.

Mapping:
  - Top-level plasmid  → sbol3.Component (type = SBO_DNA, role = SO circular region)
  - Each hierarchical_annotation whose `layer` is in {feature, module,
    submodule} → SubComponent + Feature (Range) on the parent component.
  - Every annotation's `so_role` is attached as the Component's role list.
  - Each interaction in the `interactions` payload → sbol3.Interaction with
    Participations keyed by `sbo_role`.
  - cloning_feature-layer annotations are exported as Features with
    role = SO_ENGINEERED_REGION and a custom SBOL-visual annotation key so
    viewers can distinguish them.

This exporter is a one-way converter (annotation → SBOL3); the importer is
structural (SBOL3 Component with Range features → flat list of features +
interactions) so uploaded SBOL3 files can be re-run through the annotation
pipeline.
"""

from __future__ import annotations

import io
from typing import Any, Dict, List, Optional, Tuple

import sbol3

from .so_sbo_mapping import (
    SO_ENGINEERED_REGION,
    so_role_for_feature_type,
    so_role_for_module_type,
)

# Default namespace for exported Components. Overridable per request.
DEFAULT_NAMESPACE = "https://splicify.ai/sbol/"

# SBOL3 topology URIs
SO_CIRCULAR = "http://identifiers.org/so/SO:0000988"
SO_LINEAR = "http://identifiers.org/so/SO:0000987"


def _sanitize_displayid(raw: str, idx: int) -> str:
    """SBOL3 displayIDs must match [A-Za-z_][A-Za-z0-9_]*."""
    if not raw:
        return f"feature_{idx}"
    out = []
    for ch in raw:
        if ch.isalnum() or ch == "_":
            out.append(ch)
        else:
            out.append("_")
    cleaned = "".join(out).strip("_")
    if not cleaned or not (cleaned[0].isalpha() or cleaned[0] == "_"):
        cleaned = f"f_{cleaned}" if cleaned else f"feature_{idx}"
    return cleaned[:120]


def _range_for_annotation(
    annotation: Dict[str, Any],
    parent_component: sbol3.Component,
    seq_obj: sbol3.Sequence,
    idx: int,
) -> Optional[sbol3.SequenceFeature]:
    """Create a SequenceFeature with a Range location for a single annotation."""
    start = annotation.get("start")
    end = annotation.get("end")
    if start is None or end is None:
        return None
    # SBOL3 Range is 1-based inclusive; our annotations are 0-based inclusive/exclusive
    # per internal convention. We normalize by adding 1 to start if it's zero-based.
    # We keep it simple: coerce to 1-based inclusive [start+1, end].
    start_1based = int(start) + 1
    end_1based = int(end)
    if end_1based < start_1based:
        start_1based, end_1based = end_1based, start_1based
    if end_1based < 1:
        return None

    strand = annotation.get("direction") or annotation.get("strand") or 1
    orientation = sbol3.SBOL_INLINE if strand >= 0 else sbol3.SBOL_REVERSE_COMPLEMENT

    so_role = (
        annotation.get("so_role")
        or so_role_for_module_type(annotation.get("module_type"))
        or so_role_for_feature_type(annotation.get("type"))
        or SO_ENGINEERED_REGION
    )

    disp = _sanitize_displayid(annotation.get("name") or annotation.get("module_type") or "", idx)

    rng = sbol3.Range(
        sequence=seq_obj,
        start=start_1based,
        end=end_1based,
        orientation=orientation,
    )
    sf = sbol3.SequenceFeature(
        locations=[rng],
        roles=[so_role],
        name=annotation.get("name") or annotation.get("module_type") or disp,
    )
    parent_component.features.append(sf)
    return sf


def export_annotation_to_sbol3(
    *,
    sequence: str,
    annotations: List[Dict[str, Any]],
    interactions: Optional[List[Dict[str, Any]]] = None,
    plasmid_name: str = "plasmid",
    circular: bool = True,
    namespace: str = DEFAULT_NAMESPACE,
) -> sbol3.Document:
    """Build an SBOL3 Document from an annotation pipeline output.

    Args:
        sequence: raw DNA string (ACGT / IUPAC).
        annotations: hierarchical_annotations list (after enrichment).
        interactions: the `interactions` payload from the annotation response.
        plasmid_name: used as displayId + name on the top-level Component.
        circular: True → SO_CIRCULAR, False → SO_LINEAR topology role.
        namespace: SBOL3 namespace for generated URIs.

    Returns:
        An `sbol3.Document` containing one top-level Component, one Sequence,
        plus SubComponent + Feature entries for each annotation and
        Interaction entries for each functional relationship.
    """
    sbol3.set_namespace(namespace)
    doc = sbol3.Document()

    disp = _sanitize_displayid(plasmid_name, 0)
    top = sbol3.Component(
        identity=disp,
        types=[sbol3.SBO_DNA],
        roles=[SO_CIRCULAR if circular else SO_LINEAR],
        name=plasmid_name,
    )
    doc.add(top)

    # Sequence object
    seq = sbol3.Sequence(
        identity=f"{disp}_seq",
        elements=(sequence or "").upper(),
        encoding=sbol3.IUPAC_DNA_ENCODING,
    )
    doc.add(seq)
    top.sequences = [seq]

    # Features (SequenceFeature objects on the top component)
    participant_index: Dict[Tuple[int, int, str], sbol3.SequenceFeature] = {}
    for i, ann in enumerate(annotations or []):
        sf = _range_for_annotation(ann, top, seq, i)
        if sf is None:
            continue
        key = (int(ann.get("start") or 0), int(ann.get("end") or 0), (ann.get("name") or "").lower())
        participant_index[key] = sf

    # Interactions
    for j, ix in enumerate(interactions or []):
        ix_type = ix.get("sbo_term") or sbol3.SBO_STIMULATION
        ix_id = _sanitize_displayid(ix.get("interaction_id") or f"interaction_{j}", j)

        participations: List[sbol3.Participation] = []
        for p in ix.get("participants", []):
            role = p.get("sbo_role")
            if not role:
                continue
            key = (int(p.get("start") or 0), int(p.get("end") or 0), (p.get("name") or "").lower())
            sf = participant_index.get(key)
            if sf is None:
                # External (e.g. Cre recombinase) or no matching feature — skip
                continue
            participations.append(sbol3.Participation(roles=[role], participant=sf))

        if not participations:
            continue

        interaction = sbol3.Interaction(
            types=[ix_type],
            participations=participations,
            name=ix.get("notes") or ix_id,
        )
        # attach to top component
        top.interactions.append(interaction)

    return doc


def document_to_string(doc: sbol3.Document, file_format: str = sbol3.SORTED_NTRIPLES) -> str:
    """Serialize an SBOL3 Document to a string."""
    buf = io.StringIO()
    doc.write_string(file_format=file_format)  # returns str in sbol3 >=1.1
    # sbol3.Document.write_string returns the serialization directly in recent versions
    try:
        return doc.write_string(file_format=file_format)
    except TypeError:
        # fallback older API
        doc.write(buf, file_format=file_format)
        return buf.getvalue()


def import_sbol3(source: str, *, file_format: Optional[str] = None) -> Dict[str, Any]:
    """Parse an SBOL3 document (string) into a flat annotation payload.

    Returns a dict with:
        - plasmid_name
        - sequence
        - annotations: list of flat feature dicts
        - interactions: list of interaction dicts (name + participants)

    Only reads Range locations; more complex Location types are skipped.
    """
    doc = sbol3.Document()
    if file_format is None:
        # Try to auto-detect — sbol3 accepts ttl, rdf/xml, ntriples, sorted-ntriples
        for fmt in (sbol3.SORTED_NTRIPLES, sbol3.TURTLE, sbol3.RDF_XML, sbol3.JSONLD):
            try:
                doc.read_string(source, file_format=fmt)
                break
            except Exception:
                continue
    else:
        doc.read_string(source, file_format=file_format)

    # Find the first top-level Component with a sequence
    top: Optional[sbol3.Component] = None
    for obj in doc.objects:
        if isinstance(obj, sbol3.Component):
            top = obj
            break
    if top is None:
        return {"plasmid_name": "", "sequence": "", "annotations": [], "interactions": []}

    sequence = ""
    for seq_ref in list(top.sequences):
        # In sbol3 >=1.0 Component.sequences is a list of ReferencedURIs
        seq_uri = str(seq_ref)
        seq_obj = doc.find(seq_uri)
        if isinstance(seq_obj, sbol3.Sequence) and seq_obj.elements:
            sequence = seq_obj.elements
            break

    annotations: List[Dict[str, Any]] = []
    for sf in top.features:
        if not isinstance(sf, sbol3.SequenceFeature):
            continue
        for loc in sf.locations:
            if not isinstance(loc, sbol3.Range):
                continue
            strand = 1 if loc.orientation != sbol3.SBOL_REVERSE_COMPLEMENT else -1
            annotations.append({
                "name": sf.name or "",
                "start": int(loc.start) - 1,
                "end": int(loc.end),
                "direction": strand,
                "strand": strand,
                "so_role": sf.roles[0] if sf.roles else SO_ENGINEERED_REGION,
                "source": "sbol3_import",
            })

    interactions: List[Dict[str, Any]] = []
    for ix in top.interactions:
        parts = []
        for part in ix.participations:
            participant_uri = str(part.participant)
            referent = doc.find(participant_uri)
            role = part.roles[0] if part.roles else None
            parts.append({
                "name": getattr(referent, "name", "") or "",
                "role": role,
                "sbo_role": role,
            })
        interactions.append({
            "interaction_type": ix.types[0] if ix.types else None,
            "sbo_term": ix.types[0] if ix.types else None,
            "participants": parts,
            "notes": getattr(ix, "name", "") or "",
        })

    return {
        "plasmid_name": getattr(top, "name", None) or str(top.identity),
        "sequence": sequence,
        "annotations": annotations,
        "interactions": interactions,
    }


__all__ = [
    "export_annotation_to_sbol3",
    "document_to_string",
    "import_sbol3",
    "DEFAULT_NAMESPACE",
]
