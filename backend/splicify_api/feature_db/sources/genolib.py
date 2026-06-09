"""GenoLIB SBOL v1.1 parser.

Source: `labhost_All.xml` extracted from PMC4446419 (GenoLIB 2015
supplement), licensed CC-BY-4.0.

The SBOL XML has one <Collection> per organism and one <DnaComponent>
per part. The same part (displayId) may appear in multiple Collections
when an organism uses a feature from another's KB; we dedupe by
displayId and accumulate the set of organism hosts.

Each part yields a dict shaped like the legacy KB JSON's records (so
plasmid_analyzer.KnowledgeBase + downstream callers keep working):

    {
        "displayId":      "AmpR-017",
        "name":           "AmpR-017",
        "feature_type":   "CDS",                  # from SO term lookup
        "so_resource":    "so:0000316",
        "sequence":       "atg...taa",            # nt
        "length":         861,
        "description":    "confers resistance to ampicillin...",
        "hosts":          ["Aspergillus nidulans", "Escherichia Coli", ...],
        "source": {
            "upstream": "GenoLIB",
            "supplement": "PMC4446419",
            "license": "CC-BY-4.0",
            "citation": "Wilson et al. 2016, BMC Bioinformatics",
        },
    }
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable

from ..so_ontology import so_to_type


# SBOL v1.1 namespaces
_NS_SBOL = "http://sbols.org/v1#"
_NS_RDF = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
_ABOUT = f"{{{_NS_RDF}}}about"
_RESOURCE = f"{{{_NS_RDF}}}resource"
_COLLECTION_NAME_RE = re.compile(r":([^/]+)$")


def _local(tag: str) -> str:
    """Strip namespace from a {ns}Local-name tag."""
    return tag.rsplit("}", 1)[-1]


def parse_labhost(xml_path: str | Path) -> list[dict]:
    """Parse a labhost_All.xml SBOL file into a list of part dicts.

    Iterative parsing — clears nodes after read so 2.7 MB stays under a
    few MB of RAM.
    """
    xml_path = Path(xml_path)
    parts: dict[str, dict] = {}        # displayId -> record
    current_collection: str = "?"

    for event, elem in ET.iterparse(str(xml_path), events=("start", "end")):
        tag = _local(elem.tag)

        if event == "start" and tag == "Collection":
            disp = elem.find(f"{{{_NS_SBOL}}}displayId")
            if disp is not None and disp.text:
                current_collection = disp.text.strip()
            else:
                about = elem.get(_ABOUT, "")
                m = _COLLECTION_NAME_RE.search(about)
                current_collection = m.group(1) if m else "?"

        elif event == "end" and tag == "DnaComponent":
            disp_el = elem.find(f"{{{_NS_SBOL}}}displayId")
            name_el = elem.find(f"{{{_NS_SBOL}}}name")
            desc_el = elem.find(f"{{{_NS_SBOL}}}description")
            type_el = elem.find(f"{{{_NS_RDF}}}type")
            seq_el = elem.find(
                f"{{{_NS_SBOL}}}dnaSequence/"
                f"{{{_NS_SBOL}}}DnaSequence/"
                f"{{{_NS_SBOL}}}nucleotides"
            )

            display_id = (disp_el.text or "").strip() if disp_el is not None else ""
            if not display_id:
                elem.clear()
                continue

            so_resource = ""
            if type_el is not None:
                so_resource = type_el.get(_RESOURCE, "").strip()
            feature_type = so_to_type(so_resource, default="misc_feature")

            name = (name_el.text or "").strip() if name_el is not None else display_id
            description = (desc_el.text or "").strip() if desc_el is not None else ""
            sequence = (seq_el.text or "").strip().upper().replace("\n", "").replace(" ", "") if seq_el is not None else ""

            rec = parts.get(display_id)
            if rec is None:
                rec = {
                    "displayId": display_id,
                    "name": name,
                    "feature_type": feature_type,
                    "so_resource": so_resource,
                    "sequence": sequence,
                    "length": len(sequence),
                    "description": description,
                    "hosts": [current_collection],
                    "source": {
                        "upstream": "GenoLIB",
                        "supplement": "PMC4446419",
                        "license": "CC-BY-4.0",
                        "citation": "Wilson et al. 2016, BMC Bioinformatics",
                    },
                }
                parts[display_id] = rec
            else:
                if current_collection not in rec["hosts"]:
                    rec["hosts"].append(current_collection)
                # Prefer the longer description / sequence when collections disagree
                if len(description) > len(rec["description"]):
                    rec["description"] = description
                if len(sequence) > len(rec["sequence"]):
                    rec["sequence"] = sequence
                    rec["length"] = len(sequence)
            elem.clear()

        elif event == "end" and tag == "Collection":
            elem.clear()

    return list(parts.values())


def iter_genolib_parts(xml_path: str | Path) -> Iterable[dict]:
    """Streaming variant — yields parts one at a time. Same shape as
    parse_labhost, no host accumulation across Collections."""
    yield from parse_labhost(xml_path)
