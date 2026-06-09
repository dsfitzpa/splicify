"""
SDM Demo endpoint for the AI Plasmid Design frontend.
Provides a pre-built example of deleting the His-tag from a test plasmid.
"""
from __future__ import annotations

import base64
import logging
from fastapi import APIRouter
from typing import Dict, Any

router = APIRouter()
logger = logging.getLogger(__name__)

# Simple test plasmid with His-tag
TEST_PLASMID_SEQ = """GACATTGATTATTGACTAGTTATTAATAGTAATCAATTACGGGGTCATTAGTTCATAGCC
CATATATGGAGTTCCGCGTTACATAACTTACGGTAAATGGCCCGCCTGGCTGACCGCCCA
ACGACCCCCGCCCATTGACGTCAATAATGACGTATGTTCCCATAGTAACGCCAATAGGG
ACTTTCCATTGACGTCAATGGGTGGAGTATTTACGGTAAACTGCCCACTTGGCAGTACAT
CAAGTGTATCATATGCCAAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCCGCC
TGGCATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTA
TTAGTCATCGCTATTACCATGGTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAG
CGGTTTGACTCACGGGGATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTT
TGGCACCAAAATCAACGGGACTTTCCAAAATGTCGTAACAACTCCGCCCCATTGACGCAA
ATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCTCGTTTAGTGAACCGT
CGCCACCATGCATCATCATCATCATCATATGAGTAAAGGAGAAGAACTATTTACCGGAGT
TGTACCAATTTTAGTTGAATTAGATGGAGATGTTAATGGACATAAATTTTCTGTAAGTGG
AGAAGGTGAAGGAGATGCAACTTATGGTAAATTAACATTAAAATTTATTTGTACAACTGG
TAAATTACCAGTTCCATGGCCAACATTAGTTACTACATTTTCTTATGGTGTTCAATGTTT
TTCAAGATATTCAGATCATATGAAACAACATGATTTTTTCAAAAGTGCAATGCCAGAAGG
TTATGTTCAAGAAAGAACAATTTTTTTCAAAGATGATGGTAATTACAAAACAAGAGCAGA
AGTTAAATTCGAAGGAGATACTTTAGTTAATAGAATTGAATTAAAAGGAATTGATTTTAA
AGAAGATGGAAATATTTTAAGTCATAAATTGGAATATAATTATAACTCTCATAATGTTTAT
ATTATGGCTGATAAACAAAAAAATGGAATTAAAGTTAATTTCAAAATTAGACATAATATT
GAAGACGGTTCTGTTCAATTAGCTGATCATTATCAACAAAATACTCCAATTGGTGATGGT
CCAGTTTTGTTACCAGATAATCATTATTTATCTACACAAAGTGCATTATCTAAAGACCCA
AATGAAAAAAGAGATCACATGGTTTTATTAGAATTTGTTACTGCTGCTGGAATTACACAT
GGAATGGATGAATTATACAAAGGATCCGGCGGCTCTGGAGGAAGC""".replace("\n", "").upper()

# The His-tag sequence at position 609-627
HIS_TAG_SEQ = "CATCATCATCATCATCAT"
HIS_TAG_START = 609
HIS_TAG_END = 627

# After deletion, the sequence
MUTATED_SEQ = TEST_PLASMID_SEQ[:HIS_TAG_START] + TEST_PLASMID_SEQ[HIS_TAG_END:]

# Primer CSV content
PRIMERS_CSV = """Primer Name,Sequence (5' to 3'),Length (bp),Tm (C),Notes
SDM_Forward,ATGAGTAAAGGAGAAGAAC,19,61.1,Order with standard desalting
SDM_Reverse,CATGGTGGCGACGGTTCA,18,69.1,Order with standard desalting
"""

# Protocol markdown content
PROTOCOL_MD = """# Q5 Site-Directed Mutagenesis Protocol

## Mutation Details
- Type: deletion
- Position: 609
- Original sequence: CATCATCATCATCATCAT (His-tag, 18 bp)
- New sequence: (none - deletion)
- Strategy: back_to_back

## Primers
- Forward: 5'-ATGAGTAAAGGAGAAGAAC-3' (Tm: 61.1°C)
- Reverse: 5'-CATGGTGGCGACGGTTCA-3' (Tm: 69.1°C)

## PCR Reaction Setup (25 µL)
| Component | Volume |
|-----------|--------|
| Q5 High-Fidelity 2X Master Mix | 12.5 µL |
| Forward Primer (10 µM) | 1.25 µL |
| Reverse Primer (10 µM) | 1.25 µL |
| Template DNA (1-25 ng) | 1 µL |
| Nuclease-free H2O | 9 µL |

## PCR Cycling
| Step | Temperature | Time | Cycles |
|------|-------------|------|--------|
| Initial Denaturation | 98°C | 30 sec | 1 |
| Denaturation | 98°C | 10 sec | 25 |
| Annealing | 68°C | 30 sec | 25 |
| Extension | 72°C | 30 sec/kb | 25 |
| Final Extension | 72°C | 2 min | 1 |
| Hold | 4°C | - | - |

## KLD Treatment
1. Mix 1 µL PCR product + 5 µL 2X KLD Reaction Buffer + 1 µL 10X KLD Enzyme Mix + 3 µL H2O
2. Incubate at room temperature for 5 minutes
3. Transform 5 µL into competent cells

## Expected Outcome
The His-tag (18 bp) will be removed from the plasmid, resulting in direct fusion
of the Kozak sequence to the GFP start codon.
"""


def encode_file(content: str) -> str:
    """Encode file content to base64."""
    return base64.b64encode(content.encode("utf-8")).decode("utf-8")


@router.get("/cloning/demo/sdm")
async def sdm_demo() -> Dict[str, Any]:
    """Return a pre-built SDM deletion demo."""
    logger.info("SDM demo endpoint called")
    
    response = {
        "ok": True,
        "reply": (
            "**Site-Directed Mutagenesis Demo: Delete His-tag**\n\n"
            "This demo shows Q5-style back-to-back primer design for deleting "
            f"the 6xHis-tag ({len(HIS_TAG_SEQ)} bp) from position {HIS_TAG_START}.\n\n"
            "**Strategy:** Back-to-back primers flank the deletion site. "
            "The forward primer anneals immediately downstream of the His-tag, "
            "and the reverse primer anneals immediately upstream.\n\n"
            "**Primers:**\n"
            "- Forward: `ATGAGTAAAGGAGAAGAAC` (Tm: 61.1°C)\n"
            "- Reverse: `CATGGTGGCGACGGTTCA` (Tm: 69.1°C)\n\n"
            "**Protocol:** Q5 PCR → KLD treatment (5 min) → Transform → Screen colonies"
        ),
        "sessionId": "demo-sdm-his-tag-deletion",
        "intent": "sdm_design",
        "viz": {
            "type": "design",  # Use 'design' type for compatibility with plasmid viewer
            "title": "SDM Demo: Delete His-tag",
            "sequence": MUTATED_SEQ,
            "topology": "circular",
            "total_length": len(MUTATED_SEQ),
            "mutation_type": "deletion",
            "primer_strategy": "back_to_back",
            "old_sequence": HIS_TAG_SEQ,
            "new_sequence": "",
            "annotations": [
                {
                    "name": "CMV promoter",
                    "start": 0,
                    "end": 600,
                    "direction": 1,
                    "color": "#10B981",
                    "type": "promoter",
                },
                {
                    "name": "Kozak",
                    "start": 600,
                    "end": 609,
                    "direction": 1,
                    "color": "#6B7280",
                    "type": "misc_feature",
                },
                {
                    "name": "Deletion site (18 bp removed)",
                    "start": 609,
                    "end": 610,
                    "direction": 0,
                    "color": "#EF4444",
                    "type": "sdm_mutation",
                },
                {
                    "name": "Fwd Primer (19 bp)",
                    "start": 609,
                    "end": 628,
                    "direction": 1,
                    "color": "#3B82F6",
                    "type": "primer",
                    "sequence": "ATGAGTAAAGGAGAAGAAC",
                    "tm": 61.1,
                },
                {
                    "name": "Rev Primer (18 bp)",
                    "start": 591,
                    "end": 609,
                    "direction": -1,
                    "color": "#8B5CF6",
                    "type": "primer",
                    "sequence": "CATGGTGGCGACGGTTCA",
                    "tm": 69.1,
                },
                {
                    "name": "GFP",
                    "start": 609,
                    "end": 1326,
                    "direction": 1,
                    "color": "#10B981",
                    "type": "CDS",
                },
            ],
        },
        "files": [
            {
                "fileName": "sdm_demo_primers.csv",
                "mimeType": "text/csv",
                "dataBase64": encode_file(PRIMERS_CSV),
            },
            {
                "fileName": "sdm_demo_protocol.md",
                "mimeType": "text/markdown",
                "dataBase64": encode_file(PROTOCOL_MD),
            },
        ],
    }
    
    logger.info(f"SDM demo response: viz type={response['viz']['type']}, files count={len(response['files'])}")
    return response
