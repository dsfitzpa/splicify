"""
LLM-based Module Parser for Plasmid Annotation - ENHANCED VERSION

Incorporates comprehensive heuristics from expression_module_heuristics.csv
and module_heuristics.csv for improved module identification.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path
import json
import os

# Enhanced module identification system prompt with comprehensive heuristics
MODULE_IDENTIFICATION_PROMPT = """You are a molecular biology expert specializing in plasmid functional module identification. Given annotated features from a plasmid sequence, identify all functional modules, their hierarchical relationships, and boundaries.

## Core Expression Module Types

### MAMMALIAN POL II EXPRESSION
**5' Boundaries (Promoters)**:
- Strong: CMV, CAG (CMV enhancer + chicken β-actin + intron), EF-1α + intron
- Moderate: PGK, UbC, SFFV, MSCV, RSV
- Weak: SV40, HSV TK
- Tissue-specific: hSyn, CaMKII (neuronal)
- Inducible: TRE3G (Tet-responsive)

**Enhancers/Introns** (between promoter and CDS):
- CMV enhancer (preceding non-CMV promoter)
- Chimeric/hybrid/β-globin intron
- EF-1α intron A
- Kozak sequence (at CDS start)

**Internal Elements**:
- WPRE (after CDS, before polyA) - boosts mRNA
- IRES / 2A peptides (P2A/T2A/E2A/F2A) - bicistronic
- Signal peptides: Igκ leader, tPA, hGH (secreted protein)

**3' Boundaries (polyA)**:
- bGH poly(A), SV40 poly(A), hGH poly(A), β-globin poly(A)

**Canonical Order**: [Promoter] → [optional intron] → [Kozak] → CDS → [optional WPRE] → [polyA]

### MAMMALIAN POL III (sgRNA/shRNA)
**Promoters**: U6, H1, 7SK (Pol III, for guide RNAs)
**Content**: 20bp spacer + gRNA scaffold/tracrRNA
**Terminator**: 6+ T run (TTTTTT)
**Exclusions**: NO Kozak, IRES, or polyA signals

### BACTERIAL CDS EXPRESSION
**Promoters**: T7, T5, tac, trc, lac, araBAD, tet
**Required**: RBS (Shine-Dalgarno) 5-15bp upstream of ATG
**Operators**: lac operator, tet operator (near promoter)
**Tags**: pelB (periplasmic), MBP, NusA, SUMO (solubility)
**Terminators**: T7 terminator, rrnB T1/T2, lambda t0

**Canonical Order**: [Promoter] → [operators] → RBS → CDS → [Terminator]

### INSECT CELL EXPRESSION
**Baculoviral promoters**: polyhedrin (polh), p10, p6.9, gp64, IE1
**Stable insect**: OpIE-1, OpIE-2 (no viral infection needed), Ac5
**Enhancers**: hr5 enhancer (with IE1)
**polyA**: SV40 poly(A), OpIE-2 poly(A), hsp70 poly(A)

### PLANT T-DNA EXPRESSION
**Promoters**: CaMV 35S (dicots), Ubi (monocots), Act1, NOS
**Enhancers**: TMV Ω leader, BYDV 5' UTR
**Terminators**: NOS terminator, CaMV poly(A), OCS terminator
**Boundaries**: Must be between LB and RB T-DNA repeats

### YEAST CDS EXPRESSION
**S. cerevisiae promoters**: GAL1/GAL10 (inducible), ADH1, TEF, GAP
**S. pombe**: nmt1 (thiamine-repressible)
**Pichia**: AOX1 (methanol-inducible)
**Signal**: α-factor secretion signal (secreted proteins)
**Terminators**: CYC1, ADH1, AOX1

### IVT mRNA (LNP PAYLOAD)
**Order**: T7 promoter → β-globin 5' UTR → Kozak → CDS → β-globin 3' UTR → polyA tract (80-120 A)


## Replication & Maintenance Modules

### E. COLI REPLICATION
**High-copy (500-700 copies)**: pUC ori, ColE1 ori (NO rop protein)
**Low-copy**: pSC101 ori (~5), p15A ori (~10), CloDF13 (~20-40)
**Conditional**: R6K γ ori (requires pir+ host)

### BROAD HOST-RANGE
**Replicators**: pBBR1, RSF1010 (IncQ), pVS1 (Agrobacterium), RK2/RP4 (IncP)

### BAC/F REPLICATION
**Module**: repE + sopA/B/C (mini-F replicon, single-copy)

### YEAST REPLICATION
**High-copy**: 2μ ori (20-40 copies, S. cerevisiae)
**Low-copy**: CEN/ARS (1-2 copies, chromosomal behavior)

### MAMMALIAN EPISOMAL
**SV40 ori**: Requires T-antigen (COS cells)
**EBV oriP + EBNA1**: Stable episomal maintenance in human cells


## Selection Markers

### BACTERIAL ANTIBIOTIC RESISTANCE
**Markers**: AmpR (β-lactamase), KanR/NeoR, CmR, TcR, SmR/SpR, GmR, HygR, ZeoR/BleoR
**Promoters**: Native resistance promoters (bla promoter, cat promoter)

### MAMMALIAN DRUG SELECTION
**Markers**: PuroR, NeoR/KanR (G418), HygR, BSD (blasticidin)
**Promoters**: PGK, SV40, TK, EF-1α core, EM7 (bifunctional E. coli + mammalian)
**polyA**: SV40, PGK, HSV TK poly(A)

### YEAST AUXOTROPHIC
**Markers**: URA3, LEU2, HIS3, TRP1, LYS2, ADE2

### COUNTER-SELECTION
**Markers**: ccdB (Gateway), SacB, DTA, barnase, URA3 (+ 5-FOA)


## Recombination Systems

### CRE/LOXP
**Sites**: loxP (flanking cargo = floxed)
**Variants**: lox2272, loxN, lox5171 (mutually incompatible)
**Usage**: Two direct repeats → excision; inverted → inversion

### FLP/FRT
**Sites**: FRT (48bp, flanking cargo)
**Usage**: Yeast-derived, similar to Cre/loxP

### GATEWAY (λ att)
**Sites**: attB/L/R/P + number (attL1/attL2, attR1/attR2)
**Destination vector**: attR1 → ccdB → CmR → attR2
**MultiSite**: att3-att6 for multi-fragment assembly

### INTEGRASE-MEDIATED
**Sites**: λ attP/attB, φC31 attP, HK022, phi80
**Usage**: Single-copy genomic integration


## Transposons

### SLEEPING BEAUTY
**Sites**: SB ITR (L) and SB ITR (R) flanking cargo
**Transposase**: SB100X (can be on same or separate plasmid)

### PIGGYBAC
**Sites**: PB ITR 5' and 3' TRD
**Target**: TTAA sites, seamless excision

### Tn5/Tn7
**Tn5**: Mosaic Ends (ME) in inverted orientation
**Tn7**: Tn7L and Tn7R (<15kb apart) - Bac-to-Bac signature


## Regulatory Systems

### TET-ON/OFF
**Regulators**: tTA (Tet-Off), rtTA (Tet-On), TetR, tTS
**Target**: TRE-containing promoters

### GAL4/UAS
**Regulator**: Gal4 (tissue-specific promoter)
**Target**: UAS (5X UAS → minimal promoter)
**Repressor**: Gal80/Gal80ts

### AUXIN-INDUCIBLE DEGRADATION (AID)
**Components**: OsTIR1 + AID/mini-AID tags
**Usage**: Reversible protein depletion

### DIMERIZATION
**Components**: FKBP + FRB (rapamycin-inducible)


## CRISPR Components

### CAS NUCLEASES
**SpCas9 variants**: Cas9, dCas9, HF1/2/4, eSpCas9, HypaCas9, xCas9, nickases (D10A/H840A)
**Other species**: SaCas9, St1Cas9, NmCas9, CjCas9
**Cas12**: LbCpf1, AsCpf1 (T-rich PAMs)
**Cas13**: LwCas13a, RfxCas13d (RNA targeting)

### CRISPR FUSION EFFECTORS
**CRISPRa**: dCas9-VP64/VP48/p65/Rta/p300
**CRISPRi**: dCas9-KRAB/Mxi1/DNMT3A
**Base editors**:
  - CBE: Cas9n-APOBEC/PmCDA1-UGI (C→T)
  - ABE: Cas9n-TadA* (A→G)


## Insulators & Structural Elements

**Insulators**: cHS4 (chicken β-globin HS4), gypsy insulator
**Usage**: Flank expression modules to block enhancer interference

**MCS/Polylinker & Blue-White Screening**:
- MCS: Multiple restriction sites for cloning
- **lacZα detection** (CRITICAL for blue-white screening):
  - lacZα (or "lacZ a") is the α-complementation fragment of β-galactosidase
  - Enables blue/white colony screening when MCS is within lacZα
  - Works with lacZΔM15 hosts (DH5α, XL1-Blue, TOP10)
  - **Intact lacZα** (no insert): Blue colonies on X-gal plates
  - **Disrupted lacZα** (insert present): White colonies
  - Check if lacZα reading frame is intact or disrupted
  - Module type: `lacz_alpha_screening` (functional) or `lacz_alpha_disrupted` (insert present)
- Sequencing primers: M13 fwd/rev, T7, SP6 flank MCS


## Module Hierarchy

**DO identify** (high-level architecture):
0. **Payload Expression Level** (highest): mammalian_lentiviral_expression, mammalian_aav_expression, bacterial_expression, insect_expression, plant_expression, yeast_expression, ivt_mrna_expression
1. **Payload Level**: lentiviral_payload, aav_payload, bacterial_backbone
2. **Expression Cassettes**: pol2_expression_cassette, pol3_expression_cassette, guide_expression_cassette
3. **Selection Cassettes**: bacterial_selection, mammalian_selection, yeast_selection
4. **Replication Modules**: bacterial_ori, yeast_ori, mammalian_ori
5. **Recombination Modules**: gateway_cassette, floxed_cassette, transposon
6. **Cloning/Screening Modules**: mcs_module, lacz_alpha_screening

**DO NOT identify** (handled by CDS parser):
- ❌ cds_module, protein_module, nls_module, tag_module, linker_module, gap_module


## Payload Expression Level (Highest Hierarchy)

The payload expression level is the OUTERMOST module that encompasses the entire expression system:

**mammalian_lentiviral_expression**:
- Contains: 5' LTR → lentiviral_payload → 3' LTR
- Indicates: Packaged into lentiviral particles for mammalian cell transduction
- Key markers: LTRs, Ψ, RRE, cPPT/CTS, mammalian expression cassettes

**mammalian_aav_expression**:
- Contains: 5' ITR → aav_payload → 3' ITR
- Indicates: Packaged into AAV capsids for mammalian cell transduction
- Key markers: ITRs (145bp), payload <4.7kb, mammalian expression cassettes

**bacterial_expression**:
- Contains: Bacterial promoter/RBS → CDS → terminator + bacterial_ori + selection
- Indicates: E. coli or bacterial protein expression
- Key markers: T7/lac/araBAD promoters, RBS, bacterial terminators

**insect_expression**:
- Contains: Baculoviral/insect promoters → CDS → polyA
- Indicates: Insect cell (Sf9/High Five) or Drosophila expression
- Key markers: polh, OpIE, Ac5 promoters

**plant_expression**:
- Contains: LB → plant expression cassettes → RB
- Indicates: Agrobacterium-mediated plant transformation
- Key markers: T-DNA borders, 35S/Ubi promoters, NOS terminator

**yeast_expression**:
- Contains: Yeast promoter → CDS → terminator + yeast_ori/CEN-ARS + auxotroph
- Indicates: S. cerevisiae, S. pombe, or Pichia expression
- Key markers: GAL/ADH/AOX1 promoters, CYC1 terminator, 2μ/ARS/CEN

**ivt_mrna_expression**:
- Contains: T7 → 5' UTR → ORF → 3' UTR → polyA tract
- Indicates: In vitro transcription for mRNA therapeutics/vaccines
- Key markers: T7 promoter, globin UTRs, encoded polyA (80-120 A)


## Heuristic Rules for Boundary Detection

### Expression Cassette Boundaries
**5' boundary weight**: promoter alone (0.85-0.95), promoter + enhancer (0.95-0.97)
**3' boundary weight**: polyA signal (0.85-0.95)
**Order match bonus**: Complete canonical order (0.9-0.95)

### Selection Cassette
**Required**: Promoter + resistance CDS + terminator
**Context**: Bacterial markers often in backbone; mammalian markers in payload

### Payload vs Backbone
**Lentiviral payload**: Between 5' LTR and 3' LTR (includes Ψ, RRE, cPPT)
**AAV payload**: Between 5' ITR and 3' ITR (max ~4.7kb)
**Bacterial backbone**: Origin + selection marker + structural elements (MCS, lacZ)

### Strand Direction Rules
1. Expression cassette strand = promoter strand
2. Reverse complement features (e.g. AmpR on complement) = strand -1
3. If features conflict, use promoter or CDS strand


## Output Format

Return JSON with:
```json
{
  "modules": [
    {
      "module_id": "m001",
      "module_type": "pol2_expression_cassette",
      "start": 200,
      "end": 3800,
      "strand": 1,
      "parent_module": "payload_01",
      "nested_modules": [],
      "features": ["CMV promoter", "Cas9", "WPRE", "bGH poly(A)"],
      "metadata": {
        "promoter": "CMV",
        "has_enhancer": false,
        "has_intron": false,
        "has_wpre": true,
        "polya": "bGH",
        "expression_system": "mammalian_pol2"
      }
    },
    {
      "module_id": "backbone_01",
      "module_type": "bacterial_backbone",
      "start": 11160,
      "end": 14873,
      "strand": -1,
      "parent_module": null,
      "nested_modules": ["ori_01", "sel_01"],
      "features": ["ori", "AmpR"],
      "metadata": {
        "ori_type": "pUC/ColE1",
        "copy_number": "high",
        "selection_marker": "AmpR"
      }
    }
  ]
}
```

**CRITICAL RULES**:
1. Identify module boundaries based on heuristic weights above
2. Set strand based on primary feature (promoter for cassettes)
3. Do NOT create CDS-level modules (handled separately)
4. Include metadata to explain identification logic
5. Ensure complete coverage of plasmid sequence
"""


@dataclass
class Module:
    """Represents a functional module in a plasmid."""
    module_id: str
    module_type: str
    start: int
    end: int
    strand: int = 1
    parent_module: Optional[str] = None
    nested_modules: List[str] = field(default_factory=list)
    features: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "module_id": self.module_id,
            "module_type": self.module_type,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "parent_module": self.parent_module,
            "nested_modules": self.nested_modules,
            "features": self.features,
            "metadata": self.metadata
        }


class LLMModuleParser:
    """
    Parse plasmid functional modules using LLM-based identification.

    This parser:
    1. Takes pLannotate features as input
    2. Runs CDS submodule parsing (from grammar pipeline)
    3. Sends features + submodules to LLM for module identification
    4. Returns hierarchical modules compatible with tokenizer
    """

    def __init__(self, api_key: Optional[str] = None, model: str = "gpt-4o"):
        """
        Initialize LLM module parser.

        Args:
            api_key: OpenAI API key (defaults to OPENAI_API_KEY env var)
            model: Model to use (default: gpt-4o for better biological reasoning)
        """
        self.api_key = api_key or os.getenv("OPENAI_API_KEY")
        self.model = model

        if not self.api_key:
            raise ValueError("OpenAI API key required. Set OPENAI_API_KEY or pass api_key parameter.")

    async def parse_modules(
        self,
        sequence: str,
        plannotate_features: List[Dict[str, Any]],
        cds_submodules: Optional[List[Dict[str, Any]]] = None,
        circular: bool = True
    ) -> Dict[str, Any]:
        """
        Parse functional modules from plasmid features using LLM.

        Args:
            sequence: DNA sequence
            plannotate_features: List of features from pLannotate
            cds_submodules: Optional pre-computed CDS submodules from grammar pipeline
            circular: Whether plasmid is circular

        Returns:
            Dictionary with:
                - modules: List of Module dicts
                - hierarchical_annotations: Annotations for visualization
                - summary: Statistics
        """
        seq_len = len(sequence)

        # Build context for LLM
        context = self._build_feature_context(plannotate_features, cds_submodules, seq_len)

        # Call LLM for module identification
        llm_response = await self._call_llm_for_modules(context, seq_len, circular)

        # Parse LLM response into Module objects
        modules = self._parse_llm_response(llm_response, seq_len)

        # Build hierarchical annotations for visualization
        hierarchical_annotations = self._build_annotations(modules)

        # Generate summary
        summary = self._generate_summary(modules, plannotate_features)

        return {
            "modules": [m.to_dict() for m in modules],
            "hierarchical_annotations": hierarchical_annotations,
            "summary": summary,
            "llm_raw_response": llm_response
        }

    def _build_feature_context(
        self,
        plannotate_features: List[Dict[str, Any]],
        cds_submodules: Optional[List[Dict[str, Any]]],
        seq_len: int
    ) -> str:
        """Build context string for LLM from features."""
        lines = [f"Plasmid length: {seq_len} bp\n"]
        lines.append("## Annotated Features\n")

        # Categorize features
        promoters = []
        coding = []
        terminators = []
        viral = []
        backbone = []
        regulatory = []
        recombination = []

        for f in plannotate_features:
            name = f.get('name', '').lower()
            ftype = f.get('type', '').lower()
            start = f.get('start', 0)
            end = f.get('end', 0)

            strand = f.get("strand", 1)
            strand_str = "+" if strand == 1 else "-"
            entry = f"{f.get('name')} ({start+1}..{end}, {strand_str})"

            # Categorize
            if any(x in name for x in ['promoter', 'cmv', 'cag', 'u6', 'h1', 'ef1', 'pgk', 'gal', 'tet', 't7', 'aox']):
                promoters.append(entry)
            elif any(x in name for x in ['polya', 'bgh', 'sv40 pa', 'terminator', 'nos', 'cyc1']):
                terminators.append(entry)
            elif any(x in name for x in ['ltr', 'itr', 'psi', 'rre', 'wpre', 'cppt', 'ires']):
                viral.append(entry)
            elif any(x in name for x in ['ori', 'ampr', 'kanr', 'cmr', 'neor', 'puror', 'cole1', 'puc', 'bla', 'repe']):
                backbone.append(entry)
            elif any(x in name for x in ['loxp', 'frt', 'att', 'gateway', 'itr']):
                recombination.append(entry)
            elif any(x in name for x in ['2a', 'kozak', 'rbs', 'intron', 'enhancer', 'signal', 'utr']):
                regulatory.append(entry)
            elif ftype in ['cds', 'gene'] or any(x in name for x in ['cas9', 'gfp', 'cas12', 'cas13', 'luciferase', 'cre', 'flp']):
                coding.append(entry)
            else:
                regulatory.append(entry)

        if promoters:
            lines.append(f"**Promoters**: {', '.join(promoters)}")
        if coding:
            lines.append(f"**Coding Sequences**: {', '.join(coding)}")
        if terminators:
            lines.append(f"**Terminators**: {', '.join(terminators)}")
        if viral:
            lines.append(f"**Viral Elements**: {', '.join(viral)}")
        if recombination:
            lines.append(f"**Recombination Sites**: {', '.join(recombination)}")
        if backbone:
            lines.append(f"**Backbone/Selection**: {', '.join(backbone)}")
        if regulatory:
            lines.append(f"**Regulatory Elements**: {', '.join(regulatory)}")

        # Add CDS submodules if available
        if cds_submodules:
            lines.append(f"\n## CDS Submodules ({len(cds_submodules)} identified)\n")
            lines.append("(These are already parsed - DO NOT create overlapping cds_module/protein_module/nls_module/tag_module/linker_module)")
            for sub in cds_submodules[:10]:  # Limit to first 10
                subtype = sub.get('module_type', 'unknown')
                start = sub.get('start', 0)
                end = sub.get('end', 0)
                name = sub.get('metadata', {}).get('protein_name', subtype)
                lines.append(f"- {name} ({start}..{end}): {subtype}")

        return "\n".join(lines)

    async def _call_llm_for_modules(
        self,
        context: str,
        seq_len: int,
        circular: bool
    ) -> Dict[str, Any]:
        """Call LLM to identify modules."""
        prompt = f"""{context}

Plasmid topology: {"circular" if circular else "linear"}
Sequence length: {seq_len} bp

Using the comprehensive heuristics provided in the system prompt, identify all functional modules in this plasmid.

Return a JSON object with a "modules" array. Each module must have:
- module_id: unique ID (e.g., "payload_01", "expr_01", "ori_01")
- module_type: use specific types from the heuristics (pol2_expression_cassette, pol3_expression_cassette, bacterial_selection, mammalian_selection, bacterial_ori, lentiviral_payload, etc.)
- start: start position (0-based)
- end: end position (0-based, exclusive)
- strand: 1 (forward) or -1 (reverse) based on primary feature
- parent_module: module_id of parent (null for top-level)
- nested_modules: array of child module_ids
- features: array of feature names included in this module
- metadata: dict with detailed info (promoter type, enhancers, polyA type, ori type, copy number, selection marker, etc.)

**CRITICAL**:
1. Use heuristic weights to determine boundaries
2. Do NOT create cds_module, protein_module, nls_module, tag_module, or linker_module
3. Focus on: payloads, expression cassettes, selection cassettes, origins, recombination sites
4. Include rich metadata to explain module identification logic"""

        try:
            from openai import AsyncOpenAI
            client = AsyncOpenAI(api_key=self.api_key)

            response = await client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": MODULE_IDENTIFICATION_PROMPT},
                    {"role": "user", "content": prompt}
                ],
                response_format={"type": "json_object"},
                temperature=0.2,  # Low temp for consistent biological reasoning
                max_tokens=6000  # Increased for comprehensive module lists
            )

            result = json.loads(response.choices[0].message.content)
            return result

        except Exception as e:
            print(f"[LLM Module Parser] LLM call failed: {e}")
            return {"modules": [], "error": str(e)}

    def _parse_llm_response(
        self,
        llm_response: Dict[str, Any],
        seq_len: int
    ) -> List[Module]:
        """Parse LLM JSON response into Module objects."""
        modules = []
        raw_modules = llm_response.get("modules", [])

        for idx, m in enumerate(raw_modules):
            try:
                module = Module(
                    module_id=m.get("module_id", f"m{idx:03d}"),
                    module_type=m.get("module_type", "unknown"),
                    start=m.get("start", 0),
                    end=m.get("end", seq_len),
                    strand=m.get("strand", 1),
                    parent_module=m.get("parent_module"),
                    nested_modules=m.get("nested_modules", []),
                    features=m.get("features", []),
                    metadata=m.get("metadata", {})
                )
                modules.append(module)
            except Exception as e:
                print(f"[LLM Module Parser] Failed to parse module {idx}: {e}")

        return modules

    def _build_annotations(self, modules: List[Module]) -> List[Dict[str, Any]]:
        """Build hierarchical annotations for visualization."""
        annotations = []

        # Color map for different module types
        color_map = {
            # Payload expression level (highest)
            "mammalian_lentiviral_expression": "#6A1B9A",
            "mammalian_aav_expression": "#4A148C",
            "bacterial_expression": "#BF360C",
            "insect_expression": "#E65100",
            "plant_expression": "#1B5E20",
            "yeast_expression": "#F57F17",
            "ivt_mrna_expression": "#01579B",

            # Payload/cassette level
            "pol2_expression_cassette": "#4CAF50",
            "pol3_expression_cassette": "#2196F3",
            "guide_expression_cassette": "#00BCD4",
            "lentiviral_payload": "#9C27B0",
            "aav_payload": "#673AB7",
            "bacterial_backbone": "#795548",

            # Selection/maintenance
            "bacterial_selection": "#FF9800",
            "mammalian_selection": "#FF5722",
            "yeast_selection": "#FFC107",
            "bacterial_ori": "#607D8B",
            "yeast_ori": "#9E9E9E",
            "mammalian_ori": "#BDBDBD",

            # Recombination/cloning
            "gateway_cassette": "#E91E63",
            "floxed_cassette": "#F44336",
            "transposon": "#3F51B5",
            "mcs_module": "#78909C",
            "lacz_alpha_screening": "#4DD0E1",
            "lacz_alpha_disrupted": "#B0BEC5"
        }

        for module in modules:
            annotations.append({
                "name": module.module_type.replace("_", " ").title(),
                "start": module.start,
                "end": module.end,
                "direction": module.strand,
                "color": color_map.get(module.module_type, "#9E9E9E"),
                "layer": "module",
                "module_type": module.module_type,
                "module_id": module.module_id,
                "source": "llm_parser",
                "metadata": module.metadata,
                "parent_module": module.parent_module,
                "features": module.features
            })

        return annotations

    def _generate_summary(
        self,
        modules: List[Module],
        plannotate_features: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Generate summary statistics."""
        module_types = {}
        for m in modules:
            module_types[m.module_type] = module_types.get(m.module_type, 0) + 1

        return {
            "total_modules": len(modules),
            "module_types": module_types,
            "feature_count": len(plannotate_features)
        }


# Async helper function for direct use
async def annotate_with_llm_modules(
    sequence: str,
    plannotate_rows: List[Dict[str, Any]],
    circular: bool = True,
    cds_submodules: Optional[List[Dict[str, Any]]] = None
) -> Dict[str, Any]:
    """
    Annotate plasmid with LLM-identified modules.

    Args:
        sequence: DNA sequence
        plannotate_rows: pLannotate feature rows
        circular: Whether plasmid is circular
        cds_submodules: Optional pre-computed CDS submodules

    Returns:
        Dictionary with modules, hierarchical_annotations, summary
    """
    parser = LLMModuleParser()
    return await parser.parse_modules(sequence, plannotate_rows, cds_submodules, circular)
