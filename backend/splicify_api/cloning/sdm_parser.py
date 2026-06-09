"""
SDM Mutation Parser — converts natural language mutation descriptions 
into concrete old_seq/new_seq pairs using plasmid annotations.
"""
from __future__ import annotations

import logging
logger = logging.getLogger(__name__)

import os
import json
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

# Codon table for amino acid -> nucleotide conversion (most common codons)
CODON_TABLE = {
    "Ala": "GCT", "A": "GCT",
    "Arg": "CGT", "R": "CGT",
    "Asn": "AAC", "N": "AAC",
    "Asp": "GAC", "D": "GAC",
    "Cys": "TGC", "C": "TGC",
    "Gln": "CAG", "Q": "CAG",
    "Glu": "GAG", "E": "GAG",
    "Gly": "GGC", "G": "GGC",
    "His": "CAC", "H": "CAC",
    "Ile": "ATC", "I": "ATC",
    "Leu": "CTG", "L": "CTG",
    "Lys": "AAG", "K": "AAG",
    "Met": "ATG", "M": "ATG",
    "Phe": "TTC", "F": "TTC",
    "Pro": "CCG", "P": "CCG",
    "Ser": "AGC", "S": "AGC",
    "Thr": "ACC", "T": "ACC",
    "Trp": "TGG", "W": "TGG",
    "Tyr": "TAC", "Y": "TAC",
    "Val": "GTG", "V": "GTG",
    "Stop": "TAA", "*": "TAA",
}

# Common tag sequences
TAG_SEQUENCES = {
    "his-tag": "CATCATCATCATCATCAT",
    "his6": "CATCATCATCATCATCAT",
    "6xhis": "CATCATCATCATCATCAT",
    "flag": "GATTACAAGGATGACGACGATAAG",
    "flag-tag": "GATTACAAGGATGACGACGATAAG",
    "ha": "TACCCATACGATGTTCCAGATTACGCT",
    "ha-tag": "TACCCATACGATGTTCCAGATTACGCT",
    "v5": "GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG",
    "v5-tag": "GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG",
    "myc": "GAACAAAAACTCATCTCAGAAGAGGATCTG",
    "myc-tag": "GAACAAAAACTCATCTCAGAAGAGGATCTG",
    "strep": "TGGAGCCACCCGCAGTTCGAAAAA",
    "strep-tag": "TGGAGCCACCCGCAGTTCGAAAAA",
}


@dataclass
class SDMMutationSpec:
    """Parsed mutation specification ready for primer design."""
    mutation_type: str  # "deletion" | "insertion" | "substitution"
    target_start: int   # 0-indexed position in plasmid
    target_end: int     # exclusive end
    old_sequence: str   # sequence being removed/replaced
    new_sequence: str   # sequence being inserted (empty for deletion)
    confidence: float
    reasoning: str
    feature_context: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        # Validate SDM constraints
        if self.mutation_type == "insertion" and len(self.new_sequence) > 41:
            self.warnings.append(
                f"Insertion of {len(self.new_sequence)} bp exceeds 41 bp limit; "
                "consider Gibson assembly instead"
            )
        if self.mutation_type == "substitution" and len(self.new_sequence) > 41:
            self.warnings.append(
                f"Substitution with {len(self.new_sequence)} bp exceeds 41 bp limit; "
                "consider Gibson assembly instead"
            )


class SDMMutationParser:
    """Parse natural language mutation descriptions into SDMMutationSpec."""
    
    async def parse_mutation_request(
        self,
        message: str,
        plasmid_sequence: str,
        plasmid_features: List[Dict[str, Any]],
        sdm_params: Optional[Dict[str, Any]] = None,
    ) -> SDMMutationSpec:
        """
        Parse a mutation request into concrete coordinates.
        
        If the target is described relative to a feature that isn't in the
        annotation list, returns a spec with confidence=0 indicating that
        the annotation pipeline should be run first.
        """
        sdm = sdm_params or {}
        
        # Try codon-based first: "D10A", "Y66H" etc. are most specific, and the
        # intent LLM sometimes fills target_position_start alongside codon_position.
        if sdm.get("codon_position") is not None:
            return self._parse_codon_change(plasmid_sequence, plasmid_features, sdm)

        # Try direct sequence specification
        if sdm.get("old_sequence") or sdm.get("new_sequence"):
            return self._parse_direct_specification(plasmid_sequence, sdm, message)

        # Try position-based
        if sdm.get("target_position_start") is not None:
            return self._parse_position_based(plasmid_sequence, sdm)

        # Try feature-based
        if sdm.get("target_feature_name"):
            return self._parse_feature_based(message, plasmid_sequence, plasmid_features, sdm)
        
        # Fall back to LLM parsing
        return await self._parse_with_llm(message, plasmid_sequence, plasmid_features)
    
    def _find_feature(
        self,
        feature_name: str,
        features: List[Dict[str, Any]],
    ) -> Optional[Dict[str, Any]]:
        """Find a feature by name (case-insensitive, flexible matching).

        Searches the merged feature list (GenBank features + annotation-pipeline
        modules + motifs + CDS submodules) so motif-style names like 'His-tag',
        'FLAG', 'P2A', 'NLS' resolve as long as the chat dispatch ran the
        annotation pipeline before calling us.
        """
        name_lower = feature_name.lower().strip()
        logger.info(f"Looking for feature: '{feature_name}' among {len(features)} features")

        # Common motif-name aliases — map user phrasing to canonical motif
        # names emitted by feature_motifs.fna and the rule-based detector.
        _ALIASES = {
            "histag": "his", "his-tag": "his", "6xhis": "his", "hexahis": "his",
            "polyhis": "his", "his6": "his",
            "flagtag": "flag", "flag-tag": "flag",
            "hatag": "ha", "ha-tag": "ha",
            "v5tag": "v5", "v5-tag": "v5",
            "myctag": "myc", "c-myc": "myc",
            "streptag": "strep", "strep-tag": "strep",
            "p2atag": "p2a", "t2atag": "t2a", "e2atag": "e2a", "f2atag": "f2a",
        }
        alias = _ALIASES.get(name_lower) or _ALIASES.get(name_lower.replace(" ", ""))
        
        # Normalize common names (eGFP -> egfp, EGFP -> egfp, GFP -> gfp)
        name_normalized = name_lower.replace("-", "").replace("_", "").replace(" ", "")
        
        # Exact match first
        for f in features:
            f_name = (f.get("name") or f.get("label") or f.get("gene") or "").lower()
            f_normalized = f_name.replace("-", "").replace("_", "").replace(" ", "")
            if f_name == name_lower or f_normalized == name_normalized:
                logger.info(f"Found exact match: {f_name}")
                return f
        
        # Partial match (name contains feature or feature contains name)
        for f in features:
            f_name = (f.get("name") or f.get("label") or f.get("gene") or "").lower()
            f_normalized = f_name.replace("-", "").replace("_", "").replace(" ", "")
            if name_normalized in f_normalized or f_normalized in name_normalized:
                logger.info(f"Found partial match: {f_name}")
                return f
        
        # Try matching without common prefixes/suffixes
        # e.g., "eGFP" should match "GFP", "EGFP-N1", etc.
        core_name = name_normalized.lstrip("e").rstrip("tag")
        for f in features:
            f_name = (f.get("name") or f.get("label") or f.get("gene") or "").lower()
            f_core = f_name.replace("-", "").replace("_", "").replace(" ", "").lstrip("e").rstrip("tag")
            if core_name and f_core and len(core_name) >= 2 and (core_name in f_core or f_core in core_name):
                logger.info(f"Found core match: {f_name} (core: {f_core})")
                return f

        # Alias-based motif fallback: "His-tag" → match any feature whose
        # normalized name contains "his". Same for FLAG / HA / V5 / P2A etc.
        if alias:
            for f in features:
                f_name = (f.get("name") or f.get("label") or f.get("gene") or "").lower()
                f_normalized = f_name.replace("-", "").replace("_", "").replace(" ", "")
                if alias in f_normalized:
                    logger.info(f"Found alias match for {feature_name!r} ({alias}): {f_name}")
                    return f

        logger.warning(f"Feature not found: {feature_name}")
        return None
    
    def _parse_direct_specification(
        self,
        plasmid_sequence: str,
        sdm: Dict[str, Any],
        message: str,
    ) -> SDMMutationSpec:
        """Parse when old_seq and/or new_seq are directly provided."""
        old_seq = (sdm.get("old_sequence") or "").upper().replace(" ", "")
        new_seq = (sdm.get("new_sequence") or "").upper().replace(" ", "")
        
        # Check for tag names in new_seq or message
        if not new_seq or len(new_seq) < 10:
            msg_lower = message.lower()
            for tag_name, tag_seq in TAG_SEQUENCES.items():
                if tag_name in msg_lower:
                    new_seq = tag_seq
                    break
        
        # Determine mutation type
        if not old_seq and new_seq:
            mutation_type = "insertion"
        elif old_seq and not new_seq:
            mutation_type = "deletion"
        else:
            mutation_type = "substitution"
        
        # Find position in plasmid
        if old_seq:
            start = plasmid_sequence.upper().find(old_seq)
            if start < 0:
                return SDMMutationSpec(
                    mutation_type=mutation_type,
                    target_start=-1,
                    target_end=-1,
                    old_sequence=old_seq,
                    new_sequence=new_seq,
                    confidence=0.0,
                    reasoning=f"Sequence '{old_seq[:20]}...' not found in plasmid",
                )
            end = start + len(old_seq)
        else:
            # For pure insertion, need position from sdm params
            start = sdm.get("target_position_start", 0)
            end = start
        
        return SDMMutationSpec(
            mutation_type=mutation_type,
            target_start=start,
            target_end=end,
            old_sequence=old_seq,
            new_sequence=new_seq,
            confidence=0.95,
            reasoning=f"Direct specification: {mutation_type} at position {start}",
        )
    
    def _parse_position_based(
        self,
        plasmid_sequence: str,
        sdm: Dict[str, Any],
    ) -> SDMMutationSpec:
        """Parse position-based mutation."""
        raw_start = sdm.get("target_position_start")
        raw_end = sdm.get("target_position_end")
        start = int(raw_start) if raw_start is not None else 0
        end = int(raw_end) if raw_end is not None else start
        new_seq = (sdm.get("new_sequence") or "").upper()
        mutation_type = sdm.get("mutation_type", "deletion")
        
        if end > len(plasmid_sequence):
            end = len(plasmid_sequence)
        
        old_seq = plasmid_sequence[start:end].upper() if end > start else ""
        
        return SDMMutationSpec(
            mutation_type=mutation_type,
            target_start=start,
            target_end=end,
            old_sequence=old_seq,
            new_sequence=new_seq,
            confidence=0.9,
            reasoning=f"Position-based: {mutation_type} at bp {start}-{end}",
        )
    
    def _parse_codon_change(
        self,
        plasmid_sequence: str,
        features: List[Dict[str, Any]],
        sdm: Dict[str, Any],
    ) -> SDMMutationSpec:
        """Parse codon change request (e.g., 'change codon 45 from Arg to Ala' or 'Y66H in eGFP')."""
        codon_pos = int(sdm.get("codon_position", 1)) - 1  # Convert to 0-indexed
        codon_from = sdm.get("codon_from", "")
        codon_to = sdm.get("codon_to", "")
        target_cds = sdm.get("target_feature_name", "")
        
        logger.info(f"Parsing codon change: pos={codon_pos+1}, from={codon_from}, to={codon_to}, target={target_cds}")
        
        # Find the CDS feature
        cds_start = 0
        feature_found = None
        if target_cds:
            feature_found = self._find_feature(target_cds, features)
            if feature_found:
                cds_start = feature_found.get("start", feature_found.get("location", {}).get("start", 0))
                logger.info(f"CDS start position: {cds_start}")
            else:
                # Feature not found - try to find any CDS feature
                for f in features:
                    if f.get("type") == "CDS":
                        cds_start = f.get("start", 0)
                        logger.info(f"Using first CDS feature starting at: {cds_start}")
                        break
        
        # Calculate nucleotide position
        nt_start = cds_start + (codon_pos * 3)
        nt_end = nt_start + 3
        
        # Ensure we're within bounds
        if nt_end > len(plasmid_sequence):
            return SDMMutationSpec(
                mutation_type="substitution",
                target_start=nt_start,
                target_end=nt_end,
                old_sequence="",
                new_sequence="",
                confidence=0.0,
                reasoning=f"Codon position {codon_pos+1} is out of bounds (CDS starts at {cds_start}, plasmid length {len(plasmid_sequence)})",
                feature_context=target_cds or None,
                warnings=["Position out of bounds"],
            )
        
        old_codon = plasmid_sequence[nt_start:nt_end].upper()
        
        # Look up new codon - try single letter first, then 3-letter
        new_codon = CODON_TABLE.get(codon_to.upper())
        if not new_codon and len(codon_to) >= 3:
            new_codon = CODON_TABLE.get(codon_to[:3].title())
        if not new_codon:
            new_codon = "NNN"
            logger.warning(f"Could not find codon for amino acid: {codon_to}")
        
        logger.info(f"Codon change: position {nt_start}-{nt_end}, {old_codon} -> {new_codon}")
        
        return SDMMutationSpec(
            mutation_type="substitution",
            target_start=nt_start,
            target_end=nt_end,
            old_sequence=old_codon,
            new_sequence=new_codon,
            confidence=0.85 if feature_found else 0.70,
            reasoning=f"Codon change: position {codon_pos+1}, {codon_from} ({old_codon}) -> {codon_to} ({new_codon})",
            feature_context=target_cds or None,
        )
    
    def _parse_feature_based(
        self,
        message: str,
        plasmid_sequence: str,
        features: List[Dict[str, Any]],
        sdm: Dict[str, Any],
    ) -> SDMMutationSpec:
        """Parse mutation targeting a named feature."""
        feature_name = sdm.get("target_feature_name", "")
        feature = self._find_feature(feature_name, features)
        
        if not feature:
            return SDMMutationSpec(
                mutation_type=sdm.get("mutation_type", "deletion"),
                target_start=-1,
                target_end=-1,
                old_sequence="",
                new_sequence=sdm.get("new_sequence", ""),
                confidence=0.0,
                reasoning=f"Feature '{feature_name}' not found in annotations. "
                          "Run annotation pipeline to identify feature location.",
                warnings=["NEEDS_ANNOTATION_PIPELINE"],
            )
        
        # Extract feature coordinates
        feat_start = feature.get("start", feature.get("location", {}).get("start", 0))
        feat_end = feature.get("end", feature.get("location", {}).get("end", 0))
        feat_strand = feature.get("strand", 1) or 1

        mutation_type = sdm.get("mutation_type", "deletion")
        new_seq = sdm.get("new_sequence", "")
        terminus = sdm.get("terminus")  # "N" | "C" | None

        # Handle tag insertions: look up sequence by tag name in message or description
        if mutation_type == "insertion" and (not new_seq or len(new_seq) < 4):
            search_text = (message + " " + (sdm.get("description") or "")).lower()
            for tag_name, tag_seq in TAG_SEQUENCES.items():
                if tag_name in search_text:
                    new_seq = tag_seq
                    break

        if mutation_type == "insertion":
            # Terminal insertion: place between ATG/stop codon and flanking sequence.
            # Convention: N-terminus -> immediately AFTER start codon (pos feat_start+3 for +strand)
            #             C-terminus -> immediately BEFORE stop codon (pos feat_end-3 for +strand)
            # For reverse-strand CDS, swap the offsets.
            if terminus == "C":
                if feat_strand == -1:
                    pos = feat_start + 3
                else:
                    pos = feat_end - 3
                reasoning = f"Insertion at C-terminus of '{feature_name}' (pos {pos}, just before stop codon)"
            elif terminus == "N":
                if feat_strand == -1:
                    pos = feat_end - 3
                else:
                    pos = feat_start + 3
                reasoning = f"Insertion at N-terminus of '{feature_name}' (pos {pos}, just after start codon)"
            else:
                # No terminus specified: default to immediately before the feature start
                pos = feat_start
                reasoning = f"Insertion at start of '{feature_name}' (pos {pos}); specify terminus=N|C for terminal fusion"

            return SDMMutationSpec(
                mutation_type="insertion",
                target_start=pos,
                target_end=pos,
                old_sequence="",
                new_sequence=new_seq,
                confidence=0.9 if terminus else 0.7,
                reasoning=reasoning,
                feature_context=feature_name,
            )

        # Non-insertion (deletion / substitution): span the whole feature
        old_seq = plasmid_sequence[feat_start:feat_end].upper()
        return SDMMutationSpec(
            mutation_type=mutation_type,
            target_start=feat_start,
            target_end=feat_end,
            old_sequence=old_seq,
            new_sequence=new_seq,
            confidence=0.9,
            reasoning=f"Targeting feature '{feature_name}' at positions {feat_start}-{feat_end}",
            feature_context=feature_name,
        )
    
    async def _parse_with_llm(
        self,
        message: str,
        plasmid_sequence: str,
        features: List[Dict[str, Any]],
    ) -> SDMMutationSpec:
        """Use LLM to parse ambiguous mutation requests."""
        api_key = os.getenv("OPENAI_API_KEY")
        if not api_key:
            return SDMMutationSpec(
                mutation_type="unknown",
                target_start=-1,
                target_end=-1,
                old_sequence="",
                new_sequence="",
                confidence=0.0,
                reasoning="Could not parse mutation request and no API key for LLM fallback",
            )
        
        try:
            from openai import AsyncOpenAI
            client = AsyncOpenAI(api_key=api_key)
            
            # Build feature list for context
            feature_list = "\n".join([
                f"- {f.get('name', f.get('label', 'unknown'))}: {f.get('type', 'feature')} at {f.get('start', 0)}-{f.get('end', 0)}"
                for f in features[:20]  # Limit to first 20
            ])
            
            prompt = f"""Parse this mutation request for site-directed mutagenesis.

User request: "{message}"

Available features in the plasmid:
{feature_list}

Plasmid length: {len(plasmid_sequence)} bp

Return JSON only:
{{
  "mutation_type": "deletion" | "insertion" | "substitution",
  "target_feature_name": "feature name if targeting a feature, else null",
  "target_position_start": number or null,
  "target_position_end": number or null,
  "old_sequence": "sequence to remove/replace or empty string",
  "new_sequence": "sequence to insert or empty string",
  "confidence": 0.0-1.0,
  "reasoning": "brief explanation"
}}"""
            
            response = await client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[{"role": "user", "content": prompt}],
                response_format={"type": "json_object"},
                max_tokens=500,
            )
            
            parsed = json.loads(response.choices[0].message.content)
            
            # If targeting a feature, look it up
            if parsed.get("target_feature_name"):
                feature = self._find_feature(parsed["target_feature_name"], features)
                if feature:
                    start = feature.get("start", 0)
                    end = feature.get("end", 0)
                    old_seq = plasmid_sequence[start:end].upper()
                    return SDMMutationSpec(
                        mutation_type=parsed.get("mutation_type", "deletion"),
                        target_start=start,
                        target_end=end,
                        old_sequence=old_seq if parsed.get("mutation_type") != "insertion" else "",
                        new_sequence=parsed.get("new_sequence", ""),
                        confidence=parsed.get("confidence", 0.8),
                        reasoning=parsed.get("reasoning", "LLM parsed"),
                        feature_context=parsed["target_feature_name"],
                    )
            
            return SDMMutationSpec(
                mutation_type=parsed.get("mutation_type", "unknown"),
                target_start=parsed.get("target_position_start", -1) or -1,
                target_end=parsed.get("target_position_end", -1) or -1,
                old_sequence=parsed.get("old_sequence", ""),
                new_sequence=parsed.get("new_sequence", ""),
                confidence=parsed.get("confidence", 0.5),
                reasoning=parsed.get("reasoning", "LLM parsed"),
            )
            
        except Exception as e:
            return SDMMutationSpec(
                mutation_type="unknown",
                target_start=-1,
                target_end=-1,
                old_sequence="",
                new_sequence="",
                confidence=0.0,
                reasoning=f"LLM parsing failed: {str(e)}",
            )
