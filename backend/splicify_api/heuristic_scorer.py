"""
Heuristic Scoring Engine
========================
Score module and vector calls based on heuristic rules loaded from CSV files.
"""

import csv
import re
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple, Set
from pathlib import Path

from .heuristic_motif_detector import MotifHit


# Score thresholds for filtering
MIN_CDS_SCORE = 0.70  # Minimum score to keep a CDS module
MIN_EXPRESSION_SCORE = 0.65  # Minimum score for expression modules
MAX_MODULE_COVERAGE = 0.95  # Remove modules covering >95% of plasmid


@dataclass
class HeuristicRule:
    """Parsed heuristic rule from CSV"""
    vector_type: str
    rule_id: str
    rule_type: str  # presence, presence_all, boundary, orientation, order, exclusion
    features: List[str]
    direction: str
    location_constraint: str
    weight: float
    notes: str


@dataclass
class RuleFiring:
    """Record of a rule that fired"""
    rule_id: str
    rule_type: str
    weight: float
    features_matched: List[str]
    positions: List[Tuple[int, int]]  # Start/end of matched features
    notes: str
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "rule_id": self.rule_id,
            "rule_type": self.rule_type,
            "weight": round(self.weight, 3),
            "features_matched": self.features_matched,
            "positions": self.positions,
            "notes": self.notes
        }


@dataclass 
class ModuleCall:
    """Scored module/vector call"""
    module_type: str
    start: int
    end: int
    strand: int
    score: float  # Aggregated weight from fired rules
    confidence: str  # 'high', 'medium', 'low'
    rules_fired: List[RuleFiring]
    features_included: List[str]
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "module_type": self.module_type,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "score": round(self.score, 3),
            "confidence": self.confidence,
            "rules_fired": [r.to_dict() for r in self.rules_fired],
            "features_included": self.features_included,
            "metadata": self.metadata
        }


class HeuristicScorer:
    """Score module and vector calls based on heuristic rules"""
    
    HEURISTICS_DIR = Path(__file__).parent / 'heuristics'
    
    # Feature name patterns for matching
    PROMOTER_PATTERNS = [
        r'promoter', r'CMV', r'EF-?1', r'CAG', r'PGK', r'UbC', r'SV40',
        r'T7', r'T5', r'lac', r'tac', r'araBAD', r'U6', r'H1', r'7SK'
    ]
    POLYA_PATTERNS = [
        r'poly[\s\-]?A', r'polyA', r'bGH', r'SV40.*poly', r'hGH.*poly',
        r'terminator'
    ]
    RBS_PATTERNS = [r'RBS', r'Shine[\s\-]?Dalgarno', r'ribosom']
    LTR_PATTERNS = [r'ltr', r'long\s*terminal\s*repeat']
    
    def __init__(self):
        self.cds_rules = self._load_rules('cds_module_heuristics.csv')
        self.vector_rules = self._load_rules('vector_heuristics.csv')
        self.module_rules = self._load_rules('module_heuristics.csv')
        self.expression_rules = self._load_rules('expression_module_heuristics.csv')
        
        # Build compiled patterns for feature matching
        self._promoter_re = re.compile('|'.join(self.PROMOTER_PATTERNS), re.IGNORECASE)
        self._polya_re = re.compile('|'.join(self.POLYA_PATTERNS), re.IGNORECASE)
        self._rbs_re = re.compile('|'.join(self.RBS_PATTERNS), re.IGNORECASE)
        self._ltr_re = re.compile('|'.join(self.LTR_PATTERNS), re.IGNORECASE)
    
    def _load_rules(self, filename: str) -> List[HeuristicRule]:
        """Load rules from CSV file."""
        rules = []
        filepath = self.HEURISTICS_DIR / filename
        
        if not filepath.exists():
            print(f"[HeuristicScorer] Warning: {filepath} not found")
            return rules
        
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if not row.get('rule_id'):
                        continue
                    
                    features_str = row.get('features', '')
                    features = [f.strip() for f in features_str.split('|') if f.strip()]
                    
                    weight_str = row.get('weight', '0')
                    try:
                        weight = float(weight_str)
                    except ValueError:
                        weight = 0.0
                    
                    rules.append(HeuristicRule(
                        vector_type=row.get('vector_type', ''),
                        rule_id=row.get('rule_id', ''),
                        rule_type=row.get('rule_type', 'presence'),
                        features=features,
                        direction=row.get('direction', 'forward'),
                        location_constraint=row.get('location_constraint', ''),
                        weight=weight,
                        notes=row.get('notes', '')
                    ))
            print(f"[HeuristicScorer] Loaded {len(rules)} rules from {filename}")
        except Exception as e:
            print(f"[HeuristicScorer] Error loading {filename}: {e}")
        
        return rules
    
    def _find_features_in_region(
        self,
        features: List[Dict],
        region_start: int,
        region_end: int,
        pattern: Optional[re.Pattern] = None,
        names: Optional[List[str]] = None,
        strand: Optional[int] = None
    ) -> List[Dict]:
        """Find features within a region that match criteria."""
        matches = []
        for f in features:
            f_start = f.get('start', 0)
            f_end = f.get('end', 0)
            f_name = f.get('name', '') or f.get('Feature', '') or f.get('sseqid', '')
            f_strand = f.get('strand', 1)
            
            if f_end < region_start or f_start > region_end:
                continue
            
            if strand is not None and f_strand != strand:
                continue
            
            if pattern and not pattern.search(f_name):
                continue
            
            if names:
                matched = False
                for name in names:
                    if name.lower() in f_name.lower() or f_name.lower() in name.lower():
                        matched = True
                        break
                if not matched:
                    continue
            
            matches.append(f)
        
        return matches
    
    def _find_upstream_features(
        self,
        features: List[Dict],
        position: int,
        max_distance: int,
        pattern: Optional[re.Pattern] = None,
        strand: int = 1
    ) -> List[Dict]:
        """Find features upstream of a position within max_distance."""
        if strand == 1:
            region_start = max(0, position - max_distance)
            region_end = position
        else:
            region_start = position
            region_end = position + max_distance
        
        return self._find_features_in_region(
            features, region_start, region_end,
            pattern=pattern, strand=strand
        )
    
    def _find_downstream_features(
        self,
        features: List[Dict],
        position: int,
        max_distance: int,
        pattern: Optional[re.Pattern] = None,
        strand: int = 1
    ) -> List[Dict]:
        """Find features downstream of a position within max_distance."""
        if strand == 1:
            region_start = position
            region_end = position + max_distance
        else:
            region_start = max(0, position - max_distance)
            region_end = position
        
        return self._find_features_in_region(
            features, region_start, region_end,
            pattern=pattern, strand=strand
        )
    
    def _deduplicate_overlapping_modules(
        self,
        modules: List[ModuleCall],
        overlap_threshold: float = 0.7
    ) -> List[ModuleCall]:
        """
        Remove overlapping modules on same strand, keeping highest scoring.
        
        Two modules overlap if their intersection / union > overlap_threshold.
        """
        if not modules:
            return []
        
        # Sort by score descending
        sorted_mods = sorted(modules, key=lambda m: -m.score)
        kept = []
        
        for mod in sorted_mods:
            # Check if this module significantly overlaps any already-kept module on same strand
            dominated = False
            for k in kept:
                if k.strand != mod.strand:
                    continue
                    
                # Calculate overlap
                overlap_start = max(mod.start, k.start)
                overlap_end = min(mod.end, k.end)
                
                if overlap_end <= overlap_start:
                    continue  # No overlap
                
                overlap_len = overlap_end - overlap_start
                mod_len = mod.end - mod.start
                
                # If this module is mostly contained in a higher-scoring one, skip it
                if mod_len > 0 and overlap_len / mod_len > overlap_threshold:
                    dominated = True
                    break
            
            if not dominated:
                kept.append(mod)
        
        return kept
    
    def _get_feature_name_for_orf(
        self,
        orf_start: int,
        orf_end: int,
        orf_strand: int,
        plannotate_features: List[Dict]
    ) -> Optional[str]:
        """
        Find the pLannotate feature name that best matches this ORF.
        Returns None if no good match found.
        """
        best_match = None
        best_overlap = 0
        
        for f in plannotate_features:
            f_type = (f.get('type', '') or f.get('Type', '')).lower()
            if f_type not in ('cds', 'gene', 'misc_feature'):
                continue
            
            f_start = f.get('start', 0)
            f_end = f.get('end', 0)
            f_strand = f.get('strand', 1) if 'strand' in f else (
                -1 if f.get('direction', 1) == -1 else 1
            )
            
            # Must be same strand
            if f_strand != orf_strand:
                continue
            
            # Calculate overlap
            overlap_start = max(orf_start, f_start)
            overlap_end = min(orf_end, f_end)
            
            if overlap_end > overlap_start:
                overlap_len = overlap_end - overlap_start
                orf_len = orf_end - orf_start
                
                # Good match if >70% overlap
                if orf_len > 0 and overlap_len / orf_len > 0.7:
                    if overlap_len > best_overlap:
                        best_overlap = overlap_len
                        best_match = f.get('name', '') or f.get('Feature', '')
        
        return best_match
    
    def score_cds_modules(
        self,
        motifs: Dict[str, List[MotifHit]],
        plannotate_features: List[Dict],
        sequence: str
    ) -> List[ModuleCall]:
        """
        Score potential CDS modules by evaluating rules against detected motifs and features.
        Now with deduplication and minimum score filtering.
        """
        modules = []
        seq_len = len(sequence)
        
        orfs = motifs.get('orfs', [])
        kozak_hits = motifs.get('kozak', [])
        rbs_hits = motifs.get('rbs', [])
        
        # Build index of start codon positions
        kozak_atg_positions = {
            h.metadata.get('atg_position'): h 
            for h in kozak_hits if h.metadata.get('atg_position') is not None
        }
        rbs_atg_positions = {
            h.metadata.get('atg_position'): h 
            for h in rbs_hits if h.metadata.get('atg_position') is not None
        }
        
        for orf in orfs:
            rules_fired = []
            features_included = []
            
            orf_len_aa = orf.metadata.get('length_aa', 0)
            
            # Base rule: ORF detected (CDS-01) - weight based on size
            orf_weight = 0.85 if orf_len_aa >= 100 else 0.70 if orf_len_aa >= 50 else 0.50
            rules_fired.append(RuleFiring(
                rule_id='CDS-01',
                rule_type='presence',
                weight=orf_weight,
                features_matched=[f'ORF ({orf_len_aa} aa)'],
                positions=[(orf.start, orf.end)],
                notes=f"Open reading frame detected ({orf_len_aa} amino acids)"
            ))
            
            # Check for Kozak at start
            orf_atg_pos = orf.start if orf.strand == 1 else orf.end - 3
            has_kozak = False
            if orf_atg_pos in kozak_atg_positions:
                kozak = kozak_atg_positions[orf_atg_pos]
                strength = kozak.metadata.get('strength', 'unknown')
                rule_id = {'strong': 'CDS-04', 'adequate': 'CDS-05', 'weak': 'CDS-06'}.get(strength, 'CDS-04')
                weight = {'strong': 0.95, 'adequate': 0.75, 'weak': 0.50}.get(strength, 0.5)
                
                rules_fired.append(RuleFiring(
                    rule_id=rule_id,
                    rule_type='boundary',
                    weight=weight,
                    features_matched=[f'Kozak ({strength})'],
                    positions=[(kozak.start, kozak.end)],
                    notes=f"Kozak sequence ({strength}) at translation start"
                ))
                features_included.append(f'Kozak ({strength})')
                has_kozak = True
            
            # Check for RBS at start
            has_rbs = False
            if orf_atg_pos in rbs_atg_positions:
                rbs = rbs_atg_positions[orf_atg_pos]
                strength = rbs.metadata.get('strength', 'unknown')
                rule_id = {'strong': 'CDS-07', 'adequate': 'CDS-08', 'weak': 'CDS-09'}.get(strength, 'CDS-07')
                weight = {'strong': 0.95, 'adequate': 0.75, 'weak': 0.50}.get(strength, 0.5)
                
                rules_fired.append(RuleFiring(
                    rule_id=rule_id,
                    rule_type='boundary',
                    weight=weight,
                    features_matched=[f'RBS ({strength})'],
                    positions=[(rbs.start, rbs.end)],
                    notes=f"Ribosome binding site ({strength}) upstream of start codon"
                ))
                features_included.append(f'RBS ({strength})')
                has_rbs = True
            
            # Check for pLannotate CDS/gene annotation overlap
            feature_name = self._get_feature_name_for_orf(
                orf.start, orf.end, orf.strand, plannotate_features
            )
            if feature_name:
                rules_fired.append(RuleFiring(
                    rule_id='CDS-35',
                    rule_type='presence',
                    weight=0.90,  # Boost weight for confirmed features
                    features_matched=[feature_name],
                    positions=[(orf.start, orf.end)],
                    notes=f"Matches pLannotate annotation: {feature_name}"
                ))
                features_included.append(feature_name)
            
            # Aggregate score
            score, confidence = self.aggregate_scores(rules_fired)
            
            # Filter: skip low-scoring modules
            if score < MIN_CDS_SCORE:
                continue
            
            # Use ORF boundaries directly (ATG to stop codon) - not extended to promoter/polyA
            # The ORF boundaries are the true CDS module boundaries
            module_start = orf.start
            module_end = orf.end
            
            # Filter: skip modules covering >95% of plasmid
            module_coverage = (module_end - module_start) / seq_len
            if module_coverage > MAX_MODULE_COVERAGE:
                continue
            
            modules.append(ModuleCall(
                module_type='cds_module',
                start=module_start,
                end=module_end,
                strand=orf.strand,
                score=score,
                confidence=confidence,
                rules_fired=rules_fired,
                features_included=features_included,
                metadata={
                    'orf_start': orf.start,
                    'orf_end': orf.end,
                    'orf_length_aa': orf_len_aa,
                    'feature_name': feature_name,
                    'has_kozak': has_kozak,
                    'has_rbs': has_rbs
                }
            ))
        
        # Deduplicate overlapping modules on same strand
        modules = self._deduplicate_overlapping_modules(modules, overlap_threshold=0.7)
        
        # Sort by score descending
        modules.sort(key=lambda m: -m.score)
        return modules
    
    def score_expression_modules(
        self,
        cds_modules: List[ModuleCall],
        plannotate_features: List[Dict],
        sequence: str
    ) -> List[ModuleCall]:
        """
        Score expression modules (pol2/pol3/bacterial) containing CDS modules.
        With deduplication.
        """
        modules = []
        seq_len = len(sequence)
        seen_regions = set()  # Track (start, end, strand) to avoid duplicates
        
        for cds in cds_modules:
            # Create region key for deduplication
            region_key = (cds.start, cds.end, cds.strand)
            if region_key in seen_regions:
                continue
            seen_regions.add(region_key)
            
            rules_fired = []
            module_type = 'expression_module'
            
            # Look for expression-defining features
            promoters = self._find_features_in_region(
                plannotate_features,
                cds.start, cds.end,
                pattern=self._promoter_re
            )
            
            # Classify expression type
            for prom in promoters:
                prom_name = prom.get('name', '') or prom.get('Feature', '')
                
                # Mammalian Pol II
                if re.search(r'CMV|EF-?1|CAG|PGK|UbC|SV40|SFFV|RSV', prom_name, re.IGNORECASE):
                    module_type = 'pol2_expression_mammalian'
                    rules_fired.append(RuleFiring(
                        rule_id='MPII-B5',
                        rule_type='boundary',
                        weight=0.90,
                        features_matched=[prom_name],
                        positions=[(prom.get('start', 0), prom.get('end', 0))],
                        notes="Mammalian Pol II promoter defines expression cassette"
                    ))
                    break
                
                # Mammalian Pol III
                elif re.search(r'U6|H1|7SK', prom_name, re.IGNORECASE):
                    module_type = 'pol3_expression_mammalian'
                    rules_fired.append(RuleFiring(
                        rule_id='MPIII-B5',
                        rule_type='boundary',
                        weight=0.95,
                        features_matched=[prom_name],
                        positions=[(prom.get('start', 0), prom.get('end', 0))],
                        notes="Mammalian Pol III promoter defines sgRNA/shRNA cassette"
                    ))
                    break
                
                # Bacterial
                elif re.search(r'T7|T5|lac|tac|trc|araBAD|tet', prom_name, re.IGNORECASE):
                    module_type = 'bacterial_expression'
                    rules_fired.append(RuleFiring(
                        rule_id='BAC-B5',
                        rule_type='boundary',
                        weight=0.92,
                        features_matched=[prom_name],
                        positions=[(prom.get('start', 0), prom.get('end', 0))],
                        notes="Bacterial promoter defines prokaryotic expression cassette"
                    ))
                    break
            
            if not rules_fired:
                # Unknown expression type - base it on CDS score
                rules_fired.append(RuleFiring(
                    rule_id='EXP-01',
                    rule_type='presence',
                    weight=cds.score * 0.8,
                    features_matched=['CDS module'],
                    positions=[(cds.start, cds.end)],
                    notes="Expression module inferred from CDS without clear promoter"
                ))
            
            score, confidence = self.aggregate_scores(rules_fired)
            
            # Filter: skip low-scoring expression modules
            if score < MIN_EXPRESSION_SCORE:
                continue
            
            # Filter: skip modules covering >95% of plasmid
            module_coverage = (cds.end - cds.start) / seq_len
            if module_coverage > MAX_MODULE_COVERAGE:
                continue
            
            modules.append(ModuleCall(
                module_type=module_type,
                start=cds.start,
                end=cds.end,
                strand=cds.strand,
                score=score,
                confidence=confidence,
                rules_fired=rules_fired,
                features_included=cds.features_included,
                metadata={
                    'cds_module_score': cds.score,
                    'expression_type': module_type,
                    'feature_name': cds.metadata.get('feature_name')
                }
            ))
        
        modules.sort(key=lambda m: -m.score)
        return modules
    
    def score_vector_type(
        self,
        expression_modules: List[ModuleCall],
        plannotate_features: List[Dict],
        sequence: str
    ) -> List[ModuleCall]:
        """
        Score vector type (lentiviral, AAV, transient, etc.).
        Now with proper LTR-based boundaries for lentiviral vectors.
        """
        seq_len = len(sequence)
        vector_calls = []
        
        # Build feature index by name patterns
        all_features_str = ' '.join(
            (f.get('name', '') or f.get('Feature', ''))
            for f in plannotate_features
        ).lower()
        
        # TODO: Apply boundary rules from heuristics CSV
        # Currently boundary rules (LV-01, LV-02) have location_constraint but it's not used
        ltr_5_start, ltr_3_end = None, None
        
        # Lentiviral markers
        lenti_score = 0.0
        lenti_rules = []
        
        if re.search(r'5[\'\-\s]?ltr|ltr.*5|hiv.*ltr', all_features_str):
            lenti_score += 0.95
            lenti_rules.append(RuleFiring(
                rule_id='LV-01', rule_type='boundary', weight=0.95,
                features_matched=["5' LTR"], positions=[(ltr_5_start or 0, ltr_5_start or 0)],
                notes="5' LTR marks lentiviral vector start"
            ))
        
        if re.search(r'3[\'\-\s]?ltr|ltr.*3|sin|delta.*u3', all_features_str):
            lenti_score += 0.95
            lenti_rules.append(RuleFiring(
                rule_id='LV-02', rule_type='boundary', weight=0.95,
                features_matched=["3' LTR"], positions=[(ltr_3_end or seq_len, ltr_3_end or seq_len)],
                notes="3' LTR marks lentiviral vector end"
            ))
        
        if re.search(r'psi|ψ|packaging', all_features_str):
            lenti_score += 0.90
            lenti_rules.append(RuleFiring(
                rule_id='LV-03', rule_type='presence', weight=0.90,
                features_matched=["HIV-1 Psi"], positions=[(0, 0)],
                notes="Packaging signal indicates lentiviral vector"
            ))
        
        if re.search(r'rre|rev.*response', all_features_str):
            lenti_score += 0.80
            lenti_rules.append(RuleFiring(
                rule_id='LV-04', rule_type='presence', weight=0.80,
                features_matched=["RRE"], positions=[(0, 0)],
                notes="Rev Response Element supports lentiviral identification"
            ))
        
        if re.search(r'cppt|central.*polypurine', all_features_str):
            lenti_score += 0.70
            lenti_rules.append(RuleFiring(
                rule_id='LV-05', rule_type='presence', weight=0.70,
                features_matched=["cPPT/CTS"], positions=[(0, 0)],
                notes="Central polypurine tract in lentiviral vectors"
            ))
        
        if lenti_rules:
            conf = 'high' if lenti_score >= 2.5 else 'medium' if lenti_score >= 1.5 else 'low'
            
            # Boundary detection should use rules from vector_heuristics.csv
            # Rules like LV-01 (flanks_5prime) and LV-02 (flanks_3prime) have:
            #   - features: what to look for (e.g., "5' LTR | 5' LTR (truncated)")
            #   - location_constraint: where boundary is (flanks_5prime, flanks_3prime)
            # But these aren't currently being applied to find actual positions
            lenti_start = 0
            lenti_end = seq_len
            
            vector_calls.append(ModuleCall(
                module_type='lentiviral_payload',
                start=lenti_start,
                end=lenti_end,
                strand=1,
                score=min(lenti_score / 3.0, 1.0),
                confidence=conf,
                rules_fired=lenti_rules,
                features_included=[],
                metadata={
                    'vector_class': 'lentivirus',
                    'ltr_5_start': ltr_5_start,
                    'ltr_3_end': ltr_3_end
                }
            ))
        
        # AAV markers
        aav_score = 0.0
        aav_rules = []
        
        if re.search(r'itr|inverted.*terminal', all_features_str):
            itr_count = len(re.findall(r'itr|inverted.*terminal', all_features_str))
            if itr_count >= 2:
                aav_score += 0.99
                aav_rules.append(RuleFiring(
                    rule_id='AAV-02', rule_type='presence_all', weight=0.99,
                    features_matched=["AAV2 ITR (5')", "AAV2 ITR (3')"],
                    positions=[(0, 0)],
                    notes="Two ITRs in inverted orientation diagnostic for AAV"
                ))
            else:
                aav_score += 0.98
                aav_rules.append(RuleFiring(
                    rule_id='AAV-01', rule_type='boundary', weight=0.98,
                    features_matched=["ITR"], positions=[(0, 0)],
                    notes="ITR indicates AAV vector"
                ))
        
        if lenti_rules and aav_rules:
            aav_rules.append(RuleFiring(
                rule_id='AAV-06', rule_type='exclusion', weight=-0.90,
                features_matched=["LTR/Psi elements"], positions=[(0, 0)],
                notes="LTR/Psi presence excludes AAV"
            ))
            aav_score -= 0.90
        
        if aav_rules and aav_score > 0:
            conf = 'high' if aav_score >= 0.9 else 'medium' if aav_score >= 0.5 else 'low'
            vector_calls.append(ModuleCall(
                module_type='aav_payload',
                start=0, end=seq_len, strand=1,
                score=max(0, min(aav_score, 1.0)),
                confidence=conf,
                rules_fired=aav_rules,
                features_included=[],
                metadata={'vector_class': 'aav'}
            ))
        
        # Bacterial backbone - find actual boundaries using ori and marker positions
        bac_score = 0.0
        bac_rules = []
        bac_features = []
        
        for f in plannotate_features:
            f_name = (f.get('name', '') or f.get('Feature', '')).lower()
            if re.search(r'ori|origin|cole1|pmb1|puc', f_name):
                bac_score += 0.95
                bac_features.append(f)
                bac_rules.append(RuleFiring(
                    rule_id='BR-HC-01', rule_type='boundary', weight=0.95,
                    features_matched=[f.get('name', 'ori')], 
                    positions=[(f.get('start', 0), f.get('end', 0))],
                    notes="E. coli origin of replication"
                ))
                break  # Only count once
        
        for f in plannotate_features:
            f_name = (f.get('name', '') or f.get('Feature', '')).lower()
            if re.search(r'ampr|bla|kanr|cmr|antibiotic|resistance', f_name):
                bac_score += 0.97
                bac_features.append(f)
                bac_rules.append(RuleFiring(
                    rule_id='SEL-AB-02', rule_type='presence', weight=0.97,
                    features_matched=[f.get('name', 'AmpR')],
                    positions=[(f.get('start', 0), f.get('end', 0))],
                    notes="Bacterial selection marker"
                ))
                break  # Only count once
        
        if bac_rules:
            # Calculate bacterial backbone boundaries from feature positions
            if bac_features:
                bac_start = min(f.get('start', 0) for f in bac_features)
                bac_end = max(f.get('end', seq_len) for f in bac_features)
            else:
                bac_start = 0
                bac_end = seq_len
            
            # Don't add bacterial backbone if it covers >95% of plasmid
            bac_coverage = (bac_end - bac_start) / seq_len
            if bac_coverage <= MAX_MODULE_COVERAGE:
                conf = 'high' if bac_score >= 1.5 else 'medium' if bac_score >= 0.9 else 'low'
                vector_calls.append(ModuleCall(
                    module_type='bacterial_backbone',
                    start=bac_start,
                    end=bac_end,
                    strand=1,
                    score=min(bac_score / 2.0, 1.0),
                    confidence=conf,
                    rules_fired=bac_rules,
                    features_included=[],
                    metadata={'vector_class': 'bacterial'}
                ))
        
        # If no viral elements, default to transient pDNA - but don't create full-plasmid annotation
        if not lenti_rules and not aav_rules and not vector_calls:
            vector_calls.append(ModuleCall(
                module_type='transient_pdna',
                start=0, end=seq_len, strand=1,
                score=0.6,
                confidence='medium',
                rules_fired=[RuleFiring(
                    rule_id='PD-03', rule_type='exclusion', weight=0.6,
                    features_matched=[], positions=[],
                    notes="No viral packaging elements - classified as transient pDNA"
                )],
                features_included=[],
                metadata={'vector_class': 'transient'}
            ))
        
        vector_calls.sort(key=lambda m: -m.score)
        return vector_calls
    
    def aggregate_scores(self, firings: List[RuleFiring]) -> Tuple[float, str]:
        """
        Aggregate rule firings into final score and confidence.
        """
        if not firings:
            return 0.0, 'low'
        
        total_weight = sum(f.weight for f in firings)
        max_possible = sum(f.weight for f in firings if f.weight > 0)
        
        if max_possible > 0:
            score = max(0, min(1, total_weight / max(max_possible, 1.5)))
        else:
            score = max(0, min(1, total_weight))
        
        if score >= 0.85:
            confidence = 'high'
        elif score >= 0.60:
            confidence = 'medium'
        else:
            confidence = 'low'
        
        return score, confidence
    
    def rank_modules(
        self,
        cds_modules: List[ModuleCall],
        expression_modules: List[ModuleCall],
        vector_calls: List[ModuleCall]
    ) -> Dict[str, List[ModuleCall]]:
        """
        Return ranked module calls organized by type.
        """
        return {
            'cds_modules': sorted(cds_modules, key=lambda m: -m.score),
            'expression_modules': sorted(expression_modules, key=lambda m: -m.score),
            'vector_type': sorted(vector_calls, key=lambda m: -m.score)
        }
