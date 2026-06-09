"""
Grammar-Based Mammalian Pol II Expression Cassette Detector
============================================================

Quick Win Implementation: Replaces MPII-B5a to MPII-O1 heuristic rules
with KB feature class queries for Pol II cassette detection.

Usage:
    from grammar_pol2_detector import GrammarPol2Detector

    detector = GrammarPol2Detector(features, sequence)
    cassettes = detector.detect_pol2_cassettes()
"""

import re
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, field


@dataclass
class Pol2Cassette:
    """Detected Pol II expression cassette"""
    module_type: str = 'pol2_expression_mammalian'
    start: int = 0
    end: int = 0
    strand: int = 1
    promoter_info: Dict[str, Any] = field(default_factory=dict)
    polya_info: Dict[str, Any] = field(default_factory=dict)
    components: Dict[str, Any] = field(default_factory=dict)
    confidence: float = 0.0
    score: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return {
            'module_type': self.module_type,
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'promoter': self.promoter_info,
            'polya': self.polya_info,
            'components': self.components,
            'confidence': self.confidence,
            'score': round(self.score, 3),
            'length_bp': self.end - self.start
        }


class GrammarPol2Detector:
    """
    Detect mammalian Pol II expression cassettes using KB feature classes.

    Replaces heuristic rules:
    - MPII-B5a to MPII-B5k: Promoter detection (CMV, EF1a, CAG, etc.)
    - MPII-B3a to MPII-B3h: PolyA signal detection (bGH, SV40, etc.)
    - MPII-E1 to MPII-E4: Component detection (enhancers, introns, Kozak)
    - MPII-I1 to MPII-I4: Internal elements (WPRE, IRES, 2A, signals)
    - MPII-O1: Canonical order validation
    """

    def __init__(self, features: List[Dict], sequence: str, circular: bool = True):
        """
        Initialize detector.

        Args:
            features: List of feature dicts with kb_info
            sequence: Plasmid sequence string
            circular: Whether plasmid is circular
        """
        self.features = features
        self.sequence = sequence
        self.circular = circular
        self.seq_len = len(sequence)

    def detect_pol2_cassettes(self) -> List[Pol2Cassette]:
        """
        Detect all Pol II expression cassettes.

        Returns:
            List of Pol2Cassette objects, sorted by confidence (highest first)
        """
        cassettes = []

        # Find Pol II promoters
        pol2_promoters = self._find_pol2_promoters()

        # Find polyA signals
        polya_signals = self._find_polya_signals()

        # Pair promoters with polyA signals
        for prom in pol2_promoters:
            for polya in polya_signals:
                # Validate pairing
                if not self._validate_promoter_polya_pair(prom, polya):
                    continue

                # Parse components between promoter and polyA
                components = self._parse_components(prom['end'], polya['start'], prom['strand'])

                # Classify promoter and polyA
                prom_info = self._classify_promoter(prom)
                polya_info = self._classify_polya(polya)

                # Calculate confidence and score
                confidence, score = self._calculate_confidence(prom_info, polya_info, components)

                # Build cassette
                cassette = Pol2Cassette(
                    start=prom['start'],
                    end=polya['end'],
                    strand=prom['strand'],
                    promoter_info=prom_info,
                    polya_info=polya_info,
                    components=components,
                    confidence=confidence,
                    score=score
                )

                cassettes.append(cassette)

        # Deduplicate overlapping cassettes
        cassettes = self._deduplicate_cassettes(cassettes)

        # Sort by score (highest first) - confidence is already factored into score
        cassettes.sort(key=lambda c: -c.score)

        return cassettes

    def _find_pol2_promoters(self) -> List[Dict]:
        """
        Find Pol II promoters using KB feature class.

        Replaces MPII-B5a to MPII-B5k heuristic pattern matching.

        Returns:
            List of promoter feature dicts
        """
        pol2_promoters = []

        for f in self.features:
            kb = f.get('kb_info', {})

            # Primary check: polymerase_class
            if kb.get('polymerase_class') == 'pol_ii':
                # Secondary check: host_scope (prefer mammalian)
                host_scope = kb.get('host_scope', [])
                if isinstance(host_scope, str):
                    host_scope = [host_scope]

                # Include if mammalian or no host scope specified
                if not host_scope or 'mammalian' in host_scope or 'insect' in host_scope:
                    pol2_promoters.append(f)

            # Fallback: check feature_class for promoters without polymerase_class
            elif kb.get('feature_class') == 'promoter':
                name = f.get('name', '').lower()
                # Pattern match known Pol II promoters
                if re.search(r'cmv|ef.?1|cag|pgk|ubc|sv40|sffv|rsv|tre', name):
                    pol2_promoters.append(f)

        return pol2_promoters

    def _find_polya_signals(self) -> List[Dict]:
        """
        Find polyA signals using KB feature class.

        Replaces MPII-B3a to MPII-B3h heuristic pattern matching.

        Returns:
            List of polyA signal feature dicts
        """
        polya_signals = []

        for f in self.features:
            kb = f.get('kb_info', {})

            # Primary check: feature_class
            if kb.get('feature_class') == 'polya_signal':
                polya_signals.append(f)

            # Fallback: pattern match in name
            else:
                name = f.get('name', '').lower()
                if re.search(r'poly[\s\-]?a|polya|bgh.*poly|sv40.*poly|hgh.*poly', name):
                    polya_signals.append(f)

        return polya_signals

    def _validate_promoter_polya_pair(self, promoter: Dict, polya: Dict) -> bool:
        """
        Validate that promoter and polyA can form a cassette.

        Rules:
        - Must be on same strand
        - PolyA must be downstream of promoter
        - Distance must be reasonable (< 15kb)

        Args:
            promoter: Promoter feature dict
            polya: PolyA signal feature dict

        Returns:
            True if valid pairing, False otherwise
        """
        # Same strand check
        if promoter['strand'] != polya['strand']:
            return False

        # Order check (polyA downstream of promoter)
        if promoter['strand'] == 1:
            if polya['start'] <= promoter['end']:
                return False
        else:
            if polya['end'] >= promoter['start']:
                return False

        # Distance check (max cassette size ~15kb)
        distance = abs(polya['start'] - promoter['end'])
        if distance > 15000:
            return False

        # Minimum distance check (at least 100bp for minimal CDS)
        if distance < 100:
            return False

        return True

    def _classify_promoter(self, promoter: Dict) -> Dict[str, Any]:
        """
        Classify Pol II promoter type and strength.

        Replaces MPII-B5a to MPII-B5k pattern matching with structured classification.

        Args:
            promoter: Promoter feature dict

        Returns:
            Dict with type, strength, weight, and metadata
        """
        name = promoter.get('name', '').lower()
        kb = promoter.get('kb_info', {})

        # Strong promoters
        if re.search(r'cmv|cytomegalovirus', name):
            enhanced = 'enhancer' in name or 'enhanced' in name
            return {
                'type': 'CMV',
                'full_name': promoter.get('name', 'CMV promoter'),
                'strength': 'strong',
                'weight': 0.95 if enhanced else 0.90,
                'enhanced': enhanced,
                'feature': promoter
            }

        elif re.search(r'ef.?1|ef1a', name):
            has_intron = 'intron' in name
            return {
                'type': 'EF1a',
                'full_name': promoter.get('name', 'EF1a promoter'),
                'strength': 'strong',
                'weight': 0.97 if has_intron else 0.90,
                'has_intron': has_intron,
                'feature': promoter
            }

        elif re.search(r'cag', name):
            return {
                'type': 'CAG',
                'full_name': promoter.get('name', 'CAG promoter'),
                'strength': 'strong',
                'weight': 0.97,
                'note': 'CMV enhancer + chicken β-actin promoter + intron',
                'feature': promoter
            }

        # Moderate promoters
        elif re.search(r'pgk|phosphoglycerate', name):
            return {
                'type': 'PGK',
                'full_name': promoter.get('name', 'PGK promoter'),
                'strength': 'moderate',
                'weight': 0.88,
                'feature': promoter
            }

        elif re.search(r'ubc|ubiquitin.*c', name):
            return {
                'type': 'UbC',
                'full_name': promoter.get('name', 'UbC promoter'),
                'strength': 'moderate',
                'weight': 0.88,
                'feature': promoter
            }

        elif re.search(r'sv40', name):
            return {
                'type': 'SV40',
                'full_name': promoter.get('name', 'SV40 promoter'),
                'strength': 'moderate',
                'weight': 0.85,
                'feature': promoter
            }

        # Tissue-specific
        elif re.search(r'hsyn|camkii', name):
            return {
                'type': 'neuronal',
                'full_name': promoter.get('name', 'Neuronal promoter'),
                'strength': 'strong',
                'weight': 0.90,
                'tissue': 'neuronal',
                'feature': promoter
            }

        # Inducible
        elif re.search(r'tre|tet.*response', name):
            return {
                'type': 'TRE',
                'full_name': promoter.get('name', 'TRE promoter'),
                'strength': 'strong',
                'weight': 0.90,
                'inducible': 'doxycycline',
                'feature': promoter
            }

        # Generic Pol II
        return {
            'type': 'unknown_pol2',
            'full_name': promoter.get('name', 'Unknown Pol II promoter'),
            'strength': 'weak',
            'weight': 0.70,
            'feature': promoter
        }

    def _classify_polya(self, polya: Dict) -> Dict[str, Any]:
        """
        Classify polyA signal type.

        Replaces MPII-B3a to MPII-B3h pattern matching.

        Args:
            polya: PolyA signal feature dict

        Returns:
            Dict with type, weight, and metadata
        """
        name = polya.get('name', '').lower()

        if re.search(r'bgh|bovine.*growth', name):
            return {
                'type': 'bGH',
                'full_name': polya.get('name', 'bGH poly(A)'),
                'weight': 0.90,
                'note': 'Default for AAV vectors',
                'feature': polya
            }

        elif re.search(r'sv40', name):
            return {
                'type': 'SV40',
                'full_name': polya.get('name', 'SV40 poly(A)'),
                'weight': 0.90,
                'note': 'Common in selection cassettes',
                'feature': polya
            }

        elif re.search(r'hgh|human.*growth', name):
            return {
                'type': 'hGH',
                'full_name': polya.get('name', 'hGH poly(A)'),
                'weight': 0.88,
                'feature': polya
            }

        elif re.search(r'globin', name):
            return {
                'type': 'beta_globin',
                'full_name': polya.get('name', 'β-globin poly(A)'),
                'weight': 0.88,
                'feature': polya
            }

        return {
            'type': 'generic',
            'full_name': polya.get('name', 'poly(A) signal'),
            'weight': 0.70,
            'feature': polya
        }

    def _parse_components(self, start: int, end: int, strand: int) -> Dict[str, Any]:
        """
        Parse internal components between promoter and polyA.

        Replaces MPII-E1 to MPII-I4 component detection rules.

        Args:
            start: Start position (promoter end)
            end: End position (polyA start)
            strand: Strand (1 or -1)

        Returns:
            Dict with component lists
        """
        components = {
            'enhancers': [],
            'introns': [],
            'kozak': None,
            'cds': [],
            'wpre': None,
            'ires': [],
            '2a_peptides': [],
            'signal_peptides': [],
            'tags': [],
            'nls': []
        }

        for f in self.features:
            # Skip if not in region
            if not (start <= f['start'] < end):
                continue

            # Skip if wrong strand
            if f['strand'] != strand:
                continue

            kb = f.get('kb_info', {})
            fc = kb.get('feature_class')
            sc = kb.get('subclass')
            name = f.get('name', '').lower()

            # MPII-E1: Enhancer
            if fc == 'enhancer':
                components['enhancers'].append(f)

            # MPII-E2: Intron
            elif fc == 'intron':
                components['introns'].append(f)

            # MPII-E3: Kozak
            elif fc == 'kozak':
                components['kozak'] = f

            # CDS
            elif fc == 'cds':
                components['cds'].append(f)

            # MPII-I1: WPRE
            elif 'wpre' in name or 'woodchuck' in name:
                components['wpre'] = f

            # MPII-I2: IRES
            elif fc == 'ires':
                components['ires'].append(f)

            # MPII-I3: 2A peptide
            elif fc == 'peptide' and sc == 'self_cleaving':
                components['2a_peptides'].append(f)

            # MPII-I4: Signal peptide
            elif fc == 'signal_peptide':
                components['signal_peptides'].append(f)

            # Tags
            elif fc == 'tag':
                components['tags'].append(f)

            # NLS
            elif fc == 'localization_signal':
                components['nls'].append(f)

        return components

    def _calculate_confidence(
        self,
        prom_info: Dict,
        polya_info: Dict,
        components: Dict
    ) -> Tuple[float, float]:
        """
        Calculate cassette confidence and score.

        Replaces MPII-O1 canonical order validation with structured scoring.

        Args:
            prom_info: Classified promoter info
            polya_info: Classified polyA info
            components: Parsed components dict

        Returns:
            Tuple of (confidence, score)
        """
        # Base score from promoter and polyA weights
        score = (prom_info['weight'] + polya_info['weight']) / 2

        # Bonus for complete cassette components
        if components['cds']:
            score += 0.10

        if components['kozak']:
            score += 0.05

        if components['wpre']:
            score += 0.05

        if components['introns']:
            score += 0.03

        if components['enhancers']:
            score += 0.02

        # Cap at 1.0
        score = min(score, 1.0)

        # Confidence categories
        if score >= 0.85:
            confidence = 'high'
        elif score >= 0.60:
            confidence = 'medium'
        else:
            confidence = 'low'

        return confidence, score

    def _deduplicate_cassettes(self, cassettes: List[Pol2Cassette]) -> List[Pol2Cassette]:
        """
        Remove overlapping cassettes, keeping highest scoring.

        Args:
            cassettes: List of Pol2Cassette objects

        Returns:
            Deduplicated list
        """
        if not cassettes:
            return []

        # Sort by score (highest first)
        sorted_cassettes = sorted(cassettes, key=lambda c: -c.score)

        kept = []
        for cassette in sorted_cassettes:
            # Check if significantly overlaps any kept cassette
            overlaps = False
            for kept_cassette in kept:
                if cassette.strand != kept_cassette.strand:
                    continue

                # Calculate overlap
                overlap_start = max(cassette.start, kept_cassette.start)
                overlap_end = min(cassette.end, kept_cassette.end)

                if overlap_end > overlap_start:
                    overlap_len = overlap_end - overlap_start
                    cassette_len = cassette.end - cassette.start

                    # If >70% overlap, consider it duplicate
                    if cassette_len > 0 and overlap_len / cassette_len > 0.7:
                        overlaps = True
                        break

            if not overlaps:
                kept.append(cassette)

        return kept


# Convenience function for simple usage
def detect_pol2_cassettes(features: List[Dict], sequence: str, circular: bool = True) -> List[Dict]:
    """
    Detect Pol II expression cassettes (simple interface).

    Args:
        features: List of feature dicts with kb_info
        sequence: Plasmid sequence string
        circular: Whether plasmid is circular

    Returns:
        List of cassette dicts
    """
    detector = GrammarPol2Detector(features, sequence, circular)
    cassettes = detector.detect_pol2_cassettes()
    return [c.to_dict() for c in cassettes]
