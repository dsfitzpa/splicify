"""
Module Linter - Validation and Feature Assignment

Provides:
- Module component role system
- Grammar annotation for gaps
- Module-level rule engine for validation
- Lint issue detection

This module is part of the annotation pipeline restructuring to create
cleaner separation between motif/boundary detection, validation/linting,
and inter-module analysis.
"""

from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple, Any, Set
from enum import Enum
import re
import uuid


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class LintIssue:
    """Represents a validation/lint issue"""
    issue_id: str
    rule_id: str
    severity: str  # "error", "warning", "info"
    category: str
    message: str
    feature_ids: List[str] = field(default_factory=list)
    module_ids: List[str] = field(default_factory=list)
    location: Optional[Dict[str, int]] = None
    suggestion: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ComponentRole:
    """Defines a component role within a module"""
    role: str
    required: bool
    multiple: bool  # Can there be multiple of this component?
    links_to_submodule: bool = False
    submodule_type: Optional[str] = None


@dataclass
class GapAnnotation:
    """Represents an annotated gap between features"""
    gap_type: str  # 'linker', 'tag', 'utr', 'unknown'
    start: int
    end: int
    strand: int
    upstream_feature_id: Optional[str] = None
    downstream_feature_id: Optional[str] = None
    sequence: str = ""
    description: str = ""


# =============================================================================
# MODULE COMPONENT SCHEMA
# =============================================================================

MODULE_COMPONENT_SCHEMA = {
    'pol2_expression': {
        'promoter': ComponentRole('promoter', required=True, multiple=True),
        'enhancer': ComponentRole('enhancer', required=False, multiple=True),
        'orf': ComponentRole('orf', required=True, multiple=False, links_to_submodule=True, submodule_type='transcript'),
        'polyA_signal': ComponentRole('polyA_signal', required=True, multiple=False),
    },
    'transcript': {
        'start': ComponentRole('start_codon', required=True, multiple=False),
        'cds': ComponentRole('cds', required=True, multiple=True),
        'intron': ComponentRole('intron', required=False, multiple=True, links_to_submodule=True, submodule_type='intron'),
        '2a_linker': ComponentRole('2a_linker', required=False, multiple=True, links_to_submodule=True, submodule_type='protein'),
        'ires': ComponentRole('ires', required=False, multiple=True, links_to_submodule=True, submodule_type='protein'),
        'stop': ComponentRole('stop_codon', required=True, multiple=False),
    },
    'protein': {
        'start': ComponentRole('start_codon', required=True, multiple=False),  # Required for first protein only
        'cds': ComponentRole('cds', required=True, multiple=True),
        'stop': ComponentRole('stop_codon', required=True, multiple=False),  # Required for last protein only
    },
    'pol3_expression': {
        'promoter': ComponentRole('pol3_promoter', required=True, multiple=False),
        'guide_sequence': ComponentRole('guide_sequence', required=False, multiple=False),
        'scaffold': ComponentRole('scaffold_rna', required=False, multiple=False),
        'terminator': ComponentRole('pol3_terminator', required=False, multiple=False),
    },
    'bacterial_backbone': {
        'origin': ComponentRole('origin', required=True, multiple=False),
        'marker': ComponentRole('bacterial_marker', required=False, multiple=True),
        'regulatory': ComponentRole('bacterial_regulatory', required=False, multiple=True),
    },
    'lentiviral_payload': {
        '5_ltr': ComponentRole('ltr', required=True, multiple=False),
        '3_ltr': ComponentRole('ltr', required=True, multiple=False),
        'psi': ComponentRole('psi', required=False, multiple=False),
        'rre': ComponentRole('rre', required=False, multiple=False),
        'cppt': ComponentRole('cppt', required=False, multiple=False),
        'wpre': ComponentRole('wpre', required=False, multiple=False),
    },
    'aav_payload': {
        '5_itr': ComponentRole('itr', required=True, multiple=False),
        '3_itr': ComponentRole('itr', required=True, multiple=False),
    },
}


# =============================================================================
# LENTIVIRAL INDICATORS (moved from plasmid_analyzer.py)
# =============================================================================

LENTIVIRAL_INDICATORS = [
    ('psi', ['hiv', 'packaging']),  # HIV Psi
    ('rre', []),  # RRE
    ('cppt', []), ('cts', []),  # cPPT/CTS
    ('wpre', []), ('woodchuck', []),  # WPRE
    ('tar', ['mir', 'hiv']),  # TAR element
    ('dis', ['hiv']),  # Dimerization initiation site
    ('rnai', []),  # RNAI
]


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def generate_uuid() -> str:
    """Generate a unique identifier"""
    return str(uuid.uuid4())[:8]


def _is_bacterial_promoter_by_name(name_lower: str) -> bool:
    """Check if feature name indicates a bacterial promoter."""
    if 'promoter' not in name_lower:
        return False
    bacterial_markers = ['amp', 'bla', 'kan', 'neo', 'cm', 'cat', 'tet', 'lac', 'spec', 'gen', 'strep']
    return any(marker in name_lower for marker in bacterial_markers)


# =============================================================================
# CONSTRUCT ROLE FUNCTION (moved from plasmid_analyzer.py)
# =============================================================================

def construct_role(
    feature_name: str,
    feature_type: str = "",
    description: str = "",
    kb_info: Optional[Dict[str, Any]] = None
) -> str:
    """
    Determine the functional role of a feature based on its name, type, and KB lookup.

    Returns roles like: promoter, pol2_promoter, pol3_promoter, terminator,
    polya, ori, marker, scaffold_rna, backbone_accessory, etc.
    """
    name_lower = (feature_name or "").lower()
    type_lower = (feature_type or "").lower()

    # === 1. First check for specific well-known features by name ===

    # gRNA scaffold detection (before KB lookup)
    if "scaffold" in name_lower:
        if any(x in name_lower for x in ["grna", "sgrna", "guide", "crispr"]):
            return "scaffold_rna"

    # Lentiviral cis-acting RNA elements
    for main_term, context_terms in LENTIVIRAL_INDICATORS:
        if main_term in name_lower:
            if not context_terms or any(ctx in name_lower for ctx in context_terms):
                return 'lentiviral_cis_rna'
    
    # Also check for HIV-1 stem loops and splice sites
    if name_lower.startswith('hiv'):
        if any(x in name_lower for x in ['sl1', 'sl2', 'sl3', 'sl4', 'stem', 'sd', 'sa', 'splice']):
            return 'lentiviral_cis_rna'

    # === 2. Use KnowledgeBase for authoritative lookup ===
    if kb_info:
        feature_class = kb_info.get('feature_class', '')
        polymerase_class = kb_info.get('polymerase_class', '')
        host_scope = set(kb_info.get('host_scope', []))

        # Handle promoters with polymerase class
        if feature_class == 'promoter':
            if polymerase_class == 'pol_iii':
                return 'pol3_promoter'
            elif polymerase_class == 'pol_ii':
                return 'pol2_promoter'
            elif polymerase_class == 'bacterial_or_phage':
                return 'promoter'
            else:
                if host_scope & {'mammalian', 'plant', 'yeast', 'insect'}:
                    return 'pol2_promoter'
                return 'promoter'

        # Handle terminators
        if feature_class == 'terminator':
            return 'terminator'

        # Handle polyA signals
        if feature_class == 'polyA_signal':
            return 'polya'

        # Handle origins
        if feature_class == 'origin':
            return 'origin'

        # Handle selectable markers
        if feature_class in {'selectable_marker', 'selection_marker', 'antibiotic_resistance'}:
            if host_scope & {'bacterial', 'prokaryotic'}:
                return 'bacterial_marker'
            elif host_scope & {'mammalian', 'plant', 'yeast', 'insect'}:
                return 'selection_payload'
            # Common mammalian selection markers even if host_scope unknown
            mammalian_markers = ['puror', 'puro', 'neor', 'neo', 'hygror', 'hygro', 'bleor', 'bleo', 'zeor', 'zeo', 'bsr', 'bsd']
            if any(m in (feature_name or '').lower() for m in mammalian_markers):
                return 'selection_payload'
            return 'bacterial_marker'

        # Handle gRNA scaffold from KB
        if feature_class == 'grna_scaffold':
            return 'scaffold_rna'

        # Handle enhancers
        if feature_class == 'enhancer':
            return 'enhancer'

        # Handle reporters
        # Handle nuclease (Cas9, Cas12, etc.)
        if feature_class == "nuclease":
            return "editing_payload"

        if feature_class == "reporter":
            return 'reporter_payload'

        # Handle CDS/protein coding
        if feature_class in {'cds', 'cds_payload', 'protein_coding'}:
            if any(x in (feature_name or '').lower() for x in ['cas9', 'cas12', 'crispr', 'nuclease', 'talen', 'zfn']):
                return 'editing_payload'
            if any(x in (feature_name or '').lower() for x in ['gfp', 'rfp', 'yfp', 'cfp', 'mcherry', 'luciferase', 'lacz']):
                return 'reporter_payload'
            # Check for bacterial/eukaryotic markers before returning expression_payload
            bacterial_cds_markers = ['ampr', 'kanr', 'camr', 'specr', 'tetr', 'cmr', 'genr', 'strepr', 'bla']
            eukaryotic_cds_markers = ['puror', 'neor', 'hygr', 'zeor', 'blastr', 'bleor']
            name_lower = (feature_name or '').lower()
            if any(m in name_lower for m in bacterial_cds_markers):
                return 'bacterial_marker'
            if any(m in name_lower for m in eukaryotic_cds_markers):
                return 'selection_payload'
            return 'expression_payload'
            
        # Handle LTR/ITR
        if feature_class in {'ltr', 'long_terminal_repeat'}:
            return 'ltr'
        if feature_class in {'itr', 'inverted_terminal_repeat'}:
            return 'itr'

    # === 3. Fallback: Name-based detection ===

    # Pol III promoters by name
    pol3_indicators = ['u6', 'h1', '7sk', 'u3', 'trna', 'rrna', '5s']
    if 'promoter' in name_lower or type_lower == 'promoter':
        for ind in pol3_indicators:
            if ind in name_lower:
                return 'pol3_promoter'

    # Pol II promoters by name
    pol2_indicators = ['cmv', 'ef1a', 'ef-1a', 'cag', 'pgk', 'sv40', 'ubc', 'hsyn',
                       'synapsin', 'ac5', 'actin', 'tre', 'rosa']
    if 'promoter' in name_lower or type_lower == 'promoter':
        for ind in pol2_indicators:
            if ind in name_lower:
                return 'pol2_promoter'

    # Generic promoter
    if 'promoter' in name_lower or type_lower == 'promoter':
        return 'promoter'

    # Terminators
    if 'terminator' in name_lower or 'term' in name_lower:
        return 'terminator'

    # PolyA signals
    # PolyA signals - comprehensive pattern matching
    # Match patterns like: bGH_poly(A)_signal, SV40_poly(A)_signal_(2), polyadenylation, etc.
    polya_patterns = [
        'polya', 'poly(a)', 'poly-a', 'poly_a',
        'polyadenylation', 'polyadenyl',
        'pa_signal', 'pa signal',
    ]
    if any(pat in name_lower for pat in polya_patterns):
        return 'polya'
    # Also check for "signal" combined with common polyA sources
    if 'signal' in name_lower and any(src in name_lower for src in ['bgh', 'sv40', 'hgh', 'rbglob', 'bovine']):
        return 'polya'


    # Origins
    if 'ori' in name_lower or 'origin' in name_lower:
        return 'origin'

    # Selectable markers
    bacterial_markers = ['ampr', 'kanr', 'camr', 'specr', 'tetr', 'cmr', 'genr', 'strepr']
    eukaryotic_markers = ['puror', 'neor', 'hygr', 'zeor', 'blastr', 'pac', 'neo']
    
    for ind in bacterial_markers:
        if ind in name_lower:
            return 'bacterial_marker'
    for ind in eukaryotic_markers:
        if ind in name_lower:
            return 'selection_payload'

    # gRNA scaffold
    if 'scaffold' in name_lower or 'sgrna' in name_lower:
        return 'scaffold_rna'

    # LTR/ITR
    if 'ltr' in name_lower:
        return 'ltr'
    if 'itr' in name_lower:
        return 'itr'

    # Backbone accessory features
    backbone_accessory = ['bom', 'rop', 'rom', 'f1 ori', 'flori']
    for acc in backbone_accessory:
        if acc in name_lower:
            return 'backbone_accessory'

    # Default: use feature type if available
    if type_lower:
        return type_lower

    return 'misc_feature'


# =============================================================================
# GAP ANNOTATOR CLASS
# =============================================================================

class GapAnnotator:
    """Annotate gaps between features based on context"""

    def __init__(self, sequence: str, circular: bool = True):
        self.sequence = sequence
        self.seq_len = len(sequence)
        self.circular = circular

    def annotate_linker(
        self,
        upstream_cds: Dict[str, Any],
        downstream_cds: Dict[str, Any],
        gap_seq: str
    ) -> Optional[GapAnnotation]:
        """
        In-frame gap between two CDS = linker.
        
        Check if gap is in-frame and could be a linker peptide.
        """
        # Check same strand
        if upstream_cds.get('strand', 1) != downstream_cds.get('strand', 1):
            return None
        
        strand = upstream_cds.get('strand', 1)
        gap_start = upstream_cds.get('end', 0)
        gap_end = downstream_cds.get('start', 0)
        gap_length = gap_end - gap_start
        
        # Linkers are typically short (3-90bp = 1-30 amino acids)
        if gap_length < 3 or gap_length > 90:
            return None
        
        # Check if in-frame (gap length divisible by 3)
        if gap_length % 3 != 0:
            return None
        
        return GapAnnotation(
            gap_type='linker',
            start=gap_start,
            end=gap_end,
            strand=strand,
            upstream_feature_id=upstream_cds.get('instance_id'),
            downstream_feature_id=downstream_cds.get('instance_id'),
            sequence=gap_seq,
            description=f'Linker peptide ({gap_length // 3} aa)',
        )

    def annotate_tag(
        self,
        start_or_stop: Dict[str, Any],
        cds: Dict[str, Any],
        gap_seq: str
    ) -> Optional[GapAnnotation]:
        """
        Gap between start/first CDS or last CDS/stop = tag.
        
        Epitope tags are typically at N- or C-terminus.
        """
        strand = cds.get('strand', 1)
        is_at_start = start_or_stop.get('end', 0) <= cds.get('start', 0)
        
        if is_at_start:
            gap_start = start_or_stop.get('end', 0)
            gap_end = cds.get('start', 0)
            position = 'N-terminal'
        else:
            gap_start = cds.get('end', 0)
            gap_end = start_or_stop.get('start', 0)
            position = 'C-terminal'
        
        gap_length = gap_end - gap_start
        
        # Tags are typically short (15-120bp = 5-40 amino acids)
        if gap_length < 15 or gap_length > 120:
            return None
        
        # Check if in-frame
        if gap_length % 3 != 0:
            return None
        
        return GapAnnotation(
            gap_type='tag',
            start=gap_start,
            end=gap_end,
            strand=strand,
            upstream_feature_id=start_or_stop.get('instance_id') if is_at_start else cds.get('instance_id'),
            downstream_feature_id=cds.get('instance_id') if is_at_start else start_or_stop.get('instance_id'),
            sequence=gap_seq,
            description=f'{position} tag ({gap_length // 3} aa)',
        )

    def annotate_utr(
        self,
        promoter: Dict[str, Any],
        cds: Dict[str, Any],
        gap_seq: str
    ) -> Optional[GapAnnotation]:
        """
        Gap between promoter and CDS = 5'UTR.
        """
        strand = cds.get('strand', 1)
        
        if strand == 1:
            gap_start = promoter.get('end', 0)
            gap_end = cds.get('start', 0)
        else:
            gap_start = cds.get('end', 0)
            gap_end = promoter.get('start', 0)
        
        gap_length = gap_end - gap_start
        
        # UTRs are typically 50-500bp
        if gap_length < 10 or gap_length > 1000:
            return None
        
        return GapAnnotation(
            gap_type='utr',
            start=gap_start,
            end=gap_end,
            strand=strand,
            upstream_feature_id=promoter.get('instance_id'),
            downstream_feature_id=cds.get('instance_id'),
            sequence=gap_seq,
            description=f"5' UTR ({gap_length} bp)",
        )

    def annotate_gaps_in_module(
        self,
        features: List[Dict[str, Any]]
    ) -> List[GapAnnotation]:
        """
        Annotate all gaps within a set of features (typically a module).
        """
        annotations = []
        
        # Sort features by position
        sorted_features = sorted(features, key=lambda f: f.get('start', 0))
        
        for i in range(len(sorted_features) - 1):
            upstream = sorted_features[i]
            downstream = sorted_features[i + 1]
            
            # Calculate gap
            gap_start = upstream.get('end', 0)
            gap_end = downstream.get('start', 0)
            
            if gap_end <= gap_start:
                continue  # No gap or overlap
            
            gap_seq = self.sequence[gap_start:gap_end] if gap_start < gap_end <= self.seq_len else ""
            
            # Try to annotate based on context
            upstream_role = upstream.get('role', '')
            downstream_role = downstream.get('role', '')
            
            # Check for linker between CDS
            if 'payload' in upstream_role and 'payload' in downstream_role:
                ann = self.annotate_linker(upstream, downstream, gap_seq)
                if ann:
                    annotations.append(ann)
                    continue
            
            # Check for tag
            if 'payload' in upstream_role or 'payload' in downstream_role:
                cds = upstream if 'payload' in upstream_role else downstream
                other = downstream if 'payload' in upstream_role else upstream
                ann = self.annotate_tag(other, cds, gap_seq)
                if ann:
                    annotations.append(ann)
                    continue
            
            # Check for UTR
            if 'promoter' in upstream_role and 'payload' in downstream_role:
                ann = self.annotate_utr(upstream, downstream, gap_seq)
                if ann:
                    annotations.append(ann)
        
        return annotations


# =============================================================================
# MODULE LINTER CLASS
# =============================================================================

class ModuleLinter:
    """Apply module-level validation rules"""

    def __init__(self):
        self.issues: List[LintIssue] = []

    def check_bacterial_missing_ori(
        self,
        modules: List[Dict[str, Any]],
        features: List[Dict[str, Any]]
    ) -> List[LintIssue]:
        """
        Check for AmpR promoter + resistance marker without bacterial Ori.
        """
        issues = []
        
        # Find bacterial markers
        bacterial_markers = [
            f for f in features
            if f.get('role') == 'bacterial_marker'
        ]
        
        # Find bacterial origins
        bacterial_origins = [
            f for f in features
            if f.get('role') == 'origin'
        ]
        
        # Find bacterial promoters
        bacterial_promoters = [
            f for f in features
            if _is_bacterial_promoter_by_name(f.get('feature_name', '').lower())
        ]
        
        if bacterial_markers and bacterial_promoters and not bacterial_origins:
            issues.append(LintIssue(
                issue_id=f"lint_{generate_uuid()}",
                rule_id="bacterial_missing_ori",
                severity="error",
                category="backbone",
                message="Bacterial marker cassette detected but no bacterial origin of replication found",
                feature_ids=[f.get('instance_id') for f in bacterial_markers],
                suggestion="Add a bacterial origin (ColE1, pBR322, p15A) for plasmid propagation in E. coli"
            ))
        
        return issues

    def check_lentiviral_missing_elements(
        self,
        modules: List[Dict[str, Any]],
        features: List[Dict[str, Any]]
    ) -> List[LintIssue]:
        """
        Check for lentiviral payload missing psi, rre, cppt, or wpre.
        """
        issues = []
        
        # Check if this is a lentiviral vector
        lenti_modules = [
            m for m in modules
            if 'lentiviral' in m.get('module_type', '').lower()
        ]
        
        if not lenti_modules:
            return issues
        
        # Find lentiviral elements
        feature_names = [f.get('feature_name', '').lower() for f in features]
        feature_roles = [f.get('role', '').lower() for f in features]
        
        required_elements = {
            'psi': 'Packaging signal (Psi)',
            'rre': 'Rev Response Element (RRE)',
        }
        
        recommended_elements = {
            'cppt': 'Central Polypurine Tract (cPPT)',
            'wpre': 'Woodchuck Hepatitis Virus Posttranscriptional Regulatory Element (WPRE)',
        }
        
        # Check required elements
        for elem, desc in required_elements.items():
            found = any(
                elem in name or 'lentiviral_cis_rna' in role
                for name, role in zip(feature_names, feature_roles)
            )
            if not found:
                issues.append(LintIssue(
                    issue_id=f"lint_{generate_uuid()}",
                    rule_id=f"lentiviral_missing_{elem}",
                    severity="warning",
                    category="lentiviral_packaging",
                    message=f"Lentiviral vector missing {desc}",
                    module_ids=[m.get('module_id') for m in lenti_modules],
                    suggestion=f"Consider adding {elem.upper()} for optimal lentiviral packaging"
                ))
        
        # Check recommended elements
        for elem, desc in recommended_elements.items():
            found = any(elem in name for name in feature_names)
            if not found:
                issues.append(LintIssue(
                    issue_id=f"lint_{generate_uuid()}",
                    rule_id=f"lentiviral_recommend_{elem}",
                    severity="info",
                    category="lentiviral_packaging",
                    message=f"Lentiviral vector could benefit from {desc}",
                    module_ids=[m.get('module_id') for m in lenti_modules],
                    suggestion=f"Consider adding {elem.upper()} for improved expression"
                ))
        
        return issues

    def check_module_component_requirements(
        self,
        module: Dict[str, Any],
        features: List[Dict[str, Any]],
        schema: Dict[str, ComponentRole]
    ) -> List[LintIssue]:
        """
        Validate module has all required components from schema.
        """
        issues = []
        module_type = module.get('module_type', '')
        module_id = module.get('module_id', '')
        
        # Get features in this module
        module_feature_ids = set(module.get('feature_ids', []))
        module_features = [f for f in features if f.get('instance_id') in module_feature_ids]
        
        # Check each required component
        for component_name, role_spec in schema.items():
            if not role_spec.required:
                continue
            
            # Count matching features
            matching = [
                f for f in module_features
                if f.get('role') == role_spec.role or role_spec.role in f.get('role', '')
            ]
            
            if not matching:
                issues.append(LintIssue(
                    issue_id=f"lint_{generate_uuid()}",
                    rule_id=f"missing_{component_name}",
                    severity="warning",
                    category="module_structure",
                    message=f"{module_type} module missing required {component_name}",
                    module_ids=[module_id],
                    suggestion=f"Add a {role_spec.role} element to complete the module"
                ))
        
        return issues

    def run_all_checks(
        self,
        modules: List[Dict[str, Any]],
        features: List[Dict[str, Any]]
    ) -> List[LintIssue]:
        """Run all module-level linting checks"""
        all_issues = []
        
        all_issues.extend(self.check_bacterial_missing_ori(modules, features))
        all_issues.extend(self.check_lentiviral_missing_elements(modules, features))
        
        # Check component requirements for each module
        for module in modules:
            module_type = module.get('module_type', '')
            
            # Map module type to schema
            if 'pol2' in module_type:
                schema = MODULE_COMPONENT_SCHEMA.get('pol2_expression', {})
            elif 'pol3' in module_type:
                schema = MODULE_COMPONENT_SCHEMA.get('pol3_expression', {})
            elif 'bacterial' in module_type:
                schema = MODULE_COMPONENT_SCHEMA.get('bacterial_backbone', {})
            elif 'lentiviral' in module_type:
                schema = MODULE_COMPONENT_SCHEMA.get('lentiviral_payload', {})
            elif 'aav' in module_type:
                schema = MODULE_COMPONENT_SCHEMA.get('aav_payload', {})
            else:
                continue
            
            all_issues.extend(
                self.check_module_component_requirements(module, features, schema)
            )
        
        return all_issues


# =============================================================================
# RULE ENGINE CLASS (moved from plasmid_analyzer.py)
# =============================================================================

class RuleEngine:
    """
    Rule-based validation engine for plasmid constructs.
    
    This is a compatibility layer that wraps ModuleLinter and provides
    the same interface as the original RuleEngine.
    """

    @staticmethod
    def check_orientation_violations(
        features: List[Dict[str, Any]],
        modules: List[Dict[str, Any]]
    ) -> List[LintIssue]:
        """Check for promoter/payload strand disagreement"""
        issues = []

        for module in modules:
            if 'cassette' not in module.get('module_type', ''):
                continue

            module_feature_ids = set(module.get('feature_ids', []))
            module_features = [f for f in features if f.get('instance_id') in module_feature_ids]
            promoters = [f for f in module_features if 'promoter' in f.get('role', '')]
            payloads = [f for f in module_features if 'payload' in f.get('role', '')]

            for promoter in promoters:
                for payload in payloads:
                    if promoter.get('strand', 1) != payload.get('strand', 1):
                        issues.append(LintIssue(
                            issue_id=f"lint_{generate_uuid()}",
                            rule_id="orientation_mismatch",
                            severity="error",
                            category="orientation",
                            message=f"Promoter '{promoter.get('feature_name')}' and payload '{payload.get('feature_name')}' are on opposite strands",
                            feature_ids=[promoter.get('instance_id'), payload.get('instance_id')],
                            module_ids=[module.get('module_id')],
                            suggestion="Ensure promoter and payload are on the same strand for proper transcription"
                        ))

        return issues

    @staticmethod
    def check_pol3_misuse(
        features: List[Dict[str, Any]],
        modules: List[Dict[str, Any]]
    ) -> List[LintIssue]:
        """Check for Pol III promoter driving CDS (inappropriate use)"""
        issues = []

        sorted_features = sorted(features, key=lambda f: f.get('start', 0))
        pol3_promoters = [f for f in sorted_features if f.get('role') == "pol3_promoter"]
        cds_features = [f for f in sorted_features if f.get('role') in ("expression_payload", "editing_payload", "reporter_payload")]
        guide_features = [f for f in sorted_features if f.get('role') in ("guide_rna", "scaffold_rna")]

        for prom in pol3_promoters:
            # Check if there's a guide RNA immediately after this promoter (within 500bp)
            downstream_guides = [g for g in guide_features
                               if g.get('strand', 1) == prom.get('strand', 1)
                               and g.get('start', 0) > prom.get('end', 0)
                               and g.get('start', 0) - prom.get('end', 0) < 500]

            if downstream_guides:
                continue

            # Find downstream CDS on same strand within 1000bp
            downstream_cds = [c for c in cds_features
                           if c.get('strand', 1) == prom.get('strand', 1)
                           and c.get('start', 0) > prom.get('end', 0)
                           and c.get('start', 0) - prom.get('end', 0) < 1000]
            if downstream_cds:
                cds = downstream_cds[0]
                issues.append(LintIssue(
                    issue_id=f"lint_{generate_uuid()}",
                    rule_id="pol3_driving_cds",
                    severity="error",
                    category="polymerase_misuse",
                    message=f"Pol III promoter '{prom.get('feature_name')}' appears to be driving CDS '{cds.get('feature_name')}' - Pol III cannot efficiently translate mRNA to protein",
                    feature_ids=[prom.get('instance_id'), cds.get('instance_id')],
                    suggestion="Use a Pol II promoter (CMV, EF1a, CAG) for protein expression. Pol III is only suitable for small non-coding RNAs like sgRNA or shRNA."
                ))

        return issues

    @staticmethod
    def check_frame_continuity(
        features: List[Dict[str, Any]],
        modules: List[Dict[str, Any]]
    ) -> List[LintIssue]:
        """Check for frame continuity issues: missing start, weak Kozak, premature stops"""
        issues = []

        for f in features:
            if 'payload' not in f.get('role', ''):
                continue

            # Missing start codon
            if not f.get('has_start_codon', True) and f.get('length_bp', 0) > 100:
                issues.append(LintIssue(
                    issue_id=f"lint_{generate_uuid()}",
                    rule_id="missing_start_codon",
                    severity="error",
                    category="frame_continuity",
                    message=f"CDS '{f.get('feature_name')}' lacks a start codon (ATG) at the expected position",
                    feature_ids=[f.get('instance_id')],
                    suggestion="Ensure the CDS begins with ATG"
                ))

            # Weak Kozak context
            if f.get('has_start_codon') and f.get('kozak_strength') == "weak":
                issues.append(LintIssue(
                    issue_id=f"lint_{generate_uuid()}",
                    rule_id="weak_kozak_context",
                    severity="warning",
                    category="frame_continuity",
                    message=f"CDS '{f.get('feature_name')}' has weak Kozak consensus, which may reduce translation efficiency",
                    feature_ids=[f.get('instance_id')],
                    suggestion="Consider optimizing the Kozak sequence (ideally GCCACCATGG)"
                ))

            # Premature stop codons (from feature metadata)
            if f.get('internal_stops', 0) > 0:
                issues.append(LintIssue(
                    issue_id=f"lint_{generate_uuid()}",
                    rule_id="premature_stop_codon",
                    severity="error",
                    category="frame_continuity",
                    message=f"CDS '{f.get('feature_name')}' has {f.get('internal_stops')} internal stop codon(s), resulting in truncated protein",
                    feature_ids=[f.get('instance_id')],
                    suggestion="Check for frameshift mutations or incorrect CDS boundaries"
                ))

        # Enhanced premature stop codon detection based on CDS modules
        for m in modules:
            module_type = m.get('module_type', '').lower()
            if 'cds' not in module_type:
                continue

            module_end = m.get('end', 0)
            module_start = m.get('start', 0)
            module_id = m.get('module_id', '')
            feature_ids = set(m.get('feature_ids', []))

            # Get CDS features in this module
            cds_features_in_module = [
                f for f in features
                if f.get('instance_id') in feature_ids
                and f.get('kb_feature_class', '').lower() in ('cds', 'protein_coding', 'coding_sequence')
            ]

            for cds_feat in cds_features_in_module:
                feature_end = cds_feat.get('end', 0)
                feature_name = cds_feat.get('feature_name', 'unknown')

                # Rule: CDS feature continues >3bp after CDS module end
                if feature_end > module_end + 3:
                    issues.append(LintIssue(
                        issue_id=f"lint_{generate_uuid()}",
                        rule_id="premature_stop_codon",
                        severity="error",
                        category="frame_continuity",
                        message=f"CDS '{feature_name}' annotation extends {feature_end - module_end}bp beyond detected stop codon",
                        feature_ids=[cds_feat.get('instance_id')],
                        module_ids=[module_id],
                        location={'module_end': module_end, 'feature_end': feature_end},
                        suggestion="Check for frameshift mutations, incorrect CDS boundaries, or missing intron annotations"
                    ))

            # Rule: Next feature within 100bp is also CDS (possible premature stop)
            downstream_cds = [
                f for f in features
                if f.get('start', 0) > module_end
                and f.get('start', 0) <= module_end + 100
                and f.get('kb_feature_class', '').lower() in ('cds', 'protein_coding', 'coding_sequence')
                and f.get('strand', 1) == m.get('strand', 1)  # Same strand
            ]

            if downstream_cds:
                next_cds = min(downstream_cds, key=lambda f: f.get('start', 0))
                gap = next_cds.get('start', 0) - module_end
                issues.append(LintIssue(
                    issue_id=f"lint_{generate_uuid()}",
                    rule_id="premature_stop_codon",
                    severity="warning",
                    category="frame_continuity",
                    message=f"CDS module ends {gap}bp before next CDS feature '{next_cds.get('feature_name', 'unknown')}' - possible frameshift or premature stop",
                    module_ids=[module_id],
                    feature_ids=[next_cds.get('instance_id')],
                    location={'module_end': module_end, 'next_cds_start': next_cds.get('start', 0)},
                    suggestion="Check if these CDS features should be continuous or if there's a frameshift mutation"
                ))

        return issues

    @classmethod
    def run_all_checks(
        cls,
        features: List[Dict[str, Any]],
        modules: List[Dict[str, Any]]
    ) -> List[LintIssue]:
        """Run all lint checks and return combined issues"""
        all_issues = []
        
        # Basic checks
        all_issues.extend(cls.check_orientation_violations(features, modules))
        all_issues.extend(cls.check_pol3_misuse(features, modules))
        all_issues.extend(cls.check_frame_continuity(features, modules))
        
        # Module-level checks
        linter = ModuleLinter()
        all_issues.extend(linter.run_all_checks(modules, features))
        
        return all_issues
