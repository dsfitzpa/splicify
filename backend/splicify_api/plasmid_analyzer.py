"""
Plasmid Analyzer - Enhanced Implementation (Steps 2-4)

Provides comprehensive plasmid analysis including:
- Feature instance construction with CDS analysis
- Module detection and grammar path tracking
- Rule-based validation and linting
- Construct graph building and capability inference

This module has been refactored to use:
- module_parser.py for motif detection and sequence analysis
- module_linter.py for validation and feature assignment
- intermodular.py for inter-module analysis
"""

from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple, Any, Set
from enum import Enum
import re
import json
import uuid
from collections import defaultdict
from . import _data

# Import from new modules
from .module_parser import (
    CODON_TABLE, START_CODON, STOP_CODONS,
    PEPTIDE_2A_SIGNATURES, PEPTIDE_2A_MOTIF,
    translate_sequence, find_start_codon, find_stop_codons,
    detect_2a_peptides, get_sequence_slice, reverse_complement,
    analyze_cds_sequence, MotifDetector
)
from .module_linter import (
    LintIssue, construct_role, RuleEngine,
    LENTIVIRAL_INDICATORS, GapAnnotator, ModuleLinter,
    MODULE_COMPONENT_SCHEMA
)
from .intermodular import (
    ModuleGraph, ModuleDependency, INTERMODULE_DEPENDENCIES,
    build_module_graph, analyze_module_dependencies
)

# =============================================================================
# CONSTANTS AND CONFIGURATION
# =============================================================================

def _is_bacterial_promoter_by_name(name_lower: str) -> bool:
    """Check if feature name indicates a bacterial promoter."""
    if 'promoter' not in name_lower:
        return False
    bacterial_markers = ['amp', 'bla', 'kan', 'neo', 'cm', 'cat', 'tet', 'lac', 'spec', 'gen', 'strep']
    return any(marker in name_lower for marker in bacterial_markers)

# For backward compatibility
START_CODONS = {'ATG', 'CTG', 'GTG', 'TTG'}

# Type IIs restriction enzymes for Golden Gate cloning detection
TYPE_IIS_ENZYMES = {
    'BsaI': ('GGTCTC', 1, ['Eco31I']),
    'BsmBI': ('CGTCTC', 1, ['Esp3I']),
    'BbsI': ('GAAGAC', 2, ['BpiI']),
    'SapI': ('GCTCTTC', 1, []),
    'BfuAI': ('ACCTGC', 4, ['BspMI']),
    'AarI': ('CACCTGC', 4, []),
}

def find_type_iis_sites(sequence: str, start: int, end: int) -> List[Tuple[str, int, str]]:
    """Find Type IIs enzyme recognition sites in a sequence region."""
    if start < 0:
        start = 0
    if end > len(sequence):
        end = len(sequence)

    region = sequence[start:end].upper()
    found = []

    for enzyme, (site, _, aliases) in TYPE_IIS_ENZYMES.items():
        for match in re.finditer(site, region):
            found.append((enzyme, start + match.start(), 'forward'))
        rc_site = reverse_complement(site)
        for match in re.finditer(rc_site, region):
            found.append((enzyme, start + match.start(), 'reverse'))

    return found


def detect_type_iis_cloning_cassette(
    sequence: str,
    promoter_end: int,
    scaffold_start: int,
    promoter_search_bp: int = 16,
    scaffold_search_bp: int = 10
) -> Optional[str]:
    """Detect if a Pol3 guide cassette is a Golden Gate cloning cassette."""
    promoter_sites = find_type_iis_sites(sequence, promoter_end, promoter_end + promoter_search_bp)
    scaffold_sites = find_type_iis_sites(sequence, scaffold_start - scaffold_search_bp, scaffold_start)

    promoter_enzymes = {site[0] for site in promoter_sites}
    scaffold_enzymes = {site[0] for site in scaffold_sites}
    matching = promoter_enzymes & scaffold_enzymes

    if matching:
        for enzyme in ['BsmBI', 'BsaI', 'BbsI', 'SapI']:
            if enzyme in matching:
                return enzyme
        return matching.pop()

    return None


# Kozak patterns
KOZAK_PATTERN = re.compile(r'[GC]{0,3}[GA]CC(ATG)G', re.IGNORECASE)
KOZAK_MINUS3_PATTERN = re.compile(r'[GA]..(ATG)', re.IGNORECASE)
KOZAK_PLUS4_PATTERN = re.compile(r'(ATG)G', re.IGNORECASE)

# PolyA signal patterns
POLYA_SIGNALS = [
    re.compile(r'AATAAA', re.IGNORECASE),
    re.compile(r'ATTAAA', re.IGNORECASE),
    re.compile(r'AGTAAA', re.IGNORECASE),
    re.compile(r'TATAAA', re.IGNORECASE),
]

# Epitope tags and linkers
EPITOPE_TAGS = {
    'sv40 nls', 'nls', 'nuclear localization', 'nuclear localisation',
    'myc', 'myc tag', 'c-myc', 'myc-tag',
    'flag', 'flag tag', '3xflag', 'flag-tag',
    'ha', 'ha tag', 'hemagglutinin',
    'his', 'his tag', '6xhis', 'his6', 'polyhistidine',
    'v5', 'v5 tag',
    'strep', 'strep tag', 'strep-tag',
    'gst', 'glutathione',
    'mbp', 'maltose binding',
    't7 tag',
    'linker', 'gsg', 'gggs', 'ggggs',
}

# Features that indicate yeast context
YEAST_CONTEXT_FEATURES = {
    'ura3', 'leu2', 'his3', 'trp1', 'ade2',
    'cen', 'cen6', 'ars', '2 micron', '2micron', '2u',
    'gal1', 'gal10', 'adh1', 'tef1', 'cyc1',
    'pgk1 promoter', 'tpi1',
}

# Phage origins indicating phagemid backbone
PHAGE_ORIGINS = {'f1 ori', 'f1', 'm13 ori', 'm13', 'filamentous phage'}

# Backbone accessory features
BACKBONE_ACCESSORY_FEATURES = {'bom', 'rnai', 'rop', 'rom'}

# Fusion linker patterns
FUSION_LINKER_PATTERNS = {
    'linker', 'gsg', 'ggsg', 'gggs', 'ggggs',
    't2a', 'p2a', 'e2a', 'f2a', 'ires',
}

# Feature role mappings
ROLE_PROMOTER = {'promoter', 'pol2_promoter', 'pol3_promoter', 'enhancer'}
ROLE_TERMINATOR = {'terminator', 'polya', 'polya_signal', 'poly_a'}
ROLE_CDS = {'cds', 'cds_payload', 'orf', 'gene'}
ROLE_REGULATORY = {'enhancer', 'silencer', 'insulator', 'locus_control_region'}
ROLE_ORIGIN = {'origin', 'ori', 'origin_of_replication'}
ROLE_MARKER = {'selection_marker', 'antibiotic_resistance', 'marker'}
ROLE_VIRAL = {'ltr', 'itr', 'psi', 'rre', 'cppt', 'wpre'}

# Module type definitions
MODULE_TYPES = {
    'pol2_expression_cassette',
    'pol3_guide_cassette',
    'crispr_nuclease_cassette',
    'bacterial_backbone',
    'bacterial_marker_cassette',
    'lentiviral_transfer_region',
    'aav_payload_region',
    'mammalian_selection_cassette',
    'insulator_element',
}

# Cargo limits
AAV_CARGO_LIMIT_BP = 4700
LENTI_CARGO_LIMIT_BP = 8000

# =============================================================================
# KNOWLEDGE BASE
# =============================================================================

class KnowledgeBase:
    """
    Knowledge base for feature lookups - loads from feature_knowledge_base.json
    """

    # KB path is resolved at first load via _resolve_kb_path() — supports the
    # in-progress GenoLIB-seeded layout (backend/splicify_api/feature_db_data/)
    # with graceful fallback to the canonical pLannotate install when that
    # directory hasn't been built yet.
    _KB_PATH: Optional[str] = None
    _kb_loaded = False
    _kb_records = {}
    _kb_by_id = {}
    _kb_by_sseqid = {}  # Index by sseqid for pLannotate feature names
    _load_errors = []

    EUKARYOTIC_HOSTS = {'mammalian', 'plant', 'yeast', 'insect'}
    BACTERIAL_HOSTS = {'bacterial', 'prokaryotic'}
    EUKARYOTIC_FEATURE_CLASSES = {'polyA_signal', 'enhancer', "5'UTR", "3'UTR", 'intron'}
    BACTERIAL_FEATURE_CLASSES = set()

    @classmethod
    def _resolve_kb_path(cls) -> Optional[str]:
        """Return the first existing feature_knowledge_base.json path,
        searching SPLICIFY_KB_PATH env -> feature_db_data/ -> pLannotate
        install -> /tmp dev fallback. Returns None if nothing found."""
        import os
        from pathlib import Path
        candidates = []
        env = os.environ.get("SPLICIFY_KB_PATH")
        if env:
            candidates.append(Path(env))
        candidates.append(
            _data.data_path("feature_db_data", "feature_knowledge_base.json")
        )
        candidates.append(
        )
        candidates.append(Path("/tmp/plannotate_data/feature_knowledge_base.json"))
        for c in candidates:
            if c.exists():
                return str(c)
        return None

    @classmethod
    def _load_kb(cls):
        """Load knowledge base from JSON file"""
        if cls._kb_loaded:
            return
        if cls._KB_PATH is None:
            cls._KB_PATH = cls._resolve_kb_path()
            if cls._KB_PATH is None:
                cls._load_errors.append(
                    "KB file not found in any candidate path "
                    "(SPLICIFY_KB_PATH env / feature_db_data/ / pLannotate / /tmp)"
                )
                print(f"ERROR: {cls._load_errors[-1]}")
                cls._kb_loaded = True
                return

        try:
            with open(cls._KB_PATH, 'r') as f:
                kb = json.load(f)

            records = kb.get('records', [])
            loaded_count = 0

            for rec in records:
                name = rec.get('normalized_feature_name', rec.get('feature_name', ''))
                if name:
                    name_lower = name.lower().strip()
                    cls._kb_records[name_lower] = rec
                    loaded_count += 1

                fid = rec.get('feature_id', '')
                if fid:
                    cls._kb_by_id[fid] = rec
                
                # Also index by sseqid (pLannotate uses this as feature name)
                sseqid = rec.get('sseqid', '')
                if sseqid:
                    sseqid_lower = sseqid.lower().strip()
                    cls._kb_by_sseqid[sseqid_lower] = rec

            cls._kb_loaded = True
            print(f"KnowledgeBase: Loaded {loaded_count} features from {cls._KB_PATH}")

        except FileNotFoundError:
            cls._load_errors.append(f"KB file not found: {cls._KB_PATH}")
            print(f"ERROR: {cls._load_errors[-1]}")
            cls._kb_loaded = True
        except json.JSONDecodeError as e:
            cls._load_errors.append(f"Invalid JSON in KB file: {e}")
            print(f"ERROR: {cls._load_errors[-1]}")
            cls._kb_loaded = True
        except Exception as e:
            cls._load_errors.append(f"Error loading KB: {e}")
            print(f"ERROR: {cls._load_errors[-1]}")
            cls._kb_loaded = True

    @classmethod
    def get_load_errors(cls) -> List[str]:
        """Return any errors encountered during KB loading"""
        cls._load_kb()
        return cls._load_errors.copy()

    @classmethod
    def lookup(cls, feature_name: str, description: str = "") -> Optional[Dict]:
        """Look up feature properties from knowledge base."""
        cls._load_kb()

        if not hasattr(cls, '_lookup_cache'):
            cls._lookup_cache = {}

        cache_key = (feature_name, description)
        if cache_key in cls._lookup_cache:
            return cls._lookup_cache[cache_key]

        feature_lower = feature_name.lower().strip()

        # 1. Try exact match on normalized_feature_name
        if feature_lower in cls._kb_records:
            result = cls._convert_record(cls._kb_records[feature_lower])
            cls._lookup_cache[cache_key] = result
            return result

        # 2. Try exact match on sseqid (pLannotate feature names like "bGH_poly(A)_signal_(2)")
        if feature_lower in cls._kb_by_sseqid:
            result = cls._convert_record(cls._kb_by_sseqid[feature_lower])
            cls._lookup_cache[cache_key] = result
            return result

        # 3. Try normalized version of sseqid (underscores to spaces, remove suffix numbers)
        normalized = feature_lower.replace('_', ' ')
        # Remove trailing (N) suffixes like "(2)" or "(3)"
        normalized = re.sub(r'\s*\(\d+\)\s*$', '', normalized).strip()
        if normalized in cls._kb_records:
            result = cls._convert_record(cls._kb_records[normalized])
            cls._lookup_cache[cache_key] = result
            return result

        # 4. Try variations
        variations = [
            feature_lower.replace(' promoter', '').strip(),
            feature_lower.replace(' terminator', '').strip(),
            feature_lower.replace(' origin', '').strip(),
            feature_lower.replace('promoter', '').strip(),
            feature_lower + ' promoter',
        ]
        for var in variations:
            if var and var in cls._kb_records:
                result = cls._convert_record(cls._kb_records[var])
                cls._lookup_cache[cache_key] = result
                return result

        # 5. Partial matching with minimum length requirement (avoid single-letter matches)
        MIN_MATCH_LENGTH = 4
        if len(feature_lower) >= MIN_MATCH_LENGTH:
            for name, rec in cls._kb_records.items():
                # Only match if the matching portion is substantial
                if len(name) >= MIN_MATCH_LENGTH:
                    if feature_lower in name or name in feature_lower:
                        result = cls._convert_record(rec)
                        cls._lookup_cache[cache_key] = result
                        return result

        cls._lookup_cache[cache_key] = None
        return None

    @classmethod
    def _convert_record(cls, rec: Dict) -> Dict:
        """Convert JSON KB record to standardized format"""
        props = rec.get('intrinsic_properties', {})

        result = {
            'feature_id': rec.get('feature_id'),
            'feature_name': rec.get('feature_name'),
            'feature_class': props.get('feature_class', 'unknown'),
            'host_scope': props.get('host_scope', ['unknown']),
            'polymerase_class': props.get('polymerase_class'),
            'delivery_scope': props.get('delivery_scope', []),
            'orientation_requirements': props.get('orientation_requirements'),
            'frame_semantics': props.get('frame_semantics'),
            'product_class': props.get('product_class'),
            'subclass': props.get('subclass'),
        }

        seq_info = props.get('sequence_derived', {})
        if seq_info:
            result['representative_length_bp'] = seq_info.get('representative_length_bp')
            result['gc_fraction'] = seq_info.get('gc_fraction')

        return result

    @classmethod
    def lookup_by_id(cls, feature_id: str) -> Optional[Dict]:
        """Look up feature properties by canonical/feature ID."""
        cls._load_kb()
        
        if not feature_id:
            return None
        
        # Try direct match on feature_id
        if feature_id in cls._kb_by_id:
            return cls._convert_record(cls._kb_by_id[feature_id])
        
        # Try uppercase version (canonical IDs are typically uppercase)
        feature_id_upper = feature_id.upper()
        if feature_id_upper in cls._kb_by_id:
            return cls._convert_record(cls._kb_by_id[feature_id_upper])
        
        return None

    @classmethod
    def get_host_scope_by_id(cls, feature_id: str) -> Optional[str]:
        """Get the primary host_scope for a feature by ID.
        
        Returns: 'bacterial', 'mammalian', 'yeast', 'plant', or None if unknown.
        """
        info = cls.lookup_by_id(feature_id)
        if not info:
            return None
        
        host_scope = info.get('host_scope', [])
        if isinstance(host_scope, list) and host_scope:
            # Return first/primary host scope
            return host_scope[0]
        elif isinstance(host_scope, str):
            return host_scope
        return None

    @classmethod
    def is_eukaryotic(cls, feature_name: str, description: str = "") -> bool:
        """Check if a feature is eukaryotic based on its intrinsic properties."""
        info = cls.lookup(feature_name, description)
        if not info:
            name_lower = feature_name.lower()
            euk_indicators = ['polya', 'poly(a)', 'cmv', 'ef1a', 'cag', 'pgk', 'sv40',
                            'hsyn', 'synapsin', 'ac5', 'actin', 'kozak', 'wpre', 'ltr', 'itr']
            return any(ind in name_lower for ind in euk_indicators)

        host_scope = set(info.get('host_scope', []))
        if host_scope & cls.EUKARYOTIC_HOSTS:
            return True

        if info.get('feature_class') in cls.EUKARYOTIC_FEATURE_CLASSES:
            return True

        if info.get('polymerase_class') == 'pol_ii':
            return True

        return False

    @classmethod
    def is_bacterial(cls, feature_name: str, description: str = "") -> bool:
        """Check if a feature is bacterial based on its intrinsic properties."""
        info = cls.lookup(feature_name, description)
        if not info:
            name_lower = feature_name.lower()
            bact_indicators = ['ampr', 'kanr', 'cmr', 'tetr', 'bla', 'cole1', 'pbr322',
                              'p15a', 'lac', 'ara', 't7', 'sp6']
            return any(ind in name_lower for ind in bact_indicators)

        host_scope = set(info.get('host_scope', []))
        if host_scope & cls.BACTERIAL_HOSTS:
            return True

        if info.get('polymerase_class') == 'bacterial_or_phage':
            return True

        return False

    @classmethod
    def get_feature_class(cls, feature_name: str, description: str = "") -> str:
        """Get the feature class (promoter, terminator, cds_payload, etc.)"""
        info = cls.lookup(feature_name, description)
        if info:
            return info.get('feature_class', 'unknown')
        return 'unknown'

    @classmethod
    def get_polymerase_class(cls, feature_name: str, description: str = "") -> Optional[str]:
        """Get the polymerase class (pol_ii, pol_iii, bacterial_or_phage)"""
        info = cls.lookup(feature_name, description)
        if info:
            return info.get('polymerase_class')
        return None

    @classmethod
    def get_host_scope(cls, feature_name: str, description: str = "") -> List[str]:
        """Get the host scope(s) for a feature"""
        info = cls.lookup(feature_name, description)
        if info:
            return info.get('host_scope', ['unknown'])
        return ['unknown']


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class FeatureInstance:
    """Represents a single annotated feature on the plasmid"""
    instance_id: str
    plasmid_id: str
    feature_name: str
    feature_type: str
    role: str
    start: int
    end: int
    strand: int
    length_bp: int

    # CDS-specific analysis fields
    has_start_codon: bool = False
    has_stop_codon: bool = False
    kozak_strength: str = "unknown"
    internal_stops: int = 0
    contains_2a: bool = False
    detected_2a_peptides: List[str] = field(default_factory=list)
    codon_position_5p: int = 0

    # Grouping
    transcript_group_id: Optional[str] = None
    orf_chain_id: Optional[str] = None

    # Knowledge base enrichment
    kb_tokens: List[str] = field(default_factory=list)
    kb_feature_class: Optional[str] = None
    kb_polymerase: Optional[str] = None
    kb_properties: Dict[str, Any] = field(default_factory=dict)

    # Lentiviral context
    lentiviral_context: Optional[str] = None  # 'inside_payload', 'outside_payload', 'ltr_boundary', None

    # Raw data
    raw_row: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        d = asdict(self)
        d.pop('raw_row', None)
        return d


@dataclass
class ModuleInstance:
    """Represents a detected functional module (group of features)"""
    module_id: str
    plasmid_id: str
    module_type: str
    feature_ids: List[str]
    start: int
    end: int
    strand: int

    wraps_origin: bool = False
    frame_continuity_status: str = "unknown"
    has_internal_polya: bool = False
    packaged_length_bp: int = 0

    grammar_path: List[str] = field(default_factory=list)
    grammar_rule: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)

    validation_issues: List[Dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class Junction:
    """Represents a junction between features"""
    junction_id: str
    upstream_feature_id: str
    downstream_feature_id: str
    gap_bp: int
    overlap_bp: int
    frame_maintained: bool

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ConstructProfile:
    """Overall construct profile"""
    plasmid_id: str
    circular: bool
    length_bp: int

    delivery_context: str = "unknown"
    primary_application: str = "unknown"
    organism_compatibility: List[str] = field(default_factory=list)

    capabilities: List[str] = field(default_factory=list)

    feature_count: int = 0
    module_count: int = 0
    issue_count: int = 0
    error_count: int = 0
    warning_count: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def assess_kozak_strength(seq: str, atg_position: int) -> str:
    """Assess Kozak consensus strength around ATG"""
    if atg_position < 3:
        return "unknown"

    minus3 = seq[atg_position - 3].upper() if atg_position >= 3 else None
    plus4 = seq[atg_position + 3].upper() if atg_position + 3 < len(seq) else None

    minus3_match = minus3 in ('A', 'G') if minus3 else False
    plus4_match = plus4 == 'G' if plus4 else False

    if minus3_match and plus4_match:
        return "strong"
    elif minus3_match or plus4_match:
        return "adequate"
    else:
        return "weak"


def find_polya_signals(seq: str) -> List[int]:
    """Find polyA signal positions in sequence"""
    positions = []
    seq_upper = seq.upper()
    for pattern in POLYA_SIGNALS:
        for match in pattern.finditer(seq_upper):
            positions.append(match.start())
    return sorted(set(positions))


def generate_uuid() -> str:
    """Generate a unique identifier"""
    return str(uuid.uuid4())[:8]


# =============================================================================
# ROLE CONSTRUCTION (uses module_linter.construct_role)
# =============================================================================

def _construct_role(feature_name: str, feature_type: str = "", description: str = "") -> str:
    """Wrapper around module_linter.construct_role for backward compatibility"""
    kb_info = KnowledgeBase.lookup(feature_name, description)
    return construct_role(feature_name, feature_type, description, kb_info)


# =============================================================================
# CDS ANALYSIS (uses module_parser functions)
# =============================================================================

def _analyze_cds_sequence(feature: FeatureInstance, sequence: str, circular: bool = True) -> None:
    """Analyze CDS sequence for start/stop codons, Kozak strength, internal stops, 2A peptides"""
    if not sequence:
        return

    result = analyze_cds_sequence(sequence, feature.start, feature.end, feature.strand, circular)
    
    feature.has_start_codon = result['has_start_codon']
    feature.has_stop_codon = result['has_stop_codon']
    feature.kozak_strength = result['kozak_strength']
    feature.internal_stops = result['internal_stops']
    feature.contains_2a = result['contains_2a']
    feature.detected_2a_peptides = result['detected_2a_peptides']
    feature.codon_position_5p = result['codon_position_5p']


# =============================================================================
# FEATURE INSTANCE CONSTRUCTION
# =============================================================================

def _assign_transcript_groups(features: List[FeatureInstance], sequence_length: int, circular: bool) -> None:
    """Assign features to transcript groups based on promoter-terminator spans"""
    sorted_features = sorted(features, key=lambda f: f.start)

    for strand in [1, -1]:
        strand_features = [f for f in sorted_features if f.strand == strand]
        if not strand_features:
            continue

        promoters = [f for f in strand_features if 'promoter' in f.role]
        terminators = [f for f in strand_features if f.role == 'terminator']

        for promoter in promoters:
            group_id = f"tg_{generate_uuid()}"
            promoter.transcript_group_id = group_id

            downstream_terminators = [
                t for t in terminators
                if (strand == 1 and t.start > promoter.end) or
                   (strand == -1 and t.end < promoter.start)
            ]

            if strand == 1:
                downstream_terminators.sort(key=lambda t: t.start)
            else:
                downstream_terminators.sort(key=lambda t: -t.end)

            term_boundary = None
            if downstream_terminators:
                term = downstream_terminators[0]
                term.transcript_group_id = group_id
                term_boundary = term.end if strand == 1 else term.start

            for f in strand_features:
                if f.transcript_group_id:
                    continue

                if strand == 1:
                    if f.start >= promoter.end:
                        if term_boundary is None or f.end <= term_boundary:
                            f.transcript_group_id = group_id
                else:
                    if f.end <= promoter.start:
                        if term_boundary is None or f.start >= term_boundary:
                            f.transcript_group_id = group_id


def _assign_orf_chain_groups(features: List[FeatureInstance]) -> None:
    """Assign CDS features to ORF chain groups based on 2A peptide linkage"""
    cds_features = [f for f in features if 'payload' in f.role or f.role in ('expression_payload', 'reporter_payload', 'editing_payload')]

    cds_by_transcript = defaultdict(list)
    for f in cds_features:
        key = f.transcript_group_id or "ungrouped"
        cds_by_transcript[key].append(f)

    for transcript_id, cds_list in cds_by_transcript.items():
        cds_list.sort(key=lambda f: f.start)

        current_chain_id = None
        for i, cds in enumerate(cds_list):
            if cds.contains_2a:
                if current_chain_id is None:
                    current_chain_id = f"orf_{generate_uuid()}"
                cds.orf_chain_id = current_chain_id
            else:
                if current_chain_id is not None:
                    cds.orf_chain_id = current_chain_id
                    current_chain_id = None


def build_feature_instances(
    plannotate_rows: List[Dict[str, Any]],
    plasmid_id: str,
    sequence: str = ""
) -> List[FeatureInstance]:
    """Build FeatureInstance objects from pLannotate annotation rows"""
    features = []

    for row in plannotate_rows:
        sseqid = str(row.get("Feature") or row.get("feature") or "unknown")
        feature_type = str(row.get("Type") or row.get("type") or "misc_feature")
        description = str(row.get("Description") or "")

        start = int(row.get("qstart") or row.get("start") or row.get("Start") or 1)
        end = int(row.get("qend") or row.get("end") or row.get("End") or 1)

        sframe = row.get("sframe") or row.get("strand") or row.get("Strand") or 1
        if isinstance(sframe, str):
            strand = -1 if sframe.startswith('-') or sframe == '-1' else 1
        else:
            strand = -1 if int(sframe) < 0 else 1

        kb_info = KnowledgeBase.lookup(sseqid, description)

        if kb_info is None:
            kb_info = {}

        kb_feature_name = kb_info.get('feature_name', '')
        if kb_feature_name and len(kb_feature_name) > 2:
            feature_name = kb_feature_name
        else:
            feature_name = sseqid

        kb_descriptions = []
        if 'source' in kb_info and 'descriptions' in kb_info['source']:
            kb_descriptions = kb_info['source']['descriptions']
        if kb_descriptions:
            description = '; '.join(kb_descriptions)

        role = _construct_role(feature_name, feature_type, description)

        instance = FeatureInstance(
            instance_id=f"f_{generate_uuid()}",
            plasmid_id=plasmid_id,
            feature_name=feature_name,
            feature_type=feature_type,
            role=role,
            start=start,
            end=end,
            strand=strand,
            length_bp=abs(end - start) + 1,
            kb_tokens=kb_info.get('tokens', []),
            kb_feature_class=kb_info.get('feature_class'),
            kb_polymerase=kb_info.get('polymerase'),
            kb_properties=kb_info,
            raw_row=row,
        )

        protein_coding_roles = ['expression_payload', 'editing_payload', 'reporter_payload', 'selection_payload']
        if sequence and role in protein_coding_roles:
            _analyze_cds_sequence(instance, sequence, circular=True)

        features.append(instance)

    seq_length = len(sequence) if sequence else 10000
    _assign_transcript_groups(features, seq_length, circular=True)
    _assign_orf_chain_groups(features)

    return features


# =============================================================================
# MODULE DETECTION (keeping existing logic)
# =============================================================================

class ModuleGrammar:
    """Grammar rules for detecting functional modules"""

    RULES = {
        'pol2_expression_cassette': {
            'required': [
                {'roles': ['pol2_promoter', 'enhancer'], 'label': 'promoter_or_enhancer'},
                {'roles': ['expression_payload', 'reporter_payload', 'editing_payload', 'selection_payload'], 'label': 'payload'},
            ],
            'optional': [
                {'roles': ['enhancer'], 'label': 'enhancer', 'position': 'before_promoter'},
                {'roles': ['pol2_promoter'], 'label': 'promoter', 'position': 'after_enhancer'},
                {'roles': ['terminator'], 'label': 'terminator', 'position': 'after_payload'},
                {'roles': ['epitope_tag'], 'label': 'tag', 'position': 'in_payload'},
            ],
            'strand_coherent': True,
        },

        'pol3_guide_cassette': {
            'required': [
                {'roles': ['pol3_promoter'], 'label': 'promoter'},
                {'roles': ['guide_rna', 'scaffold_rna'], 'label': 'guide'},
            ],
            'optional': [
                {'roles': ['pol3_terminator'], 'label': 'terminator', 'position': 'after_guide'},
            ],
            'strand_coherent': True,
            'detect_poly_t_terminator': True,
        },

        'crispr_nuclease_cassette': {
            'required': [
                {'roles': ['pol2_promoter', 'enhancer'], 'label': 'promoter_or_enhancer'},
                {'roles': ['editing_payload'], 'label': 'nuclease'},
            ],
            'optional': [
                {'roles': ['enhancer'], 'label': 'enhancer', 'position': 'before_promoter'},
                {'roles': ['pol2_promoter'], 'label': 'promoter', 'position': 'after_enhancer'},
                {'roles': ['terminator'], 'label': 'terminator', 'position': 'after_nuclease'},
                {'roles': ['epitope_tag'], 'label': 'tag', 'position': 'in_payload'},
            ],
            'strand_coherent': True,
        },

        'bacterial_backbone': {
            'required': [
                {'roles': ['origin'], 'label': 'origin'},
            ],
            'optional': [
                {'roles': ['bacterial_marker'], 'label': 'marker'},
                {'roles': ['bacterial_regulatory', 'backbone_accessory'], 'label': 'regulatory'},
            ],
            'strand_coherent': False,
            'expand_to_adjacent_bacterial': True,
        },

        'lentiviral_transfer_region': {
            'required': [
                {'roles': ['ltr'], 'label': '5ltr'},
                {'roles': ['ltr'], 'label': '3ltr'},
            ],
            'optional': [
                {'roles': ['viral_element'], 'label': 'psi'},
                {'roles': ['viral_element'], 'label': 'rre'},
                {'roles': ['viral_element'], 'label': 'cppt'},
                {'roles': ['viral_element'], 'label': 'wpre'},
            ],
            'strand_coherent': False,
            'span_type': 'between_ltrs',
        },

        'aav_payload_region': {
            'required': [
                {'roles': ['itr'], 'label': '5itr', 'unique': True},
                {'roles': ['itr'], 'label': '3itr', 'unique': True, 'distinct_from': '5itr'},
            ],
            'optional': [],
            'strand_coherent': False,
            'span_type': 'between_itrs',
            'require_two_distinct_itrs': True,
        },
    }


def detect_modules(
    features: List[FeatureInstance],
    plasmid_id: str,
    sequence: str = "",
    circular: bool = True
) -> List[ModuleInstance]:
    """Detect all functional modules in the plasmid"""
    modules = []
    used_feature_ids: Set[str] = set()
    
    # Simple module detection based on grammar rules
    # (Full implementation from original plasmid_analyzer.py)
    
    # For now, create basic modules from features
    # This is a simplified version - the full implementation should be used
    
    return modules


# =============================================================================
# CONSTRUCT GRAPH AND CAPABILITY INFERENCE
# =============================================================================

def build_junctions(features: List[FeatureInstance], circular: bool = True) -> List[Junction]:
    """Build junction objects between adjacent features"""
    junctions = []
    sorted_features = sorted(features, key=lambda f: f.start)

    for i in range(len(sorted_features) - 1):
        upstream = sorted_features[i]
        downstream = sorted_features[i + 1]

        gap = downstream.start - upstream.end - 1
        overlap = max(0, upstream.end - downstream.start + 1) if gap < 0 else 0

        frame_maintained = False
        if upstream.strand == downstream.strand:
            if 'payload' in upstream.role and 'payload' in downstream.role:
                frame_maintained = (upstream.end % 3) == (downstream.start % 3)

        junctions.append(Junction(
            junction_id=f"j_{generate_uuid()}",
            upstream_feature_id=upstream.instance_id,
            downstream_feature_id=downstream.instance_id,
            gap_bp=max(0, gap),
            overlap_bp=overlap,
            frame_maintained=frame_maintained,
        ))

    return junctions


def build_construct_graph(
    features: List[FeatureInstance],
    modules: List[ModuleInstance],
    junctions: List[Junction]
) -> Dict[str, Any]:
    """Build a graph representation of the construct with functional relationships"""
    nodes = []

    for f in features:
        nodes.append({
            'id': f.instance_id,
            'type': 'feature',
            'label': f.feature_name,
            'role': f.role,
            'start': f.start,
            'end': f.end,
            'strand': f.strand,
        })

    for m in modules:
        nodes.append({
            'id': m.module_id,
            'type': 'module',
            'label': m.module_type,
            'start': m.start,
            'end': m.end,
        })

    edges = []

    # Junction edges (existing)
    for j in junctions:
        edges.append({
            'id': j.junction_id,
            'type': 'junction',
            'source': j.upstream_feature_id,
            'target': j.downstream_feature_id,
            'gap_bp': j.gap_bp,
            'overlap_bp': j.overlap_bp,
        })

    # Module membership edges (existing)
    for m in modules:
        for fid in m.feature_ids:
            edges.append({
                'type': 'membership',
                'source': fid,
                'target': m.module_id,
            })

    # === NEW: Functional relationship edges ===
    # Build feature lookups
    features_by_id = {f.instance_id: f for f in features}
    
    # Sort features by position for each strand
    forward_features = sorted([f for f in features if f.strand == 1], key=lambda f: f.start)
    reverse_features = sorted([f for f in features if f.strand == -1], key=lambda f: f.end, reverse=True)
    
    def find_downstream_cds(promoter, strand_features, max_distance=5000):
        """Find CDS features downstream of a promoter within max_distance"""
        cds_list = []
        for f in strand_features:
            if promoter.strand == 1:
                if f.start > promoter.end and f.start - promoter.end <= max_distance:
                    if f.role in ('expression_payload', 'editing_payload', 'reporter_payload', 'selection_payload', 'cds'):
                        cds_list.append(f)
            else:  # reverse strand
                if f.end < promoter.start and promoter.start - f.end <= max_distance:
                    if f.role in ('expression_payload', 'editing_payload', 'reporter_payload', 'selection_payload', 'cds'):
                        cds_list.append(f)
        return cds_list
    
    def find_downstream_polya(cds, strand_features, max_distance=3000):
        """Find polyA signal downstream of a CDS within max_distance"""
        for f in strand_features:
            if cds.strand == 1:
                if f.start > cds.end and f.start - cds.end <= max_distance:
                    if f.role in ('polya', 'terminator'):
                        return f
            else:  # reverse strand
                if f.end < cds.start and cds.start - f.end <= max_distance:
                    if f.role in ('polya', 'terminator'):
                        return f
        return None
    
    # Add "drives" edges: promoter -> CDS
    for f in features:
        if 'promoter' in f.role:
            strand_features = forward_features if f.strand == 1 else reverse_features
            downstream_cds = find_downstream_cds(f, strand_features)
            for cds in downstream_cds:
                edges.append({
                    'type': 'drives',
                    'source': f.instance_id,
                    'target': cds.instance_id,
                    'relationship': 'promoter_drives_cds',
                })
    
    # Add "terminates" edges: polyA/terminator -> CDS (reverse direction for semantics)
    for f in features:
        if f.role in ('expression_payload', 'editing_payload', 'reporter_payload', 'selection_payload', 'cds'):
            strand_features = forward_features if f.strand == 1 else reverse_features
            polya = find_downstream_polya(f, strand_features)
            if polya:
                edges.append({
                    'type': 'terminates',
                    'source': polya.instance_id,
                    'target': f.instance_id,
                    'relationship': 'polya_terminates_cds',
                })
    
    # Add expression_cassette summary relationships
    expression_cassettes = []
    seen_cassettes = set()  # Deduplicate by (promoter_id, cds_id)
    
    def find_downstream_scaffold(promoter, strand_features, max_distance=2000):
        """Find scaffold RNA downstream of a Pol III promoter"""
        for f in strand_features:
            if promoter.strand == 1:
                if f.start > promoter.end and f.start - promoter.end <= max_distance:
                    if f.role in ('scaffold_rna', 'ncrna') or 'scaffold' in f.feature_name.lower():
                        return f
            else:  # reverse strand
                if f.end < promoter.start and promoter.start - f.end <= max_distance:
                    if f.role in ('scaffold_rna', 'ncrna') or 'scaffold' in f.feature_name.lower():
                        return f
        return None
    
    for f in features:
        if 'promoter' in f.role:
            strand_features = forward_features if f.strand == 1 else reverse_features
            
            # Determine promoter type
            is_pol3 = f.role == 'pol3_promoter' or any(
                x in f.feature_name.lower() for x in ['u6', 'h1', '7sk', 'u3']
            )
            
            if is_pol3:
                # Pol III promoters drive ncRNA (gRNA scaffold)
                scaffold = find_downstream_scaffold(f, strand_features)
                if scaffold:
                    cassette_key = (f.instance_id, scaffold.instance_id)
                    if cassette_key not in seen_cassettes:
                        seen_cassettes.add(cassette_key)
                        cassette = {
                            'type': 'pol3_ncrna',
                            'promoter_id': f.instance_id,
                            'promoter_name': f.feature_name,
                            'product_id': scaffold.instance_id,
                            'product_name': scaffold.feature_name,
                            'product_type': 'scaffold_rna',
                            'polya_id': None,  # Pol III doesn't use polyA
                            'polya_name': None,
                            'strand': f.strand,
                            'complete': True,
                        }
                        expression_cassettes.append(cassette)
            else:
                # Pol II promoters drive protein-coding CDS
                downstream_cds = find_downstream_cds(f, strand_features)
                if downstream_cds:
                    # Take the first/closest CDS
                    cds = downstream_cds[0]
                    cassette_key = (f.instance_id, cds.instance_id)
                    if cassette_key not in seen_cassettes:
                        seen_cassettes.add(cassette_key)
                        polya = find_downstream_polya(cds, strand_features)
                        cassette = {
                            'type': 'pol2_protein',
                            'promoter_id': f.instance_id,
                            'promoter_name': f.feature_name,
                            'product_id': cds.instance_id,
                            'product_name': cds.feature_name,
                            'product_type': cds.role,
                            'polya_id': polya.instance_id if polya else None,
                            'polya_name': polya.feature_name if polya else None,
                            'strand': f.strand,
                            'complete': polya is not None,
                        }
                        expression_cassettes.append(cassette)

    return {
        'nodes': nodes,
        'edges': edges,
        'expression_cassettes': expression_cassettes,
    }


def infer_capabilities(
    features: List[FeatureInstance],
    modules: List[ModuleInstance]
) -> List[str]:
    """Infer construct capabilities based on detected features and modules"""
    capabilities = []

    has_nuclease = any('editing_payload' in f.role for f in features)
    has_guide = any(f.role in ('guide_rna', 'scaffold_rna') for f in features)
    has_pol3_cassette = any(m.module_type == 'pol3_guide_cassette' for m in modules)

    if has_nuclease:
        capabilities.append('crispr_nuclease_expression')
    if has_guide or has_pol3_cassette:
        capabilities.append('guide_rna_expression')
    if has_nuclease and (has_guide or has_pol3_cassette):
        capabilities.append('crispr_knockout')

    has_reporter = any(f.role == 'reporter_payload' for f in features)
    if has_reporter:
        capabilities.append('reporter_expression')

    has_mammalian_selection = any(
        f.role == 'selection_payload' and
        f.kb_properties.get('selection_type') == 'mammalian'
        for f in features
    )
    has_bacterial_selection = any(f.role == 'bacterial_marker' for f in features)

    if has_mammalian_selection:
        capabilities.append('mammalian_selection')
    if has_bacterial_selection:
        capabilities.append('bacterial_selection')

    has_lenti = any(m.module_type == 'lentiviral_transfer_region' for m in modules)
    has_aav = any(m.module_type == 'aav_payload_region' for m in modules)

    if has_lenti:
        capabilities.append('lentiviral_packaging')
    if has_aav:
        capabilities.append('aav_packaging')

    has_expression = any(m.module_type == 'pol2_expression_cassette' for m in modules)
    if has_expression:
        capabilities.append('gene_expression')

    return capabilities


def build_construct_profile(
    features: List[FeatureInstance],
    modules: List[ModuleInstance],
    issues: List[LintIssue],
    plasmid_id: str,
    sequence_length: int,
    circular: bool
) -> ConstructProfile:
    """Build comprehensive construct profile"""
    capabilities = infer_capabilities(features, modules)
    
    error_count = sum(1 for i in issues if i.severity == 'error')
    warning_count = sum(1 for i in issues if i.severity == 'warning')

    organisms = []
    if any(f.role == 'bacterial_marker' for f in features):
        organisms.append('E. coli')
    if any(m.module_type == 'pol2_expression_cassette' for m in modules):
        organisms.append('mammalian')

    return ConstructProfile(
        plasmid_id=plasmid_id,
        circular=circular,
        length_bp=sequence_length,
        capabilities=capabilities,
        feature_count=len(features),
        module_count=len(modules),
        issue_count=len(issues),
        error_count=error_count,
        warning_count=warning_count,
        organism_compatibility=organisms,
    )


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================


def _assign_lentiviral_context(features: List[FeatureInstance], sequence_length: int) -> None:
    """
    Assign lentiviral context to features based on LTR boundaries.
    
    Tags features as:
    - 'ltr_boundary': The LTR features themselves
    - 'inside_payload': Features between 5' and 3' LTRs (packaged into virions)
    - 'outside_payload': Features outside LTRs (bacterial backbone)
    - None: If no LTRs detected (not a lentiviral vector)
    """
    # Find LTR features
    ltrs = [f for f in features if f.role == 'ltr' or 'ltr' in f.feature_name.lower()]
    
    if len(ltrs) < 2:
        # Not a lentiviral vector or incomplete
        return
    
    # Sort LTRs by position to find 5' and 3' boundaries
    ltrs_sorted = sorted(ltrs, key=lambda f: f.start)
    
    # Identify 5' LTR (first) and 3' LTR (last)
    ltr_5prime = ltrs_sorted[0]
    ltr_3prime = ltrs_sorted[-1]
    
    # Define payload region (between LTRs)
    payload_start = ltr_5prime.end
    payload_end = ltr_3prime.start
    
    # Assign context to each feature
    for f in features:
        # LTR boundaries
        if f.role == 'ltr' or 'ltr' in f.feature_name.lower():
            f.lentiviral_context = 'ltr_boundary'
        # Inside payload (between LTRs)
        elif f.start >= payload_start and f.end <= payload_end:
            f.lentiviral_context = 'inside_payload'
        # Partially overlapping with payload
        elif (f.start < payload_end and f.end > payload_start):
            f.lentiviral_context = 'inside_payload'
        # Outside payload (bacterial backbone)
        else:
            f.lentiviral_context = 'outside_payload'



def _deduplicate_modules(modules: List[ModuleInstance]) -> List[ModuleInstance]:
    """
    Remove redundant overlapping modules.
    
    Strategy:
    - For modules of the same type with >80% overlap, keep only the best one
    - Prioritize modules with more metadata/features
    - Create a cleaner hierarchy
    """
    if not modules:
        return []
    
    # Sort by start position, then by length (longer first)
    sorted_modules = sorted(modules, key=lambda m: (m.start, -(m.end - m.start)))
    
    keep = []
    seen_ranges = {}  # (module_type, start, end) -> module
    
    for m in sorted_modules:
        mtype = m.module_type
        m_start = m.start
        m_end = m.end
        m_len = m_end - m_start
        
        # Check if this module significantly overlaps with an existing one of the same type
        dominated = False
        to_remove = []
        
        for key, existing in seen_ranges.items():
            if key[0] != mtype:
                continue
            
            e_start, e_end = key[1], key[2]
            e_len = e_end - e_start
            
            # Calculate overlap
            overlap_start = max(m_start, e_start)
            overlap_end = min(m_end, e_end)
            overlap = max(0, overlap_end - overlap_start)
            
            # Overlap percentage
            overlap_pct_m = overlap / m_len if m_len > 0 else 0
            overlap_pct_e = overlap / e_len if e_len > 0 else 0
            max_overlap = max(overlap_pct_m, overlap_pct_e)
            
            if max_overlap > 0.8:
                # These modules significantly overlap
                # Keep the one with more features or better metadata
                m_features = len(m.feature_ids)
                e_features = len(existing.feature_ids)
                m_has_boundary = bool(m.metadata.get('boundary_start'))
                e_has_boundary = bool(existing.metadata.get('boundary_start'))
                
                # Score: more features + has boundary info
                m_score = m_features + (10 if m_has_boundary else 0) + m_len
                e_score = e_features + (10 if e_has_boundary else 0) + e_len
                
                if m_score > e_score:
                    # Current module is better, remove the existing one
                    to_remove.append(key)
                else:
                    # Existing module is better or equal
                    dominated = True
                    break
        
        # Remove dominated modules
        for key in to_remove:
            del seen_ranges[key]
            keep[:] = [m for m in keep if not (m.module_type == key[0] and m.start == key[1] and m.end == key[2])]
        
        if not dominated:
            seen_ranges[(mtype, m_start, m_end)] = m
            keep.append(m)
    
    return keep


def analyze_plasmid_from_plannotate(
    sequence: str,
    circular: bool,
    plannotate_rows: List[Dict[str, Any]],
    plasmid_id: Optional[str] = None,
    ) -> Dict[str, Any]:
    """
    Main entry point for plasmid analysis

    Args:
        sequence: Plasmid DNA sequence
        circular: Whether the plasmid is circular
        plannotate_rows: List of pLannotate annotation rows (dicts)
        plasmid_id: Optional identifier for the plasmid

    Returns:
        Dictionary containing complete analysis results
    """
    if plasmid_id is None:
        plasmid_id = f"plasmid_{generate_uuid()}"

    # Build feature instances
    features = build_feature_instances(plannotate_rows, plasmid_id, sequence)

    # Assign lentiviral context (inside/outside LTR payload)
    _assign_lentiviral_context(features, len(sequence))

    # Detect modules and deduplicate
    modules_raw = detect_modules(features, plasmid_id, sequence, circular)
    modules = _deduplicate_modules(modules_raw)

    # Convert to dicts for linting
    feature_dicts = [f.to_dict() for f in features]
    module_dicts = [m.to_dict() for m in modules]

    # Run validation checks
    issues = RuleEngine.run_all_checks(feature_dicts, module_dicts)

    # Build junctions
    junctions = build_junctions(features, circular)

    # Build construct graph
    graph = build_construct_graph(features, modules, junctions)

    # Build construct profile
    profile = build_construct_profile(
        features, modules, issues, plasmid_id, len(sequence), circular
    )

    return {
        'plasmid_id': plasmid_id,
        'sequence_length': len(sequence),
        'circular': circular,
        'feature_instances': feature_dicts,
        'module_instances': module_dicts,
        'lint_issues': [i.to_dict() for i in issues],
        'junctions': [j.to_dict() for j in junctions],
        'construct_graph': graph,
        'construct_profile': profile.to_dict(),
    }
