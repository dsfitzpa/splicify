"""
Intermodular - Module Graph and Dependencies

Provides:
- Inter-module dependency system
- Module graph builder
- Module feature annotations for GenBank
- Export for Intent Analysis workflow

This module is part of the annotation pipeline restructuring to create
cleaner separation between motif/boundary detection, validation/linting,
and inter-module analysis.
"""

from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple, Any, Set
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ModuleDependency:
    """Represents a dependency relationship between module types"""
    source_module_type: str
    target_module_type: str
    relationship: str  # 'requires', 'prefers', 'conflicts_with'
    condition: Optional[str] = None  # Optional condition for when this applies
    description: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class DependencyIssue:
    """Represents a failed dependency check"""
    source_module_id: str
    source_module_type: str
    target_module_type: str
    relationship: str
    message: str
    severity: str = "warning"  # 'error' or 'warning'

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class PreferenceWarning:
    """Represents a module preference that isn't satisfied"""
    source_module_id: str
    source_module_type: str
    preferred_module_type: str
    message: str
    suggestion: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ModuleNode:
    """Represents a module as a node in the graph"""
    module_id: str
    module_type: str
    start: int
    end: int
    strand: int
    feature_ids: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ModuleEdge:
    """Represents a relationship between modules"""
    source_id: str
    target_id: str
    relationship: str  # 'requires', 'prefers', 'conflicts_with', 'contains', 'adjacent'
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# =============================================================================
# INTER-MODULE DEPENDENCIES
# =============================================================================

INTERMODULE_DEPENDENCIES = [
    # Shuttle vector needs yeast module
    ModuleDependency(
        'shuttle_vector_backbone',
        'yeast_pol2_expression',
        'prefers',
        description='Shuttle vectors benefit from yeast expression capability'
    ),
    ModuleDependency(
        'f1_ori',
        'yeast_pol2_expression',
        'prefers',
        condition='with_shuttle_vector',
        description='F1 origin with shuttle vector typically indicates yeast/bacterial dual use'
    ),

    # Bacterial replication requirements
    ModuleDependency(
        'bacterial_marker_cassette',
        'bacterial_backbone',
        'requires',
        description='Bacterial selection marker requires bacterial origin for propagation'
    ),

    # Viral packaging requirements
    ModuleDependency(
        'lentiviral_payload',
        'lentiviral_cis_module',
        'requires',
        description='Lentiviral payload requires cis-acting elements (Psi, RRE) for packaging'
    ),
    
    # Expression cassette requirements
    ModuleDependency(
        'pol2_expression_cassette',
        'polya_signal',
        'requires',
        description='Pol II expression requires polyA signal for transcript stability'
    ),
    
    # CRISPR requirements
    ModuleDependency(
        'pol3_guide_cassette',
        'crispr_nuclease_cassette',
        'prefers',
        description='gRNA expression typically paired with nuclease expression'
    ),
    
    # AAV packaging
    ModuleDependency(
        'aav_payload',
        'aav_itr_pair',
        'requires',
        description='AAV payload requires two ITRs flanking the transgene'
    ),
]


# =============================================================================
# MODULE GRAPH CLASS
# =============================================================================

class ModuleGraph:
    """Represent plasmid as interconnected modules"""

    def __init__(
        self,
        modules: List[Dict[str, Any]],
        features: List[Dict[str, Any]],
        sequence_length: int,
        circular: bool = True
    ):
        self.modules = modules
        self.features = features
        self.seq_len = sequence_length
        self.circular = circular
        
        self.nodes: Dict[str, ModuleNode] = {}
        self.edges: List[ModuleEdge] = []
        
        # Build initial graph
        self._build_nodes()

    def _build_nodes(self):
        """Build module nodes from module list"""
        for mod in self.modules:
            node = ModuleNode(
                module_id=mod.get('module_id', ''),
                module_type=mod.get('module_type', ''),
                start=mod.get('start', 0),
                end=mod.get('end', 0),
                strand=mod.get('strand', 1),
                feature_ids=mod.get('feature_ids', []),
                metadata=mod.get('metadata', {}),
            )
            self.nodes[node.module_id] = node

    def build_graph(self):
        """Build module relationship graph"""
        # Clear existing edges
        self.edges = []
        
        # Build adjacency relationships
        self._build_adjacency_edges()
        
        # Build containment relationships
        self._build_containment_edges()
        
        # Build dependency edges
        self._build_dependency_edges()
        
        return self

    def _build_adjacency_edges(self):
        """Build edges for adjacent modules"""
        # Sort modules by start position
        sorted_modules = sorted(self.nodes.values(), key=lambda n: n.start)
        
        for i in range(len(sorted_modules) - 1):
            current = sorted_modules[i]
            next_mod = sorted_modules[i + 1]
            
            # Check if adjacent (gap < 100bp)
            gap = next_mod.start - current.end
            if 0 <= gap < 100:
                self.edges.append(ModuleEdge(
                    source_id=current.module_id,
                    target_id=next_mod.module_id,
                    relationship='adjacent',
                    metadata={'gap_bp': gap}
                ))

    def _build_containment_edges(self):
        """Build edges for modules that contain other modules"""
        for outer_id, outer in self.nodes.items():
            for inner_id, inner in self.nodes.items():
                if outer_id == inner_id:
                    continue
                
                # Check if inner is fully contained within outer
                if inner.start >= outer.start and inner.end <= outer.end:
                    self.edges.append(ModuleEdge(
                        source_id=outer_id,
                        target_id=inner_id,
                        relationship='contains',
                    ))

    def _build_dependency_edges(self):
        """Build edges based on dependency rules"""
        module_types = {n.module_type for n in self.nodes.values()}
        
        for dep in INTERMODULE_DEPENDENCIES:
            if dep.source_module_type in module_types:
                if dep.target_module_type in module_types:
                    # Find matching modules
                    source_mods = [
                        n for n in self.nodes.values()
                        if n.module_type == dep.source_module_type
                    ]
                    target_mods = [
                        n for n in self.nodes.values()
                        if n.module_type == dep.target_module_type
                    ]
                    
                    for source in source_mods:
                        for target in target_mods:
                            self.edges.append(ModuleEdge(
                                source_id=source.module_id,
                                target_id=target.module_id,
                                relationship=dep.relationship,
                                metadata={'dependency_description': dep.description}
                            ))

    def check_dependencies(self) -> List[DependencyIssue]:
        """Verify inter-module dependencies are satisfied"""
        issues = []
        module_types = {n.module_type for n in self.nodes.values()}
        
        for dep in INTERMODULE_DEPENDENCIES:
            if dep.relationship != 'requires':
                continue
            
            # Find modules that have this dependency
            source_mods = [
                n for n in self.nodes.values()
                if n.module_type == dep.source_module_type
            ]
            
            if not source_mods:
                continue
            
            # Check if target exists
            if dep.target_module_type not in module_types:
                for source in source_mods:
                    issues.append(DependencyIssue(
                        source_module_id=source.module_id,
                        source_module_type=source.module_type,
                        target_module_type=dep.target_module_type,
                        relationship='requires',
                        message=f"{source.module_type} requires {dep.target_module_type} but it was not found",
                        severity='error' if dep.relationship == 'requires' else 'warning'
                    ))
        
        return issues

    def check_preferences(self) -> List[PreferenceWarning]:
        """Check module preferences (shuttle + yeast, etc.)"""
        warnings = []
        module_types = {n.module_type for n in self.nodes.values()}
        
        for dep in INTERMODULE_DEPENDENCIES:
            if dep.relationship != 'prefers':
                continue
            
            # Find modules that have this preference
            source_mods = [
                n for n in self.nodes.values()
                if n.module_type == dep.source_module_type
            ]
            
            if not source_mods:
                continue
            
            # Check if preferred target exists
            if dep.target_module_type not in module_types:
                for source in source_mods:
                    warnings.append(PreferenceWarning(
                        source_module_id=source.module_id,
                        source_module_type=source.module_type,
                        preferred_module_type=dep.target_module_type,
                        message=f"{source.module_type} would benefit from {dep.target_module_type}",
                        suggestion=dep.description
                    ))
        
        return warnings

    def check_conflicts(self) -> List[DependencyIssue]:
        """Check for conflicting module combinations"""
        issues = []
        module_types = {n.module_type for n in self.nodes.values()}
        
        for dep in INTERMODULE_DEPENDENCIES:
            if dep.relationship != 'conflicts_with':
                continue
            
            # Check if both conflicting modules exist
            if dep.source_module_type in module_types and dep.target_module_type in module_types:
                source_mods = [
                    n for n in self.nodes.values()
                    if n.module_type == dep.source_module_type
                ]
                
                for source in source_mods:
                    issues.append(DependencyIssue(
                        source_module_id=source.module_id,
                        source_module_type=source.module_type,
                        target_module_type=dep.target_module_type,
                        relationship='conflicts_with',
                        message=f"{source.module_type} conflicts with {dep.target_module_type}",
                        severity='warning'
                    ))
        
        return issues

    def to_intent_analysis_format(self) -> Dict[str, Any]:
        """
        Export graph for Intent Analysis workflow.
        
        Returns structured representation suitable for GPT/LLM analysis.
        """
        # Organize modules by category
        expression_modules = []
        backbone_modules = []
        viral_modules = []
        other_modules = []
        
        for node in self.nodes.values():
            module_info = {
                'type': node.module_type,
                'start': node.start,
                'end': node.end,
                'length_bp': node.end - node.start,
                'feature_count': len(node.feature_ids),
            }
            
            if 'expression' in node.module_type or 'cassette' in node.module_type:
                expression_modules.append(module_info)
            elif 'backbone' in node.module_type or 'origin' in node.module_type:
                backbone_modules.append(module_info)
            elif 'lenti' in node.module_type or 'aav' in node.module_type:
                viral_modules.append(module_info)
            else:
                other_modules.append(module_info)
        
        # Build relationships summary
        relationships = []
        for edge in self.edges:
            if edge.relationship in ('requires', 'prefers', 'conflicts_with'):
                source_node = self.nodes.get(edge.source_id)
                target_node = self.nodes.get(edge.target_id)
                if source_node and target_node:
                    relationships.append({
                        'source_type': source_node.module_type,
                        'target_type': target_node.module_type,
                        'relationship': edge.relationship,
                    })
        
        # Build dependency status
        dep_issues = self.check_dependencies()
        pref_warnings = self.check_preferences()
        
        return {
            'plasmid_structure': {
                'sequence_length': self.seq_len,
                'circular': self.circular,
                'module_count': len(self.nodes),
                'expression_modules': expression_modules,
                'backbone_modules': backbone_modules,
                'viral_modules': viral_modules,
                'other_modules': other_modules,
            },
            'module_relationships': relationships,
            'dependency_status': {
                'satisfied': len(dep_issues) == 0,
                'issues': [i.to_dict() for i in dep_issues],
            },
            'preferences': {
                'warnings': [w.to_dict() for w in pref_warnings],
            },
            'summary': self._generate_summary(),
        }

    def _generate_summary(self) -> str:
        """Generate a human-readable summary of the plasmid structure"""
        lines = []
        
        # Count module types
        type_counts = {}
        for node in self.nodes.values():
            type_counts[node.module_type] = type_counts.get(node.module_type, 0) + 1
        
        if type_counts:
            lines.append("Module composition:")
            for mod_type, count in sorted(type_counts.items()):
                lines.append(f"  - {mod_type}: {count}")
        
        # Check for viral context
        if any('lenti' in t for t in type_counts.keys()):
            lines.append("\nVector type: Lentiviral")
        elif any('aav' in t for t in type_counts.keys()):
            lines.append("\nVector type: AAV")
        elif any('expression' in t for t in type_counts.keys()):
            lines.append("\nVector type: Expression plasmid")
        
        return "\n".join(lines)

    def to_dict(self) -> Dict[str, Any]:
        """Export full graph as dictionary"""
        return {
            'nodes': {k: v.to_dict() for k, v in self.nodes.items()},
            'edges': [e.to_dict() for e in self.edges],
            'sequence_length': self.seq_len,
            'circular': self.circular,
        }


# =============================================================================
# GENBANK ANNOTATION FUNCTIONS
# =============================================================================

def add_module_annotations_to_genbank(
    record: SeqRecord,
    modules: List[Dict[str, Any]]
) -> SeqRecord:
    """
    Add modules as feature annotations to final plasmid using custom 'module' type.
    
    Args:
        record: BioPython SeqRecord to modify
        modules: List of module dictionaries
        
    Returns:
        Modified SeqRecord with module features added
    """
    for module in modules:
        start = module.get('start', 0)
        end = module.get('end', 0)
        module_type = module.get('module_type', 'unknown')
        module_id = module.get('module_id', '')
        
        # Handle wrapping modules
        if module.get('wraps', False):
            # For wrapping modules, we need to create two features
            # or use a compound location
            seq_len = len(record.seq)
            
            # First segment: from start to end of sequence
            location1 = FeatureLocation(start, seq_len, strand=1)
            # Second segment: from beginning to end
            location2 = FeatureLocation(0, end, strand=1)
            
            # Create compound location
            from Bio.SeqFeature import CompoundLocation
            location = CompoundLocation([location1, location2])
        else:
            location = FeatureLocation(start, end)
        
        # Build grammar path from feature names
        grammar_path = module.get('grammar_path', [])
        if not grammar_path and 'features' in module:
            grammar_path = [f.get('feature_name', 'unknown') for f in module['features'][:5]]
        
        # Create feature
        feature = SeqFeature(
            location=location,
            type='module',  # Custom type for better identification
            qualifiers={
                'label': module_type.replace('_', ' ').title(),
                'module_type': module_type,
                'module_id': module_id,
                'components': ';'.join(grammar_path) if grammar_path else '',
                'note': f'Functional module: {module_type.replace("_", " ")}',
            }
        )
        
        # Add metadata as additional qualifiers
        metadata = module.get('metadata', {})
        if metadata:
            for key, value in metadata.items():
                if isinstance(value, (str, int, float)) and key not in feature.qualifiers:
                    feature.qualifiers[key] = str(value)
        
        record.features.append(feature)
    
    return record


def add_module_graph_to_genbank(
    record: SeqRecord,
    graph: ModuleGraph
) -> SeqRecord:
    """
    Add module graph information as a structured comment to GenBank record.
    
    This adds the graph as a comment annotation for tools that can parse it.
    """
    import json
    
    # Get intent analysis format
    intent_data = graph.to_intent_analysis_format()
    
    # Add as structured comment
    if 'structured_comment' not in record.annotations:
        record.annotations['structured_comment'] = {}
    
    record.annotations['structured_comment']['Module-Graph-Data'] = {
        'module_count': str(len(graph.nodes)),
        'relationship_count': str(len(graph.edges)),
        'dependency_issues': str(len(intent_data['dependency_status']['issues'])),
    }
    
    return record


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def build_module_graph(
    modules: List[Dict[str, Any]],
    features: List[Dict[str, Any]],
    sequence_length: int,
    circular: bool = True
) -> ModuleGraph:
    """
    Build a complete module graph from modules and features.
    
    This is the main entry point for building module graphs.
    """
    graph = ModuleGraph(modules, features, sequence_length, circular)
    graph.build_graph()
    return graph


def analyze_module_dependencies(
    modules: List[Dict[str, Any]],
    features: List[Dict[str, Any]],
    sequence_length: int,
    circular: bool = True
) -> Dict[str, Any]:
    """
    Analyze module dependencies and return comprehensive report.
    
    Returns dict with:
    - dependency_issues: List of failed required dependencies
    - preference_warnings: List of unmet preferences
    - conflict_warnings: List of conflicting module combinations
    - graph: The module graph (for further analysis)
    """
    graph = build_module_graph(modules, features, sequence_length, circular)
    
    return {
        'dependency_issues': [i.to_dict() for i in graph.check_dependencies()],
        'preference_warnings': [w.to_dict() for w in graph.check_preferences()],
        'conflict_warnings': [c.to_dict() for c in graph.check_conflicts()],
        'graph': graph.to_dict(),
        'intent_analysis': graph.to_intent_analysis_format(),
    }
