#!/usr/bin/env python3
"""
Plasmid tokenization module.

Converts hierarchical overlapping module annotations into a collapsed single-layer
representation suitable for AI training. Each position gets exactly one module path
and 0-1 features.

Usage:
    python tokenizer.py input.gb --output tokens.json --level functional
    python tokenizer.py input.gb --validate  # Check single-layer coverage
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Set

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class CollapsedSpan:
    """A non-overlapping span with single module path and optional feature."""
    start: int
    end: int
    module_path: str       # e.g., "lentiviral_payload:pol2_cassette(CAS9):cds_module"
    path_root: str         # e.g., "lentiviral_payload"
    path_leaf: str         # e.g., "cds_module"
    path_depth: int        # Number of levels in path
    payload_id: Optional[str] = None  # e.g., "CAS9"
    feature_canonical_id: Optional[str] = None
    feature_class: Optional[str] = None
    span_order: int = 0
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
    
    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass  
class TokenizationResult:
    """Result of tokenizing a plasmid."""
    plasmid_id: str
    sequence_length: int
    circular: bool
    spans: List[CollapsedSpan]
    tokens: List[str]
    validation: Dict[str, Any]
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "plasmid_id": self.plasmid_id,
            "sequence_length": self.sequence_length,
            "circular": self.circular,
            "spans": [s.to_dict() for s in self.spans],
            "tokens": self.tokens,
            "validation": self.validation,
        }


# =============================================================================
# MODULE TYPE HIERARCHY
# =============================================================================

# Define which module types are containers (can have children)
CONTAINER_TYPES = {
    "mammalian_lentiviral_expression",
    "lentiviral_payload",
    "aav_payload",
    "shuttle_vector_backbone",
    "pol2_expression_animal",
    "pol2_expression_generic",
    "pol2_expression_plant",
    "pol3_u6_sgrna_cassette",
    "pol3_expression_backbone",
    "cds_module",
    "bacterial_marker_cassette",
}

# Define module type priority (higher = more specific/preferred)
MODULE_PRIORITY = {
    # Sub-CDS level (most specific)
    "nls_module": 100,
    "linker_module": 100,
    "tag_module": 100,
    "tag": 100,
    "protein_module": 95,
    
    # CDS level
    "cds_module": 90,
    "cds_only": 85,
    
    # Selection cassettes
    "mammalian_selection_cassette": 80,
    "bacterial_marker_cassette": 80,
    
    # Expression cassettes
    "pol2_expression_animal": 70,
    "pol2_expression_generic": 70,
    "pol2_expression_plant": 70,
    "pol3_u6_sgrna_cassette": 70,
    "pol3_expression_backbone": 70,
    
    # Payload containers
    "lentiviral_payload": 60,
    "aav_payload": 60,
    
    # Vector-level
    "mammalian_lentiviral_expression": 50,
    "mammalian_replication_module": 50,
    "polya_only": 45,
    
    # Backbone elements
    "bacterial_origin": 40,
    "bacterial_backbone": 35,
    "shuttle_vector_backbone": 35,
    
    # Gaps (lowest priority)
    "inter_module_gap": 10,
    "gap": 10,
}

# Path name normalization
PATH_NAME_MAP = {
    "mammalian_lentiviral_expression": "lenti_vector",
    "lentiviral_payload": "lenti_payload",
    "pol2_expression_animal": "pol2_cassette",
    "pol2_expression_generic": "pol2_cassette",
    "pol2_expression_plant": "pol2_cassette",
    "pol3_u6_sgrna_cassette": "pol3_cassette",
    "pol3_expression_backbone": "pol3_cassette",
    "mammalian_selection_cassette": "selection",
    "bacterial_marker_cassette": "bac_marker",
    "bacterial_origin": "bac_origin",
    "mammalian_replication_module": "mam_ori",
    "polya_only": "polya",
    "cds_module": "cds",
    "nls_module": "nls",
    "linker_module": "linker",
    "tag_module": "tag",
    "tag": "tag",
    "protein_module": "protein",
}


# =============================================================================
# TOKEN VOCABULARY
# =============================================================================

def normalize_canonical_id(canonical_id: str) -> str:
    """Normalize canonical_id to token-friendly format."""
    if not canonical_id:
        return "UNKNOWN"
    
    # Remove common prefixes
    cid = canonical_id.upper()
    for prefix in ["CDS_", "PROMOTER_", "POLYA_", "ORI_", "MARKER_", "LENTI_", "MISC_"]:
        if cid.startswith(prefix):
            cid = cid[len(prefix):]
            break
    
    # Clean up
    cid = re.sub(r"[^A-Z0-9_]", "_", cid)
    cid = re.sub(r"_+", "_", cid).strip("_")
    
    return cid or "UNKNOWN"


def map_feature_to_token(
    canonical_id: Optional[str],
    feature_class: Optional[str],
    feature_name: Optional[str] = None,
) -> Optional[str]:
    """Map a feature to its token representation."""
    if not canonical_id and not feature_name:
        return None
    
    cid = (canonical_id or "").upper()
    fclass = (feature_class or "").lower()
    fname = (feature_name or "").lower()
    
    # Promoters
    if "promoter" in cid.lower() or fclass == "promoter" or "promoter" in fname:
        name = normalize_canonical_id(canonical_id) if canonical_id else fname.upper()
        return f"<PROMOTER:{name}>"
    
    # PolyA signals
    if "polya" in cid.lower() or "poly(a)" in fname or "polya" in fname:
        name = normalize_canonical_id(canonical_id) if canonical_id else "GENERIC"
        return f"<POLYA:{name}>"
    
    # Origins
    if cid.startswith("ORI_") or fclass == "origin" or "ori" in fname:
        name = normalize_canonical_id(canonical_id) if canonical_id else "GENERIC"
        return f"<ORI:{name}>"
    
    # Selection markers
    if cid.startswith("MARKER_") or "marker" in fclass:
        name = normalize_canonical_id(canonical_id) if canonical_id else "GENERIC"
        return f"<MARKER:{name}>"
    
    # NLS
    if "nls" in cid.lower() or "nls" in fname:
        name = normalize_canonical_id(canonical_id) if canonical_id else "GENERIC"
        return f"<NLS:{name}>"
    
    # 2A linkers
    if any(x in cid for x in ["P2A", "T2A", "E2A", "F2A"]) or "2a" in fname:
        name = normalize_canonical_id(canonical_id) if canonical_id else "GENERIC"
        return f"<LINKER:{name}>"
    
    # Tags
    if any(x in cid for x in ["FLAG", "HA", "MYC", "V5", "HIS", "STREP"]):
        name = normalize_canonical_id(canonical_id) if canonical_id else "GENERIC"
        return f"<TAG:{name}>"
    
    # Lentiviral elements
    if any(x in cid for x in ["LTR", "PSI", "RRE", "CPPT", "WPRE", "GAG", "ENV"]):
        name = normalize_canonical_id(canonical_id) if canonical_id else "GENERIC"
        return f"<LENTI:{name}>"
    
    # CDS/proteins (default for CDS_* prefixes)
    if cid.startswith("CDS_") or fclass in ("cds", "cds_payload", "nuclease"):
        name = normalize_canonical_id(canonical_id) if canonical_id else "GENERIC"
        return f"<CDS:{name}>"
    
    # Generic feature
    if canonical_id:
        return f"<FEATURE:{normalize_canonical_id(canonical_id)}>"
    
    return None


# =============================================================================
# MODULE COLLAPSING ALGORITHM
# =============================================================================

def _extract_payload_from_label(label: str) -> Optional[str]:
    """Extract payload ID from module label like 'cds module (CDS_CAS9)'."""
    match = re.search(r"\(([^)]+)\)", label)
    if match:
        return normalize_canonical_id(match.group(1))
    return None


def _normalize_module_type(module_type: str) -> str:
    """Normalize module type to standard format."""
    # Remove parenthetical suffixes
    base = re.sub(r"\s*\([^)]*\)", "", module_type).strip()
    # Convert to snake_case
    base = re.sub(r"[\s-]+", "_", base.lower())
    # Remove common suffixes
    base = re.sub(r"_cassette$", "", base) if "selection" not in base and "marker" not in base else base
    return base


def build_containment_tree(modules: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """
    Build a containment tree from modules.
    
    Returns dict mapping module_id -> {"module": module, "parent": parent_module or None, "children": []}
    """
    tree: Dict[str, Dict[str, Any]] = {}
    
    # Initialize tree nodes
    for mod in modules:
        mod_id = mod.get("module_id") or mod.get("id") or f"{mod['start']}_{mod['end']}"
        tree[mod_id] = {
            "module": mod,
            "parent": None,
            "children": [],
            "id": mod_id,
        }
    
    # Sort by span size (largest first) to find parents
    sorted_mods = sorted(modules, key=lambda m: -(m["end"] - m["start"]))
    
    # For each module, find smallest containing module as parent
    for mod in sorted_mods:
        mod_id = mod.get("module_id") or mod.get("id") or f"{mod['start']}_{mod['end']}"
        mod_start = mod["start"]
        mod_end = mod["end"]
        
        best_parent = None
        best_parent_size = float("inf")
        
        for candidate in sorted_mods:
            if candidate is mod:
                continue
            
            cand_start = candidate["start"]
            cand_end = candidate["end"]
            cand_size = cand_end - cand_start
            
            # Check if candidate contains this module
            if cand_start <= mod_start and cand_end >= mod_end and cand_size > (mod_end - mod_start):
                if cand_size < best_parent_size:
                    best_parent = candidate
                    best_parent_size = cand_size
        
        if best_parent:
            parent_id = best_parent.get("module_id") or best_parent.get("id") or f"{best_parent['start']}_{best_parent['end']}"
            tree[mod_id]["parent"] = tree[parent_id]
            tree[parent_id]["children"].append(tree[mod_id])
    
    return tree


def build_ancestry_path(
    module: Dict[str, Any],
    tree: Dict[str, Dict[str, Any]],
    max_depth: int = 4,
) -> str:
    """Build colon-separated path from root to this module."""
    mod_id = module.get("module_id") or module.get("id") or f"{module['start']}_{module['end']}"
    
    if mod_id not in tree:
        # Not in tree, use module type directly
        mod_type = _normalize_module_type(module.get("module_type", "unknown"))
        payload = _extract_payload_from_label(module.get("label", ""))
        if payload:
            return f"{mod_type}({payload})"
        return mod_type
    
    # Build path from root to leaf
    path_parts = []
    node = tree[mod_id]
    
    while node and len(path_parts) < max_depth:
        mod = node["module"]
        mod_type = _normalize_module_type(mod.get("module_type", "unknown"))
        
        # Get display name from PATH_NAME_MAP or use normalized type
        display_name = PATH_NAME_MAP.get(mod_type, mod_type)
        
        # Extract payload ID
        payload = mod.get("payload_id") or _extract_payload_from_label(mod.get("label", ""))
        
        if payload:
            path_parts.append(f"{display_name}({payload})"  )
        else:
            path_parts.append(display_name)
        
        node = node.get("parent")
    
    # Reverse to get root-first order
    path_parts.reverse()
    
    return ":".join(path_parts)


def get_module_priority(module: Dict[str, Any]) -> int:
    """Get priority score for a module (higher = more specific)."""
    mod_type = _normalize_module_type(module.get("module_type", ""))
    return MODULE_PRIORITY.get(mod_type, 50)


def collapse_modules_to_spans(
    modules: List[Dict[str, Any]],
    features: List[Dict[str, Any]],
    sequence_length: int,
) -> List[CollapsedSpan]:
    """
    Convert overlapping hierarchical modules to non-overlapping spans.
    
    Each position gets exactly one module path and 0-1 features.
    """
    if not modules:
        # No modules - entire sequence is a gap
        return [CollapsedSpan(
            start=0,
            end=sequence_length,
            module_path="gap:uncharacterized",
            path_root="gap",
            path_leaf="uncharacterized",
            path_depth=1,
            span_order=0,
        )]
    
    # Build containment tree
    tree = build_containment_tree(modules)
    
    # Find all boundary points
    boundaries: Set[int] = {0, sequence_length}
    for mod in modules:
        boundaries.add(mod["start"])
        boundaries.add(mod["end"])
    
    sorted_boundaries = sorted(boundaries)
    
    # For each segment, find the deepest (most specific) covering module
    spans: List[CollapsedSpan] = []
    
    for i in range(len(sorted_boundaries) - 1):
        seg_start = sorted_boundaries[i]
        seg_end = sorted_boundaries[i + 1]
        
        if seg_start >= seg_end:
            continue
        
        # Find all modules covering this segment
        covering = [
            m for m in modules
            if m["start"] <= seg_start and m["end"] >= seg_end
        ]
        
        if not covering:
            # Gap region
            # Determine gap context from adjacent modules
            prev_mod = None
            next_mod = None
            for m in modules:
                if m["end"] == seg_start:
                    prev_mod = m
                if m["start"] == seg_end:
                    next_mod = m
            
            gap_context = "inter_module"
            if prev_mod and next_mod:
                gap_context = "between_cassettes"
            elif seg_start == 0 or seg_end == sequence_length:
                gap_context = "terminal"
            
            spans.append(CollapsedSpan(
                start=seg_start,
                end=seg_end,
                module_path=f"gap:{gap_context}",
                path_root="gap",
                path_leaf=gap_context,
                path_depth=1,
                span_order=len(spans),
            ))
            continue
        
        # Find deepest (most specific) module - use priority then smallest span
        best_mod = max(covering, key=lambda m: (get_module_priority(m), -(m["end"] - m["start"])))
        
        # Build ancestry path
        path = build_ancestry_path(best_mod, tree)
        path_parts = path.split(":")
        
        # Extract payload from path if present
        payload = None
        if "(" in path:
            match = re.search(r"\(([^)]+)\)", path)
            if match:
                payload = match.group(1)
        
        # Find most specific feature in this segment
        seg_features = [
            f for f in features
            if f.get("start", 0) <= seg_start and f.get("end", 0) >= seg_end
        ]
        
        # Prefer smallest feature that covers the segment
        primary_feature = None
        if seg_features:
            primary_feature = min(seg_features, key=lambda f: f.get("end", 0) - f.get("start", 0))
        
        spans.append(CollapsedSpan(
            start=seg_start,
            end=seg_end,
            module_path=path,
            path_root=path_parts[0].split("(")[0] if path_parts else "unknown",
            path_leaf=path_parts[-1].split("(")[0] if path_parts else "unknown",
            path_depth=len(path_parts),
            payload_id=payload,
            feature_canonical_id=primary_feature.get("canonical_id") if primary_feature else None,
            feature_class=primary_feature.get("feature_class") or primary_feature.get("canonical_type") if primary_feature else None,
            span_order=len(spans),
        ))
    
    # Merge adjacent spans with identical module paths
    merged = _merge_adjacent_spans(spans)
    
    # Reassign span orders
    for i, span in enumerate(merged):
        span.span_order = i
    
    return merged


def _merge_adjacent_spans(spans: List[CollapsedSpan]) -> List[CollapsedSpan]:
    """Merge adjacent spans with identical module paths."""
    if not spans:
        return []
    
    merged = [spans[0]]
    
    for span in spans[1:]:
        prev = merged[-1]
        
        # Check if can merge (same path and adjacent)
        if (span.module_path == prev.module_path and 
            span.start == prev.end and
            span.feature_canonical_id == prev.feature_canonical_id):
            # Merge by extending previous span
            merged[-1] = CollapsedSpan(
                start=prev.start,
                end=span.end,
                module_path=prev.module_path,
                path_root=prev.path_root,
                path_leaf=prev.path_leaf,
                path_depth=prev.path_depth,
                payload_id=prev.payload_id,
                feature_canonical_id=prev.feature_canonical_id,
                feature_class=prev.feature_class,
                span_order=prev.span_order,
            )
        else:
            merged.append(span)
    
    return merged


# =============================================================================
# TOKENIZATION
# =============================================================================

def tokenize_spans(
    spans: List[CollapsedSpan],
    level: str = "functional",
) -> List[str]:
    """Convert collapsed spans to token sequence."""
    tokens = []
    
    for span in spans:
        # Module-level token
        tokens.append(f"<MODULE:{span.module_path}>")
        
        # Feature-level token (if not a gap and level is functional or full)
        if level in ("functional", "full") and span.feature_canonical_id:
            feat_token = map_feature_to_token(
                span.feature_canonical_id,
                span.feature_class,
            )
            if feat_token:
                tokens.append(feat_token)
    
    return tokens


# =============================================================================
# VALIDATION
# =============================================================================

def validate_collapsed_spans(
    spans: List[CollapsedSpan],
    sequence_length: int,
) -> Dict[str, Any]:
    """Validate that collapsed spans provide complete non-overlapping coverage."""
    issues = []
    
    # Check complete coverage
    covered = set()
    for span in spans:
        for pos in range(span.start, span.end):
            if pos in covered:
                issues.append(f"Position {pos} covered multiple times")
            covered.add(pos)
    
    # Check for gaps in coverage
    for pos in range(sequence_length):
        if pos not in covered:
            issues.append(f"Position {pos} not covered")
    
    # Check adjacency
    for i in range(len(spans) - 1):
        if spans[i].end != spans[i + 1].start:
            issues.append(f"Gap between spans {i} and {i+1}: {spans[i].end} != {spans[i+1].start}")
    
    # Check single module per span
    for span in spans:
        if ":" not in span.module_path and "gap" not in span.module_path:
            # Single level path is fine
            pass
    
    return {
        "valid": len(issues) == 0,
        "coverage_complete": len(covered) == sequence_length,
        "span_count": len(spans),
        "issues": issues[:10],  # Limit to first 10 issues
    }


# =============================================================================
# MAIN TOKENIZATION FUNCTION
# =============================================================================

def tokenize_plasmid(
    modules: List[Dict[str, Any]],
    features: List[Dict[str, Any]],
    sequence_length: int,
    plasmid_id: str = "unknown",
    circular: bool = True,
    level: str = "functional",
) -> TokenizationResult:
    """
    Main tokenization function.
    
    Args:
        modules: List of module dicts from annotation pipeline
        features: List of feature dicts from annotation pipeline
        sequence_length: Total sequence length
        plasmid_id: Identifier for the plasmid
        circular: Whether plasmid is circular
        level: Tokenization level ('module', 'functional', 'full')
    
    Returns:
        TokenizationResult with spans, tokens, and validation info
    """
    # Collapse modules to spans
    spans = collapse_modules_to_spans(modules, features, sequence_length)
    
    # Generate tokens
    tokens = tokenize_spans(spans, level)
    
    # Validate
    validation = validate_collapsed_spans(spans, sequence_length)
    
    return TokenizationResult(
        plasmid_id=plasmid_id,
        sequence_length=sequence_length,
        circular=circular,
        spans=spans,
        tokens=tokens,
        validation=validation,
    )


# =============================================================================
# CLI INTERFACE
# =============================================================================

def _parse_genbank_modules(gb_path: Path) -> Tuple[List[Dict], List[Dict], int, str]:
    """Parse modules and features from annotated GenBank file."""
    from Bio import SeqIO
    
    record = SeqIO.read(str(gb_path), "genbank")
    sequence_length = len(record.seq)
    plasmid_id = record.id or gb_path.stem
    
    modules = []
    features = []
    
    for feat in record.features:
        if feat.type == "source":
            continue
        
        label = feat.qualifiers.get("label", [""])[0]
        
        # Determine if this is a module annotation
        is_module = False
        module_type = None
        
        # Check for module patterns in label
        module_patterns = [
            "lentiviral", "expression", "cassette", "payload", "backbone",
            "origin", "marker", "selection", "cds module", "nls module",
            "linker module", "tag module", "pol2", "pol3", "polya",
            "inter-module gap", "gap", "replication module",
        ]
        
        label_lower = label.lower()
        for pattern in module_patterns:
            if pattern in label_lower:
                is_module = True
                module_type = label_lower
                break
        
        start = int(feat.location.start)
        end = int(feat.location.end)
        strand = 1 if feat.location.strand != -1 else -1
        
        if is_module:
            modules.append({
                "module_id": f"{start}_{end}",
                "module_type": module_type,
                "label": label,
                "start": start,
                "end": end,
                "strand": strand,
                "payload_id": _extract_payload_from_label(label),
            })
        else:
            # Regular feature
            features.append({
                "feature_id": f"{start}_{end}_{label}",
                "canonical_id": label.replace(" ", "_").upper() if label else None,
                "canonical_type": feat.type,
                "feature_class": feat.type,
                "label": label,
                "start": start,
                "end": end,
                "strand": strand,
            })
    
    return modules, features, sequence_length, plasmid_id


def main():
    parser = argparse.ArgumentParser(
        description="Plasmid tokenization tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python tokenizer.py input.gb
  python tokenizer.py input.gb --output tokens.json
  python tokenizer.py input.gb --validate
  python tokenizer.py input.gb --format table
        """,
    )
    parser.add_argument("input", help="GenBank file path")
    parser.add_argument("--output", "-o", help="Output JSON file")
    parser.add_argument("--level", choices=["module", "functional", "full"], default="functional",
                       help="Tokenization level")
    parser.add_argument("--format", choices=["json", "tokens", "table"], default="json",
                       help="Output format")
    parser.add_argument("--validate", action="store_true", help="Validate single-layer coverage")
    
    args = parser.parse_args()
    
    # Parse input file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: File not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    
    try:
        modules, features, seq_len, plasmid_id = _parse_genbank_modules(input_path)
    except Exception as e:
        print(f"Error parsing GenBank file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Parsed {len(modules)} modules and {len(features)} features from {input_path.name}", file=sys.stderr)
    
    # Tokenize
    result = tokenize_plasmid(
        modules=modules,
        features=features,
        sequence_length=seq_len,
        plasmid_id=plasmid_id,
        level=args.level,
    )
    
    # Validate if requested
    if args.validate:
        print(f"\nValidation Results:", file=sys.stderr)
        print(f"  Valid: {result.validation['valid']}", file=sys.stderr)
        print(f"  Coverage complete: {result.validation['coverage_complete']}", file=sys.stderr)
        print(f"  Span count: {result.validation['span_count']}", file=sys.stderr)
        if result.validation["issues"]:
            print(f"  Issues:", file=sys.stderr)
            for issue in result.validation["issues"]:
                print(f"    - {issue}", file=sys.stderr)
    
    # Output
    if args.format == "json":
        output = json.dumps(result.to_dict(), indent=2)
    elif args.format == "tokens":
        output = "\n".join(result.tokens)
    elif args.format == "table":
        lines = ["Start\tEnd\tLength\tModule Path\tFeature"]
        for span in result.spans:
            lines.append(f"{span.start}\t{span.end}\t{span.length}\t{span.module_path}\t{span.feature_canonical_id or '-'}")
        output = "\n".join(lines)
    
    if args.output:
        Path(args.output).write_text(output)
        print(f"Output written to {args.output}", file=sys.stderr)
    else:
        print(output)


if __name__ == "__main__":
    main()
