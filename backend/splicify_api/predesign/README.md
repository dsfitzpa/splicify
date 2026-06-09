# Unified Pre-Design Architecture

A standardized layer for plasmid cloning workflow design that sits between intent parsing and operator selection.

## Quick Start

```python
from splicify_api.predesign import (
    DesignRequest, PartSpecification, InputSource,
    PartResolver, TargetPlasmidBuilder, CloningRouter
)

# 1. Create design request
request = DesignRequest(
    session_id="sess_001",
    user_message="Design Gibson assembly for CMV + eGFP + bGH",
    parts=[
        PartSpecification(
            name="CMV",
            source=InputSource.DIRECT_SEQUENCE,
            sequence="ATGC..." * 100,
            role="promoter"
        ),
        # ... more parts
    ]
)

# 2. Resolve parts
resolver = PartResolver()
resolved = await resolver.resolve_all(request.parts)

# 3. Build target
builder = TargetPlasmidBuilder()
target = builder.build_from_parts(resolved, topology="circular")

# 4. Route to workflow
router = CloningRouter()
candidates = await router.route(target, resolved, objective="balanced")

# 5. Get recommendation
best = candidates[0]
print(f"Recommended: {best.method.value}")
```

## Components

### 1. `design_request.py` - Input Structures

**InputSource** - Part source types:
- `DIRECT_SEQUENCE` - User-provided DNA
- `UPLOADED_FILE` - GenBank/FASTA file
- `FEATURE_NAME` - Lookup by name (e.g., "CMV")
- `FEATURE_PATTERN` - Search by description
- `HOMOLOGY_DERIVED` - Detected from homology

**PartSpecification** - Part before resolution:
```python
PartSpecification(
    name="CMV promoter",
    source=InputSource.FEATURE_NAME,
    feature_name="CMV",
    role="promoter"
)
```

**DesignRequest** - Complete request:
- `session_id`: Session identifier
- `user_message`: Original message
- `parts`: List of PartSpecification
- `target`: Optional target specification
- `suggested_workflow`: User preference

### 2. `part_resolver.py` - Part Resolution

**ResolvedPart** - Part with sequence:
- `sequence`: DNA sequence (uppercase)
- `canonical_id`: Standard identifier
- `confidence`: Resolution confidence (0-1)
- `to_module_dict()`: Convert for operators

**PartResolver** - Multi-source resolution:
```python
resolver = PartResolver()
resolved = await resolver.resolve_all(parts, context)
```

Strategies:
- ✅ Direct sequences
- ✅ GenBank/FASTA files
- 🔄 Feature name lookup (Phase 2)
- 🔄 Pattern matching (Phase 2)
- 🔄 Homology detection (Phase 6)

### 3. `target_builder.py` - Target Assembly

**TargetPlasmid** - Complete target:
- `sequence`: Assembled sequence
- `parts`: Ordered resolved parts
- `restriction_sites`: RE site map
- `type_iis_sites`: Type IIS site map
- `get_part_junctions()`: Junction info

**TargetPlasmidBuilder** - Assembly:
```python
builder = TargetPlasmidBuilder()
target = builder.build_from_parts(
    parts=resolved,
    assembly_order="listed",
    topology="circular"
)
```

Features:
- Sequence concatenation
- RE site scanning (14 enzymes)
- Type IIS scanning (BsaI, BsmBI, BbsI, SapI)
- Homology detection
- GC content calculation

### 4. `cloning_router.py` - Workflow Selection

**WorkflowMethod** - Supported workflows:
- `GIBSON` - Gibson assembly
- `GOLDEN_GATE` - Golden Gate assembly
- `RESTRICTION` - Restriction cloning
- `SDM` - Site-directed mutagenesis
- `GATEWAY` - Gateway cloning

**WorkflowCandidate** - Evaluated option:
- `method`: Workflow type
- `compatible`: Compatibility flag
- `total_cost_usd`: Estimated cost
- `total_calendar_days`: Time estimate
- `overall_risk_score`: Risk (0-1)

**CloningRouter** - Evaluation:
```python
router = CloningRouter()
candidates = await router.route(
    target=target,
    parts=resolved,
    objective="balanced"  # or "cost", "time", "risk"
)
```

Compatibility constraints:
- **Gibson:** 2-6 fragments, min 100bp/fragment
- **Golden Gate:** 2-32 fragments
- **Restriction:** 1-3 fragments, unique RE sites
- **SDM:** Requires anchor plasmid

### 5. `knowledge_base.py` - Feature Lookup

**PlannotateKnowledgeBase** - pLannotate integration:
```python
from splicify_api.predesign import get_knowledge_base

kb = get_knowledge_base()
kb.load()

# Search by name
results = await kb.search_feature("CMV")

# Search by pattern
results = await kb.search_pattern("strong promoter")
```

Configuration:
- Set `PLANNOTATE_DATA_DIR` environment variable
- Or place files in `/tmp/plannotate_data/`

Required files:
- `feature_knowledge_base.json`
- `swissprot_knowledge_base.json`

## Validation

```bash
# Quick validation
python3 validate_phase1.py

# Comprehensive test
python3 validate_unified_predesign.py

# Interactive demo
python3 demo_predesign.py
```

## Integration with Existing Code

### Operators (No Changes Required)

All existing operators remain unchanged:
- `GibsonOperator`
- `GoldenGateOperator`
- `RestrictionOperator`
- `SDMOperator`

The router calls `operator.evaluate()` to get metrics.

### chat.py Integration (Phase 5)

Insert at line ~375:
```python
# NEW: Pre-design phase
design_request = await _build_design_request(...)
resolved_parts = await PartResolver().resolve_all(...)
target = TargetPlasmidBuilder().build_from_parts(...)
candidates = await CloningRouter().route(...)

# Delegate to existing operator
best = candidates[0]
if best.method == WorkflowMethod.GIBSON:
    return await _execute_gibson(target, resolved_parts, ...)
```

## Architecture Benefits

✅ **Consistent Input Handling** - All workflows use same structures
✅ **Flexible Part Sources** - Sequences, files, names, patterns
✅ **Explicit Targets** - Clear representation of desired outcome
✅ **Smart Routing** - Automatic workflow selection
✅ **Backward Compatible** - Existing code unchanged
✅ **Future-Proof** - Easy to add new workflows

## File Structure

```
predesign/
├── __init__.py              # Module exports
├── design_request.py        # DesignRequest, PartSpecification, etc.
├── part_resolver.py         # PartResolver, ResolvedPart
├── knowledge_base.py        # PlannotateKnowledgeBase
├── target_builder.py        # TargetPlasmid, TargetPlasmidBuilder
├── cloning_router.py        # CloningRouter, WorkflowCandidate
└── README.md                # This file
```

## Status

**Implemented (Phases 1-4):**
- ✅ Core data structures
- ✅ Direct sequence resolution
- ✅ File parsing (GenBank/FASTA)
- ✅ Target assembly and annotation
- ✅ Workflow compatibility checking
- ✅ Multi-objective ranking

**Pending (Phases 2, 5-6):**
- 🔄 Knowledge base integration (files not yet deployed)
- 🔄 Feature name lookup
- 🔄 Pattern matching
- 🔄 chat.py integration
- 🔄 Operator delegation (using placeholder metrics)
- 🔄 Homology-based assembly
- 🔄 NCBI fetching
- 🔄 Synthesis cost estimation

## Documentation

- **API Guide:** See `PREDESIGN_IMPLEMENTATION.md`
- **Examples:** Run `demo_predesign.py`
- **Tests:** See `tests/predesign/`

## Support

For questions or issues:
1. Check `PREDESIGN_IMPLEMENTATION.md` for detailed documentation
2. Run demo scripts to see examples
3. Review test files for usage patterns

---

**Version:** 1.0 (Phases 1-4)
**Date:** April 9, 2026
**Status:** Core components complete, ready for Phase 5 integration
