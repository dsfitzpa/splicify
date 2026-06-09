import re
import pathlib
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
import time
import tempfile
import os
import sys
from pathlib import Path
from . import _data
from io import StringIO
from fastapi import APIRouter, HTTPException, UploadFile, File, Form
import json

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from .hierarchical_annotator import annotate_hierarchy_from_plannotate, annotate_hierarchy_from_plannotate_v2
from .grammar_pol2_detector import detect_pol2_cassettes

# Create router

# Load knowledge base once at module level
_KB_CACHE = None
_KB_RECORDS_CACHE = None

def _resolve_kb_path():
    """Find feature_knowledge_base.json. Search order:
      1. SPLICIFY_KB_PATH env var (explicit override)
      2. backend/splicify_api/feature_db_data/feature_knowledge_base.json
         (the GenoLIB-seeded location built by feature_db/build_reference_db.py)
         (legacy pLannotate install — kept for portability; may not exist).
    Returns None when nothing is found so callers can degrade gracefully
    instead of raising at module load."""
    import os
    candidates = []
    env = os.environ.get("SPLICIFY_KB_PATH")
    if env:
        candidates.append(Path(env))
    candidates.append(
        _data.data_path("feature_db_data", "feature_knowledge_base.json")
    )
    candidates.append(
    )
    for c in candidates:
        if c.exists():
            return c
    return None


def _load_kb():
    """Load knowledge base once and cache it (sseqid-indexed). Returns
    {} when no KB file is available so _lookup_kb_feature degrades to
    "name lookup miss" instead of raising FileNotFoundError. The latter
    used to be swallowed by _plannotate_rows_to_seqviz_annotations's
    bare-except, silently dropping every non-swissprot row."""
    global _KB_CACHE
    if _KB_CACHE is not None:
        return _KB_CACHE
    import json
    kb_path = _resolve_kb_path()
    if kb_path is None:
        print("[KB] WARN: feature_knowledge_base.json not found in any candidate path; KB lookups will return None.")
        _KB_CACHE = {}
        return _KB_CACHE
    print(f"[KB] loading {kb_path}")
    with open(kb_path, "r") as f:
        kb_data = json.load(f)
    _KB_CACHE = {}
    if "records" in kb_data:
        for record in kb_data["records"]:
            sseqid = record.get("sseqid", "")
            if sseqid:
                _KB_CACHE[sseqid] = record
    print(f"[KB] indexed {len(_KB_CACHE)} records by sseqid")
    return _KB_CACHE


def _load_kb_records():
    """Return the raw KB record list (same payload as _load_kb but unindexed).

    chat.py iterates over the full record list when matching against
    feature_name / normalized_feature_name / sseqid / alternative_names,
    which the sseqid-indexed dict does not support cleanly.
    """
    global _KB_RECORDS_CACHE
    if _KB_RECORDS_CACHE is None:
        import json
        try:
            with open(kb_path, "r") as f:
                kb_data = json.load(f)
            _KB_RECORDS_CACHE = list(kb_data.get("records", []))
        except Exception:
            _KB_RECORDS_CACHE = []
    return _KB_RECORDS_CACHE

def _lookup_kb_feature(sseqid):
    """Look up a feature in the knowledge base by sseqid"""
    kb = _load_kb()
    return kb.get(sseqid, {})

_SWISSPROT_KB_CACHE = None

def _resolve_swissprot_kb_path():
    """Find the SwissProt KB. Prefer the short-form lookup table built
    by feature_db/sources/swissprot_curate.py (~15 MB, dict keyed by
    entry_name) over the full curated KB (~110 MB)."""
    import os
    candidates = []
    env = os.environ.get("SPLICIFY_SWISSPROT_KB_PATH")
    if env:
        candidates.append(Path(env))
    fdd = _data.data_path("feature_db_data")
    candidates.append(fdd / "swissprot_short_kb.json")
    candidates.append(fdd / "swissprot_curated_kb.json")
    # Legacy pLannotate fallback — kept for portability; may not exist.
    for c in candidates:
        if c.exists():
            return c
    return None


def _load_swissprot_kb():
    """Load SwissProt KB once and cache it. Returns {} when no file is
    available so _lookup_swissprot degrades to "miss" instead of
    raising. Handles three input shapes:
      (a) the short KB: flat dict keyed by entry_name with display fields
      (b) the curated KB: {"records": [{"sseqid": ..., "intrinsic_properties": {...}}, ...]}
      (c) the legacy pLannotate KB: {"records": [...]} flat fields"""
    global _SWISSPROT_KB_CACHE
    if _SWISSPROT_KB_CACHE is not None:
        return _SWISSPROT_KB_CACHE
    import json
    kb_path = _resolve_swissprot_kb_path()
    if kb_path is None:
        print("[KB] WARN: swissprot KB not found; swissprot lookups will return None.")
        _SWISSPROT_KB_CACHE = {}
        return _SWISSPROT_KB_CACHE
    print(f"[KB] loading swissprot from {kb_path}")
    with open(kb_path, "r") as f:
        kb_data = json.load(f)
    out = {}
    if isinstance(kb_data, dict) and "records" not in kb_data:
        # Short-form dict (entry_name -> fields) — most common
        for entry_name, rec in kb_data.items():
            out[entry_name] = rec
    elif isinstance(kb_data, dict) and "records" in kb_data:
        # records[] shape — flatten intrinsic_properties up to top level
        for record in kb_data["records"]:
            eid = record.get("sseqid") or record.get("normalized_feature_name")
            if not eid:
                continue
            props = record.get("intrinsic_properties", {}) or {}
            descs = ((record.get("source") or {}).get("descriptions") or [])
            full_desc = descs[0] if descs else ""
            protein_name = full_desc
            for marker in (" OS=", " OX=", " GN=", " PE="):
                if marker in protein_name:
                    protein_name = protein_name.split(marker, 1)[0].strip()
            out[eid] = {
                "entry_name": eid,
                "gene_name": record.get("gene_name") or props.get("GN", "") or "",
                "protein_name": record.get("protein_name") or protein_name,
                "organism": record.get("organism") or props.get("OS", "") or "",
                "taxonomy_id": record.get("taxonomy_id") or props.get("OX", "") or "",
                "protein_existence": record.get("protein_existence") or props.get("PE", "") or "",
            }
    _SWISSPROT_KB_CACHE = out
    print(f"[KB] indexed {len(out)} swissprot entries")
    return _SWISSPROT_KB_CACHE

def _lookup_swissprot(sseqid):
    """
    Look up a SwissProt entry by sseqid (entry_name like FRSA_ECO24).
    Returns dict with: protein_name, gene_name, organism, etc.
    """
    kb = _load_swissprot_kb()
    entry = kb.get(sseqid, {})
    if entry:
        return {
            'feature_name': entry.get('gene_name') or entry.get('protein_name') or sseqid,
            'protein_name': entry.get('protein_name', ''),
            'gene_name': entry.get('gene_name', ''),
            'organism': entry.get('organism', ''),
            'taxonomy_id': entry.get('taxonomy_id', ''),
            'protein_existence': entry.get('protein_existence', ''),
            'description': entry.get('description', ''),
            'entry_name': entry.get('entry_name', sseqid),
        }
    return {}



def _deduplicate_annotations(df):
    # DEDUP: Added logging
    print(f"DEDUP INPUT: {len(df)} features" if not df.empty else "DEDUP INPUT: empty")
    """
    Additional deduplication to remove:
    1. Exact duplicates (same position and label)
    2. Overlapping synonymous features (same position, different names for same gene)
    3. High-overlap features (>80% overlap) - keep the longer/higher-scoring one
    """
    if df.empty:
        return df
    
    df = df.copy()
    qlen = int(df["qlen"].iloc[0]) if "qlen" in df.columns else 10000
    
    def get_positions(row):
        qstart = int(row.get("qstart", 0))
        qend = int(row.get("qend", 0))
        if qstart <= qend:
            return set(range(qstart, qend + 1))
        else:
            return set(range(qstart, qlen)) | set(range(0, qend + 1))
    
    # Create position key (start, end, strand)
    df["_pos_key"] = df.apply(
        lambda r: (int(r.get("qstart", 0)), int(r.get("qend", 0)), int(r.get("sframe", 1) > 0)), 
        axis=1
    )
    
    # Create label key (lowercase, stripped)
    df["_label_key"] = df["Feature"].apply(lambda x: str(x or "").lower().strip())
    
    # Sort by score descending to prefer higher-scoring features
    df = df.sort_values(by=["percmatch", "score"], ascending=[False, False])
    
    # Remove exact duplicates (same position AND same label)
    df = df.drop_duplicates(subset=["_pos_key", "_label_key"], keep="first")
    
    # For features at the exact same position with different labels,
    # keep only the highest-scoring one (they are likely synonyms)
    seen_exact_positions = set()
    indices_after_exact = []
    
    for idx in df.index:
        pos_key = df.loc[idx, "_pos_key"]
        if pos_key not in seen_exact_positions:
            seen_exact_positions.add(pos_key)
            indices_after_exact.append(idx)
    
    df = df.loc[indices_after_exact]
    print(f"DEDUP AFTER POSITION: {len(df)} features")
    
    # Handle high-overlap features (>80% overlap) with same label - keep longest
    keep_indices = []
    skipped = set()
    
    for i, idx1 in enumerate(df.index):
        if idx1 in skipped:
            continue
        
        row1 = df.loc[idx1]
        pos1 = get_positions(row1)
        label1 = row1["_label_key"]
        
        should_keep = True
        for j, idx2 in enumerate(df.index):
            if i >= j or idx2 in skipped:
                continue
            
            row2 = df.loc[idx2]
            label2 = row2["_label_key"]
            
            if label1 != label2:
                continue
            
            pos2 = get_positions(row2)
            
            overlap = pos1 & pos2
            if not pos1 or not pos2:
                continue
            
            overlap_frac1 = len(overlap) / len(pos1)
            overlap_frac2 = len(overlap) / len(pos2)
            max_overlap = max(overlap_frac1, overlap_frac2)
            
            if max_overlap > 0.8:
                if len(pos2) > len(pos1):
                    should_keep = False
                    break
                else:
                    skipped.add(idx2)
        
        if should_keep:
            keep_indices.append(idx1)
    
    df = df.loc[keep_indices]
    df = df.drop(columns=["_pos_key", "_label_key"], errors="ignore")
    
    return df.reset_index(drop=True)

router = APIRouter(prefix="/plannotate", tags=["plannotate"])

# Try to import pLannotate
PLANNOTATE_AVAILABLE = False
PLANNOTATE_ERROR = None

# 2026-05-14: replaced "from plannotate.annotate import annotate" with the
# splicify-owned wrapper. The wrapper still uses pLannotate under the hood
# (so its BLAST/scoring code is reused) but searches the feature_db_data/
# tree built by scripts/build_feature_db_data.py. This lets us own the
# reference data without forking pLannotate.
try:
    from .feature_annotator import annotate
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    PLANNOTATE_AVAILABLE = True
    print("✅ splicify feature_annotator loaded (backed by feature_db_data/)")
except ImportError as e1:
    if plannotate_path.exists():
        sys.path.insert(0, str(plannotate_path))
        try:
            from .feature_annotator import annotate
            from Bio import SeqIO
            from Bio.SeqRecord import SeqRecord
            from Bio.SeqFeature import SeqFeature, FeatureLocation
            PLANNOTATE_AVAILABLE = True
            print(f"✅ splicify feature_annotator loaded (after adding {plannotate_path} to sys.path)")
        except ImportError as e2:
            PLANNOTATE_ERROR = f"Found pLannotate repo but import failed: {e2}"
            print(f"❌ {PLANNOTATE_ERROR}")
    else:
        PLANNOTATE_ERROR = "pLannotate not installed"
        print(f"❌ {PLANNOTATE_ERROR}")


# Request/Response models
class AnnotateRequest(BaseModel):
    gb_text: str
    session_id: Optional[str] = "unknown"
    options: Optional[Dict[str, Any]] = {
        "linear": False,
        "detailed": True,
        "batch_size": 100,
        "file_name": "output"
    }


class Annotation(BaseModel):
    type: str
    location: str
    name: Optional[str] = None
    qualifiers: Dict[str, str] = {}


class AnnotateResponse(BaseModel):
    ok: bool
    annotated_gb: Optional[str] = None
    annotation_count: Optional[int] = 0
    annotations: Optional[List[Annotation]] = []
    module_annotations: Optional[List[Dict]] = []
    modules: Optional[List[Dict]] = []
    error: Optional[str] = None
    details: Optional[str] = None


class HealthResponse(BaseModel):
    ok: bool
    plannotate_available: bool
    error: Optional[str] = None
    version: str = "1.0"


@router.get("/health", response_model=HealthResponse)
async def health_check():
    """Check if pLannotate is available"""
    return HealthResponse(
        ok=True,
        plannotate_available=PLANNOTATE_AVAILABLE,
        error=PLANNOTATE_ERROR if not PLANNOTATE_AVAILABLE else None
    )


@router.post("/annotate_genbank", response_model=AnnotateResponse)
async def annotate_genbank(request: AnnotateRequest):
    """
    Annotate a GenBank file using pLannotate
    
    pLannotate's annotate() returns a pandas DataFrame with columns:
    - qstart, qend: positions
    - sframe: strand (1 or -1)
    - Feature: feature name
    - Type: feature type (CDS, misc_feature, etc.)
    - Description: feature description
    """
    
    if not PLANNOTATE_AVAILABLE:
        raise HTTPException(
            status_code=500,
            detail={
                "error": "pLannotate not available",
                "details": PLANNOTATE_ERROR
            }
        )
    
    if not request.gb_text:
        raise HTTPException(
            status_code=400,
            detail="Missing gb_text parameter"
        )
    
    try:
        # Extract options
        options = request.options or {}
        linear = options.get('linear', False)
        detailed = options.get('detailed', True)
        
        # Parse the input GenBank to get the original SeqRecord
        
        # Fix common LOCUS line formatting issues
        gb_text_fixed = request.gb_text
        lines = gb_text_fixed.split('\n')

        # Fix LOCUS line
        if lines and lines[0].startswith('LOCUS'):
            parts = lines[0].split()
            if len(parts) >= 7:
                name = parts[1]
                length = parts[2]
                mol_type = parts[3] if len(parts) > 3 else 'DNA'
                topology = parts[4] if len(parts) > 4 else 'linear'
                division = parts[5] if len(parts) > 5 else 'SYN'
                date = parts[6] if len(parts) > 6 else ''
                lines[0] = f"LOCUS       {name:<16} {length:>7} bp    {mol_type:<6}  {topology:<8} {division} {date}"

        # Fix ORIGIN section - add line numbers if missing
        origin_idx = None
        for i, line in enumerate(lines):
            if line.startswith('ORIGIN'):
                origin_idx = i
                break

        if origin_idx is not None:
            # Collect all sequence lines after ORIGIN
            seq_lines = []
            i = origin_idx + 1
            while i < len(lines) and not lines[i].startswith('//'):
                line = lines[i].strip()
                if line and not line.startswith('ORIGIN'):
                    # Remove line numbers (they're at the start, digits + whitespace)
                    # Example: "        1 agaagtacag catcggcct" or "1 agaagtacag"
                    import re
                    # Remove leading digits and whitespace
                    cleaned = re.sub(r'^\s*\d+\s*', '', line)
                    # Remove ALL whitespace to get just sequence
                    seq_only = ''.join(cleaned.split())
                    # Only keep valid DNA characters
                    seq_only = re.sub(r'[^ACGTNacgtn]', '', seq_only)
                    if seq_only:
                        seq_lines.append(seq_only)
                i += 1
            
            # Rebuild with proper numbering
            full_seq = ''.join(seq_lines)
            new_origin_lines = ['ORIGIN      ']
            
            for pos in range(0, len(full_seq), 60):
                chunk = full_seq[pos:pos+60]
                # Add spaces every 10 bases
                spaced = ' '.join([chunk[i:i+10] for i in range(0, len(chunk), 10)])
                line_num = pos + 1
                new_origin_lines.append(f"{line_num:>9} {spaced}")
            
            # Replace origin section
            end_idx = i
            lines = lines[:origin_idx] + new_origin_lines + lines[end_idx:]

        gb_text_fixed = '\n'.join(lines)

        # Parse the input GenBank to get the original SeqRecord
        gb_io = StringIO(gb_text_fixed)
        original_record = SeqIO.read(gb_io, "genbank")
        
        print(f"🔬 Annotating GenBank (linear={linear}, detailed={detailed})...")

        # annotate() expects a raw DNA sequence string (not a file path)
        dna_sequence = str(original_record.seq)
        df_annotations = annotate(
            inSeq=dna_sequence,
            linear=linear,
            is_detailed=detailed
        )
        
        # Additional deduplication for exact duplicates and high-overlap synonyms
        if not df_annotations.empty:
            df_annotations = _deduplicate_annotations(df_annotations)
        
        # Additional deduplication for exact duplicates and high-overlap synonyms
        if not df_annotations.empty:
            df_annotations = _deduplicate_annotations(df_annotations)

        # Create new SeqRecord with pLannotate features added
        rec_annotations = original_record.annotations.copy()
        if 'molecule_type' not in rec_annotations:
            rec_annotations['molecule_type'] = 'DNA'
        if 'topology' not in rec_annotations:
            rec_annotations['topology'] = 'circular' if not linear else 'linear'

        new_record = SeqRecord(
            seq=original_record.seq,
            id=original_record.id,
            name=original_record.name,
            description=original_record.description,
            annotations=rec_annotations
        )

        # Keep original features for display in GenBank output
        new_record.features = list(original_record.features)

        # Add pLannotate annotations from DataFrame
        annotation_count = 0
        annotations = []

        if not df_annotations.empty:
            for idx, row in df_annotations.iterrows():
                try:
                    start = int(row.get('qstart', 0))
                    end = int(row.get('qend', 0))
                    strand = 1 if row.get('sframe', 1) > 0 else -1
                    feature_type = str(row.get('Type', 'misc_feature'))
                    sseqid = str(row.get('sseqid', ''))
                    blast_description = str(row.get('Description', ''))
                    
                    # Look up full KB entry for better naming and descriptions
                    db_name = str(row.get('db', ''))
                    kb_entry = None
                    swiss_entry = None
                    
                    if db_name == 'swissprot':
                        # Use SwissProt KB for swissprot hits
                        swiss_entry = _lookup_swissprot(sseqid)
                        if swiss_entry:
                            feature_name = swiss_entry.get('feature_name', sseqid)
                        else:
                            feature_name = csv_feature or sseqid
                    else:
                        kb_entry = _lookup_kb_feature(sseqid)
                        feature_name = kb_entry.get('feature_name', sseqid)  # Use KB name if available
                    
                    # Build comprehensive description
                    descriptions = []
                    if swiss_entry:
                        # Build description from SwissProt metadata
                        desc_parts = []
                        if swiss_entry.get('protein_name'):
                            desc_parts.append(swiss_entry['protein_name'])
                        if swiss_entry.get('organism'):
                            desc_parts.append('Organism: ' + swiss_entry['organism'])
                        if swiss_entry.get('gene_name'):
                            desc_parts.append('Gene: ' + swiss_entry['gene_name'])
                        if desc_parts:
                            descriptions.append('; '.join(desc_parts))
                    elif kb_entry and 'source' in kb_entry and 'descriptions' in kb_entry['source']:
                        descriptions.extend(kb_entry['source']['descriptions'])
                    elif blast_description:
                        descriptions.append(blast_description)
                    
                    description = '; '.join(descriptions) if descriptions else ''

                    location = FeatureLocation(start - 1, end, strand=strand)
                    feature = SeqFeature(
                        location=location,
                        type=feature_type,
                        qualifiers={
                            'label': feature_name,  # Now using KB feature_name
                            'note': description,    # Full KB descriptions
                            'description': description,  # Also as 'description' qualifier
                            'source': 'pLannotate',
                            'sseqid': sseqid       # Keep sseqid for reference
                        }
                    )

                    new_record.features.append(feature)
                    annotation_count += 1

                    annotations.append({
                        "type": feature_type,
                        "location": f"{start}..{end}",
                        "name": feature_name,  # Now using KB feature_name
                        "description": description,  # Top-level description for frontend display
                        "qualifiers": {
                            "label": feature_name,
                            "note": description,  # Full KB descriptions
                            "description": description,  # Also in qualifiers for GenBank compatibility
                            "sseqid": sseqid
                            # KB metadata available internally but not exposed in response
                        }
                    })
                except Exception as feature_err:
                    print(f"Warning: Could not add feature {idx}: {feature_err}")
                    continue

        print(f"✅ Added {annotation_count} pLannotate annotations")

        # Generate hierarchical annotations BEFORE writing GenBank
        module_annotations = []
        modules = []
        try:
            # Use pLannotate rows, or fall back to original file features
            if not df_annotations.empty:
                plannotate_rows = df_annotations.to_dict(orient="records")
            else:
                # Convert original GenBank features to pLannotate-style rows
                plannotate_rows = []
                for f in original_record.features:
                    if f.type == "source":
                        continue
                    name = f.qualifiers.get("label", f.qualifiers.get("gene", ["unknown"]))[0]
                    plannotate_rows.append({
                        "qstart": int(f.location.start) + 1,
                        "qend": int(f.location.end),
                        "sframe": 1 if f.location.strand == 1 else -1,
                        "Feature": name,
                        "Type": f.type,
                        "sseqid": name,
                        "Description": f.qualifiers.get("note", [""])[0] if "note" in f.qualifiers else "",
                    })
                print(f"Using {len(plannotate_rows)} original features (pLannotate returned 0)")
            
            hierarchy = annotate_hierarchy_from_plannotate_v2(
                sequence=dna_sequence,
                circular=not linear,
                plannotate_rows=plannotate_rows,
            )
            module_annotations = hierarchy.get("module_annotations", [])
            modules = hierarchy.get("modules", [])
            print(f"✅ Generated {len(module_annotations)} hierarchical module annotations")
            
            # Remove existing module annotations from input to avoid duplicates
            new_record.features = [
                f for f in new_record.features
                if not any(x in f.qualifiers.get('label', [''])[0].lower() 
                          for x in ['module', 'inter-module gap'])
            ]

            # Add module annotations as features to the GenBank record
            from Bio.SeqFeature import SeqFeature, FeatureLocation
            detected_feature_count = 0
            for mod_ann in module_annotations:
                mod_start = mod_ann.get("start", 0)
                mod_end = mod_ann.get("end", 0)
                mod_name = mod_ann.get("name", "unknown module")
                mod_type = mod_ann.get("module_type", "unknown")
                mod_strand = mod_ann.get("direction", 1)
                mod_metadata = mod_ann.get("metadata", {})
                
                # Create feature location (BioPython uses 0-based coordinates)
                location = FeatureLocation(mod_start, mod_end, strand=mod_strand)
                
                # Create misc_feature for module
                feature = SeqFeature(
                    location=location,
                    type="misc_feature",
                    qualifiers={
                        "label": mod_name,
                        "note": f"Functional module: {mod_type}",
                        "module_type": mod_type,
                    }
                )
                new_record.features.append(feature)
                
                # If module has a detected feature (tag/linker), also add it as a separate annotation
                detected_name = mod_metadata.get("detected_feature_name")
                canonical_id = mod_metadata.get("canonical_id")
                if detected_name and mod_type in ("tag_module", "flexible_linker_module"):
                    detected_feature = SeqFeature(
                        location=location,
                        type="misc_feature",
                        qualifiers={
                            "label": detected_name,
                            "note": f"Detected by sequence match",
                            "canonical_id": canonical_id or "",
                        }
                    )
                    new_record.features.append(detected_feature)
                    detected_feature_count += 1
            
            print(f"✅ Added {len(module_annotations)} module features to GenBank")
            if detected_feature_count > 0:
                print(f"✅ Added {detected_feature_count} detected tag/linker features")

            # Apply CDS boundary corrections and remove filtered features
            cds_filtered = hierarchy.get("cds_filtered_features", [])
            cds_corrections = hierarchy.get("cds_boundary_corrections", [])

            if cds_filtered or cds_corrections:
                print(f"[CDS Resolution] Filtered features: {len(cds_filtered)}, Boundary corrections: {len(cds_corrections)}")

                # Build set of filtered feature names/positions for removal
                filtered_set = set()
                for ff in cds_filtered:
                    # Match by name and approximate position
                    filtered_set.add((ff["name"], ff["start"], ff["end"]))

                # Build dict of boundary corrections: (name, orig_start, orig_end) -> (new_start, new_end)
                correction_map = {}
                for bc in cds_corrections:
                    key = (bc["name"], bc["original_start"], bc["original_end"])
                    correction_map[key] = (bc["corrected_start"], bc["corrected_end"])

                # Filter and correct features in new_record
                filtered_features = []
                for feat in new_record.features:
                    if feat.type == "source":
                        filtered_features.append(feat)
                        continue

                    label = feat.qualifiers.get("label", [""])[0]
                    # BioPython uses 0-based start
                    feat_start = int(feat.location.start)
                    feat_end = int(feat.location.end)

                    # Check if this feature should be removed (filtered)
                    should_remove = False
                    for ff_name, ff_start, ff_end in filtered_set:
                        # Match by name and position (allowing 3bp tolerance)
                        if label == ff_name:
                            if abs(feat_start - ff_start) <= 3 and abs(feat_end - ff_end) <= 3:
                                should_remove = True
                                print(f"  [REMOVE] {label} at {feat_start+1}..{feat_end} (filtered as >90% covered)")
                                break

                    if should_remove:
                        continue

                    # Check if this feature needs boundary correction
                    for (bc_name, bc_orig_start, bc_orig_end), (bc_new_start, bc_new_end) in correction_map.items():
                        if label == bc_name:
                            # Match by approximate position
                            if abs(feat_start - bc_orig_start) <= 3 and abs(feat_end - bc_orig_end) <= 3:
                                # Apply correction
                                from Bio.SeqFeature import FeatureLocation
                                new_location = FeatureLocation(
                                    bc_new_start,  # Already 0-based
                                    bc_new_end,
                                    strand=feat.location.strand
                                )
                                feat.location = new_location
                                print(f"  [CORRECT] {label}: {feat_start+1}..{feat_end} -> {bc_new_start+1}..{bc_new_end}")
                                break

                    filtered_features.append(feat)

                # Also filter out mammalian_selection_cassette features that overlap with protein submodules
                # Get protein submodule ranges from module_annotations
                protein_submodule_ranges = [
                    (ann.get('start', 0), ann.get('end', 0))
                    for ann in module_annotations
                    if ann.get('module_type') == 'protein_module'
                    and ann.get('metadata', {}).get('module_family') == 'cds_submodule'
                ]

                def overlaps_protein_submodule(feat_start, feat_end):
                    feat_len = feat_end - feat_start
                    if feat_len <= 0:
                        return False
                    for ps_start, ps_end in protein_submodule_ranges:
                        overlap_start = max(feat_start, ps_start)
                        overlap_end = min(feat_end, ps_end)
                        overlap = max(0, overlap_end - overlap_start)
                        if overlap / feat_len > 0.8:
                            return True
                    return False

                final_features = []
                for feat in filtered_features:
                    label = feat.qualifiers.get('label', [''])[0].lower()
                    if 'selection cassette' in label or 'marker_puro' in label or 'marker_neo' in label:
                        feat_start = int(feat.location.start)
                        feat_end = int(feat.location.end)
                        if overlaps_protein_submodule(feat_start, feat_end):
                            print(f"  [REMOVE] {feat.qualifiers.get('label', [''])[0]} (selection cassette overlaps protein submodule)")
                            continue
                    final_features.append(feat)

                new_record.features = final_features
                print(f"[CDS Resolution] Final feature count: {len(new_record.features)}")
            
        except Exception as hier_err:
            print(f"⚠️  Hierarchical annotation failed: {hier_err}")
            import traceback
            traceback.print_exc()

        # Convert SeqRecord back to GenBank text (now includes module features)
        output = StringIO()
        SeqIO.write([new_record], output, "genbank")
        annotated_gb = output.getvalue()

        return AnnotateResponse(
            ok=True,
            annotated_gb=annotated_gb,
            annotation_count=annotation_count,
            annotations=annotations[:20],
            module_annotations=module_annotations,  # Add hierarchical annotations
            modules=modules
        )

    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        print(f"❌ Annotation error: {error_details}")
        
        raise HTTPException(
            status_code=500,
            detail={
                "error": str(e),
                "traceback": error_details
            }
        )
@router.post("/annotate_genbank_file", response_model=AnnotateResponse)
async def annotate_genbank_file(
    file: UploadFile = File(...),
    session_id: str = Form("unknown"),
    options_json: str = Form('{"linear": false, "batch_size": 100, "file_name": "output"}'),
):
    """
    Accepts a .gb/.gbk upload as multipart/form-data (n8n binary upload),
    converts bytes -> text, and annotates with pLannotate.
    """
    if not PLANNOTATE_AVAILABLE:
        raise HTTPException(status_code=500, detail={"error": "pLannotate not available", "details": PLANNOTATE_ERROR})

    # Decode uploaded file bytes -> text
    raw = await file.read()
    try:
        gb_text = raw.decode("utf-8")
    except UnicodeDecodeError:
        # GenBank is almost always UTF-8/ASCII; fall back safely
        gb_text = raw.decode("latin-1", errors="replace")

    # Parse options
    try:
        options = json.loads(options_json) if options_json else {}
    except Exception:
        options = {}

    # Reuse the same logic as annotate_genbank by calling it directly
    req = AnnotateRequest(gb_text=gb_text, session_id=session_id, options=options)
    return await annotate_genbank(req)


# ---------------------------------------------------------------------------
# Annotate a raw DNA sequence — returns SeqViz-compatible annotation list
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Analyze plasmid intent — LLM-powered biological purpose summary
# ---------------------------------------------------------------------------

class AnnotationItem(BaseModel):
    name: str
    start: int
    end: int
    role: Optional[str] = None
    color: Optional[str] = None


class AnalyzeIntentRequest(BaseModel):
    annotations: List[AnnotationItem]
    title: str = "Plasmid"
    sequence_length: int = 0
    circular: bool = True

class AnalyzePlasmidRequest(BaseModel):
    """Request model for comprehensive plasmid analysis"""
    annotations: List[AnnotationItem] = []  # Optional - will be generated from sequence if not provided
    sequence: str = ""
    title: str = "Plasmid"
    sequence_length: int = 0
    circular: bool = True


def _format_plasmid_analysis(result: Dict[str, Any], title: str, gpt_description: str = "") -> str:
    """Format the plasmid analysis result into a human-readable report"""
    lines = []

    # Header
    lines.append(f"# Plasmid Analysis: {title}")
    lines.append("")

    # GPT-Generated Overview (if available)
    if gpt_description:
        lines.append("## Overview")
        lines.append(gpt_description)
        lines.append("")

    # Intent/Purpose
    profile = result.get("construct_profile", {})
    primary_intent = profile.get("primary_intent")
    if primary_intent:
        intent_map = {
            "mammalian_transient_expression": "Mammalian Transient Expression Vector",
            "lentiviral_delivery": "Lentiviral Transfer/Delivery Vector",
            "aav_delivery": "AAV Gene Delivery Vector",
            "crispr_genome_editing": "CRISPR Genome Editing Construct",
            "bacterial_expression": "Bacterial Expression Vector",
        }
        lines.append(f"## Intended Use")
        lines.append(f"**{intent_map.get(primary_intent, primary_intent.replace('_', ' ').title())}**")
        lines.append("")

    # Get features early (needed for module components)
    features = result.get("feature_instances", [])
    
    # Functional Modules
    modules = result.get("module_instances", [])
    if modules:
        lines.append("## Functional Modules")
        for mod in modules:
            mod_type = mod.get("module_type", "unknown")
            mod_name = mod_type.replace("_", " ").title()
            grammar = mod.get("grammar_path", [])

            lines.append(f"### {mod_name}")
            
            # Show actual feature names as components
            mod_feature_ids = mod.get("feature_ids", [])
            mod_features = [f for f in features if f.get("instance_id") in mod_feature_ids]
            if mod_features:
                feature_names = [f.get("feature_name", "unknown") for f in mod_features[:5]]  # Limit to first 5
                if len(mod_features) > 5:
                    feature_names.append(f"...({len(mod_features) - 5} more)")
                lines.append(f"Components: {', '.join(feature_names)}")

            # Add module-specific details
            if "expression_cassette" in mod_type:
                lines.append("Function: Drives expression of protein payload")
            elif "guide_cassette" in mod_type:
                lines.append("Function: Expresses sgRNA for CRISPR targeting")
            elif "nuclease_cassette" in mod_type:
                lines.append("Function: Expresses Cas9/nuclease for genome editing")
            elif "backbone" in mod_type:
                lines.append("Function: Enables plasmid replication and selection")
            elif "lentiviral" in mod_type:
                lines.append("Function: Lentiviral packaging and transfer region")
            elif "aav" in mod_type:
                lines.append("Function: AAV payload packaging region")

            lines.append("")

    # Replication System
    origins = [f for f in features if f.get("role") == "origin"]
    markers = [f for f in features if f.get("role") == "bacterial_marker"]

    if origins or markers:
        lines.append("## Replication & Selection")
        if origins:
            ori_names = [f.get("feature_name", "unknown") for f in origins]
            lines.append(f"**Origin of replication:** {', '.join(ori_names)}")
        if markers:
            marker_names = [f.get("feature_name", "unknown") for f in markers]
            lines.append(f"**Selection marker:** {', '.join(marker_names)}")
        lines.append("")

    # Payload System
    # Only show actual protein-coding payloads, not RNA elements
    protein_coding_roles = ["expression_payload", "editing_payload", "reporter_payload", "selection_payload"]
    payloads = [f for f in features if f.get("role") in protein_coding_roles]
    if payloads:
        lines.append("## Payload")
        for p in payloads:
            name = p.get("feature_name", "unknown")
            role = p.get("role", "").replace("_", " ").title()
            has_start = p.get("has_start_codon")
            kozak = p.get("kozak_strength")

            lines.append(f"**{name}** ({role})")
            if has_start is not None:
                lines.append(f"  - Start codon: {'Present' if has_start else 'Missing ⚠️'}")
            if kozak:
                lines.append(f"  - Kozak context: {kozak.title()}")
        lines.append("")

    # Host System
    promoters = [f for f in features if "promoter" in f.get("role", "")]
    if promoters:
        lines.append("## Expression System")
        pol2 = [f for f in promoters if f.get("role") == "pol2_promoter"]
        pol3 = [f for f in promoters if f.get("role") == "pol3_promoter"]

        if pol2:
            names = [f.get("feature_name", "") for f in pol2]
            lines.append(f"**Pol II promoters (protein expression):** {', '.join(names)}")
        if pol3:
            names = [f.get("feature_name", "") for f in pol3]
            lines.append(f"**Pol III promoters (small RNA expression):** {', '.join(names)}")
        lines.append("")

    # Errors and Warnings
    lint_issues = result.get("lint_issues", [])
    errors = [i for i in lint_issues if i.get("severity") == "error"]
    warnings = [i for i in lint_issues if i.get("severity") == "warning"]

    if errors or warnings:
        lines.append("## Design Issues")

        if errors:
            lines.append("### ❌ Errors")
            for e in errors:
                lines.append(f"- **{e.get('rule_id', 'unknown')}**: {e.get('message', '')}")
                if e.get("suggestion"):
                    lines.append(f"  - *Suggestion: {e.get('suggestion')}*")
            lines.append("")

        if warnings:
            lines.append("### ⚠️ Warnings")
            for w in warnings:
                lines.append(f"- **{w.get('rule_id', 'unknown')}**: {w.get('message', '')}")
            lines.append("")
    else:
        lines.append("## Design Validation")
        lines.append("✅ No design issues detected. Plasmid appears to be correctly constructed.")
        lines.append("")

    return "\n".join(lines)


@router.post("/analyze_plasmid")
async def analyze_plasmid_endpoint(request: AnalyzePlasmidRequest):
    """
    Comprehensive plasmid analysis using the enhanced KB-backed analyzer.

    Returns a structured analysis including:
    - Intended use/purpose
    - Functional modules detected
    - Replication and selection systems
    - Payload information
    - Expression system details
    - Design errors and warnings
    """
    from .plasmid_analyzer import analyze_plasmid_from_plannotate

    if not request.sequence:
        return {"ok": False, "error": "Sequence required for analysis"}

    # Run pLannotate directly on sequence to get fresh annotations
    # This ensures we only use pLannotate features, not original GenBank features
    try:
        df_annotations = annotate(
            inSeq=request.sequence,
            linear=not request.circular,
            is_detailed=True
        )
        plannotate_rows = df_annotations.to_dict(orient="records") if not df_annotations.empty else []
    except Exception as e:
        import traceback
        return {"ok": False, "error": f"pLannotate failed: {str(e)}", "traceback": traceback.format_exc()}

    try:
        # Run the analyzer
        result = analyze_plasmid_from_plannotate(
            sequence=request.sequence,
            circular=request.circular,
            plannotate_rows=plannotate_rows,
            plasmid_id=request.title.replace(" ", "_")[:20]
        )

        # Get GPT description for the analysis
        gpt_description = ""
        try:
            import os
            api_key = os.getenv("OPENAI_API_KEY")
            if api_key and api_key != "your-openai-api-key-here":
                from openai import OpenAI
                client = OpenAI(api_key=api_key)

                # Build feature summary for GPT
                features = result.get("feature_instances", [])
                feature_summary = []
                for f in features[:30]:  # Limit to 30 features
                    name = f.get("feature_name", f.get("name", "unknown"))
                    role = f.get("role", "unknown")
                    start = f.get("start", 0)
                    end = f.get("end", 0)
                    feature_summary.append(f"- {name} ({role}): {start}-{end}")

                prompt = f"""Analyze this plasmid and write a concise 2-3 paragraph description of its purpose and function.

Plasmid: {request.title}
Sequence length: {len(request.sequence)} bp
Topology: {"circular" if request.circular else "linear"}

Features:
{chr(10).join(feature_summary)}

Write a clear, informative description for a molecular biologist. Focus on:
1. What type of vector this is and its primary purpose
2. Key expression cassettes and what they produce
3. Notable functional elements (selection markers, viral components, etc.)"""

                response = client.chat.completions.create(
                    model="gpt-4o-mini",
                    messages=[{"role": "user", "content": prompt}],
                    max_tokens=600,
                    temperature=0.3,
                )
                gpt_description = response.choices[0].message.content.strip()
        except Exception as e:
            print(f"[DEBUG] GPT analysis in analyze_plasmid failed: {e}")
            gpt_description = ""

        # Format the analysis
        analysis_text = _format_plasmid_analysis(result, request.title, gpt_description)

        # Build module graph data for download
        module_graph_data = {
            "plasmid_title": request.title,
            "sequence_length": len(request.sequence),
            "circular": request.circular,
            "modules": result.get("module_instances", []),
            "features": result.get("feature_instances", []),
            "construct_graph": result.get("construct_graph", {}),
            "lint_issues": result.get("lint_issues", []),
        }

        return {
            "ok": True,
            "analysis": analysis_text,
            "modules": [m.get("module_type") for m in result.get("module_instances", [])],
            "errors": len([i for i in result.get("lint_issues", []) if i.get("severity") == "error"]),
            "warnings": len([i for i in result.get("lint_issues", []) if i.get("severity") == "warning"]),
            "module_graph": module_graph_data,  # Include module graph for download
            "raw_result": result  # Include raw data for debugging
        }
    except Exception as e:
        import traceback
        return {"ok": False, "error": str(e), "traceback": traceback.format_exc()}




@router.post("/analyze_intent")
async def analyze_intent_endpoint(request: AnalyzeIntentRequest):
    """
    Generate an LLM-powered natural-language analysis of what a plasmid is
    designed to do, based on its annotation map.

    Input:  {annotations, title, sequence_length, circular}
    Output: {ok, analysis: string}
    """
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        return {"ok": False, "error": "OPENAI_API_KEY not configured"}

    if not request.annotations:
        return {"ok": False, "error": "No annotations provided"}

    # Sort annotations by start position for a coherent map description
    sorted_anns = sorted(request.annotations, key=lambda a: a.start)

    topology = "circular" if request.circular else "linear"
    length_str = f"{request.sequence_length:,} bp" if request.sequence_length else "unknown length"

    feature_lines = "\n".join(
        f"  - {a.name} ({a.start}..{a.end}{', role: ' + a.role if a.role else ''})"
        for a in sorted_anns
    )

    prompt = (
        f"You are a molecular biology expert. Analyze the following plasmid map and "
        f"write a clear, concise explanation of its biological intent.\n\n"
        f"Plasmid: {request.title}\n"
        f"Topology: {topology}\n"
        f"Total length: {length_str}\n\n"
        f"Features (in order):\n{feature_lines}\n\n"
        f"Write 2–3 paragraphs covering:\n"
        f"1. What type of vector this is and its primary purpose\n"
        f"2. How the expression cassette(s) work — what proteins or RNAs are produced "
        f"and what drives their expression\n"
        f"3. Any notable functional features (selection, viral elements, cloning sites, "
        f"delivery mechanism)\n\n"
        f"Be specific about the elements present. Write in plain English for a molecular "
        f"biologist. Do not repeat the feature list verbatim — synthesise it into a "
        f"coherent narrative."
    )

    try:
        from openai import AsyncOpenAI
        client = AsyncOpenAI(api_key=api_key)
        response = await client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[{"role": "user", "content": prompt}],
            max_tokens=600,
            temperature=0.3,
        )
        analysis = response.choices[0].message.content.strip()
        return {"ok": True, "analysis": analysis}
    except Exception as e:
        return {"ok": False, "error": str(e)}


class AnnotateSequenceRequest(BaseModel):
    sequence: str
    circular: bool = True
    detailed: bool = True
    hierarchical: bool = False



def _filter_hierarchical_dups(plannotate_annotations, hierarchical):
    """Drop hierarchical-layer features that exactly duplicate a Step 1
    plannotate annotation (same lowercase-name + same [start,end]).
    Avoids the same feature appearing twice in the frontend export when
    Step 4 echoes a Step 1 hit."""
    seen = set()
    for a in plannotate_annotations or []:
        seen.add((str(a.get("name","")).lower().strip(),
                  int(a.get("start",0)), int(a.get("end",0))))
    out = []
    for h in hierarchical or []:
        key = (str(h.get("name","")).lower().strip(),
               int(h.get("start",0)), int(h.get("end",0)))
        if key in seen:
            continue
        out.append(h)
    return out

def _dedup_cross_tier_annotations(annotations):
    """Two passes:
      1. Within each annotation, drop strand-duplicates: same display name +
         same [start,end] but opposite strand. (Symmetric palindromic features
         like ori / f1 ori / SV40 ori sometimes match both strands when blastn
         is run with --word_size 11; the user sees them twice.)
      2. Cross-tier dedup: drop a swissprot row if it overlaps >=70 % with
         any higher-priority hit from feature_protein / feature_reference /
         fpbase. Preserves the canonical GenoLIB name on shared CDSes.
    Lower-priority dbs:
      swissprot > fpbase > feature_motifs > feature_protein > feature_reference"""
    if not annotations:
        return annotations

    # ---- Pass 1: strand-duplicate dedup ---------------------------
    seen = {}
    kept = []
    for a in annotations:
        key = (a.get("name", "").lower().strip(),
               int(a.get("start", 0)), int(a.get("end", 0)))
        if key in seen:
            prev = seen[key]
            # Prefer + strand when in doubt
            if a.get("direction", 1) > prev.get("direction", 1):
                kept[kept.index(prev)] = a
                seen[key] = a
            continue
        seen[key] = a
        kept.append(a)

    # ---- Pass 2: cross-tier dedup (swissprot only) ----------------
    # Track covered intervals per "tier of authority"
    auth_intervals = []  # list of (start, end) from higher-priority dbs
    for a in kept:
        db = a.get("db", "")
        if db in ("feature_protein", "feature_reference", "feature_motifs", "fpbase"):
            auth_intervals.append((int(a.get("start", 0)), int(a.get("end", 0))))

    def _overlap_frac(s1, e1, s2, e2):
        if e1 <= s2 or e2 <= s1:
            return 0.0
        inter = max(0, min(e1, e2) - max(s1, s2))
        smaller = min(e1 - s1, e2 - s2) or 1
        return inter / smaller

    final = []
    for a in kept:
        if a.get("db") != "swissprot":
            final.append(a)
            continue
        s, e = int(a.get("start", 0)), int(a.get("end", 0))
        dominated = any(_overlap_frac(s, e, ai, ae) >= 0.70 for (ai, ae) in auth_intervals)
        if dominated:
            continue
        final.append(a)

    return final


def _plannotate_rows_to_seqviz_annotations(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    annotations: List[Dict[str, Any]] = []
    for row in rows:
        try:
            qstart = int(row.get("qstart", 1))
            qend = int(row.get("qend", 0))
            raw_frame = row.get("sframe", 1)
            sframe = int(raw_frame) if raw_frame is not None and str(raw_frame) not in ("", "nan") else 1
            direction = 1 if sframe > 0 else -1

            # KB lookup for better naming and descriptions.
            # The Feature column comes from _get_details CSV merge — it can
            # be the description text (e.g. swissprot.csv writes the full
            # protein description there). The real lookup key is sseqid.
            sseqid = str(row.get("sseqid", "") or row.get("Feature", "unknown"))
            csv_feature = str(row.get("Feature", "") or sseqid)
            description = str(row.get("Description", "") or "")
            db_name = str(row.get("db", ""))
            
            kb_entry = None
            swiss_entry = None
            kb_data = None  # Full KB metadata for gene card
            
            if db_name == 'swissprot':
                swiss_entry = _lookup_swissprot(sseqid)
                if swiss_entry:
                    # Display name preference: gene_name (short, user-readable)
                    # -> protein_name -> entry_name. Examples:
                    #   "bla"     instead of "AMPR_ECOLI"
                    #   "EGFP"    instead of "EGFP_AEQVI"
                    gene = (swiss_entry.get('gene_name') or '').strip()
                    pname = (swiss_entry.get('protein_name') or '').strip()
                    organism = (swiss_entry.get('organism') or '').strip()
                    if gene:
                        feature_name = gene
                    elif pname:
                        feature_name = pname[:60]
                    else:
                        feature_name = csv_feature or sseqid
                    desc_parts = []
                    if pname:
                        desc_parts.append(pname)
                    if organism:
                        desc_parts.append(organism)
                    desc_parts.append(f"UniProt: {sseqid}")
                    description = ' — '.join(desc_parts)
                    kb_data = {
                        "source_type": "swissprot",
                        "protein_name": pname,
                        "gene_name": gene,
                        "organism": organism,
                        "taxonomy_id": swiss_entry.get('taxonomy_id'),
                        "protein_existence": swiss_entry.get('protein_existence'),
                        "entry_name": swiss_entry.get('entry_name', sseqid),
                    }
                else:
                    feature_name = csv_feature or sseqid
            else:
                # Use regular KB for other databases
                kb_entry = _lookup_kb_feature(sseqid)
                kb_feature_name = kb_entry.get('feature_name', '') if kb_entry else ''
                if kb_feature_name and len(kb_feature_name) > 2:
                    feature_name = kb_feature_name
                else:
                    feature_name = csv_feature or sseqid
                    
                # Get KB descriptions if available
                if kb_entry and 'source' in kb_entry and 'descriptions' in kb_entry['source']:
                    kb_descriptions = kb_entry['source']['descriptions']
                    if kb_descriptions:
                        description = '; '.join(kb_descriptions)
                
                # Build kb_data for gene card (Feature KB format)
                if kb_entry:
                    props = kb_entry.get('intrinsic_properties', {})
                    source = kb_entry.get('source', {})
                    kb_data = {
                        "source_type": "feature_kb",
                        "feature_id": kb_entry.get('feature_id'),
                        "feature_name": kb_entry.get('feature_name'),
                        "feature_type": kb_entry.get('feature_type'),
                        "feature_class": props.get('feature_class'),
                        "subclass": props.get('subclass'),
                        "host_scope": props.get('host_scope', []),
                        "delivery_scope": props.get('delivery_scope', []),
                        "descriptions": source.get('descriptions', []),
                        "annotation_source": source.get('annotation_source'),
                        "polymerase_class": props.get('polymerase_class'),
                        "orientation_requirements": props.get('orientation_requirements'),
                        "frame_semantics": props.get('frame_semantics'),
                    }

            annotations.append({
                "name": feature_name,
                "type": str(row.get("Type", "misc_feature")),
                "start": max(0, qstart - 1),
                "end": qend,
                "direction": direction,
                "color": "#7C3AED",
                "source": "pLannotate",
                "layer": "feature",
                "description": description,
                "sseqid": sseqid,  # Keep original sseqid for reference
                "db": db_name,     # Database source
                "kb_data": kb_data,  # Full KB metadata for gene card
            })
        except Exception:
            continue
    return annotations


@router.post("/annotate_sequence")
async def annotate_sequence_endpoint(request: AnnotateSequenceRequest):
    """
    Annotate a raw DNA sequence with pLannotate.
    Returns SeqViz-compatible annotations (0-based start, 0-based end-exclusive).
    """
    if not PLANNOTATE_AVAILABLE:
        raise HTTPException(
            status_code=500,
            detail={"error": "pLannotate not available", "details": PLANNOTATE_ERROR},
        )

    if not request.sequence:
        return {"ok": False, "annotations": [], "error": "No sequence provided"}

    if request.hierarchical:
        return await annotate_sequence_with_hierarchy_endpoint(request)

    try:
        linear = not request.circular
        # strip whitespace — same Biopython 1.86 SeqIO.write assertion
        # that broke annotate_sequence_llm_endpoint hits this call site too.
        request_sequence_clean = re.sub(r'\s+', '', request.sequence)
        df = annotate(inSeq=request_sequence_clean, linear=linear, is_detailed=request.detailed)
        # Deduplicate annotations
        if df is not None and not df.empty:
            df = _deduplicate_annotations(df)
        rows = df.to_dict(orient="records") if df is not None and not df.empty else []
        annotations = _plannotate_rows_to_seqviz_annotations(rows)

        return {
            "ok": True,
            "annotations": annotations,
            "summary": {
                "plannotate_feature_count": len(annotations),
                "detailed": request.detailed,
                "source_endpoint": "annotate_sequence",
            },
        }

    except Exception as e:
        import traceback
        return {
            "ok": False,
            "annotations": [],
            "error": str(e),
            "traceback": traceback.format_exc(),
        }




async def _get_gpt_intent_analysis(features: list, graph: dict, sequence_length: int, circular: bool) -> dict:
    """
    Use GPT-4o-mini to generate natural language analysis of plasmid intent.
    Returns dict with 'gpt_analysis' or 'error'/'fallback' if unavailable.
    """
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key or not features:
        return {"fallback": "OpenAI API key not configured or no features detected"}
    
    try:
        # Build annotations list for GPT analysis (limit to first 30)
        features_for_gpt = features[:30]
        gpt_annotations = [
            {
                "name": f.get("feature_name", "unknown"),
                "start": f.get("start", 0),
                "end": f.get("end", 0),
                "role": f.get("role"),
                "lentiviral_context": f.get("lentiviral_context"),
            }
            for f in features_for_gpt
        ]
        
        # Build expression cassettes summary
        cassettes = graph.get("expression_cassettes", [])
        cassette_summary = ""
        if cassettes:
            cassette_lines = []
            for c in cassettes[:10]:
                prom = c.get("promoter_name", "?")
                prod = c.get("product_name", "?")
                prod_type = c.get("product_type", "?")
                polya = c.get("polya_name", "None")
                ctype = c.get("type", "unknown")
                cassette_lines.append(f"  - [{ctype}] {prom} -> {prod} ({prod_type}) -> {polya}")
            cassette_summary = "\n\nExpression cassettes detected:\n" + "\n".join(cassette_lines)
        
        # Sort annotations by position
        sorted_anns = sorted(gpt_annotations, key=lambda a: a.get("start", 0))
        
        topology = "circular" if circular else "linear"
        length_str = f"{sequence_length:,} bp"
        
        feature_lines = "\n".join(
            f"  - {a['name']} ({a.get('start')}..{a.get('end')}, role: {a.get('role', 'unknown')}, context: {a.get('lentiviral_context', 'N/A')})"
            for a in sorted_anns
        )
        
        prompt = (
            f"You are a molecular biology expert. Analyze the following plasmid map and "
            f"write a clear, concise explanation of its biological intent.\n\n"
            f"Plasmid: Annotated Plasmid\n"
            f"Topology: {topology}\n"
            f"Total length: {length_str}\n\n"
            f"Features (in order):\n{feature_lines}"
            f"{cassette_summary}\n\n"
            f"Write 2-3 paragraphs covering:\n"
            f"1. What type of vector this is and its primary purpose\n"
            f"2. How the expression cassette(s) work - what proteins or RNAs are produced "
            f"and what drives their expression\n"
            f"3. Any notable functional features (selection, viral elements, delivery mechanism)\n\n"
            f"Be specific about the elements present. Write in plain English for a molecular "
            f"biologist."
        )
        
        from openai import AsyncOpenAI
        client = AsyncOpenAI(api_key=api_key)
        response = await client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[{"role": "user", "content": prompt}],
            max_tokens=800,
            temperature=0.3,
        )
        gpt_analysis = response.choices[0].message.content.strip()
        
        return {
            "gpt_analysis": gpt_analysis,
            "model": "gpt-4o-mini",
            "features_analyzed": len(gpt_annotations),
            "cassettes_detected": len(cassettes),
        }
    except Exception as e:
        print(f"[DEBUG] GPT analysis failed: {e}")
        return {"error": str(e), "fallback": "GPT analysis unavailable"}


@router.post("/annotate_sequence_with_hierarchy")
async def annotate_sequence_with_hierarchy_endpoint(request: AnnotateSequenceRequest):
    """
    Annotate a raw DNA sequence with pLannotate and a higher-order grammar layer.
    Returns a combined SeqViz-compatible annotation list plus module metadata.
    """
    if not PLANNOTATE_AVAILABLE:
        raise HTTPException(
            status_code=500,
            detail={"error": "pLannotate not available", "details": PLANNOTATE_ERROR},
        )

    if not request.sequence:
        return {"ok": False, "annotations": [], "error": "No sequence provided"}

    try:
        print("[DEBUG] Step 1: Starting annotation...")
        linear = not request.circular
        # strip whitespace — same Biopython 1.86 SeqIO.write assertion
        # that broke annotate_sequence_llm_endpoint hits this call site too.
        request_sequence_clean = re.sub(r'\s+', '', request.sequence)
        df = annotate(inSeq=request_sequence_clean, linear=linear, is_detailed=request.detailed)
        print(f"[DEBUG] Step 2: pLannotate returned df with shape: {df.shape if df is not None else 'None'}")
        
        # Deduplicate annotations
        if df is not None and not df.empty:
            df = _deduplicate_annotations(df)
            print(f"[DEBUG] Step 2b: After dedup: {df.shape}")
        
        rows = df.to_dict(orient="records") if df is not None and not df.empty else []
        print(f"[DEBUG] Step 3: Converted to {len(rows)} rows")
        
        plannotate_annotations = _plannotate_rows_to_seqviz_annotations(rows)
        plannotate_annotations = _dedup_cross_tier_annotations(plannotate_annotations)
        print(f"[DEBUG] Step 4: Created {len(plannotate_annotations)} plannotate annotations")
        
        hierarchy = annotate_hierarchy_from_plannotate_v2(
            sequence=request.sequence,
            circular=request.circular,
            plannotate_rows=rows,
            )
        print(f"[DEBUG] Step 5: Hierarchy created successfully")
        
        # Apply boundary corrections to plannotate annotations
        boundary_corrections = hierarchy.get("boundary_corrections", [])
        if boundary_corrections:
            correction_map = {
                corr["feature_name"].lower(): corr
                for corr in boundary_corrections
            }
            for ann in plannotate_annotations:
                ann_name = (ann.get("name") or "").lower()
                if ann_name in correction_map:
                    corr = correction_map[ann_name]
                    # Update to corrected boundaries (SeqViz uses 0-based start)
                    ann["start"] = corr["corrected_start"]
                    ann["end"] = corr["corrected_end"]
                    ann["boundary_corrected"] = True
                    print(f"[DEBUG] Corrected {ann_name}: {corr['original_start']}-{corr['original_end']} -> {corr['corrected_start']}-{corr['corrected_end']}")
            print(f"[DEBUG] Applied {len(boundary_corrections)} boundary corrections")
        
        module_anns = hierarchy.get("module_annotations", [])
        gap_anns = hierarchy.get("gap_annotations", [])
        motif_anns = hierarchy.get("motif_annotations", [])
        print(f"[DEBUG] Hierarchical: {len(module_anns)} modules, {len(gap_anns)} gaps, {len(motif_anns)} motifs")
        if motif_anns:
            print(f"[DEBUG] Sample motifs: {motif_anns[:3]}")
        hierarchical_annotations = [
            *module_anns,
            *gap_anns,
            *motif_anns,
        ]
        combined = plannotate_annotations  # hierarchical merged by frontend

        # Note: GPT intent analysis is NOT run automatically during annotation
        # Use the separate /analyze_plasmid endpoint to get GPT-powered analysis

        return {
            "ok": True,
            "annotations": combined,
            "plannotate_annotations": plannotate_annotations,
            "all_hierarchical": hierarchical_annotations,
            "modules": hierarchy.get("modules", []),
            "feature_instances": hierarchy.get("feature_instances", []),
            "module_instances": hierarchy.get("module_instances", []),
            "junctions": hierarchy.get("junctions", []),
            "build_profiles": hierarchy.get("build_profiles", []),
            "construct_graph": hierarchy.get("construct_graph", {}),
            "module_graph": hierarchy.get("module_graph", {}),  # For download
            "capabilities": hierarchy.get("capabilities", []),
            "rule_findings": hierarchy.get("rule_findings", []),
            "boundary_corrections": boundary_corrections,
            "summary": {
                "plannotate_feature_count": len(plannotate_annotations),
                "boundary_corrections_applied": len(boundary_corrections),
                "hierarchical_annotation_count": len(hierarchical_annotations),
                "module_count": hierarchy.get("module_count", 0),
                "features_inferred": hierarchy.get("features_inferred", 0),
                "feature_instance_count": hierarchy.get("analyzer_summary", {}).get("feature_instance_count", 0),
                "module_instance_count": hierarchy.get("analyzer_summary", {}).get("module_instance_count", 0),
                "junction_count": hierarchy.get("analyzer_summary", {}).get("junction_count", 0),
                "build_profile_count": hierarchy.get("analyzer_summary", {}).get("build_profile_count", 0),
                "capability_count": hierarchy.get("analyzer_summary", {}).get("capability_count", 0),
                "rule_finding_count": hierarchy.get("analyzer_summary", {}).get("rule_finding_count", 0),
                "detailed": request.detailed,
                "source_endpoint": "annotate_sequence_with_hierarchy",
            },
        }
    except Exception as e:
        import traceback
        import sys
        exc_type, exc_value, exc_tb = sys.exc_info()
        print(f"❌ ERROR in annotate_sequence_with_hierarchy:")
        print(f"   Type: {exc_type.__name__}")
        print(f"   Message: {str(e)}")
        print(f"   File: {exc_tb.tb_frame.f_code.co_filename}:{exc_tb.tb_lineno}")
        print("   Full traceback:")
        traceback.print_exc()
        import traceback
        return {
            "ok": False,
            "annotations": [],
            "error": str(e),
            "traceback": traceback.format_exc(),
        }


# =============================================================================

# =============================================================================
# HEURISTIC ANNOTATION ENDPOINT
# =============================================================================

from .heuristic_motif_detector import HeuristicMotifDetector, MotifHit
from .heuristic_scorer import HeuristicScorer, ModuleCall, RuleFiring


@router.post("/annotate_sequence_heuristic")
async def annotate_sequence_heuristic_endpoint(request: AnnotateSequenceRequest):
    """
    Annotate plasmid sequence using heuristic-based module detection.
    
    Pipeline:
    1. Run pLannotate for feature annotation (reuse existing)
    2. Run heuristic motif detection (ORFs, Kozak, RBS)
    3. Score CDS modules using heuristic rules
    4. Score expression modules using heuristic rules
    5. Score vector type using heuristic rules
    6. Return annotated features with scores and rules fired
    """
    if not PLANNOTATE_AVAILABLE:
        raise HTTPException(
            status_code=500,
            detail={"error": "pLannotate not available", "details": PLANNOTATE_ERROR},
        )

    if not request.sequence:
        return {"ok": False, "annotations": [], "error": "No sequence provided"}

    try:
        sequence = request.sequence.upper()
        circular = request.circular if request.circular is not None else True
        linear = not circular
        
        print("[HEURISTIC] Step 1: Running pLannotate...")
        # Step 1: Run pLannotate (reuse existing function)
        # Strip whitespace from the input sequence — modern Biopython's
        # SeqIO.write asserts no newlines or spaces inside a FASTA record
        # (Bio/SeqIO/FastaIO.py write_record), which pLannotate trips on
        # for any GenBank-derived input where the sequence string still
        # carries its source-line newlines. Without this strip, annotate()
        # raises AssertionError, the endpoint catches it and silently
        # returns 0 features — every annotation-driven flow then breaks.
        sequence = re.sub(r'\s+', '', sequence)
        df = annotate(inSeq=sequence, linear=linear, is_detailed=request.detailed)
        
        # Deduplicate annotations
        if df is not None and not df.empty:
            df = _deduplicate_annotations(df)
        
        rows = df.to_dict(orient="records") if df is not None and not df.empty else []
        plannotate_annotations = _plannotate_rows_to_seqviz_annotations(rows)
        plannotate_annotations = _dedup_cross_tier_annotations(plannotate_annotations)
        print(f"[HEURISTIC] Step 1 complete: {len(plannotate_annotations)} pLannotate features")
        
        # Build feature list for scoring
        plannotate_features = []
        for ann in plannotate_annotations:
            plannotate_features.append({
                'name': ann.get('name', ''),
                'Feature': ann.get('name', ''),
                'type': ann.get('type', ''),
                'Type': ann.get('type', ''),
                'start': ann.get('start', 0),
                'end': ann.get('end', 0),
                'strand': ann.get('direction', 1),
                'sseqid': ann.get('sseqid', ''),
                'description': ann.get('description', ''),
            })
        
        # Step 2: Run heuristic motif detection
        print("[HEURISTIC] Step 2: Running motif detection...")
        motif_detector = HeuristicMotifDetector(sequence, circular)
        motifs = motif_detector.detect_all_motifs(plannotate_features)
        print(f"[HEURISTIC] Step 2 complete: {len(motifs.get('orfs', []))} ORFs, {len(motifs.get('kozak', []))} Kozak, {len(motifs.get('rbs', []))} RBS")
        
        # Step 3: Score modules using heuristics
        print("[HEURISTIC] Step 3: Scoring modules...")
        scorer = HeuristicScorer()
        cds_modules = scorer.score_cds_modules(motifs, plannotate_features, sequence)
        expression_modules = scorer.score_expression_modules(cds_modules, plannotate_features, sequence)
        vector_calls = scorer.score_vector_type(expression_modules, plannotate_features, sequence)
        print(f"[HEURISTIC] Step 3 complete: {len(cds_modules)} CDS modules, {len(expression_modules)} expression modules, {len(vector_calls)} vector calls")
        
        # Step 4: Rank and organize results
        ranked = scorer.rank_modules(cds_modules, expression_modules, vector_calls)
        
        # Step 5: Build nested module structure
        nested_modules = _build_nested_modules_heuristic(
            ranked['cds_modules'],
            ranked['expression_modules'],
            ranked['vector_type'],
            len(sequence)
        )
        
        # Step 6: Build hierarchical annotations for visualization
        hierarchical_annotations = []
        
        # Add CDS modules as annotations
        for m in ranked['cds_modules']:
            hierarchical_annotations.append({
                'name': m.metadata.get('feature_name') or f"CDS Module ({m.metadata.get('orf_length_aa', 0)} aa)",
                'type': 'module',
                'start': m.start,
                'end': m.end,
                'direction': m.strand,
                'color': '#4CAF50' if m.confidence == 'high' else '#FFC107' if m.confidence == 'medium' else '#F44336',
                'layer': 'module',
                'module_type': m.module_type,
                'score': m.score,
                'confidence': m.confidence,
                'rules_fired': [r.to_dict() for r in m.rules_fired],
                'metadata': m.metadata,
            })
        
        # Add expression modules
        for m in ranked['expression_modules']:
            exp_type = m.module_type.replace('_', ' ').title()
            hierarchical_annotations.append({
                'name': exp_type,
                'type': 'module',
                'start': m.start,
                'end': m.end,
                'direction': m.strand,
                'color': '#2196F3' if m.confidence == 'high' else '#03A9F4' if m.confidence == 'medium' else '#B3E5FC',
                'layer': 'module',
                'module_type': m.module_type,
                'score': m.score,
                'confidence': m.confidence,
                'rules_fired': [r.to_dict() for r in m.rules_fired],
                'metadata': m.metadata,
            })
        
        # Add vector type annotations
        for m in ranked['vector_type']:
            if m.score > 0.5:
                vec_type = m.module_type.replace('_', ' ').title()
                hierarchical_annotations.append({
                    'name': vec_type,
                    'type': 'vector',
                    'start': m.start,
                    'end': m.end,
                    'direction': m.strand,
                    'color': '#9C27B0' if m.confidence == 'high' else '#CE93D8',
                    'layer': 'module',
                    'module_type': m.module_type,
                    'score': m.score,
                    'confidence': m.confidence,
                    'rules_fired': [r.to_dict() for r in m.rules_fired],
                    'metadata': m.metadata,
                })
        
        # Determine top vector type
        top_vector = ranked['vector_type'][0].module_type if ranked['vector_type'] else "unknown"
        
        return {
            "ok": True,
            "method": "heuristic",
            "plannotate_annotations": plannotate_annotations,
            "annotations": plannotate_annotations,  # For compatibility
            "heuristic_motifs": {
                "orfs": [m.to_dict() for m in motifs.get('orfs', [])],
                "kozak": [m.to_dict() for m in motifs.get('kozak', [])],
                "rbs": [m.to_dict() for m in motifs.get('rbs', [])],
                "start_codons": [m.to_dict() for m in motifs.get('start_codons', [])],
                "stop_codons": [m.to_dict() for m in motifs.get('stop_codons', [])],
            },
            "cds_modules": [m.to_dict() for m in ranked['cds_modules']],
            "expression_modules": [m.to_dict() for m in ranked['expression_modules']],
            "vector_type": [m.to_dict() for m in ranked['vector_type']],
            "nested_modules": nested_modules,
            "hierarchical_annotations": all_hierarchical,
            "summary": {
                "plannotate_feature_count": len(plannotate_annotations),
                "orf_count": len(motifs.get('orfs', [])),
                "kozak_count": len(motifs.get('kozak', [])),
                "rbs_count": len(motifs.get('rbs', [])),
                "cds_module_count": len(ranked['cds_modules']),
                "expression_module_count": len(ranked['expression_modules']),
                "top_vector_type": top_vector,
                "source_endpoint": "annotate_sequence_heuristic"
            }
        }
        
    except Exception as e:
        import traceback
        import sys
        exc_type, exc_value, exc_tb = sys.exc_info()
        print(f"❌ ERROR in annotate_sequence_heuristic:")
        print(f"   Type: {exc_type.__name__}")
        print(f"   Message: {str(e)}")
        traceback.print_exc()
        return {
            "ok": False,
            "annotations": [],
            "error": str(e),
            "traceback": traceback.format_exc(),
        }


def _build_nested_modules_heuristic(
    cds_modules: list,
    expression_modules: list,
    vector_calls: list,
    seq_len: int
) -> list:
    """
    Build nested module structure following hierarchical nesting rules.
    
    Nesting hierarchy:
    1. Vector level (lentiviral_payload, aav_payload, etc.)
    2. Expression module level (pol2_expression, bacterial_expression)
    3. CDS module level (cds_module with features)
    """
    nested = []
    
    # Build all modules with level
    all_modules = []
    for m in cds_modules:
        all_modules.append({
            'level': 3,
            'data': m.to_dict(),
            'children': []
        })
    
    for m in expression_modules:
        all_modules.append({
            'level': 2,
            'data': m.to_dict(),
            'children': []
        })
    
    for m in vector_calls:
        if m.score > 0.5:
            all_modules.append({
                'level': 1,
                'data': m.to_dict(),
                'children': []
            })
    
    # Sort by level then start position
    all_modules.sort(key=lambda x: (x['level'], x['data']['start']))
    
    # Nest children inside parents
    for parent in all_modules:
        for potential_child in all_modules:
            if parent['level'] >= potential_child['level']:
                continue
            
            p_start = parent['data']['start']
            p_end = parent['data']['end']
            c_start = potential_child['data']['start']
            c_end = potential_child['data']['end']
            
            # Check if child is contained within parent
            if c_start >= p_start and c_end <= p_end:
                parent['children'].append(potential_child['data'])
    
    # Collect top-level modules (not nested)
    nested_ids = set()
    for item in all_modules:
        for child in item['children']:
            # Use module_type + start + end as identifier
            key = f"{child['module_type']}_{child['start']}_{child['end']}"
            nested_ids.add(key)
    
    for item in all_modules:
        key = f"{item['data']['module_type']}_{item['data']['start']}_{item['data']['end']}"
        if key not in nested_ids:
            result = item['data'].copy()
            result['children'] = item['children']
            result['layer'] = 'module'
            nested.append(result)
    
    return nested


# ---------------------------------------------------------------------------
# Import annotations from a user-supplied list (CSV-derived)
# ---------------------------------------------------------------------------

class ImportAnnotationEntry(BaseModel):
    name: str
    sequence: str
    type: Optional[str] = None
    location: Optional[str] = None
    length: Optional[str] = None
    description: Optional[str] = None


class ImportAnnotationsRequest(BaseModel):
    sequence: str
    circular: bool = False
    entries: List[ImportAnnotationEntry]
    max_mismatches: int = 0
    default_color: str = "#84B0DC"


def _revcomp(seq: str) -> str:
    table = str.maketrans("ATGCNatgcn", "TACGNtacgn")
    return seq.translate(table)[::-1]


def _scan_exact(haystack: str, needle: str) -> List[int]:
    hits: List[int] = []
    if not needle:
        return hits
    start = 0
    while True:
        i = haystack.find(needle, start)
        if i < 0:
            break
        hits.append(i)
        start = i + 1
    return hits


def _scan_mismatch(haystack: str, needle: str, max_mm: int) -> List[int]:
    if max_mm <= 0 or not needle or len(needle) > len(haystack):
        return []
    hits: List[int] = []
    n = len(needle)
    for i in range(len(haystack) - n + 1):
        mm = 0
        for j in range(n):
            if haystack[i + j] != needle[j]:
                mm += 1
                if mm > max_mm:
                    break
        if mm <= max_mm:
            hits.append(i)
    return hits


@router.post("/import_annotations")
async def import_annotations_endpoint(request: ImportAnnotationsRequest):
    """Map user-supplied (name, sequence) entries onto a plasmid by sequence
    identity. Returns SeqViz-shaped annotations for every match found on
    either strand (with circular wrap-around when applicable). Optional CSV
    metadata (type, location, length, description) is passed through to the
    annotation's gene-card popup via the same `kb_data` shape the frontend
    already renders for KB-derived features.
    """
    import re

    seq_raw = (request.sequence or "").upper()
    seq = re.sub(r"[^ATGCN]", "", seq_raw)
    if not seq:
        raise HTTPException(status_code=400, detail="Empty plasmid sequence")
    seq_len = len(seq)

    # Wrap-extended haystack for circular plasmids: searches that span the
    # origin still resolve via modulo arithmetic on the hit index.
    haystack_fwd = seq + (seq[: max(0, seq_len - 1)] if request.circular else "")
    haystack_rev = _revcomp(haystack_fwd)

    annotations: List[Dict[str, Any]] = []
    unmatched: List[Dict[str, Any]] = []
    n_input = len(request.entries)
    n_matched = 0

    valid_chars = set("ATGCN")

    for entry in request.entries:
        name = (entry.name or "").strip()
        if not name:
            unmatched.append({"name": "", "sequence": entry.sequence,
                              "reason": "missing name"})
            continue

        q_raw = (entry.sequence or "").upper().strip()
        q = "".join(ch for ch in q_raw if ch.isalpha())
        if not q:
            unmatched.append({"name": name, "sequence": entry.sequence,
                              "reason": "empty sequence"})
            continue
        bad = [ch for ch in q if ch not in valid_chars]
        if bad:
            unmatched.append({"name": name, "sequence": q,
                              "reason": f"non-ATGCN characters: {''.join(sorted(set(bad)))}"})
            continue
        if len(q) > seq_len:
            unmatched.append({"name": name, "sequence": q,
                              "reason": "query longer than plasmid"})
            continue

        q_rc = _revcomp(q)

        fwd_hits = _scan_exact(haystack_fwd, q)
        rev_hits = _scan_exact(haystack_fwd, q_rc)

        # Mismatch fallback only if no exact hits on either strand.
        if not fwd_hits and not rev_hits and request.max_mismatches > 0:
            fwd_hits = _scan_mismatch(haystack_fwd, q, request.max_mismatches)
            rev_hits = _scan_mismatch(haystack_fwd, q_rc, request.max_mismatches)

        # Normalise circular wrap: any hit i in [seq_len, len(haystack)) is a
        # wrap-around hit whose linear start is i % seq_len.
        def _make_ann(start_in_haystack: int, direction: int) -> Dict[str, Any]:
            start = start_in_haystack % seq_len
            end = (start + len(q)) % seq_len
            if end == 0:
                end = seq_len
            descriptions: List[str] = []
            if entry.description:
                descriptions.append(entry.description)
            if entry.location:
                descriptions.append(f"location: {entry.location}")
            if entry.length:
                descriptions.append(f"length: {entry.length}")
            ann: Dict[str, Any] = {
                "name": name,
                "start": start,
                "end": end,
                "direction": direction,
                "strand": direction,
                "color": request.default_color or "#84B0DC",
                "type": entry.type or "misc_feature",
                "layer": "feature",
                "description": entry.description or "",
                "source": "csv_import",
            }
            ann["kb_data"] = {
                "source_type": "feature_kb",
                "feature_name": name,
                "feature_type": entry.type or "misc_feature",
                "annotation_source": "CSV import",
                "descriptions": descriptions or [f"Imported from CSV: {name}"],
            }
            return ann

        # Emit at most ONE annotation per CSV entry. Forward strand wins ties;
        # reverse-strand hits are used only if there is no forward match. This
        # honors the "one feature per row" rule and avoids creating duplicate
        # annotations of differing types for the same logical entry.
        n_for_entry = 0
        if fwd_hits:
            annotations.append(_make_ann(fwd_hits[0], 1))
            n_for_entry = 1
        elif rev_hits:
            # Reverse-strand hits are still indexed against the forward-strand
            # haystack (we revcomp'd the query, not the plasmid), so the start
            # position is already the genomic start of the feature.
            annotations.append(_make_ann(rev_hits[0], -1))
            n_for_entry = 1

        if n_for_entry == 0:
            unmatched.append({"name": name, "sequence": q,
                              "reason": "no match on either strand"})
        else:
            n_matched += 1

    return {
        "ok": True,
        "annotations": annotations,
        "unmatched": unmatched,
        "summary": {
            "n_input": n_input,
            "n_matched": n_matched,
            "n_annotations": len(annotations),
            "n_unmatched": len(unmatched),
        },
    }


# ---------------------------------------------------------------------------
# Design primers (alias under /plannotate so reverse proxies that whitelist
# only the /plannotate/* path can reach it). Delegates to pcr.design_primers.
# ---------------------------------------------------------------------------
from .pcr import design_primers as _pcr_design_primers, PrimerRequest as _PcrPrimerRequest


@router.post("/design_primers")
def design_primers_alias(req: _PcrPrimerRequest):
    return _pcr_design_primers(req)


# ---------------------------------------------------------------------------
# Design CRISPR sgRNAs (deterministic Python re-implementation of the
# crisprVerse plasmid-scoped design path — pure Python, no R, no genome
# off-target search, no Bioconductor dependency).
# ---------------------------------------------------------------------------
from .guide_designer import design_guides as _design_guides


class GuideDesignRequest(BaseModel):
    sequence: str
    region_start: int        # 1-indexed inclusive
    region_end: int          # 1-indexed inclusive
    pam: str = "NGG"
    guide_length: int = 20
    pam_position: str = "3prime"   # "3prime" (Cas9) | "5prime" (Cas12a)
    max_guides: int = 50
    min_score: float = 0.0
    score_method: str = "doench2014"   # "doench2014" | "heuristic"


@router.post("/design_guides")
def design_guides_endpoint(req: GuideDesignRequest):
    return _design_guides(
        sequence=req.sequence,
        region_start=req.region_start,
        region_end=req.region_end,
        pam=req.pam,
        guide_length=req.guide_length,
        pam_position=req.pam_position,
        max_guides=req.max_guides,
        min_score=req.min_score,
        score_method=req.score_method,
    )

# ---------------------------------------------------------------------------
# Design pegRNAs for prime editing (easy_prime port: PE3 XGBoost model,
# Anzalone 2019 pegRNA architecture spacer-scaffold-RTT-PBS).
# ---------------------------------------------------------------------------
from .pegrna_designer import design_pegrnas as _design_pegrnas


class PegRNADesignRequest(BaseModel):
    sequence: str
    edit_start: int   # 1-indexed inclusive
    edit_end: int     # 1-indexed inclusive
    alt: str = ""
    edit_type: str = "substitution"   # substitution | insertion | deletion
    n_results: int = 3
    use_pe3: bool = True


@router.post("/design_pegrnas")
def design_pegrnas_endpoint(req: PegRNADesignRequest):
    return _design_pegrnas(
        sequence=req.sequence,
        edit_start_1based=req.edit_start,
        edit_end_1based=req.edit_end,
        alt=req.alt,
        edit_type=req.edit_type,
        n_results=req.n_results,
        use_pe3=req.use_pe3,
    )



# TOKENIZATION ENDPOINTS

@router.post("/annotate_pol2_grammar")
async def annotate_pol2_grammar_endpoint(request: AnnotateSequenceRequest):
    """
    Annotate plasmid using grammar-based Pol II cassette detection.
    
    This endpoint uses KB feature class queries instead of heuristic rules
    for detecting mammalian Pol II expression cassettes.
    
    Returns:
    - Detected Pol II cassettes with confidence scores
    - Component breakdown (promoters, polyA, enhancers, introns, etc.)
    - Compatible with existing frontend
    """
    if not PLANNOTATE_AVAILABLE:
        raise HTTPException(
            status_code=500,
            detail={"error": "pLannotate not available", "details": PLANNOTATE_ERROR},
        )

    if not request.sequence:
        return {"ok": False, "error": "No sequence provided"}

    try:
        sequence = request.sequence.upper()
        circular = request.circular if request.circular is not None else True
        linear = not circular
        
        print("[GRAMMAR] Step 1: Running pLannotate...")
        # Step 1: Run pLannotate
        # Strip whitespace from the input sequence — modern Biopython's
        # SeqIO.write asserts no newlines or spaces inside a FASTA record
        # (Bio/SeqIO/FastaIO.py write_record), which pLannotate trips on
        # for any GenBank-derived input where the sequence string still
        # carries its source-line newlines. Without this strip, annotate()
        # raises AssertionError, the endpoint catches it and silently
        # returns 0 features — every annotation-driven flow then breaks.
        sequence = re.sub(r'\s+', '', sequence)
        df = annotate(inSeq=sequence, linear=linear, is_detailed=request.detailed)
        
        # Deduplicate annotations
        if df is not None and not df.empty:
            df = _deduplicate_annotations(df)
        
        rows = df.to_dict(orient="records") if df is not None and not df.empty else []
        plannotate_annotations = _plannotate_rows_to_seqviz_annotations(rows)
        plannotate_annotations = _dedup_cross_tier_annotations(plannotate_annotations)
        print(f"[GRAMMAR] pLannotate found {len(plannotate_annotations)} features")
        
        # Step 2: Enrich with KB
        print("[GRAMMAR] Step 2: Enriching with KB...")
        plannotate_features = []
        for ann in plannotate_annotations:
            kb_info = _lookup_kb_full(ann.get('name', ''))
            feature = {
                'name': ann.get('name', ''),
                'start': ann.get('start', 0),
                'end': ann.get('end', 0),
                'strand': ann.get('direction', 1),
                'type': ann.get('type', 'misc_feature'),
                'kb_info': kb_info
            }
            plannotate_features.append(feature)
        
        print(f"[GRAMMAR] KB enrichment complete")
        
        # Step 3: Detect Pol II cassettes using grammar
        print("[GRAMMAR] Step 3: Detecting Pol II cassettes...")
        pol2_cassettes = detect_pol2_cassettes(plannotate_features, sequence, circular)
        print(f"[GRAMMAR] Detected {len(pol2_cassettes)} Pol II cassette(s)")
        
        # Step 4: Format response
        hierarchical_annotations = []
        for cassette in pol2_cassettes:
            # Add module annotation
            module_ann = {
                'name': f"{cassette['promoter']['type']} expression cassette",
                'start': cassette['start'],
                'end': cassette['end'],
                'direction': cassette['strand'],
                'type': 'module',
                'layer': 'module',
                'module_type': cassette['module_type'],
                'confidence': cassette['confidence'],
                'score': cassette['score'],
                'promoter_type': cassette['promoter']['type'],
                'polya_type': cassette['polya']['type'],
                'components': {
                    'cds_count': len(cassette['components']['cds']),
                    'has_kozak': cassette['components']['kozak'] is not None,
                    'has_wpre': cassette['components']['wpre'] is not None,
                    'intron_count': len(cassette['components']['introns']),
                    '2a_count': len(cassette['components']['2a_peptides'])
                }
            }
            hierarchical_annotations.append(module_ann)
        
        return {
            "ok": True,
            "method": "grammar_pol2",
            "annotations": plannotate_annotations,
            "hierarchical_annotations": all_hierarchical,
            "pol2_cassettes": pol2_cassettes,
            "summary": {
                "plannotate_feature_count": len(plannotate_annotations),
                "pol2_cassette_count": len(pol2_cassettes),
                "detection_method": "KB grammar-based",
                "total_annotations": len(plannotate_annotations) + len(hierarchical_annotations)
            }
        }
        
    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        print(f"[GRAMMAR] Error: {error_details}")
        return {
            "ok": False,
            "error": str(e),
            "details": error_details
        }


# =============================================================================

from .tokenizer import (
    tokenize_plasmid,
    CollapsedSpan,
    TokenizationResult,
    map_feature_to_token,
)


class TokenizationRequest(BaseModel):
    """Request model for tokenization endpoint."""
    sequence: str
    circular: bool = True
    plasmid_id: str = "unknown"
    level: str = "functional"  # 'module', 'functional', 'full'


class TokenizationResponse(BaseModel):
    """Response model for tokenization endpoint."""
    ok: bool
    plasmid_id: str
    sequence_length: int
    circular: bool
    spans: List[Dict[str, Any]]
    tokens: List[str]
    validation: Dict[str, Any]
    error: Optional[str] = None


@router.post("/tokenize", response_model=TokenizationResponse)
async def tokenize_sequence(request: TokenizationRequest):
    """
    Tokenize a plasmid sequence.
    
    This endpoint:
    1. Runs pLannotate annotation
    2. Extracts hierarchical modules  
    3. Collapses to single-layer representation
    4. Returns token sequence
    
    Args:
        request: TokenizationRequest with sequence and options
        
    Returns:
        TokenizationResponse with collapsed spans and tokens
    """
    try:
        # Run annotation pipeline first
        sequence = request.sequence.upper()
        seq_len = len(sequence)
        
        # Run pLannotate
        db_dir = _get_db_dir()
        _, df = run_plannotate(
            sequence=sequence,
            is_linear=not request.circular,
            yaml_file=str(db_dir / "databases.yaml"),
            detailed=True,
        )
        
        # Deduplicate
        if df is not None and not df.empty:
            df = _deduplicate_annotations(df)
        
        rows = df.to_dict(orient="records") if df is not None and not df.empty else []
        
        # Get hierarchical modules
        hierarchy = annotate_hierarchy_from_plannotate_v2(
            sequence=sequence,
            circular=request.circular,
            plannotate_rows=rows,
            )
        
        # Extract modules and features for tokenization
        modules = []
        features = []
        
        # Parse module annotations into tokenizer format
        for ann in hierarchy.get("module_annotations", []):
            label = ann.get("name", "")
            modules.append({
                "module_id": f"{ann.get('start', 0)}_{ann.get('end', 0)}",
                "module_type": label.lower(),
                "label": label,
                "start": ann.get("start", 0),
                "end": ann.get("end", 0),
                "strand": 1,
                "payload_id": None,
            })
        
        # Parse plannotate annotations as features
        for row in rows:
            features.append({
                "feature_id": f"{row.get('qstart', 0)}_{row.get('qend', 0)}",
                "canonical_id": row.get("sseqid", ""),
                "canonical_type": row.get("type", ""),
                "feature_class": row.get("feature_class", ""),
                "label": row.get("Feature", ""),
                "start": row.get("qstart", 0),
                "end": row.get("qend", 0),
                "strand": 1 if row.get("strand", 1) >= 0 else -1,
            })
        
        # Run tokenization
        result = tokenize_plasmid(
            modules=modules,
            features=features,
            sequence_length=seq_len,
            plasmid_id=request.plasmid_id,
            circular=request.circular,
            level=request.level,
        )
        
        return TokenizationResponse(
            ok=True,
            plasmid_id=result.plasmid_id,
            sequence_length=result.sequence_length,
            circular=result.circular,
            spans=[s.to_dict() for s in result.spans],
            tokens=result.tokens,
            validation=result.validation,
        )
        
    except Exception as e:
        import traceback
        return TokenizationResponse(
            ok=False,
            plasmid_id=request.plasmid_id,
            sequence_length=len(request.sequence),
            circular=request.circular,
            spans=[],
            tokens=[],
            validation={"valid": False, "issues": [str(e)]},
            error=traceback.format_exc(),
        )


@router.post("/tokenize_annotated")
async def tokenize_from_annotations(
    modules: List[Dict[str, Any]],
    features: List[Dict[str, Any]],
    sequence_length: int,
    plasmid_id: str = "unknown",
    circular: bool = True,
    level: str = "functional",
):
    """
    Tokenize from pre-computed annotations.
    
    Use this endpoint when you already have module annotations
    (e.g., from /annotate_sequence_with_hierarchy).
    """
    try:
        result = tokenize_plasmid(
            modules=modules,
            features=features,
            sequence_length=sequence_length,
            plasmid_id=plasmid_id,
            circular=circular,
            level=level,
        )
        
        return {
            "ok": True,
            "plasmid_id": result.plasmid_id,
            "sequence_length": result.sequence_length,
            "circular": result.circular,
            "spans": [s.to_dict() for s in result.spans],
            "tokens": result.tokens,
            "validation": result.validation,
        }
        
    except Exception as e:
        import traceback
        return {
            "ok": False,
            "error": str(e),
            "traceback": traceback.format_exc(),
        }


@router.get("/tokens/vocabulary")
async def get_token_vocabulary():
    """
    Get the current token vocabulary.
    
    Returns information about token classes and common tokens.
    """
    from .tokenizer import (
        MODULE_PRIORITY,
        PATH_NAME_MAP,
    )
    
    return {
        "ok": True,
        "module_types": list(MODULE_PRIORITY.keys()),
        "path_name_map": PATH_NAME_MAP,
        "token_classes": [
            "MODULE", "GAP", "PROMOTER", "CDS", "POLYA", "ORI",
            "MARKER", "NLS", "LINKER", "TAG", "LENTI", "FEATURE",
        ],
    }


# =============================================================================
# LLM-BASED MODULE ANNOTATION ENDPOINT
# =============================================================================

@router.post("/annotate_sequence_llm")
async def annotate_sequence_llm_endpoint(request: AnnotateSequenceRequest):
    """
    Annotate plasmid using hybrid LLM + grammar CDS parsing pipeline.

    Pipeline:
    1. Run pLannotate for feature annotation
    2. Run CDS submodule parsing (grammar pipeline logic)
    3. Send features + submodules to LLM for module identification
    4. Return hierarchical modules + collapsed representation
    """
    if not PLANNOTATE_AVAILABLE:
        raise HTTPException(
            status_code=500,
            detail={"error": "pLannotate not available", "details": PLANNOTATE_ERROR},
        )

    if not request.sequence:
        return {"ok": False, "annotations": [], "error": "No sequence provided"}

    try:
        sequence = request.sequence.upper()
        circular = request.circular if request.circular is not None else True
        linear = not circular

        print("[LLM] Step 1: Running pLannotate...")
        # Step 1: Run pLannotate
        # Strip whitespace from the input sequence — modern Biopython's
        # SeqIO.write asserts no newlines or spaces inside a FASTA record
        # (Bio/SeqIO/FastaIO.py write_record), which pLannotate trips on
        # for any GenBank-derived input where the sequence string still
        # carries its source-line newlines. Without this strip, annotate()
        # raises AssertionError, the endpoint catches it and silently
        # returns 0 features — every annotation-driven flow then breaks.
        sequence = re.sub(r'\s+', '', sequence)
        df = annotate(inSeq=sequence, linear=linear, is_detailed=request.detailed)

        if df is not None and not df.empty:
            df = _deduplicate_annotations(df)

        rows = df.to_dict(orient="records") if df is not None and not df.empty else []
        plannotate_annotations = _plannotate_rows_to_seqviz_annotations(rows)
        plannotate_annotations = _dedup_cross_tier_annotations(plannotate_annotations)
        print(f"[LLM] Step 1 complete: {len(plannotate_annotations)} pLannotate features")

        # Step 2: Detect ORFs and parse CDS submodules
        print("[LLM] Step 2: Detecting ORFs (>150aa with start/stop codons)...")
        from .Module_Library_gb.module_extractor import resolve_cds_submodules, Module, sha256_text, seq_hash
        from .orf_finder import find_orfs
        
        # Find ORFs >150aa with ATG start and stop codons
        detected_orfs = find_orfs(sequence, min_aa_length=150)
        print(f"[LLM] Found {len(detected_orfs)} ORFs >150aa")
        
        # Parse each ORF into CDS submodules
        cds_submodules_list = []
        orf_modules = []
        
        for orf_idx, orf in enumerate(detected_orfs):
            cds_start = orf['start']
            cds_end = orf['end']
            cds_strand = orf['strand']
            cds_seq = sequence[cds_start:cds_end]
            cds_length = cds_end - cds_start
            
            # Create Module object for this ORF
            cds_module = Module(
                id=f"orf_{orf_idx}_{cds_start}_{cds_end}",
                plasmid_id="temp_plasmid",
                module_type="cds_module",
                payload_id=None,
                start=cds_start,
                end=cds_end,
                wraps=False,
                length=cds_length,
                sequence=cds_seq,
                seq_hash=sha256_text(cds_seq)[:24],
                end_inferred=False,
                metadata={'strand': cds_strand, 'orf_detected': True, 'aa_length': orf['aa_length']},
                features=[]
            )
            
            orf_modules.append({
                'start': cds_start,
                'end': cds_end,
                'strand': cds_strand,
                'aa_length': orf['aa_length']
            })
            
            try:
                # Convert plannotate annotations to Feature-like objects
                class SimpleFeature:
                    def __init__(self, d):
                        self.start = d['start']
                        self.end = d['end']
                        self.name = d.get('name', '')
                        self.canonical_id = d.get('sseqid', '')
                        self.canonical_type = d.get('type', '')
                        self.kb_feature_class = d.get('kb_data', {}).get('feature_class', '') if d.get('kb_data') else ''
                
                simple_features = [SimpleFeature(a) for a in plannotate_annotations]
                
                # Run CDS submodule resolution
                result = resolve_cds_submodules(
                    cds_module, sequence, simple_features,
                    "temp_plasmid", Module, sha256_text, seq_hash
                )
                
                # Convert Module objects to dicts with proper strand
                for submod in result["submodules"]:
                    cds_submodules_list.append({
                        "module_type": submod.module_type,
                        "payload_id": submod.payload_id,
                        "start": submod.start,
                        "end": submod.end,
                        "strand": cds_strand,  # Preserve ORF strand
                        "metadata": submod.metadata
                    })
            except Exception as e:
                print(f"[WARN] CDS submodule resolution failed for ORF {cds_start}-{cds_end}: {e}")
        
        print(f"[LLM] Step 2 complete: {len(cds_submodules_list)} CDS submodules from {len(detected_orfs)} ORFs")

        # Step 3: Deterministic rule-based module detection
        # (REPLACES former LLM-based step — rule schema is authoritative).
        print("[ANN] Step 3a: RuleBasedModuleDetector...")
        from .rule_based_module_detector import RuleBasedModuleDetector
        _heuristics_csv = str(
            pathlib.Path(__file__).parent / "heuristics" / "module_heuristics.csv"
        )
        _rule_detector = RuleBasedModuleDetector(_heuristics_csv)
        rule_based_modules = _rule_detector.detect_modules(plannotate_annotations, sequence)
        print(f"[ANN] Step 3a complete: {len(rule_based_modules)} rule-based modules")

        print("[ANN] Step 3b: MammalianPol2 detector...")
        from .mammalian_pol2_detector import detect_mammalian_pol2_cassettes
        pol2_cassettes, _filtered_orfs = detect_mammalian_pol2_cassettes(
            features=plannotate_annotations,
            sequence=sequence,
            orfs=orf_modules,
            circular=circular,
        )
        print(f"[ANN] Step 3b complete: {len(pol2_cassettes)} Pol II cassettes")

        print("[ANN] Step 3c: interaction_builder...")
        from .interaction_builder import build_interactions
        interactions = build_interactions(
            features=plannotate_annotations,
            rule_based_modules=rule_based_modules,
            cds_submodules=cds_submodules_list,
            mammalian_pol2_cassettes=pol2_cassettes,
        )
        print(f"[ANN] Step 3c complete: {len(interactions)} interactions")

        # Combined modules list (rule-based + Pol II cassettes).
        modules = list(rule_based_modules) + list(pol2_cassettes)

        # ----------------------------------------------------------------
        # Build module-level hierarchical_annotations from the deterministic
        # detector outputs. No LLM, no fallback.
        # ----------------------------------------------------------------
        def _module_color(mt: str) -> str:
            return {
                "mammalian_pol2_expression_cassette": "#4CAF50",
                "pol2_expression_cassette":           "#4CAF50",
                "upstream_regulatory_module":         "#8BC34A",
                "downstream_regulatory_module":       "#CDDC39",
                "cds_module":                         "#673AB7",
                "bacterial_ori":                      "#607D8B",
                "bacterial_selection_cassette":       "#FF9800",
                "bacterial_selection":                "#FF9800",
                "bacterial_backbone":                 "#795548",
                "pol3_expression_cassette":           "#3F51B5",
                "lentiviral_payload":                 "#E91E63",
                "lentiviral_upstream_regulatory":     "#03A9F4",
                "lentiviral_downstream_regulatory":   "#9C27B0",
                "lac_alpha_blue_white_module":        "#FFEB3B",
                "gateway_dest_cassette":              "#00BCD4",
                "floxed_region":                      "#F44336",
                "frt_flanked_region":                 "#E91E63",
            }.get(mt, "#607D8B")

        def _human_module_name(mt: str) -> str:
            return mt.replace("_", " ").title()

        hierarchical_annotations = []
        for mod in rule_based_modules:
            hierarchical_annotations.append({
                "name":            mod.get("name") or _human_module_name(mod.get("module_type", "module")),
                "start":           mod.get("start"),
                "end":             mod.get("end"),
                "direction":       mod.get("strand", 1),
                "color":           _module_color(mod.get("module_type", "")),
                "layer":           "module",
                "module_type":     mod.get("module_type"),
                "module_id":       mod.get("module_id"),
                "source":          "rule_based_detector",
                "metadata":        mod.get("metadata", {}),
                "parent_module":   mod.get("parent_module"),
                "features":        mod.get("features", []),
                "nested_modules":  mod.get("nested_modules", []),
            })

        # Flatten each Pol II cassette into 4 hierarchical entries
        # (cassette + upstream + cds_module + downstream).
        for cidx, cas in enumerate(pol2_cassettes):
            cassette_id = f"pol2_cassette_{cidx + 1:02d}"
            hierarchical_annotations.append({
                "name":         "Pol2 Expression Cassette",
                "start":        cas.get("start"),
                "end":          cas.get("end"),
                "direction":    cas.get("strand", 1),
                "color":        _module_color("pol2_expression_cassette"),
                "layer":        "module",
                "module_type":  "pol2_expression_cassette",
                "module_id":    cassette_id,
                "source":       "mammalian_pol2_detector",
                "metadata":     {"detection_method": cas.get("detection_method"), "weight": cas.get("weight")},
                "features":     [],
            })
            ur = cas.get("upstream_regulatory") or {}
            if ur.get("start") is not None:
                hierarchical_annotations.append({
                    "name":          "Upstream Regulatory",
                    "start":         ur.get("start"),
                    "end":           ur.get("end"),
                    "direction":     ur.get("strand", 1),
                    "color":         _module_color("upstream_regulatory_module"),
                    "layer":         "module",
                    "module_type":   "upstream_regulatory_module",
                    "module_id":     f"{cassette_id}_ur",
                    "parent_module": cassette_id,
                    "source":        "mammalian_pol2_detector",
                    "metadata":      {"primary_promoter": ur.get("primary_promoter"), "components": ur.get("components", [])},
                    "features":      ur.get("components", []),
                })
            cds = cas.get("cds_module") or {}
            if cds.get("start") is not None:
                hierarchical_annotations.append({
                    "name":          f"CDS Module ({cds.get('aa_length', '?')} aa)",
                    "start":         cds.get("start"),
                    "end":           cds.get("end"),
                    "direction":     cds.get("strand", 1),
                    "color":         _module_color("cds_module"),
                    "layer":         "module",
                    "module_type":   "cds_module",
                    "module_id":     f"{cassette_id}_cds",
                    "parent_module": cassette_id,
                    "source":        "mammalian_pol2_detector",
                    "metadata":      {k: cds.get(k) for k in ("kozak_start", "kozak_strength", "aa_length", "stop_codon", "exons", "introns")},
                    "features":      [],
                })
            dr = cas.get("downstream_regulatory") or {}
            if dr.get("start") is not None:
                hierarchical_annotations.append({
                    "name":          "Downstream Regulatory",
                    "start":         dr.get("start"),
                    "end":           dr.get("end"),
                    "direction":     dr.get("strand", 1),
                    "color":         _module_color("downstream_regulatory_module"),
                    "layer":         "module",
                    "module_type":   "downstream_regulatory_module",
                    "module_id":     f"{cassette_id}_dr",
                    "parent_module": cassette_id,
                    "source":        "mammalian_pol2_detector",
                    "metadata":      {"polya_signal": dr.get("polya_signal"), "components": dr.get("components", [])},
                    "features":      dr.get("components", []),
                })

        # Step 4: Convert CDS submodules to hierarchical annotations
        print(f"[LLM] Step 4: Converting {len(cds_submodules_list)} CDS submodules to annotations...")
        
        cds_hierarchical = []
        color_map = {
            "protein_module": "#9C27B0",
            "nls_module": "#FF5722",
            "tag_module": "#00BCD4",
            "linker_module": "#FF9800",
            "flexible_linker_module": "#FF9800",
            "gap_module": "#9E9E9E",
        }

        # Filter ORFs: drop those that have NO same-direction overlapping
        # plannotate feature, OR for which same-direction gap_module
        # submodules cover >50% of the ORF's range. Survivors are replaced
        # with a 'translation' annotation carrying per-AA metadata.
        _CODON_TABLE = {
            'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
            'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
            'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
            'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
            'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
            'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
            'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
        }
        def _rev_comp(d):
            comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
            return ''.join(comp.get(b, 'N') for b in reversed(d))
        def _translate(nt_seq, strand):
            seq = nt_seq if strand == 1 else _rev_comp(nt_seq)
            aas = []
            for i in range(0, len(seq) - 2, 3):
                aas.append(_CODON_TABLE.get(seq[i:i+3].upper(), 'X'))
            return ''.join(aas)
        def _overlaps(a_start, a_end, b_start, b_end):
            return a_start < b_end and b_start < a_end

        surviving_orfs = []
        for orf in orf_modules:
            o_start, o_end = orf["start"], orf["end"]
            o_strand = orf.get("strand", 1)
            o_len = max(o_end - o_start, 1)

            # (a) same-direction overlapping plannotate feature?
            has_same_dir_feat = any(
                _overlaps(o_start, o_end, int(a.get("start", 0)), int(a.get("end", 0)))
                and int(a.get("direction", a.get("strand", 1)) or 1) == o_strand
                and (a.get("type", "") or "").lower() != "source"
                for a in plannotate_annotations
            )

            # (b) gap_module coverage from cds_submodules_list, same strand,
            # within the ORF range.
            gap_bp = 0
            for sub in cds_submodules_list:
                if sub.get("module_type") != "gap_module":
                    continue
                if sub.get("strand", 1) != o_strand:
                    continue
                g_start = max(int(sub.get("start", 0)), o_start)
                g_end = min(int(sub.get("end", 0)), o_end)
                if g_end > g_start:
                    gap_bp += g_end - g_start
            gap_pct = gap_bp / o_len

            if not has_same_dir_feat or gap_pct > 0.50:
                continue
            surviving_orfs.append(orf)

        print(f"[LLM] CDS ORF filter: {len(surviving_orfs)}/{len(orf_modules)} ORFs survived "
              f"(dropped {len(orf_modules) - len(surviving_orfs)})")

        # Replace surviving CDS ORFs with translation annotations.
        for orf in surviving_orfs:
            o_start, o_end = orf["start"], orf["end"]
            o_strand = orf.get("strand", 1)
            nt_seq = sequence[o_start:o_end]
            aa_seq = _translate(nt_seq, o_strand)
            # strip the trailing stop codon for display
            if aa_seq.endswith('*'):
                aa_seq = aa_seq[:-1]

            # Build feature_regions: for each plannotate annotation overlapping
            # this ORF in same direction, project its nucleotide range onto
            # 1-based AA coordinates within the ORF.
            feature_regions = []
            for a in plannotate_annotations:
                if int(a.get("direction", a.get("strand", 1)) or 1) != o_strand:
                    continue
                a_start = int(a.get("start", 0))
                a_end = int(a.get("end", 0))
                if not _overlaps(o_start, o_end, a_start, a_end):
                    continue
                if (a.get("type", "") or "").lower() == "source":
                    continue
                # nt offset within ORF, projected to aa index (1-based)
                if o_strand == 1:
                    nt_off_start = max(0, a_start - o_start)
                    nt_off_end = min(o_end - o_start, a_end - o_start)
                else:
                    nt_off_start = max(0, o_end - a_end)
                    nt_off_end = min(o_end - o_start, o_end - a_start)
                aa_start_idx = nt_off_start // 3 + 1
                aa_end_idx = (nt_off_end + 2) // 3
                if aa_end_idx < aa_start_idx:
                    continue
                feature_regions.append({
                    "name": a.get("name", "?"),
                    "aa_start": aa_start_idx,
                    "aa_end": aa_end_idx,
                    "feature_type": a.get("type", a.get("feature_type", "misc_feature")),
                })

            cds_hierarchical.append({
                "name": f"Translation ({len(aa_seq)} aa)",
                "start": o_start,
                "end": o_end,
                "direction": o_strand,
                "color": "#673AB7",
                "layer": "translation",
                "module_type": "translation",
                "source": "orf_detection",
                "metadata": {
                    "aa_length": len(aa_seq),
                    "aa_sequence": aa_seq,
                    "feature_regions": feature_regions,
                    "orf_detected": True,
                },
                "payload_id": None,
            })
        
        for sub in cds_submodules_list:
            module_type = sub.get("module_type", "cds_module")
            payload_id = sub.get("payload_id")
            md = sub.get("metadata", {}) or {}
            # Prefer the specific feature name (Cas9, PuroR-001, P2A, FLAG, ...)
            # over the generic module_type title. Falls back to module_type
            # for cds_submodules without a named feature (true gaps).
            display_name = (
                payload_id
                or md.get("detected_2a")
                or md.get("detected_nls")
                or md.get("detected_tag")
                or md.get("detected_linker")
                or module_type.replace("_", " ").title()
            )
            cds_hierarchical.append({
                "name": display_name,
                "start": sub["start"],
                "end": sub["end"],
                "direction": sub.get("strand", 1),
                "color": color_map.get(module_type, "#607D8B"),
                "layer": "module",
                "module_type": module_type,
                "source": "cds_submodule_parser",
                "metadata": md,
                "payload_id": payload_id,
                "module_family": "cds_submodule"
            })
        
        # Step 2.75: Scan for cloning features (restriction II/IIs, Gateway att,
        # PCR feasibility warnings). Fast, hidden-by-default in the viewer.
        step275_start = time.time()
        print("[LLM] Step 2.75: Scanning cloning features (restriction II/IIs, Gateway att, PCR warnings)...")
        from .cloning_feature_annotator import (
            scan_cloning_features,
            cloning_features_to_hierarchical,
        )
        cloning_scan = scan_cloning_features(sequence)
        cloning_hierarchical = cloning_features_to_hierarchical(cloning_scan.features)
        step275_time = time.time() - step275_start
        print(f"[LLM] Step 2.75 complete: {len(cloning_scan.features)} cloning features "
              f"({sum(1 for f in cloning_scan.features if f.feature_family == 'restriction_site_II')} Type II, "
              f"{sum(1 for f in cloning_scan.features if f.feature_family == 'restriction_site_IIs')} Type IIs, "
              f"{sum(1 for f in cloning_scan.features if f.feature_family == 'gateway_att')} Gateway, "
              f"{sum(1 for f in cloning_scan.features if f.feature_family == 'primer_design_warning')} PCR warnings); "
              f"{len(cloning_scan.non_cutters)} non-cutters")
        print(f"[TIMING] Step 2.75 (cloning features): {step275_time:.2f}s")

        # Merge LLM modules with CDS submodules and cloning features
        all_hierarchical = hierarchical_annotations + cds_hierarchical + cloning_hierarchical
        all_hierarchical = _filter_hierarchical_dups(plannotate_annotations, all_hierarchical)
        print(f"[LLM] Step 4 complete: {len(all_hierarchical)} total annotations ({len(hierarchical_annotations)} LLM + {len(cds_hierarchical)} CDS + {len(cloning_hierarchical)} cloning)")

        return {
            "ok": True,
            "method": "rule_based",
            "plannotate_annotations": plannotate_annotations,
            "annotations": plannotate_annotations,
            "cds_submodules": cds_submodules_list,
            "modules": modules,
            "interactions": interactions,
            "hierarchical_annotations": all_hierarchical,
            "cloning_features": cloning_scan.to_dict(),
            "summary": {
                "plannotate_feature_count": len(plannotate_annotations),
                "cds_submodule_count": len(cds_submodules_list),
                "module_count": len(modules),
                "rule_based_module_count": len(rule_based_modules),
                "pol2_cassette_count": len(pol2_cassettes),
                "interaction_count": len(interactions),
                "cloning_feature_count": len(cloning_scan.features),
                "non_cutter_count": len(cloning_scan.non_cutters),
                "source_endpoint": "annotate_sequence_llm",
            },
        }

    except Exception as e:
        import traceback
        import sys
        exc_type, exc_value, exc_tb = sys.exc_info()
        print(f"❌ ERROR in annotate_sequence_llm:")
        print(f"   Type: {exc_type.__name__}")
        print(f"   Message: {str(e)}")
        traceback.print_exc()
        return {
            "ok": False,
            "annotations": [],
            "error": str(e),
            "traceback": traceback.format_exc(),
        }
