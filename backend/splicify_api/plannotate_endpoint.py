"""
pLannotate API Endpoint - Updated for GitHub Installation
"""

from flask import Blueprint, request, jsonify
import tempfile
import os
import sys
from pathlib import Path

# Create blueprint
plannotate_bp = Blueprint('plannotate', __name__)

# Try to import pLannotate
PLANNOTATE_AVAILABLE = False
PLANNOTATE_ERROR = None

# Method 1: Try pip-installed version
try:
    from plannotate import annotate
    PLANNOTATE_AVAILABLE = True
    print("✅ pLannotate loaded from pip installation")
except ImportError as e1:
    # Method 2: Try to import from git clone location
    if plannotate_path.exists():
        sys.path.insert(0, str(plannotate_path))
        try:
            from plannotate import annotate
            PLANNOTATE_AVAILABLE = True
            print(f"✅ pLannotate loaded from {plannotate_path}")
        except ImportError as e2:
            PLANNOTATE_ERROR = f"Found pLannotate repo but import failed: {e2}"
            print(f"❌ {PLANNOTATE_ERROR}")
    else:
        PLANNOTATE_ERROR = f"pLannotate not found. Install with: cd /root/python-libraries && git clone https://github.com/mmcguffi/pLannotate.git && cd pLannotate && pip install -e . --break-system-packages"
        print(f"❌ {PLANNOTATE_ERROR}")


@plannotate_bp.route('/annotate_genbank', methods=['POST'])
def annotate_genbank():
    """
    Annotate a GenBank file using pLannotate
    
    Request JSON:
    {
        "gb_text": "LOCUS ...",
        "session_id": "abc123",
        "options": {
            "linear": false,
            "batch_size": 100,
            "file_name": "output"
        }
    }
    
    Response JSON:
    {
        "ok": true,
        "annotated_gb": "LOCUS ... [with annotations]",
        "annotation_count": 15,
        "annotations": [
            {"type": "promoter", "name": "T7 promoter", ...},
            ...
        ]
    }
    """
    
    if not PLANNOTATE_AVAILABLE:
        return jsonify({
            "ok": False,
            "error": "pLannotate not available",
            "details": PLANNOTATE_ERROR
        }), 500
    
    try:
        data = request.get_json()
        
        gb_text = data.get('gb_text')
        session_id = data.get('session_id', 'unknown')
        options = data.get('options', {})
        
        if not gb_text:
            return jsonify({
                "ok": False,
                "error": "Missing gb_text parameter"
            }), 400
        
        # Write to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gb', delete=False) as tmp_in:
            tmp_in.write(gb_text)
            tmp_in_path = tmp_in.name
        
        # Create output directory
        output_dir = tempfile.mkdtemp()
        
        try:
            # Extract options
            linear = options.get('linear', False)
            batch_size = options.get('batch_size', 100)
            file_name = options.get('file_name', f'gibson_{session_id}')
            
            print(f"🔬 Annotating {file_name} (linear={linear})...")
            
            # Call pLannotate
            # Note: Function signature may vary by version
            # Check: https://github.com/mmcguffi/pLannotate
            try:
                # Try newer API (with all params)
                annotate(
                    file_path=tmp_in_path,
                    output_path=output_dir,
                    linear=linear,
                    batch_size=batch_size,
                    html=False,
                    csv=False,
                    file_name=file_name
                )
            except TypeError:
                # Try older API (fewer params)
                annotate(
                    file_path=tmp_in_path,
                    output_path=output_dir,
                    linear=linear
                )
            
            # Find the output file (pLannotate may use .gbk or .gb)
            output_gb_path = None
            for ext in ['.gbk', '.gb', '.genbank']:
                candidate = Path(output_dir) / f"{file_name}{ext}"
                if candidate.exists():
                    output_gb_path = candidate
                    break
            
            # If not found with file_name, look for any .gbk file
            if not output_gb_path:
                gbk_files = list(Path(output_dir).glob("*.gbk"))
                if gbk_files:
                    output_gb_path = gbk_files[0]
            
            if not output_gb_path or not output_gb_path.exists():
                # List what files were created
                files = list(Path(output_dir).glob("*"))
                return jsonify({
                    "ok": False,
                    "error": "pLannotate did not generate output file",
                    "output_dir": output_dir,
                    "files_created": [str(f) for f in files]
                }), 500
            
            # Read the annotated GenBank file
            with open(output_gb_path, 'r') as f:
                annotated_gb = f.read()
            
            # Parse annotations
            annotation_count = 0
            annotations = []
            
            in_features = False
            current_feature = None
            
            for line in annotated_gb.split('\n'):
                if line.startswith('FEATURES'):
                    in_features = True
                    continue
                elif line.startswith('ORIGIN') or line.startswith('//'):
                    in_features = False
                    if current_feature:
                        annotations.append(current_feature)
                    continue
                
                if not in_features:
                    continue
                
                # Feature line: "     CDS             123..456"
                if line.startswith('     ') and not line.strip().startswith('/'):
                    # Save previous feature
                    if current_feature:
                        annotations.append(current_feature)
                        annotation_count += 1
                    
                    # Parse new feature
                    parts = line.strip().split(None, 1)
                    if len(parts) >= 2:
                        feature_type = parts[0]
                        location = parts[1]
                        current_feature = {
                            "type": feature_type,
                            "location": location,
                            "qualifiers": {}
                        }
                    else:
                        current_feature = None
                
                # Qualifier line: "                     /label="T7 promoter""
                elif line.strip().startswith('/') and current_feature:
                    qualifier_line = line.strip()[1:]  # Remove /
                    if '=' in qualifier_line:
                        key, value = qualifier_line.split('=', 1)
                        # Remove quotes
                        value = value.strip('"')
                        current_feature['qualifiers'][key] = value
                        
                        # Store name from common qualifiers
                        if key in ['label', 'note', 'product', 'gene']:
                            current_feature['name'] = value
            
            # Add last feature
            if current_feature:
                annotations.append(current_feature)
                annotation_count += 1
            
            print(f"✅ Added {annotation_count} annotations")
            
            return jsonify({
                "ok": True,
                "annotated_gb": annotated_gb,
                "annotation_count": annotation_count,
                "annotations": annotations[:20]  # Return first 20
            })
            
        finally:
            # Clean up temp files
            try:
                os.unlink(tmp_in_path)
                import shutil
                shutil.rmtree(output_dir)
            except Exception as cleanup_error:
                print(f"Warning: Cleanup failed: {cleanup_error}")
    
    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        print(f"❌ Annotation error: {error_details}")
        
        return jsonify({
            "ok": False,
            "error": str(e),
            "traceback": error_details
        }), 500


@plannotate_bp.route('/health', methods=['GET'])
def health_check():
    """Check if pLannotate is available"""
    return jsonify({
        "ok": True,
        "plannotate_available": PLANNOTATE_AVAILABLE,
        "error": PLANNOTATE_ERROR if not PLANNOTATE_AVAILABLE else None,
        "version": "1.0"
    })


# Register blueprint in main app.py:
# from plannotate_endpoint import plannotate_bp
# app.register_blueprint(plannotate_bp, url_prefix='/plannotate')
