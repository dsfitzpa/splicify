#!/usr/bin/env python3
"""
Update script to replace the LLM endpoint's Step 2 with ORF detection
"""

# Read the current file
with open('plannotate_router.py', 'r') as f:
    lines = f.readlines()

# Find the start and end of Step 2
start_line = None
end_line = None

for i, line in enumerate(lines):
    if '# Step 2: Run CDS submodule parsing' in line:
        start_line = i
    if start_line and '# Step 3: Use LLM to identify functional modules' in line:
        end_line = i
        break

if not start_line or not end_line:
    print('Could not find Step 2 section')
    exit(1)

print(f'Found Step 2: lines {start_line+1} to {end_line}')

# New Step 2 implementation
new_step2 = '''        # Step 2: Detect ORFs and parse CDS submodules
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
                        "start": submod.start,
                        "end": submod.end,
                        "strand": cds_strand,  # Preserve ORF strand
                        "metadata": submod.metadata
                    })
            except Exception as e:
                print(f"[WARN] CDS submodule resolution failed for ORF {cds_start}-{cds_end}: {e}")
        
        print(f"[LLM] Step 2 complete: {len(cds_submodules_list)} CDS submodules from {len(detected_orfs)} ORFs")

'''

# Replace lines
new_lines = lines[:start_line] + [new_step2] + lines[end_line:]

# Write back
with open('plannotate_router.py', 'w') as f:
    f.writelines(new_lines)

print('✓ Updated plannotate_router.py with ORF detection')
