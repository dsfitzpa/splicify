#!/usr/bin/env python3
"""Add CDS submodules to hierarchical annotations"""

with open('plannotate_router.py', 'r') as f:
    lines = f.readlines()

# Find the line with 'print(f"[LLM] Step 3 complete'
insert_line = None
for i, line in enumerate(lines):
    if 'print(f"[LLM] Step 3 complete:' in line:
        insert_line = i + 1
        break

if not insert_line:
    print('Could not find insertion point')
    exit(1)

# Insert code to convert CDS submodules to hierarchical annotations
new_code = '''
        # Step 4: Convert CDS submodules to hierarchical annotations
        print(f"[LLM] Step 4: Converting {len(cds_submodules_list)} CDS submodules to annotations...")
        
        cds_hierarchical = []
        color_map = {
            "protein_module": "#9C27B0",
            "nls_module": "#FF5722",
            "tag_module": "#00BCD4",
            "linker_module": "#FF9800",
            "gap_module": "#9E9E9E"
        }
        
        for sub in cds_submodules_list:
            module_type = sub.get("module_type", "cds_module")
            cds_hierarchical.append({
                "name": module_type.replace("_", " ").title(),
                "start": sub["start"],
                "end": sub["end"],
                "direction": sub.get("strand", 1),
                "color": color_map.get(module_type, "#607D8B"),
                "layer": "module",
                "module_type": module_type,
                "source": "cds_submodule_parser",
                "metadata": sub.get("metadata", {}),
                "payload_id": None,
                "module_family": "cds_submodule"
            })
        
        # Merge LLM modules with CDS submodules
        all_hierarchical = hierarchical_annotations + cds_hierarchical
        print(f"[LLM] Step 4 complete: {len(all_hierarchical)} total annotations ({len(hierarchical_annotations)} LLM + {len(cds_hierarchical)} CDS)")
'''

lines.insert(insert_line, new_code)

# Update the return statement to use all_hierarchical
for i, line in enumerate(lines):
    if '"hierarchical_annotations": hierarchical_annotations,' in line:
        lines[i] = line.replace('hierarchical_annotations', 'all_hierarchical')
        break

with open('plannotate_router.py', 'w') as f:
    f.writelines(lines)

print('✓ Added CDS submodule to hierarchical annotation conversion')
