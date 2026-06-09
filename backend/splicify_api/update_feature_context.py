#!/usr/bin/env python3
"""Update feature context to include strand information"""

with open('llm_module_parser.py', 'r') as f:
    lines = f.readlines()

# Find and replace the entry line
for i, line in enumerate(lines):
    if 'entry = f"{f.get(\'name\')} ({start+1}..{end})"' in line:
        # Replace with strand-aware version
        lines[i] = '            strand = f.get("strand", 1)\n'
        lines.insert(i+1, '            strand_str = "+" if strand == 1 else "-"\n')
        lines.insert(i+2, '            entry = f"{f.get(\'name\')} ({start+1}..{end}, {strand_str})"\n')
        lines.insert(i+3, '\n')
        break

with open('llm_module_parser.py', 'w') as f:
    f.writelines(lines)

print('✓ Updated feature context to include strand')
