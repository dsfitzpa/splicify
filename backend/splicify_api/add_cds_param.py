#!/usr/bin/env python3
"""Add cds_submodules parameter to annotate_with_llm_modules"""

# Update llm_module_parser.py
with open('llm_module_parser.py', 'r') as f:
    content = f.read()

# Update function signature
old_sig = '''async def annotate_with_llm_modules(
    sequence: str,
    plannotate_rows: List[Dict[str, Any]],
    circular: bool = True,
    api_key: Optional[str] = None
) -> Dict[str, Any]:'''

new_sig = '''async def annotate_with_llm_modules(
    sequence: str,
    plannotate_rows: List[Dict[str, Any]],
    circular: bool = True,
    cds_submodules: Optional[List[Dict[str, Any]]] = None,
    api_key: Optional[str] = None
) -> Dict[str, Any]:'''

content = content.replace(old_sig, new_sig)

# Update the parse_modules call to pass cds_submodules
old_parse = '''        modules_result = await parser.parse_modules(
            sequence=sequence,
            plannotate_features=plannotate_features,
            cds_submodules=None,
            circular=circular
        )'''

new_parse = '''        modules_result = await parser.parse_modules(
            sequence=sequence,
            plannotate_features=plannotate_features,
            cds_submodules=cds_submodules,
            circular=circular
        )'''

content = content.replace(old_parse, new_parse)

with open('llm_module_parser.py', 'w') as f:
    f.write(content)

print('✓ Added cds_submodules parameter to annotate_with_llm_modules')

# Update plannotate_router.py to pass cds_submodules
with open('plannotate_router.py', 'r') as f:
    content = f.read()

old_call = '''        llm_result = await annotate_with_llm_modules(
            sequence=sequence,
            plannotate_rows=rows,
            circular=circular
        )'''

new_call = '''        llm_result = await annotate_with_llm_modules(
            sequence=sequence,
            plannotate_rows=rows,
            circular=circular,
            cds_submodules=cds_submodules_list
        )'''

content = content.replace(old_call, new_call)

with open('plannotate_router.py', 'w') as f:
    f.write(content)

print('✓ Updated plannotate_router.py to pass cds_submodules')
