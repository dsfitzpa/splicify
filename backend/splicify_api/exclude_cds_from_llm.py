#!/usr/bin/env python3
"""Update LLM prompt to exclude CDS-level modules"""

with open('llm_module_parser.py', 'r') as f:
    content = f.read()

# Add exclusion rule to the prompt
old_hierarchy = '''## Module Hierarchy (Outermost to Innermost)

1. **Payload Level**: lentiviral_payload, aav_payload, bacterial_backbone
2. **Expression Level**: pol2_expression_cassette, pol3_expression_cassette, guide_expression_cassette
3. **CDS Level**: cds_module
4. **Protein Level**: protein_submodule, nls_module, tag_module, linker_module'''

new_hierarchy = '''## Module Hierarchy (Outermost to Innermost)

**LLM should identify ONLY these levels:**
1. **Payload Level**: lentiviral_payload, aav_payload, bacterial_backbone
2. **Expression Level**: pol2_expression_cassette, pol3_expression_cassette, guide_expression_cassette
3. **Selection Cassette Level**: selection_cassette (bacterial/mammalian resistance markers)

**DO NOT identify these (handled by CDS submodule parser):**
- ❌ cds_module
- ❌ protein_module / protein_submodule
- ❌ nls_module
- ❌ tag_module
- ❌ linker_module
- ❌ gap_module

These are detected by ORF finder and CDS submodule parsing step.'''

content = content.replace(old_hierarchy, new_hierarchy)

# Add critical instruction
old_critical = '''**CRITICAL**: 
- Ensure every position 1..sequence_length is covered by at least one module
- Set strand based on the primary feature (promoter for expression cassettes, CDS for protein modules)
- Reverse complement features (like AmpR) should have strand=-1'''

new_critical = '''**CRITICAL**: 
- Ensure every position 1..sequence_length is covered by at least one module
- Set strand based on the primary feature (promoter for expression cassettes)
- Reverse complement features (like AmpR) should have strand=-1
- **DO NOT create cds_module, protein_module, nls_module, tag_module, or linker_module**
- Focus on: payloads, expression cassettes, selection cassettes, and backbone
- CDS-level modules are handled separately by ORF detection + submodule parsing'''

content = content.replace(old_critical, new_critical)

with open('llm_module_parser.py', 'w') as f:
    f.write(content)

print('✓ Updated LLM prompt to exclude CDS modules')
