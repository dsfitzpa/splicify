#!/usr/bin/env python3
"""Update LLM prompt with strand detection rules"""

with open('llm_module_parser.py', 'r') as f:
    content = f.read()

# Find the MODULE_IDENTIFICATION_PROMPT
import re

# Add strand rule to the prompt after the Protein Expression section
old_rules = '''**Protein Expression**:
- 2A peptides (P2A/T2A/E2A/F2A) produce SEPARATE proteins via ribosome skipping
- IRES produces SEPARATE proteins via internal ribosome entry
- No separator = ONE fusion protein

**Promoter Classes**:'''

new_rules = '''**Protein Expression**:
- 2A peptides (P2A/T2A/E2A/F2A) produce SEPARATE proteins via ribosome skipping
- IRES produces SEPARATE proteins via internal ribosome entry
- No separator = ONE fusion protein

**Strand Direction**:
- Use the strand of the PRIMARY feature in each module
- Forward strand: +1, Reverse strand: -1
- Promoter strand determines cassette strand
- CDS strand determines protein module strand
- If features conflict, use the promoter/CDS strand (most important)

**Promoter Classes**:'''

content = content.replace(old_rules, new_rules)

# Also update the output format example to show strand usage
old_example = '''      "start": 200,
      "end": 3800,
      "strand": 1,'''

new_example = '''      "start": 200,
      "end": 3800,
      "strand": 1,  // Use promoter strand (CMV is forward)'''

content = content.replace(old_example, new_example)

# Add a reminder in the CRITICAL section
old_critical = '''**CRITICAL**: Ensure every position 1..sequence_length is covered by at least one module.'''

new_critical = '''**CRITICAL**: 
- Ensure every position 1..sequence_length is covered by at least one module
- Set strand based on the primary feature (promoter for expression cassettes, CDS for protein modules)
- Reverse complement features (like AmpR) should have strand=-1'''

content = content.replace(old_critical, new_critical)

with open('llm_module_parser.py', 'w') as f:
    f.write(content)

print('✓ Updated LLM prompt with strand detection rules')
