#!/usr/bin/env python3
"""Add dotenv loading to main.py"""

with open('main.py', 'r') as f:
    lines = f.readlines()

# Check if load_dotenv is already there
if any('load_dotenv' in line for line in lines):
    print('✓ load_dotenv already present')
else:
    # Add after 'from __future__ import annotations'
    for i, line in enumerate(lines):
        if 'from __future__ import annotations' in line:
            # Insert blank line and dotenv loading
            lines.insert(i+1, '\n')
            lines.insert(i+2, '# Load environment variables from .env file\n')
            lines.insert(i+3, 'from dotenv import load_dotenv\n')
            lines.insert(i+4, 'from pathlib import Path\n')
            lines.insert(i+5, 'env_path = Path(__file__).parent.parent / ".env"\n')
            lines.insert(i+6, 'load_dotenv(dotenv_path=env_path)\n')
            lines.insert(i+7, '\n')
            break
    
    with open('main.py', 'w') as f:
        f.writelines(lines)
    
    print('✓ Added load_dotenv to main.py')
