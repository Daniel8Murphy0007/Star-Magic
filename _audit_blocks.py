#!/usr/bin/env python3
"""Re-audit session blocks with proper quote handling."""
import re

t = open('CondensedPhysics4.py', encoding='utf-8').read()
blocks = re.findall(r'(_SESSION_\w+_CLASSES)\s*=\s*\[([^\]]*)\]', t)
for name, content in blocks:
    entries = re.findall(r'''["'](\w+)["']''', content)
    status = f'{len(entries)} entries' if entries else 'EMPTY'
    print(f'  {name}: {status}')
