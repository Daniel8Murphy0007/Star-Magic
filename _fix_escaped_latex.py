#!/usr/bin/env python3
"""Fix broken LaTeX escaping patterns that break xelatex compilation."""
import glob, os, re

files = sorted(glob.glob('whitepapers/PAPER_*.md'))
fixed_count = 0

for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    original = c

    # Pattern 1: \\$\\ ... \\$\\$ (double-escaped display math)
    # These should be $$...$$
    c = re.sub(r'\\\$\\([^$]+)\\\$\\\$', r'$$\1$$', c)

    # Pattern 2: \\$\{ ... \$\$ (another variant)
    c = re.sub(r'\\\$\\\{([^}]+)\}\s*\\\$\\\$', r'$$\1$$', c)

    # Pattern 3: Remaining \$\{ ... } patterns
    c = re.sub(r'\\\$\\\{([^}]+)\}', r'$\1$', c)

    # Pattern 4: Standalone \\$\\ at start of line
    c = re.sub(r'^\\\$\\', '$$', c, flags=re.MULTILINE)

    # Pattern 5: \\$ at end of line
    c = re.sub(r'\\\$\\\$$', '$$', c, flags=re.MULTILINE)

    if c != original:
        with open(f, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(c)
        fixed_count += 1
        print(f'Fixed: {os.path.basename(f)}')

print(f'\nTotal: {fixed_count} files fixed')
