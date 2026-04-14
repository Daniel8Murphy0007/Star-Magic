#!/usr/bin/env python3
"""Fix remaining backtick-in-math patterns."""
import glob, os, re

files = sorted(glob.glob('whitepapers/PAPER_*.md'))
fixed = 0

for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    original = c

    # Fix `word_` and `word` patterns back to \word_ and \word
    # inside $$ and $ blocks
    def fix_math_block(match):
        block = match.group(0)
        # Replace `word` -> \word (any alphabetic word in backticks)
        block = re.sub(r'`([a-zA-Z_]+)`', r'\\\1', block)
        return block

    # Fix $$ display math blocks
    c = re.sub(r'\$\$(.*?)\$\$', fix_math_block, c, flags=re.DOTALL)

    # Fix $ inline math (be careful not to match across line boundaries)
    # Only fix single-line inline math
    def fix_inline_math(match):
        block = match.group(0)
        if '`' in block:
            block = re.sub(r'`([a-zA-Z_]+)`', r'\\\1', block)
        return block

    c = re.sub(r'\$([^$\n]+)\$', fix_inline_math, c)

    if c != original:
        with open(f, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(c)
        fixed += 1

print(f'Fixed backtick-in-math in {fixed} papers')
