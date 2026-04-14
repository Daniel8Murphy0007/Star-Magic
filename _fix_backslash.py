#!/usr/bin/env python3
"""Fix backslash-filename issues that break xelatex."""
import glob, os, re

files = sorted(glob.glob('whitepapers/PAPER_*.md'))
fixed = 0
for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    original = c

    # Fix: \alidate -> validate (was \validate with \v consumed as LaTeX)
    c = c.replace('\\alidate', 'validate')

    # Fix bare backslash before filenames of any extension
    # \word.py\ -> `word.py`, \source27.cpp\ -> `source27.cpp`
    c = re.sub(r'\\([a-zA-Z_][a-zA-Z0-9_]*\.[a-z]+)\\', r'`\1`', c)
    
    # Also fix \source without extension (bare LaTeX-breaking backslash)
    c = re.sub(r'\\(source\d+)', r'`\1`', c)
    
    # Fix \function_name\ patterns (function/method names)
    # \compute_Ug1_SOURCE4\ -> `compute_Ug1_SOURCE4`
    c = re.sub(r'\\([a-zA-Z_][a-zA-Z0-9_]*(?:\(\))?)\\', r'`\1`', c)

    if c != original:
        with open(f, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(c)
        fixed += 1

print(f'Fixed in {fixed} papers')
