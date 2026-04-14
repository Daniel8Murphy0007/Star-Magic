#!/usr/bin/env python3
"""Replace Unicode characters that break xelatex with safe alternatives."""
import glob, os, re

files = sorted(glob.glob('whitepapers/PAPER_*.md'))
fixed = 0

# Unicode replacements (char -> safe alternative)
REPLACEMENTS = {
    '\u2713': 'PASS',    # ✓ check mark -> PASS
    '\u2714': 'PASS',    # ✔ heavy check mark -> PASS  
    '\u2717': 'FAIL',    # ✗ ballot x -> FAIL
    '\u2718': 'FAIL',    # ✘ heavy ballot x -> FAIL
    '\u207B': '-',        # ⁻ superscript minus -> -
    '\u2070': '0',        # ⁰ superscript zero
    '\u00B9': '1',        # ¹ superscript one
    '\u00B2': '2',        # ² superscript two
    '\u00B3': '3',        # ³ superscript three
    '\u2074': '4',        # ⁴ superscript four
    '\u2075': '5',        # ⁵ superscript five
    '\u2076': '6',        # ⁶ superscript six
    '\u2077': '7',        # ⁷ superscript seven
    '\u2078': '8',        # ⁸ superscript eight
    '\u2079': '9',        # ⁹ superscript nine
}

for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    original = c

    # Replace "PASS ✓" with just "PASS"
    c = c.replace('PASS \u2713', 'PASS')
    c = c.replace('PASS \u2714', 'PASS')

    # Replace standalone check marks
    for old, new in REPLACEMENTS.items():
        if old in c:
            c = c.replace(old, new)

    # Replace Greek letters in non-math context with LaTeX equivalents
    # κ in prose -> leave as-is (handled by font), but in headers remove
    # Actually, for xelatex with the right font, Greek should be fine
    # The issue is specific to checkmarks and superscript chars

    if c != original:
        with open(f, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(c)
        fixed += 1

print(f'Replaced Unicode in {fixed} papers')
