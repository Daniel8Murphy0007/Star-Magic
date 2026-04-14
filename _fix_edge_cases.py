#!/usr/bin/env python3
"""Fix remaining 34 edge cases with backtick/text in LaTeX math blocks."""
import glob, os, re

files = sorted(glob.glob('whitepapers/PAPER_*.md'))
fixed = 0

for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    original = c

    # Fix 1: `word_digits` and `Word_digits` in math → \word_digits
    c = re.sub(r'(?<=\$)`([A-Za-z_]\w*)`(?=[^$]*\$)', r'\\\1', c)

    # Fix 2: $$ blocks containing Markdown headers/tables → unwrap to prose
    def unwrap_bad_math(match):
        block = match.group(1)
        # If block contains Markdown headers ## or tables | or long prose
        if ('###' in block or
            block.count('|') > 4 or
            block.count('\\\\') > 5 and '---' in block):
            # This is probably a code-block that was wrongly converted to $$
            # Convert back to text/code
            lines = block.strip().split('\n')
            cleaned = []
            for line in lines:
                line = line.replace('\\\\', '')
                line = re.sub(r'^\s*& ?', '', line)
                cleaned.append(line)
            return '\n'.join(cleaned)
        return match.group(0)

    c = re.sub(r'\$\$\n\\begin\{aligned\}\n(.*?)\n\\end\{aligned\}\n\$\$',
               unwrap_bad_math, c, flags=re.DOTALL)

    # Fix 3: Remaining backtick-number patterns: `times100` -> \times 100
    c = re.sub(r'`times(\d+)`', r'\\times \1', c)
    c = re.sub(r'`right\)', r'\\right)', c)
    c = re.sub(r'`left\(', r'\\left(', c)

    if c != original:
        with open(f, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(c)
        fixed += 1

print(f'Fixed {fixed} papers')
