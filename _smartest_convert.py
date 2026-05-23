#!/usr/bin/env python3
"""Smartest conversion: protect math mode, convert text mode only"""

import re

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.tex', 'r', encoding='utf-8') as f:
    content = f.read()

# Extract body only
m = re.search(r'\\begin\{document\}(.*?)\\end\{document\}', content, re.DOTALL)
if m:
    content = m.group(1)
else:
    content = ''

content = content.replace(r'\maketitle', '')

# === EXTRACT AND PROTECT MATH MODE ===
# Find all $...$ and $$...$$ blocks and replace with placeholders
math_blocks = {}
math_counter = 0

# Protect $$...$$ first (before single $)
def protect_double_math(m):
    global math_counter
    key = f'__MATH_DOUBLE_{math_counter}__'
    math_blocks[key] = m.group(0)
    math_counter += 1
    return key

def protect_single_math(m):
    global math_counter
    key = f'__MATH_SINGLE_{math_counter}__'
    math_blocks[key] = m.group(0)
    math_counter += 1
    return key

content = re.sub(r'\$\$[^$]*\$\$', protect_double_math, content, flags=re.DOTALL)
content = re.sub(r'\$[^$]*\$', protect_single_math, content)

# === NOW PROCESS REMAINING LaTeX ===

# Remove bibliography
content = re.sub(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', '', content, flags=re.DOTALL)

# Remove citations/refs
content = re.sub(r'~?\\cite\{[^}]*\}', '', content)
content = re.sub(r'~?\\eqref\{[^}]*\}', '', content)
content = re.sub(r'~?\\label\{[^}]*\}', '', content)
content = re.sub(r'~?\\ref\{[^}]*\}', '', content)

# === CONVERT SECTIONS/HEADINGS ===
content = re.sub(r'\\begin\{abstract\}', '## Abstract', content)
content = re.sub(r'\\end\{abstract\}', '', content)
content = re.sub(r'\\section\{([^}]+)\}', r'## \1', content)
content = re.sub(r'\\subsection\{([^}]+)\}', r'### \1', content)
content = re.sub(r'\\subsubsection\{([^}]+)\}', r'#### \1', content)

# === CONVERT TABLES ===
def convert_table(m):
    """Convert LaTeX tabular to markdown"""
    table = m.group(1)
    # Clean up table markup
    table = table.replace(r'\toprule', '')
    table = table.replace(r'\midrule', '')
    table = table.replace(r'\bottomrule', '')
    table = table.replace(r'\hline', '')
    
    # Split by \\ to get rows
    rows = table.split('\\\\')
    rows = [r.strip() for r in rows if r.strip()]
    
    if not rows:
        return ''
    
    # Split cells by &
    cells_per_row = []
    for row in rows:
        cells = [c.strip() for c in row.split('&')]
        cells = [c for c in cells if c]
        if cells:
            cells_per_row.append(cells)
    
    if not cells_per_row:
        return ''
    
    # Get max cols
    max_cols = max(len(row) for row in cells_per_row) if cells_per_row else 0
    
    # Build table
    lines = []
    for i, row in enumerate(cells_per_row):
        # Pad to max_cols
        row = row + [''] * (max_cols - len(row))
        lines.append('| ' + ' | '.join(row[:max_cols]) + ' |')
        # Add separator after header
        if i == 0:
            lines.append('| ' + ' | '.join(['---'] * max_cols) + ' |')
    
    return '\n' + '\n'.join(lines) + '\n'

content = re.sub(r'\\begin\{center\}\s*\\begin\{tabular\}[^}]*\}(.*?)\\end\{tabular\}\s*\\end\{center\}', 
                 convert_table, content, flags=re.DOTALL)

# === CONVERT EQUATION ENVIRONMENTS ===
content = re.sub(r'\\begin\{equation\}(.*?)\\end\{equation\}', r'$$\1$$', content, flags=re.DOTALL)
content = re.sub(r'\\begin\{align\*?\}(.*?)\\end\{align\*?\}', r'$$\1$$', content, flags=re.DOTALL)

# Clean up equation content (remove & and \\)
def clean_equation(m):
    eq = m.group(1)
    eq = eq.replace(' &', '')
    eq = eq.replace('&', '')
    eq = eq.replace('\\\\', '\n')
    return '$$' + eq.strip() + '$$'

content = re.sub(r'\$\$([^$]*)\$\$', clean_equation, content, flags=re.DOTALL)

# === LISTS ===
content = re.sub(r'\\begin\{itemize\}', '', content)
content = re.sub(r'\\end\{itemize\}', '', content)
content = re.sub(r'\\item\s+', '- ', content)

content = re.sub(r'\\begin\{enumerate\}', '', content)
content = re.sub(r'\\end\{enumerate\}', '', content)

# === TEXT FORMATTING ===
content = re.sub(r'\\textbf\{([^}]*)\}', r'**\1**', content)
content = re.sub(r'\\textit\{([^}]*)\}', r'*\1*', content)
content = re.sub(r'\\emph\{([^}]*)\}', r'*\1*', content)
content = re.sub(r'\\text\{([^}]*)\}', r'\1', content)

# === SPECIAL CHARS ===
content = content.replace(r'\textdegree', '°')
content = content.replace(r'\ ', ' ')
content = content.replace('~', ' ')
content = content.replace(r'---', '—')
content = content.replace(r'--', '–')
content = content.replace(r'\$', '$')

# === REMOVE REMAINING LaTeX COMMANDS (except protected math) ===
content = re.sub(r'\\boxed\{([^}]*)\}', r'\1', content)
content = re.sub(r'\\[a-zA-Z]+\{[^}]*\}', '', content)
content = re.sub(r'\\[a-zA-Z]+', '', content)

# === REMOVE COMMENTS ===
content = re.sub(r'^%.*$', '', content, flags=re.MULTILINE)

# === RESTORE MATH BLOCKS ===
for key, math_block in math_blocks.items():
    content = content.replace(key, math_block)

# === FINAL CLEANUP ===
content = re.sub(r'\n\n\n+', '\n\n', content)
content = '\n'.join(line.rstrip() for line in content.split('\n'))
content = content.strip()

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.write(content)

print('[SMARTEST CONVERT] Math-mode protected conversion complete')
