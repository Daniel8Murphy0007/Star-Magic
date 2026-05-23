#!/usr/bin/env python3
"""Final smartest: protect ALL math before ANY processing"""

import re

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.tex', 'r', encoding='utf-8') as f:
    content = f.read()

# Extract body only
m = re.search(r'\\begin\{document\}(.*?)\\end\{document\}', content, re.DOTALL)
if m:
    content = m.group(1)

# === STEP 1: PROTECT ALL MATH MODE BEFORE ANYTHING ELSE ===
protected = {}
counter = 0

# Pattern to match \(...\) escaped math
def protect_escaped_math(m):
    global counter
    key = f'__PROTECT_MATH_{counter:04d}__'
    protected[key] = m.group(0)
    counter += 1
    return key

# Pattern to match $...$ single dollar
def protect_single_math(m):
    global counter
    key = f'__PROTECT_MATH_{counter:04d}__'
    protected[key] = m.group(0)
    counter += 1
    return key

# Pattern to match $$...$$ double dollar
def protect_double_math(m):
    global counter
    key = f'__PROTECT_MATH_{counter:04d}__'
    protected[key] = m.group(0)
    counter += 1
    return key

# Protect in order: escaped, then double dollar, then single dollar
content = re.sub(r'\\\([^)]*\\\)', protect_escaped_math, content, flags=re.DOTALL)
content = re.sub(r'\$\$[^$]+\$\$', protect_double_math, content, flags=re.DOTALL)
content = re.sub(r'\$[^\$]+\$', protect_single_math, content)

# === STEP 2: PROCESS LATEX ===

# Remove \maketitle
content = content.replace('\\maketitle', '')

# Remove bibliography
content = re.sub(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', '', content, flags=re.DOTALL)

# Remove citations/refs/labels
content = re.sub(r'~?\\cite\{[^}]*\}', '', content)
content = re.sub(r'~?\\eqref\{[^}]*\}', '', content)
content = re.sub(r'~?\\label\{[^}]*\}', '', content)
content = re.sub(r'~?\\ref\{[^}]*\}', '', content)

# === CONVERT SECTIONS ===
content = re.sub(r'\\begin\{abstract\}', '## Abstract', content)
content = re.sub(r'\\end\{abstract\}', '', content)
content = re.sub(r'\\section\{', '## {', content)
content = re.sub(r'\\subsection\{', '### {', content)
content = re.sub(r'\\subsubsection\{', '#### {', content)

# Fix heading braces
content = re.sub(r'^## \{([^}]+)\}', r'## \1', content, flags=re.MULTILINE)
content = re.sub(r'^### \{([^}]+)\}', r'### \1', content, flags=re.MULTILINE)
content = re.sub(r'^#### \{([^}]+)\}', r'#### \1', content, flags=re.MULTILINE)

# === CONVERT TABLES ===
def convert_tabular(m):
    tabular_content = m.group(1)
    # Remove LaTeX table directives
    tabular_content = tabular_content.replace('\\toprule', '')
    tabular_content = tabular_content.replace('\\midrule', '')
    tabular_content = tabular_content.replace('\\bottomrule', '')
    tabular_content = tabular_content.replace('\\hline', '')
    
    # Split rows by \\
    rows = tabular_content.split('\\\\')
    rows = [r.strip() for r in rows if r.strip()]
    
    if not rows:
        return ''
    
    # Process each row
    table_rows = []
    for row in rows:
        cells = [c.strip() for c in row.split('&')]
        cells = [c for c in cells if c]
        if cells:
            table_rows.append(cells)
    
    if not table_rows:
        return ''
    
    # Build markdown
    lines = []
    max_cols = max(len(r) for r in table_rows) if table_rows else 0
    
    for i, row in enumerate(table_rows):
        # Pad row
        while len(row) < max_cols:
            row.append('')
        lines.append('| ' + ' | '.join(row[:max_cols]) + ' |')
        # Add separator after first row
        if i == 0:
            lines.append('| ' + ' | '.join(['---'] * max_cols) + ' |')
    
    return '\n' + '\n'.join(lines) + '\n'

content = re.sub(r'\\begin\{center\}\s*\\begin\{tabular\}[^}]*\}(.*?)\\end\{tabular\}\s*\\end\{center\}', 
                 convert_tabular, content, flags=re.DOTALL)

# === CONVERT EQUATION ENVIRONMENTS ===
content = re.sub(r'\\begin\{equation\*?\}(.*?)\\end\{equation\*?\}', r'$$\1$$', content, flags=re.DOTALL)
content = re.sub(r'\\begin\{align\*?\}(.*?)\\end\{align\*?\}', r'$$\1$$', content, flags=re.DOTALL)

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

# === SPECIAL CHARS ===
content = content.replace('~', ' ')
content = content.replace(r'\ ', ' ')
content = content.replace('\\textdegree', '°')
content = content.replace('---', '—')
content = content.replace('--', '–')

# === REMOVE REMAINING LaTeX ===
# Remove \boxed, \center, etc.
content = re.sub(r'\\boxed\{([^}]*)\}', r'\1', content)
content = re.sub(r'\\begin\{center\}', '', content)
content = re.sub(r'\\end\{center\}', '', content)

# Remove other commands
content = re.sub(r'\\[a-zA-Z]+\{[^}]*\}', '', content)
content = re.sub(r'\\[a-zA-Z]+', '', content)

# Remove comment lines
content = re.sub(r'^%.*$', '', content, flags=re.MULTILINE)

# === RESTORE MATH BLOCKS ===
for key, math_block in protected.items():
    content = content.replace(key, math_block)

# === FINAL CLEANUP ===
content = re.sub(r'\n\n\n+', '\n\n', content)
content = '\n'.join(line.rstrip() for line in content.split('\n'))
content = content.strip()

# Add title
content = '# SCm - Rossi E-Cat Variants Unified\n\n' + content

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.write(content)

print('[FINAL SMARTEST] Math-protected conversion complete with title')
