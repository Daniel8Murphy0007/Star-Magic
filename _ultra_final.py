#!/usr/bin/env python3
"""Ultra-final: protect EVERYTHING including equation blocks"""

import re

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.tex', 'r', encoding='utf-8') as f:
    content = f.read()

# Extract body
m = re.search(r'\\begin\{document\}(.*?)\\end\{document\}', content, re.DOTALL)
if m:
    content = m.group(1)

# === PROTECT ALL MATH/EQUATION BLOCKS FIRST ===
protected = {}
counter = 0

def protect(m):
    global counter
    key = f'___PROT_{counter:05d}___'
    protected[key] = m.group(0)
    counter += 1
    return key

# Protect equation blocks
content = re.sub(r'\\begin\{equation\*?\}.*?\\end\{equation\*?\}', protect, content, flags=re.DOTALL)
# Protect align blocks
content = re.sub(r'\\begin\{align\*?\}.*?\\end\{align\*?\}', protect, content, flags=re.DOTALL)
# Protect $$...$$ blocks
content = re.sub(r'\$\$[^$]*\$\$', protect, content, flags=re.DOTALL)
# Protect $...$ blocks
content = re.sub(r'\$[^\$]*\$', protect, content)

# === NOW SAFE TO PROCESS LATEX ===

content = content.replace('\\maketitle', '')

# Remove bibliography
content = re.sub(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', '', content, flags=re.DOTALL)

# Remove citations
content = re.sub(r'~?\\(cite|eqref|label|ref)\{[^}]*\}', '', content)

# === CONVERT SECTIONS ===
content = re.sub(r'\\begin\{abstract\}', '## Abstract', content)
content = re.sub(r'\\end\{abstract\}', '', content)
content = re.sub(r'\\section\{([^}]*)\}', r'## \1', content)
content = re.sub(r'\\subsection\{([^}]*)\}', r'### \1', content)
content = re.sub(r'\\subsubsection\{([^}]*)\}', r'#### \1', content)

# === CONVERT TABLES ===
def convert_table(m):
    tab = m.group(1)
    tab = tab.replace('\\toprule', '').replace('\\midrule', '').replace('\\bottomrule', '').replace('\\hline', '')
    rows = [r.strip() for r in tab.split('\\\\') if r.strip()]
    if not rows:
        return ''
    table_rows = [[c.strip() for c in r.split('&') if c.strip()] for r in rows]
    if not table_rows:
        return ''
    max_cols = max(len(r) for r in table_rows)
    lines = []
    for i, row in enumerate(table_rows):
        while len(row) < max_cols:
            row.append('')
        lines.append('| ' + ' | '.join(row[:max_cols]) + ' |')
        if i == 0:
            lines.append('| ' + ' | '.join(['---'] * max_cols) + ' |')
    return '\n' + '\n'.join(lines) + '\n'

content = re.sub(r'\\begin\{center\}\s*\\begin\{tabular\}[^}]*\}(.*?)\\end\{tabular\}\s*\\end\{center\}', 
                 convert_table, content, flags=re.DOTALL)

# === LISTS ===
content = re.sub(r'\\begin\{(itemize|enumerate)\}', '', content)
content = re.sub(r'\\end\{(itemize|enumerate)\}', '', content)
content = re.sub(r'\\item\s+', '- ', content)

# === TEXT FORMAT ===
content = re.sub(r'\\textbf\{([^}]*)\}', r'**\1**', content)
content = re.sub(r'\\textit\{([^}]*)\}', r'*\1*', content)
content = re.sub(r'\\emph\{([^}]*)\}', r'*\1*', content)

# === SPECIAL CHARS ===
content = content.replace('~', ' ')
content = content.replace(r'\ ', ' ')
content = content.replace('\\textdegree', '°')
content = content.replace('---', '—')
content = content.replace('--', '–')

# === REMOVE OTHER LaTeX (safe now that math is protected) ===
content = re.sub(r'\\boxed\{([^}]*)\}', r'\1', content)
content = re.sub(r'\\(begin|end)\{(center|tabular)\}', '', content)
content = re.sub(r'\\[a-zA-Z]+\{[^}]*\}', '', content)
content = re.sub(r'\\[a-zA-Z]+\[[^\]]*\]', '', content)
content = re.sub(r'\\[a-zA-Z]+', '', content)

# Remove comments
content = re.sub(r'^%.*$', '', content, flags=re.MULTILINE)

# === RESTORE PROTECTED BLOCKS ===
for key, block in protected.items():
    content = content.replace(key, block)

# === CLEANUP ===
content = re.sub(r'\n\n\n+', '\n\n', content)
content = '\n'.join(line.rstrip() for line in content.split('\n'))
content = content.strip()

# Add title
content = '# Rossi E-Cat Variants Under SCm\n\n' + content

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.write(content)

print('[ULTRA FINAL] Conversion with full protection complete')
