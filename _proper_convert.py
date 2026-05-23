#!/usr/bin/env python3
"""Proper LaTeX to Markdown conversion with table support"""

import re

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.tex', 'r', encoding='utf-8') as f:
    content = f.read()

# Extract body
content = re.sub(r'^.*?\\begin\{document\}', '', content, flags=re.DOTALL)
content = re.sub(r'\\end\{document\}.*$', '', content, flags=re.DOTALL)

# === REMOVE BIBLIOGRAPHY ===
content = re.sub(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', '', content, flags=re.DOTALL)

# === REMOVE CITATIONS/REFS ===
content = re.sub(r'~?\\cite\{[^}]*\}', '', content)
content = re.sub(r'~?\\eqref\{[^}]*\}', '', content)
content = re.sub(r'~?\\label\{[^}]*\}', '', content)
content = re.sub(r'~?\\ref\{[^}]*\}', '', content)

# === CONVERT TABULAR TABLES TO MARKDOWN ===
def convert_table(m):
    """Convert LaTeX tabular to markdown table"""
    table_content = m.group(1)
    
    # Split by \\ to get rows
    rows = table_content.split(r'\\')
    rows = [row.strip() for row in rows if row.strip()]
    
    if not rows:
        return ''
    
    # Remove \toprule, \midrule, \bottomrule, \hline
    rows = [re.sub(r'\\(toprule|midrule|bottomrule|hline)', '', row).strip() for row in rows]
    rows = [row for row in rows if row]
    
    if not rows:
        return ''
    
    # Split cells by &
    table_rows = []
    for row in rows:
        cells = [cell.strip() for cell in row.split('&')]
        cells = [cell for cell in cells if cell]  # Remove empty cells
        if cells:
            table_rows.append(cells)
    
    if not table_rows:
        return ''
    
    # Build markdown table
    # Header is first row
    header = table_rows[0]
    body = table_rows[1:] if len(table_rows) > 1 else []
    
    # Create markdown
    md_lines = []
    md_lines.append('| ' + ' | '.join(header) + ' |')
    md_lines.append('| ' + ' | '.join(['---'] * len(header)) + ' |')
    for row in body:
        # Pad row to match header length
        while len(row) < len(header):
            row.append('')
        md_lines.append('| ' + ' | '.join(row[:len(header)]) + ' |')
    
    return '\n' + '\n'.join(md_lines) + '\n'

# Convert tables
content = re.sub(r'\\begin\{center\}\s*\\begin\{tabular\}[^}]*\}(.*?)\\end\{tabular\}\s*\\end\{center\}', 
                 convert_table, content, flags=re.DOTALL)

# === CONVERT SECTIONS ===
content = content.replace(r'\maketitle', '')
content = re.sub(r'\\begin\{abstract\}', '## Abstract', content)
content = re.sub(r'\\end\{abstract\}', '', content)
content = re.sub(r'\\section\{([^}]+)\}', r'## \1', content)
content = re.sub(r'\\subsection\{([^}]+)\}', r'### \1', content)
content = re.sub(r'\\subsubsection\{([^}]+)\}', r'#### \1', content)

# === CONVERT EQUATIONS ===
content = re.sub(r'\\begin\{equation\}', '$$', content)
content = re.sub(r'\\end\{equation\}', '$$', content)
content = re.sub(r'\\begin\{align\*?\}', '$$', content)
content = re.sub(r'\\end\{align\*?\}', '$$', content)

# Within $$...$$, clean up
parts = re.split(r'(\$\$[^$]*\$\$)', content, flags=re.DOTALL)
for i in range(1, len(parts), 2):
    if parts[i].startswith('$$'):
        math = parts[i][2:-2]
        math = math.replace(' &', '')
        math = math.replace('&', '')
        math = math.replace('\\\\', '\n')
        math = re.sub(r'\\boxed\{([^}]*)\}', r'\1', math)
        math = math.strip()
        parts[i] = '$$' + math + '$$'
content = ''.join(parts)

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
content = content.replace(r'\textdegree', '°')
content = content.replace(r'---', '—')
content = content.replace(r'--', '–')
content = content.replace(r'\$', '$')  # Escaped $

# === REMOVE REMAINING LaTeX ===
content = re.sub(r'\\[a-zA-Z]+\{[^}]*\}', '', content)
content = re.sub(r'\\[a-zA-Z]+', '', content)

# === REMOVE COMMENTS & DIRECTIVES ===
content = re.sub(r'\\usepackage[^}]*\}', '', content)
content = re.sub(r'\\geometry[^}]*\}', '', content)
content = re.sub(r'\\hypersetup[^}]*\}', '', content)
content = re.sub(r'\\title[^}]*\}', '', content)
content = re.sub(r'\\author[^}]*\}', '', content)
content = re.sub(r'\\date[^}]*\}', '', content)
content = re.sub(r'^%.*$', '', content, flags=re.MULTILINE)

# === CLEANUP ===
content = re.sub(r'\n\n\n+', '\n\n', content)
content = '\n'.join(line.rstrip() for line in content.split('\n'))
content = content.strip()

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.write(content)

print('[PROPER CONVERT] LaTeX to Markdown complete with table support')
