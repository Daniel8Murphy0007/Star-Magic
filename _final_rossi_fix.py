#!/usr/bin/env python3
"""Read source .tex and convert to proper Markdown with robust parsing"""

import re

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.tex', 'r', encoding='utf-8') as f:
    content = f.read()

# Extract document body only
start = content.find(r'\begin{document}')
end = content.find(r'\end{document}')
if start >= 0 and end >= 0:
    content = content[start + len(r'\begin{document}'):end]

# Remove \maketitle, \begin{abstract}...\end{abstract}
content = content.replace(r'\maketitle', '')
content = re.sub(r'\\begin\{abstract\}', '## Abstract', content)
content = re.sub(r'\\end\{abstract\}', '', content)

def clean_equation(eq):
    """Clean equation content"""
    eq = re.sub(r'\\label\{[^}]*\}', '', eq)
    eq = eq.replace('\\\\', '\n')  # line breaks in align
    eq = eq.replace('&', '\n')     # alignment in align
    eq = eq.strip()
    return eq

# Handle sections - properly handle nested braces
def replace_sections(text):
    r"""Replace \section, \subsection, etc., handling nested braces"""
    text = re.sub(r'\\section\{([^}]+)\}', r'## \1', text)
    # Manual replacement for subsection
    while r'\subsection{' in text:
        start_idx = text.find(r'\subsection{')
        if start_idx == -1:
            break
        brace_start = start_idx + len(r'\subsection')
        depth = 0
        end_idx = -1
        for i in range(brace_start, len(text)):
            if text[i] == '{':
                depth += 1
            elif text[i] == '}':
                depth -= 1
                if depth == 0:
                    end_idx = i
                    break
        if end_idx > 0:
            content_inside = text[brace_start+1:end_idx]
            text = text[:start_idx] + '### ' + content_inside + text[end_idx+1:]
    
    # subsubsection
    text = re.sub(r'\\subsubsection\{([^}]+)\}', r'#### \1', text)
    return text

content = replace_sections(content)

# Remove bibliography sections
content = re.sub(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', '', content, flags=re.DOTALL)

# Remove \cite{}, \eqref{}, \label{}, \ref{}
content = re.sub(r'~?\\cite\{[^}]*\}', '', content)
content = re.sub(r'~?\\eqref\{[^}]*\}', '', content)
content = re.sub(r'~?\\label\{[^}]*\}', '', content)
content = re.sub(r'~?\\ref\{[^}]*\}', '', content)

# Convert equation/align environments to $$...$$ blocks
content = re.sub(r'\\begin\{equation\}(.*?)\\end\{equation\}', lambda m: '$$' + clean_equation(m.group(1)) + '$$', content, flags=re.DOTALL)
content = re.sub(r'\\begin\{align\*?\}(.*?)\\end\{align\*?\}', lambda m: '$$' + clean_equation(m.group(1)) + '$$', content, flags=re.DOTALL)

# Convert itemize/enumerate
def convert_itemize(m):
    items = re.findall(r'\\item\s+(.*?)(?=\\item|$)', m.group(1), flags=re.DOTALL)
    return '\n'.join(f'- {item.strip()}' for item in items) + '\n'

def convert_enumerate(m):
    items = re.findall(r'\\item\s+(.*?)(?=\\item|$)', m.group(1), flags=re.DOTALL)
    return '\n'.join(f'{i}. {item.strip()}' for i, item in enumerate(items, 1)) + '\n'

content = re.sub(r'\\begin\{itemize\}(.*?)\\end\{itemize\}', convert_itemize, content, flags=re.DOTALL)
content = re.sub(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', convert_enumerate, content, flags=re.DOTALL)

# Text formatting
content = re.sub(r'\\textbf\{([^}]*)\}', r'**\1**', content)
content = re.sub(r'\\textit\{([^}]*)\}', r'*\1*', content)
content = re.sub(r'\\emph\{([^}]*)\}', r'*\1*', content)

# Special characters
content = content.replace(r'\ ', ' ')
content = content.replace('~', ' ')
content = content.replace(r'\textdegree', '°')
content = content.replace(r'---', '—')
content = content.replace(r'--', '–')

# Remove remaining LaTeX commands (preserve math mode content)
parts = re.split(r'(\$\$.*?\$\$|\$[^$]*\$)', content, flags=re.DOTALL)
result = []
for i, part in enumerate(parts):
    if i % 2 == 1:  # Math or equation environment
        result.append(part)
    else:  # Text
        # Remove LaTeX commands
        part = re.sub(r'\\[a-zA-Z]+\{[^}]*\}', '', part)  # \cmd{...}
        part = re.sub(r'\\[a-zA-Z]+', '', part)  # \cmd
        result.append(part)
content = ''.join(result)

# Remove comments
content = re.sub(r'^%.*$', '', content, flags=re.MULTILINE)

# Clean up whitespace
content = re.sub(r'\n\n\n+', '\n\n', content)
lines = [line.rstrip() for line in content.split('\n')]
content = '\n'.join(lines).strip()

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.write(content)

print('[FINAL FIX] Robust Markdown conversion complete')
