#!/usr/bin/env python3
"""Fresh conversion of Rossi .tex to clean .md"""

import re

# Read original .tex file
with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.tex', 'r', encoding='utf-8') as f:
    content = f.read()

# Convert LaTeX structures to Markdown
# Abstract
content = re.sub(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', 
                 lambda m: '## Abstract\n\n' + m.group(1).strip() + '\n\n', 
                 content, flags=re.DOTALL)

# Equation environments - convert to $$...$$
content = re.sub(r'\\begin\{equation\*?\}(.*?)\\end\{equation\*?\}',
                 lambda m: '\n$$' + m.group(1).strip() + '$$\n',
                 content, flags=re.DOTALL)

# Align environments - convert to $$...$$
content = re.sub(r'\\begin\{align\*?\}(.*?)\\end\{align\*?\}',
                 lambda m: '\n$$\n' + re.sub(r'\\\\', '\n', m.group(1).strip()) + '\n$$\n',
                 content, flags=re.DOTALL)

# Remove \label{} and \eqref{}
content = re.sub(r'\\label\{[^}]*\}', '', content)
content = re.sub(r'\\eqref\{[^}]*\}', '', content)
content = re.sub(r'\\ref\{[^}]*\}', '', content)

# Remove bibliography
content = re.sub(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', '', content, flags=re.DOTALL)

# Remove \cite{}
content = re.sub(r'\\cite\{[^}]*\}', '', content)

# Remove \citep{} and \citet{}
content = re.sub(r'\\citet?\{[^}]*\}', '', content)

# Section headers - already start with #, but remove any LaTeX
content = re.sub(r'\\section\{(.*?)\}', r'## \1\n', content)
content = re.sub(r'\\subsection\{(.*?)\}', r'### \1\n', content)
content = re.sub(r'\\subsubsection\{(.*?)\}', r'#### \1\n', content)

# Abstract fix if not caught
content = re.sub(r'\\begin\{abstract\}', '', content)
content = re.sub(r'\\end\{abstract\}', '', content)

# Remove all other LaTeX commands at start of lines
content = re.sub(r'^\\[a-zA-Z]+', '', content, flags=re.MULTILINE)

# Remove }{, {}, etc. from formatting
content = re.sub(r'\{\}', '', content)

# Clean up multiple spaces
content = re.sub(r'  +', ' ', content)

# Clean up multiple newlines
content = re.sub(r'\n\n\n+', '\n\n', content)

# Write clean .md file
with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.write(content)

print('[CONVERTED] Fresh .md file from .tex')
