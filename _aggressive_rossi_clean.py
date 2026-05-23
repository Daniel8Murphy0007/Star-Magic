#!/usr/bin/env python3
"""Aggressive cleaning to remove ALL LaTeX from Rossi .md"""

import re

# Read current .md
with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Completely reconstruct from .tex source
with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.tex', 'r', encoding='utf-8') as f:
    tex = f.read()

# Remove everything before \begin{document}
tex = re.sub(r'^.*?\\begin\{document\}', '', tex, flags=re.DOTALL)
# Remove everything after \end{document}
tex = re.sub(r'\\end\{document\}.*$', '', tex, flags=re.DOTALL)

# === STRUCTURAL CONVERSIONS ===

# Abstract - before any other sections
tex = re.sub(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', 
             r'\n## Abstract\n\n\1\n\n', 
             tex, flags=re.DOTALL)

# Sections
tex = re.sub(r'\\section\{(.*?)\}', r'\n## \1\n', tex)
tex = re.sub(r'\\subsection\{(.*?)\}', r'\n### \1\n', tex)
tex = re.sub(r'\\subsubsection\{(.*?)\}', r'\n#### \1\n', tex)

# === EQUATION CONVERSIONS ===

# equation* and equation
tex = re.sub(r'\\begin\{equation\*?\}(.*?)\\end\{equation\*?\}',
             lambda m: '\n$$\n' + m.group(1).strip() + '\n$$\n',
             tex, flags=re.DOTALL)

# align* and align
tex = re.sub(r'\\begin\{align\*?\}(.*?)\\end\{align\*?\}',
             lambda m: '\n$$\n' + re.sub(r'\\\\', '\n', m.group(1).strip()) + '\n$$\n',
             tex, flags=re.DOTALL)

# multline and similar
tex = re.sub(r'\\begin\{multline\*?\}(.*?)\\end\{multline\*?\}',
             lambda m: '\n$$\n' + m.group(1).strip() + '\n$$\n',
             tex, flags=re.DOTALL)

# === ENVIRONMENT CONVERSIONS ===

# Lists - itemize
tex = re.sub(r'\\begin\{itemize\}(.*?)\\end\{itemize\}',
             lambda m: '\n' + re.sub(r'\\item\s+', '- ', m.group(1)) + '\n',
             tex, flags=re.DOTALL)

# Lists - enumerate
tex = re.sub(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}',
             lambda m: '\n' + re.sub(r'\\item\s+', '1. ', m.group(1)) + '\n',
             tex, flags=re.DOTALL)

# Tables - convert tabular to Markdown pipes
# This is complex so just remove for now
tex = re.sub(r'\\begin\{table\}.*?\\end\{table\}', '', tex, flags=re.DOTALL)
tex = re.sub(r'\\begin\{center\}(.*?)\\end\{center\}', r'\1', tex, flags=re.DOTALL)
tex = re.sub(r'\\begin\{tabular\}.*?\\end\{tabular\}', '', tex, flags=re.DOTALL)

# Code blocks
tex = re.sub(r'\\begin\{lstlisting\}.*?\[language=([^\]]+)\](.*?)\\end\{lstlisting\}',
             lambda m: '\n```' + m.group(1) + '\n' + m.group(2) + '\n```\n',
             tex, flags=re.DOTALL)

# === CITATION REMOVAL ===

# Remove bibliography
tex = re.sub(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', '', tex, flags=re.DOTALL)

# Remove all citations
tex = re.sub(r'\\cite[ptPN]*\{[^}]*\}', '', tex)
tex = re.sub(r'\\eqref\{[^}]*\}', '', tex)
tex = re.sub(r'\\ref\{[^}]*\}', '', tex)
tex = re.sub(r'\\label\{[^}]*\}', '', tex)

# === TEXT FORMATTING ===

# Bold and italic
tex = re.sub(r'\\textbf\{(.*?)\}', r'**\1**', tex)
tex = re.sub(r'\\textit\{(.*?)\}', r'*\1*', tex)
tex = re.sub(r'\\emph\{(.*?)\}', r'*\1*', tex)
tex = re.sub(r'\\text\{(.*?)\}', r'\1', tex)
tex = re.sub(r'\\texttt\{(.*?)\}', r'`\1`', tex)

# === SPECIAL CHARACTERS ===

tex = tex.replace(r'\%', '%')
tex = tex.replace(r'\$', '$')
tex = tex.replace(r'\&', '&')
tex = tex.replace(r'\#', '#')
tex = tex.replace(r'\_', '_')
tex = tex.replace(r'\textdegree', '°')
tex = tex.replace(r'\degree', '°')
tex = tex.replace(r'~', ' ')  # tilde becomes space
tex = tex.replace(r'\,', ' ')  # thin space

# === UNICODE FIXES ===

tex = tex.replace('─', '-')  # U+2500
tex = tex.replace('│', '|')  # U+2502
tex = tex.replace('├', '+')  # U+251C
tex = tex.replace('└', '|')  # U+2514
tex = tex.replace('ΓöÇ', '-')  # mangled unicode

# === REMOVE REMAINING LaTeX COMMANDS ===

# Remove comment lines
tex = re.sub(r'^%.*$', '', tex, flags=re.MULTILINE)

# Remove orphaned braces
tex = re.sub(r'\{\}', '', tex)
tex = re.sub(r'\{([^}]*?)\}(?![a-z_])', r'\1', tex)  # Remove unnecessary braces

# Remove stray backslashes
tex = re.sub(r'\\([^a-zA-Z0-9])', r'\1', tex)  # Escaped special chars
tex = re.sub(r'\\$', '', tex)  # Trailing backslash

# Clean up double spaces
tex = re.sub(r'  +', ' ', tex)

# Clean up excessive newlines
tex = re.sub(r'\n\n\n+', '\n\n', tex)

# Write clean version
with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.write(tex)

print('[AGGRESSIVE] LaTeX removed, file cleaned')
