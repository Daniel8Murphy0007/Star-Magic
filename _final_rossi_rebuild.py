#!/usr/bin/env python3
"""Final rebuild of Rossi .md - extract body, clean aggressively"""

import re

# Read .tex
with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.tex', 'r', encoding='utf-8') as f:
    tex = f.read()

# Remove preamble
tex = re.sub(r'^.*?\\begin\{document\}', '', tex, flags=re.DOTALL)
tex = re.sub(r'\\end\{document\}.*$', '', tex, flags=re.DOTALL)

# Remove \maketitle
tex = tex.replace(r'\maketitle', '')

# === STRIP ALL REMAINING LaTeX COMMANDS ===
# This is aggressive but necessary for pdflatex compatibility

# Remove \textbf, \textit, etc. - keep content
tex = re.sub(r'\\text[a-z]+\{(.*?)\}', r'\1', tex)
tex = re.sub(r'\\emph\{(.*?)\}', r'\1', tex)

# Remove \frac{a}{b} - replace with a/b
tex = re.sub(r'\\frac\{([^}]*)\}\{([^}]*)\}', r'\1/\2', tex)

# Remove ALL other backslash commands (except those in $$...$$)
# Split by $$...$$
parts = re.split(r'(\$\$.*?\$\$)', tex, flags=re.DOTALL)
result = []
for i, part in enumerate(parts):
    if i % 2 == 1:  # This is math mode
        result.append(part)
    else:  # This is regular text
        # Remove backslash commands
        part = re.sub(r'\\[a-zA-Z]+\{[^}]*\}', '', part)  # \cmd{...}
        part = re.sub(r'\\[a-zA-Z]+\[[^\]]*\]', '', part)  # \cmd[...]
        part = re.sub(r'\\[a-zA-Z]+', '', part)  # \cmd
        result.append(part)
tex = ''.join(result)

# === FIX MATH MODE ===

# Within $$...$$, fix common issues
# Remove \\ from equations
tex = re.sub(r'(\$\$.*?)\\\\\n', r'\1\n', tex, flags=re.DOTALL)

# Replace & with newline in equations
tex = re.sub(r'(\$\$.*?)&', r'\1\n', tex, flags=re.DOTALL)

# Fix broken subscripts/superscripts
tex = re.sub(r'_([0-9]+)', r'_{\1}', tex)  # Convert _123 to _{123}
tex = re.sub(r'\^([0-9]+)', r'^{\1}', tex)  # Convert ^123 to ^{123}

# === SPECIAL CHARACTER FIXES ===

tex = tex.replace('~', ' ')
tex = tex.replace('--', '–')  # en-dash
tex = tex.replace('---', '—')  # em-dash

# === REMOVE BIBLIOGRAPHY ===

tex = re.sub(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', '', tex, flags=re.DOTALL)

# === CLEAN UP ===

# Remove empty braces
tex = re.sub(r'\{\}', '', tex)

# Remove orphaned braces
tex = re.sub(r'(?<![\\$_^])(\{)([^}]*)(\})', r'\2', tex)

# Multiple spaces
tex = re.sub(r'  +', ' ', tex)

# Multiple newlines
tex = re.sub(r'\n\n\n+', '\n\n', tex)

# Strip leading/trailing spaces on lines
lines = tex.split('\n')
lines = [line.rstrip() for line in lines]
tex = '\n'.join(lines)

# Write
with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.write(tex)

print('[FINAL REBUILD] Complete aggressive LaTeX removal')
