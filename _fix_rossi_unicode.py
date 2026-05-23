#!/usr/bin/env python3
"""Fix Unicode and LaTeX issues in Rossi .md file"""

import re

# Read file
with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'r', encoding='utf-8') as f:
    content = f.read()

# Fix Unicode characters
content = content.replace('─', '-')  # U+2500 box drawing line
content = content.replace('│', '|')  # U+2502 box drawing vertical
content = content.replace('├', '|')  # U+251C box drawing
content = content.replace('└', '|')  # U+2514 box drawing
content = content.replace('─', '-')  # Additional pass
content = content.replace('ΓöÇ', '-')  # The exact error from pdflatex
content = content.replace('´', "'")   # accent characters
content = content.replace('`', "'")   # backticks

# Fix LaTeX commands
content = content.replace(r'\textdegree', '°')  # degree symbol
content = content.replace(r'\degree', '°')     # alternate
content = content.replace(r'\%', '%')          # percent
content = content.replace(r'\$', '$')          # dollar (but keep math mode $)
content = content.replace(r'\{', '{')          # braces
content = content.replace(r'\}', '}')

# Remove remaining common LaTeX commands
content = re.sub(r'\\emph\{(.*?)\}', r'*\1*', content)  # emphasis
content = re.sub(r'\\textbf\{(.*?)\}', r'**\1**', content)  # bold
content = re.sub(r'\\textit\{(.*?)\}', r'*\1*', content)  # italic
content = re.sub(r'\\text\{(.*?)\}', r'\1', content)  # text mode

# Remove comment lines that only contain dashes
content = re.sub(r'^% *-+.*$', '', content, flags=re.MULTILINE)
content = re.sub(r'^% *[─┼├┤└┘┐┌│─═╪╫╬]+.*$', '', content, flags=re.MULTILINE)

# Remove stray backslashes
content = re.sub(r'\\(?![a-zA-Z])', '', content)

# Clean up spacing
content = re.sub(r'\n\n\n+', '\n\n', content)

# Write clean file
with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.write(content)

print('[FIXED] Unicode and LaTeX commands removed')
