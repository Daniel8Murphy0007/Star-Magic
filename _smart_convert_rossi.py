#!/usr/bin/env python3
"""Smart conversion: proper LaTeX to Markdown"""

import re

def convert_itemize(content):
    """Convert itemize to markdown bullets"""
    content = re.sub(r'\\item\s+', '- ', content)
    return content

def convert_enumerate(content):
    """Convert enumerate to markdown numbered"""
    items = re.findall(r'\\item\s+(.*?)(?=\\item|$)', content, flags=re.DOTALL)
    result = []
    for i, item in enumerate(items, 1):
        result.append(f'{i}. {item.strip()}')
    return '\n'.join(result) + '\n'

def clean_math_block(math):
    """Clean math - remove \\\\, &, preserve subscripts/superscripts"""
    # Remove equation labels
    math = re.sub(r'\\label\{[^}]*\}', '', math)
    # Convert \\ and & to newlines (for align blocks)
    math = math.replace('\\\\', '\n')
    math = math.replace('&', '\n')
    return f'$${math}$$'

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.tex', 'r', encoding='utf-8') as f:
    tex = f.read()

# Remove preamble
tex = re.sub(r'^.*?\\begin\{document\}', '', tex, flags=re.DOTALL)
tex = re.sub(r'\\end\{document\}.*$', '', tex, flags=re.DOTALL)

# ===== REMOVE BIBLIOGRAPHY =====
tex = re.sub(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', '', tex, flags=re.DOTALL)

# ===== REMOVE ALL \cite, \eqref, \label, \ref =====
tex = re.sub(r'~?\\cite\{[^}]*\}', '', tex)
tex = re.sub(r'~?\\eqref\{[^}]*\}', '', tex)
tex = re.sub(r'~?\\label\{[^}]*\}', '', tex)
tex = re.sub(r'~?\\ref\{[^}]*\}', '', tex)

# ===== CONVERT DOCUMENT COMMANDS =====
tex = tex.replace(r'\maketitle', '')
tex = re.sub(r'\\section\{([^}]*)\}', r'## \1', tex)
tex = re.sub(r'\\subsection\{([^}]*)\}', r'### \1', tex)
tex = re.sub(r'\\subsubsection\{([^}]*)\}', r'#### \1', tex)

# ===== CONVERT ENVIRONMENTS =====
# abstract
tex = re.sub(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', r'## Abstract\n\n\1', tex, flags=re.DOTALL)

# itemize
tex = re.sub(r'\\begin\{itemize\}(.*?)\\end\{itemize\}', lambda m: convert_itemize(m.group(1)), tex, flags=re.DOTALL)

# enumerate
tex = re.sub(r'\\begin\{enumerate\}(.*?)\\end\{enumerate\}', lambda m: convert_enumerate(m.group(1)), tex, flags=re.DOTALL)

# equations - keep as $$...$$ blocks
tex = re.sub(r'\\begin\{equation\}(.*?)\\end\{equation\}', r'$$\1$$', tex, flags=re.DOTALL)
tex = re.sub(r'\\begin\{align\}(.*?)\\end\{align\}', r'$$\1$$', tex, flags=re.DOTALL)
tex = re.sub(r'\\begin\{align\*\}(.*?)\\end\{align\*\}', r'$$\1$$', tex, flags=re.DOTALL)
tex = re.sub(r'\$\$(.*?)\$\$', lambda m: clean_math_block(m.group(1)), tex, flags=re.DOTALL)

# lstlisting
tex = re.sub(r'\\begin\{lstlisting\}\[language=([^\]]*)\](.*?)\\end\{lstlisting\}', 
             r'```\1\n\2\n```', tex, flags=re.DOTALL)
tex = re.sub(r'\\begin\{lstlisting\}(.*?)\\end\{lstlisting\}', 
             r'```\n\1\n```', tex, flags=re.DOTALL)

# ===== TEXT FORMATTING =====
# \textbf, \textit, \emph
tex = re.sub(r'\\textbf\{([^}]*)\}', r'**\1**', tex)
tex = re.sub(r'\\textit\{([^}]*)\}', r'*\1*', tex)
tex = re.sub(r'\\emph\{([^}]*)\}', r'*\1*', tex)
tex = re.sub(r'\\text\{([^}]*)\}', r'\1', tex)

# ===== SPECIAL CHARACTERS =====
tex = tex.replace('~', ' ')  # non-breaking space
tex = tex.replace(r'\ ', ' ')  # escaped space
tex = tex.replace(r'\textdegree', '°')
tex = tex.replace(r'---', '—')  # em-dash
tex = tex.replace(r'--', '–')   # en-dash

# ===== REMOVE REMAINING LaTeX COMMANDS (but not subscripts!) =====
# Remove \(...), \[...], \(...\)
tex = re.sub(r'\\\(', '$', tex)
tex = re.sub(r'\\\)', '$', tex)
tex = re.sub(r'\\\[', '$$', tex)
tex = re.sub(r'\\\]', '$$', tex)

# Remove backslash commands that take args but NOT subscripts/superscripts
tex = re.sub(r'\\(frac|sqrt|left|right|big[lr]?|Text[a-z]*)\{', r'{', tex)

# Remove remaining backslash commands (but carefully)
# \usepackage, \includegraphics, etc. - remove whole command
tex = re.sub(r'\\(usepackage|documentclass|geometry|hypersetup|title|author|date|maketitle|pagestyle|setlength|addtolength)\{[^}]*\}', '', tex)
tex = re.sub(r'\\(setlength|addtolength)\{[^}]*\}\{[^}]*\}', '', tex)

# Remove other backslash commands (except those in $$...$$)
parts = re.split(r'(\$\$.*?\$\$)', tex, flags=re.DOTALL)
result_parts = []
for i, part in enumerate(parts):
    if i % 2 == 1:  # Math mode
        result_parts.append(part)
    else:  # Regular text
        # Remove one-off commands
        part = re.sub(r'\\(hline|midrule|toprule|bottomrule|centering|small|large|Large|tiny|scriptsize|footnotesize)', '', part)
        # Remove comment lines
        part = re.sub(r'^%.*$', '', part, flags=re.MULTILINE)
        result_parts.append(part)
tex = ''.join(result_parts)

# ===== CLEAN UP =====
# Remove multiple empty lines
tex = re.sub(r'\n\n\n+', '\n\n', tex)

# Remove trailing spaces
lines = tex.split('\n')
lines = [line.rstrip() for line in lines]
tex = '\n'.join(lines)

# Strip leading/trailing
tex = tex.strip()

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.write(tex)

print('[SMART CONVERT] LaTeX properly converted to Markdown')
