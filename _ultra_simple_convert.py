#!/usr/bin/env python3
"""Ultra-simple conversion: just remove LaTeX structures, keep content"""

import re

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.tex', 'r', encoding='utf-8') as f:
    content = f.read()

# Extract body
content = re.sub(r'^.*?\\begin\{document\}', '', content, flags=re.DOTALL)
content = re.sub(r'\\end\{document\}.*$', '', content, flags=re.DOTALL)

# Remove bibliography
content = re.sub(r'\\begin\{thebibliography\}.*?\\end\{thebibliography\}', '', content, flags=re.DOTALL)

# Remove citations/refs
content = re.sub(r'~?\\cite\{[^}]*\}', '', content)
content = re.sub(r'~?\\eqref\{[^}]*\}', '', content)
content = re.sub(r'~?\\label\{[^}]*\}', '', content)
content = re.sub(r'~?\\ref\{[^}]*\}', '', content)

# Remove \maketitle
content = content.replace(r'\maketitle', '')

# Convert sections
content = re.sub(r'\\begin\{abstract\}', '## Abstract', content)
content = re.sub(r'\\end\{abstract\}', '', content)
content = re.sub(r'\\section\{', '## {', content)
content = re.sub(r'\}(\s*%)', r'}\1', content)  # Preserve structure
content = re.sub(r'\\subsection\{', '### {', content)
content = re.sub(r'\\subsubsection\{', '#### {', content)

# Now fix the braces for headings - convert ### {text} to ### text
content = re.sub(r'^### \{([^}]*)\}', r'### \1', content, flags=re.MULTILINE)
content = re.sub(r'^## \{([^}]*)\}', r'## \1', content, flags=re.MULTILINE)
content = re.sub(r'^#### \{([^}]*)\}', r'#### \1', content, flags=re.MULTILINE)

# Convert equation/align blocks - SIMPLE approach: just keep content as-is
# Convert \begin{align}...\end{align} to $$...$$
content = re.sub(r'\\begin\{equation\}', '$$', content)
content = re.sub(r'\\end\{equation\}', '$$', content)
content = re.sub(r'\\begin\{align\*?\}', '$$', content)
content = re.sub(r'\\end\{align\*?\}', '$$', content)

# Within $$...$$, remove alignment markers
parts = re.split(r'(\$\$[^$]*\$\$)', content, flags=re.DOTALL)
for i in range(1, len(parts), 2):  # Only process math blocks
    if parts[i].startswith('$$'):
        # Remove & alignment markers
        math_content = parts[i][2:-2]  # Strip $$ markers
        math_content = math_content.replace(' &', '')  # Remove space-ampersand
        math_content = math_content.replace('&', '')  # Remove ampersand
        math_content = math_content.replace('\\\\', '  \n')  # Convert \\ to newline
        parts[i] = '$$' + math_content + '$$'
content = ''.join(parts)

# Convert itemize/enumerate
content = re.sub(r'\\begin\{itemize\}', '', content)
content = re.sub(r'\\end\{itemize\}', '', content)
content = re.sub(r'\\item\s+', '- ', content)

content = re.sub(r'\\begin\{enumerate\}', '', content)
content = re.sub(r'\\end\{enumerate\}', '', content)

# Text formatting
content = re.sub(r'\\textbf\{([^}]*)\}', r'**\1**', content)
content = re.sub(r'\\textit\{([^}]*)\}', r'*\1*', content)
content = re.sub(r'\\emph\{([^}]*)\}', r'*\1*', content)

# Special characters
content = content.replace('~', ' ')
content = content.replace(r'\ ', ' ')
content = content.replace(r'\textdegree', '°')
content = content.replace(r'---', '—')
content = content.replace(r'--', '–')

# Remove remaining LaTeX commands
content = re.sub(r'\\[a-zA-Z]+\{[^}]*\}', '', content)  # \cmd{...}
content = re.sub(r'\\[a-zA-Z]+', '', content)  # \cmd

# Remove comment lines
content = re.sub(r'^%.*$', '', content, flags=re.MULTILINE)

# Clean whitespace
content = re.sub(r'\n\n\n+', '\n\n', content)
content = '\n'.join(line.rstrip() for line in content.split('\n'))
content = content.strip()

with open('whitepapers/SCm_Rossi_ECat_Variants_Unified.md', 'w', encoding='utf-8') as f:
    f.write(content)

print('[ULTRA SIMPLE] Conversion complete')
