#!/usr/bin/env python3
"""Remove bibliography sections from SCm .md files"""

files = [
    "whitepapers/SCm_PonsFleischmann_Derivation.md",
    "whitepapers/SCm_Rossi_ECat_Variants_Unified.md"
]

for fpath in files:
    try:
        with open(fpath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Remove \begin{thebibliography}...\end{thebibliography}
        if '\\begin{thebibliography}' in content:
            start = content.find('\\begin{thebibliography}')
            end = content.find('\\end{thebibliography}')
            if start >= 0 and end >= 0:
                # Keep everything before start, and after end + len('...')
                before = content[:start].rstrip()
                after = content[end + len('\\end{thebibliography}'):].lstrip()
                content = before + "\n" + after
        
        # Remove \cite{...} references - replace with empty or simple text
        import re
        content = re.sub(r'~?\\cite\{[^}]*\}', '', content)
        content = re.sub(r'~?\\eqref\{[^}]*\}', '', content)
        
        with open(fpath, 'w', encoding='utf-8') as f:
            f.write(content)
        
        print(f"[OK] {fpath}")
    except Exception as e:
        print(f"[ERROR] {fpath}: {e}")
