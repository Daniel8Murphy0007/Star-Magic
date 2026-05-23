#!/usr/bin/env python3
"""Final cleanup of all SCm .md files - remove remaining LaTeX and Unicode issues"""

import re

files = [
    "whitepapers/SCm_Holmlid_Parkhomov_PonsFleischmann_Upgrade.md",
    "whitepapers/SCm_Holmlid_Rossi_Parkhomov_Validation.md",
    "whitepapers/SCm_Mizuno_LENR_Transmutation.md",
    "whitepapers/SCm_PonsFleischmann_Derivation.md",
    "whitepapers/SCm_Rossi_ECat_Variants_Unified.md"
]

for fpath in files:
    try:
        with open(fpath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        original = content
        
        # Remove any remaining \equation with labels
        content = re.sub(r'\\begin\{equation\}.*?\\end\{equation\}', '', content, flags=re.DOTALL)
        
        # Remove any remaining \align with labels
        content = re.sub(r'\\begin\{align\}.*?\\end\{align\}', '', content, flags=re.DOTALL)
        
        # Remove any remaining itemize/enumerate
        content = re.sub(r'\\begin\{itemize\}.*?\\end\{itemize\}', '', content, flags=re.DOTALL)
        content = re.sub(r'\\begin\{enumerate\}.*?\\end\{enumerate\}', '', content, flags=re.DOTALL)
        
        # Remove remaining \item commands
        content = re.sub(r'\\item\s+', '- ', content)
        
        # Remove \cite and \eqref references
        content = re.sub(r'~?\\cite\{[^}]*\}', '', content)
        content = re.sub(r'~?\\eqref\{[^}]*\}', '', content)
        
        # Remove remaining LaTeX commands that might cause issues
        content = re.sub(r'\\textdegree', 'degree', content)
        content = re.sub(r'\\textdash', '-', content)
        content = re.sub(r'\\newblock\s+', '', content)
        
        # Remove problematic Unicode characters (the line dash character U+2500)
        content = content.replace('\u2500', '-')
        
        # Clean up excessive newlines
        content = re.sub(r'\n\n\n+', '\n\n', content)
        
        if content != original:
            with open(fpath, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"[CLEANED] {fpath}")
        else:
            print(f"[OK] {fpath}")
    except Exception as e:
        print(f"[ERROR] {fpath}: {e}")
