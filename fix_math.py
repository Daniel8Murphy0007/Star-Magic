"""Fix LaTeX math delimiters in upgraded papers 1109-1125.
Removes inner $...$ from display math $$...$$ blocks."""

import re
import glob

def fix_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    original = content
    
    # Fix 1: $$$content$$$ -> $$content$$  (triple dollar from $$ + $content$ + $$)
    # Pattern: $$$ at start of line, $$$ at end
    content = re.sub(r'^\$\$\$(.+?)\$\$\$', r'$$\1$$', content, flags=re.MULTILINE)
    
    # Fix 2: $$\boxed{$content$}$$ -> $$\boxed{content}$$
    content = re.sub(r'\$\$\\boxed\{\$(.+?)\$\}\$\$', r'$$\\boxed{\1}$$', content, flags=re.MULTILINE)
    
    # Fix 3: $$\boxed{$content$ (text)}$$ -> $$\boxed{content}$$\n\n(text)
    content = re.sub(
        r'\$\$\\boxed\{\$(.+?)\$\s*(\([^)]+\))\}\$\$',
        r'$$\\boxed{\1}$$\n\n\2',
        content, flags=re.MULTILINE
    )
    
    if content != original:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    return False

fixed = 0
for n in range(1109, 1126):
    files = glob.glob(f"whitepapers/PAPER_{n}*.md")
    for fp in files:
        if fix_file(fp):
            print(f"  FIXED: {fp}")
            fixed += 1
        else:
            print(f"  OK:    {fp}")

print(f"\n{fixed} files fixed.")
