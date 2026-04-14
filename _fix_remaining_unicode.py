"""Replace remaining problematic Unicode: U+2299 (⊙) and others."""
import os

replacements = {
    '\u2299': 'M_sun',   # ⊙ circled dot operator / solar symbol
}

count = 0
for fn in sorted(os.listdir('whitepapers')):
    if not fn.endswith('.md'):
        continue
    path = os.path.join('whitepapers', fn)
    with open(path, encoding='utf-8') as f:
        content = f.read()
    orig = content
    for char, repl in replacements.items():
        content = content.replace(char, repl)
    if content != orig:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        count += 1

print(f"Fixed remaining Unicode in {count} papers")
