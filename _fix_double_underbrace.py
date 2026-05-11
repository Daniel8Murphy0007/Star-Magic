"""Fix double \\underbrace typo: \\underbrace{\\underbrace{X}_{label} -> \\underbrace{X}_{label}"""
import re, glob, pathlib
PAT = re.compile(r'\\underbrace\{\\underbrace\{')
count = 0
for fp in glob.glob('whitepapers/PAPER_*.md'):
    p = pathlib.Path(fp)
    s = p.read_text(encoding='utf-8')
    if '\\underbrace{\\underbrace{' not in s:
        continue
    new = PAT.sub(r'\\underbrace{', s)
    if new != s:
        p.write_text(new, encoding='utf-8')
        count += 1
        print(f'  fixed {p.name}')
print(f'Updated {count} files.')
