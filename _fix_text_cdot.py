"""Bulk-fix invalid \\text{X\\cdot Y} → \\text{X}\\cdot\\text{Y} across whitepapers."""
import pathlib, re

ROOT = pathlib.Path('whitepapers')

# Match \text{ ... \cdot ... } and split. Handle nested-free single-level only.
TEXT_CDOT = re.compile(r'\\text\{([^{}]*?)\\cdot\s*([^{}]*?)\}')

def split_once(m):
    a = m.group(1).strip()
    b = m.group(2).strip()
    # If either side empty, leave alone
    if not a or not b:
        return m.group(0)
    return f'\\text{{{a}}}\\cdot\\text{{{b}}}'

count_files = 0
count_subs = 0
for p in ROOT.glob('PAPER_*.md'):
    s = p.read_text(encoding='utf-8')
    new_s, n = TEXT_CDOT.subn(split_once, s)
    # Repeat to catch chained \cdot
    while n > 0:
        s2 = new_s
        new_s, n = TEXT_CDOT.subn(split_once, s2)
        if new_s == s2:
            break
    if new_s != s:
        p.write_text(new_s, encoding='utf-8')
        count_files += 1
        count_subs += 1
        print(f'FIXED {p.name}')

print(f'\nTotal files modified: {count_files}')
