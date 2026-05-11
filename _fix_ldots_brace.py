"""Fix \\ldots} missing closing brace → \\ldots\\}"""
import pathlib, re
ROOT = pathlib.Path('whitepapers')
# Pattern: \{...\ldots} where the } is meant to be escaped \}
# Match: \{ ... \ldots\s*\} but only when the } is bare (not preceded by \)
PAT = re.compile(r'\\\{([^{}]*?\\ldots)\}')
n_total = 0
for p in ROOT.glob('PAPER_*.md'):
    s = p.read_text(encoding='utf-8')
    new_s, n = PAT.subn(lambda m: '\\{' + m.group(1) + '\\}', s)
    if n:
        p.write_text(new_s, encoding='utf-8')
        print(f'{p.name}: {n}')
        n_total += n
print(f'TOTAL fixes: {n_total}')
