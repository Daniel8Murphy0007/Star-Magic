"""Fix the problematic Unicode supplement block comment in all 17 merged files,
then regenerate their PDFs."""

import sys, re
sys.path.insert(0, '.')
from pathlib import Path
from generate_pdfs import generate_pdf

WHITEPAPERS = Path('whitepapers')

# The bad Unicode comment (written by the merge script with f-string CANONICAL values)
BAD_COMMENT = (
    '*Merged from companion derivation file. Canonical UQFF constants: '
    '\u03ba=5.0e-4 day^-1, [SSq]=0.57, '
    '\u03b2_i=0.603, \u03c1_SCm=7.09e-37 J/m^3*'
)
GOOD_COMMENT = (
    '*Merged from companion derivation file. '
    'Canonical UQFF constants: kappa=5.0e-4/day, [SSq]=0.57, '
    'beta\\_i=0.603, rho\\_SCm=7.09e-37 J/m3.*'
)

fixes = 0
for n in range(1109, 1126):
    candidates = sorted(WHITEPAPERS.glob(f'PAPER_{n}_*.md'))
    if not candidates:
        print(f'  PAPER_{n}: no named file found')
        continue
    path = candidates[0]
    text = path.read_text(encoding='utf-8', errors='replace')
    if BAD_COMMENT in text:
        text = text.replace(BAD_COMMENT, GOOD_COMMENT)
        path.write_text(text, encoding='utf-8')
        fixes += 1
        print(f'  PAPER_{n}: fixed comment in {path.name}')
    else:
        print(f'  PAPER_{n}: comment not found (already fixed or different content)')

print(f'\nFixed {fixes} files. Now regenerating PDFs...\n')

fails = []
for n in range(1109, 1126):
    candidates = sorted(WHITEPAPERS.glob(f'PAPER_{n}_*.md'))
    if not candidates:
        print(f'  PAPER_{n}: no named file, skip')
        continue
    p = candidates[0]
    result = generate_pdf(p)
    success = result[2] if isinstance(result, tuple) and len(result) >= 3 else bool(result)
    detail  = result[3] if isinstance(result, tuple) and len(result) >= 4 else ''
    err     = result[4] if isinstance(result, tuple) and len(result) >= 5 else ''
    if success:
        print(f'  PAPER_{n}: OK  {detail}  -> {p.stem}.pdf')
    else:
        print(f'  PAPER_{n}: FAIL  {err[-300:]}')
        fails.append(p.name)

print(f'\nDone. Failures: {len(fails)}')
for f in fails:
    print(f'  FAIL: {f}')
