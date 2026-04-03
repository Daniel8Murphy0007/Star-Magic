#!/usr/bin/env python3
"""
Sample the heaviest UFFFD papers and analyze what characters were stripped.
Also check the superscript/subscript coverage in pdf_header.tex.
Also check NO_H1 papers to see which format style they use.
"""
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import glob, os, re, collections, json
pnum = re.compile(r'PAPER_(\d+)')
all_files = {int(pnum.search(os.path.basename(f)).group(1)): f
             for f in glob.glob('whitepapers/PAPER_*.md')
             if pnum.search(os.path.basename(f))}

with open('glyph_audit_results.json') as f:
    data = json.load(f)

# ── 1. Sample UFFFD context ──────────────────────────────────────────────────
heavy_ufffd = sorted([(int(k), v['UFFFD']) for k, v in data.items() if 'UFFFD' in v],
                     key=lambda x: -x[1])
print('=== UFFFD context samples (top 5 papers) ===')
for num, cnt in heavy_ufffd[:5]:
    fp = all_files.get(num)
    if not fp:
        continue
    with open(fp, 'r', encoding='utf-8', errors='replace') as f:
        text = f.read()
    lines = text.splitlines()
    print(f'\nPAPER_{num:04d} ({cnt} UFFDs) — first 8 UFFFD lines:')
    shown = 0
    for i, line in enumerate(lines, 1):
        if '\ufffd' in line:
            print(f'  L{i:4d}: {line.rstrip().encode("ascii","replace").decode()[:120]}')
            shown += 1
            if shown >= 8:
                break

# ── 2. What chars are MISSING from pdf_header.tex newunicodechar? ─────────────
print('\n=== Characters currently NOT covered by pdf_header.tex ===')

# Read what's currently in pdf_header.tex
with open('pdf_header.tex', 'r', encoding='utf-8') as f:
    header = f.read()

# Extract which codepoints are already handled
covered = set(re.findall(r'\\newunicodechar\{(.)\}', header))
print(f'Currently covered: {len(covered)} chars: ' + 
      ' '.join(f'U+{ord(c):04X}({c})' for c in sorted(covered, key=ord)))

# Now scan all files to find the full set of non-ASCII chars that appear
char_freq = collections.Counter()
for fp in all_files.values():
    with open(fp, 'r', encoding='utf-8', errors='replace') as f:
        text = f.read()
    for ch in text:
        o = ord(ch)
        if o > 127 and o != 0xfffd:
            char_freq[ch] += 1

print(f'\nTotal distinct non-ASCII chars in corpus: {len(char_freq)}')
print('\nTop 50 most frequent non-ASCII chars (not in pdf_header.tex):')
for ch, cnt in char_freq.most_common(60):
    if ch not in covered:
        cat = 'sup' if ord(ch) in range(0x2070,0x20a0) or ord(ch) in (0xb2,0xb3,0xb9) else \
              'sub' if ord(ch) in range(0x2080,0x20a0) else \
              'grk' if ord(ch) in range(0x0391,0x03ca) else \
              'math' if ord(ch) in range(0x2190,0x2400) else \
              'lat' if ord(ch) in range(0x00c0,0x0100) else 'oth'
        print(f'  U+{ord(ch):04X} {repr(ch).encode("ascii","replace").decode():10s} [{cat:4s}] count={cnt:6d}')

# ── 3. NO_H1 papers — what format style do they use? (first line types) ──────
print('\n=== NO_H1 papers — first non-blank line patterns ===')
no_h1 = sorted([int(k) for k, v in data.items() if 'NO_H1' in v])
first_line_types = collections.Counter()
for num in no_h1:
    fp = all_files.get(num)
    if not fp:
        continue
    with open(fp, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    first = next((l.strip() for l in lines if l.strip()), '')
    if first.startswith('**Title:**'):
        first_line_types['**Title:**'] += 1
    elif first.startswith('**'):
        first_line_types['**other**'] += 1
    elif first.startswith('<!--'):
        first_line_types['<!-- comment'] += 1
    else:
        first_line_types[repr(first[:30])] += 1

print('First-line patterns in NO_H1 papers:')
for pat, cnt in first_line_types.most_common():
    print(f'  {cnt:4d}  {pat}')
