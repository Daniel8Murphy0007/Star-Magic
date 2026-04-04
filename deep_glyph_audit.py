#!/usr/bin/env python3
"""
deep_glyph_audit.py
Scan all whitepaper markdown files for glyph corruption patterns.
Produces a categorized report with per-file counts.
"""
import glob, os, re, collections, json

WP_DIR = 'whitepapers'
pnum = re.compile(r'PAPER_(\d+)')

# ── Pattern definitions ──────────────────────────────────────────────────────
# Each entry: (category_name, regex_or_bytes, description)
PATTERNS = [
    # 1. Unicode replacement character U+FFFD (shows as □ or ? in PDFs)
    ('UFFFD',         re.compile(r'\ufffd'),
     'U+FFFD replacement char (shows as □/? in PDF)'),

    # 2. Literal question mark standing in for Greek/math chars
    #    Common: κ→?, β→?, ρ→?, α→?, Ω→?, θ→?, Δ→?, σ→?, π→?, λ→?, μ→?
    #    Detected as: standalone ? surrounded by math context
    ('BARE_QMARK',    re.compile(r'(?<!\S)\?(?!\?|\s*PASS|\s*ACCEPT|\s*FAIL|\s*=|\d)|\b\?\s*=\s*[0-9]|\?\s*[_^({]|[_^(]?\s*\?[\s,)]'),
     'Bare ? in math/formula context (likely stripped Greek letter)'),

    # 3. Superscript Unicode outside math (U+2070-U+209F, U+00B2/B3/B9)
    ('SUPERSCRIPT_U', re.compile(r'[\u2070-\u209f\u00b2\u00b3\u00b9]'),
     'Unicode superscript chars (may not render in all fonts)'),

    # 4. Subscript Unicode (U+2080-U+209F not caught above, U+1D00-U+1DBF)
    ('SUBSCRIPT_U',   re.compile(r'[\u2080-\u209f\u1d62\u1d63\u1d64\u1d65\u1d6e\u1d70\u2099]'),
     'Unicode subscript chars'),

    # 5. Solar symbol ☉ U+2609
    ('SOLAR',         re.compile(r'\u2609'),
     'Solar symbol ☉ (U+2609) — handled by pdf_header.tex newunicodechar'),

    # 6. Checkmark ✓ U+2713
    ('CHECKMARK',     re.compile(r'\u2713'),
     'Checkmark ✓ (U+2713) — handled by pdf_header.tex'),

    # 7. Arrow/math symbols not in DejaVu Serif
    ('MATH_UNICODE',  re.compile(r'[\u2190-\u21ff\u2200-\u22ff\u2300-\u23ff]'),
     'Math/arrow Unicode (U+2190-23FF) — may need font fallback'),

    # 8. Greek letters in running text (outside $$) — fine in xelatex with DejaVu but check
    ('GREEK_TEXT',    re.compile(r'[^\$\`]{0,3}[\u0391-\u03c9]{1,3}[^\$\`]{0,3}'),
     'Greek letters in text context (not in math block)'),

    # 9. Degree sign, multiplication × etc.
    ('LATIN_SUPP',    re.compile(r'[\u00b0\u00d7\u00f7\u00b1\u00b7]'),
     'Latin-1 supplement chars (°×÷±·)'),

    # 10. En-dash / em-dash / smart quotes
    ('DASHES_QUOTES', re.compile(r'[\u2013\u2014\u2018\u2019\u201c\u201d]'),
     'En/em-dash or smart quotes'),

    # 11. Narrow no-break space U+202F / non-breaking space U+00A0
    ('NBSP',          re.compile(r'[\u00a0\u202f\u2009\u200b]'),
     'Non-breaking/thin/zero-width spaces'),

    # 12. Missing H1 title (first non-blank line is not # heading)
    ('NO_H1',         None,   # handled specially
     'No H1 (#) title line — paper ID missing from PDF header'),

    # 13. Literal \n or \t that ended up in text (from bad script injection)
    ('CTRL_IN_TEXT',  re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f]'),
     'Raw control characters (should have been fixed)'),
]

results = {}  # paper_num -> {category: count, ...}
all_files = sorted(glob.glob(os.path.join(WP_DIR, 'PAPER_*.md')))

for fpath in all_files:
    m = pnum.search(os.path.basename(fpath))
    if not m:
        continue
    num = int(m.group(1))
    
    try:
        with open(fpath, 'r', encoding='utf-8', errors='replace') as f:
            text = f.read()
            lines = text.splitlines()
    except Exception as e:
        results[num] = {'READ_ERROR': 1}
        continue

    rec = {}

    for cat, pat, desc in PATTERNS:
        if cat == 'NO_H1':
            # Check if first non-blank line starts with '# '
            first = next((l for l in lines if l.strip()), '')
            if not first.startswith('# '):
                rec['NO_H1'] = 1
        elif pat is not None:
            count = len(pat.findall(text))
            if count:
                rec[cat] = count

    if rec:
        results[num] = rec

# ── Summary stats ─────────────────────────────────────────────────────────────
total_papers = len(all_files)
category_totals = collections.Counter()
category_paper_counts = collections.Counter()

for num, rec in results.items():
    for cat, cnt in rec.items():
        category_totals[cat] += cnt
        category_paper_counts[cat] += 1

print(f'Total papers scanned: {total_papers}')
print(f'Papers with at least one issue: {len(results)}')
print()
print(f'{"Category":<20} {"Papers affected":>16} {"Total occurrences":>18}  Description')
print('-' * 95)
for cat, desc in [(c, d) for c, _, d in PATTERNS]:
    p = category_paper_counts.get(cat, 0)
    t = category_totals.get(cat, 0)
    if p:
        print(f'{cat:<20} {p:>16,} {t:>18,}  {desc[:55]}')

# ── Top-10 worst papers ────────────────────────────────────────────────────────
print()
print('Top 20 papers with most issues (by total issue count):')
sorted_results = sorted(results.items(), key=lambda kv: sum(kv[1].values()), reverse=True)
for num, rec in sorted_results[:20]:
    total = sum(rec.values())
    cats = ', '.join(f'{k}:{v}' for k, v in sorted(rec.items(), key=lambda x: -x[1])[:4])
    print(f'  PAPER_{num:04d}  total={total:>5}  [{cats}]')

# Save full results to JSON for later processing
with open('glyph_audit_results.json', 'w', encoding='utf-8') as f:
    json.dump({str(k): v for k, v in sorted(results.items())}, f, indent=2)
print()
print('Full results saved to glyph_audit_results.json')
