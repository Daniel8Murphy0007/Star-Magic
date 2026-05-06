"""
Fix formatting in PAPER_1109-1125:
1. Remove bare `---` lines that bracket the short ## References section
   (pandoc simple_table parser creates a broken 1-col table at \real{0.0556} width)
2. Ensure wide tables (5+ columns) have proportional separators that give
   a minimum of 10% per column by redistributing dash counts.
Regenerates all 17 PDFs.
"""
import re
from pathlib import Path
import sys
sys.path.insert(0, '.')
from generate_pdfs import generate_pdf

# -----------------------------------------------------------------------
# Fix 1: Remove bare `---` surrounding ## References sections
# Pattern: blank lines, then ---, then ## References / list items, then ---
# Replace: just keep ## References and the list items, remove the `---` fences
# -----------------------------------------------------------------------

# Matches the `---` (alone on a line) immediately before ## References
REF_BEFORE = re.compile(r'\n---\n(## References\b)', re.MULTILINE)
# Matches the `---` (alone on a line) right after the references list ends
# (followed by blank line then ## Session or ## Supplementary or end-of-section)
REF_AFTER  = re.compile(r'(\n(?:- .+\n)+)\n---\n(\n## (?:Session|Supplementary))', re.MULTILINE)

# -----------------------------------------------------------------------
# Fix 2: For any pipe table separator where any column < 12% of total,
# redistribute so no column is narrower than 12% of total.
# -----------------------------------------------------------------------

def is_separator(line):
    s = line.strip()
    if not (s.startswith('|') and s.endswith('|')): return False
    cells = s[1:-1].split('|')
    return all(re.fullmatch(r'[\s\-:]+', c) for c in cells)

def get_dash_counts(sep_line):
    s = sep_line.strip()[1:-1].split('|')
    return [len(c.strip()) for c in s]

def rebuild_separator(sep_line, new_counts):
    """Rebuild separator using same alignment markers."""
    old_cells = sep_line.strip()[1:-1].split('|')
    new_cells = []
    for old, n in zip(old_cells, new_counts):
        t = old.strip()
        if t.startswith(':') and t.endswith(':'):
            new_cells.append(':' + '-'*(n-2) + ':')
        elif t.startswith(':'):
            new_cells.append(':' + '-'*(n-1))
        elif t.endswith(':'):
            new_cells.append('-'*(n-1) + ':')
        else:
            new_cells.append('-'*n)
    return '| ' + ' | '.join(new_cells) + ' |'

MIN_PCT = 0.12   # minimum 12% per column

def fix_wide_table_separators(text):
    """
    For pipe tables with 4+ columns where any column < MIN_PCT of total,
    redistribute dashes so each column gets at least MIN_PCT of total.
    """
    lines = text.split('\n')
    changed = 0
    i = 0
    while i < len(lines):
        if is_separator(lines[i]):
            counts = get_dash_counts(lines[i])
            if len(counts) < 4:
                i += 1
                continue
            total = sum(counts)
            fractions = [c/total for c in counts]
            if any(f < MIN_PCT for f in fractions):
                # Boost thin columns to MIN_PCT, scale others down proportionally
                # Use iterative approach: boost all thin cols, renormalize fat cols
                new_counts = list(counts)
                min_count = int(total * MIN_PCT)
                # First pass: set minimums
                for j in range(len(new_counts)):
                    if new_counts[j] < min_count:
                        new_counts[j] = min_count
                # Scale total: keep proportions of non-boosted columns
                # Recalculate proportions
                new_total = sum(new_counts)
                # If total grew too much, scale everything proportionally
                # so the sum stays the same as original to avoid changing page layout
                if new_total != total:
                    scale = total / new_total
                    new_counts = [max(min_count, int(c * scale)) for c in new_counts]
                    # Adjust last col to match total exactly
                    diff = total - sum(new_counts)
                    new_counts[-1] += diff
                # Only apply if truly changed
                if new_counts != counts:
                    lines[i] = rebuild_separator(lines[i], new_counts)
                    changed += 1
        i += 1
    return '\n'.join(lines), changed

# -----------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------
results = []

for n in range(1109, 1126):
    candidates = sorted(Path('whitepapers').glob(f'PAPER_{n}_*.md'))
    if not candidates:
        continue
    p = candidates[0]
    text = p.read_text(encoding='utf-8')
    fixes = []

    # Fix 1: remove `---` before ## References
    new_text, cnt = REF_BEFORE.subn(r'\n\1', text)
    if cnt:
        fixes.append(f'removed {cnt} --- before References')
    text = new_text

    # Fix 1b: remove `---` after references list
    new_text, cnt = REF_AFTER.subn(r'\1\2', text)
    if cnt:
        fixes.append(f'removed {cnt} --- after References list')
    text = new_text

    # Fix 2: fix wide table proportions
    text, wide_fixes = fix_wide_table_separators(text)
    if wide_fixes:
        fixes.append(f'redistributed {wide_fixes} wide-table separator(s)')

    if fixes:
        p.write_text(text, encoding='utf-8')
        fix_str = ', '.join(fixes)
        print(f'{p.name}: {fix_str}')
        result = generate_pdf(p)
        ok = 'OK' if result[2] else 'FAIL'
        detail = result[3] or result[4] or ''
        print(f'  PDF: {ok}  {detail}')
        results.append((p.name, ok))
    else:
        print(f'{p.name}: nothing to fix')

print(f'\n=== {sum(1 for _,s in results if s=="OK")}/{len(results)} PDFs regenerated ===')
