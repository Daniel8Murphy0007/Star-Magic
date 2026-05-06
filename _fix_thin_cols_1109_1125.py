"""
Fix thin table columns in PAPER_1109-1125.
Each Markdown pipe table separator row (|---|---| line) determines
relative column widths in pdflatex output. We rewrite separator rows
so each column has width = max(15, len(widest_cell_in_that_column)).
"""
import re
from pathlib import Path
import sys
sys.path.insert(0, '.')
from generate_pdfs import generate_pdf

SEP_LINE = re.compile(r'^\|[\s\-:|]+\|[\s\-:||\s]*$')

def is_separator(line):
    stripped = line.strip()
    if not stripped.startswith('|') or not stripped.endswith('|'):
        return False
    cells = stripped[1:-1].split('|')
    return all(re.fullmatch(r'[\s\-:]+', c) for c in cells)

def parse_table(lines, sep_idx):
    """Return (header_idx, sep_idx, data_start, data_end) for a table block."""
    header_idx = sep_idx - 1
    if header_idx < 0 or not lines[header_idx].strip().startswith('|'):
        return None
    data_end = sep_idx + 1
    while data_end < len(lines) and lines[data_end].strip().startswith('|'):
        data_end += 1
    return header_idx, sep_idx, data_end

def cell_content_widths(row):
    stripped = row.strip()
    if not (stripped.startswith('|') and stripped.endswith('|')):
        return []
    return [len(c.strip()) for c in stripped[1:-1].split('|')]

def fix_tables_in_text(text):
    """Fix all thin separator rows in the text. Returns (new_text, fix_count)."""
    lines = text.split('\n')
    fix_count = 0
    i = 0
    while i < len(lines):
        if is_separator(lines[i]):
            sep_idx = i
            result = parse_table(lines, sep_idx)
            if result is None:
                i += 1
                continue
            header_idx, sep_idx, data_end = result

            # Collect all rows in this table (header + data)
            all_rows = [lines[header_idx]] + lines[sep_idx+1:data_end]

            # Get column count from separator
            sep_cells = lines[sep_idx].strip()[1:-1].split('|')
            ncols = len(sep_cells)

            # Find max content width per column across all rows
            max_widths = [0] * ncols
            for row in all_rows:
                widths = cell_content_widths(row)
                for col in range(min(len(widths), ncols)):
                    max_widths[col] = max(max_widths[col], widths[col])

            # Ensure minimum of 15 dashes per column
            target_widths = [max(15, w) for w in max_widths]

            # Build new separator row preserving alignment markers
            new_sep_cols = []
            for col, (old_cell, tw) in enumerate(zip(sep_cells, target_widths)):
                old = old_cell.strip()
                if old.startswith(':') and old.endswith(':'):
                    new_sep_cols.append(':' + '-' * (tw - 2) + ':')
                elif old.startswith(':'):
                    new_sep_cols.append(':' + '-' * (tw - 1))
                elif old.endswith(':'):
                    new_sep_cols.append('-' * (tw - 1) + ':')
                else:
                    new_sep_cols.append('-' * tw)

            new_sep = '| ' + ' | '.join(new_sep_cols) + ' |'
            if new_sep.strip() != lines[sep_idx].strip():
                lines[sep_idx] = new_sep
                fix_count += 1
            i = data_end
        else:
            i += 1
    return '\n'.join(lines), fix_count

results = []
total_fixes = 0

for n in range(1109, 1126):
    candidates = sorted(Path('whitepapers').glob(f'PAPER_{n}_*.md'))
    if not candidates:
        continue
    p = candidates[0]
    text = p.read_text(encoding='utf-8')
    new_text, fixes = fix_tables_in_text(text)
    if fixes > 0:
        p.write_text(new_text, encoding='utf-8')
        total_fixes += fixes
        print(f'{p.name}: fixed {fixes} table separator(s)')
        result = generate_pdf(p)
        ok = 'OK' if result[2] else 'FAIL'
        print(f'  PDF: {ok}  {result[3] or result[4] or ""}')
        results.append((p.name, ok))
    else:
        print(f'{p.name}: no thin separators found')

print(f'\n=== DONE: {total_fixes} separator rows fixed, {sum(1 for _,s in results if s=="OK")}/{len(results)} PDFs regenerated OK ===')
