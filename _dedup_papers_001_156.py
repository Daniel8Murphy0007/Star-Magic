"""
_dedup_papers_001_156.py
Remove duplicate document copies from PAPER_001-156 whitepapers.

Pattern: each duplicated file contains 2-4 identical copies of the full document,
created when upgrade scripts appended the entire document as a suffix.

Strategy:
  1. Detect duplication: count occurrences of '## Abstract'
  2. Find the start of the SECOND copy (the heading line '# PAPER_' just before
     the second '## Abstract')
  3. Check if the Appendix section appears ONLY in a later copy (not first copy)
  4. Keep: first copy + any Appendix unique to later copies
  5. Fix garbled title lines (PowerShell/Python format strings in title)
  6. Write back only if changed
"""

import glob, re, os, sys

WHITEPAPER_DIR = "whitepapers"

# Regex to find document title lines to use as copy-boundary markers
TITLE_RE  = re.compile(r'^#\s+(?:"PAPER_\{.*?\}"[^#]*#\s*)?PAPER_\d+\w*[:\s]', re.MULTILINE)
HEADER_RE = re.compile(
    r'#\s+(?:"[^"]*"\s*-f\s*[^\n#]+#\s*)?PAPER_(\d+\w*)',
    re.MULTILINE
)

def find_copy_starts(content):
    """Return character positions of each document header block."""
    # A copy starts when we see a heading '# PAPER_NNN...' followed eventually
    # by '## Abstract'.  Find by scanning backwards from each Abstract occurrence.
    abstract_positions = [m.start() for m in re.finditer(r'^## Abstract$', content, re.MULTILINE)]
    if len(abstract_positions) <= 1:
        return []   # no duplication

    boundaries = []
    for abs_pos in abstract_positions:
        # Walk backwards to find the nearest # PAPER_ heading
        snippet = content[:abs_pos]
        m = None
        for candidate in re.finditer(r'^#\s', snippet, re.MULTILINE):
            m = candidate
        if m:
            # Refine: the line must contain PAPER_
            line_start = m.start()
            line_end   = snippet.find('\n', line_start)
            title_line = snippet[line_start:line_end] if line_end >= 0 else snippet[line_start:]
            boundaries.append(line_start)
        else:
            boundaries.append(0)
    return boundaries


def fix_garbled_title(title_line):
    """Remove PowerShell / Python format artifacts from title line."""
    # Pattern 1: #  "PAPER_{0:D3}" -f [int]# PAPER_001: ...
    fixed = re.sub(r'^#\s+"PAPER_\{[^}]+\}"[^#]*#\s+', '# ', title_line)
    # Pattern 2: .Groups[1].Value : -> ignore (appears inside table cells only)
    return fixed.rstrip()


def extract_appendix(text):
    """Return the Appendix block if it exists, else ''."""
    m = re.search(r'\n---\n\n## Appendix', text)
    if m:
        return text[m.start():]
    m = re.search(r'\n## Appendix', text)
    if m:
        return '\n---\n' + text[m.start():]
    return ''


def deduplicate(filepath):
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        original = f.read()

    abstract_count = len(re.findall(r'^## Abstract$', original, re.MULTILINE))
    if abstract_count <= 1:
        return False, "clean (no duplication)"

    # Find copy boundaries
    boundaries = find_copy_starts(original)
    if not boundaries or len(boundaries) < 2:
        return False, f"could not split ({abstract_count} abstracts, {len(boundaries)} boundaries)"

    # First copy: from content[0] to boundaries[1]
    first_copy = original[:boundaries[1]].rstrip()

    # Check if Appendix exists in any later copy but not in first copy
    appendix_in_first = bool(re.search(r'^## Appendix', first_copy, re.MULTILINE) or
                              re.search(r'^### A\.', first_copy, re.MULTILINE))
    appendix_later = ''
    if not appendix_in_first:
        # Look for appendix in remainder
        remainder = original[boundaries[1]:]
        appendix_later = extract_appendix(remainder)

    # Fix garbled title line (first line of first_copy)
    lines = first_copy.split('\n')
    if lines and lines[0].startswith('#'):
        clean_title = fix_garbled_title(lines[0])
        if clean_title != lines[0]:
            lines[0] = clean_title
    first_copy = '\n'.join(lines)

    # Remove any garbled inline artifacts: '.Groups[1].Value : ...' in table cells
    # These appear as table cell content: '| ... |.Groups[1].Value : PAPER name\n'
    first_copy = re.sub(r'\|\.Groups\[1\]\.Value\s*:[^\n]*', '|', first_copy)
    first_copy = re.sub(r'\.Groups\[1\]\.Value\s*:[^\n]*\n', '\n', first_copy)

    # Build clean content
    clean_content = first_copy.rstrip()
    if appendix_later:
        clean_content += '\n' + appendix_later.rstrip()
    clean_content += '\n'

    if clean_content == original:
        return False, "no change needed"

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(clean_content)

    removed = len(original) - len(clean_content)
    return True, f"OK — removed {abstract_count-1} duplicate(s), -{removed:,} chars"


def main():
    files = sorted(glob.glob(f"{WHITEPAPER_DIR}/PAPER_*.md"))
    target = []
    for f in files:
        bn = os.path.basename(f)
        m = re.match(r'PAPER_(\d+)', bn)
        if m and 1 <= int(m.group(1)) <= 156:
            target.append(f)

    print(f"Scanning {len(target)} files in PAPER_001-156 range...\n")

    changed = 0
    errors  = 0
    skipped = 0

    for fpath in target:
        name = os.path.basename(fpath)
        try:
            ok, msg = deduplicate(fpath)
            if ok:
                changed += 1
                print(f"  [FIXED] {name}: {msg}")
            elif "clean" in msg or "no change" in msg:
                skipped += 1
            else:
                print(f"  [SKIP]  {name}: {msg}")
                skipped += 1
        except Exception as e:
            errors += 1
            print(f"  [ERR]   {name}: {e}")

    print(f"\n{'='*60}")
    print(f"  Fixed  : {changed}")
    print(f"  Skipped: {skipped}")
    print(f"  Errors : {errors}")
    print(f"  Total  : {len(target)}")
    print(f"{'='*60}")
    return errors


if __name__ == "__main__":
    sys.exit(main())
