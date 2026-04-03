"""
Corrected deduplication script for PAPER_001-156 whitepapers.
Algorithm: find second ## Abstract as the copy boundary, strip inter-copy artifacts,
preserve the UQFF Production Framework appendix if only in later copies.
"""
import re
import glob
import os

GARBLED_PATTERNS = [
    r'^\s*\.Groups\[1\]\.Value\s*$',
    r'^\s*"PAPER_\{0:D3\}"\s*-f\s*\$n\s*$',
    r'^\s*"PAPER_\{0:D3\}"\s*-f\s*\[int\]\s*$',
    r'^---\s*$',
]

def is_garbage_line(line):
    stripped = line.strip()
    if stripped == '':
        return True
    for pat in GARBLED_PATTERNS:
        if re.match(pat, line):
            return True
    return False

def strip_trailing_garbage(text):
    """Remove inter-copy artifact lines from the end of copy 1 content."""
    lines = text.split('\n')
    while lines and is_garbage_line(lines[-1]):
        lines.pop()
    return '\n'.join(lines)

def fix_garbled_title(text):
    """Fix garbled PowerShell/Python title line at the start of a document."""
    def _fix(m):
        # e.g.  #  "PAPER_{0:D3}" -f [int]# PAPER_001: ...
        # or    #  "PAPER_{0:D3}" -f $n# PAPER_001: ...
        line = m.group(0)
        # Extract the real title after the garbled prefix
        real = re.search(r'(# PAPER_[^\n]+)', line)
        if real:
            return '# ' + real.group(1)[2:]  # strip leading '# '
        return line
    text = re.sub(r'^#\s+"PAPER_\{0:D3\}"[^\n]*# PAPER_', '# PAPER_', text, flags=re.MULTILINE)
    text = re.sub(r'^#\s+"PAPER_\{0:D3\}"[^\n]+\n', '', text, flags=re.MULTILINE)
    return text

def extract_prod_appendix(text):
    """Extract UQFF Production Framework Reference appendix from text."""
    m = re.search(r'(## Appendix: UQFF Production Framework Reference.*)', text, re.DOTALL)
    return m.group(1).rstrip() if m else None

def deduplicate(filepath):
    content = open(filepath, 'r', encoding='utf-8', errors='replace').read()
    orig_len = len(content)

    abstracts = [m.start() for m in re.finditer(r'^## Abstract$', content, re.MULTILINE)]
    if len(abstracts) < 2:
        return (True, f"Already clean (abstracts={len(abstracts)})")

    n_dups = len(abstracts) - 1

    # Keep from start of file to just before second Abstract
    first_copy = content[:abstracts[1]]

    # Strip inter-copy garbage from end of first copy
    first_copy = strip_trailing_garbage(first_copy)

    # Fix garbled title line at the start
    first_copy = fix_garbled_title(first_copy)

    # Remove .Groups[1].Value and "PAPER_{0:D3}" artifacts from within the body
    first_copy = re.sub(r'^\.Groups\[1\]\.Value\s*$\n?', '', first_copy, flags=re.MULTILINE)
    first_copy = re.sub(r'^\s*"PAPER_\{0:D3\}"\s*-f\s*\$n\s*$\n?', '', first_copy, flags=re.MULTILINE)

    # Check if first copy already has UQFF Production Framework appendix
    has_prod = '## Appendix: UQFF Production Framework Reference' in first_copy

    # Extract it from the full document if needed
    prod_appendix = ''
    if not has_prod:
        pa = extract_prod_appendix(content)
        if pa:
            prod_appendix = '\n\n---\n\n' + pa

    result = first_copy.rstrip() + prod_appendix + '\n'
    new_len = len(result)

    if result == content:
        return (True, f"No change needed")

    open(filepath, 'w', encoding='utf-8').write(result)
    return (True, f"OK — removed {n_dups} duplicate(s), {orig_len - new_len:+,d} chars")


if __name__ == '__main__':
    files = sorted(glob.glob('whitepapers/PAPER_*.md'))
    target = [f for f in files
              if re.match(r'PAPER_(\d+)', os.path.basename(f))
              and 1 <= int(re.match(r'PAPER_(\d+)', os.path.basename(f)).group(1)) <= 156]

    fixed = 0
    skipped = 0
    errors = 0
    for f in target:
        ok, msg = deduplicate(f)
        label = '[FIXED]' if 'removed' in msg else '[SKIP]' if 'clean' in msg or 'No change' in msg else '[ERR]'
        if label == '[FIXED]':
            fixed += 1
        elif label == '[SKIP]':
            skipped += 1
        else:
            errors += 1
        print(f"  {label} {os.path.basename(f)}: {msg}")

    print(f"\nTotal: {fixed} fixed, {skipped} already clean, {errors} errors out of {len(target)} files")
