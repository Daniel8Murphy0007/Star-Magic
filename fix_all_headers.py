"""
fix_all_headers.py — Fix H1 headings across all whitepapers.

Actions:
  1. Papers whose first non-blank line is NOT a `# PAPER...` heading but has
     a `**Title:**` line → inject `# PAPER_XXX: {title}` at top.
  2. Papers whose first non-blank line starts with `#` but not `# PAPER` 
     (non-standard formats like `# Paper #24`, `# Whitepaper #27`, 
     `# PAPER #32b`) → normalise to `# PAPER_XXX: {title}`.
  3. Papers with PowerShell template corruption at top → strip corrupt lines
     then apply rule 1/2 as appropriate.
"""

import sys, glob, os, re
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# Regex patterns
PNUM_RE  = re.compile(r'PAPER[_\s#]*0*(\d+)', re.IGNORECASE)
H1_OK    = re.compile(r'^# PAPER[_\s#]', re.IGNORECASE)
ANY_H1   = re.compile(r'^#\s', re.IGNORECASE)
TITLE_RE = re.compile(r'^\*\*[Tt]itle[:：]\*\*\s*(.*)')
CORRUPT  = re.compile(r'PAPER_\{0:D3\}|\"PAPER_\{0|\[int\]# P', re.IGNORECASE)

changed = []
skipped = []
errors  = []

all_files = sorted(glob.glob('whitepapers/PAPER_*.md'))

for fpath in all_files:
    fname = os.path.basename(fpath)
    m = PNUM_RE.search(fname)
    if not m:
        errors.append(f'Cannot parse number: {fname}')
        continue
    num = int(m.group(1))
    num_str = f'{num:03d}'

    with open(fpath, encoding='utf-8', errors='replace') as fh:
        content = fh.read()

    lines = content.splitlines(keepends=True)

    # ── Step 1: strip corrupt PowerShell template lines ────────────────────
    stripped = False
    clean_lines = []
    for line in lines:
        if CORRUPT.search(line):
            stripped = True
        else:
            clean_lines.append(line)
    if stripped:
        lines = clean_lines
        content = ''.join(lines)

    # Re-check first non-blank line
    first_nb = next((l.rstrip() for l in lines if l.strip()), '')

    # ── Already correct ────────────────────────────────────────────────────
    if H1_OK.match(first_nb) and not stripped:
        skipped.append(fname)
        continue

    # ── Extract title ──────────────────────────────────────────────────────
    title = None

    # Case A: non-standard H1 exists (starts with # but wrong format)
    if ANY_H1.match(first_nb):
        # Strip leading # signs and common prefixes
        raw = re.sub(r'^#+\s*', '', first_nb)
        # Remove patterns like "Paper #24:", "PAPER #32b —", "Whitepaper #27 —"
        raw = re.sub(r'^(paper|whitepaper)\s*[#\s]*\d+[a-z]?\s*[:\-–—]\s*', '', raw, flags=re.IGNORECASE).strip()
        # Also try extracting from **Title:** if raw is empty/unhelpful
        if not raw:
            for line in lines:
                tm = TITLE_RE.match(line.strip())
                if tm:
                    raw = tm.group(1).strip()
                    break
        title = raw if raw else f'PAPER_{num_str}'
        new_h1 = f'# PAPER_{num_str}: {title}\n'
        # Replace first line
        new_lines = []
        replaced = False
        for line in lines:
            if not replaced and line.rstrip() == first_nb:
                new_lines.append(new_h1)
                replaced = True
            else:
                new_lines.append(line)
        new_content = ''.join(new_lines)

    else:
        # Case B: no H1 at all — find **Title:** and prepend
        for line in lines:
            tm = TITLE_RE.match(line.strip())
            if tm:
                title = tm.group(1).strip()
                break
        if not title:
            # Fall back to filename-based title
            title = f'PAPER_{num_str}'
        new_h1 = f'# PAPER_{num_str}: {title}\n\n'
        new_content = new_h1 + ''.join(lines)

    # Write back
    try:
        with open(fpath, 'w', encoding='utf-8') as fh:
            fh.write(new_content)
        changed.append(fname)
    except Exception as e:
        errors.append(f'{fname}: {e}')

print(f'Fixed   : {len(changed)}')
print(f'Already OK: {len(skipped)}')
print(f'Errors  : {len(errors)}')
if errors:
    for e in errors:
        print(f'  ERROR: {e}')
print()
print('Fixed papers:')
for f in changed[:20]:
    print(f'  {f}')
if len(changed) > 20:
    print(f'  ... and {len(changed)-20} more')
