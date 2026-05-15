"""
Fix corrupted ? characters in 7 Batch-3 MD files using their pristine git originals.
Strategy:
  - Load original (clean unicode) and current (corrupted with ?)
  - For each ? site in current, search original for the matching context window
    (with ? replaced by . regex wildcard); if a unique match found, extract the
    unicode char that sat there and substitute.
  - Bulk replacements for repeating patterns are also applied.
"""
import re, os, subprocess
from collections import Counter

FILES = {
 'PAPER_009b_Aether_String_TRZ_Damping_GW': 'ca2c552d',
 'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF': 'ca2c552d',
 'PAPER_016b_White_Dwarf_Foreground_UQFF': 'ca2c552d',
 'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime': '98b86da4',
 'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density': 'eae91218',
}
# PAPER_016 (3 ?s = real question sentences) and PAPER_018 (1 ? = real) are skipped.

# Small literal map applied AFTER per-site reconstruction. Anything still left.
GLOBAL_MAP = {
    # Common Greek/math ASCII fallbacks unique to this corpus
    # (kept empty - rely on originals; only used for context replacement)
}

def load_original(stem, sha):
    r = subprocess.run(['git', 'show', f'{sha}:whitepapers/{stem}.md'],
                       capture_output=True)
    return r.stdout.decode('utf-8', errors='replace')

def fix_one(stem, sha):
    cur_path = f'whitepapers/{stem}.md'
    cur = open(cur_path, encoding='utf-8').read()
    orig = load_original(stem, sha)

    # Find every '?' position in current
    positions = [m.start() for m in re.finditer(r'\?', cur)]
    cur_chars = list(cur)
    fixed = 0
    skipped = 0
    log = []

    for pos in positions:
        # Build a context window (15 chars before, 15 after) excluding nearby ?
        L, R = max(0, pos - 15), min(len(cur), pos + 16)
        before = cur[L:pos]
        after = cur[pos+1:R]
        # Strip other ? in the window (we want a clean anchor)
        if '?' in before or '?' in after:
            # shrink window to avoid other ?s
            before = before.split('?')[-1]
            after = after.split('?')[0]
        if len(before) < 4 and len(after) < 4:
            skipped += 1
            continue
        # Build regex: escape, then look in original
        pat = re.escape(before) + r'(.)' + re.escape(after)
        matches = re.findall(pat, orig)
        if not matches:
            # Try shorter window
            b2 = before[-8:] if len(before) >= 8 else before
            a2 = after[:8] if len(after) >= 8 else after
            if len(b2) < 3 and len(a2) < 3:
                skipped += 1; continue
            pat = re.escape(b2) + r'(.)' + re.escape(a2)
            matches = re.findall(pat, orig)
        if not matches:
            skipped += 1
            log.append(f'  NOMATCH  ...{before!r} ? {after!r}...')
            continue
        # Pick most common candidate
        cand = Counter(matches).most_common(1)[0][0]
        if cand == '?':
            skipped += 1
            log.append(f'  CIRCULAR ...{before!r} ? {after!r}... -> ?')
            continue
        cur_chars[pos] = cand
        fixed += 1

    new = ''.join(cur_chars)
    return new, fixed, skipped, log

if __name__ == '__main__':
    for stem, sha in FILES.items():
        new, fixed, skipped, log = fix_one(stem, sha)
        print(f'\n=== {stem} ===')
        print(f'  fixed={fixed}  skipped={skipped}  remaining_?={new.count("?")}')
        for line in log[:20]:
            print(line)
        # Don't write yet - dry run mode
        out = f'_fixed_md/{stem}.md'
        os.makedirs('_fixed_md', exist_ok=True)
        open(out, 'w', encoding='utf-8').write(new)
