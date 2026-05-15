"""
Reconstruct ? characters in CURRENT files by aligning to pre-corruption
versions from git history.

For each file:
  1. Walk git history to find pre/post corruption commits.
  2. Pull the PRE version (before ? corruption).
  3. For each line in CURRENT containing '?':
     - Strip markdown formatting differences.
     - Search PRE for the line whose normalized form best matches.
     - Extract the unicode chars at ? positions.
  4. Apply substitutions to current MD; write to _fixed_md/.
"""
import os, re, subprocess, difflib
from collections import Counter
import unicodedata

# (stem, pre_corruption_sha) - sha BEFORE the corruption that introduced ?'s
# Determined by walking history; first commit where ?-count jumps significantly.
TARGETS = [
 ('PAPER_009b_Aether_String_TRZ_Damping_GW',                        None),
 ('PAPER_015b_Multiband_GW_LISA_LIGO_UQFF',                          None),
 ('PAPER_016b_White_Dwarf_Foreground_UQFF',                          None),
 ('PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime',                 None),
 ('PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density', None),
]

def find_pre_corruption_sha(stem):
    """Walk history; return SHA of last commit BEFORE first significant ?-jump (>=10)."""
    path = f'whitepapers/{stem}.md'
    r = subprocess.run(['git','log','--reverse','--format=%H','--', path],
                       capture_output=True, text=True)
    shas = r.stdout.strip().split('\n')
    prev_q = 0
    prev_sha = None
    for sha in shas:
        r2 = subprocess.run(['git','show', f'{sha}:{path}'], capture_output=True)
        txt = r2.stdout.decode('utf-8', errors='replace')
        qc = txt.count('?')
        if qc - prev_q >= 5:  # significant jump
            return prev_sha, sha
        prev_q = qc
        prev_sha = sha
    return None, None

def load_at(sha, path):
    if sha is None:
        return None
    r = subprocess.run(['git','show', f'{sha}:{path}'], capture_output=True)
    return r.stdout.decode('utf-8', errors='replace')

def normalize(s):
    """Normalize a line for fuzzy matching: collapse whitespace, lowercase."""
    return re.sub(r'\s+', ' ', s).strip().lower()

def fix_file(stem):
    path = f'whitepapers/{stem}.md'
    pre_sha, post_sha = find_pre_corruption_sha(stem)
    cur = open(path, encoding='utf-8').read()
    pre = load_at(pre_sha, path) if pre_sha else None
    print(f'\n=== {stem} ===')
    print(f'  pre_sha={pre_sha}  post_sha={post_sha}')
    if not pre:
        print('  NO PRE - cannot fix')
        return cur, 0
    print(f'  PRE chars={len(pre):6} ?-count={pre.count("?"):3}')
    print(f'  CUR chars={len(cur):6} ?-count={cur.count("?"):3}')

    # Build a list of PRE lines for alignment
    pre_lines = pre.split('\n')
    pre_norm_idx = {}  # normalized form -> line index (first occurrence)
    for i, line in enumerate(pre_lines):
        n = normalize(line)
        if n and n not in pre_norm_idx:
            pre_norm_idx[n] = i

    cur_lines = cur.split('\n')
    fixed_lines = []
    total_fixed = 0
    total_unfixed = 0
    unfixed_log = []

    for li, line in enumerate(cur_lines):
        if '?' not in line:
            fixed_lines.append(line)
            continue
        # Search pre for matching line
        n = normalize(line.replace('?', ''))  # remove ?s for matching
        # Make a regex from line where ? matches any single non-whitespace
        # Actually use difflib to find best matching pre line
        cur_norm = normalize(line)
        # Replace ? with . regex wildcard
        regex_line = re.escape(cur_norm).replace(r'\?', '.')
        candidates = []
        for i, pline in enumerate(pre_lines):
            pn = normalize(pline)
            if not pn or len(pn) < 5:
                continue
            if re.search(regex_line, pn):
                candidates.append((i, pline, pn))
        if not candidates:
            # Try short windows around each ?
            new_line = line
            for m in list(re.finditer(r'\?', line))[::-1]:
                pos = m.start()
                # Get window
                before = line[max(0,pos-12):pos]
                after = line[pos+1:pos+13]
                # Skip if both windows trivially small
                if len(before) < 3 and len(after) < 3:
                    continue
                # Build regex
                pat = re.escape(before) + r'(.{0,4})' + re.escape(after)
                ms = re.findall(pat, pre)
                if ms:
                    cand = Counter(ms).most_common(1)[0][0]
                    if cand and '?' not in cand and cand != '':
                        new_line = new_line[:pos] + cand + new_line[pos+1:]
                        total_fixed += 1
                        continue
                total_unfixed += 1
                unfixed_log.append(f'    line{li}: ...{before!r} ? {after!r}...')
            fixed_lines.append(new_line)
        else:
            # Use best match (most similar in length)
            best = min(candidates, key=lambda x: abs(len(x[2]) - len(cur_norm)))
            pre_line = best[1]
            # Align word-by-word: find ? positions in cur and lookup chars in pre
            new_line = line
            # Use regex with capture groups
            cur_pattern = re.escape(cur_norm)
            cur_pattern = cur_pattern.replace(r'\?', r'(.)')
            m = re.search(cur_pattern, normalize(pre_line))
            if m:
                replacements = list(m.groups())
                # Replace ? in original (case-preserved) line
                idx = 0
                for q_match in re.finditer(r'\?', line):
                    if idx < len(replacements):
                        rep = replacements[idx]
                        if rep and '?' not in rep:
                            qp = q_match.start()
                            new_line = new_line[:qp] + rep + new_line[qp+1:]
                            total_fixed += 1
                        idx += 1
            fixed_lines.append(new_line)

    new = '\n'.join(fixed_lines)
    print(f'  Fixed: {total_fixed}  Unfixed: {total_unfixed}  Remaining ?: {new.count("?")}')
    if unfixed_log:
        print('  Unfixed sites:')
        for l in unfixed_log[:15]:
            print(l)
    os.makedirs('_fixed_md', exist_ok=True)
    open(f'_fixed_md/{stem}.md', 'w', encoding='utf-8').write(new)
    return new, total_fixed

if __name__ == '__main__':
    for stem, _ in TARGETS:
        fix_file(stem)
