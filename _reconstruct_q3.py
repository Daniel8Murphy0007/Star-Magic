"""Per-line git blame walk to resolve remaining ? sites.

For each line in CURRENT containing ?:
  1. git blame -L <line>,<line> -- <path>  -> get blame commit C
  2. Walk commits that touched <path> from C backward; for each commit,
     get the line text in that commit (using same line number after rename
     via git log --follow not needed since path stable).
  3. Find a commit whose text at the same/nearby line has NO ? AND
     matches the context (before+after of the ?).
  4. Extract the substituted unicode char(s).
"""
import os, re, subprocess
from collections import Counter

TARGETS = [
 'PAPER_009b_Aether_String_TRZ_Damping_GW',
 'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF',
 'PAPER_016b_White_Dwarf_Foreground_UQFF',
 'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime',
 'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density',
]

LATEX_TO_UNI = {
    r'\times': '×', r'\cdot': '·', r'\rightarrow': '→', r'\to': '→',
    r'\nabla': '∇', r'\partial': '∂', r'\propto': '∝', r'\approx': '≈',
    r'\langle': '⟨', r'\rangle': '⟩', r'\div': '÷', r'\pm': '±',
    r'\le': '≤', r'\leq': '≤', r'\ge': '≥', r'\geq': '≥',
    r'\neq': '≠', r'\equiv': '≡', r'\infty': '∞',
    r'\alpha': 'α', r'\beta': 'β', r'\gamma': 'γ', r'\Gamma': 'Γ',
    r'\delta': 'δ', r'\Delta': 'Δ', r'\epsilon': 'ε', r'\eta': 'η',
    r'\theta': 'θ', r'\Theta': 'Θ', r'\kappa': 'κ', r'\lambda': 'λ',
    r'\Lambda': 'Λ', r'\mu': 'μ', r'\nu': 'ν', r'\xi': 'ξ', r'\pi': 'π',
    r'\rho': 'ρ', r'\sigma': 'σ', r'\Sigma': 'Σ', r'\tau': 'τ',
    r'\phi': 'φ', r'\Phi': 'Φ', r'\chi': 'χ', r'\psi': 'ψ',
    r'\omega': 'ω', r'\Omega': 'Ω',
}

def latex_to_unicode(s):
    s = re.sub(r'\$([^$]*)\$', r'\1', s)
    for cmd, uni in LATEX_TO_UNI.items():
        s = s.replace(cmd, uni)
    return s

def norm(s):
    return re.sub(r'\s+', ' ', latex_to_unicode(s)).strip()

def get_line_history(path, lineno):
    """git log -L <n>,<n>:<path> returns successive versions of the line."""
    r = subprocess.run(
        ['git','log','--no-patch','--format=%H','-L',
         f'{lineno},{lineno}:{path}'],
        capture_output=True, text=True, encoding='utf-8', errors='replace')
    # Parse hashes
    shas = [l.strip() for l in r.stdout.splitlines() if re.fullmatch(r'[0-9a-f]{40}', l.strip())]
    return shas

def file_at(sha, path):
    r = subprocess.run(['git','show', f'{sha}:{path}'], capture_output=True)
    return r.stdout.decode('utf-8', errors='replace')

def try_resolve(current_line, qpos, path, lineno):
    """Try to find a pre-corruption version of this exact site."""
    # Build context
    cur = current_line
    before_full = cur[:qpos]
    after_full = cur[qpos+1:]
    # Normalize for matching
    b_norm = norm(before_full)
    a_norm = norm(after_full)
    # Walk line history
    shas = get_line_history(path, lineno)
    for sha in shas:
        try:
            content = file_at(sha, path)
        except Exception:
            continue
        plines = content.split('\n')
        # search within +/-3 of lineno
        lo = max(0, lineno-4); hi = min(len(plines), lineno+3)
        for pl in plines[lo:hi]:
            pn = norm(pl)
            if '?' in pn:
                continue
            # Try suffix/prefix anchors
            for wb in (20, 12, 8, 5):
                for wa in (20, 12, 8, 5):
                    b = b_norm[-wb:] if len(b_norm) >= wb else b_norm
                    a = a_norm[:wa] if len(a_norm) >= wa else a_norm
                    if len(b) < 3 and len(a) < 3:
                        continue
                    pat = re.escape(b) + r'(.{1,4}?)' + re.escape(a)
                    m = re.search(pat, pn)
                    if m:
                        cand = m.group(1)
                        if cand and '?' not in cand:
                            return cand, sha
    return None, None

def fix_file(stem):
    path = f'whitepapers/{stem}.md'
    fixed_path = f'_fixed_md/{stem}.md'
    if os.path.exists(fixed_path):
        cur = open(fixed_path, encoding='utf-8').read()
    else:
        cur = open(path, encoding='utf-8').read()
    lines = cur.split('\n')
    fixed = 0
    unfixed = 0
    unfixed_log = []
    new_lines = []
    for li, line in enumerate(lines, start=1):
        if '?' not in line:
            new_lines.append(line); continue
        new_line = line
        for m in list(re.finditer(r'\?', line))[::-1]:
            pos = m.start()
            rep, sha = try_resolve(line, pos, path, li)
            if rep:
                new_line = new_line[:pos] + rep + new_line[pos+1:]
                fixed += 1
            else:
                unfixed += 1
                if len(unfixed_log) < 30:
                    b = line[max(0,pos-15):pos]; a = line[pos+1:pos+16]
                    unfixed_log.append(f'    L{li:4}: ...{b!r} ? {a!r}...')
        new_lines.append(new_line)
    new = '\n'.join(new_lines)
    print(f'\n=== {stem} ===  fixed={fixed} unfixed={unfixed} remaining?={new.count("?")}')
    for l in unfixed_log:
        print(l)
    os.makedirs('_fixed_md', exist_ok=True)
    open(fixed_path, 'w', encoding='utf-8').write(new)

if __name__ == '__main__':
    for s in TARGETS:
        fix_file(s)
