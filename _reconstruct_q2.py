"""Smarter ? reconstructor.

Key insight: between PRE (107906c7) and CURRENT, a cleanup pipeline:
  - replaced unicode \xD7 with \\times (wrapped in $...$)
  - replaced \cdot, \rightarrow, \nabla etc. with LaTeX commands
  - some bytes were normalized to literal Latin characters (ø = U+00F8 etc.)

So normalization for matching must STRIP these LaTeX wrappers from CURRENT.
"""
import os, re, subprocess
from collections import Counter

PRE_SHA = '107906c7'
TARGETS = [
 'PAPER_009b_Aether_String_TRZ_Damping_GW',
 'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF',
 'PAPER_016b_White_Dwarf_Foreground_UQFF',
 'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime',
 'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density',
]

# LaTeX commands that came in to REPLACE a unicode char that got corrupted
# Map of LaTeX -> probable original unicode (for context normalization)
LATEX_TO_UNI = {
    r'\times': '×', r'\cdot': '·', r'\rightarrow': '→', r'\to': '→',
    r'\nabla': '∇', r'\partial': '∂', r'\propto': '∝', r'\approx': '≈',
    r'\langle': '⟨', r'\rangle': '⟩', r'\div': '÷', r'\pm': '±',
    r'\sim': '~', r'\le': '≤', r'\leq': '≤', r'\ge': '≥', r'\geq': '≥',
    r'\neq': '≠', r'\equiv': '≡', r'\infty': '∞',
    r'\alpha': 'α', r'\beta': 'β', r'\gamma': 'γ', r'\Gamma': 'Γ',
    r'\delta': 'δ', r'\Delta': 'Δ', r'\epsilon': 'ε', r'\zeta': 'ζ',
    r'\eta': 'η', r'\theta': 'θ', r'\Theta': 'Θ', r'\kappa': 'κ',
    r'\lambda': 'λ', r'\Lambda': 'Λ', r'\mu': 'μ', r'\nu': 'ν',
    r'\xi': 'ξ', r'\pi': 'π', r'\Pi': 'Π', r'\rho': 'ρ', r'\sigma': 'σ',
    r'\Sigma': 'Σ', r'\tau': 'τ', r'\phi': 'φ', r'\Phi': 'Φ',
    r'\chi': 'χ', r'\psi': 'ψ', r'\Psi': 'Ψ', r'\omega': 'ω', r'\Omega': 'Ω',
}

def latex_to_unicode(s):
    """Replace $\\times$ and \\command with unicode equivalent."""
    # Strip $...$ math wrappers (keep inner)
    s = re.sub(r'\$([^$]*)\$', r'\1', s)
    # Replace LaTeX commands
    for cmd, uni in LATEX_TO_UNI.items():
        s = s.replace(cmd, uni)
    # Other cleanups: ø is sometimes ÷ or ·, leave for now
    return s

def normalize(s, keep_q=True):
    s = latex_to_unicode(s)
    s = re.sub(r'\s+', ' ', s).strip()
    if not keep_q:
        s = s.replace('?', '')
    return s

def fix_file(stem):
    path = f'whitepapers/{stem}.md'
    cur = open(path, encoding='utf-8').read()
    pre = subprocess.run(['git','show', f'{PRE_SHA}:{path}'],
                         capture_output=True).stdout.decode('utf-8', errors='replace')
    cur_lines = cur.split('\n')
    pre_lines = pre.split('\n')
    pre_norm_all = [normalize(l) for l in pre_lines]
    pre_full_norm = '\n'.join(pre_norm_all)

    fixed_lines = []
    fixed_count = 0
    unfixed_count = 0
    unfixed_log = []

    for li, line in enumerate(cur_lines):
        if '?' not in line:
            fixed_lines.append(line)
            continue
        new_line = line
        # Iterate ? positions right-to-left (to keep indices valid)
        for m in list(re.finditer(r'\?', line))[::-1]:
            pos = m.start()
            # Try expanding context windows
            replacement = None
            for win in (15, 10, 7, 5, 4):
                before = line[max(0,pos-win):pos]
                after = line[pos+1:pos+1+win]
                # Normalize each
                b_norm = latex_to_unicode(before)
                a_norm = latex_to_unicode(after)
                # Whitespace normalize
                b_norm = re.sub(r'\s+', ' ', b_norm)
                a_norm = re.sub(r'\s+', ' ', a_norm)
                if len(b_norm) < 2 and len(a_norm) < 2:
                    continue
                # Build regex
                pat = re.escape(b_norm) + r'(.{1,4}?)' + re.escape(a_norm)
                ms = re.findall(pat, pre_full_norm)
                if ms:
                    cnt = Counter(ms)
                    cand, freq = cnt.most_common(1)[0]
                    # Prefer non-? non-empty unique answers
                    if cand and '?' not in cand and freq >= 1:
                        # If multiple candidates, only use if dominant
                        if freq == sum(cnt.values()) or freq >= 2*max(c for k,c in cnt.items() if k!=cand) if len(cnt)>1 else True:
                            replacement = cand
                            break
            if replacement is not None:
                new_line = new_line[:pos] + replacement + new_line[pos+1:]
                fixed_count += 1
            else:
                unfixed_count += 1
                if len(unfixed_log) < 30:
                    b = line[max(0,pos-15):pos]
                    a = line[pos+1:pos+16]
                    unfixed_log.append(f'    L{li:4}: ...{b!r} ? {a!r}...')
        fixed_lines.append(new_line)

    new = '\n'.join(fixed_lines)
    print(f'\n=== {stem} ===')
    print(f'  Fixed: {fixed_count}  Unfixed: {unfixed_count}  Remaining ?: {new.count("?")}')
    for l in unfixed_log[:25]:
        print(l)
    os.makedirs('_fixed_md', exist_ok=True)
    open(f'_fixed_md/{stem}.md', 'w', encoding='utf-8').write(new)

if __name__ == '__main__':
    for s in TARGETS:
        fix_file(s)
