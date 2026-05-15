"""Resolve remaining ? sites using ONLY the CP1/CP2/CP3/CP4 pipeline as source.

The 100 unfixed sites are mostly equation fragments and calibrated numerical
values. CP1-4 contain the canonical Python implementations with proper unicode.
"""
import os, re, sys
from collections import Counter

CP_FILES = [
    'CondensedPhysics.py',
    'CondensedPhysics2.py',
    'CondensedPhysics3.py',
    'CondensedPhysics4.py',
    'CondensedPhysics_OutputData.py',
    'CondensedPhysics_InputData.py',
    'CondensedPhysics_Validation.py',
    'CondensedPhysicsAggregator.py',
]

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
    r'\leq': '≤', r'\le': '≤', r'\geq': '≥', r'\ge': '≥',
    r'\neq': '≠', r'\equiv': '≡', r'\infty': '∞',
    r'\alpha':'α',r'\beta':'β',r'\gamma':'γ',r'\Gamma':'Γ',
    r'\delta':'δ',r'\Delta':'Δ',r'\epsilon':'ε',r'\eta':'η',
    r'\theta':'θ',r'\Theta':'Θ',r'\kappa':'κ',r'\lambda':'λ',
    r'\Lambda':'Λ',r'\mu':'μ',r'\nu':'ν',r'\xi':'ξ',r'\pi':'π',
    r'\rho':'ρ',r'\sigma':'σ',r'\Sigma':'Σ',r'\tau':'τ',
    r'\phi':'φ',r'\Phi':'Φ',r'\chi':'χ',r'\psi':'ψ',
    r'\omega':'ω',r'\Omega':'Ω',
}
def latex_to_unicode(s):
    s = re.sub(r'\$([^$]*)\$', r'\1', s)
    for cmd, uni in sorted(LATEX_TO_UNI.items(), key=lambda x: -len(x[0])):
        s = s.replace(cmd, uni)
    return s
def norm(s):
    return re.sub(r'\s+', ' ', latex_to_unicode(s)).strip()

print('Loading CP pipeline corpus...', flush=True)
CORPUS = []
for f in CP_FILES:
    if not os.path.exists(f):
        print(f'  MISSING: {f}', flush=True)
        continue
    with open(f, encoding='utf-8', errors='replace') as fh:
        txt = fh.read()
    CORPUS.append((f, norm(txt)))
    print(f'  {f}: {len(txt):>10} chars  ?-count={txt.count("?")}', flush=True)
print(f'Corpus: {len(CORPUS)} files\n', flush=True)

def resolve_site(line, qpos):
    before = line[:qpos]; after = line[qpos+1:]
    b_n = norm(before); a_n = norm(after)
    # Strip cp437-mangled trailing junk from before: keep last printable chunk
    candidates = Counter()
    for wb in (40, 25, 15, 10, 6, 4):
        for wa in (40, 25, 15, 10, 6, 4):
            b = b_n[-wb:] if len(b_n) >= wb else b_n
            a = a_n[:wa] if len(a_n) >= wa else a_n
            if '?' in b or '?' in a:
                continue
            if len(b) < 3 and len(a) < 3:
                continue
            pat = re.escape(b) + r'(.{1,6}?)' + re.escape(a)
            try:
                cre = re.compile(pat)
            except re.error:
                continue
            for path, txt in CORPUS:
                for m in cre.finditer(txt):
                    cand = m.group(1)
                    if cand and '?' not in cand:
                        candidates[cand] += 1
            if candidates:
                break
        if candidates:
            break
    if not candidates:
        return None, None
    top, freq = candidates.most_common(1)[0]
    return top, freq

def fix_file(stem):
    path = f'whitepapers/{stem}.md'
    fixed_path = f'_fixed_md/{stem}.md'
    src = fixed_path if os.path.exists(fixed_path) else path
    cur = open(src, encoding='utf-8').read()
    lines = cur.split('\n')
    fixed = 0; unfixed = 0; sites = []
    new_lines = []
    for li, line in enumerate(lines, 1):
        if '?' not in line:
            new_lines.append(line); continue
        new_line = line
        for m in list(re.finditer(r'\?', line))[::-1]:
            pos = m.start()
            rep, freq = resolve_site(line, pos)
            b = line[max(0,pos-12):pos]; a = line[pos+1:pos+13]
            if rep:
                new_line = new_line[:pos] + rep + new_line[pos+1:]
                fixed += 1
                sites.append(f'  +L{li:4} q{pos:4} freq={freq:3} rep={rep!r:6}  ctx=...{b[-12:]!r}?{a[:12]!r}')
            else:
                unfixed += 1
                sites.append(f'  -L{li:4} q{pos:4} UNRESOLVED            ctx=...{b[-12:]!r}?{a[:12]!r}')
        new_lines.append(new_line)
    new = '\n'.join(new_lines)
    print(f'\n=== {stem} ===  fixed={fixed} unfixed={unfixed} remaining?={new.count("?")}', flush=True)
    for s in sites:
        print(s, flush=True)
    open(fixed_path,'w',encoding='utf-8').write(new)

if __name__ == '__main__':
    for s in TARGETS:
        fix_file(s)
