"""Resolve remaining ? sites by searching the ENTIRE codebase for the surrounding
context. The same equations (e.g. ∂N/∂t = ∇(D∇N) - ∇·(v N) + Q) and numerical
values (4.1499 × 10²¹) appear in many other files (whitepapers, .cpp, .py, .js,
.tex, .md) with proper unicode.

For each unfixed ? site:
  1. Extract before/after window (LaTeX-normalized to unicode equivalents).
  2. Grep the whole codebase (excluding the corrupted files and pdf/ binaries).
  3. Find files whose normalized content contains the same context with a
     non-? char between.
  4. Tally candidates; pick the dominant one.
"""
import os, re, glob
from collections import Counter

TARGETS = [
 'PAPER_009b_Aether_String_TRZ_Damping_GW',
 'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF',
 'PAPER_016b_White_Dwarf_Foreground_UQFF',
 'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime',
 'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density',
]
TARGET_PATHS = {f'whitepapers/{t}.md' for t in TARGETS}

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
    # Superscript LaTeX patterns
    r'^{-1}': '⁻¹', r'^{-2}': '⁻²', r'^{-3}': '⁻³',
    r'^{-5}': '⁻⁵', r'^{-6}': '⁻⁶', r'^{-7}': '⁻⁷', r'^{-8}': '⁻⁸',
    r'^{-9}': '⁻⁹', r'^{-10}': '⁻¹⁰',
}

def latex_to_unicode(s):
    s = re.sub(r'\$([^$]*)\$', r'\1', s)
    for cmd, uni in sorted(LATEX_TO_UNI.items(), key=lambda x: -len(x[0])):
        s = s.replace(cmd, uni)
    return s

def norm(s):
    return re.sub(r'\s+', ' ', latex_to_unicode(s)).strip()

# Build corpus: all text files, EXCLUDE the 5 corrupted files
EXTS = ('.md','.txt','.cpp','.h','.py','.js','.tex','.json','.csv')
EXCLUDE_DIRS = {'.git','pdf','pdf_backup_pandoc_2026-05-12','node_modules',
                'build','build_msvc','__pycache__','dist','out','_fixed_md',
                '_orig_md','pdf_backup','.vs','.vscode'}

def load_corpus():
    import sys
    corpus = []
    skipped_corrupt = 0
    print('Loading corpus...', flush=True)
    for root, dirs, files in os.walk('.'):
        dirs[:] = [d for d in dirs if d not in EXCLUDE_DIRS and not d.startswith('.')]
        for f in files:
            if not f.endswith(EXTS):
                continue
            p = os.path.relpath(os.path.join(root, f)).replace('\\','/')
            if p in TARGET_PATHS:
                continue
            if p.startswith('_') or '/backup' in p.lower():
                continue
            try:
                with open(p, encoding='utf-8', errors='replace') as fh:
                    txt = fh.read()
            except Exception:
                continue
            if '?' in txt and txt.count('?') > 50:
                skipped_corrupt += 1
                continue
            corpus.append((p, norm(txt)))
            if len(corpus) % 200 == 0:
                print(f'  loaded {len(corpus)} files...', flush=True)
    print(f'  Corpus loaded: {len(corpus)} files (skipped {skipped_corrupt} as ?-corrupt)', flush=True)
    return corpus

CORPUS = None

def resolve_site(line, qpos):
    global CORPUS
    if CORPUS is None:
        CORPUS = load_corpus()
    before = line[:qpos]
    after = line[qpos+1:]
    b_n = norm(before)
    a_n = norm(after)
    # Strip mangled cp437 noise from the windows: drop chars outside printable ASCII+common math
    # Actually keep all - just use shorter windows from clean tail
    candidates = Counter()
    for wb in (25, 18, 12, 8, 5):
        for wa in (25, 18, 12, 8, 5):
            b = b_n[-wb:] if len(b_n) >= wb else b_n
            a = a_n[:wa] if len(a_n) >= wa else a_n
            if len(b) < 3 and len(a) < 3:
                continue
            # Skip if b/a contain '?' themselves (recursive)
            if '?' in b or '?' in a:
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
        return None
    top, freq = candidates.most_common(1)[0]
    # Require dominance: top must be >= 2x runner-up OR unique
    if len(candidates) == 1 or freq >= 2 * sum(c for k,c in candidates.items() if k != top):
        return top
    # Otherwise pick by length=1 preference (single char most likely)
    short = [(k,v) for k,v in candidates.items() if len(k)==1]
    if short:
        short.sort(key=lambda x:-x[1])
        if short[0][1] >= 2:
            return short[0][0]
    return top  # best effort

def fix_file(stem):
    path = f'whitepapers/{stem}.md'
    fixed_path = f'_fixed_md/{stem}.md'
    src = fixed_path if os.path.exists(fixed_path) else path
    cur = open(src, encoding='utf-8').read()
    lines = cur.split('\n')
    fixed = 0; unfixed = 0; log = []
    new_lines = []
    for li, line in enumerate(lines, 1):
        if '?' not in line:
            new_lines.append(line); continue
        new_line = line
        for m in list(re.finditer(r'\?', line))[::-1]:
            pos = m.start()
            rep = resolve_site(line, pos)
            if rep:
                new_line = new_line[:pos] + rep + new_line[pos+1:]
                fixed += 1
                if len(log) < 50:
                    b = line[max(0,pos-15):pos]; a = line[pos+1:pos+16]
                    log.append(f'  L{li:4} +"{rep}"  ctx=...{b[-15:]!r}?{a[:15]!r}')
            else:
                unfixed += 1
        new_lines.append(new_line)
    new = '\n'.join(new_lines)
    print(f'\n=== {stem} ===  fixed={fixed} unfixed={unfixed} remaining?={new.count("?")}')
    for l in log[:30]:
        print(l)
    open(fixed_path,'w',encoding='utf-8').write(new)

if __name__ == '__main__':
    for s in TARGETS:
        fix_file(s)
