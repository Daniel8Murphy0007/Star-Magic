# fix_latex2.py - Fix pre-existing LaTeX errors in whitepapers/*.md
import re
from pathlib import Path

ROOT = Path(__file__).parent
WP_DIR = ROOT / "whitepapers"

GREEK = [
    'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta',
    'iota', 'kappa', 'lambda', 'mu', 'nu', 'xi', 'pi', 'rho', 'sigma',
    'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega',
    'Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Zeta', 'Eta', 'Theta',
    'Kappa', 'Lambda', 'Mu', 'Nu', 'Xi', 'Pi', 'Rho', 'Sigma',
    'Tau', 'Upsilon', 'Phi', 'Chi', 'Psi', 'Omega',
]

MATH_FUNCS = [
    'sin', 'cos', 'tan', 'cot', 'sec', 'csc', 'sinh', 'cosh', 'tanh',
    'arcsin', 'arccos', 'arctan', 'ln', 'log', 'exp', 'min', 'max',
    'lim', 'det', 'gcd', 'ker', 'arg', 'deg', 'dim', 'to', 'left', 'right',
    'frac', 'sqrt', 'cdot', 'ell',
]

UPPERCASE_CMDS = [
    'Delta', 'Sigma', 'Gamma', 'Lambda', 'Omega', 'Theta', 'Pi',
    'Phi', 'Psi', 'Xi', 'Upsilon', 'Nabla',
]

LOWERCASE_CMDS = [
    'delta', 'sigma', 'gamma', 'lambda', 'omega', 'theta', 'pi',
    'phi', 'psi', 'xi', 'upsilon', 'alpha', 'beta', 'kappa',
    'mu', 'nu', 'rho', 'tau', 'epsilon', 'eta', 'iota', 'zeta',
    'nabla',
]

ALL_CMD_LOWER = set(g.lower() for g in GREEK) | set(f.lower() for f in MATH_FUNCS)


def fix_content(text):
    original = text

    # 0. \rm → \mathrm (deprecated font switch in math mode, must come first)
    # {\rm word} → {\mathrm{word}}, {\rm word,word} → {\mathrm{word,word}}
    text = re.sub(r'\{\\rm\b([^}]*)\}',
                  lambda m: '{\\mathrm{' + m.group(1).strip() + '}}', text)
    # \rm word (not in braces) → \mathrm{word}
    text = re.sub(r'\\rm\s+(\w[\w,\-]*)',
                  lambda m: '\\mathrm{' + m.group(1) + '}', text)

    # 0b. \left\arrow → just \arrow (\left must precede a bracket/delimiter, not arrows)
    text = re.sub(r'\\left(\\(?:rightarrow|leftarrow|leftrightarrow|Rightarrow|Leftarrow|Leftrightarrow|to|gets|mapsto))',
                  lambda m: m.group(1), text)

    # 1. cdotBig/Bigl/Bigr/Bigg etc. - missing backslash
    text = re.sub(r'\\cdot(Bigg[lr]?|Bigg?|Big[lr]?)',
                  lambda m: '\\cdot\\' + m.group(1), text)

    # 2. cdot/times followed by any word
    BIG_CMDS = {'Bigl', 'Bigr', 'Bigg', 'Biggl', 'Biggr', 'Big',
                'Left', 'Right', 'Frac', 'Sqrt', 'Text', 'Mathrm'}

    def fix_cdot(m):
        following = m.group(1)
        if following.lower() in ALL_CMD_LOWER or following in BIG_CMDS:
            return '\\cdot\\' + following
        return '\\cdot ' + following
    text = re.sub(r'\\cdot([A-Za-z]+)', fix_cdot, text)

    def fix_times(m):
        following = m.group(1)
        if following.lower() in ALL_CMD_LOWER or following in BIG_CMDS:
            return '\\times\\' + following
        return '\\times ' + following
    text = re.sub(r'\\times([A-Za-z]+)', fix_times, text)

    # 3. Math func + Greek without backslash: sintheta, cosomega, lnsigma etc.
    # Use (?![a-zA-Z]) instead of \b so we match even when followed by _ digit etc.
    math_func_names = [
        'sin', 'cos', 'tan', 'cot', 'sec', 'csc', 'sinh', 'cosh', 'tanh',
        'arcsin', 'arccos', 'arctan', 'ln', 'log', 'exp', 'min', 'max',
        'lim', 'det', 'gcd', 'ker', 'ell',
    ]
    for func in math_func_names:
        for greek in GREEK:
            text = re.sub('\\\\' + func + greek + r'(?![a-zA-Z])',
                          lambda m, f=func, g=greek: '\\' + f + '\\' + g, text)

    # 4. ell+greek: ellnu -> ell\nu
    for greek in GREEK:
        text = re.sub('\\\\ell' + greek + r'(?![a-zA-Z])',
                      lambda m, g=greek: '\\ell\\' + g, text)

    # 5. to+units: tofb -> to\,\text{fb}
    text = re.sub(r'\\to(fb|pb|nb|ub|ab)\b',
                  lambda m: '\\to\\,\\text{' + m.group(1) + '}', text)
    text = re.sub(r'\\min(left)\b', r'\\min\\left', text)

    # 6. Uppercase cmd + uppercase: DeltaX -> Delta X or Delta\X if known cmd
    def make_upper_repl(c):
        def r(m):
            f = m.group(1)
            if f.lower() in ALL_CMD_LOWER:
                return '\\' + c + '\\' + f
            return '\\' + c + ' ' + f
        return r
    for cmd in UPPERCASE_CMDS:
        text = re.sub(r'\\' + cmd + r'([A-Z][a-zA-Z]*)', make_upper_repl(cmd), text)

    # 7. Lowercase cmd + uppercase: deltaT -> delta T or delta\T if known
    def make_lower_repl(c):
        def r(m):
            f = m.group(1)
            if f.lower() in ALL_CMD_LOWER:
                return '\\' + c + '\\' + f
            return '\\' + c + ' ' + f
        return r
    for cmd in LOWERCASE_CMDS:
        text = re.sub(r'\\' + cmd + r'([A-Z][a-zA-Z]*)', make_lower_repl(cmd), text)

    # 8. Cmd + math func: Deltalog -> Delta\log, deltalog -> delta\log
    all_cmds = UPPERCASE_CMDS + LOWERCASE_CMDS
    extra_funcs = [
        'log', 'ln', 'sin', 'cos', 'tan', 'exp', 'max', 'min', 'det',
        'lim', 'sup', 'inf', 'arg', 'deg', 'dim', 'gcd',
        'sec', 'csc', 'cot', 'sinh', 'cosh', 'tanh',
        'arcsin', 'arccos', 'arctan',
    ]
    for cmd in all_cmds:
        for func in extra_funcs:
            text = re.sub('\\\\' + cmd + func + r'(?![a-zA-Z])',
                          lambda m, c=cmd, f=func: '\\' + c + '\\' + f, text)

    # 9. word_$\Cmd$ -> $word_{\Cmd}$ (subscript before math span)
    # e.g. a_$\Lambda$ -> $a_{\Lambda}$, F_U_$\text{Bi}$ -> $F_{U_{\text{Bi}}}$
    def fix_sub(m):
        word, sub, inner = m.group(1), m.group(2), m.group(3)
        inner_stripped = inner[1:-1]  # strip outer dollar signs
        return '$' + word + sub + '{' + inner_stripped + '}$'
    text = re.sub(r'([A-Za-z][A-Za-z0-9]*)([_\^])(\$\\[A-Za-z][^$\n]{0,50}\$)',
                  fix_sub, text)

    # 10. Unicode Greek/math chars → LaTeX (chars missed by upgrade_latex.py)
    UNICODE_TO_LATEX = {
        # C1 control chars stored as W1252 surrogates (causes pandoc YAML failure)
        '\u0085': r'\ldots',   # U+0085 (NEXT LINE) used as ellipsis
        '\u0096': '\u2013',    # U+0096 → en-dash (safe for pandoc YAML)
        '\u0097': '\u2014',    # U+0097 → em-dash (safe for pandoc YAML)
        '\u0098': r'\approx',  # U+0098 (SPA) used as ≈
        # Greek lowercase
        'μ': r'\mu', 'µ': r'\mu', 'η': r'\eta', 'α': r'\alpha', 'β': r'\beta',
        'γ': r'\gamma', 'δ': r'\delta', 'ε': r'\epsilon', 'ζ': r'\zeta',
        'θ': r'\theta', 'ι': r'\iota', 'κ': r'\kappa', 'λ': r'\lambda',
        'ν': r'\nu', 'ξ': r'\xi', 'π': r'\pi', 'ρ': r'\rho', 'σ': r'\sigma',
        'τ': r'\tau', 'υ': r'\upsilon', 'φ': r'\phi', 'χ': r'\chi', 'ψ': r'\psi',
        'ω': r'\omega', 'Γ': r'\Gamma', 'Δ': r'\Delta', 'Θ': r'\Theta',
        'Λ': r'\Lambda', 'Ξ': r'\Xi', 'Π': r'\Pi', 'Σ': r'\Sigma',
        'Υ': r'\Upsilon', 'Φ': r'\Phi', 'Ψ': r'\Psi', 'Ω': r'\Omega',
        '∑': r'\sum', '∏': r'\prod', '∫': r'\int', '∂': r'\partial',
        '∇': r'\nabla', '∞': r'\infty', '≈': r'\approx', '≤': r'\leq',
        '≥': r'\geq', '≠': r'\neq', '±': r'\pm', '×': r'\times', '·': r'\cdot',
    }
    for uc, latex in UNICODE_TO_LATEX.items():
        text = text.replace(uc, latex)

    # 11. log10 → \log_{10} (both plain and after \)
    text = re.sub(r'\\log10\b', r'\\log_{10}', text)           # \log10 → \log_{10}
    text = re.sub(r'(?<![a-zA-Z\\])log10\b', r'\\log_{10}', text)  # log10 → \log_{10}

    # 11. Double/triple subscripts x_y_z → x_{y\_z} (fix TeX "double subscript" error)
    # Apply globally - these patterns only appear as math-mode identifiers in these files
    def fix_multi_sub(m):
        full = m.group(0)
        base = m.group(1)
        subs = re.findall(r'_([A-Za-z0-9]+)', full[1:])
        if len(subs) < 2:
            return full
        return base + '_{' + '\\_'.join(subs) + '}'
    # Match: base char + 2 or more _word parts (not already inside braces)
    text = re.sub(r'(?<![{\\])([A-Za-z0-9])(?:_[A-Za-z0-9]+){2,}', fix_multi_sub, text)

    n_changed = sum(1 for a, b in zip(original.split('\n'), text.split('\n')) if a != b)
    return text, n_changed


def process_files():
    md_files = sorted(WP_DIR.glob("PAPER_*.md"))
    fixed, changed_list = 0, []

    for i, md_path in enumerate(md_files):
        try:
            content = md_path.read_text(encoding='utf-8')
        except UnicodeDecodeError:
            try:
                content = md_path.read_text(encoding='latin-1')
            except Exception:
                print(f"  SKIP: {md_path.name}")
                continue

        new_content, changes = fix_content(content)
        if new_content != content:
            md_path.write_text(new_content, encoding='utf-8')
            fixed += 1
            changed_list.append(str(md_path))

        if (i + 1) % 300 == 0:
            print(f"  Progress: {i+1}/{len(md_files)} ({fixed} fixed)")

    print(f"\nDone: {fixed} files fixed")
    (ROOT / "_fixed_files2.txt").write_text('\n'.join(changed_list), encoding='utf-8')
    return changed_list


if __name__ == '__main__':
    print(f"Scanning {WP_DIR}...")
    process_files()
