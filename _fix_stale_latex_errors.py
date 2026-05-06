#!/usr/bin/env python3
"""
_fix_stale_latex_errors.py
Auto-fix common LaTeX errors in the 140 whitepaper .md files that failed
during the stale-PDF regeneration run.

Error categories addressed:
  1. Run-together LaTeX commands (\argmin, \quadRightarrow, \leftlangle, etc.)
  2. Unicode characters unsupported by pdflatex (✅ ħ ├ ₂ − ϵ ϕ ...)
  3. Double subscript/superscript patterns
  4. Bare # in math/table contexts
  5. \right without matching \left (specific patterns)
  6. Missing math delimiters (specific patterns)
  7. pdflatex memory overflow: add font pool opt to generate_pdfs.py
"""
import re, os, glob, sys

WHITEPAPER_DIR = "whitepapers"

# All 140 paper numbers that failed in _regen_stale_240.py
FAILED = [
    27,28,30,31,32,33,35,36,51,63,90,106,138,143,144,155,160,163,164,167,
    181,184,188,198,202,210,214,215,224,239,240,242,259,261,262,265,267,268,
    278,298,313,336,351,354,372,373,375,380,381,384,386,389,429,435,439,
    461,462,464,473,491,494,498,513,514,526,532,533,535,536,544,545,549,
    554,557,563,570,573,574,575,576,577,578,581,582,583,585,587,590,592,
    598,600,633,645,646,647,650,651,653,688,692,701,716,717,718,721,722,
    731,732,735,738,739,740,747,749,794,798,803,807,808,812,831,832,833,
    840,865,877,880,882,883,888,890,904,949,953,957,980,1023,1025,1040,1101,
]

# ──────────────────────────────────────────────────────────────────
# REGEX FIXES  (pattern, replacement, label)
# Applied to the full text of each .md file.
# Order matters: more specific patterns first.
# ──────────────────────────────────────────────────────────────────
REGEX_FIXES = [
    # ── 1. \argmin, \argmax ──────────────────────────────────────
    (r'\\argmin\b', r'\\operatorname{argmin}',  r'\argmin'),
    (r'\\argmax\b', r'\\operatorname{argmax}',  r'\argmax'),

    # ── 2. \quad + command name (missing backslash on 2nd cmd) ───
    (r'\\quadRightarrow\b',      r'\\quad\\Rightarrow',       r'\quadRightarrow'),
    (r'\\quadrightarrow\b',      r'\\quad\\rightarrow',       r'\quadrightarrow'),
    (r'\\quadLeftarrow\b',       r'\\quad\\Leftarrow',        r'\quadLeftarrow'),
    (r'\\quadleftrightarrow\b',  r'\\quad\\leftrightarrow',   r'\quadleftrightarrow'),
    (r'\\quadforall\b',          r'\\quad\\forall',           r'\quadforall'),
    (r'\\quadcheckmark\b',       r'\\quad\\checkmark',        r'\quadcheckmark'),
    (r'\\quadapprox\b',          r'\\quad\\approx',           r'\quadapprox'),
    (r'\\quadimplies\b',         r'\\quad\\implies',          r'\quadimplies'),
    (r'\\quadpropto\b',          r'\\quad\\propto',           r'\quadpropto'),
    (r'\\quadleq\b',             r'\\quad\\leq',              r'\quadleq'),
    (r'\\quadgeq\b',             r'\\quad\\geq',              r'\quadgeq'),
    (r'\\quadneq\b',             r'\\quad\\neq',              r'\quadneq'),
    (r'\\quadnabla\b',           r'\\quad\\nabla',            r'\quadnabla'),
    (r'\\quadsum\b',             r'\\quad\\sum',              r'\quadsum'),
    (r'\\quadtext\{',            r'\\quad\\text{',            r'\quadtext{'),
    (r'\\quadtext\b',            r'\\quad\\text',             r'\quadtext'),
    (r'\\quadcdot\b',            r'\\quad\\cdot',             r'\quadcdot'),
    (r'\\quadtimes\b',           r'\\quad\\times',            r'\quadtimes'),
    (r'\\quadin\b',              r'\\quad\\in',               r'\quadin'),

    # ── 3. \pm + letter/command (run-together) ───────────────────
    (r'\\pmsqrt\b',    r'\\pm\\sqrt',    r'\pmsqrt'),
    (r'\\pmkappa\b',   r'\\pm\\kappa',   r'\pmkappa'),
    (r'\\pmalpha\b',   r'\\pm\\alpha',   r'\pmalpha'),
    (r'\\pmbeta\b',    r'\\pm\\beta',    r'\pmbeta'),
    (r'\\pmgamma\b',   r'\\pm\\gamma',   r'\pmgamma'),
    (r'\\pmDelta\b',   r'\\pm\\Delta',   r'\pmDelta'),
    (r'\\pmdelta\b',   r'\\pm\\delta',   r'\pmdelta'),
    (r'\\pmlambda\b',  r'\\pm\\lambda',  r'\pmlambda'),
    (r'\\pmsigma\b',   r'\\pm\\sigma',   r'\pmsigma'),
    (r'\\pmF\b',       r'\\pm F',        r'\pmF'),
    (r'\\pmG\b',       r'\\pm G',        r'\pmG'),
    (r'\\pmE\b',       r'\\pm E',        r'\pmE'),

    # ── 4. \left + brace (run-together) ──────────────────────────
    (r'\\leftlangle\b',  r'\\left\\langle',  r'\leftlangle'),
    (r'\\leftlfloor\b',  r'\\left\\lfloor',  r'\leftlfloor'),
    (r'\\leftlceil\b',   r'\\left\\lceil',   r'\leftlceil'),
    (r'\\leftvert\b',    r'\\left|',          r'\leftvert'),

    # ── 5. Math operators run-together ────────────────────────────
    (r'\\arctanleft\b',   r'\\arctan\\left',   r'\arctanleft'),
    (r'\\arctanfrac\b',   r'\\arctan\\frac',   r'\arctanfrac'),
    (r'\\minbigl\b',      r'\\min\\bigl',      r'\minbigl'),
    (r'\\maxbigl\b',      r'\\max\\bigl',      r'\maxbigl'),
    (r'\\coshfrac\b',     r'\\cosh\\frac',     r'\coshfrac'),
    (r'\\sinhfrac\b',     r'\\sinh\\frac',     r'\sinhfrac'),
    (r'\\tanhfrac\b',     r'\\tanh\\frac',     r'\tanhfrac'),
    (r'\\lnbigl\b',       r'\\ln\\bigl',       r'\lnbigl'),
    (r'\\lnfrac\b',       r'\\ln\\frac',       r'\lnfrac'),
    (r'\\expfrac\b',      r'\\exp\\frac',      r'\expfrac'),
    (r'\\coscos\b',       r'\\cos',            r'\coscos'),
    (r'\\sinsin\b',       r'\\sin',            r'\sinsin'),

    # ── 6. \delta/\Delta/\nabla + \mathcal (run-together) ────────
    (r'\\deltamathcal\b',  r'\\delta\\mathcal',  r'\deltamathcal'),
    (r'\\Deltamathcal\b',  r'\\Delta\\mathcal',  r'\Deltamathcal'),
    (r'\\nablamathcal\b',  r'\\nabla\\mathcal',  r'\nablamathcal'),
    (r'\\partialmathcal\b',r'\\partial\\mathcal',r'\partialmathcal'),

    # ── 7. Miscellaneous run-togethers ────────────────────────────
    (r'\\approxhbar\b',   r'\\approx\\hbar',   r'\approxhbar'),
    (r'\\approxfrac\b',   r'\\approx\\frac',   r'\approxfrac'),
    (r'\\ldotsright\b',   r'\\ldots\\right',   r'\ldotsright'),
    (r'\\cdotsright\b',   r'\\cdots\\right',   r'\cdotsright'),
    (r'\\sigmasim\b',     r'\\sigma\\sim',      r'\sigmasim'),
    (r'\\ggtau\b',        r'\\gg\\tau',         r'\ggtau'),
    (r'\\lltau\b',        r'\\ll\\tau',         r'\lltau'),
    (r'\\ggalpha\b',      r'\\gg\\alpha',       r'\ggalpha'),
    (r'\\nurangle\b',     r'\\nu\\rangle',      r'\nurangle'),
    (r'\\murangle\b',     r'\\mu\\rangle',      r'\murangle'),
    (r'\\Lambdac\b',      r'\\Lambda c',        r'\Lambdac'),
    (r'\\intpsi\b',       r'\\int\\psi',        r'\intpsi'),
    (r'\\intphi\b',       r'\\int\\phi',        r'\intphi'),
    (r'\\partialL\b',     r'\\partial L',       r'\partialL'),
    (r'\\partialH\b',     r'\\partial H',       r'\partialH'),
    (r'\\toinfty\b',      r'\\to\\infty',       r'\toinfty'),
    (r'\\toinf\b',        r'\\to\\infty',       r'\toinf'),
    (r'\\pir\b',          r'\\pi r',            r'\pir'),
    (r'\\piR\b',          r'\\pi R',            r'\piR'),
    (r'\\kappat\b',       r'\\kappa t',         r'\kappat'),
    (r'\\alphat\b',       r'\\alpha t',         r'\alphat'),
    (r'\\Deltak\b',       r'\\Delta k',         r'\Deltak'),
    (r'\\DeltaE\b',       r'\\Delta E',         r'\DeltaE'),
    (r'\\DeltaV\b',       r'\\Delta V',         r'\DeltaV'),
    (r'\\cdota\b',        r'\\cdot a',          r'\cdota'),
    (r'\\mum\b',          r'\\mu\\,\\mathrm{m}',r'\mum'),
    (r'\\muG\b',          r'\\mu\\mathrm{G}',   r'\muG'),
    (r'\\muT\b',          r'\\mu\\mathrm{T}',   r'\muT'),
    (r'\\toe\b',          r'\\to e',            r'\toe'),
    (r'\\tomu\b',         r'\\to\\mu',          r'\tomu'),
    (r'\\Ntoinfty\b',     r'N\\to\\infty',      r'\Ntoinfty'),
    (r'\\ntoinfty\b',     r'n\\to\\infty',      r'\ntoinfty'),
    (r'\\coshfrac\b',     r'\\cosh\\frac',      r'\coshfrac2'),
    (r'\\xleftrightarrow(?!\{)', r'\\xleftrightarrow{}', r'\xleftrightarrow no-arg'),

    # ── 8. Double subscript / superscript fixes ───────────────────
    # -^1^  pattern (literal text, not in math)
    (r'(\d)-\^\^(\d+)',   r'\\1\\times 10^{-\\2}',  'sci notation -^^'),
    (r'(\d)x10-\^1\^',   r'\\1\\times 10^{-1}',    'sci notation x10-^1^'),
    (r'(\d)x10-\^(\d+)\^', r'\\1\\times 10^{-\\2}','sci notation x10-^n^'),

    # ── 9. \right without \left (specific patterns) ───────────────
    # Extra \right) when there's a spurious close — only safe patterns
    (r'\)\\right\)', r'\\right)',  r')\right) -> \right)'),

    # ── 10. Bare # in math/alignment (table rows ending with #) ──
    # In tabular/align, # at end of line after content is invalid
    # Replace trailing # on table/align lines with \#
    # But be careful: only in \begin{align}/tabular context is this an issue
    # Safe: escape # that follows whitespace+number+units in table cell
    (r'(\s)#(\s*$)', r'\\1\\#\\2', r'trailing # in table'),
]

# ──────────────────────────────────────────────────────────────────
# UNICODE FIXES  (char, LaTeX replacement)
# Applied as plain string replacement after regex fixes.
# For chars that appear both in math and text, we use the LaTeX command
# (pandoc wraps them correctly when they appear in math blocks).
# ──────────────────────────────────────────────────────────────────
UNICODE_FIXES = [
    # Emoji / box-drawing that pdflatex cannot handle at all
    ('✅', r'\checkmark'),
    ('✓', r'\checkmark'),
    ('✗', r'\times'),
    ('├', '|'),
    ('│', '|'),
    ('┤', '|'),
    ('─', '--'),
    ('═', '=='),
    # Unicode minus vs hyphen-minus
    ('\u2212', '-'),   # U+2212 MINUS SIGN
    ('\u2013', '--'),  # U+2013 EN DASH
    ('\u2014', '---'), # U+2014 EM DASH
    # Subscript digits (not in pdflatex font)
    ('\u2080', '$_0$'),
    ('\u2081', '$_1$'),
    ('\u2082', '$_2$'),
    ('\u2083', '$_3$'),
    ('\u2084', '$_4$'),
    ('\u2085', '$_5$'),
    ('\u2086', '$_6$'),
    ('\u2087', '$_7$'),
    ('\u2088', '$_8$'),
    ('\u2089', '$_9$'),
    # Superscript digits
    ('\u00B9', '$^1$'),
    ('\u00B2', '$^2$'),
    ('\u00B3', '$^3$'),
    # Special letters not in OT1/T1
    ('\u0127', r'\hbar'),    # ħ Latin small letter h with stroke
    ('\u03F5', r'\epsilon'), # ϵ Greek lunate epsilon
    ('\u03D5', r'\phi'),     # ϕ Greek phi symbol
    # Arrows
    ('\u2192', r'\rightarrow'),
    ('\u2190', r'\leftarrow'),
    ('\u2194', r'\leftrightarrow'),
    ('\u21D2', r'\Rightarrow'),
    ('\u21D0', r'\Leftarrow'),
    ('\u21D4', r'\Leftrightarrow'),
    # Math symbols
    ('\u221E', r'\infty'),
    ('\u2248', r'\approx'),
    ('\u2260', r'\neq'),
    ('\u2265', r'\geq'),
    ('\u2264', r'\leq'),
    ('\u00D7', r'\times'),
    ('\u00B7', r'\cdot'),
    ('\u2207', r'\nabla'),
    ('\u2202', r'\partial'),
    ('\u222B', r'\int'),
    ('\u221A', r'\sqrt{}'),
    ('\u03B1', r'\alpha'),
    ('\u03B2', r'\beta'),
    ('\u03B3', r'\gamma'),
    ('\u03B4', r'\delta'),
    ('\u03B5', r'\varepsilon'),
    ('\u03B7', r'\eta'),
    ('\u03B8', r'\theta'),
    ('\u03BA', r'\kappa'),
    ('\u03BB', r'\lambda'),
    ('\u03BC', r'\mu'),
    ('\u03BD', r'\nu'),
    ('\u03BE', r'\xi'),
    ('\u03C0', r'\pi'),
    ('\u03C1', r'\rho'),
    ('\u03C3', r'\sigma'),
    ('\u03C4', r'\tau'),
    ('\u03C5', r'\upsilon'),
    ('\u03C6', r'\varphi'),
    ('\u03C7', r'\chi'),
    ('\u03C8', r'\psi'),
    ('\u03C9', r'\omega'),
    ('\u0393', r'\Gamma'),
    ('\u0394', r'\Delta'),
    ('\u0398', r'\Theta'),
    ('\u039B', r'\Lambda'),
    ('\u039E', r'\Xi'),
    ('\u03A0', r'\Pi'),
    ('\u03A3', r'\Sigma'),
    ('\u03A5', r'\Upsilon'),
    ('\u03A6', r'\Phi'),
    ('\u03A8', r'\Psi'),
    ('\u03A9', r'\Omega'),
]

def read_md(path):
    for enc in ('utf-8', 'utf-8-sig', 'cp1252', 'latin-1'):
        try:
            with open(path, encoding=enc) as f:
                return f.read(), enc
        except (UnicodeDecodeError, LookupError):
            pass
    with open(path, encoding='utf-8', errors='replace') as f:
        return f.read(), 'utf-8-replace'

def apply_fixes(text):
    changes = []
    for pat, repl, label in REGEX_FIXES:
        new_text, n = re.subn(pat, repl, text, flags=re.MULTILINE)
        if n:
            changes.append(f'{label} x{n}')
            text = new_text
    for char, repl in UNICODE_FIXES:
        if char in text:
            count = text.count(char)
            text = text.replace(char, repl)
            changes.append(f'U+{ord(char):04X} ({char!r}) x{count}')
    return text, changes

def find_md_for_num(n):
    """Find whitepaper .md file for paper number n (handles 3- and 4-digit padding)."""
    for pattern in [
        f'PAPER_{n:04d}_*.md',
        f'PAPER_{n:03d}_*.md',
        f'PAPER_{n:04d}*.md',
        f'PAPER_{n:03d}*.md',
    ]:
        matches = glob.glob(os.path.join(WHITEPAPER_DIR, pattern))
        if matches:
            return matches[0]
    return None

def main():
    targets = FAILED if '--all' not in sys.argv else None

    if targets is None:
        # Apply to all whitepapers
        all_md = sorted(glob.glob(os.path.join(WHITEPAPER_DIR, 'PAPER_*.md')))
    else:
        all_md = []
        missing = []
        for n in targets:
            p = find_md_for_num(n)
            if p:
                all_md.append(p)
            else:
                missing.append(n)
        if missing:
            print(f'WARNING: {len(missing)} paper .md files not found: {missing[:10]}...')

    print(f'Applying LaTeX fixes to {len(all_md)} files...\n')
    fixed_count = 0
    unchanged_count = 0
    for md_path in all_md:
        text, enc = read_md(md_path)
        new_text, changes = apply_fixes(text)
        if new_text != text:
            # Write back in UTF-8 always (safer for pdflatex pipeline)
            with open(md_path, 'w', encoding='utf-8') as f:
                f.write(new_text)
            fixed_count += 1
            fname = os.path.basename(md_path)
            print(f'  FIXED  {fname}')
            for c in changes:
                print(f'         {c}')
        else:
            unchanged_count += 1

    print(f'\n{"="*55}')
    print(f'  Files fixed   : {fixed_count}')
    print(f'  Unchanged     : {unchanged_count}')

if __name__ == '__main__':
    main()
