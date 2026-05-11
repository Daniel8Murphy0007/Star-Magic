"""
Targeted shared-typo fix across whitepapers.
Replaces the exact buggy boilerplate strings that appear identically in
many files. NO builds are run here.
"""
import pathlib, re

WP = pathlib.Path("whitepapers")

# Literal replacements (exact-match substring)
LITERAL = [
    (r"\tanh!\left(\frac{t - t_{\mathrm{sat}}}{\tau _{\mathrm{BSH}}}\right\right)",
     r"\tanh\!\left(\frac{t - t_{\mathrm{sat}}}{\tau_{\mathrm{BSH}}}\right)\right)"),
    (r"\tanh!\left(\frac{t - t_{\mathrm{sat}}}{\tau_{\mathrm{BSH}}}\right\right)",
     r"\tanh\!\left(\frac{t - t_{\mathrm{sat}}}{\tau_{\mathrm{BSH}}}\right)\right)"),
]

# Regex replacements
REGEX = [
    # \cos!\left -> \cos\!\left (function-name followed by ! and \left/\Bigl/\bigl/\Big/\big)
    (re.compile(r"\\(cos|sin|tan|tanh|sinh|cosh|exp|ln|log|sec|csc|cot)!\\(left|right|Bigl|Bigr|bigl|bigr|Big|big|biggl|biggr|Biggl|Biggr)\b"),
     r"\\\1\\!\\\2"),
    # \right\right -> \right)\right  (missing closing delimiter, optional ) follows)
    (re.compile(r"\\right\\right(\))?"), r"\\right)\\right)"),
]

n_files = 0
for md in WP.glob("PAPER_*.md"):
    s = md.read_text(encoding="utf-8")
    orig = s
    for old, new in LITERAL:
        s = s.replace(old, new)
    for pat, rep in REGEX:
        s = pat.sub(rep, s)
    if s != orig:
        md.write_text(s, encoding="utf-8")
        n_files += 1

print(f"Updated {n_files} files.")

