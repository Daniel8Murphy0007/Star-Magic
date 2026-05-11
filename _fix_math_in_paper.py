"""
Focused math cleanup for a single whitepaper markdown file.
Applies safe transformations only inside $$...$$ display blocks
and $...$ inline math spans.

Usage:
    python _fix_math_in_paper.py whitepapers/PAPER_027_xxx.md
"""
import re, sys, pathlib

# (pattern, replacement) — applied INSIDE math segments only
MATH_FIXES = [
    # cmd-without-space: \pit_n -> \pi t_n etc.
    (re.compile(r"\\pi(?=[A-Za-z_])"), r"\\pi "),
    (re.compile(r"\\Delta(?=[A-Za-z_])"), r"\\Delta "),
    (re.compile(r"\\Sigma(?=[A-Za-z_])"), r"\\Sigma "),
    (re.compile(r"\\Gamma(?=[A-Za-z_])"), r"\\Gamma "),
    (re.compile(r"\\Lambda(?=[A-Za-z_])"), r"\\Lambda "),
    (re.compile(r"\\Phi(?=[A-Za-z_])"), r"\\Phi "),
    (re.compile(r"\\Omega(?=[A-Za-z_])"), r"\\Omega "),
    (re.compile(r"\\Theta(?=[A-Za-z_])"), r"\\Theta "),
    (re.compile(r"\\Psi(?=[A-Za-z_])"), r"\\Psi "),
    (re.compile(r"\\alpha(?=[A-Za-z_])"), r"\\alpha "),
    (re.compile(r"\\beta(?=[A-Za-z_])"), r"\\beta "),
    (re.compile(r"\\gamma(?=[A-Za-z_])"), r"\\gamma "),
    (re.compile(r"\\delta(?=[A-Za-z_])"), r"\\delta "),
    (re.compile(r"\\rho(?=[A-Za-z_])"), r"\\rho "),
    (re.compile(r"\\sigma(?=[A-Za-z_])"), r"\\sigma "),
    (re.compile(r"\\mu(?=[A-Za-z_])"), r"\\mu "),
    (re.compile(r"\\nu(?=[A-Za-z_])"), r"\\nu "),
    (re.compile(r"\\tau(?=[A-Za-z_])"), r"\\tau "),
    (re.compile(r"\\phi(?=[A-Za-z_])"), r"\\phi "),
    (re.compile(r"\\theta(?=[A-Za-z_])"), r"\\theta "),
    (re.compile(r"\\omega(?=[A-Za-z_])"), r"\\omega "),
    (re.compile(r"\\kappa(?=[A-Za-z_])"), r"\\kappa "),
    (re.compile(r"\\lambda(?=[A-Za-z_])"), r"\\lambda "),
    (re.compile(r"\\epsilon(?=[A-Za-z_])"), r"\\epsilon "),
    (re.compile(r"\\eta(?=[A-Za-z_])"), r"\\eta "),
    (re.compile(r"\\chi(?=[A-Za-z_])"), r"\\chi "),
    (re.compile(r"\\xi(?=[A-Za-z_])"), r"\\xi "),
    (re.compile(r"\\zeta(?=[A-Za-z_])"), r"\\zeta "),
    (re.compile(r"\\hbar(?=[A-Za-z_])"), r"\\hbar "),
    (re.compile(r"\\nabla(?=[A-Za-z_])"), r"\\nabla "),
    (re.compile(r"\\to(?=[A-Za-z_])"), r"\\to "),
    (re.compile(r"\\approx(?=[A-Za-z_])"), r"\\approx "),
    (re.compile(r"\\sim(?=[A-Za-z_])"), r"\\sim "),
    (re.compile(r"\\times(?=[A-Za-z_])"), r"\\times "),
    (re.compile(r"\\cdot(?=[A-Za-z_])"), r"\\cdot "),
    (re.compile(r"\\propto(?=[A-Za-z_])"), r"\\propto "),
    (re.compile(r"\\partial(?=[A-Za-z_])"), r"\\partial "),
    (re.compile(r"\\infty(?=[A-Za-z_])"), r"\\infty "),
    (re.compile(r"\\pm(?=[A-Za-z_])"), r"\\pm "),
    (re.compile(r"\\mp(?=[A-Za-z_])"), r"\\mp "),
    (re.compile(r"\\leq(?=[A-Za-z_])"), r"\\leq "),
    (re.compile(r"\\geq(?=[A-Za-z_])"), r"\\geq "),
    (re.compile(r"\\neq(?=[A-Za-z_])"), r"\\neq "),
    (re.compile(r"\\rightarrow(?=[A-Za-z_])"), r"\\rightarrow "),
    (re.compile(r"\\leftarrow(?=[A-Za-z_])"), r"\\leftarrow "),
    # bare function names -> backslashed
    (re.compile(r"(?<![A-Za-z\\])(cos|sin|tan|exp|ln|log|sec|csc|cot|tanh|sinh|cosh)\("), r"\\\1("),
    # 10-N -> 10^{-N} (one or two digit exponent)
    (re.compile(r"10-(\d{1,3})\b"), r"10^{-\1}"),
    # c2  -> c^2  (only when followed by ) or space, not letter)
    (re.compile(r"\bc2(?=[^A-Za-z0-9_]|$)"), r"c^2"),
    # m_X2 -> m_X^2  (greek/latin token then 2)
    (re.compile(r"(m_[A-Za-z]|m_\\[A-Za-z]+|G_F|m_W|m_Z|m_B|m_t|m_b|m_e|m_\\nu)2\b"), r"\1^2"),
    # subscript followed by underscore tokens stay—do not touch
]

INLINE_MATH = re.compile(r"(\$\$.*?\$\$|\$[^\n$]+?\$)", re.DOTALL)

def fix_math_segment(s):
    for pat, rep in MATH_FIXES:
        s = pat.sub(rep, s)
    return s

def fix_text(text):
    out = []
    last = 0
    for m in INLINE_MATH.finditer(text):
        out.append(text[last:m.start()])
        out.append(fix_math_segment(m.group(0)))
        last = m.end()
    out.append(text[last:])
    return "".join(out)

def main():
    if len(sys.argv) != 2:
        print("Usage: python _fix_math_in_paper.py <path>")
        sys.exit(2)
    p = pathlib.Path(sys.argv[1])
    src = p.read_text(encoding="utf-8")
    fixed = fix_text(src)
    if fixed == src:
        print("no change")
        return
    p.write_text(fixed, encoding="utf-8")
    print(f"updated {p}")

if __name__ == "__main__":
    main()
