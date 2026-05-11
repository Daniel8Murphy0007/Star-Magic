"""
Clean problematic unicode characters that pdflatex can't render even in
code/verbatim blocks. Converts to ASCII equivalents safe for plain text.
For math mode, conversions still parse correctly.

Usage: python _clean_unicode_paper.py path/to/PAPER_xxx.md
"""
import sys, pathlib, re, unicodedata

# Combining-accent letter -> ASCII replacement (after-decomposition)
# Letter + U+0304 (macron) -> letter+"bar"
# Letter + U+0302 (circumflex) -> letter+"hat"
# Letter + U+0307 (dot above) -> letter+"dot"
# Letter + U+0303 (tilde) -> letter+"tilde"

ACCENT_MAP = {
    "\u0304": "bar",
    "\u0302": "hat",
    "\u0307": "dot",
    "\u0303": "tilde",
    "\u0308": "ddot",
    "\u030A": "ring",
    "\u030C": "check",
    "\u0301": "acute",
    "\u0300": "grave",
}

# Superscript / subscript unicode -> ASCII forms
SUPER_SUB = {
    "⁰":"^0","¹":"^1","²":"^2","³":"^3","⁴":"^4","⁵":"^5","⁶":"^6","⁷":"^7","⁸":"^8","⁹":"^9",
    "⁺":"^+","⁻":"^-","⁽":"^(","⁾":"^)","ⁿ":"^n",
    "₀":"_0","₁":"_1","₂":"_2","₃":"_3","₄":"_4","₅":"_5","₆":"_6","₇":"_7","₈":"_8","₉":"_9",
    "₊":"_+","₋":"_-","₍":"_(","₎":"_)",
    "ₐ":"_a","ₑ":"_e","ₒ":"_o","ₓ":"_x","ₕ":"_h","ₖ":"_k","ₗ":"_l","ₘ":"_m","ₙ":"_n",
    "ₚ":"_p","ₛ":"_s","ₜ":"_t","ᵢ":"_i","ⱼ":"_j",
}

# Other unicode that often appears in code blocks and breaks pdflatex
MISC = {
    "ℓ": "l",   # script l in code/text (in math, \ell would be better)
    "Ṁ": "Mdot",
    "ħ": "hbar",
    "↕": "<->",
    "↔": "<->",
    "∎": "QED",
    "→": "->",
    "←": "<-",
    "⇒": "=>",
    "⇐": "<=",
    "≈": "~",
    "≃": "~=",
    "≠": "!=",
    "≤": "<=",
    "≥": ">=",
    "≡": "==",
    "∼": "~",
    "∝": "prop",
    "∞": "inf",
    "·": "*",
    "×": "x",   # only in code/text contexts (we already remap $\times$ in regen)
    "±": "+/-",
    "∓": "-/+",
    "∂": "d",
    "∇": "grad",
    "∫": "int",
    "∮": "oint",
    "∑": "sum",
    "∏": "prod",
    "√": "sqrt",
    "∈": "in",
    "⊂": "subset",
    "∪": "union",
    "∩": "intersect",
    "∅": "empty",
    "∀": "forall",
    "∃": "exists",
    "°": "deg",
    "ℏ": "hbar",
    "ℝ": "R",
    "ℂ": "C",
    "ℕ": "N",
    "ℤ": "Z",
    "⟨": "<",
    "⟩": ">",
    "–": "-",
    "—": "--",
    "…": "...",
    "α":"alpha","β":"beta","γ":"gamma","δ":"delta","ε":"epsilon","ζ":"zeta",
    "η":"eta","θ":"theta","ι":"iota","κ":"kappa","λ":"lambda","μ":"mu",
    "ν":"nu","ξ":"xi","π":"pi","ρ":"rho","σ":"sigma","τ":"tau","υ":"upsilon",
    "φ":"phi","χ":"chi","ψ":"psi","ω":"omega","ϑ":"theta","ϕ":"phi","ϵ":"epsilon",
    "Γ":"Gamma","Δ":"Delta","Θ":"Theta","Λ":"Lambda","Ξ":"Xi","Π":"Pi","Σ":"Sigma",
    "Υ":"Upsilon","Φ":"Phi","Ψ":"Psi","Ω":"Omega",
}

# Regions to clean unconditionally: code blocks, indented code, tables.
# For everything else (markdown body, math), we preserve LaTeX-compatible chars
# already handled by _pdf_unicode_header.tex.

CODE_BLOCK = re.compile(r"```.*?```", re.DOTALL)

def clean_combining_accents_in_codeblocks(text):
    def cb(m):
        block = m.group(0)
        # Decompose then re-pair letter+accent
        d = unicodedata.normalize("NFD", block)
        out_chars = []
        i = 0
        while i < len(d):
            ch = d[i]
            if i + 1 < len(d) and d[i+1] in ACCENT_MAP:
                out_chars.append(f"{ch}{ACCENT_MAP[d[i+1]]}")
                i += 2
            else:
                out_chars.append(ch)
                i += 1
        cleaned = "".join(out_chars)
        # Now replace remaining mapped unicode
        for u, a in {**SUPER_SUB, **MISC}.items():
            cleaned = cleaned.replace(u, a)
        return cleaned
    return CODE_BLOCK.sub(cb, text)

def clean_combining_accents_in_text(text):
    """Outside code blocks: convert letter+combining-accent to wrapped LaTeX.
    Wraps in $...$ when emitted outside existing math mode.
    """
    LATEX_ACCENT = {
        "\u0304": "bar",
        "\u0302": "hat",
        "\u0307": "dot",
        "\u0303": "tilde",
        "\u0308": "ddot",
    }
    d = unicodedata.normalize("NFD", text)
    out = []
    i = 0
    in_math = False
    while i < len(d):
        ch = d[i]
        if ch == "$":
            # Track math state, handling $$ as toggle of display math
            if i + 1 < len(d) and d[i+1] == "$":
                in_math = not in_math
                out.append("$$")
                i += 2
                continue
            in_math = not in_math
            out.append("$")
            i += 1
            continue
        if i + 1 < len(d) and d[i+1] in LATEX_ACCENT and ch.isalpha():
            name = LATEX_ACCENT[d[i+1]]
            if in_math:
                out.append(f"\\{name}{{{ch}}}")
            else:
                out.append(f"${{\\{name} {{{ch}}}}}$")
            i += 2
        else:
            out.append(ch)
            i += 1
    return "".join(out)

def main():
    if len(sys.argv) != 2:
        print("usage: _clean_unicode_paper.py <path>")
        sys.exit(2)
    p = pathlib.Path(sys.argv[1])
    src = p.read_text(encoding="utf-8")
    # Process code blocks first (ASCII replacement)
    step1 = clean_combining_accents_in_codeblocks(src)
    # Then process remaining combining accents in body (LaTeX form)
    step2 = clean_combining_accents_in_text(step1)
    if step2 == src:
        print("no change")
        return
    p.write_text(step2, encoding="utf-8")
    print(f"updated {p}")

if __name__ == "__main__":
    main()
