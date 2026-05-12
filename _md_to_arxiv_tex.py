"""
_md_to_arxiv_tex.py — Star-Magic markdown -> arXiv-canonical LaTeX.
NO pandoc. NO xelatex. NO reportlab. Direct text->LaTeX, then pdflatex.

v2 (2026-05-11): proper title/metadata block, bold-across-math, hard line
breaks, LaTeX-aware escaping (don't double-escape \_, \{, \}; don't escape
braces or backslashes — author already writes inline LaTeX in prose).
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

# ----------------------------------------------------------------------------
# Preamble — matches whitepapers/PAPER_495 and Manuscript 1_12Feb2026
# ----------------------------------------------------------------------------
PREAMBLE = r"""\documentclass[12pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,amsfonts,amsthm}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage{booktabs}
\usepackage{array}
\usepackage{longtable}
\usepackage{multicol}
\usepackage{geometry}
\usepackage{float}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{enumitem}
\usepackage{textcomp}
\usepackage{xcolor}
\usepackage{tabularx}
\geometry{margin=1in}
\usepackage{url}
\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue,citecolor=blue,breaklinks=true}
\setlength{\parskip}{0.5em}
\setlength{\parindent}{0pt}
% Prevent overfull \hbox warnings (margins blown by long URLs/identifiers)
\sloppy
\emergencystretch=3em
\setcounter{secnumdepth}{3}
\setcounter{tocdepth}{3}
% Allow URL-style break points inside long identifiers like MAIN_1_CoAnQi
\makeatletter
\def\UrlBreaks{\do\.\do\@\do\\\do\/\do\!\do\_\do\?\do\&\do\<\do\>\do\%\do\-\do\+\do\=\do\:\do\;\do\,}
\makeatother
"""

TYPEOUT = r"\typeout{get arXiv to do 4 passes: Label(s) may have changed. Rerun}"

# ----------------------------------------------------------------------------
# Unicode -> LaTeX
# ----------------------------------------------------------------------------
UNICODE_MAP = {
    "α": r"$\alpha$", "β": r"$\beta$", "γ": r"$\gamma$", "δ": r"$\delta$",
    "ε": r"$\varepsilon$", "ζ": r"$\zeta$", "η": r"$\eta$", "θ": r"$\theta$",
    "ι": r"$\iota$", "κ": r"$\kappa$", "λ": r"$\lambda$", "μ": r"$\mu$",
    "ν": r"$\nu$", "ξ": r"$\xi$", "π": r"$\pi$",
    "ρ": r"$\rho$", "σ": r"$\sigma$", "τ": r"$\tau$", "υ": r"$\upsilon$",
    "φ": r"$\phi$", "χ": r"$\chi$", "ψ": r"$\psi$", "ω": r"$\omega$",
    "Γ": r"$\Gamma$", "Δ": r"$\Delta$", "Θ": r"$\Theta$",
    "Λ": r"$\Lambda$", "Ξ": r"$\Xi$", "Π": r"$\Pi$",
    "Σ": r"$\Sigma$", "Υ": r"$\Upsilon$",
    "Φ": r"$\Phi$", "Ψ": r"$\Psi$", "Ω": r"$\Omega$",
    "∝": r"$\propto$", "∓": r"$\mp$", "±": r"$\pm$",
    "≳": r"$\gtrsim$", "≲": r"$\lesssim$", "≈": r"$\approx$",
    "≠": r"$\neq$", "≤": r"$\leq$", "≥": r"$\geq$", "≡": r"$\equiv$",
    "∞": r"$\infty$", "∂": r"$\partial$", "∇": r"$\nabla$",
    "∈": r"$\in$", "∉": r"$\notin$", "⊂": r"$\subset$", "⊃": r"$\supset$",
    "∑": r"$\sum$", "∏": r"$\prod$", "∫": r"$\int$",
    "→": r"$\to$", "←": r"$\leftarrow$", "↔": r"$\leftrightarrow$",
    "⇒": r"$\Rightarrow$", "⇐": r"$\Leftarrow$", "⇔": r"$\Leftrightarrow$",
    "·": r"$\cdot$", "×": r"$\times$", "÷": r"$\div$",
    "°": r"$^{\circ}$", "′": r"$'$", "″": r"$''$",
    "—": r"---", "–": r"--", "…": r"\ldots{}",
    "“": r"``", "”": r"''", "‘": r"`", "’": r"'",
    "•": r"\textbullet{}",
    "⊙": r"$\odot$", "⊕": r"$\oplus$", "⊗": r"$\otimes$",
    "ℏ": r"$\hbar$", "ℓ": r"$\ell$",
    "□": r"$\Box$", "ℳ": r"$\mathcal{M}$",
    "ℰ": r"$\mathcal{E}$", "ℋ": r"$\mathcal{H}$",
    "ℒ": r"$\mathcal{L}$", "ℐ": r"$\mathcal{I}$",
    "ℝ": r"$\mathbb{R}$", "ℕ": r"$\mathbb{N}$",
    "ℤ": r"$\mathbb{Z}$", "ℂ": r"$\mathbb{C}$", "ℚ": r"$\mathbb{Q}$",
    "∼": r"$\sim$", "≅": r"$\cong$",
    "²": r"$^{2}$", "³": r"$^{3}$", "¹": r"$^{1}$", "⁰": r"$^{0}$",
    "⁴": r"$^{4}$", "⁵": r"$^{5}$", "⁶": r"$^{6}$", "⁷": r"$^{7}$",
    "⁸": r"$^{8}$", "⁹": r"$^{9}$",
    "₀": r"$_{0}$", "₁": r"$_{1}$", "₂": r"$_{2}$", "₃": r"$_{3}$",
    "₄": r"$_{4}$", "₅": r"$_{5}$", "₆": r"$_{6}$", "₇": r"$_{7}$",
    "₈": r"$_{8}$", "₉": r"$_{9}$",
    # Modifier letter sub/superscripts (Phonetic Extensions block U+1D00-U+1D7F, Latin Sub/Superscripts U+2080-U+209F)
    "ᵃ": r"$^{a}$", "ᵇ": r"$^{b}$", "ᶜ": r"$^{c}$", "ᵈ": r"$^{d}$",
    "ᵉ": r"$^{e}$", "ᶠ": r"$^{f}$", "ᵍ": r"$^{g}$", "ʰ": r"$^{h}$",
    "ⁱ": r"$^{i}$", "ʲ": r"$^{j}$", "ᵏ": r"$^{k}$", "ˡ": r"$^{l}$",
    "ᵐ": r"$^{m}$", "ⁿ": r"$^{n}$", "ᵒ": r"$^{o}$", "ᵖ": r"$^{p}$",
    "ʳ": r"$^{r}$", "ˢ": r"$^{s}$", "ᵗ": r"$^{t}$", "ᵘ": r"$^{u}$",
    "ᵛ": r"$^{v}$", "ʷ": r"$^{w}$", "ˣ": r"$^{x}$", "ʸ": r"$^{y}$", "ᶻ": r"$^{z}$",
    "ₐ": r"$_{a}$", "ₑ": r"$_{e}$", "ₕ": r"$_{h}$", "ᵢ": r"$_{i}$",
    "ⱼ": r"$_{j}$", "ₖ": r"$_{k}$", "ₗ": r"$_{l}$", "ₘ": r"$_{m}$",
    "ₙ": r"$_{n}$", "ₒ": r"$_{o}$", "ₚ": r"$_{p}$", "ᵣ": r"$_{r}$",
    "ₛ": r"$_{s}$", "ₜ": r"$_{t}$", "ᵤ": r"$_{u}$", "ᵥ": r"$_{v}$", "ₓ": r"$_{x}$",
    "⁺": r"$^{+}$", "⁻": r"$^{-}$", "⁼": r"$^{=}$", "⁽": r"$^{(}$", "⁾": r"$^{)}$",
    "₊": r"$_{+}$", "₋": r"$_{-}$", "₌": r"$_{=}$", "₍": r"$_{(}$", "₎": r"$_{)}$",
    "Å": r"\AA{}", "ø": r"\o{}", "Ø": r"\O{}", "ß": r"\ss{}",
    "£": r"\pounds{}", "˜": r"$\sim$",
    "✅": r"$\checkmark$", "✓": r"$\checkmark$", "❌": r"$\times$",
    "⚠": r"!", "★": r"$\star$", "☆": r"$\star$",
    "§": r"\S{}", "¶": r"\P{}", "†": r"$\dagger$", "‡": r"$\ddagger$",
    "©": r"\textcopyright{}", "®": r"\textregistered{}", "™": r"\texttrademark{}",
    "\u00a0": " ",
    "\ufeff": "",
    "\ufe0f": "",
    "\ufe0e": "",
    "\u200b": "", "\u200c": "", "\u200d": "",
}


def map_unicode(s: str) -> str:
    for u, repl in UNICODE_MAP.items():
        if u in s:
            s = s.replace(u, repl)
    return s


# Inside math contexts, Latin lookalikes typed instead of Greek letters
# (e.g. "ß" used for Greek β) need math-mode mappings. \ss{} is text-mode
# only and triggers "Undefined control sequence" inside math.
MATH_UNICODE_MAP = {
    "\u00df": r"\beta",     # ß typo -> \beta
    "\u03b1": r"\alpha",
    "\u03b2": r"\beta",
    "\u03b3": r"\gamma",
    "\u03b4": r"\delta",
    "\u03b5": r"\epsilon",
    "\u03b6": r"\zeta",
    "\u03b7": r"\eta",
    "\u03b8": r"\theta",
    "\u03ba": r"\kappa",
    "\u03bb": r"\lambda",
    "\u03bc": r"\mu",
    "\u03bd": r"\nu",
    "\u03be": r"\xi",
    "\u03c0": r"\pi",
    "\u03c1": r"\rho",
    "\u03c3": r"\sigma",
    "\u03c4": r"\tau",
    "\u03c6": r"\phi",
    "\u03c7": r"\chi",
    "\u03c8": r"\psi",
    "\u03c9": r"\omega",
    "\u0394": r"\Delta",
    "\u0398": r"\Theta",
    "\u039b": r"\Lambda",
    "\u03a0": r"\Pi",
    "\u03a3": r"\Sigma",
    "\u03a6": r"\Phi",
    "\u03a8": r"\Psi",
    "\u03a9": r"\Omega",
    "\u2207": r"\nabla",
    "\u2202": r"\partial",
    "\u221e": r"\infty",
    "\u00b1": r"\pm",
    "\u2213": r"\mp",
    "\u00d7": r"\times",
    "\u00f7": r"\div",
    "\u2265": r"\geq",
    "\u2264": r"\leq",
    "\u2260": r"\neq",
    "\u2248": r"\approx",
    "\u2261": r"\equiv",
    "\u221d": r"\propto",
    "\u22c5": r"\cdot",
    "\u2192": r"\to",
    "\u21d2": r"\Rightarrow",
    "\u210f": r"\hbar",
    "\u2113": r"\ell",
    "\u25a1": r"\Box",
    "\u2133": r"\mathcal{M}",
    "\u2130": r"\mathcal{E}",
    "\u210b": r"\mathcal{H}",
    "\u2112": r"\mathcal{L}",
    "\u2110": r"\mathcal{I}",
    "\u211d": r"\mathbb{R}",
    "\u2115": r"\mathbb{N}",
    "\u2124": r"\mathbb{Z}",
    "\u2102": r"\mathbb{C}",
    "\u211a": r"\mathbb{Q}",
    "\u2299": r"\odot",
    "\u2295": r"\oplus",
    "\u2297": r"\otimes",
    "\u223c": r"\sim",
    "\u2245": r"\cong",
}


def map_unicode_math(s: str) -> str:
    for u, repl in MATH_UNICODE_MAP.items():
        if u in s:
            # Add trailing space if next char is a letter (so \beta_x not \betax)
            # NOTE: pass replacement as a lambda to avoid re.sub interpreting
            # backslash sequences (e.g. \B) as template escapes.
            spaced = repl + " "
            s = re.sub(re.escape(u) + r"(?=[A-Za-z])", lambda _m, r=spaced: r, s)
            s = s.replace(u, repl)
    return s


# Unicode characters inside verbatim blocks can't be rendered by T1 fonts.
# Replace common math sub/superscripts and Greek letters with ASCII so the
# verbatim still compiles. (Verbatim already loses Markdown fidelity; ASCII
# is acceptable.)
_VERBATIM_ASCII_MAP = {
    # Subscripts U+2080-U+2089
    "\u2080": "_0", "\u2081": "_1", "\u2082": "_2", "\u2083": "_3", "\u2084": "_4",
    "\u2085": "_5", "\u2086": "_6", "\u2087": "_7", "\u2088": "_8", "\u2089": "_9",
    "\u208a": "_+", "\u208b": "_-", "\u208c": "_=", "\u208d": "_(", "\u208e": "_)",
    # Superscripts U+2070-U+207F
    "\u2070": "^0", "\u00b9": "^1", "\u00b2": "^2", "\u00b3": "^3",
    "\u2074": "^4", "\u2075": "^5", "\u2076": "^6", "\u2077": "^7",
    "\u2078": "^8", "\u2079": "^9",
    "\u207a": "^+", "\u207b": "^-", "\u207c": "^=", "\u207d": "^(", "\u207e": "^)",
    # Greek lowercase (common in code comments)
    "\u03b1": "alpha", "\u03b2": "beta", "\u03b3": "gamma", "\u03b4": "delta",
    "\u03b5": "epsilon", "\u03b6": "zeta", "\u03b7": "eta", "\u03b8": "theta",
    "\u03ba": "kappa", "\u03bb": "lambda", "\u03bc": "mu", "\u03bd": "nu",
    "\u03be": "xi", "\u03c0": "pi", "\u03c1": "rho", "\u03c3": "sigma",
    "\u03c4": "tau", "\u03c6": "phi", "\u03c7": "chi", "\u03c8": "psi",
    "\u03c9": "omega",
    "\u0394": "Delta", "\u0398": "Theta", "\u039b": "Lambda", "\u03a0": "Pi",
    "\u03a3": "Sigma", "\u03a6": "Phi", "\u03a8": "Psi", "\u03a9": "Omega",
    # Operators
    "\u00d7": "x", "\u00f7": "/", "\u00b1": "+/-",
    "\u2265": ">=", "\u2264": "<=", "\u2260": "!=", "\u2248": "~",
    "\u2192": "->", "\u21d2": "=>", "\u22c5": "*", "\u221e": "inf",
    "\u00b0": "deg", "\u03bc": "u",
    # Misc
    "\u00df": "ss", "\u2026": "...", "\u2013": "-", "\u2014": "--",
    "\u201c": '"', "\u201d": '"', "\u2018": "'", "\u2019": "'",
    # Math-letter unicode that can appear in code spans
    "\u210f": "hbar", "\u2113": "l",
    "\u25a1": "Box", "\u2133": "M", "\u2130": "E",
    "\u210b": "H", "\u2112": "L", "\u2110": "I",
    "\u211d": "R", "\u2115": "N", "\u2124": "Z", "\u2102": "C", "\u211a": "Q",
    "\u2299": "(.)", "\u2295": "(+)", "\u2297": "(x)",
    "\u223c": "~", "\u2245": "~=",
    # Strip zero-width / variation selectors
    "\ufeff": "", "\ufe0f": "", "\ufe0e": "",
    "\u200b": "", "\u200c": "", "\u200d": "", "\u00a0": " ",
}


def sanitize_verbatim(s: str) -> str:
    for u, repl in _VERBATIM_ASCII_MAP.items():
        if u in s:
            s = s.replace(u, repl)
    # Final pass: replace any remaining non-ASCII char with '?' so pdflatex
    # can compile. (Pure ASCII verbatim is required by T1 font encoding.)
    return s.encode("ascii", errors="replace").decode("ascii")


# ----------------------------------------------------------------------------
# LaTeX-aware escape:
#   - Pass through `\X` (already-escaped or LaTeX command)
#   - Leave `{` `}` untouched (LaTeX grouping; author writes inline LaTeX)
#   - Escape only: _ % # & ~ $ (where not already preceded by \)
# ----------------------------------------------------------------------------
def escape_text(s: str) -> str:
    out = []
    i = 0
    n = len(s)
    while i < n:
        ch = s[i]
        if ch == "\\" and i + 1 < n:
            out.append(s[i:i + 2])
            i += 2
        elif ch == "_":
            out.append(r"\_")
            i += 1
        elif ch == "%":
            out.append(r"\%")
            i += 1
        elif ch == "#":
            out.append(r"\#")
            i += 1
        elif ch == "&":
            out.append(r"\&")
            i += 1
        elif ch == "~":
            out.append(r"\textasciitilde{}")
            i += 1
        elif ch == "$":
            out.append(r"\$")
            i += 1
        else:
            out.append(ch)
            i += 1
    return "".join(out)


# ----------------------------------------------------------------------------
# Inline processing: stash math/code, run markdown subs, escape, restore.
# Placeholders are alphanumeric so they survive escape unchanged.
# ----------------------------------------------------------------------------
# Text-mode LaTeX commands we should NOT auto-mathify.
TEXT_MODE_CMDS = {
    "textbf", "textit", "emph", "texttt", "verb", "section", "subsection",
    "subsubsection", "paragraph", "subparagraph", "href", "url", "title",
    "author", "date", "maketitle", "begin", "end", "item", "label", "ref",
    "cite", "footnote", "noindent", "indent", "newpage", "clearpage",
    "newline", "linebreak", "smallskip", "medskip", "bigskip", "hrule",
    "hline", "toprule", "midrule", "bottomrule", "endhead", "tabularnewline",
    "pounds", "AA", "aa", "o", "O", "ss", "S", "P", "ldots", "textbullet",
    "textcopyright", "textregistered", "texttrademark", "textasciitilde",
    "textasciicircum", "rule", "flushleft", "flushright", "center", "left",
    "right", "vspace", "hspace", "tableofcontents", "appendix", "chapter",
    "part", "today", "thepage", "pagestyle", "thispagestyle",
    "renewcommand", "newcommand", "providecommand", "setlength", "addtolength",
    "documentclass", "usepackage", "hypersetup", "typeout", "input", "include",
    "checkmark", "dagger", "ddagger", "S", "P",
}

# Math-leak detector: after inline rendering, find `\command` that aren't
# text-mode and aren't already inside math (`$...$` or `\(...\)` or `\[...\]`).
_MATH_LEAK_RE = re.compile(r"\\([a-zA-Z]+)")


def _is_in_math(text: str, pos: int) -> bool:
    """True if position pos in text is inside a $...$ or \\(...\\) span."""
    # Count unescaped $ before pos. Odd = inside inline math.
    dollars = 0
    i = 0
    while i < pos:
        ch = text[i]
        if ch == "\\" and i + 1 < pos:
            i += 2
            continue
        if ch == "$":
            # Skip $$ as display delimiters (count as 0)
            if i + 1 < pos and text[i + 1] == "$":
                i += 2
                continue
            dollars += 1
        i += 1
    return dollars % 2 == 1


def _escape_bare_carets(s: str) -> str:
    """Escape any remaining bare ^ outside inline math as literal `\\^{}`."""
    out = []
    for i, ch in enumerate(s):
        if ch == "^" and not _is_in_math(s, i):
            out.append(r"\^{}")
        else:
            out.append(ch)
    return "".join(out)


def wrap_math_leaks(s: str) -> str:
    """Wrap stray LaTeX math commands in inline math to prevent compile errors.

    Only handles standalone tokens — `\\command{arg}` is wrapped along with
    the immediately following braced group if present. Heuristic, but safe
    for malformed source markdown.
    """
    # Pre-pass: split known math commands that the author wrote without a
    # separator before a trailing letter (e.g. `\cdotU`, `\timesM`, `\Deltax`).
    # The tokenizer below treats `\cdotU` as the single (undefined) command
    # "cdotU", so we must insert a space here BEFORE tokenizing.
    # Longer alternates must come first; Python regex matches left-to-right.
    _KNOWN_MATH_CMDS = (
        r"cdots|cdot|ldots|times|div|pm|mp|leq|geq|neq|approx|equiv|propto|"
        r"mapsto|gets|to|nabla|partial|infty|forall|exists|notin|in|"
        r"subset|supset|cup|cap|wedge|vee|land|lor|neg|hbar|ell|"
        r"varepsilon|epsilon|vartheta|theta|varphi|phi|varpi|pi|varrho|rho|"
        r"varsigma|sigma|upsilon|omega|alpha|beta|gamma|delta|zeta|eta|"
        r"iota|kappa|lambda|nu|xi|tau|chi|psi|mu|"
        r"Alpha|Beta|Gamma|Delta|Epsilon|Zeta|Eta|Theta|Iota|Kappa|Lambda|"
        r"Mu|Nu|Xi|Pi|Rho|Sigma|Tau|Upsilon|Phi|Chi|Psi|Omega|"
        r"sinh|cosh|tanh|sin|cos|tan|cot|sec|csc|log|ln|lg|exp|"
        r"iiint|iint|oint|int|sum|prod|lim|"
        r"det|dim|ker|deg|gcd"
    )
    _split_re = re.compile(r"\\(" + _KNOWN_MATH_CMDS + r")(?=[A-Za-z])")
    s = _split_re.sub(lambda m: "\\" + m.group(1) + " ", s)
    # Pre-pass 2: globally wrap math-only commands followed by an optional
    # sub/superscript chain into inline math, so they survive even inside
    # text-mode wrappers like \textbf{...} (which the tokenizer below would
    # otherwise emit verbatim without recursing). Handles both raw `_{..}`
    # and the post-escape form `\_{..}`.
    _MATH_LEAK_RE = re.compile(
        r"\\(" + _KNOWN_MATH_CMDS + r")"
        r"((?:\\?[_^](?:\{[^{}]*\}|[A-Za-z0-9]))*)"
    )

    def _wrap_leak(m, src=s):
        if _is_in_math(src, m.start()):
            return m.group(0)
        return "$\\" + m.group(1) + m.group(2) + "$"

    s = _MATH_LEAK_RE.sub(_wrap_leak, s)
    out = []
    i = 0
    n = len(s)
    while i < n:
        if s[i] == "\\" and i + 1 < n and s[i + 1].isalpha():
            j = i + 1
            while j < n and s[j].isalpha():
                j += 1
            name = s[i + 1:j]
            # Trailing braced arg
            k = j
            depth = 0
            if k < n and s[k] == "{":
                depth = 1
                k += 1
                while k < n and depth > 0:
                    if s[k] == "\\" and k + 1 < n:
                        k += 2
                        continue
                    if s[k] == "{":
                        depth += 1
                    elif s[k] == "}":
                        depth -= 1
                    k += 1
            # Also consume trailing sub/superscript chains: _{..}, ^{..},
            # _x, ^x (so `\log_{10}` wraps as a single math span, not just
            # `\log`). Also handle the post-escape form `\_{..}` `\^{..}`
            # produced by escape_text running before wrap_math_leaks.
            while k < n:
                # `_` or `^` directly
                if s[k] in ("_", "^"):
                    k += 1
                # escaped `\_` or `\^`
                elif s[k] == "\\" and k + 1 < n and s[k + 1] in ("_", "^"):
                    k += 2
                else:
                    break
                if k < n and s[k] == "{":
                    depth = 1
                    k += 1
                    while k < n and depth > 0:
                        if s[k] == "\\" and k + 1 < n:
                            k += 2
                            continue
                        if s[k] == "{":
                            depth += 1
                        elif s[k] == "}":
                            depth -= 1
                        k += 1
                elif k < n and (s[k].isalnum() or s[k] == "\\"):
                    if s[k] == "\\" and k + 1 < n and s[k + 1].isalpha():
                        k += 2
                        while k < n and s[k].isalpha():
                            k += 1
                    else:
                        k += 1
            token = s[i:k]
            if name in TEXT_MODE_CMDS or _is_in_math(s, i):
                out.append(token)
            else:
                out.append("$" + token + "$")
            i = k
        else:
            out.append(s[i])
            i += 1
    return "".join(out)


def render_inline(line: str) -> str:
    stash: dict[str, str] = {}

    def make_token(prefix: str, value: str) -> str:
        key = f"ZZ{prefix}{len(stash):04d}ZZ"
        stash[key] = value
        return key

    # Stash inline code FIRST (before any unicode mapping; code is literal)
    line = re.sub(r"`([^`\n]+)`", lambda m: make_token("CODE", m.group(1)), line)
    # Authors sometimes write adjacent inline math like "$\Delta$$\alpha$"
    # without separation. The block-level $$ display handler already ran in
    # process_lines (we only see inline content here), so we can safely
    # split any remaining "$$" into "$ $" to make each side a clean inline
    # math span before stashing.
    line = line.replace("$$", "$ $")
    # Stash display math $$...$$ (rare inline cases that survive the split,
    # e.g. authored deliberately as inline display)
    line = re.sub(r"\$\$(.+?)\$\$",
                  lambda m: make_token("DMATH", m.group(1)), line, flags=re.DOTALL)
    # Stash inline math $...$
    line = re.sub(r"(?<!\$)\$([^\$]+?)\$(?!\$)",
                  lambda m: make_token("MATH", m.group(1)), line, flags=re.DOTALL)

    # Now apply text-mode unicode mapping to the remaining prose only.
    line = map_unicode(line)

    # map_unicode may have introduced new inline math like `$\hbar$` for
    # symbols (□, ℳ, ℏ, etc.). Stash those too so escape_text doesn't
    # mangle the `$` into `\$` (which would leave \mathcal in text mode).
    line = re.sub(r"(?<!\$)\$([^\$]+?)\$(?!\$)",
                  lambda m: make_token("MATH", m.group(1)), line, flags=re.DOTALL)

    # Wrap bare ASCII scientific notation like `10^10`, `10^-4`, `x^{n}` so
    # they survive text mode. Done AFTER math stashing so existing math
    # expressions are protected. Wrap into a MATH placeholder directly.
    def caret_sub(m):
        content = m.group(1) + "^{" + m.group(2).strip("{}") + "}"
        return make_token("MATH", content)

    line = re.sub(r"([A-Za-z0-9])\^(\{[^}]+\}|-?\d+|[A-Za-z])",
                  caret_sub, line)

    # Auto-fix half-text / half-math leaks that authors write in prose.
    # Done AFTER math stashing so existing math expressions are protected.
    # 1) "MM_sun" / "M_sun" / "M_\\odot" (typo) -> $M_\odot$
    line = re.sub(r"\bMM[_\\]+sun\b",
                  lambda m: make_token("MATH", r"M_{\odot}"), line)
    line = re.sub(r"\bM[_\\]+sun\b",
                  lambda m: make_token("MATH", r"M_{\odot}"), line)
    # 2) Text-mode negative exponents like "10-22" / "10-4" -> $10^{-N}$
    line = re.sub(r"\b10-(\d+)\b",
                  lambda m: make_token("MATH", f"10^{{-{m.group(1)}}}"), line)
    # 3) Unit-with-negative-exponent like "day-1", "s-1", "Hz-1"
    line = re.sub(
        r"\b(day|days|s|sec|Hz|m|cm|mm|km|kg|g|yr|year|K|eV|keV|MeV|GeV|J|mol|nm|pc|kpc|Mpc)-(\d+)\b",
        lambda m: make_token("MATH",
                              r"\mathrm{" + m.group(1) + "}^{-" + m.group(2) + "}"),
        line)

    # Links [text](url)
    def link_sub(m):
        text = m.group(1)
        url = m.group(2)
        url_tok = make_token("URL", url)
        return r"\href{" + url_tok + r"}{" + text + r"}"
    line = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", link_sub, line)

    # Bold **...** and __...__
    line = re.sub(r"\*\*([^*\n][^*\n]*?)\*\*",
                  lambda m: r"\textbf{" + m.group(1) + r"}", line)
    line = re.sub(r"__([^_\n][^_\n]*?)__",
                  lambda m: r"\textbf{" + m.group(1) + r"}", line)
    # Italic *...*
    line = re.sub(r"(?<![*\w])\*([^*\n]+?)\*(?!\*)",
                  lambda m: r"\emph{" + m.group(1) + r"}", line)

    # Escape LaTeX specials in remaining prose
    line = escape_text(line)

    # Restore tokens
    def restore(m):
        key = m.group(0)
        v = stash[key]
        if key.startswith("ZZCODE"):
            # Always use \texttt so it works inside \emph/\textbf (\verb is
            # fragile and breaks inside command arguments). Inside the code
            # body: literal braces -> \{ \}, literal backslash -> \textbackslash\
            # (no braces, since braces in code are usually literal). Then
            # escape_text handles _ % # & ~ $.
            # First sanitize unicode to ASCII (T1 font can't render arbitrary
            # unicode inside \texttt -- same constraint as verbatim).
            v = sanitize_verbatim(v)
            v_e = v.replace("{", r"\{").replace("}", r"\}")
            v_e = v_e.replace("\\", r"\textbackslash ")
            v_e = escape_text(v_e)
            return r"\texttt{" + v_e + r"}"
        if key.startswith("ZZDMATH"):
            return f"\\({map_unicode_math(v)}\\)"
        if key.startswith("ZZMATH"):
            return f"${map_unicode_math(v)}$"
        if key.startswith("ZZURL"):
            return v
        return v

    line = re.sub(r"ZZ(?:CODE|DMATH|MATH|URL)\d{4}ZZ", restore, line)
    line = wrap_math_leaks(line)
    line = _escape_bare_carets(line)
    return line


# ----------------------------------------------------------------------------
# Block parsing
# ----------------------------------------------------------------------------
def parse_frontmatter(lines: list[str]) -> tuple[dict, list[str]]:
    meta: dict[str, str] = {}
    if not lines:
        return meta, lines
    if lines[0].startswith("\ufeff"):
        lines = [lines[0].lstrip("\ufeff")] + lines[1:]
    if not lines or lines[0].strip() != "---":
        return meta, lines
    end = None
    for i in range(1, len(lines)):
        if lines[i].strip() == "---":
            end = i
            break
    if end is None:
        return meta, lines
    for ln in lines[1:end]:
        if ":" in ln:
            k, _, v = ln.partition(":")
            meta[k.strip()] = v.strip().strip('"').strip("'")
    return meta, lines[end + 1:]


def render_table(rows: list[str]) -> str:
    def cells(line: str) -> list[str]:
        line = line.strip()
        if line.startswith("|"):
            line = line[1:]
        if line.endswith("|"):
            line = line[:-1]
        return [c.strip() for c in line.split("|")]
    if len(rows) < 2:
        return ""
    header = cells(rows[0])
    ncols = len(header)
    # Build all body rows first so we can decide column strategy.
    body: list[list[str]] = []
    for r in rows[2:]:
        c = cells(r)
        if len(c) < ncols:
            c += [""] * (ncols - len(c))
        elif len(c) > ncols:
            c = c[:ncols]
        body.append(c)

    # Estimate max cell width per column to pick fixed-width vs flexible.
    col_max = [len(h) for h in header]
    for r in body:
        for j, v in enumerate(r):
            if len(v) > col_max[j]:
                col_max[j] = len(v)

    # Use tabularx with X columns when any cell is wide enough to overflow.
    # X columns wrap automatically inside \linewidth.
    use_tabularx = any(w > 18 for w in col_max) or ncols >= 4

    out = []

    def fmt(cs):
        rendered = [render_inline(c) for c in cs]
        rendered = [("{}" + r) if r.lstrip().startswith("[") else r for r in rendered]
        return " & ".join(rendered) + r" \tabularnewline"

    if use_tabularx:
        # First column left-aligned ragged-right, rest equal X
        # so labels stay readable and prose wraps.
        col_spec = ">{\\raggedright\\arraybackslash}X" * ncols
        out.append(r"\begin{table}[H]")
        out.append(r"\centering")
        out.append(r"\small")
        out.append(r"\begin{tabularx}{\linewidth}{@{}" + col_spec + r"@{}}")
        out.append(r"\toprule")
        out.append(fmt(header))
        out.append(r"\midrule")
        for r in body:
            out.append(fmt(r))
        out.append(r"\bottomrule")
        out.append(r"\end{tabularx}")
        out.append(r"\end{table}")
        return "\n".join(out)

    # Narrow tables: longtable for page-spanning.
    out.append(r"\begin{center}")
    out.append(r"\begin{longtable}{@{}" + "l" * ncols + r"@{}}")
    out.append(r"\toprule")
    out.append(fmt(header))
    out.append(r"\midrule")
    out.append(r"\endhead")
    for r in body:
        out.append(fmt(r))
    out.append(r"\bottomrule")
    out.append(r"\end{longtable}")
    out.append(r"\end{center}")
    return "\n".join(out)


META_RE = re.compile(r"^\*\*([^*:]+?):\*\*\s*(.*)$")


def consume_metadata_block(lines: list[str], i: int) -> tuple[int, list[tuple[str, str]]]:
    items: list[tuple[str, str]] = []
    n = len(lines)
    while i < n:
        stripped = lines[i].strip()
        m = META_RE.match(stripped)
        if not m:
            break
        field = m.group(1).strip()
        value = m.group(2).strip()
        j = i + 1
        while j < n:
            nxt = lines[j]
            if not nxt.strip():
                break
            if META_RE.match(nxt.strip()):
                break
            if nxt.strip().startswith("#") or nxt.strip().startswith("---") \
                    or nxt.strip().startswith("|"):
                break
            if re.match(r"^\s*[-*+]\s+|^\s*\d+\.\s+", nxt):
                break
            value += " " + nxt.strip()
            j += 1
        items.append((field, value))
        i = j
    return i, items


_MONTHS = {"01":"January","02":"February","03":"March","04":"April",
           "05":"May","06":"June","07":"July","08":"August",
           "09":"September","10":"October","11":"November","12":"December"}


def _format_date(d: str) -> str:
    """YYYY-MM-DD -> 'Month D, YYYY'. Leaves other formats untouched."""
    m = re.match(r"^(\d{4})-(\d{2})-(\d{2})$", d.strip())
    if not m:
        return d
    y, mo, dd = m.group(1), m.group(2), m.group(3)
    month = _MONTHS.get(mo, mo)
    return f"{month} {int(dd)}, {y}"


_TITLE_FIELDS = {"author", "date", "session", "framework", "source",
                 "cross-links", "cross links", "affiliation", "email",
                 "contact"}


def convert(md: str, paper_id: str = "") -> str:
    lines = md.splitlines()
    meta, lines = parse_frontmatter(lines)

    # Use first H1 as title (typically "PAPER_NNN: Title").
    title = ""
    h1_idx = None
    for idx, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("# ") and not s.startswith("## "):
            title = s[2:].strip()
            h1_idx = idx
            break
        if s:
            break
    if not title:
        title = meta.get("title", paper_id)

    # Pull leading metadata block (after H1) BEFORE writing preamble so we
    # can use Author/Date/Affiliation/Email in \maketitle. Field values not
    # used by title (Session, Framework, Source, Cross-links) are dropped
    # because the canonical arXiv format keeps them out of \maketitle.
    start = (h1_idx + 1) if h1_idx is not None else 0
    while start < len(lines) and not lines[start].strip():
        start += 1
    pre_meta: dict[str, str] = {}
    if start < len(lines) and META_RE.match(lines[start].strip()):
        new_start, items = consume_metadata_block(lines, start)
        for k, v in items:
            pre_meta[k.strip().lower()] = v.strip()
        start = new_start

    author = pre_meta.get("author") or meta.get("author", "Daniel T. Murphy")
    affiliation = pre_meta.get("affiliation") or meta.get("affiliation",
                                                          "Independent Research")
    email = pre_meta.get("email") or meta.get("email",
                                              "daniel8murphy0007@github.com")
    # Date MUST come from frontmatter only (stable, single source of truth).
    # Do NOT use body-level "**Date:**" metadata which can drift.
    date = meta.get("date") or pre_meta.get("date", "")
    date_fmt = _format_date(date) if date else ""

    # Build two-line bold title: split on first ': ' if present.
    if ": " in title:
        prefix, _, rest = title.partition(": ")
        title_tex = (r"\textbf{" + render_inline(prefix) + r": \\ "
                     + render_inline(rest) + r"}")
    else:
        title_tex = r"\textbf{" + render_inline(title) + r"}"

    out: list[str] = [PREAMBLE]
    out.append(r"\title{" + title_tex + r"}")
    # Author block with superscript affiliation + email (arXiv canonical).
    out.append(r"\author{")
    out.append(escape_text(author) + r"\textsuperscript{1} \\")
    out.append(r"\small \textsuperscript{1}" + escape_text(affiliation) + r" \\")
    out.append(r"\small \texttt{" + escape_text(email) + r"}")
    out.append(r"}")
    if date_fmt:
        out.append(r"\date{" + escape_text(date_fmt) + r"}")
    out.append(r"\begin{document}")
    out.append(r"\maketitle")
    out.append("")

    # ------------------------------------------------------------------
    # Pull References block to END of document. Everything else stays in
    # source order (session updates, calibration notes, SM anchors are
    # NOT appendices -- they are timeline updates and must appear inline).
    # ------------------------------------------------------------------
    def _is_refs_heading(h2_text: str) -> bool:
        t = h2_text.lower().strip()
        t = re.sub(r"^\d+(?:\.\d+)*\.?\s+", "", t)
        return (t.startswith("references") or t.startswith("bibliography")
                or t.startswith("works cited") or t.startswith("citations"))

    body_lines: list[str] = []
    refs_lines: list[str] = []
    in_refs = False
    n0 = len(lines)
    for idx in range(start, n0):
        ln = lines[idx]
        s = ln.strip()
        m_h2 = re.match(r"^##\s+(.+)$", s)
        if m_h2:
            in_refs = _is_refs_heading(m_h2.group(1))
            # Skip the references heading itself; we synthesize sections at end
            if in_refs:
                continue
        elif in_refs and re.match(r"^#{1,6}\s+", s):
            # any further heading inside refs zone ends refs collection
            in_refs = False
            body_lines.append(ln)
            continue
        if in_refs:
            refs_lines.append(ln)
        else:
            body_lines.append(ln)

    lines = body_lines
    n = len(lines)
    i = 0

    para: list[tuple[str, bool]] = []

    def flush_para():
        if not para:
            return
        # Join paragraph FIRST so multi-line inline math/code/links work,
        # then render once. Hard-breaks become literal " \\ " inside the
        # joined string and survive (the placeholder protection in
        # render_inline ignores backslash-only tokens).
        joined_parts = []
        for idx, (txt, hb) in enumerate(para):
            joined_parts.append(txt)
            if hb and idx < len(para) - 1:
                joined_parts.append("ZZBREAKZZ")
        joined = " ".join(joined_parts)

        # Auto-promote a standalone single-$ math expression to a numbered
        # display equation. Body text often uses "$expr$" alone on a line
        # as a pseudo-equation; render it as proper \begin{equation}.
        stripped_j = joined.strip()
        m_disp = re.fullmatch(r"\$([^$].*?[^$])\$", stripped_j)
        if m_disp and "$" not in m_disp.group(1):
            out.append(r"\begin{equation}")
            out.append(m_disp.group(1))
            out.append(r"\end{equation}")
            out.append("")
            para.clear()
            return

        rendered = render_inline(joined)
        rendered = rendered.replace("ZZBREAKZZ", r"\\ ")
        out.append(rendered.strip())
        out.append("")
        para.clear()

    n = len(lines)
    while i < n:
        raw = lines[i]
        stripped = raw.strip()

        if not stripped:
            flush_para()
            i += 1
            continue

        if re.fullmatch(r"-{3,}|\*{3,}|_{3,}", stripped):
            flush_para()
            out.append(r"\medskip\hrule\medskip")
            out.append("")
            i += 1
            continue

        m = re.match(r"^```(\w*)\s*$", raw)
        if m:
            flush_para()
            i += 1
            code = []
            while i < n and not re.match(r"^```\s*$", lines[i]):
                code.append(lines[i])
                i += 1
            i += 1
            out.append(r"\begin{verbatim}")
            out.extend(sanitize_verbatim(c) for c in code)
            out.append(r"\end{verbatim}")
            out.append("")
            continue

        if stripped.startswith("$$"):
            flush_para()
            if stripped.count("$$") >= 2 and len(stripped) > 4:
                content = stripped[2:].rsplit("$$", 1)[0]
                out.append(r"\begin{equation}")
                out.append(map_unicode_math(content))
                out.append(r"\end{equation}")
                out.append("")
                i += 1
                continue
            first = stripped[2:]
            buf = [map_unicode_math(first)] if first.strip() else []
            i += 1
            while i < n:
                ln = lines[i]
                if "$$" in ln:
                    pre = ln.split("$$", 1)[0]
                    if pre.strip():
                        buf.append(map_unicode_math(pre))
                    i += 1
                    break
                # Skip blank lines inside math env (LaTeX forbids them).
                if ln.strip():
                    buf.append(map_unicode_math(ln))
                i += 1
            out.append(r"\begin{equation}")
            out.extend(buf)
            out.append(r"\end{equation}")
            out.append("")
            continue

        m = re.match(r"^(#{1,6})\s+(.*)$", stripped)
        if m:
            flush_para()
            level = len(m.group(1))
            text = re.sub(r"\s*#+\s*$", "", m.group(2).strip())
            # Strip leading numeric prefix like "1. " / "1.2 " so we don't
            # double-number when using numbered \section / \subsection.
            text_clean = re.sub(r"^\d+(?:\.\d+)*\.?\s+", "", text)
            r = render_inline(text_clean)

            # H2 "Abstract" -> abstract environment. Slurp until next H2 and
            # render as one block, lifting any trailing "**Keywords:** ..."
            # line into the abstract per arXiv convention.
            if level == 2 and text_clean.strip().lower() == "abstract":
                i += 1
                abs_lines: list[str] = []
                while i < n:
                    nxt = lines[i].strip()
                    if re.match(r"^#{1,2}\s+", nxt):
                        break
                    if re.fullmatch(r"-{3,}|\*{3,}|_{3,}", nxt):
                        break
                    abs_lines.append(lines[i])
                    i += 1
                # Find optional Keywords line at end
                keywords = ""
                kept: list[str] = []
                for ln in abs_lines:
                    mk = re.match(r"^\s*\*\*Keywords?:\*\*\s*(.+)$", ln,
                                  re.IGNORECASE)
                    if mk:
                        keywords = mk.group(1).strip()
                    else:
                        kept.append(ln)
                abs_text = render_inline(" ".join(
                    l.strip() for l in kept if l.strip()))
                out.append(r"\begin{abstract}")
                out.append(abs_text)
                if keywords:
                    out.append("")
                    out.append(r"\textbf{Keywords:} " + render_inline(keywords))
                out.append(r"\end{abstract}")
                out.append("")
                continue

            if level == 1:
                out.append(r"\section{" + r + r"}")
            elif level == 2:
                out.append(r"\section{" + r + r"}")
            elif level == 3:
                out.append(r"\subsection{" + r + r"}")
            elif level == 4:
                out.append(r"\subsubsection{" + r + r"}")
            else:
                out.append(r"\paragraph{" + r + r"}")
            out.append("")
            i += 1
            continue

        if stripped.startswith("|") and i + 1 < n \
                and re.match(r"^\s*\|?[\s\-:|]+\|?\s*$", lines[i + 1]) \
                and "-" in lines[i + 1]:
            flush_para()
            tbl = [raw, lines[i + 1]]
            i += 2
            while i < n and lines[i].strip().startswith("|"):
                tbl.append(lines[i])
                i += 1
            out.append(render_table(tbl))
            out.append("")
            continue

        if META_RE.match(stripped) and not para:
            j, items = consume_metadata_block(lines, i)
            if len(items) >= 2:
                out.append(r"\begin{flushleft}")
                for k, v in items:
                    out.append(r"\textbf{" + escape_text(k) + r":} "
                               + render_inline(v) + r" \\")
                out.append(r"\end{flushleft}")
                out.append("")
                i = j
                continue

        if re.match(r"^\s*[-*+]\s+", raw):
            flush_para()
            out.append(r"\begin{itemize}")
            while i < n and re.match(r"^\s*[-*+]\s+", lines[i]):
                item = re.sub(r"^\s*[-*+]\s+", "", lines[i])
                j = i + 1
                while j < n and lines[j].startswith("  ") and lines[j].strip() \
                        and not re.match(r"^\s*[-*+]\s+|^\s*\d+\.\s+", lines[j].lstrip()):
                    item += " " + lines[j].strip()
                    j += 1
                out.append(r"  \item " + render_inline(item))
                i = j
            out.append(r"\end{itemize}")
            out.append("")
            continue

        if re.match(r"^\s*\d+\.\s+", raw):
            flush_para()
            out.append(r"\begin{enumerate}")
            while i < n and re.match(r"^\s*\d+\.\s+", lines[i]):
                item = re.sub(r"^\s*\d+\.\s+", "", lines[i])
                j = i + 1
                while j < n and lines[j].startswith("  ") and lines[j].strip() \
                        and not re.match(r"^\s*[-*+]\s+|^\s*\d+\.\s+", lines[j].lstrip()):
                    item += " " + lines[j].strip()
                    j += 1
                out.append(r"  \item " + render_inline(item))
                i = j
            out.append(r"\end{enumerate}")
            out.append("")
            continue

        if stripped.startswith(">"):
            flush_para()
            out.append(r"\begin{quote}")
            while i < n and lines[i].strip().startswith(">"):
                out.append(render_inline(lines[i].strip().lstrip(">").lstrip()))
                i += 1
            out.append(r"\end{quote}")
            out.append("")
            continue

        # Plain paragraph line with hard-break detection
        hardbreak = raw.rstrip("\n").endswith("  ")
        para.append((stripped, hardbreak))
        i += 1

    flush_para()
    out.append("")

    # ------------------------------------------------------------------
    # Split References block into TWO end sections:
    #   - \section*{Citations} = external validation (arXiv, DOI, NASA, papers)
    #   - \section*{References} = code base (.py, .cpp, .h, scripts, repo files)
    # Per author convention: "References = code base, Citations = external."
    # ------------------------------------------------------------------
    if refs_lines:
        raw_items: list[tuple[str, str]] = []
        for ln in refs_lines:
            s = ln.strip()
            if not s:
                continue
            if s.startswith("#") or s.startswith("<!--"):
                continue
            if re.fullmatch(r"-{3,}|\*{3,}|_{3,}", s):
                continue
            m_ref = re.match(r"^(\d+[a-zA-Z]?)\.\s+(.+)$", s)
            if m_ref:
                raw_items.append((m_ref.group(1), m_ref.group(2)))
            elif raw_items and (ln.startswith(" ") or ln.startswith("\t")):
                k, v = raw_items[-1]
                raw_items[-1] = (k, v + " " + s)

        # Classify each entry: code file refs vs external citations.
        code_pat = re.compile(
            r"\.(py|cpp|cxx|cc|c|h|hpp|js|ts|ps1|sh|bat|cmake|json|csv|md|txt|yaml|yml)\b"
            r"|MAIN_\{?1|source\d+|repository|namespace|module",
            re.IGNORECASE)
        ext_pat = re.compile(
            r"arXiv:|doi:|Phys\.|ApJ|ApJL|MNRAS|Nature|Science|NASA|ESA|LIGO|"
            r"Virgo|Observatory|Collaboration|et al\.|https?://",
            re.IGNORECASE)

        citations: list[tuple[str, str]] = []
        references: list[tuple[str, str]] = []
        for k, v in raw_items:
            is_code = bool(code_pat.search(v))
            is_ext = bool(ext_pat.search(v))
            # Code-file beats external if both match (the file itself wins).
            if is_code and not is_ext:
                references.append((k, v))
            elif is_ext and not is_code:
                citations.append((k, v))
            elif is_code:
                references.append((k, v))
            else:
                citations.append((k, v))

        if citations:
            out.append(r"\section*{Citations}")
            out.append(r"\addcontentsline{toc}{section}{Citations}")
            out.append(r"\begin{thebibliography}{99}")
            for k, v in citations:
                out.append(r"\bibitem{cite" + k + "} " + render_inline(v))
            out.append(r"\end{thebibliography}")
            out.append("")

        if references:
            out.append(r"\section*{References}")
            out.append(r"\addcontentsline{toc}{section}{References}")
            out.append(r"\begin{itemize}")
            for k, v in references:
                out.append(r"  \item[" + k + r".] " + render_inline(v))
            out.append(r"\end{itemize}")
            out.append("")

    out.append(TYPEOUT)
    out.append(r"\end{document}")
    return "\n".join(out) + "\n"


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    for path in argv[1:]:
        p = Path(path)
        if not p.exists():
            print(f"SKIP (missing): {p}")
            continue
        md = p.read_text(encoding="utf-8")
        tex = convert(md, paper_id=p.stem)
        out_path = p.with_suffix(".tex")
        out_path.write_text(tex, encoding="utf-8")
        print(f"OK: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
