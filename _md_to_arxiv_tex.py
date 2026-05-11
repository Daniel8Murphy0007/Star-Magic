"""
_md_to_arxiv_tex.py
Star-Magic UQFF whitepapers: convert markdown -> arXiv-canonical LaTeX.

NO pandoc. NO xelatex. NO reportlab. Direct text->LaTeX with pdflatex output.
Template matches whitepapers/PAPER_495_Cosmic_Quantum_Egg_Theory.tex and
Manuscript 1_12Feb2026/uqff_production_arxiv.tex.

Usage:
    python _md_to_arxiv_tex.py whitepapers/PAPER_NNN_*.md [...]
Outputs PAPER_NNN_*.tex alongside the .md file.
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

# ----------------------------------------------------------------------------
# LaTeX preamble (arXiv-canonical, matches PAPER_495 + uqff_production_arxiv)
# ----------------------------------------------------------------------------
PREAMBLE = r"""\documentclass[11pt]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath, amssymb, amsthm}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage{booktabs}
\usepackage{array}
\usepackage{longtable}
\usepackage[margin=1in]{geometry}
\usepackage{enumitem}
\usepackage{textcomp}
\usepackage{xcolor}
\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue,citecolor=blue}
\setlength{\parskip}{0.5em}
\setlength{\parindent}{0pt}
"""

TYPEOUT = r"\typeout{get arXiv to do 4 passes: Label(s) may have changed. Rerun}"

# ----------------------------------------------------------------------------
# Special character escaping for plain text (NOT inside math or verbatim)
# ----------------------------------------------------------------------------
LATEX_ESCAPES = [
    ("\\", r"\textbackslash{}"),  # must be first
    ("&", r"\&"),
    ("%", r"\%"),
    ("$", r"\$"),
    ("#", r"\#"),
    ("_", r"\_"),
    ("{", r"\{"),
    ("}", r"\}"),
    ("~", r"\textasciitilde{}"),
    ("^", r"\textasciicircum{}"),
]

# Unicode -> LaTeX command map (common physics/Greek)
UNICODE_MAP = {
    "α": r"$\alpha$", "β": r"$\beta$", "γ": r"$\gamma$", "δ": r"$\delta$",
    "ε": r"$\epsilon$", "ζ": r"$\zeta$", "η": r"$\eta$", "θ": r"$\theta$",
    "ι": r"$\iota$", "κ": r"$\kappa$", "λ": r"$\lambda$", "μ": r"$\mu$",
    "ν": r"$\nu$", "ξ": r"$\xi$", "ο": r"o", "π": r"$\pi$",
    "ρ": r"$\rho$", "σ": r"$\sigma$", "τ": r"$\tau$", "υ": r"$\upsilon$",
    "φ": r"$\phi$", "χ": r"$\chi$", "ψ": r"$\psi$", "ω": r"$\omega$",
    "Α": r"A", "Β": r"B", "Γ": r"$\Gamma$", "Δ": r"$\Delta$",
    "Ε": r"E", "Ζ": r"Z", "Η": r"H", "Θ": r"$\Theta$",
    "Ι": r"I", "Κ": r"K", "Λ": r"$\Lambda$", "Μ": r"M",
    "Ν": r"N", "Ξ": r"$\Xi$", "Ο": r"O", "Π": r"$\Pi$",
    "Ρ": r"P", "Σ": r"$\Sigma$", "Τ": r"T", "Υ": r"$\Upsilon$",
    "Φ": r"$\Phi$", "Χ": r"X", "Ψ": r"$\Psi$", "Ω": r"$\Omega$",
    "∝": r"$\propto$", "∓": r"$\mp$", "±": r"$\pm$",
    "≳": r"$\gtrsim$", "≲": r"$\lesssim$", "≈": r"$\approx$",
    "≠": r"$\neq$", "≤": r"$\leq$", "≥": r"$\geq$", "≡": r"$\equiv$",
    "∞": r"$\infty$", "∂": r"$\partial$", "∇": r"$\nabla$",
    "∈": r"$\in$", "∉": r"$\notin$", "⊂": r"$\subset$", "⊃": r"$\supset$",
    "∑": r"$\sum$", "∏": r"$\prod$", "∫": r"$\int$",
    "→": r"$\to$", "←": r"$\leftarrow$", "↔": r"$\leftrightarrow$",
    "⇒": r"$\Rightarrow$", "⇐": r"$\Leftarrow$", "⇔": r"$\Leftrightarrow$",
    "·": r"$\cdot$", "×": r"$\times$", "÷": r"$\div$",
    "°": r"$^\circ$", "′": r"'", "″": r"''",
    "—": r"---", "–": r"--", "…": r"\ldots{}",
    "“": r"``", "”": r"''", "‘": r"`", "’": r"'",
    "•": r"\textbullet{}", "·": r"$\cdot$",
    "⊙": r"$\odot$", "⊕": r"$\oplus$", "⊗": r"$\otimes$",
    "ℏ": r"$\hbar$", "ℓ": r"$\ell$", "ℝ": r"$\mathbb{R}$", "ℕ": r"$\mathbb{N}$",
    "ℤ": r"$\mathbb{Z}$", "ℂ": r"$\mathbb{C}$", "ℚ": r"$\mathbb{Q}$",
    "√": r"$\sqrt{}$", "∼": r"$\sim$", "≅": r"$\cong$",
    "²": r"$^2$", "³": r"$^3$", "¹": r"$^1$", "⁰": r"$^0$",
    "⁴": r"$^4$", "⁵": r"$^5$", "⁶": r"$^6$", "⁷": r"$^7$",
    "⁸": r"$^8$", "⁹": r"$^9$",
    "₀": r"$_0$", "₁": r"$_1$", "₂": r"$_2$", "₃": r"$_3$",
    "₄": r"$_4$", "₅": r"$_5$", "₆": r"$_6$", "₇": r"$_7$",
    "₈": r"$_8$", "₉": r"$_9$",
    "Å": r"\AA{}", "ø": r"\o{}", "Ø": r"\O{}",
    "ß": r"\ss{}", "£": r"\pounds{}", "¥": r"Y",
    "˜": r"$\sim$",
    "✅": r"$\checkmark$", "❌": r"$\times$", "✓": r"$\checkmark$",
    "⚠": r"!", "★": r"$\star$", "☆": r"$\star$",
    "\ufeff": r"",  # BOM
    "©": r"\textcopyright{}", "®": r"\textregistered{}", "™": r"\texttrademark{}",
    "§": r"\S{}", "¶": r"\P{}", "†": r"\dag{}", "‡": r"\ddag{}",
    "\u00a0": r" ",  # non-breaking space
}


def escape_text(s: str) -> str:
    """Escape LaTeX specials in plain text segments only."""
    # Apply unicode mappings first
    for u, repl in UNICODE_MAP.items():
        s = s.replace(u, repl)
    # Then escape LaTeX specials
    for ch, repl in LATEX_ESCAPES:
        s = s.replace(ch, repl)
    return s


# ----------------------------------------------------------------------------
# Tokenizer: split line into (kind, text) segments where kind is text/math/code
# ----------------------------------------------------------------------------
INLINE_RE = re.compile(
    r"(\$[^\$\n]+\$)"         # inline math $...$
    r"|(`[^`\n]+`)"           # inline code `...`
)


def render_inline(line: str) -> str:
    """Process inline markdown: bold, italic, links, escapes (math/code preserved)."""
    parts = []
    last = 0
    for m in INLINE_RE.finditer(line):
        if m.start() > last:
            parts.append(("text", line[last:m.start()]))
        if m.group(1):
            parts.append(("math", m.group(1)))
        elif m.group(2):
            parts.append(("code", m.group(2)[1:-1]))
        last = m.end()
    if last < len(line):
        parts.append(("text", line[last:]))

    out = []
    for kind, txt in parts:
        if kind == "math":
            out.append(txt)
        elif kind == "code":
            # \verb|...|  pick a delimiter not in text
            for d in "|!@^+":
                if d not in txt:
                    out.append(f"\\verb{d}{txt}{d}")
                    break
            else:
                out.append("\\texttt{" + escape_text(txt) + "}")
        else:
            # text segment: process markdown bold/italic/links, then escape
            t = txt
            # links [text](url)
            t = re.sub(
                r"\[([^\]]+)\]\(([^)]+)\)",
                lambda m: r"\href{" + m.group(2).replace("\\", "/").replace("%", r"\%") + "}{" + escape_text(m.group(1)) + "}",
                t,
            )
            # bold **x** and __x__
            placeholders = {}
            def stash(s):
                k = f"\x00BOLD{len(placeholders)}\x00"
                placeholders[k] = s
                return k
            t = re.sub(r"\*\*([^\*\n]+)\*\*", lambda m: stash(r"\textbf{" + escape_text(m.group(1)) + "}"), t)
            t = re.sub(r"__([^_\n]+)__", lambda m: stash(r"\textbf{" + escape_text(m.group(1)) + "}"), t)
            # italic *x* (avoid matching ** leftovers) and _x_ (be cautious)
            t = re.sub(r"(?<![\*\w])\*([^\*\n]+)\*(?!\*)", lambda m: stash(r"\emph{" + escape_text(m.group(1)) + "}"), t)
            # Now escape the remainder
            t = escape_text(t)
            for k, v in placeholders.items():
                t = t.replace(escape_text(k), v).replace(k, v)
            out.append(t)
    return "".join(out)


# ----------------------------------------------------------------------------
# Block-level parser
# ----------------------------------------------------------------------------
def parse_table(rows: list[str]) -> str:
    """Convert pipe table to longtable."""
    # rows[0]: header; rows[1]: separator; rows[2:]: data
    def cells(line):
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
    out = []
    out.append(r"\begin{center}")
    out.append(r"\begin{longtable}{" + "l" * ncols + "}")
    out.append(r"\toprule")
    def fmt_row(cs):
        rendered = [render_inline(x) for x in cs]
        # Protect cells that start with '[' so following row's '[' doesn't become
        # an optional arg to \tabularnewline.
        rendered = [("{}" + r) if r.lstrip().startswith("[") else r for r in rendered]
        return " & ".join(rendered) + r" \tabularnewline"
    out.append(fmt_row(header))
    out.append(r"\midrule")
    out.append(r"\endhead")
    for r in rows[2:]:
        c = cells(r)
        if len(c) < ncols:
            c += [""] * (ncols - len(c))
        elif len(c) > ncols:
            c = c[:ncols]
        out.append(fmt_row(c))
    out.append(r"\bottomrule")
    out.append(r"\end{longtable}")
    out.append(r"\end{center}")
    return "\n".join(out)


def parse_frontmatter(lines: list[str]) -> tuple[dict, list[str]]:
    """Extract YAML frontmatter if present."""
    meta = {}
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


def convert(md: str, paper_id: str = "") -> str:
    lines = md.splitlines()
    meta, lines = parse_frontmatter(lines)

    title = meta.get("title", paper_id or "Untitled")
    author = meta.get("author", "Daniel T. Murphy")
    date = meta.get("date", "")

    out = [PREAMBLE]
    out.append(r"\title{" + escape_text(title) + "}")
    out.append(r"\author{" + escape_text(author) + "}")
    if date:
        out.append(r"\date{" + escape_text(date) + "}")
    out.append(r"\begin{document}")
    out.append(r"\maketitle")

    i = 0
    n = len(lines)
    in_para = []

    def flush_para():
        if in_para:
            txt = " ".join(in_para).strip()
            if txt:
                out.append(render_inline(txt))
                out.append("")
            in_para.clear()

    # Strip leading H1 if it duplicates the title
    # (papers typically begin with "# PAPER_NNN: ...")
    first_h1_skipped = False

    while i < n:
        line = lines[i]
        stripped = line.strip()

        # Blank line ends paragraph
        if not stripped:
            flush_para()
            i += 1
            continue

        # Horizontal rule
        if re.fullmatch(r"-{3,}|\*{3,}|_{3,}", stripped):
            flush_para()
            out.append(r"\noindent\rule{\textwidth}{0.4pt}")
            out.append("")
            i += 1
            continue

        # Fenced code block ```...```
        m = re.match(r"^```(\w*)\s*$", line)
        if m:
            flush_para()
            lang = m.group(1)
            i += 1
            code_lines = []
            while i < n and not re.match(r"^```\s*$", lines[i]):
                code_lines.append(lines[i])
                i += 1
            i += 1  # skip closing ```
            out.append(r"\begin{verbatim}")
            out.extend(code_lines)
            out.append(r"\end{verbatim}")
            out.append("")
            continue

        # Display math $$...$$
        if stripped.startswith("$$"):
            flush_para()
            # collect until closing $$
            buf = [line]
            # Check if it closes on same line: $$...$$
            if stripped.count("$$") >= 2 and len(stripped) > 4:
                out.append(line)
                out.append("")
                i += 1
                continue
            i += 1
            while i < n:
                buf.append(lines[i])
                if "$$" in lines[i]:
                    i += 1
                    break
                i += 1
            out.append("\n".join(buf))
            out.append("")
            continue

        # Headings
        m = re.match(r"^(#{1,6})\s+(.*)$", stripped)
        if m:
            flush_para()
            level = len(m.group(1))
            text = m.group(2).strip()
            # Strip trailing trailing punctuation like ## Foo ##
            text = re.sub(r"\s*#+\s*$", "", text)
            rendered = render_inline(text)
            if level == 1:
                if not first_h1_skipped:
                    # First H1 is title duplicate; skip
                    first_h1_skipped = True
                    i += 1
                    continue
                out.append(r"\section*{" + rendered + "}")
            elif level == 2:
                out.append(r"\section*{" + rendered + "}")
            elif level == 3:
                out.append(r"\subsection*{" + rendered + "}")
            elif level == 4:
                out.append(r"\subsubsection*{" + rendered + "}")
            else:
                out.append(r"\paragraph{" + rendered + "}")
            out.append("")
            i += 1
            continue

        # Table block: line starts with | and next line is separator
        if stripped.startswith("|") and i + 1 < n and re.match(r"^\s*\|?[\s\-:|]+\|?\s*$", lines[i + 1]) and "-" in lines[i + 1]:
            flush_para()
            tbl_lines = [line, lines[i + 1]]
            i += 2
            while i < n and lines[i].strip().startswith("|"):
                tbl_lines.append(lines[i])
                i += 1
            out.append(parse_table(tbl_lines))
            out.append("")
            continue

        # Bullet list
        if re.match(r"^\s*[\-\*\+]\s+", line):
            flush_para()
            out.append(r"\begin{itemize}")
            while i < n and re.match(r"^\s*[\-\*\+]\s+", lines[i]):
                item = re.sub(r"^\s*[\-\*\+]\s+", "", lines[i])
                # continuation lines (indented further)
                j = i + 1
                while j < n and lines[j].startswith("  ") and lines[j].strip() and not re.match(r"^\s*[\-\*\+\d]", lines[j].lstrip()):
                    item += " " + lines[j].strip()
                    j += 1
                out.append(r"  \item " + render_inline(item))
                i = j
            out.append(r"\end{itemize}")
            out.append("")
            continue

        # Numbered list
        if re.match(r"^\s*\d+\.\s+", line):
            flush_para()
            out.append(r"\begin{enumerate}")
            while i < n and re.match(r"^\s*\d+\.\s+", lines[i]):
                item = re.sub(r"^\s*\d+\.\s+", "", lines[i])
                j = i + 1
                while j < n and lines[j].startswith("  ") and lines[j].strip() and not re.match(r"^\s*[\-\*\+\d]", lines[j].lstrip()):
                    item += " " + lines[j].strip()
                    j += 1
                out.append(r"  \item " + render_inline(item))
                i = j
            out.append(r"\end{enumerate}")
            out.append("")
            continue

        # Blockquote
        if stripped.startswith(">"):
            flush_para()
            out.append(r"\begin{quote}")
            while i < n and lines[i].strip().startswith(">"):
                out.append(render_inline(lines[i].strip().lstrip(">").lstrip()))
                i += 1
            out.append(r"\end{quote}")
            out.append("")
            continue

        # Default: accumulate into paragraph
        in_para.append(stripped)
        i += 1

    flush_para()
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
        paper_id = p.stem
        tex = convert(md, paper_id=paper_id)
        out = p.with_suffix(".tex")
        out.write_text(tex, encoding="utf-8")
        print(f"OK: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
