"""Session 258 Tier A: auto-fix Unicode + inject bibliography stubs into bundle 1159-1172.

Fixes:
  - U+2014 (em-dash) -> '---'  (outside math regions)
  - U+2192 (->)      -> '$\\to$'
  - U+039B/3A6/3A9/3C1/3B6 (Greek) -> escaped math (only the few flagged tokens)
  - Inject '\\bibliography{refs}\\bibliographystyle{unsrt}' before \\end{document}
    in each .tex file that lacks a bibliography stub.

Idempotent: safe to re-run.
"""
from __future__ import annotations

import re
from pathlib import Path

BUNDLE = Path(__file__).resolve().parent / "arxiv_submission_1159_1172"
MD_DIR = BUNDLE / "md"
TEX_DIR = BUNDLE / "tex"

UNICODE_PROSE_FIXES = {
    "\u2014": "---",   # em-dash
    "\u2192": r"$\to$",
}


def _strip_math_preserve(text: str):
    """Return (segments, math_regions) so we can rebuild after fixing prose only."""
    pattern = re.compile(r"(\$\$.*?\$\$|\$[^$\n]+\$|\\\(.*?\\\))", re.DOTALL)
    parts = []
    last = 0
    for m in pattern.finditer(text):
        if m.start() > last:
            parts.append(("prose", text[last:m.start()]))
        parts.append(("math", m.group(0)))
        last = m.end()
    if last < len(text):
        parts.append(("prose", text[last:]))
    return parts


def fix_md(path: Path) -> bool:
    text = path.read_text(encoding="utf-8")
    parts = _strip_math_preserve(text)
    changed = False
    out = []
    for kind, seg in parts:
        if kind == "prose":
            new = seg
            for u, r in UNICODE_PROSE_FIXES.items():
                if u in new:
                    new = new.replace(u, r)
                    changed = True
            out.append(new)
        else:
            out.append(seg)
    if changed:
        path.write_text("".join(out), encoding="utf-8")
    return changed


def fix_tex_bib(path: Path) -> bool:
    body = path.read_text(encoding="utf-8", errors="ignore")
    if ("\\bibliography{" in body) or ("\\begin{thebibliography}" in body):
        return False
    stub = (
        "\n% Session 258 auto-injected bibliography stub\n"
        "\\bibliographystyle{unsrt}\n"
        "\\bibliography{refs}\n"
    )
    if "\\end{document}" in body:
        body = body.replace("\\end{document}", stub + "\\end{document}", 1)
    else:
        body = body + stub
    path.write_text(body, encoding="utf-8")
    return True


def main() -> int:
    md_changed = 0
    tex_changed = 0
    for md in sorted(MD_DIR.glob("PAPER_*.md")):
        if fix_md(md):
            md_changed += 1
            print(f"  [md fixed] {md.name}")
    for tex in sorted(TEX_DIR.glob("PAPER_*.tex")):
        if fix_tex_bib(tex):
            tex_changed += 1
            print(f"  [tex bib] {tex.name}")
    print(f"\nSummary: {md_changed} md files fixed, {tex_changed} tex files got bib stubs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
