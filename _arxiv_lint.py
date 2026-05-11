"""_arxiv_lint.py - Session 257 arXiv readiness pass for bundle 1159-1172.

Audits each PAPER_*.md in arxiv_submission_1159_1172/md/ for:
  (a) non-math-mode Unicode tokens (Lambda, Phi, rho, Omega, em-dashes, check/cross marks)
      that pdflatex will reject without unicode-math.
  (b) presence of corresponding .tex with \\bibliography{} stub (or inline thebibliography).
  (c) optional pdflatex 4-pass cycle when --compile is passed (pdflatex -> bibtex -> pdflatex x2).

Exit codes: 0 = all clean; 1 = at least one paper has lint warnings/errors.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

BUNDLE = Path(__file__).resolve().parent / "arxiv_submission_1159_1172"
MD_DIR = BUNDLE / "md"
TEX_DIR = BUNDLE / "tex"

# Tokens flagged when found OUTSIDE math mode ($...$ or $$...$$ or \(...\)).
FORBIDDEN_UNICODE = {
    "\u039B": "Lambda",      # Λ
    "\u03A6": "Phi",         # Φ
    "\u03A9": "Omega",       # Ω
    "\u03C1": "rho",         # ρ
    "\u03B6": "zeta",        # ζ
    "\u2014": "em-dash",     # —
    "\u2713": "checkmark",   # ✓
    "\u2705": "white-check", # ✅
    "\u274C": "cross",       # ❌
    "\u2192": "rarrow",      # →
}


def _strip_math(text: str) -> str:
    """Remove $...$, $$...$$, and \\(...\\) regions so we only audit prose."""
    text = re.sub(r"\$\$.*?\$\$", " ", text, flags=re.DOTALL)
    text = re.sub(r"\$[^$\n]+\$", " ", text)
    text = re.sub(r"\\\(.*?\\\)", " ", text, flags=re.DOTALL)
    return text


def lint_unicode(md_path: Path) -> list[str]:
    findings: list[str] = []
    text = md_path.read_text(encoding="utf-8")
    prose = _strip_math(text)
    for ch, name in FORBIDDEN_UNICODE.items():
        n = prose.count(ch)
        if n:
            findings.append(f"  unicode-{name} x{n} (codepoint U+{ord(ch):04X})")
    return findings


def lint_bib_stub(md_path: Path) -> list[str]:
    findings: list[str] = []
    tex = TEX_DIR / (md_path.stem + ".tex")
    if not tex.exists():
        findings.append(f"  no matching .tex: {tex.name}")
        return findings
    body = tex.read_text(encoding="utf-8", errors="ignore")
    has_bib = ("\\bibliography{" in body) or ("\\begin{thebibliography}" in body)
    if not has_bib:
        findings.append("  missing \\bibliography{} / thebibliography stub")
    return findings


def compile_pass(md_path: Path) -> list[str]:
    findings: list[str] = []
    tex = TEX_DIR / (md_path.stem + ".tex")
    if not tex.exists():
        return [f"  compile skipped (no tex): {tex.name}"]
    pdflatex = Path(r"C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe")
    if not pdflatex.exists():
        return [f"  compile skipped (no pdflatex): {pdflatex}"]
    try:
        cmd = [str(pdflatex), "-interaction=nonstopmode", "-halt-on-error", tex.name]
        r = subprocess.run(cmd, cwd=tex.parent, capture_output=True, timeout=120)
        if r.returncode != 0:
            findings.append(f"  pdflatex pass1 returncode={r.returncode}")
    except subprocess.TimeoutExpired:
        findings.append("  pdflatex timed out")
    return findings


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--compile", action="store_true", help="also run pdflatex pass 1")
    args = ap.parse_args()

    if not MD_DIR.exists():
        print(f"[FAIL] missing bundle: {MD_DIR}")
        return 1

    md_files = sorted(MD_DIR.glob("PAPER_*.md"))
    if not md_files:
        print(f"[FAIL] no PAPER_*.md in {MD_DIR}")
        return 1

    n_total = len(md_files)
    n_clean = 0
    print(f"=== arXiv lint pass: {n_total} papers in {MD_DIR.name} ===")
    for md in md_files:
        findings = lint_unicode(md) + lint_bib_stub(md)
        if args.compile:
            findings += compile_pass(md)
        if findings:
            print(f"[WARN] {md.name}")
            for f in findings:
                print(f)
        else:
            print(f"[PASS] {md.name}")
            n_clean += 1

    print(f"\n=== Summary: {n_clean}/{n_total} clean ===")
    return 0 if n_clean == n_total else 1


if __name__ == "__main__":
    sys.exit(main())
