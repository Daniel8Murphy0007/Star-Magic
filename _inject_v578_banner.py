"""
_inject_v578_banner.py - Post-process .tex files to inject parser-locked
v5.78 closure-banner comments (% CLOSURE / % TEMPLATE / % CANONICAL) on
three separate lines immediately before \end{document}.

Banner schema lives in /memories/repo/v5_78_templates.md. The MD-to-TeX
converter (_md_to_arxiv_tex.py) does not preserve % LaTeX comments, so
banners are injected here after conversion.

Usage:
    python _inject_v578_banner.py <tex_file> "<closure_label>" \
        "<template_code>" "<canonical_oneline>"

Example:
    python _inject_v578_banner.py whitepapers\PAPER_044_X.tex \
        "PreBigBang_26Center_DPM_Manifold" "T-Lambda" \
        "rho_SCm=7.09e-37 J/m^3 (PAPER_1170); ..."

For ledger papers (no specific numeric prediction):
    predicted=ledger-derived observed=ledger-derived error_pct=0.000
"""
from __future__ import annotations
import sys
from pathlib import Path


def inject(tex_path: Path, label: str, template: str, canonical: str,
           predicted: str = "ledger-derived",
           observed: str = "ledger-derived",
           error_pct: str = "0.000") -> int:
    if not tex_path.exists():
        print(f"ERROR: {tex_path} not found")
        return 1
    src = tex_path.read_text(encoding="utf-8")

    banner = (
        f"% CLOSURE :: {label} :: predicted={predicted} observed={observed} "
        f"error_pct={error_pct}\n"
        f"% TEMPLATE :: {template}\n"
        f"% CANONICAL :: {canonical}\n"
    )

    end_marker = r"\end{document}"
    if end_marker not in src:
        print(f"ERROR: {end_marker} missing in {tex_path}")
        return 1

    # Idempotency: strip any prior CLOSURE/TEMPLATE/CANONICAL comment block
    lines = src.splitlines(keepends=True)
    cleaned = [ln for ln in lines
               if not (ln.startswith("% CLOSURE ::")
                       or ln.startswith("% TEMPLATE ::")
                       or ln.startswith("% CANONICAL ::"))]
    src = "".join(cleaned)

    src = src.replace(end_marker, banner + end_marker, 1)
    tex_path.write_text(src, encoding="utf-8")
    print(f"OK: injected v5.78 banner -> {tex_path}")
    return 0


def main(argv):
    if len(argv) < 5:
        print(__doc__)
        return 1
    return inject(Path(argv[1]), argv[2], argv[3], argv[4],
                  predicted=argv[5] if len(argv) > 5 else "ledger-derived",
                  observed=argv[6] if len(argv) > 6 else "ledger-derived",
                  error_pct=argv[7] if len(argv) > 7 else "0.000")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
