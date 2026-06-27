# scm_latex_exporter.py
# Restored 2026-06-26 from pre-9c1c7083 pdf/scm_latex_exporter.py.
# Commit 9c1c7083 deleted pdf/scm_latex_exporter.py but did NOT move it to root.
# This is the ONLY pdf/.py file from that commit without a root counterpart.
# OUT_DIR adjusted to "pdf/" so generated .tex/.pdf still land in pdf/ folder.
from scm_vacuum_manifold import export_all_to_latex, monte_carlo_fubi_i, progress_metric
import subprocess
import os

# Output to pdf/ folder (relative to repo root where this script lives now)
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pdf")
if not os.path.isdir(OUT_DIR):
    os.makedirs(OUT_DIR, exist_ok=True)

def generate_whitepaper_section():
    latex = export_all_to_latex()
    mean, std, rng = monte_carlo_fubi_i()
    mc_line = (
        f"$\\langle F_{{U,Bi,i}} \\rangle = {mean:.3e}$ N, "
        f"$\\sigma = {std:.3e}$ N, "
        f"90\\%\\ CI $= [{rng[0]:.3e},\\ {rng[1]:.3e}]$ N. "
        f"Progress: {progress_metric}\\%."
    )

    content = r"""
\documentclass{article}
\usepackage{amsmath, amssymb}
\title{SCm Vacuum Manifold --- Primordial First-Principle Derivation}
\author{Daniel T. Murphy --- Clean Thread 27FEB2026\_A.docx}
\date{\today}

\begin{document}
\maketitle

\section{SCm Vacuum Manifold (Pre-Gravity First Principle)}
The primordial ``matter'' is the superconducting vacuum manifold:
\[
\rho_{\text{vac,SCm}} = """ + latex['rho_scm'] + r"""
\]

Negative-time modulation:
\[
\cos(\pi t_n) = """ + latex['cos_pi_tn'] + r"""
\]

Buoyancy force (long-form integral):
\[
""" + latex['F_U_Bi_i'] + r"""
\]

26D Vacuum Density Series:
\[
""" + latex['VDS'] + r"""
\]

Phonon activation at 1.25 THz:
\[
\omega_{SCm} = 1.25\ \text{THz}
\]

Holmlid Rydberg bridge:
\[
""" + latex.get('Holmlid', r'E_{Holmlid} = 630\ \text{eV}') + r"""
\]

\section{Numerical Validation (Monte-Carlo on F\_U\_Bi\_i)}
""" + mc_line + r"""

\section{Refined Proof Fragment Outline (Riemann + Cosmogenesis)}
\begin{enumerate}
\item SCm Vacuum Manifold (1 page)
\item Primordial Split \& Cosmic Quantum Egg (0.5 page)
\item Holmlid Rydberg Bridge (0.5 page)
\item Riemann Closure via SCm (1 page)
\item Numerical Validation (0.5--1 page)
\end{enumerate}

\end{document}
"""
    tex_path = os.path.join(OUT_DIR, "SCm_Whitepaper_Section.tex")
    with open(tex_path, "w", encoding="utf-8") as f:
        f.write(content)

    # Compile to PDF inside pdf/ so output lands there
    subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-output-directory", OUT_DIR,
         "SCm_Whitepaper_Section.tex"],
        cwd=OUT_DIR
    )
    pdf_path = os.path.join(OUT_DIR, "SCm_Whitepaper_Section.pdf")
    if os.path.exists(pdf_path):
        print(f"[OK] PDF ready: {pdf_path}")
    else:
        print("[FAIL] PDF generation failed -- check SCm_Whitepaper_Section.log")

if __name__ == "__main__":
    generate_whitepaper_section()
