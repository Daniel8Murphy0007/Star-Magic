"""
build_papers_599_601.py — Session 158 PDF build
Generates PDFs for PAPER_599–601 using generate_pdfs.py genpdf pattern.
Run from Star-Magic workspace root.
"""
import importlib, importlib.util, sys, os

def load_genpdf():
    spec = importlib.util.spec_from_file_location("generate_pdfs",
                                                   os.path.join(os.path.dirname(__file__),
                                                                "generate_pdfs.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

from pathlib import Path

PAPERS = [
    "whitepapers/PAPER_599_UQFF_BSD_Conjecture_Rank_Cohomology.md",
    "whitepapers/PAPER_600_UQFF_Hodge_Conjecture_Algebraic_Cycles.md",
    "whitepapers/PAPER_601_UQFF_Magnetic_Gateway_Cosmic_Flux.md",
]

def main():
    genpdf = load_genpdf()
    styles = genpdf.make_styles()
    out_dir = Path("pdf")
    out_dir.mkdir(exist_ok=True)
    ok = 0
    for md_path_str in PAPERS:
        md_path = Path(md_path_str)
        try:
            genpdf.md_file_to_pdf(md_path, out_dir, styles)
            pdf_name = md_path.stem + ".pdf"
            pdf_path = out_dir / pdf_name
            size = pdf_path.stat().st_size if pdf_path.exists() else 0
            print(f"  [OK]  {pdf_path}  ({size/1024:.1f} KB)")
            ok += 1
        except Exception as e:
            print(f"  [FAIL] {md_path}: {e}")
    print(f"\n{ok}/{len(PAPERS)} PDFs generated successfully.")
    if ok == len(PAPERS):
        print("Session 158 PDF build complete.")

if __name__ == "__main__":
    main()
