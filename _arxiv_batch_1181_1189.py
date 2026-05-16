"""Non-interactive arXiv bundle builder for PAPER_1181-1189.

Per-paper steps:
  1. Copy .tex (or convert .md -> .tex via pandoc if no .tex exists)
  2. Strip LaTeX comments
  3. pdflatex pass x2 (in-place)
  4. tar -> PAPER_NNNN.tar inside per-paper subdir

Final step: zip all per-paper folders + PDFs into arxiv_submission_1181_1189.zip
"""
import os, re, shutil, subprocess, sys, tarfile, zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
WP   = ROOT / "whitepapers"
PDF  = ROOT / "pdf"
BUNDLE = ROOT / "arxiv_submission_1181_1189"
PDFLATEX = r"C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
PANDOC   = Path(os.environ.get("LOCALAPPDATA",""), "Pandoc", "pandoc.exe")

PAPERS = {
    1181: "PAPER_1181_UQFF_Grand_Unification_S266_S295_Thirty_Closures",
    1182: "PAPER_1182_UQFF_Millennium_Prize_Unified_Proof_Set",
    1183: "PAPER_1183_UQFF_Aggressive_Paradox_Unified_Proof_Set",
    1184: "PAPER_1184_UQFF_Open_Problems_Unified_Proof_Set",
    1185: "PAPER_1185_UQFF_Quantum_Gravity_Unified_Proof_Set",
    1186: "PAPER_1186_UQFF_Standard_Model_Unified_Proof_Set",
    1187: "PAPER_1187_UQFF_Cosmological_Tensions_Unified_Proof_Set",
    1188: "PAPER_1188_UQFF_Number_Theory_Frontier_Set",
    1189: "PAPER_1189_UQFF_Chemistry_Atomic_Unified_Proof_Set",
}

# Reset bundle dir
if BUNDLE.exists(): shutil.rmtree(BUNDLE)
BUNDLE.mkdir()

def strip_tex_comments(tex_path: Path):
    out_lines = []
    for ln in tex_path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True):
        if "%" in ln:
            parts = re.split(r"(?<!\\)%", ln, 1)
            cleaned = parts[0].rstrip() + "\n"
            if cleaned.strip(): out_lines.append(cleaned)
        else:
            out_lines.append(ln)
    tex_path.write_text("".join(out_lines), encoding="utf-8")

def md_to_tex(md_path: Path, tex_path: Path):
    if not PANDOC.exists():
        raise FileNotFoundError(f"pandoc not found at {PANDOC}")
    subprocess.run([str(PANDOC), str(md_path), "-s", "-o", str(tex_path),
                    "--from=markdown", "--to=latex"], check=True)

results = []
for pid, stem in PAPERS.items():
    print(f"\n=== PAPER_{pid} :: {stem} ===")
    pdir = BUNDLE / f"PAPER_{pid}"
    pdir.mkdir()

    src_tex = WP / f"{stem}.tex"
    src_md  = WP / f"{stem}.md"
    dst_tex = pdir / "main.tex"

    if src_tex.exists():
        shutil.copy2(src_tex, dst_tex)
        print(f"  copied .tex ({src_tex.stat().st_size} bytes)")
    elif src_md.exists():
        try:
            md_to_tex(src_md, dst_tex)
            print(f"  converted .md -> .tex via pandoc")
        except Exception as e:
            print(f"  PANDOC FAIL: {e}")
            results.append((pid, "PANDOC_FAIL", str(e)))
            continue
    else:
        print(f"  SKIP - neither .tex nor .md present")
        results.append((pid, "MISSING_SOURCE", ""))
        continue

    strip_tex_comments(dst_tex)

    # Compile twice for refs
    ok = True
    for pass_idx in (1, 2):
        try:
            r = subprocess.run([PDFLATEX, "-interaction=nonstopmode", "-halt-on-error", "main.tex"],
                               cwd=str(pdir), capture_output=True, text=True, timeout=120)
            if r.returncode != 0 and pass_idx == 2:
                print(f"  pdflatex pass {pass_idx} returncode={r.returncode}")
                ok = False
        except Exception as e:
            print(f"  pdflatex EXCEPTION pass {pass_idx}: {e}")
            ok = False; break

    pdf_out = pdir / "main.pdf"
    if pdf_out.exists():
        # Rename to canonical paper name
        canon_pdf = pdir / f"PAPER_{pid}.pdf"
        shutil.copy2(pdf_out, canon_pdf)
        print(f"  PDF: {canon_pdf.stat().st_size} bytes")
    else:
        print(f"  WARN: no main.pdf produced")
        ok = False

    # Clean aux files but keep main.tex + main.pdf
    for ext in (".aux", ".log", ".out", ".toc"):
        for f in pdir.glob(f"*{ext}"):
            f.unlink()

    # Tar the per-paper source
    tar_path = pdir / f"PAPER_{pid}.tar"
    with tarfile.open(tar_path, "w") as tf:
        tf.add(str(dst_tex), arcname="main.tex")
        bbl = pdir / "main.bbl"
        if bbl.exists(): tf.add(str(bbl), arcname="main.bbl")
    print(f"  TAR: {tar_path.stat().st_size} bytes")

    results.append((pid, "OK" if ok else "PDF_WARN", f"{tar_path.stat().st_size}B"))

# Zip everything
zip_path = ROOT / "arxiv_submission_1181_1189.zip"
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for pid in PAPERS:
        pdir = BUNDLE / f"PAPER_{pid}"
        if not pdir.exists(): continue
        for f in pdir.iterdir():
            if f.suffix in (".tex", ".tar", ".pdf"):
                zf.write(f, arcname=f"PAPER_{pid}/{f.name}")
print(f"\nBundle: {zip_path}  ({zip_path.stat().st_size} bytes)")

print("\n=== SUMMARY ===")
for pid, status, note in results:
    print(f"  PAPER_{pid}: {status}  {note}")
