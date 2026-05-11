"""
_build_batch.py — Run md->tex conversion and pdflatex x2 for a list of papers.
Backs up old PDFs to pdf_backup_pandoc_<date>/, replaces with new pdflatex PDFs.
Usage:  python _build_batch.py <stem1> <stem2> ...
   where stem is filename without extension (e.g. PAPER_001_GW170817_UQFF_Damping_Analysis)
"""
from __future__ import annotations
import shutil
import subprocess
import sys
from pathlib import Path
from datetime import date

ROOT = Path(__file__).resolve().parent
WHITEPAPERS = ROOT / "whitepapers"
PDF_DIR = ROOT / "pdf"
BACKUP = ROOT / f"pdf_backup_pandoc_{date.today().isoformat()}"
PDFLATEX = r"C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
CONVERTER = ROOT / "_md_to_arxiv_tex.py"


def run(cmd, cwd=None):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, encoding="utf-8", errors="replace")


def build(stem: str) -> dict:
    md = WHITEPAPERS / f"{stem}.md"
    tex = WHITEPAPERS / f"{stem}.tex"
    pdf_src = WHITEPAPERS / f"{stem}.pdf"
    pdf_dst = PDF_DIR / f"{stem}.pdf"

    result = {"stem": stem, "ok": False, "tex_bytes": 0, "pdf_bytes": 0, "pages": 0, "error": ""}

    if not md.exists():
        result["error"] = f"missing md: {md}"
        return result

    # 1. Convert md -> tex
    r = run([sys.executable, str(CONVERTER), str(md)])
    if r.returncode != 0:
        result["error"] = f"convert failed: {r.stderr[:500]}"
        return result
    result["tex_bytes"] = tex.stat().st_size

    # 2. pdflatex x 2
    for i in (1, 2):
        r = run([PDFLATEX, "-interaction=nonstopmode", "-halt-on-error", f"{stem}.tex"], cwd=WHITEPAPERS)
        if not pdf_src.exists():
            # extract last error
            log = WHITEPAPERS / f"{stem}.log"
            err = ""
            if log.exists():
                txt = log.read_text(encoding="utf-8", errors="replace")
                # last !-prefixed error
                for line in reversed(txt.splitlines()):
                    if line.startswith("!"):
                        err = line
                        break
            result["error"] = f"pdflatex pass {i} produced no PDF: {err}"
            return result

    # 3. Backup old pandoc PDF
    BACKUP.mkdir(exist_ok=True)
    if pdf_dst.exists():
        shutil.copy2(pdf_dst, BACKUP / pdf_dst.name)

    # 4. Replace
    shutil.move(str(pdf_src), str(pdf_dst))
    result["pdf_bytes"] = pdf_dst.stat().st_size

    # 5. count pages by inspecting %%EOF / /Type /Page in the PDF (cheap heuristic)
    try:
        with open(pdf_dst, "rb") as f:
            data = f.read()
        result["pages"] = data.count(b"/Type /Page ") + data.count(b"/Type/Page ")
        if result["pages"] == 0:
            # fallback: count /Type /Page lines without trailing space
            result["pages"] = data.count(b"/Type /Page")
    except Exception:
        pass

    # 6. clean aux files
    for ext in (".aux", ".log", ".out", ".toc"):
        p = WHITEPAPERS / f"{stem}{ext}"
        if p.exists():
            p.unlink()

    result["ok"] = True
    return result


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    results = []
    for stem in argv[1:]:
        print(f">>> {stem}")
        r = build(stem)
        results.append(r)
        if r["ok"]:
            print(f"    OK tex={r['tex_bytes']}B pdf={r['pdf_bytes']}B pages~={r['pages']}")
        else:
            print(f"    FAIL {r['error']}")
    print()
    print("=" * 70)
    okc = sum(1 for r in results if r["ok"])
    print(f"Summary: {okc}/{len(results)} OK")
    for r in results:
        if not r["ok"]:
            print(f"  FAIL {r['stem']}: {r['error']}")
    return 0 if okc == len(results) else 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
