#!/usr/bin/env python3
"""
generate_pdfs.py
Generate PDFs from all whitepapers/PAPER_*.md using pandoc + pdflatex (arXiv approved).

NOTE: pdflatex is the ONLY approved engine (matches arXiv submission standard).
      xelatex is NOT used here. All .md sources use arXiv LaTeX math syntax ($...$).

Preprocessing handles:
  - Non-UTF-8 encoded files (cp1252 fallback)
  - U+FFFD replacement chars from earlier processing
  - YAML metadata block parse errors (disabled)
  - Problematic superscript/special Unicode in text

Usage:
    python generate_pdfs.py             # all whitepapers
    python generate_pdfs.py 633 642     # specific range (inclusive)
    python generate_pdfs.py --new       # only papers with no existing PDF
"""

import os, sys, glob, subprocess, time, concurrent.futures, re, tempfile

WHITEPAPER_DIR = "whitepapers"
PDF_DIR        = "pdf"
HEADER_FILE    = "pdf_header.tex"
WORKERS        = 2
TIMEOUT_SEC    = 180

BASE_CMD = [
    "pandoc",
    "--pdf-engine=pdflatex",  # arXiv approved engine
    "-V", "geometry:a4paper,top=0.75in,bottom=0.75in,left=0.75in,right=0.75in",
    "-V", "fontsize=11pt",
    "-V", "documentclass=article",
    "--pdf-engine-opt=-interaction=nonstopmode",
    "--pdf-engine-opt=-pool-size=10000000",   # fix: tiletter memory overflow
    "--pdf-engine-opt=-extra-mem-top=10000000",
    "--from=markdown-yaml_metadata_block-raw_tex+smart",
    "--standalone",
    "--wrap=none",
]

_PNUM = re.compile(r"PAPER_(\d+)")

# Characters that cause "invalid character" in LaTeX
# U+FFFD: replacement char from errors='replace' encoding
# U+202F: narrow no-break space
# U+0000-U+001F: control chars (except tab/newline)
_BAD_CHARS = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f\ufffd\u202f]')

def read_as_utf8(path):
    """Read file, trying UTF-8 first then cp1252/latin-1."""
    for enc in ('utf-8', 'utf-8-sig', 'cp1252', 'latin-1'):
        try:
            with open(path, encoding=enc) as f:
                return f.read()
        except (UnicodeDecodeError, LookupError):
            continue
    # Last resort: raw bytes with replace
    with open(path, encoding='utf-8', errors='replace') as f:
        return f.read()

def preprocess(path, aggressive=False):
    """Return path to a clean UTF-8 temp file ready for pandoc."""
    text = read_as_utf8(path)
    # Remove problematic control/replacement characters
    text = _BAD_CHARS.sub('', text)
    if aggressive:
        # Strip ALL non-ASCII outside of math blocks (for stubborn papers)
        # Split on $$ and $ delimiters to preserve math, strip text sections
        # Simple heuristic: strip non-ASCII that's not already a safe Unicode letter
        # For xelatex + DejaVu: keep Greek U+0370-03FF, math U+2200-22FF
        KEEP = set(range(0x0370, 0x0400)) | set(range(0x2200, 0x2300))
        cleaned = []
        for ch in text:
            o = ord(ch)
            if o < 128 or o in KEEP:
                cleaned.append(ch)
            elif o == 0x00B7:   # middle dot
                cleaned.append('*')
            elif o in (0x2013, 0x2014):   # en/em dash
                cleaned.append('--')
            elif o in (0x2018, 0x2019, 0x201C, 0x201D):  # smart quotes
                cleaned.append('"')
            elif o == 0x2248:   # approximately
                cleaned.append('~')
            elif o == 0x00D7:   # multiplication sign
                cleaned.append('x')
            elif o >= 0x00B0 and o <= 0x00FF:  # latin supplement
                cleaned.append(ch)  # keep - DejaVu Serif has these
            # else: drop (removes U+FFFD, replacement chars, etc.)
        text = ''.join(cleaned)
    # Write temp UTF-8 file
    tmp = tempfile.NamedTemporaryFile(
        mode='w', encoding='utf-8', suffix='.md', delete=False
    )
    tmp.write(text)
    tmp.close()
    return tmp.name

def paper_num(fname):
    m = _PNUM.search(fname)
    return int(m.group(1)) if m else 0

def generate_pdf(md_path):
    fname    = os.path.basename(md_path)
    pdf_name = fname.replace(".md", ".pdf")
    pdf_path = os.path.join(PDF_DIR, pdf_name)

    tmp_path = None
    try:
        tmp_path = preprocess(md_path)
        cmd = BASE_CMD + [tmp_path, "-o", pdf_path]
        r = subprocess.run(cmd, capture_output=True, timeout=TIMEOUT_SEC)
        if r.returncode == 0 and os.path.exists(pdf_path):
            return (paper_num(fname), fname, True, f"{os.path.getsize(pdf_path)//1024}KB (pdflatex)", None)
        # Try aggressive preprocessing (strip non-ASCII causing pdflatex issues)
        os.unlink(tmp_path)
        tmp_path = preprocess(md_path, aggressive=True)
        cmd_agg = BASE_CMD + [tmp_path, "-o", pdf_path]
        r_agg = subprocess.run(cmd_agg, capture_output=True, timeout=TIMEOUT_SEC)
        if r_agg.returncode == 0 and os.path.exists(pdf_path):
            return (paper_num(fname), fname, True, f"{os.path.getsize(pdf_path)//1024}KB (pdflatex-agg)", None)
        err = r.stderr.decode("utf-8", errors="replace")[-300:]
        return (paper_num(fname), fname, False, None, err)
    except subprocess.TimeoutExpired:
        return (paper_num(fname), fname, False, None, "TIMEOUT")
    except Exception as e:
        return (paper_num(fname), fname, False, None, str(e))
    finally:
        if tmp_path and os.path.exists(tmp_path):
            try: os.unlink(tmp_path)
            except: pass

def main():
    os.makedirs(PDF_DIR, exist_ok=True)
    all_files = sorted(
        glob.glob(os.path.join(WHITEPAPER_DIR, "PAPER_*.md")),
        key=lambda f: paper_num(os.path.basename(f))
    )
    only_new   = "--new" in sys.argv
    range_args = [a for a in sys.argv[1:] if a.isdigit()]
    if len(range_args) == 2:
        lo, hi = int(range_args[0]), int(range_args[1])
        all_files = [f for f in all_files if lo <= paper_num(os.path.basename(f)) <= hi]
        print(f"Range: PAPER_{lo}-{hi} -> {len(all_files)} files")
    elif only_new:
        all_files = [f for f in all_files
                     if not os.path.exists(os.path.join(PDF_DIR, os.path.basename(f).replace(".md",".pdf")))]
        print(f"New-only: {len(all_files)} files")

    total = len(all_files)
    if not total:
        print("No files to process.")
        return

    print(f"\nGenerating {total} PDFs -> {PDF_DIR}/  (workers={WORKERS})\n")
    t0 = time.time()
    done, failures = 0, []

    with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(generate_pdf, f): f for f in all_files}
        for fut in concurrent.futures.as_completed(futs):
            pnum, fname, ok, info, err = fut.result()
            done += 1
            elapsed = time.time() - t0
            rate = done / elapsed if elapsed > 0 else 1
            eta  = int((total - done) / rate)
            if ok:
                print(f"  [{done:3d}/{total}] PAPER_{pnum:04d}  {info}  ETA {eta}s")
            else:
                failures.append((pnum, fname, err))
                print(f"  [{done:3d}/{total}] PAPER_{pnum:04d}  FAIL  {fname[:55]}")
                if err: print(f"            {err[:100]}")

    elapsed = int(time.time() - t0)
    pdf_count = len(glob.glob(os.path.join(PDF_DIR, "PAPER_*.pdf")))
    print(f"\n{'='*55}")
    print(f"  generate_pdfs.py - Session 164 complete")
    print(f"{'='*55}")
    print(f"  Processed : {total}  |  OK : {total-len(failures)}  |  Fail : {len(failures)}")
    print(f"  PDFs in {PDF_DIR}/  : {pdf_count}")
    print(f"  Time       : {elapsed//60}m {elapsed%60}s")
    print(f"{'='*55}")
    if failures:
        print("\nFAILURES:")
        for pnum, fname, err in sorted(failures):
            print(f"  PAPER_{pnum:04d}  {err[:150] if err else ''}")

if __name__ == "__main__":
    main()
