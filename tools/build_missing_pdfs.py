#!/usr/bin/env python3
"""
build_missing_pdfs.py — Build PDFs for every whitepaper .md that lacks a corresponding .pdf.

Usage:
    python tools/build_missing_pdfs.py

Runs pandoc + xelatex over each md in whitepapers/ that has no matching pdf/*.pdf.
Writes progress every 10 papers to build_pdfs_progress.log.

Estimated wall time on typical desktop: 3-8 hours (674 papers x avg ~15s each).
Safe to interrupt; already-built PDFs are skipped on re-run.

Requirements: pandoc, xelatex (TeX Live), DejaVu Serif font.
On Windows: install via miktex + pandoc from official sites.
"""
import subprocess, os, sys, time, re

WHITEPAPERS = 'whitepapers'
PDF_OUT = 'pdf'
PROGRESS = 'build_pdfs_progress.log'
DETAIL = 'build_pdfs_detail.log'

def find_missing():
    have_pdf = set()
    for d in [PDF_OUT, WHITEPAPERS]:
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            m = re.match(r'PAPER_(\d+)_', f)
            if m and f.endswith('.pdf'):
                have_pdf.add(int(m.group(1)))
    missing = []
    for f in sorted(os.listdir(WHITEPAPERS)):
        m = re.match(r'PAPER_(\d+)_', f)
        if not (m and f.endswith('.md') and '_ASCII' not in f and '_RESERVED' not in f):
            continue
        n = int(m.group(1))
        if n not in have_pdf:
            missing.append((n, f))
    return missing

def build_one(md_name):
    src = os.path.join(WHITEPAPERS, md_name)
    dst = os.path.join(PDF_OUT, md_name[:-3] + '.pdf')
    if os.path.exists(dst):
        return 'skip', 0
    try:
        r = subprocess.run(['pandoc', src, '-o', dst,
            '--pdf-engine=xelatex',
            '-V', 'mainfont=DejaVu Serif',
            '-V', 'monofont=DejaVu Sans Mono',
            '-V', 'sansfont=DejaVu Sans',
            '-V', 'geometry:margin=1in',
            '-V', 'colorlinks=true',
            '-V', 'linkcolor=blue',
            '-V', 'urlcolor=blue'],
            capture_output=True, text=True, timeout=180)
        if r.returncode == 0 and os.path.exists(dst) and os.path.getsize(dst) > 5000:
            return 'ok', os.path.getsize(dst)
        return 'fail', r.stderr[:500]
    except subprocess.TimeoutExpired:
        return 'timeout', 0
    except Exception as e:
        return 'exc', str(e)

def main():
    missing = find_missing()
    print(f"Papers with missing PDF: {len(missing)}")
    if not missing:
        return
    os.makedirs(PDF_OUT, exist_ok=True)
    ok = fail = skip = 0
    t0 = time.time()
    with open(PROGRESS, 'w') as pf, open(DETAIL, 'w') as df:
        pf.write(f"Started at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        pf.write(f"Papers to build: {len(missing)}\n\n")
        pf.flush()
        for i, (n, md_name) in enumerate(missing, 1):
            status, info = build_one(md_name)
            if status == 'ok':
                ok += 1
            elif status == 'skip':
                skip += 1
            else:
                fail += 1
                df.write(f"{status} PAPER_{n}: {info}\n")
                df.flush()
            if i % 10 == 0 or i == len(missing):
                elapsed = time.time() - t0
                rate = i / elapsed
                eta_min = (len(missing) - i) / rate / 60 if rate > 0 else 0
                msg = f"[{i}/{len(missing)}] ok={ok} fail={fail} skip={skip} elapsed={elapsed/60:.1f}min rate={rate:.2f}/s eta={eta_min:.0f}min"
                print(msg, flush=True)
                pf.write(msg + "\n")
                pf.flush()
        pf.write(f"\nDone at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        pf.write(f"Final: ok={ok} fail={fail} skip={skip}\n")

if __name__ == '__main__':
    main()
