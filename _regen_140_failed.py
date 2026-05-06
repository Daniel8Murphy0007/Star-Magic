#!/usr/bin/env python3
"""
_regen_140_failed.py
Regenerate PDFs for the 140 papers that failed in the stale-PDF run.
Must be run AFTER _fix_stale_latex_errors.py has been applied.
"""
import os, glob, sys, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from generate_pdfs import generate_pdf, PDF_DIR, WHITEPAPER_DIR, paper_num

WORKERS = 2  # keep low to avoid pdflatex OOM

FAILED = [
    27,28,30,31,32,33,35,36,51,63,90,106,138,143,144,155,160,163,164,167,
    181,184,188,198,202,210,214,215,224,239,240,242,259,261,262,265,267,268,
    278,298,313,336,351,354,372,373,375,380,381,384,386,389,429,435,439,
    461,462,464,473,491,494,498,513,514,526,532,533,535,536,544,545,549,
    554,557,563,570,573,574,575,576,577,578,581,582,583,585,587,590,592,
    598,600,633,645,646,647,650,651,653,688,692,701,716,717,718,721,722,
    731,732,735,738,739,740,747,749,794,798,803,807,808,812,831,832,833,
    840,865,877,880,882,883,888,890,904,949,953,957,980,1023,1025,1040,1101,
]

def find_md(n):
    for pat in [f'PAPER_{n:04d}_*.md', f'PAPER_{n:03d}_*.md',
                f'PAPER_{n:04d}*.md', f'PAPER_{n:03d}*.md']:
        m = glob.glob(os.path.join(WHITEPAPER_DIR, pat))
        if m:
            return m[0]
    return None

def main():
    targets = []
    missing = []
    for n in FAILED:
        p = find_md(n)
        if p:
            targets.append(p)
        else:
            missing.append(n)
    if missing:
        print(f"WARNING: no .md found for: {missing}")

    print(f"Regenerating {len(targets)} PDFs with {WORKERS} workers...\n")
    ok, fail = [], []
    t0 = time.time()

    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(generate_pdf, p): p for p in targets}
        for fut in as_completed(futures):
            md_path = futures[fut]
            num = paper_num(md_path)
            try:
                result = fut.result()
                if result:
                    ok.append(num)
                    print(f"  OK   {num:4d}  {os.path.basename(md_path)}")
                else:
                    fail.append(num)
                    print(f"  FAIL {num:4d}  {os.path.basename(md_path)}")
            except Exception as e:
                fail.append(num)
                print(f"  ERR  {num:4d}  {e}")

    elapsed = time.time() - t0
    print(f"\n{'='*55}")
    print(f"  OK   : {len(ok)}")
    print(f"  FAIL : {len(fail)}")
    print(f"  Time : {elapsed:.0f}s")
    if fail:
        print(f"\nStill failing: {sorted(fail)}")

if __name__ == '__main__':
    main()
