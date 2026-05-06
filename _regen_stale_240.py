#!/usr/bin/env python3
"""
_regen_stale_240.py
Regenerate the 240 whitepapers whose PDFs are stale (identified by _audit_stale_pdfs.py).
Uses generate_pdfs.py infrastructure (pdflatex, arXiv approved).
"""
import os, sys, glob, time, concurrent.futures
sys.path.insert(0, os.path.dirname(__file__))
from generate_pdfs import generate_pdf, PDF_DIR, WHITEPAPER_DIR, paper_num

STALE = [
    16,
    27, 28,
    30, 31, 32, 33, 34, 35, 36,
    46, 51, 63,
    69, 70, 71,
    90, 92, 94, 106, 119, 138,
    143, 144,
    155,
    158, 159, 160,
    163, 164,
    167, 181, 184, 188,
    197, 198,
    202, 210,
    214, 215,
    217, 222,
    224, 225,
    228,
    239, 240,
    242, 248,
    250, 251,
    255, 257, 259,
    261, 262,
    265,
    267, 268,
    273, 278,
    297, 298, 299,
    301, 304, 313, 332,
    335, 336,
    346, 347,
    349, 350, 351,
    354,
    359, 360,
    370,
    372, 373,
    375,
    380, 381,
    384, 386, 389, 429, 435, 439,
    461, 462,
    464, 473, 477, 491, 494, 498,
    513, 514,
    526,
    532, 533,
    535, 536,
    544, 545,
    549, 554, 557, 563, 570,
    573, 574, 575, 576, 577, 578,
    580, 581, 582, 583,
    585, 587,
    590, 591, 592,
    598, 600, 633,
    645, 646, 647,
    650, 651,
    653,
    660, 661, 662, 663, 664, 665,
    667, 668,
    670, 671, 672, 673,
    688, 689, 690, 691, 692, 693, 694, 695, 696, 697,
    698, 699, 700, 701, 702, 703, 704, 705, 706, 707,
    708, 709, 710, 711, 712, 713, 714, 715, 716, 717,
    718, 719, 720, 721, 722, 723, 724, 725, 726, 727,
    728, 729, 730, 731, 732, 733, 734, 735,
    738, 739, 740,
    747, 749, 794, 798,
    800, 801,
    803,
    807, 808,
    812,
    831, 832, 833,
    835, 838, 840, 865, 877,
    879, 880, 881, 882, 883,
    888, 890, 892, 898, 904, 949, 953, 957,
    979, 980, 981,
    985, 999, 1001, 1023, 1025, 1040, 1101,
]

def main():
    stale_set = set(STALE)
    all_md = sorted(
        glob.glob(os.path.join(WHITEPAPER_DIR, "PAPER_*.md")),
        key=lambda f: paper_num(os.path.basename(f))
    )
    targets = [f for f in all_md if paper_num(os.path.basename(f)) in stale_set]

    missing_nums = stale_set - {paper_num(os.path.basename(f)) for f in all_md}
    if missing_nums:
        print(f"WARNING: {len(missing_nums)} paper(s) not found in whitepapers/:")
        for n in sorted(missing_nums):
            print(f"  PAPER_{n:04d}.md")

    total = len(targets)
    print(f"Regenerating {total} stale PDFs  (engine: pdflatex / arXiv approved)")
    print(f"Output -> {PDF_DIR}/")
    print()

    os.makedirs(PDF_DIR, exist_ok=True)
    t0 = time.time()
    done, failures = 0, []

    WORKERS = 4
    with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(generate_pdf, f): f for f in targets}
        for fut in concurrent.futures.as_completed(futs):
            pnum, fname, ok, info, err = fut.result()
            done += 1
            elapsed = time.time() - t0
            rate = done / elapsed if elapsed > 0 else 1
            eta = int((total - done) / rate)
            if ok:
                print(f"  [{done:3d}/{total}] PAPER_{pnum:04d}  OK  {info}  ETA {eta}s")
            else:
                failures.append((pnum, fname, err))
                print(f"  [{done:3d}/{total}] PAPER_{pnum:04d}  FAIL  {(err or '')[:100]}")

    elapsed = int(time.time() - t0)
    print(f"\n{'='*55}")
    print(f"  Stale PDF regeneration complete")
    print(f"  Processed: {total}  OK: {total-len(failures)}  Fail: {len(failures)}")
    print(f"  Time: {elapsed}s")
    if failures:
        print(f"\nFailed papers:")
        for pnum, fname, err in sorted(failures):
            print(f"  PAPER_{pnum:04d}  {(err or '')[:120]}")

if __name__ == "__main__":
    main()
