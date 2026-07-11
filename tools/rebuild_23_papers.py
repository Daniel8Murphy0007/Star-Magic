#!/usr/bin/env python3
"""
rebuild_23_papers.py — Rebuild the 23 papers where F:\Aetheric was replaced with F:/Aetheric.
Excludes PAPER_1218 and PAPER_1801 which have separate LaTeX issues deferred to v5.58.1.

Usage:
    python tools/rebuild_23_papers.py

Wall time: ~7 minutes on Windows (23 papers x ~15s avg).
Idempotent: deletes stale PDF first so pandoc overwrites cleanly.
"""
import subprocess, os, glob, time

FIXED_IDS = [1472, 1473, 1474, 1475, 1476, 1477, 1478, 1479, 1480,
             1481, 1482, 1483, 1484, 1485, 1486, 1487, 1488, 1489, 1490,
             1491, 1492, 1493, 1498]

ok = fail = 0
t0 = time.time()
for i, paper_id in enumerate(FIXED_IDS, 1):
    matches = glob.glob(f'whitepapers/PAPER_{paper_id}_*.md')
    if not matches:
        print(f"[{i}/{len(FIXED_IDS)}] PAPER_{paper_id}: no md found")
        fail += 1
        continue
    md_name = os.path.basename(matches[0])
    src = matches[0]
    dst = os.path.join('pdf', md_name[:-3] + '.pdf')
    # Delete stale pdf if any
    if os.path.exists(dst):
        os.remove(dst)
    print(f"[{i}/{len(FIXED_IDS)}] PAPER_{paper_id}: building...", flush=True)
    try:
        r = subprocess.run(['pandoc', src, '-o', dst,
            '--pdf-engine=xelatex',
            '-V', 'mainfont=DejaVu Serif',
            '-V', 'monofont=DejaVu Sans Mono',
            '-V', 'sansfont=DejaVu Sans',
            '-V', 'geometry:margin=1in'],
            capture_output=True, timeout=90)
        if r.returncode == 0 and os.path.exists(dst) and os.path.getsize(dst) > 5000:
            ok += 1
            print(f"    OK ({os.path.getsize(dst)} bytes)")
        else:
            fail += 1
            print(f"    FAIL rc={r.returncode}")
    except Exception as e:
        fail += 1
        print(f"    EXC: {e}")

elapsed = time.time() - t0
print(f"\nDone: ok={ok}/{len(FIXED_IDS)} fail={fail} elapsed={elapsed/60:.1f}min")
