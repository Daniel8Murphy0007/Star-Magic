#!/usr/bin/env python3
"""Rebuild PAPER_1493 (Aetheric Propulsion second backslash fix)."""
import subprocess, os, glob
matches = glob.glob('whitepapers/PAPER_1493_*.md')
if not matches:
    raise SystemExit("PAPER_1493 md not found")
src = matches[0]
dst = os.path.join('pdf', os.path.basename(src)[:-3] + '.pdf')
if os.path.exists(dst):
    os.remove(dst)
print(f"Building {src}...")
r = subprocess.run(['pandoc', src, '-o', dst,
    '--pdf-engine=xelatex',
    '-V', 'mainfont=DejaVu Serif',
    '-V', 'monofont=DejaVu Sans Mono',
    '-V', 'sansfont=DejaVu Sans',
    '-V', 'geometry:margin=1in'],
    capture_output=True, timeout=90)
if r.returncode == 0 and os.path.exists(dst) and os.path.getsize(dst) > 5000:
    print(f"OK: {os.path.getsize(dst)} bytes")
else:
    print(f"FAIL rc={r.returncode}")
    try:
        print(r.stderr.decode('utf-8', errors='replace')[:500])
    except Exception:
        pass
