#!/usr/bin/env python3
"""Check header format of sample whitepapers."""
import glob, os

files = sorted(glob.glob('whitepapers/PAPER_*.md'))
# Sample papers at regular intervals
samples = [files[i] for i in [50, 100, 150, 200, 250, 300, 400]]
for wp in samples:
    with open(wp, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    print(f'=== {os.path.basename(wp)} ===')
    for i, line in enumerate(lines[:6], start=1):
        print(f'{i}: {repr(line.rstrip())[:110]}')
    print()
