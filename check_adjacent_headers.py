#!/usr/bin/env python3
"""Check papers adjacent to 072, 114, 129, 137 for header format."""
import glob, os

def show_header(wp, n=3):
    try:
        with open(wp, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        print(f'=== {os.path.basename(wp)} ===')
        for i, line in enumerate(lines[:n], start=1):
            print(f'{i}: {repr(line.rstrip())[:110]}')
        print()
    except FileNotFoundError:
        pass

all_files = {int(f.split('PAPER_')[1].split('_')[0]): f
             for f in glob.glob('whitepapers/PAPER_*.md')}

for target in [72, 114, 129, 137]:
    print(f'--- Papers around {target} ---')
    for n in [target-2, target-1, target, target+1, target+2]:
        if n in all_files:
            show_header(all_files[n])
