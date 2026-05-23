#!/usr/bin/env python3
"""Count wiring details in sessions 389-787"""
import os
import re
from collections import defaultdict

session_files = defaultdict(list)
session_papers = defaultdict(list)

for root, dirs, files in os.walk('.'):
    for f in files:
        # Match session files 389-787
        match = re.search(r'_session(\d+)', f)
        if match:
            sess_num = int(match.group(1))
            if 389 <= sess_num <= 787:
                session_files[sess_num].append(f)
        
        # Match PAPER files (late range)
        paper_match = re.search(r'PAPER_(\d+)', f)
        if paper_match:
            paper_num = int(paper_match.group(1))
            if paper_num >= 700:  # Late papers
                session_papers[paper_num].append(f)

total_session_files = sum(len(v) for v in session_files.values())
unique_sessions = len(session_files)

print(f"Sessions 389-787:\n")
print(f"  Total session files: {total_session_files}")
print(f"  Unique session numbers: {unique_sessions}")
if session_files:
    print(f"  Range of sessions found: {min(session_files.keys())} to {max(session_files.keys())}")
print()

file_types = defaultdict(int)
for sess, files in session_files.items():
    for f in files:
        if f.endswith('.py'):
            file_types['Python scripts'] += 1
        elif f.endswith('.txt'):
            file_types['Text outputs'] += 1
        elif f.endswith('.json'):
            file_types['JSON'] += 1
        elif f.endswith('.md'):
            file_types['Markdown'] += 1
        elif f.endswith('.csv'):
            file_types['CSV'] += 1
        else:
            ext = f.split('.')[-1]
            file_types[f'Other ({ext})'] += 1

print("Breakdown by file type:")
for ftype, count in sorted(file_types.items(), key=lambda x: -x[1]):
    print(f"  {ftype:30s}: {count:4d}")

print()
print(f"Late papers (PAPER_700+):")
print(f"  Count: {len(session_papers)}")

print()
print("=" * 60)
print("TOTAL ADDITIONAL WIRING DETAILS (Sessions 389-787):")
print("=" * 60)
print(f"  Session-specific files:  {total_session_files}")
print(f"  Late framework papers:   {len(session_papers)}")
print(f"  ────────────────────────────────────")
print(f"  SUBTOTAL (Sessions):     {total_session_files + len(session_papers)}")
print()
print("COMBINED WITH SESSIONS 201-388:")
print(f"  Sessions 201-388 count:  82 (from WIRING_INVENTORY.md)")
print(f"  Sessions 389-787 count:  {total_session_files + len(session_papers)}")
print(f"  ────────────────────────────────────")
print(f"  TOTAL ALL SESSIONS:      {82 + total_session_files + len(session_papers)}")

# Sample the sessions found
if session_files:
    print()
    print("Sample sessions found:")
    for sess in sorted(session_files.keys())[:10]:
        files = session_files[sess]
        print(f"  Session {sess:3d}: {len(files):2d} files - {', '.join(files[:3])}")
    if len(session_files) > 10:
        print(f"  ... ({unique_sessions - 10} more sessions)")
