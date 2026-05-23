#!/usr/bin/env python3
"""Show summary of all 630 constants by session."""

import csv
import re

scripts = []
total = 0

with open('MASTER_LEDGER_BY_SCRIPT.csv', 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['script']:
            count = int(row['count'])
            scripts.append({
                'script': row['script'],
                'count': count,
                'status': row['status_breakdown']
            })
            total += count

print('CONSTANTS BY SESSION - COMPLETE INVENTORY')
print('=' * 110)
print()
print(f'TOTAL: {total} constants across {len(scripts)} session scripts')
print()
print('TOP 20 SESSIONS BY CONSTANT COUNT:')
print('-' * 110)

for i, s in enumerate(sorted(scripts, key=lambda x: x['count'], reverse=True)[:20], 1):
    print(f"{i:2d}. {s['script']:<55} count={s['count']:>2}  {s['status']}")

print()
print('=' * 110)
print('SESSIONS BY ERA:')
print('-' * 110)

# Extract session numbers
def get_session_number(script_name):
    m = re.search(r'_session(\d+)', script_name)
    return int(m.group(1)) if m else 0

by_session = {get_session_number(s['script']): s for s in scripts if get_session_number(s['script']) > 0}

eras = [
    ('FOUNDATION (201-210)', 201, 210),
    ('EARLY (211-260)', 211, 260),
    ('GROWTH (262-268)', 262, 268),
    ('EXPANSION (270-279)', 270, 279),
    ('CALIBRATION (280-304)', 280, 304),
    ('SPECIALIZED (305-785)', 305, 785),
]

for era_name, start, end in eras:
    count = sum(s['count'] for num, s in by_session.items() if start <= num <= end)
    num_sessions = len([num for num in by_session.keys() if start <= num <= end])
    print(f"  {era_name:<30} {count:>4} constants in {num_sessions:>3} sessions")

print()
print('=' * 110)
print('STATUS SUMMARY:')
print('-' * 110)

status_totals = {}
for s in scripts:
    statuses = s['status'].split(';')
    for status_pair in statuses:
        if '=' in status_pair:
            status, count_str = status_pair.split('=')
            status_totals[status] = status_totals.get(status, 0) + int(count_str)

for status in sorted(status_totals.keys()):
    count = status_totals[status]
    pct = 100 * count / total if total > 0 else 0
    print(f"  {status:<30} {count:>4} constants ({pct:>5.1f}%)")

print()
print(f"  {'TOTAL':<30} {total:>4} constants (100.0%)")
print('=' * 110)
