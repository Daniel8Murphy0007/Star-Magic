#!/usr/bin/env python3
"""Generate quick-reference inventory of all 630 constants by domain and session."""

import csv

constants_by_domain = {}
total_count = 0
status_counts = {}

with open('MASTER_LEDGER_BY_SCRIPT.csv', 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if not row.get('script_id'):
            continue
        
        domain = row.get('physics_domain', 'Unknown')
        status = row.get('validation_status', 'Unknown')
        count = int(row.get('constant_count', 0))
        
        if domain not in constants_by_domain:
            constants_by_domain[domain] = []
        constants_by_domain[domain].append({
            'script': row.get('script_id'),
            'count': count,
            'status': status
        })
        
        status_counts[status] = status_counts.get(status, 0) + count
        total_count += count

print("=" * 90)
print("MASTER CONSTANTS INVENTORY - BY DOMAIN & SESSION")
print("=" * 90)
print()

for domain in sorted(constants_by_domain.keys()):
    scripts = constants_by_domain[domain]
    domain_total = sum(s['count'] for s in scripts)
    print(f"{domain.upper()} ({domain_total} total across {len(scripts)} sessions)")
    print("-" * 90)
    
    # Show first few sessions in each domain
    for script in sorted(scripts, key=lambda x: x['script'])[:5]:
        print(f"  {script['script']:<45} count={script['count']:>3}  status={script['status']}")
    
    if len(scripts) > 5:
        print(f"  ... and {len(scripts)-5} more sessions in this domain ...")
    print()

print("=" * 90)
print("STATUS SUMMARY (630 TOTAL)")
print("=" * 90)
for status in sorted(status_counts.keys()):
    count = status_counts[status]
    pct = 100 * count / total_count
    print(f"{status:<30} {count:>4} constants ({pct:>5.1f}%)")

print()
print(f"TOTAL: {total_count} constants across {len(constants_by_domain)} physics domains")
print("=" * 90)
