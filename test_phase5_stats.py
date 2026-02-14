#!/usr/bin/env python3
"""Quick stats check for Phase 5 database update"""

from SymbolicDB import SymbolicDatabase

db = SymbolicDatabase('uqff_equations.db')
stats = db.get_statistics()

print("\n" + "="*80)
print("SYMBOLIC DATABASE STATISTICS (POST PHASE 5)")
print("="*80)
print(f"Total equations: {stats['total_equations']}")
print(f"Self-expanding: {stats['self_expanding_count']}")
print(f"Categories: {len(stats['by_category'])}")
print(f"Sources: {len(stats['by_source'])}")

print("\nBy Category:")
for cat, count in sorted(stats['by_category'].items()):
    print(f"  • {cat}: {count}")

print("\nBy Source:")
for src, count in sorted(stats['by_source'].items()):
    if src and src.startswith('source') and int(src[6:-4]) >= 52:
        print(f"  ★ {src}: {count} (Phase 5)")
    elif src:
        print(f"  • {src}: {count}")

print("\n" + "="*80)
print("Phase 5 Equations:")
print("="*80)
phase5_sources = [f'source{i}.cpp' for i in range(52, 66)]
phase5_eqs = []
for src in phase5_sources:
    eqs = db.query_by_source(src)
    phase5_eqs.extend(eqs)

for i, eq in enumerate(sorted(phase5_eqs), 1):
    print(f"{i:2d}. {eq}")

print(f"\n✓ Total Phase 5 equations: {len(phase5_eqs)}")
print("="*80)
