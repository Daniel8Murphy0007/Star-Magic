#!/usr/bin/env python3
"""
Pre-commit duplicate checker for CondensedPhysics.py (CP1), CondensedPhysics2.py (CP2),
CondensedPhysics3.py (CP3), and CondensedPhysics4.py (CP4).

WIRED PARALLEL: CP3/CP4 hooks now checked exactly like CP1/CP2 (DYNAMIC _SIMULTANEOUS_CALLING safety).
Called by .git/hooks/pre-commit. Prevents import collisions across the full CP1-4 pipeline.

Exit codes:
  0 = No duplicates found (commit allowed)
  1 = Duplicates found (commit blocked)
  2 = File error (commit allowed with warning)
"""

import re
import sys
import os

def check_duplicates():
    """Check for duplicate Calculator class names across CP1/CP2/CP3/CP4 (parallel wiring)."""
    
    # Get script directory (repo root)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cp_paths = {
        'CP1': os.path.join(script_dir, 'CondensedPhysics.py'),
        'CP2': os.path.join(script_dir, 'CondensedPhysics2.py'),
        'CP3': os.path.join(script_dir, 'CondensedPhysics3.py'),
        'CP4': os.path.join(script_dir, 'CondensedPhysics4.py'),
    }
    
    cp_classes = {}
    for label, path in cp_paths.items():
        if not os.path.exists(path):
            print(f"Warning: {path} not found")
            # Continue; missing CP3/CP4 is non-fatal for older check but we still run what exists
            cp_classes[label] = set()
            continue
        try:
            with open(path, 'r', encoding='utf-8', errors='replace') as f:
                content = f.read()
            cp_classes[label] = set(re.findall(r'^class\s+(\w+)', content, re.MULTILINE))
        except Exception as e:
            print(f"Warning: Error reading {path}: {e}")
            cp_classes[label] = set()
    
    # Report header (parallel to old CP1/CP2)
    print("=" * 70)
    print("CP1/CP2/CP3/CP4 DUPLICATE CHECK (parallel hooks for DYNAMIC_SIMULTANEOUS_CALLING)")
    print("=" * 70)
    for label in ['CP1', 'CP2', 'CP3', 'CP4']:
        print(f"{label} classes (unique): {len(cp_classes.get(label, set()))}")
    
    # All pairwise checks (parallel wiring)
    pairs = [('CP1','CP2'), ('CP1','CP3'), ('CP1','CP4'), ('CP2','CP3'), ('CP2','CP4'), ('CP3','CP4')]
    any_duplicates = False
    all_dup_details = {}
    
    for a, b in pairs:
        dups = cp_classes.get(a, set()) & cp_classes.get(b, set())
        if dups:
            any_duplicates = True
            all_dup_details[f"{a}/{b}"] = dups
    
    print()
    if any_duplicates:
        print("ERROR: Duplicate Calculator classes found across CP layers!")
        print()
        for pair, dups in all_dup_details.items():
            print(f"[{pair}] IDENTICAL NAME classes (import conflict risk):")
            for d in sorted(dups):
                print(f"  - {d}")
            print()
        print("RESOLUTION OPTIONS (apply to the higher CP layer):")
        print("  1. REMOVE if identical physics (already in lower CP)")
        print("  2. RENAME with layer suffix (e.g. ...CP4, ...Orb, ...Simul)")
        print("  3. Move to proper CP layer per CONDENSEDPHYSICS_ARCHITECTURE_REFRESH.md")
        print()
        print("Example: FooBarCalculator in CP3 -> FooBarCP3Calculator or FooBarSimultaneousCalculator")
        print()
        return 1  # Block commit
    else:
        print("[OK] No duplicate class names between any CP1/CP2/CP3/CP4 pair")
        print("     (full parallel hook surface ready for CP3/CP4 dynamic simultaneous use)")
        print()
        return 0  # Allow commit


if __name__ == '__main__':
    sys.exit(check_duplicates())
