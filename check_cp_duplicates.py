#!/usr/bin/env python3
"""
Pre-commit duplicate checker for CondensedPhysics.py and CondensedPhysics2.py

This script is called by the git pre-commit hook to prevent committing
duplicate Calculator classes between CP1 and CP2.

Exit codes:
  0 = No duplicates found (commit allowed)
  1 = Duplicates found (commit blocked)
  2 = File error (commit allowed with warning)
"""

import re
import sys
import os

def check_duplicates():
    """Check for duplicate Calculator class names between CP1 and CP2."""
    
    # Get script directory (repo root)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cp1_path = os.path.join(script_dir, 'CondensedPhysics.py')
    cp2_path = os.path.join(script_dir, 'CondensedPhysics2.py')
    
    # Check files exist
    if not os.path.exists(cp1_path):
        print(f"Warning: {cp1_path} not found")
        return 2
    if not os.path.exists(cp2_path):
        print(f"Warning: {cp2_path} not found")
        return 2
    
    try:
        # Extract class names from CP1
        with open(cp1_path, 'r', encoding='utf-8', errors='replace') as f:
            cp1 = f.read()
        cp1_classes = set(re.findall(r'^class\s+(\w+)', cp1, re.MULTILINE))
        
        # Extract class names from CP2
        with open(cp2_path, 'r', encoding='utf-8', errors='replace') as f:
            cp2 = f.read()
        cp2_classes = set(re.findall(r'^class\s+(\w+)', cp2, re.MULTILINE))
        
    except Exception as e:
        print(f"Warning: Error reading files: {e}")
        return 2
    
    # Find duplicates
    duplicates = cp1_classes & cp2_classes
    
    # Report results
    print("=" * 60)
    print("CP1/CP2 DUPLICATE CHECK")
    print("=" * 60)
    print(f"CP1 classes (unique): {len(cp1_classes)}")
    print(f"CP2 classes (unique): {len(cp2_classes)}")
    print(f"Duplicate class names: {len(duplicates)}")
    print()
    
    if duplicates:
        print("ERROR: Duplicate Calculator classes found!")
        print()
        print("IDENTICAL NAME classes (will cause import conflicts):")
        for d in sorted(duplicates):
            print(f"  - {d}")
        print()
        print("RESOLUTION OPTIONS:")
        print("  1. REMOVE from CP2 if identical physics (already in CP1)")
        print("  2. RENAME in CP2 if different physics (add Orb# suffix)")
        print()
        print("Example rename: HydrogenResonanceCalculator -> HydrogenResonanceOrb36Calculator")
        print()
        return 1  # Block commit
    else:
        print("[OK] No duplicate class names between CP1 and CP2")
        print()
        return 0  # Allow commit


if __name__ == '__main__':
    sys.exit(check_duplicates())
