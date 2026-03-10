#!/usr/bin/env python3
"""
Fix duplicate classes in CondensedPhysics2.py:
  - Remove 4 true duplicates (same physics, keep first/better instance)
  - Rename 9 different-physics classes (rename second instance)
"""
import re

with open('CondensedPhysics2.py', 'r', encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

print(f"Total lines before: {len(lines)}")

# Find all top-level class/def positions (col 0)
top_level_starts = [
    i for i, line in enumerate(lines)
    if re.match(r'^(?:class|def)\s+\w+', line)
]

def find_class_occurrences(name):
    """Returns list of (start_idx, end_idx) for each top-level occurrence."""
    occurrences = []
    for i, line in enumerate(lines):
        if re.match(rf'^class\s+{re.escape(name)}\b', line):
            # end = start of next top-level definition
            end = len(lines)
            for tl in top_level_starts:
                if tl > i:
                    end = tl
                    break
            occurrences.append((i, end))
    return occurrences

# ── 4 TRUE DUPLICATES: remove second occurrence ──────────────────────────────
to_remove = [
    'EnergyDensityDecayCalculator',         # keep L≈12818, remove L≈14175
    'NegativeTimeUPCalculator',             # keep L≈10444, remove L≈12155
    'UnitConversionCalculator',             # keep L≈31869, remove L≈33172
    'UniversalPermanenceCalculator',        # keep L≈7366,  remove L≈12450
]

# ── 9 DIFFERENT PHYSICS: rename second occurrence ────────────────────────────
to_rename = [
    ('CycleDynamicsCalculator',                  'CycleDynamicsOrb33Calculator'),
    ('ErrorPropagationCalculator',               'ErrorPropagationFormulasCalculator'),
    ('InterferenceFactorCalculator',             'InterferenceFactorComplexCalculator'),
    ('OutflowPressureCalculator',                'OutflowRamPressureCalculator'),
    ('SelfSimilarQuotientCalculator',            'SSqExponentScalingCalculator'),
    ('StressEnergyTensorCalculator',             'VacuumStressEnergyTensorCalculator'),
    ('TimeReversalZoneCalculator',               'TimeReversalZoneDynamicCalculator'),
    ('UniversalPermanenceEquationCalculator',    'UniversalPermanenceMultiplicativeCalculator'),
    ('YangMillsMassGapCalculator',               'YangMillsStringSpectrumCalculator'),
]

# ── Validate before applying ──────────────────────────────────────────────────
print("\n--- Verification pass ---")
lines_to_remove = set()
line_renames = {}  # line_idx -> (old_name, new_name)

ok = True
for name in to_remove:
    occs = find_class_occurrences(name)
    if len(occs) < 2:
        print(f"  ERROR: {name} has {len(occs)} occurrence(s), expected 2")
        ok = False
    else:
        start, end = occs[1]
        print(f"  REMOVE  {name}  (lines {start+1}–{end})")
        for i in range(start, end):
            lines_to_remove.add(i)

for old_name, new_name in to_rename:
    occs = find_class_occurrences(old_name)
    if len(occs) < 2:
        print(f"  ERROR: {old_name} has {len(occs)} occurrence(s), expected 2")
        ok = False
    else:
        start, _ = occs[1]
        print(f"  RENAME  {old_name}  (line {start+1}) → {new_name}")
        line_renames[start] = (old_name, new_name)

if not ok:
    print("\nAborted — fix errors above before proceeding.")
    exit(1)

# ── Apply changes ─────────────────────────────────────────────────────────────
new_lines = []
removed = 0
renamed = 0

for i, line in enumerate(lines):
    if i in lines_to_remove:
        removed += 1
        continue
    if i in line_renames:
        old_name, new_name = line_renames[i]
        line = re.sub(rf'\b{re.escape(old_name)}\b', new_name, line)
        renamed += 1
    new_lines.append(line)

with open('CondensedPhysics2.py', 'w', encoding='utf-8') as f:
    f.writelines(new_lines)

print(f"\n--- Results ---")
print(f"Lines removed: {removed}")
print(f"Classes renamed: {renamed}")
print(f"Total lines after: {len(new_lines)}")

# ── Final count verification ──────────────────────────────────────────────────
content = ''.join(new_lines)
all_classes = re.findall(r'^class\s+(\w+)', content, re.MULTILINE)
from collections import Counter
counts = Counter(all_classes)
remaining_dupes = {k: v for k, v in counts.items() if v > 1}
print(f"Total class definitions: {len(all_classes)}")
print(f"Unique class names: {len(counts)}")
print(f"Remaining duplicates: {len(remaining_dupes)}")
if remaining_dupes:
    for k, v in sorted(remaining_dupes.items()):
        print(f"  {k} ({v}x)")
else:
    print("[OK] No duplicate class names in CP2")
