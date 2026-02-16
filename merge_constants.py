#!/usr/bin/env python3
"""
Merge extracted_constants.py into QCalc.py CONSTANTS dict
- Identifies NEW constants not already in QCalc.py
- Generates merge code block
- Updates QCalc.py directly
"""
import re
import json

# Read extracted constants
with open('extracted_constants.py', 'r', encoding='utf-8') as f:
    content = f.read()
    # Extract the dict
    match = re.search(r'EXTRACTED_CONSTANTS\s*=\s*\{(.+)\}', content, re.DOTALL)
    if match:
        dict_content = '{' + match.group(1) + '}'
        # Safely evaluate
        extracted = eval(dict_content)
    else:
        raise ValueError("Could not parse EXTRACTED_CONSTANTS")

print(f"Loaded {len(extracted)} extracted constants")

# Read QCalc.py CONSTANTS
with open('QCalc.py', 'r', encoding='utf-8') as f:
    qcalc_content = f.read()

# Extract existing CONSTANTS dict keys
existing_keys = set()
# Find CONSTANTS = { ... } block
const_match = re.search(r"CONSTANTS\s*=\s*\{", qcalc_content)
if const_match:
    start = const_match.end()
    # Find matching closing brace
    brace_count = 1
    i = start
    while brace_count > 0 and i < len(qcalc_content):
        if qcalc_content[i] == '{':
            brace_count += 1
        elif qcalc_content[i] == '}':
            brace_count -= 1
        i += 1
    const_block = qcalc_content[const_match.start():i]
    
    # Extract keys using regex
    key_pattern = re.compile(r"'([a-zA-Z_][a-zA-Z0-9_]*)':")
    existing_keys = set(key_pattern.findall(const_block))
    print(f"Found {len(existing_keys)} existing keys in QCalc.py CONSTANTS")

# Find NEW constants (not in QCalc.py)
new_constants = {}
duplicate_count = 0
for key, value in extracted.items():
    if key not in existing_keys:
        new_constants[key] = value
    else:
        duplicate_count += 1

print(f"New constants to add: {len(new_constants)}")
print(f"Skipped duplicates: {duplicate_count}")

# Generate merge block
merge_block = """
    # ═══════════════════════════════════════════════════════════════════════════
    # EXTRACTED FROM source*.js (163 files) - Auto-merged
    # ═══════════════════════════════════════════════════════════════════════════
"""

# Group by category (based on naming patterns)
categories = {
    'Magnetic': [], 'Mass': [], 'Distance': [], 'Time': [],
    'Cosmological': [], 'Quantum': [], 'Coupling': [], 'Velocity': [],
    'Energy': [], 'Frequency': [], 'Temperature': [], 'Other': []
}

for key, value in sorted(new_constants.items()):
    added = False
    if key.startswith(('B_', 'B0', 'B_')) or 'magnetic' in key.lower():
        categories['Magnetic'].append((key, value))
        added = True
    elif key.startswith(('M_', 'mass', 'm_')) or 'mass' in key.lower():
        categories['Mass'].append((key, value))
        added = True
    elif key.startswith(('d_', 'r_', 'dg', 'distance', 'radius')):
        categories['Distance'].append((key, value))
        added = True
    elif key.startswith(('t_', 'tau_', 'T_', 'dt', 'period', 'age')):
        categories['Time'].append((key, value))
        added = True
    elif key.startswith(('H0', 'Hz', 'Lambda', 'Omega', 'z_')):
        categories['Cosmological'].append((key, value))
        added = True
    elif key.startswith(('Delta', 'delta', 'epsilon', 'integral')):
        categories['Quantum'].append((key, value))
        added = True
    elif key.startswith(('k_', 'k1', 'k2', 'k3', 'k4', 'kappa')):
        categories['Coupling'].append((key, value))
        added = True
    elif key.startswith(('v_', 'V_', 'velocity')):
        categories['Velocity'].append((key, value))
        added = True
    elif key.startswith(('E_', 'L_', 'P_', 'energy')):
        categories['Energy'].append((key, value))
        added = True
    elif key.startswith(('f_', 'f', 'omega', 'freq')):
        categories['Frequency'].append((key, value))
        added = True
    elif key.startswith(('T_', 'temp')):
        categories['Temperature'].append((key, value))
        added = True
    if not added:
        categories['Other'].append((key, value))

# Build merge block
for cat_name, items in categories.items():
    if items:
        merge_block += f"\n    # {cat_name}\n"
        for key, value in items:
            # Format value
            if isinstance(value, float):
                if abs(value) > 1e6 or (abs(value) < 1e-3 and value != 0):
                    val_str = f"{value:.6e}"
                else:
                    val_str = str(value)
            else:
                val_str = str(value)
            merge_block += f"    '{key}': {val_str},\n"

# Find insertion point (before closing brace of CONSTANTS)
# Look for the last entry before }
insert_pattern = re.compile(r"(CONSTANTS\s*=\s*\{.+?)(^\})", re.MULTILINE | re.DOTALL)
match = insert_pattern.search(qcalc_content)

if match:
    # Find the position of the closing brace
    const_start = qcalc_content.find("CONSTANTS = {")
    brace_count = 0
    end_pos = const_start
    
    for i, c in enumerate(qcalc_content[const_start:], const_start):
        if c == '{':
            brace_count += 1
        elif c == '}':
            brace_count -= 1
            if brace_count == 0:
                end_pos = i
                break
    
    # Insert before the closing brace
    new_content = qcalc_content[:end_pos] + merge_block + qcalc_content[end_pos:]
    
    # Write back
    with open('QCalc.py', 'w', encoding='utf-8') as f:
        f.write(new_content)
    
    print(f"\n✓ Merged {len(new_constants)} new constants into QCalc.py")
    print(f"  Insertion point: line ~{qcalc_content[:end_pos].count(chr(10))}")
else:
    print("ERROR: Could not find CONSTANTS closing brace")
    
# Summary
print("\nSummary by category:")
for cat, items in sorted(categories.items(), key=lambda x: -len(x[1])):
    if items:
        print(f"  {cat}: {len(items)}")
