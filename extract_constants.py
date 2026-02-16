#!/usr/bin/env python3
"""
Extract all constants from JS source files and generate CONSTANTS dict for QCalc.py
"""

import re
import glob
import json

def extract_constants():
    constants = {}
    
    # Scan all JS files for constant definitions
    for filepath in sorted(glob.glob('source*.js')):
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
        
        # Pattern 1: this.NAME = VALUE;  // comment
        pattern1 = r'this\.(\w+)\s*=\s*([0-9eE.+\-]+(?:\s*\*\s*[0-9eE.+\-]+)?)\s*;?\s*(?://\s*(.+))?'
        for match in re.finditer(pattern1, content):
            name = match.group(1)
            value = match.group(2).strip()
            comment = match.group(3) or ''
            if name not in constants:
                try:
                    v = eval(value)
                    if isinstance(v, (int, float)) and v != 0:
                        constants[name] = {'value': value, 'comment': comment.strip()[:60], 'source': filepath}
                except:
                    pass
        
        # Pattern 2: const NAME = VALUE;
        pattern2 = r'const\s+(\w+)\s*=\s*([0-9eE.+\-]+(?:\s*\*\s*[0-9eE.+\-]+)?)\s*;'
        for match in re.finditer(pattern2, content):
            name = match.group(1)
            value = match.group(2).strip()
            if name not in constants:
                try:
                    v = eval(value)
                    if isinstance(v, (int, float)) and v != 0:
                        constants[name] = {'value': value, 'comment': '', 'source': filepath}
                except:
                    pass
    
    return constants

def generate_constants_dict(constants):
    """Generate Python CONSTANTS dict code"""
    lines = []
    lines.append("# =============================================================================")
    lines.append("# CONSTANTS - Extracted from source*.js (163 files)")
    lines.append("# =============================================================================")
    lines.append("")
    lines.append("EXTRACTED_CONSTANTS = {")
    
    # Group by category based on naming patterns
    categories = {
        'Fundamental': [],
        'Magnetic': [],
        'Vacuum': [],
        'Time': [],
        'Mass': [],
        'Distance': [],
        'Cosmological': [],
        'Quantum': [],
        'Other': []
    }
    
    for name, data in sorted(constants.items()):
        value = data['value']
        comment = data['comment']
        
        # Categorize
        if name in ('G', 'c', 'c_light', 'hbar', 'q_charge', 'proton_mass', 'PI'):
            cat = 'Fundamental'
        elif 'B' in name.upper() and ('FIELD' in name.upper() or 'CRIT' in name.upper() or name.startswith('B')):
            cat = 'Magnetic'
        elif 'vac' in name.lower() or 'rho' in name.lower():
            cat = 'Vacuum'
        elif 'tau' in name.lower() or 't_' in name.lower() or 'Hubble' in name:
            cat = 'Time'
        elif name.startswith('M') and not name.startswith('Math'):
            cat = 'Mass'
        elif name.startswith('r') or name.startswith('d') and 'delta' not in name.lower():
            cat = 'Distance'
        elif 'Lambda' in name or 'H0' in name or 'Hz' in name or 'Omega' in name:
            cat = 'Cosmological'
        elif 'hbar' in name.lower() or 'delta_x' in name.lower() or 'psi' in name.lower():
            cat = 'Quantum'
        else:
            cat = 'Other'
        
        categories[cat].append((name, value, comment))
    
    for cat, items in categories.items():
        if items:
            lines.append(f"    # --- {cat} ---")
            for name, value, comment in sorted(items):
                if comment:
                    lines.append(f"    '{name}': {value},  # {comment}")
                else:
                    lines.append(f"    '{name}': {value},")
            lines.append("")
    
    lines.append("}")
    return '\n'.join(lines)

def main():
    print("Extracting constants from source*.js files...")
    constants = extract_constants()
    print(f"Found {len(constants)} unique constants")
    
    # Save to JSON
    with open('extracted_constants.json', 'w') as f:
        json.dump(constants, f, indent=2)
    print("Saved to extracted_constants.json")
    
    # Generate Python code
    code = generate_constants_dict(constants)
    with open('extracted_constants.py', 'w', encoding='utf-8') as f:
        f.write(code)
    print(f"Generated extracted_constants.py ({len(code):,} chars)")
    
    # Print summary
    print("\nSummary by category:")
    categories = {}
    for name, data in constants.items():
        if name in ('G', 'c', 'c_light', 'hbar', 'q_charge', 'proton_mass', 'PI'):
            cat = 'Fundamental'
        elif 'B' in name.upper() and len(name) <= 10:
            cat = 'Magnetic'
        elif 'vac' in name.lower() or 'rho' in name.lower():
            cat = 'Vacuum'
        elif 'tau' in name.lower() or 't_' in name.lower():
            cat = 'Time'
        elif name.startswith('M') and not name.startswith('Math'):
            cat = 'Mass'
        else:
            cat = 'Other'
        categories[cat] = categories.get(cat, 0) + 1
    
    for cat, count in sorted(categories.items(), key=lambda x: -x[1]):
        print(f"  {cat}: {count}")

if __name__ == '__main__':
    main()
