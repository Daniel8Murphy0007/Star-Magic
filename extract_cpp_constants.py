#!/usr/bin/env python3
"""
Extract all constants from MAIN_1_CoAnQi.cpp and compare with QCalc.py
to identify missing constants that need to be integrated.
"""

import re
import os

def extract_cpp_constants(filepath):
    """Extract constants from MAIN_1_CoAnQi.cpp using multiple patterns."""
    constants = {}
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # Pattern 1: constants["NAME"] = value;
    pattern1 = r'constants\s*\[\s*"([A-Za-z_][A-Za-z0-9_]*)"\s*\]\s*=\s*([^;]+);'
    for match in re.finditer(pattern1, content):
        name = match.group(1)
        value_str = match.group(2).strip()
        try:
            # Try to evaluate numeric value
            value = eval(value_str.replace('constants["PI"]', '3.141592653589793'))
            constants[name] = value
        except:
            pass
    
    # Pattern 2: double NAME = value; or const double NAME = value;
    pattern2 = r'(?:const\s+)?double\s+([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([\d\.\-eE+]+)\s*;'
    for match in re.finditer(pattern2, content):
        name = match.group(1)
        value_str = match.group(2)
        try:
            value = float(value_str)
            constants[name] = value
        except:
            pass
    
    # Pattern 3: {\"NAME\", value} or {"NAME", value} in initializer lists
    pattern3 = r'\{\s*"([A-Za-z_][A-Za-z0-9_]*)"\s*,\s*([\d\.\-eE+]+)\s*\}'
    for match in re.finditer(pattern3, content):
        name = match.group(1)
        value_str = match.group(2)
        try:
            value = float(value_str)
            constants[name] = value
        except:
            pass
    
    # Pattern 4: #define NAME value
    pattern4 = r'#define\s+([A-Z_][A-Z0-9_]*)\s+([\d\.\-eE+]+)'
    for match in re.finditer(pattern4, content):
        name = match.group(1)
        value_str = match.group(2)
        try:
            value = float(value_str)
            constants[name] = value
        except:
            pass
    
    # Pattern 5: Parameter structs - double field = value;
    pattern5 = r'(\w+)\s*=\s*([\d\.\-eE+]+)\s*(?:,|;|\})'
    for match in re.finditer(pattern5, content):
        name = match.group(1)
        value_str = match.group(2)
        # Filter out common non-constants
        if name in ['line', 'col', 'i', 'j', 'k', 'n', 'idx', 'size', 'count']:
            continue
        try:
            value = float(value_str)
            # Only add if looks like a physics constant (not loop variable)
            if len(name) > 1 and (name[0].isupper() or '_' in name or name in ['c', 'G', 'h', 'k', 'q']):
                constants[name] = value
        except:
            pass
    
    return constants

def load_qcalc_constants():
    """Load existing QCalc.py constants."""
    import QCalc
    return QCalc.CONSTANTS.copy()

def main():
    cpp_file = "MAIN_1_CoAnQi.cpp"
    
    print("="*80)
    print("EXTRACTING CONSTANTS FROM MAIN_1_CoAnQi.cpp")
    print("="*80)
    
    # Extract C++ constants
    cpp_constants = extract_cpp_constants(cpp_file)
    print(f"\nExtracted from C++: {len(cpp_constants)} constants")
    
    # Load QCalc constants
    qcalc_constants = load_qcalc_constants()
    print(f"Existing in QCalc.py: {len(qcalc_constants)} constants")
    
    # Find missing (in C++ but not in QCalc)
    cpp_keys = set(cpp_constants.keys())
    qcalc_keys = set(qcalc_constants.keys())
    
    missing = cpp_keys - qcalc_keys
    print(f"\nMISSING FROM QCalc.py: {len(missing)} constants")
    
    # Save extracted constants
    with open('cpp_extracted_constants.py', 'w', encoding='utf-8') as f:
        f.write('"""Constants extracted from MAIN_1_CoAnQi.cpp"""\n\n')
        f.write('CPP_EXTRACTED_CONSTANTS = {\n')
        for name in sorted(cpp_constants.keys()):
            value = cpp_constants[name]
            f.write(f'    "{name}": {value},\n')
        f.write('}\n')
    
    # Save missing constants
    with open('missing_constants.py', 'w', encoding='utf-8') as f:
        f.write('"""Constants from MAIN_1_CoAnQi.cpp missing in QCalc.py"""\n\n')
        f.write('MISSING_CONSTANTS = {\n')
        for name in sorted(missing):
            value = cpp_constants[name]
            f.write(f'    "{name}": {value},\n')
        f.write('}\n')
    
    print(f"\nSaved to: cpp_extracted_constants.py ({len(cpp_constants)} total)")
    print(f"Saved to: missing_constants.py ({len(missing)} missing)")
    
    if missing:
        print("\n" + "="*80)
        print("MISSING CONSTANTS (first 50):")
        print("="*80)
        for i, name in enumerate(sorted(missing)[:50], 1):
            print(f"  {i:3}. {name:30} = {cpp_constants[name]}")
        if len(missing) > 50:
            print(f"  ... and {len(missing) - 50} more")
    
    return len(missing)

if __name__ == "__main__":
    main()
