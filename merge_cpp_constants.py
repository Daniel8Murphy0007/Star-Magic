#!/usr/bin/env python3
"""
Merge missing constants from MAIN_1_CoAnQi.cpp into QCalc.py CONSTANTS dict.
"""

import re

def main():
    # Load missing constants
    from missing_constants import MISSING_CONSTANTS
    
    # Filter to ASCII-only valid Python identifiers
    valid_constants = {}
    for name, value in MISSING_CONSTANTS.items():
        if name.isascii() and name.isidentifier():
            valid_constants[name] = value
    
    print(f"Loading {len(valid_constants)} valid missing constants to merge...")
    print(f"(Skipped {len(MISSING_CONSTANTS) - len(valid_constants)} invalid names)")
    
    # Read QCalc.py
    with open('QCalc.py', 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Find the line with 'yaw': -90.0 and insert after it
    insert_pattern = r"('yaw':\s*-90\.0,\n)"
    
    if "# CPP_EXTRACTED_CONSTANTS" in content:
        print("Constants already merged (CPP_EXTRACTED_CONSTANTS marker found)")
        return
    
    # Build the new constants block
    new_block = "\n    # =========================================================================\n"
    new_block += "    # CPP_EXTRACTED_CONSTANTS (from MAIN_1_CoAnQi.cpp - 573 Wolfram entities)\n"
    new_block += "    # =========================================================================\n"
    
    for name in sorted(valid_constants.keys()):
        value = valid_constants[name]
        new_block += f'    "{name}": {value},\n'
    
    # Replace pattern with pattern + new block
    modified = re.sub(insert_pattern, r'\1' + new_block, content)
    
    if modified == content:
        print("ERROR: Could not find insertion point")
        return
    
    # Write back
    with open('QCalc.py', 'w', encoding='utf-8') as f:
        f.write(modified)
    
    print(f"Successfully merged {len(valid_constants)} constants into QCalc.py")
    
    # Verify
    import importlib
    import QCalc
    importlib.reload(QCalc)
    print(f"\nVerification: QCalc.CONSTANTS now has {len(QCalc.CONSTANTS)} constants")

if __name__ == "__main__":
    main()
