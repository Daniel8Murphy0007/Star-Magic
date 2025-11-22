#!/usr/bin/env python3
"""
Fix integration errors:
1. Comment out Qt dependencies (QDialog, QWidget)
2. Remove duplicate class definitions
3. Keep all physics intact
"""

import re
from pathlib import Path

def fix_integration_file():
    filepath = Path("MAIN_1_CoAnQi.cpp")
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # Track class definitions we've seen
    seen_classes = {}
    lines = content.split('\n')
    fixed_lines = []
    skip_until_close = False
    brace_depth = 0
    
    for i, line in enumerate(lines):
        # Check for Qt includes and dependencies
        if any(qt in line for qt in ['#include <Q', 'QDialog', 'QWidget', 'QPushButton', 'QLabel']):
            fixed_lines.append('// [Qt dependency commented out] ' + line)
            continue
        
        # Check for class definition start
        class_match = re.match(r'^(class\s+(\w+))\s*(?::\s*public\s+\w+)?\s*{?\s*$', line.strip())
        if class_match:
            class_name = class_match.group(2)
            
            # If we've seen this exact class before, comment it out
            if class_name in seen_classes:
                fixed_lines.append(f'// [Duplicate class definition] {line}')
                skip_until_close = True
                brace_depth = 0
                continue
            else:
                seen_classes[class_name] = i
                fixed_lines.append(line)
                if '{' in line:
                    brace_depth = 1
                continue
        
        # Track braces for duplicate class skipping
        if skip_until_close:
            if '{' in line:
                brace_depth += line.count('{')
            if '}' in line:
                brace_depth -= line.count('}')
                
            fixed_lines.append(f'// [Duplicate] {line}')
            
            if brace_depth == 0 and '}' in line:
                skip_until_close = False
            continue
        
        # Keep all other lines
        fixed_lines.append(line)
    
    # Write fixed content
    fixed_content = '\n'.join(fixed_lines)
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(fixed_content)
    
    print(f"Fixed {len(lines)} lines")
    print(f"Commented out {len([l for l in fixed_lines if l.startswith('//')])} duplicate/Qt lines")
    print(f"Preserved {len([l for l in fixed_lines if not l.startswith('//')])} physics lines")

if __name__ == "__main__":
    fix_integration_file()
