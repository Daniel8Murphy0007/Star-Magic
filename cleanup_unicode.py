#!/usr/bin/env python3
with open('MAIN_1_repaired.cpp', 'r', encoding='utf-8', errors='replace') as f:
    content = f.read()
    
# Aggressive Unicode replacement
content = content.replace('\u2018', "'")  # Left single quotation mark
content = content.replace('\u2019', "'")  # Right single quotation mark (apostrophe)
content = content.replace('\u201c', '"')  # Left double quotation mark
content = content.replace('\u201d', '"')  # Right double quotation mark
content = content.replace('\u00b3', '^3')  # Superscript three
content = content.replace('\u00b2', '^2')  # Superscript two
content = content.replace('\u00b9', '^1')  # Superscript one
content = content.replace('\u00b0', ' degrees')  # Degree sign

with open('MAIN_1_repaired.cpp', 'w', encoding='utf-8') as f:
    f.write(content)
    
print("✓ Additional Unicode cleanup complete")
