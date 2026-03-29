#!/usr/bin/env python3
"""
Fix misplaced default constructors: remove ones placed outside public section,
re-insert them right before the class's existing parameterized constructor (which is in public:).
"""
import re

BRIDGE_FILE = 'wolfram_sources_bridge.cpp'

with open(BRIDGE_FILE, 'r', encoding='utf-8') as f:
    content = f.read()

original_len = len(content)

# These classes got their default ctor placed wrongly (before the private: section)
CLASSES = [
    'CompressedExpansionTerm',
    'CompressedSuperAdjTerm',
    'CompressedFluidTerm',
    'CompressedPerturbationTerm',
    'ResonanceADPMTerm',
    'ResonanceATHzTerm',
    'ResonanceAvacDiffTerm',
    'ResonanceASuperFreqTerm',
    'ResonanceAAetherResTerm',
    'ResonanceUg4iTerm',
    'ResonanceAQuantumFreqTerm',
    'ResonanceAAetherFreqTerm',
    'ResonanceAFluidFreqTerm',
    'ResonanceAExpFreqTerm',
    'ResonanceFTRZTerm',
    'ResonanceWormholeTerm',
]

fixed = 0
for class_name in CLASSES:
    # Remove the misplaced default ctor (before private: section — after class {)
    bad_pattern = re.compile(
        r'(class\s+' + re.escape(class_name) + r'\s*:[^{]+\{)\s*\n\s*' 
        + re.escape(class_name) + r'\(\)\s*=\s*default;\s*\n',
        re.MULTILINE
    )
    content_new = bad_pattern.sub(r'\1\n', content)
    if content_new != content:
        content = content_new
    # Now find the parameterized constructor in public section and insert before it
    # Pattern: the constructor call with params inside public: section
    pub_ctor_pattern = re.compile(
        r'(public:\s*\n(?:.*\n)*?)'              # public: followed by any lines
        r'(\s+' + re.escape(class_name) + r'\([^)]+\))',  # the ctor with params
        re.MULTILINE
    )
    # Simpler approach: find the ctor line directly and insert default ctor before it
    ctor_line_pat = re.compile(
        r'(\s{4})(' + re.escape(class_name) + r'\([^)]+\))',
        re.MULTILINE
    )
    m = ctor_line_pat.search(content)
    if m:
        # Check it's the class definition ctor (not a call)
        pos = m.start()
        preceding = content[max(0, pos-600):pos]
        if f'class {class_name}' in preceding:
            # Check if default ctor already right before it
            before = content[max(0, pos-80):pos]
            if f'{class_name}() = default' not in before:
                indent = m.group(1)
                insert_text = f'{indent}{class_name}() = default;\n'
                content = content[:pos] + insert_text + content[pos:]
                fixed += 1
                print(f'[OK] Fixed {class_name}')
            else:
                print(f'[SKIP] {class_name} already has default ctor in right place')
        else:
            print(f'[SKIP] {class_name} - ctor match not in class definition context')
    else:
        print(f'[WARN] {class_name} - could not find parameterized ctor')

with open(BRIDGE_FILE, 'w', encoding='utf-8') as f:
    f.write(content)

print(f'\nFixed {fixed} classes')
print(f'Written {BRIDGE_FILE} ({original_len} -> {len(content)} chars)')
