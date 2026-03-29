#!/usr/bin/env python3
"""
Fix abstract classes that set 'description' in constructor but don't override getDescription().
Adds: std::string getDescription() const override { return description; }
just before the closing brace of each affected class.

Also handles sw7comp::CompressedBaseTerm which has no default constructor issue
(it uses a non-default constructor, so std::make_unique<T>() fails).
"""
import re

BRIDGE_FILE = 'wolfram_sources_bridge.cpp'

# Classes to fix - all missing getDescription() override
CLASSES_MISSING_GETDESC = [
    'YAMLConfigLoaderTerm',
    'ArchiveMediaManagerTerm',
    'CompressedBaseTerm',
    'CompressedCosmTerm',
    'CompressedEnvTerm',
    'CompressedExpansionTerm',
    'CompressedFluidTerm',
    'CompressedPerturbationTerm',
    'CompressedQuantumTerm',
    'CompressedSuperAdjTerm',
    'CompressedUgSumTerm',
    'ResonanceAAetherResTerm',
    'ResonanceADPMTerm',
    'ResonanceAQuantumFreqTerm',
    'ResonanceASuperFreqTerm',
    'ResonanceATHzTerm',
    'ResonanceAvacDiffTerm',
    'ResonanceUg4iTerm',
]

with open(BRIDGE_FILE, 'r', encoding='utf-8') as f:
    content = f.read()

original_len = len(content)
fixed_count = 0

for class_name in CLASSES_MISSING_GETDESC:
    # Find the class definition
    class_pattern = re.compile(
        rf'(class\s+{re.escape(class_name)}\s*:[^\{{]+\{{)',
        re.MULTILINE
    )
    m = class_pattern.search(content)
    if not m:
        print(f'[SKIP] class {class_name} not found')
        continue

    class_start = m.start()
    # Find the matching closing brace of this class
    brace_depth = 0
    i = m.start(0)
    in_string = False
    class_end = -1
    while i < len(content):
        ch = content[i]
        if ch == '{':
            brace_depth += 1
        elif ch == '}':
            brace_depth -= 1
            if brace_depth == 0:
                class_end = i
                break
        i += 1

    if class_end == -1:
        print(f'[SKIP] could not find closing brace for {class_name}')
        continue

    class_body = content[class_start:class_end + 1]

    # Check if getDescription is already there
    if 'getDescription()' in class_body:
        print(f'[SKIP] {class_name} already has getDescription()')
        continue

    # Check if getName is there
    if 'getName()' not in class_body:
        print(f'[WARN] {class_name} is also missing getName() - adding both')
        # Find the closing brace position and insert before it
        insert_at = class_end
        insert_text = '\n    std::string getName() const override { return "' + class_name + '"; }\n    std::string getDescription() const override { return description; }\n'
    else:
        # Just add getDescription
        insert_at = class_end
        insert_text = '\n    std::string getDescription() const override { return description; }\n'

    content = content[:insert_at] + insert_text + content[insert_at:]
    fixed_count += 1
    print(f'[OK] Fixed {class_name}')

# Also fix CompressedBaseTerm default constructor issue if needed
# sw7comp::CompressedBaseTerm has error C2512 (no default constructor)
# Check if it has a constructor with required args
comp_base_match = re.search(r'class\s+CompressedBaseTerm[^{]+\{([^}]+(?:\{[^}]*\}[^}]*)*)', content)
if comp_base_match:
    body_preview = comp_base_match.group(1)[:400]
    if 'CompressedBaseTerm()' not in body_preview and 'CompressedBaseTerm(' in body_preview:
        # Has non-default constructor but no default - add one
        # Find the constructor
        ctor_match = re.search(r'([ \t]+CompressedBaseTerm\([^)]+\))', content)
        if ctor_match:
            indent = '    '
            insert_pos = ctor_match.start()
            default_ctor = f'{indent}CompressedBaseTerm() = default;\n'
            if 'CompressedBaseTerm() = default' not in content[max(0,insert_pos-200):insert_pos+200]:
                content = content[:insert_pos] + default_ctor + content[insert_pos:]
                print(f'[OK] Added default constructor for CompressedBaseTerm')

with open(BRIDGE_FILE, 'w', encoding='utf-8') as f:
    f.write(content)

print(f'\nTotal fixes: {fixed_count}')
print(f'Written {BRIDGE_FILE} ({original_len} -> {len(content)} chars)')
