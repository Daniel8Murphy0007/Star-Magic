#!/usr/bin/env python3
"""
Add missing getName() pure virtual implementation to ALL sw* bridge classes
that are abstract because they omit getName().

Strategy: Find every getDescription() override in the bridge file that is NOT
preceded by a getName() within the same class body, and insert getName() before it.
Works across all namespaces.
"""
import re

BRIDGE_FILE = 'wolfram_sources_bridge.cpp'

with open(BRIDGE_FILE, 'r', encoding='utf-8') as f:
    content = f.read()

original_len = len(content)
fixed_count = 0

# For each occurrence of `std::string getDescription() const override {`
# check if `getName()` appears within the ~1000 chars before it (class body scope).
# If not, insert a getName() just before it.

# We process in reverse to preserve offsets.
desc_pattern = re.compile(
    r'([ \t]+)(std::string getDescription\(\) const override \{)',
    re.MULTILINE
)

positions = [(m.start(), m.group(1)) for m in desc_pattern.finditer(content)]

# Build new content by collecting (insertion_offset, insertion_text) pairs
insertions = []
for pos, indent in positions:
    preceding = content[max(0, pos - 1200):pos]
    if 'getName()' in preceding:
        continue  # Already has getName
    
    # Try to extract class name from preceding context
    class_match = re.search(r'class\s+(\w+)\s*:', preceding)
    # Also check for namespace context
    ns_match = re.search(r'namespace\s+(\w+)\s*\{', preceding)
    class_name = class_match.group(1) if class_match else f'BridgeTerm_{fixed_count:04d}'
    ns_name = ns_match.group(1) if ns_match else 'unknown'
    
    insertion = f'{indent}std::string getName() const override {{ return "{ns_name}::{class_name}"; }}\n\n'
    insertions.append((pos, insertion))
    fixed_count += 1

# Apply insertions in reverse order to preserve offsets
insertions.sort(key=lambda x: x[0], reverse=True)
for pos, insertion in insertions:
    content = content[:pos] + insertion + content[pos:]

with open(BRIDGE_FILE, 'w', encoding='utf-8') as f:
    f.write(content)

print(f'Total getName() added: {fixed_count}')
print(f'Written {BRIDGE_FILE} ({original_len} -> {len(content)} chars)')
