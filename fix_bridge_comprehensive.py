#!/usr/bin/env python3
"""
Comprehensive fix for wolfram_sources_bridge.cpp:
1. Add 'getDescription() const override { return description; }' to all PhysicsTerm subclasses
   that set 'description' in constructor but don't override getDescription().
2. Add '= default' constructor to all PhysicsTerm subclasses that lack a default constructor
   (needed for std::make_unique<T>() calls with no args).
"""
import re

BRIDGE_FILE = 'wolfram_sources_bridge.cpp'

with open(BRIDGE_FILE, 'r', encoding='utf-8') as f:
    content = f.read()

original_len = len(content)

# -----------------------------------------------------------------------
# Find all PhysicsTerm subclass bodies
# -----------------------------------------------------------------------
# Pattern to find class declarations inheriting from PhysicsTerm (directly or ::PhysicsTerm)
class_pattern = re.compile(
    r'class\s+(\w+)\s*:\s*public\s*::?PhysicsTerm\b',
    re.MULTILINE
)

def get_class_body(content, class_start):
    """Returns (body_start, body_end) as character offsets for the class body."""
    # Find opening brace
    i = class_start
    while i < len(content) and content[i] != '{':
        i += 1
    if i >= len(content):
        return None, None
    brace_depth = 0
    start = i
    while i < len(content):
        if content[i] == '{':
            brace_depth += 1
        elif content[i] == '}':
            brace_depth -= 1
            if brace_depth == 0:
                return start, i
        i += 1
    return None, None

# Collect all classes with their body ranges
classes = []
for m in class_pattern.finditer(content):
    name = m.group(1)
    body_start, body_end = get_class_body(content, m.start())
    if body_start is not None:
        classes.append({
            'name': name,
            'class_start': m.start(),
            'body_start': body_start,
            'body_end': body_end,
        })

# -----------------------------------------------------------------------
# Determine what each class needs
# -----------------------------------------------------------------------
needs_getdesc = []     # (class_name, insert_before_closing_brace_offset)
needs_default_ctor = []

for cls in classes:
    body = content[cls['body_start']:cls['body_end'] + 1]
    name = cls['name']
    closing = cls['body_end']
    
    # Check for getDescription override
    has_getdesc = 'getDescription()' in body
    sets_description = 'description =' in body
    
    if not has_getdesc and sets_description:
        needs_getdesc.append((name, closing))
    
    # Check for default constructor (needed for make_unique<T>())
    # A class needs a default ctor if it ONLY has parameterized constructors
    # Pattern: has constructor declaration(s) but none with empty params
    ctor_pattern = re.compile(rf'\b{re.escape(name)}\s*\(([^)]*)\)')
    ctor_matches = ctor_pattern.findall(body)
    
    if ctor_matches:
        # Filter out destructor (starts with ~)
        has_empty_ctor = any(
            p.strip() == '' or                          # no params
            all(                                         # all params have defaults
                '=' in param
                for param in p.split(',')
                if param.strip()
            )
            for p in ctor_matches
        )
        if not has_empty_ctor:
            needs_default_ctor.append((name, cls['body_start']))  # insert after {

print(f'Classes needing getDescription(): {len(needs_getdesc)}')
for name, _ in needs_getdesc:
    print(f'  {name}')

print(f'\nClasses needing default constructor: {len(needs_default_ctor)}')
for name, _ in needs_default_ctor:
    print(f'  {name}')

# -----------------------------------------------------------------------
# Apply fixes in reverse order (to preserve offsets)
# -----------------------------------------------------------------------
# Collect all insertions as (offset, text)
insertions = []

for name, closing_brace_pos in needs_getdesc:
    insert_text = '\n    std::string getDescription() const override { return description; }\n'
    insertions.append((closing_brace_pos, insert_text))

for name, open_brace_pos in needs_default_ctor:
    # Insert default ctor right after '{'
    insert_text = f'\n    {name}() = default;'
    insertions.append((open_brace_pos + 1, insert_text))

# Sort by offset descending to apply without offset shifts
insertions.sort(key=lambda x: x[0], reverse=True)

# Deduplicate (same offset shouldn't be inserted twice)
seen_offsets = set()
unique_insertions = []
for offset, text in insertions:
    if offset not in seen_offsets:
        seen_offsets.add(offset)
        unique_insertions.append((offset, text))

applied = 0
for offset, text in unique_insertions:
    # Check if text is already there (avoid duplicates on re-run)
    surrounding = content[max(0, offset-50):min(len(content), offset+100)]
    check = text.strip().split('\n')[0].strip() if '\n' in text else text.strip()
    if check and check in content[max(0, offset-200):min(len(content), offset+200)]:
        continue  # Already there
    content = content[:offset] + text + content[offset:]
    applied += 1

print(f'\nApplied {applied} fixes')
print(f'Written {BRIDGE_FILE} ({original_len} -> {len(content)} chars)')

with open(BRIDGE_FILE, 'w', encoding='utf-8') as f:
    f.write(content)
