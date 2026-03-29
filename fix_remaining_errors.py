#!/usr/bin/env python3
"""Fix remaining 2 error categories before build:
1. Add #undef PI / #undef P guards around wolfram_sources_bridge.cpp include
2. Comment out registrations for undefined hash-suffixed _84A767D3 classes
"""
import re

MAIN_FILE = 'MAIN_1_CoAnQi.cpp'

with open(MAIN_FILE, 'r', encoding='utf-8-sig') as f:
    content = f.read()

original_len = len(content)

# -----------------------------------------------------------------------
# Fix 1: Add #undef PI / #undef P guards around bridge include
# -----------------------------------------------------------------------
OLD_INCLUDE = '#include "wolfram_sources_bridge.cpp"'
NEW_INCLUDE = (
    '// Temporarily undefine macros that conflict with wolfram_sources_bridge.cpp identifiers\n'
    '#ifdef PI\n'
    '#  undef PI\n'
    '#endif\n'
    '#ifdef P\n'
    '#  undef P\n'
    '#endif\n'
    '#include "wolfram_sources_bridge.cpp"\n'
    '// Restore PI macro after bridge include\n'
    '#ifndef PI\n'
    '#  define PI 3.141592653589793\n'
    '#endif'
)

if OLD_INCLUDE in content:
    content = content.replace(OLD_INCLUDE, NEW_INCLUDE, 1)
    print(f'[OK] Added #undef PI/#undef P guards around bridge include')
else:
    print('[SKIP] wolfram_sources_bridge.cpp include not found as plain string (may already be guarded)')

# -----------------------------------------------------------------------
# Fix 2: Comment out _84A767D3 hash-suffixed class registrations
# -----------------------------------------------------------------------
# Match the whole registerPhysicsTerm call spanning 1-3 lines for _84A767D3 classes
pattern = re.compile(
    r'([ \t]*core\.registerPhysicsTerm\("[^"]+",\s*\n?\s*std::make_unique<[A-Za-z0-9:_]+_84A767D3>\(\)[^;]*;)',
    re.MULTILINE
)

def comment_out(m):
    lines = m.group(1).split('\n')
    return '\n'.join('    // REMOVED (undefined class): ' + l.lstrip() for l in lines)

new_content, count = pattern.subn(comment_out, content)
print(f'[OK] Commented out {count} _84A767D3 class registrations')
content = new_content

# -----------------------------------------------------------------------
# Write back
# -----------------------------------------------------------------------
with open(MAIN_FILE, 'w', encoding='utf-8') as f:
    f.write(content)

print(f'Written {MAIN_FILE} ({original_len} -> {len(content)} chars)')
