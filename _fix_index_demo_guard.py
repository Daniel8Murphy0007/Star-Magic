"""
Add require.main === module guard around demo/execution code in index.js.
Also add try/catch around require('./sourceNN.js') blocks in exports.

Demo code runs from line 11483 to ~24002 (before exports block).
Exports block starts at line 24003.

This prevents the demo analysis from running when index.js is require()'d as a library.
"""

with open('index.js', encoding='utf-8') as f:
    content = f.read()

lines = content.splitlines(keepends=True)
total = len(lines)
print(f'Total lines: {total}')

# Find demo start: "// Enhanced Demonstration with Predefined Systems"
demo_start_idx = None
for i, line in enumerate(lines):
    if '// Enhanced Demonstration with Predefined Systems' in line:
        demo_start_idx = i
        print(f'Demo start at line {i+1}: {line.rstrip()[:80]}')
        break

# Find exports block start: "// Export all modules"
exports_start_idx = None
for i, line in enumerate(lines):
    if '// Export all modules' in line and i > 20000:
        exports_start_idx = i
        print(f'Exports start at line {i+1}: {line.rstrip()[:80]}')
        break

if demo_start_idx is None or exports_start_idx is None:
    print('ERROR: could not find markers')
    exit(1)

# Verify the exports block line
print(f'\nLines around exports block:')
for i in range(exports_start_idx, min(exports_start_idx+5, total)):
    print(f'  {i+1}: {lines[i].rstrip()[:80]}')

print(f'\nLines around demo start:')
for i in range(max(0, demo_start_idx-2), min(demo_start_idx+4, total)):
    print(f'  {i+1}: {lines[i].rstrip()[:80]}')

# Build new content:
# - Lines 0..demo_start_idx-1: unchanged
# - Insert: "if (require.main === module) {\n"
# - Lines demo_start_idx..exports_start_idx-1: indented + blank line separators
# - Insert: "}\n\n"
# - Lines exports_start_idx..end: unchanged

# Actually indenting 12520 lines is risky - just wrap with guard block WITHOUT indenting
# Node.js doesn't require indentation for correctness; the guard works without it.

before = ''.join(lines[:demo_start_idx])
demo_block = ''.join(lines[demo_start_idx:exports_start_idx])
after = ''.join(lines[exports_start_idx:])

new_content = (
    before +
    '// === DEMO/ANALYSIS BLOCK: only runs when executed directly, not when require()\'d ===\n'
    'if (require.main === module) {\n' +
    demo_block +
    '} // end require.main guard\n\n' +
    after
)

# Also fix the require() calls in exports block - wrap each in try/catch
# Pattern: const { XYZ } = require('./sourceNN.js');
#          module.exports.XYZ = XYZ;
# Becomes: try { const { XYZ } = require('./sourceNN.js'); module.exports.XYZ = XYZ; } catch(e) { /* optional fallback */ }

import re

def wrap_require_exports(text):
    """Wrap each source require+export pair in try/catch."""
    # Pattern: two consecutive lines inside the if(typeof module...) block
    #   const { NAME } = require('./sourceNN.js');
    #   module.exports.NAME = NAME;
    # Optional blank line between them
    pat = re.compile(
        r'(    const \{ (\w+) \} = require\(\'\.\/source\d+\.js\'\);\n)'
        r'(    module\.exports\.\2 = \2;\n)',
        re.MULTILINE
    )
    
    def replacer(m):
        name = m.group(2)
        return (
            f'    try {{ {m.group(1).strip()} {m.group(3).strip()} }}'
            f' catch(e) {{ /* {name} module not found */ }}\n'
        )
    
    result, n = pat.subn(replacer, text)
    print(f'Wrapped {n} require/export pairs in try/catch')
    return result

new_content = wrap_require_exports(new_content)

# Syntax check via node --check by writing to temp
with open('_index_tmp.js', 'w', encoding='utf-8') as f:
    f.write(new_content)

import subprocess
result = subprocess.run(['node', '--check', '_index_tmp.js'], capture_output=True, text=True)
if result.returncode == 0:
    print('\nSyntax check PASSED')
    with open('index.js', 'w', encoding='utf-8') as f:
        f.write(new_content)
    new_lines = new_content.count('\n')
    print(f'Written: {new_lines} lines')
else:
    print(f'\nSyntax check FAILED:')
    print(result.stderr[:500])
    import os
    os.remove('_index_tmp.js')
    exit(1)

import os
os.remove('_index_tmp.js')
print('Done.')
