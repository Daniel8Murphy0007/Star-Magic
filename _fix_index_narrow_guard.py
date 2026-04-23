"""
Fix the require.main guard in index.js:
1. Remove the overly-broad guard (lines 11483-11484 and line 24005) that incorrectly
   wrapped class definitions needed by the exports block
2. Apply a narrow guard ONLY around the actual demo execution code (lines 11483-11632)

The demo execution code ends after the 'breakthroughs' block (around line 11632).
After that come function/class definitions that must be at module scope.
"""

with open('index.js', encoding='utf-8') as f:
    content = f.read()

lines = content.splitlines(keepends=True)

total = len(lines)
print(f'Total lines: {total}')

# Find markers in current (broken) state
demo_comment_idx = None  # line with "DEMO/ANALYSIS BLOCK"
guard_open_idx = None    # line with "if (require.main === module) {"
guard_close_idx = None   # line with "} // end require.main guard"
demo_exec_end_idx = None # line AFTER last demo execution statement

for i, line in enumerate(lines):
    stripped = line.strip()
    if 'DEMO/ANALYSIS BLOCK' in line and demo_comment_idx is None:
        demo_comment_idx = i
    if 'if (require.main === module) {' in line and guard_open_idx is None:
        guard_open_idx = i
    if '} // end require.main guard' in line:
        guard_close_idx = i
    # Find the last console.log/forEach in the demo execution block
    # The demo exec block ends after line ~11632 (the else clause for breakthroughs)
    # We identify it as: the line matching "No breakthrough thresholds exceeded"
    if 'No breakthrough thresholds exceeded' in line:
        demo_exec_end_idx = i  # last actual execution statement

print(f'demo_comment_idx: {demo_comment_idx+1}')
print(f'guard_open_idx:   {guard_open_idx+1}')
print(f'guard_close_idx:  {guard_close_idx+1}')
print(f'demo_exec_end_idx: {demo_exec_end_idx+1}')

# Verify we can find the closing brace of the if/else at demo_exec_end_idx
# It should be on the same line or shortly after
for i in range(demo_exec_end_idx, min(demo_exec_end_idx+5, total)):
    print(f'  {i+1}: {lines[i].rstrip()[:80]}')

if None in (demo_comment_idx, guard_open_idx, guard_close_idx, demo_exec_end_idx):
    print('ERROR: could not find all markers')
    exit(1)

# Find the closing brace of the if(breakthroughs.length > 0) else block
# Starting from demo_exec_end_idx, find the next closing "}"
close_brace_after_exec = None
for i in range(demo_exec_end_idx, min(demo_exec_end_idx+10, total)):
    if lines[i].strip() == '}':
        close_brace_after_exec = i
        break

print(f'close_brace_after_exec: {close_brace_after_exec+1 if close_brace_after_exec else "NOT FOUND"}')
if close_brace_after_exec:
    for i in range(close_brace_after_exec, min(close_brace_after_exec+5, total)):
        print(f'  {i+1}: {lines[i].rstrip()[:80]}')

# Build new content:
# - Lines 0..demo_comment_idx-1: unchanged (before demo block)
# - New narrow guard: wrap ONLY lines demo_comment_idx..close_brace_after_exec
#   with "if (require.main === module) { ... }"
# - Lines close_brace_after_exec+1..guard_close_idx-1: unchanged (class definitions - NO guard)
# - Line guard_close_idx: SKIP (remove the old broad guard close)
# - Lines guard_close_idx+1..end: unchanged (exports block)

# The broad guard_open was at guard_open_idx. We need to:
# 1. Replace line guard_open_idx with just a blank line (remove the broad guard open)
# 2. Replace line guard_close_idx with just a blank line (remove the broad guard close)
# 3. Insert narrow guard around demo_comment_idx..close_brace_after_exec

# Easier: rebuild in one pass
before = lines[:demo_comment_idx]
demo_exec_block = lines[demo_comment_idx:close_brace_after_exec+1]  # includes last closing }
class_defs_block = lines[close_brace_after_exec+1:guard_close_idx]  # class defs outside guard
after = lines[guard_close_idx+1:]  # exports block

# Remove the old "if (require.main === module) {" from demo_exec_block
# (it's at guard_open_idx which is demo_comment_idx+1 in current state)
relative_guard_idx = guard_open_idx - demo_comment_idx
print(f'\nRelative guard open position in demo block: {relative_guard_idx}')
print(f'  Line: {demo_exec_block[relative_guard_idx].rstrip()[:80]}')

# Remove the guard open line from demo_exec_block and replace with nothing
demo_exec_block_clean = [l for i, l in enumerate(demo_exec_block) if i != relative_guard_idx]

print(f'demo_exec_block: {len(demo_exec_block)} -> {len(demo_exec_block_clean)} lines (removed guard open)')
print(f'class_defs_block: {len(class_defs_block)} lines')
print(f'after (exports): {len(after)} lines')

# Build final content
new_lines = (
    before +
    ['// === DEMO/ANALYSIS: runs only when executed directly (node index.js) ===\n',
     'if (require.main === module) {\n'] +
    demo_exec_block_clean +
    ['} // end demo guard\n',
     '\n'] +
    class_defs_block +
    after
)

new_content = ''.join(new_lines)

# Syntax check
with open('_index_tmp2.js', 'w', encoding='utf-8') as f:
    f.write(new_content)

import subprocess
result = subprocess.run(['node', '--check', '_index_tmp2.js'], capture_output=True, text=True)
if result.returncode == 0:
    print(f'\nSyntax check PASSED')
    with open('index.js', 'w', encoding='utf-8') as f:
        f.write(new_content)
    print(f'Written: {len(new_lines)} lines')
else:
    print(f'\nSyntax check FAILED:')
    print(result.stderr[:800])
    import os; os.remove('_index_tmp2.js')
    exit(1)

import os; os.remove('_index_tmp2.js')
print('Done.')
