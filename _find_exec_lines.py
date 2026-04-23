"""
Find all top-level executable (non-definition) statements after line 11482 in index.js.
These are the ones that run on require() and should be wrapped in require.main guard.
"""
import re

with open('index.js', encoding='utf-8') as f:
    lines = f.readlines()

# Find ALL top-level executable statements (not class/function definitions) after line 11482
# Also find the blocks: console.log, systemsToAnalyze.forEach, etc.

exec_lines = []
for i, line in enumerate(lines[11482:], 11483):
    stripped = line.lstrip()
    indent_len = len(line) - len(stripped)
    if indent_len == 0 and stripped.strip():
        # Top-level non-comment
        if stripped.startswith('//') or stripped.startswith('/*') or stripped.startswith('*'):
            continue
        # Skip definitions
        if re.match(r'^(class |function |const |let |var |if |try |catch |finally |async function |"use strict|module\.exports|} )', stripped):
            continue
        # Skip single closing brace
        if stripped.strip() in ('}', '};', '});'):
            continue
        exec_lines.append((i, line.rstrip()[:100]))

print(f'Top-level exec statements after 11482: {len(exec_lines)}')
for ln, text in exec_lines[:40]:
    print(f'  {ln}: {text}')
