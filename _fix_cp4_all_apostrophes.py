"""
Fix all lines in CondensedPhysics4.py where a single-quoted Python string value
(in a dict) contains unescaped apostrophes that cause SyntaxError.

Strategy:
  For each line matching: <indent>'key'<whitespace>:<whitespace>'value_with_apostrophe',
  convert the value from single-quoted to double-quoted string.
  
  We use a greedy regex that can match across internal apostrophes.
"""
import re

with open('CondensedPhysics4.py', 'r', encoding='utf-8', errors='replace') as f:
    content = f.read()

lines = content.splitlines(keepends=True)
fixed = 0

# Pattern: an indented dict entry with a single-quoted VALUE that contains internal apostrophes
# We detect: line has a ':' followed by spaces then a single-quoted string that contains '
# 
# Specifically we look for the pattern where between the opening ' and closing ',
# there is at least one ' followed by a non-' character (i.e., an internal apostrophe)
#
# We match the FULL intended string value by being greedy and taking everything
# from the first ' after ': ' to the last ' before the final ',

pattern = re.compile(
    r"^(\s+'[^']+'\s*:\s+)'((?:[^'\\]|\\.|'(?=\S))*(?:[a-zA-Z\d\ufffd\u00b2\u00b7])'[^\s,][^']*|[^']*[a-zA-Z\d]'[/(= ][^']*)',(\s*)\n",
    re.MULTILINE
)

# Simpler: just look for lines where the value string APPEARS to contain an internal '
# The greedy approach: find all occurrences of  ': '...','  where ... contains a '
# 
# Revised approach: scan line by line and detect lines where Python would see
# the string as ending prematurely.
# 
# Detection: line matches  'key': 'value',\n  but value contains unescaped '
# We use: match the outer structure, extract the raw value, check if it has '

line_pattern = re.compile(
    r"^(\s+'[^'\n]+'\s*:\s+)'(.+)',(\s*)\n?$"
)

new_lines = []
for i, line in enumerate(lines):
    m = line_pattern.match(line)
    if m:
        key_part = m.group(1)   # '  'keyname':   '
        val_inner = m.group(2)  # the greedy-matched inner content
        tail = m.group(3)       # trailing whitespace
        
        # Check if this value actually contains an apostrophe (meaning it was broken)
        if "'" in val_inner:
            # The inner content contains apostrophes - this was a broken string
            # Re-wrap in double quotes, escaping any existing double quotes
            safe_inner = val_inner.replace('"', '\\"')
            new_line = key_part + '"' + safe_inner + '",' + '\n'
            new_lines.append(new_line)
            if new_line != line:
                fixed += 1
                if fixed <= 20:
                    print(f'Fixed line {i+1}: {repr(line[:80])} -> {repr(new_line[:80])}')
            continue
    new_lines.append(line)

print(f'\nTotal fixed: {fixed}')

with open('CondensedPhysics4.py', 'w', encoding='utf-8', errors='replace') as f:
    f.writelines(new_lines)
print('File written.')
