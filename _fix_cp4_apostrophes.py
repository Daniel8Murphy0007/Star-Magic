import re

with open('CondensedPhysics4.py', 'r', encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

fixed_count = 0
# Lines we already know about and need to fix
# Pattern: dict value is a single-quoted string containing unescaped apostrophes
# We look for lines like:  'key': 'text with e'(... or e'/...
# Strategy: detect lines where between ': ' and the apparent end ',\n'
# there are internal single quotes followed by ( or / or space-then-letter

target_lines = [11257, 11258]  # 1-indexed

for ln in target_lines:
    idx = ln - 1
    line = lines[idx]
    print(f'Before [{ln}]: {repr(line[:120])}')
    
    # Find the dict value part: after ': '
    # The value is everything from ': ' to end, which is a broken single-quoted string
    # We want to re-wrap the FULL INTENDED value in double quotes
    
    # Approach: find the pattern '  'value',  where value may contain '
    # The full intended value starts after ': ' and ends before ',\n' or '}\n'
    m = re.match(r"^(\s+'[^']+'\s*:\s+)'(.+)',(\s*)$", line)
    if m:
        key_part = m.group(1)
        val_inner = m.group(2)
        tail = m.group(3)
        # val_inner may still have issues, but at least we wrap in double quotes
        new_line = key_part + '"' + val_inner + '",' + tail + '\n'
        lines[idx] = new_line
        fixed_count += 1
        print(f'After  [{ln}]: {repr(new_line[:120])}')
    else:
        # Manual approach: find ': ' position and replace single-quoted value
        # with double-quoted version by taking everything from ': ' to last ',\n'
        colon_pos = line.find(": '")
        if colon_pos != -1:
            prefix = line[:colon_pos + 2]  # up to and including ': '
            rest = line[colon_pos + 3:]    # after the opening '
            # Find the last ', at end (stripping \n)
            stripped = rest.rstrip()
            if stripped.endswith("',"):
                inner = stripped[:-2]  # remove trailing ',
            elif stripped.endswith("'"):
                inner = stripped[:-1]
            else:
                inner = stripped
            new_line = prefix[:-1] + '"' + inner + '",\n'
            lines[idx] = new_line
            fixed_count += 1
            print(f'After  [{ln}]: {repr(new_line[:120])}')
        else:
            print(f'  Could not fix line {ln}')

print(f'\nFixed {fixed_count} lines')

with open('CondensedPhysics4.py', 'w', encoding='utf-8', errors='replace') as f:
    f.writelines(lines)
print('File written')
