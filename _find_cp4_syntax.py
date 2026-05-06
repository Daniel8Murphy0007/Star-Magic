import re
with open('CondensedPhysics4.py', 'r', encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

# Find lines where a single-quoted string value has unescaped apostrophes
bad = []
for i, line in enumerate(lines):
    # Pattern: a dict value like 'some text with letter' followed by ( or / — breaks Python string
    if re.search(r"': '.*[a-zA-Z]'[/(]", line):
        bad.append(i+1)
    # Also catch: 'text e'(r) pattern
    elif re.search(r"'[^']*[a-zA-Z]'[/(][^']*'", line) and "': '" in line:
        bad.append(i+1)

print(f'Lines with potential unescaped apostrophes: {len(bad)}')
for ln in bad:
    print(ln, repr(lines[ln-1][:150]))
