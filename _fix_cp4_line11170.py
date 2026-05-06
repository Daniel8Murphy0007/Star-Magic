with open('CondensedPhysics4.py', 'r', encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

old = lines[11169]
print('old:', repr(old))

# The value starts with 'G and ends with aa)', with unescaped ' inside
# Replace the single-quoted string value with a double-quoted one
import re
new = re.sub(
    r"'(G\?_\{[^}]+\} = -e'[^']+)'",
    lambda m: '"' + m.group(1).replace("'", "'") + '"',
    old
)
# Simpler: just fix the specific pattern
new = old.replace(
    "'G\ufffd_{\ufffd\ufffd} = -e'/(2A_rr),  G\ufffd_{ar} = e'/(2A_aa)'",
    '"G\ufffd_{\ufffd\ufffd} = -e\'/(2A_rr),  G\ufffd_{ar} = e\'/(2A_aa)"'
)
if old == new:
    # Try alternative: just escape the apostrophes
    # Find the problematic substring and replace
    prefix = "                'Christoffels': "
    if old.startswith(prefix):
        rest = old[len(prefix):]
        # rest should be: 'G?_... = -e'/(2A_rr),  G?_{ar} = e'/(2A_aa)',\n
        # Fix by changing to double-quoted string
        fixed_val = rest.strip().rstrip(',\n')
        # fixed_val = 'G?_... = -e'/(2A_rr),  G?_{ar} = e'/(2A_aa)'
        # Remove outer single quotes and re-wrap with double quotes
        inner = fixed_val[1:-1]  # strip leading and trailing '
        new = prefix + '"' + inner + '",\n'

print('new:', repr(new))
if old == new:
    print('NO CHANGE')
else:
    lines[11169] = new
    with open('CondensedPhysics4.py', 'w', encoding='utf-8', errors='replace') as f:
        f.writelines(lines)
    print('File patched OK')
