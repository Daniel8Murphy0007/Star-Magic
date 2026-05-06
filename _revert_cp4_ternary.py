with open('CondensedPhysics4.py', 'r', encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

# Revert lines 9423 and 11530 (1-indexed) back to their original ternary form
# Current (wrong):  'key': "FIRST_PART' if CONDITION else 'SECOND_PART",
# Original:         'key': 'FIRST_PART' if CONDITION else 'SECOND_PART',

def revert_ternary(line):
    """Revert a mistakenly-fixed ternary dict value back to original."""
    import re
    # Pattern: 'key': "value1' if COND else 'value2",
    m = re.match(r"^(\s+'[^']+'\s*:\s+)\"(.+)'( if .+ else )'(.+)\",(\s*)$", line.rstrip('\n'))
    if m:
        key = m.group(1)
        v1 = m.group(2)
        middle = m.group(3)
        v2 = m.group(4)
        tail = m.group(5)
        return key + "'" + v1 + "'" + middle + "'" + v2 + "'," + tail + "\n"
    return None

for idx in [9422, 11529]:  # 0-indexed
    original = lines[idx]
    fixed = revert_ternary(original)
    if fixed:
        print(f'Reverted line {idx+1}:')
        print(f'  Before: {repr(original[:100])}')
        print(f'  After:  {repr(fixed[:100])}')
        lines[idx] = fixed
    else:
        print(f'Could not revert line {idx+1}: {repr(original[:100])}')

with open('CondensedPhysics4.py', 'w', encoding='utf-8', errors='replace') as f:
    f.writelines(lines)
print('Done')
