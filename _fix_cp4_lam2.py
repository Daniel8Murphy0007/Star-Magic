"""Fix remaining ? variable usages in CP4 — replace ? as Python identifier with lam."""
with open('CondensedPhysics4.py', 'r', encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

# Fix line 10687: "YM_gap": ?,  -> "YM_gap": lam,
idx = 10686
old = lines[idx]
new = old.replace('"YM_gap":         ?,', '"YM_gap":         lam,')
if old != new:
    print(f'Fixed line 10687: {repr(old.strip())} -> {repr(new.strip())}')
    lines[idx] = new
else:
    print(f'Line 10687 not matched: {repr(old)}')

# Also scan for any other ? used as identifiers (not in strings/comments)
# Look for patterns like: ?, or ?. or ?: or = ? or (? at start of expression
import re
other_issues = []
for i, line in enumerate(lines):
    # Skip if it's already been fixed
    if i == 10686:
        continue
    # Look for ? used as a Python identifier (preceded and followed by non-string chars)
    # This is tricky - let's just check for ?, and ?: patterns that indicate dict values
    if re.search(r'(?<!["\'])\?[,)]', line) and not line.strip().startswith('#'):
        other_issues.append((i+1, repr(line.strip()[:100])))

if other_issues:
    print(f'\nOther potential ? identifier usages ({len(other_issues)}):')
    for ln, r in other_issues[:20]:
        print(f'  {ln}: {r}')
else:
    print('\nNo other ? identifier issues found')

with open('CondensedPhysics4.py', 'w', encoding='utf-8', errors='replace') as f:
    f.writelines(lines)
print('Done')
