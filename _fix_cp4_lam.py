"""Fix CondensedPhysics4.py line 10643 where U+FFFD replacement char is used as variable name.
Context: lam12 = P/3, lam3 = 2P/3, ? = lam12  (? is likely lambda, renamed to 'lam')
Also fix the two f-string usages of this variable on lines 10657 and 10668.
"""
with open('CondensedPhysics4.py', 'r', encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

FFFD = '?'  # literal ASCII question mark (not U+FFFD)

# Fix line 10643 (0-indexed: 10642): replace variable assignment
idx = 10642
old = lines[idx]
new = old.replace(FFFD + '      = lam12', 'lam      = lam12')
if old != new:
    print(f'Fixed line 10643: {repr(old.strip())} -> {repr(new.strip())}')
    lines[idx] = new
else:
    print(f'Line 10643 not matched: {repr(old)}')

# Fix line 10657 (0-indexed: 10656): f-string {?:.2e} -> {lam:.2e}
idx = 10656
old = lines[idx]
# In the f-string, {?:.2e} should become {lam:.2e}
new = old.replace('{' + FFFD + ':.2e}', '{lam:.2e}')
if old != new:
    print(f'Fixed line 10657: {repr(old.strip())} -> {repr(new.strip())}')
    lines[idx] = new
else:
    print(f'Line 10657 not matched: {repr(old)}')

# Fix line 10668 (0-indexed: 10667): f-string {?:.3e} -> {lam:.3e}  
idx = 10667
old = lines[idx]
new = old.replace('{' + FFFD + ':.3e}', '{lam:.3e}')
if old != new:
    print(f'Fixed line 10668: {repr(old.strip())} -> {repr(new.strip())}')
    lines[idx] = new
else:
    print(f'Line 10668 not matched: {repr(old)}')

with open('CondensedPhysics4.py', 'w', encoding='utf-8', errors='replace') as f:
    f.writelines(lines)
print('Done')
