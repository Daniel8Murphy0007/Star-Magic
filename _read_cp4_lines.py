with open('CondensedPhysics4.py', 'r', encoding='utf-8', errors='replace') as f:
    lines = f.readlines()
# Print lines 11254-11260 (0-indexed: 11253-11259)
for i in range(11253, 11261):
    print(i+1, repr(lines[i]))
