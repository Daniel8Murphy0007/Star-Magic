with open('CondensedPhysics4.py', 'r', encoding='utf-8', errors='replace') as f:
    lines = f.readlines()
for i in range(10682, 10695):
    print(i+1, repr(lines[i]))
