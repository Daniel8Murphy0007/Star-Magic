import re
t = open('CondensedPhysics2.py', 'r', encoding='utf-8').read()
classes = re.findall(r'^class \w+.*:', t, re.MULTILINE)
print(f'CP2 classes: {len(classes)}')
print(f'CP2 lines: {len(t.splitlines())}')
