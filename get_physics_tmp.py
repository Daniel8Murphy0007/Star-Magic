import re

def get_methods(fname):
    with open(fname, encoding='utf-8') as f:
        content = f.read()
    # Get all method names from cpp
    methods = re.findall(r'^\w[\w<>:, ]+\s+\w+::(compute\w+|updateVariable|getEquation|computeFBi)\s*\(', content, re.MULTILINE)
    return methods, content

def get_systems(content):
    return re.findall(r'system == "([^"]+)"', content)

print('=== UQFFBuoyancyModule ===')
methods, content = get_methods('UQFFBuoyancyModule.cpp')
print('Methods:', methods[:30])
print('Systems:', get_systems(content))
print()

print('=== UQFFBuoyancyAstroModule ===')
methods, content = get_methods('UQFFBuoyancyAstroModule.cpp')
print('Methods:', methods[:30])
print('Systems:', get_systems(content))
# Get computeFBi logic
idx = content.find('computeFBi')
print()
print('computeFBi (first 1500 chars):')
print(content[idx:idx+1500])
print()

print('=== UQFFBuoyancyCNBModule ===')
methods, content3 = get_methods('UQFFBuoyancyCNBModule.cpp')
print('Methods:', methods[:30])
print('Systems:', get_systems(content3))
# Get computeFBi logic
idx = content3.find('computeFBi')
print()
print('computeFBi CNB (first 2500 chars):')
print(content3[idx:idx+2500])
