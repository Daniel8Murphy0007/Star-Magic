#!/usr/bin/env python
"""Quick analysis script for CondensedPhysics.py document coverage."""
import re

with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
    content = f.read()
    lines = content.split('\n')

print("=" * 80)
print("CONDENSEDPHYSICS.PY DEEP ANALYSIS REPORT")
print("=" * 80)

# 1. Basic counts
classes = re.findall(r'^class (\w+)', content, re.MULTILINE)
compute_methods = re.findall(r'def (compute\w*)', content)
validate_methods = re.findall(r'def (validate\w*|run_tests?\w*)', content)

print(f"\n1. BASIC METRICS:")
print(f"   Total Lines: {len(lines)}")
print(f"   Total Classes: {len(classes)}")
print(f"   Compute Methods: {len(compute_methods)}")
print(f"   Validate/Test Methods: {len(validate_methods)}")

# 2. Document references
docs_found = set()
for match in re.finditer(r'Document (\d+)', content):
    docs_found.add(int(match.group(1)))

print(f"\n2. DOCUMENT REFERENCES:")
print(f"   Found: {sorted(docs_found)}")
missing = set(range(1, 30)) - docs_found
print(f"   Missing (1-29): {sorted(missing) if missing else 'None'}")

# 3. UQFF Model mapping (Documents 21-29)
print(f"\n3. UQFF MODEL TO DOCUMENT MAPPING:")
for i, line in enumerate(lines):
    if line.strip().startswith('class UQFF') and 'Model:' in line:
        # Get next 5 lines for docstring
        docstring = '\n'.join(lines[i:i+10])
        doc_match = re.search(r'Document (\d+)', docstring)
        doc_num = doc_match.group(1) if doc_match else 'NOT FOUND'
        class_name = re.search(r'class (\w+)', line).group(1)
        print(f"   Line {i+1}: {class_name} -> Document {doc_num}")

# 4. Global instances
print(f"\n4. GLOBAL MODEL INSTANCES:")
global_instances = re.findall(r'^([A-Z_]+_MODEL)\s*=\s*(\w+)\(\)', content, re.MULTILINE)
for inst, cls in global_instances:
    print(f"   {inst} = {cls}()")

# 5. Equation constants (self.X = numeric)
eq_constants = re.findall(r'self\.(\w+)\s*=\s*([\d\.e\-\+]+)\s*#?\s*(.*?)$', content, re.MULTILINE)
print(f"\n5. EQUATION CONSTANTS (self.x = value):")
print(f"   Total: {len(eq_constants)}")

# 6. Key physics terms
key_terms = ['SSq', 'beta_i', 'kappa', 'k_eta', 'rho_vac', 'Ug1', 'Ug2', 'Ug3', 'Ug4', 'H_SCm', 'U_UA']
print(f"\n6. KEY PHYSICS TERMS PRESENCE:")
for term in key_terms:
    count = len(re.findall(rf'self\.{term}', content))
    print(f"   self.{term}: {count} occurrences")

# 7. Master equations
print(f"\n7. MASTER EQUATION CLASSES:")
master_classes = [c for c in classes if 'Master' in c or 'Triadic' in c]
for c in master_classes:
    print(f"   {c}")

print("\n" + "=" * 80)
print("END OF ANALYSIS")
print("=" * 80)
