"""
Show the actual long-form equations computed by:
  (A) A QC-correct calculator (uses _RHO_VAC_SCM)
  (B) A hardcoded calculator (uses SCm_density=1e15 or RHO_SCM=9.47e-27 etc.)
Pull the full compute() method body for representative examples from each CP file.
"""
import re

def extract_method(text, start_line, name='compute', lines_max=80):
    """Extract method starting at start_line (1-indexed), return method text."""
    lines = text.splitlines()
    idx = start_line - 1
    # Find the def line
    while idx < len(lines) and f'def {name}' not in lines[idx]:
        idx += 1
    if idx >= len(lines):
        return f"[method {name} not found near line {start_line}]"
    method_lines = [lines[idx]]
    indent = len(lines[idx]) - len(lines[idx].lstrip())
    idx += 1
    count = 1
    while idx < len(lines) and count < lines_max:
        line = lines[idx]
        if line.strip() == '':
            method_lines.append(line)
        elif len(line) - len(line.lstrip()) <= indent and line.strip().startswith('def '):
            break
        else:
            method_lines.append(line)
        idx += 1
        count += 1
    if count >= lines_max:
        method_lines.append("    ... [truncated]")
    return '\n'.join(method_lines)

def find_class_containing_line(text, lineno):
    """Find class name containing a given line number."""
    lines = text.splitlines()
    cls = None
    for i in range(min(lineno-1, len(lines)-1), -1, -1):
        m = re.match(r'^class (\w+)', lines[i])
        if m:
            return m.group(1)
    return None

print("=" * 72)
print("ACTUAL EQUATIONS IN CP FILES — QC-DERIVED vs HARDCODED")
print("=" * 72)

# ---- CP1: QC-correct example ----
with open('CondensedPhysics.py', 'r', encoding='utf-8', errors='replace') as f:
    cp1 = f.read()
cp1_lines = cp1.splitlines()

print("\n\n>>> CP1 — CORRECT: module uses _RHO_VAC_SCM (line 805, rho_ua = _RHO_VAC_UA)")
print("Class context:", find_class_containing_line(cp1, 805))
# Find the compute method of that class
cls_name = find_class_containing_line(cp1, 805)
if cls_name:
    # Find where that class defines compute()
    class_start = None
    for i, ln in enumerate(cp1_lines):
        if f'class {cls_name}' in ln:
            class_start = i
            break
    if class_start:
        # Find compute in that class
        for i in range(class_start, min(class_start+300, len(cp1_lines))):
            if 'def compute(' in cp1_lines[i]:
                print(f"\n[Lines ~{i+1}+]")
                print(extract_method(cp1, i+1, name='compute', lines_max=60))
                break

print("\n\n>>> CP1 — HARDCODED: SCm_density=1e15 (line 17190)")
cls_name2 = find_class_containing_line(cp1, 17190)
print("Class:", cls_name2)
if cls_name2:
    class_start2 = None
    for i, ln in enumerate(cp1_lines):
        if f'class {cls_name2}' in ln:
            class_start2 = i
            break
    if class_start2:
        for i in range(class_start2, min(class_start2+400, len(cp1_lines))):
            if 'def compute(' in cp1_lines[i]:
                print(f"\n[Lines ~{i+1}+]")
                print(extract_method(cp1, i+1, name='compute', lines_max=70))
                break

# ---- CP4: hardcoded example ----
with open('CondensedPhysics4.py', 'r', encoding='utf-8', errors='replace') as f:
    cp4 = f.read()
cp4_lines = cp4.splitlines()

print("\n\n>>> CP4 — HARDCODED: RHO_SCM=9.47e-27 kg/m³ (line 5632)")
cls_name4 = find_class_containing_line(cp4, 5632)
print("Class:", cls_name4)
if cls_name4:
    class_start4 = None
    for i, ln in enumerate(cp4_lines):
        if f'class {cls_name4}' in ln:
            class_start4 = i
            break
    if class_start4:
        for i in range(class_start4, min(class_start4+500, len(cp4_lines))):
            if 'def compute(' in cp4_lines[i]:
                print(f"\n[Lines ~{i+1}+]")
                print(extract_method(cp4, i+1, name='compute', lines_max=80))
                break

print("\n\n>>> CP4 — HARDCODED: RHO_SCM=1.60e19 J/m³ (line 7812)")
cls_name4b = find_class_containing_line(cp4, 7812)
print("Class:", cls_name4b)
if cls_name4b:
    class_start4b = None
    for i, ln in enumerate(cp4_lines):
        if f'class {cls_name4b}' in ln:
            class_start4b = i
            break
    if class_start4b:
        for i in range(class_start4b, min(class_start4b+500, len(cp4_lines))):
            if 'def compute(' in cp4_lines[i]:
                print(f"\n[Lines ~{i+1}+]")
                print(extract_method(cp4, i+1, name='compute', lines_max=80))
                break

# ---- CP3: module-level hardcoded RHO_VAC_A=1e-23 ----
with open('CondensedPhysics3.py', 'r', encoding='utf-8', errors='replace') as f:
    cp3 = f.read()
cp3_lines = cp3.splitlines()

print("\n\n>>> CP3 — HARDCODED local RHO_SCM=1e-26 (line 13137)")
cls_name3 = find_class_containing_line(cp3, 13137)
print("Class:", cls_name3)
if cls_name3:
    class_start3 = None
    for i, ln in enumerate(cp3_lines):
        if f'class {cls_name3}' in ln:
            class_start3 = i
            break
    if class_start3:
        for i in range(class_start3, min(class_start3+400, len(cp3_lines))):
            if 'def compute(' in cp3_lines[i]:
                print(f"\n[Lines ~{i+1}+]")
                print(extract_method(cp3, i+1, name='compute', lines_max=80))
                break
