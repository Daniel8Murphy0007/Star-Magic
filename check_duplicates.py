"""Check for duplicate Calculator classes between CondensedPhysics.py and CondensedPhysics2.py"""
import re
import hashlib

def get_class_content(text, class_name):
    """Extract full class definition"""
    pattern = rf'^(class\s+{class_name}.*?)(?=^class |\Z)'
    match = re.search(pattern, text, re.MULTILINE | re.DOTALL)
    if match:
        return match.group(1).strip()
    return None

def get_docstring(content):
    """Extract docstring from class"""
    match = re.search(r'"""(.+?)"""', content, re.DOTALL)
    return match.group(1).strip()[:100] if match else ''

# Load files
with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
    cp1 = f.read()
with open('CondensedPhysics2.py', 'r', encoding='utf-8') as f:
    cp2 = f.read()

# Extract class names
cp1_classes = set(re.findall(r'^class\s+(\w+Calculator)', cp1, re.MULTILINE))
cp2_classes = set(re.findall(r'^class\s+(\w+Calculator)', cp2, re.MULTILINE))
duplicates = sorted(cp1_classes & cp2_classes)

print(f"CP1 classes: {len(cp1_classes)}")
print(f"CP2 classes: {len(cp2_classes)}")
print(f"Duplicates: {len(duplicates)}")
print()

identical_count = 0
different_count = 0

for dup in duplicates:
    cp1_content = get_class_content(cp1, dup)
    cp2_content = get_class_content(cp2, dup)
    
    if cp1_content and cp2_content:
        cp1_hash = hashlib.md5(cp1_content.encode()).hexdigest()[:8]
        cp2_hash = hashlib.md5(cp2_content.encode()).hexdigest()[:8]
        identical = cp1_hash == cp2_hash
        
        if identical:
            identical_count += 1
            status = "IDENTICAL"
        else:
            different_count += 1
            status = "DIFFERENT"
        
        cp1_doc = get_docstring(cp1_content)
        cp2_doc = get_docstring(cp2_content)
        
        print(f"{dup}: {status}")
        print(f"  CP1 hash: {cp1_hash}")
        print(f"  CP2 hash: {cp2_hash}")
        if not identical:
            print(f"  CP1 docstring: {cp1_doc[:80]}...")
            print(f"  CP2 docstring: {cp2_doc[:80]}...")
        print()

print("=" * 60)
print(f"IDENTICAL implementations: {identical_count}")
print(f"DIFFERENT implementations (potential new physics): {different_count}")
print()

if identical_count > 0:
    print("RECOMMENDATION: Remove or rename identical duplicates in CP2")
