import re

files = {
    'scm': 'scm_vacuum_manifold.py',
    'ua':  'ua_vacuum_manifold.py',
    'dpm': 'dpm_vacuum_manifold.py',
}

contents = {}
for key, fname in files.items():
    with open(fname, 'r', encoding='utf-8') as f:
        contents[key] = f.read()

print("=== Hardcoded 7.09e-37 / 7.09e-36 occurrences ===")
for key, txt in contents.items():
    hits = re.findall(r'7\.09e-3[67]', txt)
    print(f"  {key}: {len(hits)} hits")

print("\n=== derive_from_quantum_chain occurrences ===")
for key, txt in contents.items():
    print(f"  {key}: {txt.count('derive_from_quantum_chain')}")

print("\n=== Key exports present in each file ===")
for key, fname in files.items():
    txt = contents[key]
    checks = ['RHO_VAC_SCM', 'RHO_VAC_UA', 'DPM_DENSITY_RATIO', 'derive_from_quantum_chain']
    for c in checks:
        print(f"  {key}: {c!r} present = {c in txt}")
