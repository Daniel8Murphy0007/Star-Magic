#!/usr/bin/env python3
"""Quick audit of the 4 QCalc files."""
import re, os

files = {
    'QCalc_API.py': open('QCalc_API.py', 'r', encoding='utf-8', errors='replace').read(),
    'QCalc_Performance.py': open('QCalc_Performance.py', 'r', encoding='utf-8', errors='replace').read(),
    'QCalc_Advanced.py': open('QCalc_Advanced.py', 'r', encoding='utf-8', errors='replace').read(),
    'QCalc_Wolfram_Extensions.py': open('QCalc_Wolfram_Extensions.py', 'r', encoding='utf-8', errors='replace').read(),
}

for fn, content in files.items():
    lines = content.splitlines()
    print(f"\n{'='*60}")
    print(f"FILE: {fn}  ({len(lines)} lines)")
    print(f"{'='*60}")

    # Encoding corruption
    corruption = len(re.findall(r'[â\x80-\x9f\xc2-\xc3]', content))
    print(f"  Encoding artifacts: {corruption}")

    # Import issues
    imports = [l for l in lines if l.startswith('from ') or l.startswith('import ')]
    print(f"  Imports: {len(imports)}")
    for imp in imports:
        print(f"    {imp}")

    # UQFF constant checks
    for pat, desc in [
        (r'beta_i\s*=\s*(0\.5|0\.6)[^0-9]', 'WRONG beta_i (not 0.603)'),
        (r'eta\s*=\s*1e-23', 'WRONG eta (should be 1e-22)'),
        (r'gamma\s*=\s*0\.0001', 'WRONG gamma (should be 5e-5)'),
        (r'kappa\s*=\s*0\.001', 'WRONG kappa (should be 5e-4)'),
        (r'GM\s*/\s*r\s*\*\*\s*2', 'Newtonian GM/r^2 as base'),
    ]:
        hits = [(i+1, l) for i,l in enumerate(lines) if re.search(pat, l)]
        if hits:
            print(f"  [ISSUE] {desc}:")
            for n,l in hits[:5]:
                print(f"    {n}: {l.strip()}")

    # Check for missing cos(pi*t_n)
    trz = len([l for l in lines if re.search(r'cos.*pi.*t_n|np\.cos.*np\.pi.*t_n', l)])
    ereact = len([l for l in lines if 'E_react' in l])
    print(f"  cos(pi*t_n) occurrences: {trz}")
    print(f"  E_react occurrences: {ereact}")

# Check missing dpm_helpers
print(f"\n  dpm_helpers.py exists: {os.path.exists('dpm_helpers.py')}")
