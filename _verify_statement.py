import re, subprocess, sys

print("=" * 65)
print("FULL PIPELINE STATE VERIFICATION")
print("=" * 65)

# 1. Repo file presence
import os
manifolds = ['scm_vacuum_manifold.py', 'ua_vacuum_manifold.py', 'dpm_vacuum_manifold.py']
print("\n--- 1. Manifold files in repo ---")
for f in manifolds:
    exists = os.path.isfile(f)
    size   = os.path.getsize(f) if exists else 0
    print(f"  {f}: exists={exists}  size={size:,} bytes")

# 2. Hardcoded constants in manifolds
print("\n--- 2. Hardcoded 7.09e-37/36 hits in manifold files ---")
for f in manifolds:
    with open(f, 'r', encoding='utf-8') as fh:
        txt = fh.read()
    hits = re.findall(r'7\.09e-3[67][^\n]*', txt)
    qc   = txt.count('derive_from_quantum_chain')
    print(f"  {f}: hardcoded_hits={len(hits)}  derive_QC_calls={qc}")
    for h in hits:
        print(f"    -> {h.strip()[:90]}")

# 3. CP file imports — full-file grep (not just first 8000 chars)
print("\n--- 3. Import block presence in CP files (full-file scan) ---")
cp_files = {
    'CP1': 'CondensedPhysics.py',
    'CP2': 'CondensedPhysics2.py',
    'CP3': 'CondensedPhysics3.py',
    'CP4': 'CondensedPhysics4.py',
}
for label, fname in cp_files.items():
    with open(fname, 'r', encoding='utf-8', errors='replace') as fh:
        txt = fh.read()
    scm_imp = 'from scm_vacuum_manifold import' in txt
    ua_imp  = 'from ua_vacuum_manifold import' in txt
    dpm_imp = 'from dpm_vacuum_manifold import' in txt
    ua_flag = '_UA_MANIFOLD_LOADED' in txt
    dpm_flag= '_DPM_MANIFOLD_LOADED' in txt
    # Find line number of UA import if present
    ua_line = None
    if ua_imp:
        for i, line in enumerate(txt.splitlines(), 1):
            if 'from ua_vacuum_manifold import' in line:
                ua_line = i
                break
    print(f"  {label}: scm={scm_imp}  ua={ua_imp}(L{ua_line})  dpm={dpm_imp}  ua_flag={ua_flag}  dpm_flag={dpm_flag}")

# 4. Runtime import test
print("\n--- 4. Runtime import test (all 3 manifolds) ---")
res = subprocess.run(
    [sys.executable, '-X', 'utf8', '-c',
     'from scm_vacuum_manifold import derive_from_quantum_chain, RHO_VAC_SCM;'
     'from ua_vacuum_manifold  import ua_dpm_total_density, DPM_DENSITY_RATIO;'
     'from dpm_vacuum_manifold import ELEMENT, E_CRACK, DPM_DENSITY_RATIO as dr;'
     'rho,_=derive_from_quantum_chain();'
     'print(f"SCm RHO_VAC_SCM={RHO_VAC_SCM:.4e}");'
     'print(f"UA  total_density={ua_dpm_total_density(0):.4e}  ratio={DPM_DENSITY_RATIO}");'
     'print(f"DPM E_CRACK={E_CRACK:.4e}  ratio={dr:.4f}  Z1={ELEMENT[1].symbol}");'
     'print("All 3 manifolds importable: PASS")'],
    capture_output=True, text=True
)
print(res.stdout.strip())
if res.returncode != 0:
    print("STDERR:", res.stderr[:300])

print("\n--- 5. scm hardcoded value context ---")
with open('scm_vacuum_manifold.py', 'r', encoding='utf-8') as fh:
    scm_lines = fh.readlines()
for i, ln in enumerate(scm_lines):
    if '7.09e-3' in ln:
        start = max(0, i-2)
        end   = min(len(scm_lines), i+3)
        print(f"  Line {i+1}: {ln.rstrip()}")
        for j in range(start, end):
            print(f"    {j+1}: {scm_lines[j].rstrip()}")
        print()

print("\n--- 6. dpm hardcoded values (first 5) ---")
with open('dpm_vacuum_manifold.py', 'r', encoding='utf-8') as fh:
    dpm_lines = fh.readlines()
count = 0
for i, ln in enumerate(dpm_lines):
    if re.search(r'7\.09e-3[67]', ln):
        print(f"  Line {i+1}: {ln.rstrip()}")
        count += 1
        if count >= 5:
            print("  ...")
            break
