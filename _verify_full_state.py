import re, ast

# --- Check CP files for UA/DPM import blocks ---
cp_files = {
    'CP1': 'CondensedPhysics.py',
    'CP2': 'CondensedPhysics2.py',
    'CP3': 'CondensedPhysics3.py',
    'CP4': 'CondensedPhysics4.py',
}

print("=== UA/DPM import blocks in CP files ===")
for key, fname in cp_files.items():
    with open(fname, 'r', encoding='utf-8', errors='replace') as f:
        head = f.read(8000)  # read first ~8000 chars
    ua_loaded  = '_UA_MANIFOLD_LOADED' in head
    dpm_loaded = '_DPM_MANIFOLD_LOADED' in head
    ua_import  = 'from ua_vacuum_manifold import' in head
    dpm_import = 'from dpm_vacuum_manifold import' in head
    scm_import = 'from scm_vacuum_manifold import' in head
    print(f"  {key}: scm_import={scm_import}  ua_import={ua_import}  dpm_import={dpm_import}  ua_loaded_flag={ua_loaded}  dpm_loaded_flag={dpm_loaded}")

print("\n=== Syntax check all 4 CP files ===")
import py_compile
for key, fname in cp_files.items():
    try:
        py_compile.compile(fname, doraise=True)
        print(f"  {key}: syntax OK")
    except py_compile.PyCompileError as e:
        print(f"  {key}: SYNTAX ERROR - {e}")

print("\n=== Runtime import test ===")
import subprocess, sys
result = subprocess.run(
    [sys.executable, '-X', 'utf8', '-c',
     'from ua_vacuum_manifold import ua_dpm_total_density, DPM_DENSITY_RATIO; '
     'from dpm_vacuum_manifold import ELEMENT, E_CRACK, DPM_DENSITY_RATIO as dr; '
     'from scm_vacuum_manifold import derive_from_quantum_chain, RHO_VAC_SCM; '
     'rho,_=derive_from_quantum_chain(); '
     'print(f"SCm RHO_VAC_SCM={RHO_VAC_SCM:.4e}"); '
     'print(f"UA total={ua_dpm_total_density(0.0):.4e}"); '
     'print(f"UA DPM_DENSITY_RATIO={DPM_DENSITY_RATIO}"); '
     'print(f"DPM E_CRACK={E_CRACK:.4e}"); '
     'print(f"DPM ratio={dr:.4f}"); '
     'print("All 3 manifolds importable OK")'],
    capture_output=True, text=True
)
print(result.stdout)
if result.returncode != 0:
    print("STDERR:", result.stderr[:400])
