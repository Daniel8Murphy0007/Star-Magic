"""
Scan CP1-CP4 for ALL density constant assignments, both:
  (a) module-level  — these SHOULD use _RHO_VAC_SCM / Quantum Chain
  (b) class-internal (self.x = ..., local x = ...) — these may be hardcoded
Report every unique value found and whether it traces back to QC or is hardcoded.
"""
import re, os

CP_FILES = {
    'CP1': 'CondensedPhysics.py',
    'CP2': 'CondensedPhysics2.py',
    'CP3': 'CondensedPhysics3.py',
    'CP4': 'CondensedPhysics4.py',
}

# Patterns for vacuum density constants — any float that looks like a density constant
# We care about: rho_vac, RHO_VAC, rho_SCm, rho_UA, SCm_DENSITY, UA_DENSITY, etc.
DENSITY_VAR_PAT = re.compile(
    r'(?:rho_vac|RHO_VAC|rho_SCm|rho_UA|SCm_DENSITY|UA_DENSITY|SCm_INF|UA_VAL|'
    r'RHO_SCM|RHO_UA|vacuum_density|vac_density)',
    re.IGNORECASE
)

# Hardcoded float pattern
HARDCODED_PAT = re.compile(r'=\s*([\d.]+e[+-]?\d+|\d+\.\d+(?:e[+-]?\d+)?)')

# Reference to module-level QC variable
QC_REF_PAT = re.compile(r'=\s*(_RHO_VAC_SCM|_RHO_VAC_UA|_UA_MANIFOLD_LOADED|_DPM_MANIFOLD_LOADED'
                         r'|derive_from_quantum_chain|_ua_layer_density|_ua_dpm_total|RHO_VAC_SCM'
                         r'|RHO_VAC_UA|scm_vacuum_manifold)')

for label, fname in CP_FILES.items():
    print(f"\n{'='*65}")
    print(f"{label}: {fname}")
    print(f"{'='*65}")

    with open(fname, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    hardcoded_class = []   # inside class __init__ or methods
    hardcoded_module = []  # module-level
    qc_refs = []           # points to QC value
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        if not DENSITY_VAR_PAT.search(stripped):
            continue
        if '=' not in stripped:
            continue
        
        lineno = i + 1
        is_class = line.startswith('        ') or line.startswith('    ')  # indented = class/method
        is_module = not line.startswith(' ')

        if QC_REF_PAT.search(stripped):
            qc_refs.append((lineno, stripped[:100]))
        elif HARDCODED_PAT.search(stripped):
            val_m = HARDCODED_PAT.search(stripped)
            val = val_m.group(1) if val_m else '?'
            if is_module:
                hardcoded_module.append((lineno, val, stripped[:100]))
            else:
                hardcoded_class.append((lineno, val, stripped[:100]))

    print(f"  QC-derived references: {len(qc_refs)}")
    for ln, s in qc_refs[:10]:
        print(f"    L{ln}: {s}")

    print(f"\n  HARDCODED module-level: {len(hardcoded_module)}")
    for ln, val, s in hardcoded_module[:10]:
        print(f"    L{ln} val={val}: {s}")

    print(f"\n  HARDCODED class-internal: {len(hardcoded_class)}")
    for ln, val, s in hardcoded_class[:20]:
        print(f"    L{ln} val={val}: {s}")
    if len(hardcoded_class) > 20:
        # Group by unique values
        vals = {}
        for ln, val, s in hardcoded_class:
            vals.setdefault(val, []).append(ln)
        print(f"  ... ({len(hardcoded_class)} total class-internal hardcoded)")
        print(f"  Unique values found:")
        for v, lns in sorted(vals.items(), key=lambda x: len(x[1]), reverse=True):
            print(f"    {v} ({len(lns)} occurrences, first L{lns[0]})")
