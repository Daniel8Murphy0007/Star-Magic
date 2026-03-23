"""
Extract system-specific parameters from each module's constructor.
Outputs a summary of all unique physics parameters per system.
"""
import re

FNAME = 'grok_share_bdfb3a05b06.txt'

with open(FNAME, encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

# Map of class name → constructor content
sep_lines = []
for i, line in enumerate(lines):
    if line.strip() == 'cpp':
        sep_lines.append(i)

# Pair them
pairs = []
for k in range(0, len(sep_lines) - 1, 2):
    h_sep = sep_lines[k]
    c_sep = sep_lines[k + 1]
    end_idx = sep_lines[k + 2] if k + 2 < len(sep_lines) else len(lines)
    pairs.append((h_sep, c_sep, end_idx))

def get_classname(lines_slice):
    for line in lines_slice[:10]:
        m = re.match(r'class\s+(\w+)', line.strip())
        if m:
            return m.group(1)
    return None

# Dedup map
class_registry = {}
for idx, (h_sep, c_sep, end_idx) in enumerate(pairs):
    h_lines = lines[h_sep+1:c_sep]
    c_lines = lines[c_sep+1:end_idx]
    classname = get_classname(h_lines) or get_classname(c_lines)
    if classname:
        size = len(c_lines)
        if classname not in class_registry or size > class_registry[classname][1]:
            class_registry[classname] = (idx, size, c_sep, end_idx)

PARAM_KEYS = ['M', 'r', 'L_X', 'B0', 'omega0', 'omega_LENR', 'E_cm_astro', 
              'F_neutrino', 'f_TRZ', 't', 'rho_gas', 'z', 'k_rel', 'k_LENR',
              'x2', 'T', 'n_H', 'omega_H', 'k_DE', 'k_act', 'F0']
SYSTEM_PARAMS = ['M', 'r', 'L_X', 'B0', 'omega0', 't', 'rho_gas', 'E_cm_astro', 'z']

results = {}
for classname, (idx, size, c_sep, end_idx) in sorted(class_registry.items(), key=lambda x: x[1][0]):
    c_lines = lines[c_sep+1:end_idx]
    params = {}
    for line in c_lines:
        m = re.match(r'\s*variables\["(\w+)"\]\s*=\s*\{([^}]+)\}', line)
        if m:
            key = m.group(1)
            val = m.group(2).split(',')[0].strip()
            params[key] = val
    results[classname] = params

print("=== SYSTEM PARAMETER TABLE ===\n")
print(f"{'Class':<35} {'M(kg)':<14} {'r(m)':<14} {'omega0':<12} {'L_X(W)':<14} {'E_cm_astro':<14} {'t(s)':<14}")
print("-" * 120)
for classname, params in results.items():
    M = params.get('M', 'N/A')
    r = params.get('r', 'N/A')
    om0 = params.get('omega0', 'N/A')
    LX = params.get('L_X', 'N/A')
    Ecm = params.get('E_cm_astro', 'N/A')
    t = params.get('t', 'N/A')
    print(f"{classname:<35} {M:<14} {r:<14} {om0:<12} {LX:<14} {Ecm:<14} {t:<14}")

print("\n\n=== ALL SYSTEM-SPECIFIC VARIABLES ===\n")
# Universal params that appear in every module (base constants)
UNIVERSAL = {'G', 'c', 'hbar', 'q', 'pi', 'm_e', 'mu_B', 'g_Lande', 'k_B', 'mu0',
             'rho_vac_UA', 'DPM_momentum', 'DPM_gravity', 'DPM_stability',
             'k_LENR', 'omega_LENR', 'k_act', 'phi', 'k_DE', 'k_neutron', 'sigma_n',
             'k_rel', 'F_neutrino', 'F0', 'x2', 'beta_i', 'V_infl_UA', 'rho_vac_A',
             'a_universal', 'lambda_i', 'rho_vac_SCm', 'omega_s', 'f_TRZ', 't_scale',
             'V', 'B0', 'theta', 'omega0'}

for classname, params in results.items():
    unique_params = {k: v for k, v in params.items() if k not in UNIVERSAL}
    print(f"\n[{classname}]")
    # Always show M, r, system-level params  
    for pk in SYSTEM_PARAMS:
        if pk in params:
            print(f"  {pk} = {params[pk]}")
    if unique_params:
        print(f"  UNIQUE params: {unique_params}")
    # E_cm value (varies per system)
    if 'E_cm' in params:
        print(f"  E_cm = {params['E_cm']} J")
    if 'omega0' in params:
        print(f"  omega0 = {params['omega0']}")
    if 't' in params:
        print(f"  t_default = {params['t']} s")
    if 'L_X' in params:
        print(f"  L_X = {params['L_X']} W")
