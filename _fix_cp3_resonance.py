"""Fix CP3 CoAnQiModularResonanceMUGECalculator to implement all 12 PAPER_371 terms."""
import re

filepath = r'c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics3.py'

with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
    content = f.read()

# Find the old compute method - search for unique marker
marker = "dataset.get('M', M_BH_SGR)"
idx = content.find(marker)
if idx < 0:
    print("ERROR: Could not find marker in file")
    exit(1)

print(f"Found marker at position {idx}")

# Find the start of this compute method (go back to find "def compute")
method_start = content.rfind("    def compute(self, dataset: dict) -> dict:", 0, idx)
if method_start < 0:
    print("ERROR: Could not find method start")
    exit(1)

# Find the end - next class definition or next method at same indent level
# Look for the simulation_set closing brace + next class
end_marker = "'simulation_set': {'g_res_vs_omega': 'omega from 1e-6 to 1e3 rad/s'},"
end_idx = content.find(end_marker, method_start)
if end_idx < 0:
    print("ERROR: Could not find end marker")
    exit(1)

# Find the closing of the return block - go past the end_marker to find the closing braces
rest = content[end_idx:]
# Find the "}" closing the return dict, then the class boundary
closing = rest.find("\n        }")
if closing < 0:
    print("ERROR: Could not find closing brace")
    exit(1)

method_end = end_idx + closing + len("\n        }")

print(f"Method from {method_start} to {method_end}")
print(f"Old method length: {method_end - method_start}")
print("Old method first 100 chars:", repr(content[method_start:method_start+100]))

new_method = '''    def compute(self, dataset: dict) -> dict:
        # System parameters
        I       = dataset.get('I', 1e21)
        A_area  = dataset.get('A', 3.142e8)
        w1      = dataset.get('omega1', 1e-3)
        w2      = dataset.get('omega2', 0.0)
        Vsys    = dataset.get('Vsys', 4.189e12)
        vexp    = dataset.get('vexp', 1e3)
        t       = dataset.get('t', 3.799e10)
        ffluid  = dataset.get('ffluid', 1.269e-14)
        kappa   = dataset.get('kappa', 0.0005 / 86400)
        r       = dataset.get('r', 1e4)
        f_worm  = dataset.get('f_worm', 1e-10)
        b_worm  = dataset.get('b_worm', 1e6)
        t_n     = dataset.get('t_n', 0.0)

        # Resonance constants (PAPER_371 section 3)
        fDPM    = 1e12;  fTHz    = 1e12
        Evac_nb = 7.09e-36;  Evac_ISM = 7.09e-37;  Delta_E = 6.381e-36
        Fsuper  = 6.287e-19; UA_SCM  = 10.0
        omega_i = 1e-8;  k4r     = 1.0;  freact  = 1e10
        fq      = 1.445e-17;  fAe     = 1.576e-35
        fosc    = 4.57e14
        fTRZ    = 0.1;  H_z     = 2.270e-18;  c_res   = 3e8

        # All 12 terms per PAPER_371 section 2.1-2.12
        FDPM    = I * A_area * (w1 - w2)
        aDPM    = FDPM * fDPM * Evac_nb * c_res * Vsys
        aTHz    = fTHz * Evac_nb * vexp * aDPM / Evac_ISM / c_res
        avd     = Delta_E * vexp**2 * aDPM / Evac_nb / c_res**2
        asf     = Fsuper * fTHz * aDPM / Evac_nb / c_res
        aar     = UA_SCM * omega_i * fTHz * aDPM * (1.0 + fTRZ)
        Ereact  = 1e46 * math.exp(-kappa * t)
        Ug4i    = k4r * Ereact * freact * aDPM / Evac_nb * c_res
        aqf     = fq * Evac_nb * aDPM / Evac_ISM / c_res
        aAf     = fAe * Evac_nb * aDPM / Evac_ISM / c_res
        afl     = ffluid * Evac_nb * Vsys / Evac_ISM / c_res
        Osc     = fosc * math.cos(2.0 * math.pi * fosc * t)
        aexp    = 2.0 * math.pi * H_z * t * Evac_nb * aDPM / Evac_ISM / c_res
        a_worm  = f_worm * Evac_nb / (b_worm**2 + r**2) if r > 0 else f_worm * Evac_nb

        g_res = (aDPM + aTHz + avd + asf + aar + Ug4i + aqf + aAf
                 + afl + Osc + aexp + fTRZ + a_worm)
        eqs = {
            'g_resonance_MUGE': f'{g_res:.4e} m/s^2',
            'aDPM': f'{aDPM:.4e}', 'aTHz': f'{aTHz:.4e}',
            'avac_diff': f'{avd:.4e}', 'asuper_freq': f'{asf:.4e}',
            'aaether_res': f'{aar:.4e}', 'Ug4i': f'{Ug4i:.4e}',
            'aquantum_freq': f'{aqf:.4e}', 'aAether_freq': f'{aAf:.4e}',
            'afluid_freq': f'{afl:.4e}', 'Osc_term': f'{Osc:.4e}',
            'aexp_freq': f'{aexp:.4e}', 'fTRZ': f'{fTRZ}',
            'a_wormhole_13th': f'{a_worm:.4e}',
            'CoAnQi_plugin': 'ResonanceMUGEPlugin (13-term)',
        }
        return {
            'primary_equations': eqs,
            'available_equations': [
                'g_res = aDPM + aTHz + avac_diff + asuper_freq + aaether_res + Ug4i'
                ' + aquantum_freq + aAether_freq + afluid_freq + Osc + aexp + fTRZ + a_worm',
                'All 12 core terms per PAPER_371 + Morris-Thorne wormhole (13th)',
            ],
            'simulation_set': {
                'SGR1745': {'I': 1e21, 'A': 3.142e8, 'omega1': 1e-3, 'Vsys': 4.189e12},
                'SagA':    {'I': 1e23, 'A': 2.813e30, 'omega1': 1e-5, 'Vsys': 3.552e45},
                'Pillars': {'I': 1e21, 'A': 2.813e32, 'omega1': 1e-3, 'Vsys': 3.552e48},
            },
        }'''

new_content = content[:method_start] + new_method + content[method_end:]

with open(filepath, 'w', encoding='utf-8') as f:
    f.write(new_content)

print("SUCCESS: CP3 CoAnQiModularResonanceMUGECalculator updated with all 12 PAPER_371 terms")
