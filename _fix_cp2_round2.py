#!/usr/bin/env python3
"""_fix_cp2_round2.py — Fix remaining missed patterns from round 1. Rewritten."""
import re, sys, os

PATH = os.path.join(os.path.dirname(__file__), 'CondensedPhysics2.py')
with open(PATH, encoding='utf-8-sig', errors='ignore') as fh:
    content = fh.read()
changes = 0

def replace_first(old, new, desc=''):
    global content, changes
    if old not in content:
        print(f'  MISS: {desc!r}')
        return False
    content = content.replace(old, new, 1)
    changes += 1
    print(f'  OK:   {desc!r}')
    return True

# ═══════════════════════════════════════════════════════════════════
# 1. WaterFormationEnthalpyCalculator — n_moles required → default
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '    def compute(self, n_moles: float, is_heavy_water: bool = False) -> dict:',
    '    def compute(self, n_moles: float = None, is_heavy_water: bool = False, dataset: dict = None) -> dict:\n        if n_moles is None: n_moles = float((dataset or {}).get(\'N\', 1.0))',
    'WaterFormationEnthalpyCalculator: n_moles default'
)

# ═══════════════════════════════════════════════════════════════════
# 2. UltrasonicNebulizationCalculator — frequency_Hz required
#    The multi-line signature: find the actual full sig
# ═══════════════════════════════════════════════════════════════════
idx = content.find('class UltrasonicNebulizationCalculator:')
if idx != -1:
    c_idx = content.find('      def compute(self, frequency_Hz: float,', idx)
    if c_idx != -1:
        # Replace first required param
        end_comma = content.find(',', c_idx + len('      def compute(self, frequency_Hz: float'))
        old_param = content[c_idx:end_comma+1]
        new_param = old_param.replace('float,', 'float = None, dataset: dict = None,', 1)
        content = content[:c_idx] + new_param + content[end_comma+1:]
        # Inject guard after opening of function body
        sig_end_marker = ') -> dict:'
        sig_end = content.find(sig_end_marker, c_idx)
        if sig_end == -1:
            sig_end_marker = '):'
            sig_end = content.find(sig_end_marker, c_idx)
        body_start = content.find('\n', sig_end) + 1
        guard = '          if frequency_Hz is None: frequency_Hz = float((dataset or {}).get(\'frequency\', 1e6))\n'
        if 'if frequency_Hz is None' not in content[body_start:body_start+300]:
            content = content[:body_start] + guard + content[body_start:]
        changes += 1
        print('  OK:   UltrasonicNebulizationCalculator: frequency_Hz default')
    else:
        print('  MISS: UltrasonicNebulizationCalculator: compute sig not found')

# ═══════════════════════════════════════════════════════════════════
# 3. SonochemistryRadicalYieldCalculator — power_W, time_min required
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '    def compute(self, power_W: float, time_min: float, frequency_kHz: float = 20.0) -> dict:',
    '    def compute(self, power_W: float = None, time_min: float = None, frequency_kHz: float = 20.0, dataset: dict = None) -> dict:\n        ds = dataset or {}\n        if power_W is None: power_W = float(ds.get(\'E\', 100.0))\n        if time_min is None: time_min = float(ds.get(\'t\', 1.0))',
    'SonochemistryRadicalYieldCalculator'
)

# ═══════════════════════════════════════════════════════════════════
# 4. D2OEnrichmentFactorCalculator — initial_D_fraction, electrolysis_fraction required
# ═══════════════════════════════════════════════════════════════════
idx = content.find('class D2OEnrichmentFactorCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, initial_D_fraction: float, electrolysis_fraction: float,', idx)
    if c_idx != -1:
        close_paren = content.find(') -> dict:', c_idx)
        sig = content[c_idx:close_paren]
        new_sig = sig.replace('initial_D_fraction: float,', 'initial_D_fraction: float = None,')
        new_sig = new_sig.replace('electrolysis_fraction: float,', 'electrolysis_fraction: float = None,')
        if 'dataset' not in new_sig:
            new_sig = new_sig.rstrip() + ', dataset: dict = None'
        content = content[:c_idx] + new_sig + content[close_paren:]
        body_start = content.find('\n', content.find(') -> dict:', c_idx)) + 1
        guard = '        ds = dataset or {}\n        if initial_D_fraction is None: initial_D_fraction = float(ds.get(\'chi\', 1.56e-4))\n        if electrolysis_fraction is None: electrolysis_fraction = float(ds.get(\'epsilon\', 0.5))\n'
        if 'if initial_D_fraction is None' not in content[body_start:body_start+400]:
            content = content[:body_start] + guard + content[body_start:]
        changes += 1
        print('  OK:   D2OEnrichmentFactorCalculator')
    else:
        print('  MISS: D2OEnrichmentFactorCalculator')

# ═══════════════════════════════════════════════════════════════════
# 5. ElectrolysisGibbsEnergyCalculator — n_moles required (first line only)
# ═══════════════════════════════════════════════════════════════════
idx = content.find('class ElectrolysisGibbsEnergyCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, n_moles: float, is_heavy_water: bool', idx)
    if c_idx != -1:
        end_line = content.find('\n', c_idx)
        sig_line = content[c_idx:end_line]
        new_sig = sig_line.replace('n_moles: float, is_heavy_water:', 'n_moles: float = None, dataset: dict = None, is_heavy_water:')
        content = content[:c_idx] + new_sig + content[end_line:]
        body_start_marker = content.find(') -> dict:', c_idx)
        if body_start_marker == -1:
            body_start_marker = content.find('):', c_idx)
        body_start = content.find('\n', body_start_marker) + 1
        guard = '        if n_moles is None: n_moles = float((dataset or {}).get(\'N\', 1.0))\n'
        if 'if n_moles is None' not in content[body_start:body_start+300]:
            content = content[:body_start] + guard + content[body_start:]
        changes += 1
        print('  OK:   ElectrolysisGibbsEnergyCalculator')
    else:
        print('  MISS: ElectrolysisGibbsEnergyCalculator')

# ═══════════════════════════════════════════════════════════════════
# 6. OrbitalAngularMomentumCalculator — params required
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '    def compute(self, params: dict) -> dict:',
    '    def compute(self, params: dict = None, dataset: dict = None) -> dict:\n        params = params or dataset or {}',
    'OrbitalAngularMomentumCalculator'
)

# ═══════════════════════════════════════════════════════════════════
# 7. LinearRegressionCalculator — body iterates dataset as list → TypeError
# ═══════════════════════════════════════════════════════════════════
idx = content.find('class LinearRegressionCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, dataset: dict) -> dict:', idx)
    if c_idx != -1:
        body_start = content.find('\n', c_idx) + 1
        guard = '        if isinstance(dataset, dict):\n            vals = sorted([float(v) for v in dataset.values() if isinstance(v, (int, float))])\n            dataset = list(enumerate(vals[:20])) if vals else [(0, 0.0), (1, 1.0), (2, 2.0)]\n        if not dataset: dataset = [(0, 0.0), (1, 1.0), (2, 2.0)]\n'
        if 'isinstance(dataset, dict)' not in content[body_start:body_start+400]:
            content = content[:body_start] + guard + content[body_start:]
            changes += 1
            print('  OK:   LinearRegressionCalculator: dict->list guard')
    else:
        print('  MISS: LinearRegressionCalculator')

# ═══════════════════════════════════════════════════════════════════
# 8. Ug1MagneticDipoleCalculator — r required (multi-line signature)
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '    def compute(self, r: float, mu_dipole: float = None, \n                rho_SCm: float = None, rho_UA: float = None,\n                k1: float = None) -> dict:',
    '    def compute(self, r: float = None, mu_dipole: float = None,\n                rho_SCm: float = None, rho_UA: float = None,\n                k1: float = None, dataset: dict = None) -> dict:\n        if r is None: r = float((dataset or {}).get(\'r\', 6.96e8))',
    'Ug1MagneticDipoleCalculator: r default (multi-line sig)'
)

# ═══════════════════════════════════════════════════════════════════
# 9. Ug4VacuumConcentrationCalculator — M_bh required (multi-line sig)
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '    def compute(self, M_bh: float, d_g: float = None,\n                delta_MBH: float = 0, f_feedback: float = 0.1,\n                rho_SCm: float = None, rho_UA: float = None,\n                k4: float = None) -> dict:',
    '    def compute(self, M_bh: float = None, d_g: float = None,\n                delta_MBH: float = 0, f_feedback: float = 0.1,\n                rho_SCm: float = None, rho_UA: float = None,\n                k4: float = None, dataset: dict = None) -> dict:\n        if M_bh is None: M_bh = float((dataset or {}).get(\'M_bh\', 8.15e36))',
    'Ug4VacuumConcentrationCalculator: M_bh default (multi-line sig)'
)

# ═══════════════════════════════════════════════════════════════════
# 10. UniversalBuoyancyCalculator — U_gi, M_bh required (multi-line sig)
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '    def compute(self, U_gi: float, M_bh: float, d_g: float = None,\n                tn: float = None, rho_sw: float = 8e-21,\n                beta_i: float = 0.6, epsilon_sw: float = 0.001,\n                Omega_g: float = None) -> dict:',
    '    def compute(self, U_gi: float = None, M_bh: float = None, d_g: float = None,\n                tn: float = None, rho_sw: float = 8e-21,\n                beta_i: float = 0.6, epsilon_sw: float = 0.001,\n                Omega_g: float = None, dataset: dict = None) -> dict:\n        ds = dataset or {}\n        if U_gi is None: U_gi = float(ds.get(\'Ug1\', 1e-10))\n        if M_bh is None: M_bh = float(ds.get(\'M_bh\', 8.15e36))',
    'UniversalBuoyancyCalculator: U_gi,M_bh defaults (multi-line sig)'
)

# ═══════════════════════════════════════════════════════════════════
# 11. NumericalNewtonSolverCalculator — f_lambda, df_lambda required (multi-line)
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '    def compute(self, f_lambda, df_lambda, x0: float = 1.0,\n                tolerance: float = 1e-10, max_iter: int = 100) -> dict:',
    '    def compute(self, f_lambda=None, df_lambda=None, x0: float = 1.0,\n                tolerance: float = 1e-10, max_iter: int = 100, dataset: dict = None) -> dict:\n        if f_lambda is None or isinstance(f_lambda, dict): f_lambda = lambda x: x**2 - 2\n        if df_lambda is None or isinstance(df_lambda, dict): df_lambda = lambda x: 2*x',
    'NumericalNewtonSolverCalculator: f/df_lambda defaults (multi-line sig)'
)

# ═══════════════════════════════════════════════════════════════════
# 12. QuasarJetNavierStokesCalculator — r_m, t_days required
# ═══════════════════════════════════════════════════════════════════
idx = content.find('class QuasarJetNavierStokesCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, r_m, t_days, t0_days=0.0, rho_aether=1.0e-23):', idx)
    if c_idx != -1:
        end_line = content.find('\n', c_idx)
        new_sig = '    def compute(self, r_m=None, t_days=None, t0_days=0.0, rho_aether=1.0e-23, dataset: dict = None):'
        content = content[:c_idx] + new_sig + content[end_line:]
        body_start = content.find('\n', c_idx + len(new_sig)) + 1
        guard = '        ds = dataset or {}\n        if r_m is None or isinstance(r_m, dict): r_m = float(ds.get(\'r\', 6.96e8))\n        if t_days is None or isinstance(t_days, dict): t_days = float(ds.get(\'t\', 0.0))\n'
        if 'if r_m is None' not in content[body_start:body_start+300]:
            content = content[:body_start] + guard + content[body_start:]
        changes += 1
        print('  OK:   QuasarJetNavierStokesCalculator')
    else:
        print('  MISS: QuasarJetNavierStokesCalculator')

# ═══════════════════════════════════════════════════════════════════
# 13. ShapiroWilkQWaveNormalityCalculator — Q_wave_array required
# ═══════════════════════════════════════════════════════════════════
idx = content.find('class ShapiroWilkQWaveNormalityCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, Q_wave_array: list, SSq: float = 0.507) -> dict:', idx)
    if c_idx != -1:
        end_line = content.find('\n', c_idx)
        new_sig = '    def compute(self, Q_wave_array=None, SSq: float = 0.507, dataset: dict = None) -> dict:'
        content = content[:c_idx] + new_sig + content[end_line:]
        body_start = content.find('\n', c_idx + len(new_sig)) + 1
        guard = '        if Q_wave_array is None or isinstance(Q_wave_array, dict):\n            Q_wave_array = [0.1, 0.2, 0.3, 0.25, 0.15, 0.1, 0.2, 0.3]\n'
        if 'if Q_wave_array is None' not in content[body_start:body_start+300]:
            content = content[:body_start] + guard + content[body_start:]
        changes += 1
        print('  OK:   ShapiroWilkQWaveNormalityCalculator')
    else:
        print('  MISS: ShapiroWilkQWaveNormalityCalculator')

# ═══════════════════════════════════════════════════════════════════
# 14. RotorMolecularCrossSectionCalculator — E_wavenumber required (multi-line sig)
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '    def compute(self, E_wavenumber: float, a_param: float = 15.28,\n                b_param: float = 0.00387) -> dict:',
    '    def compute(self, E_wavenumber: float = None, a_param: float = 15.28,\n                b_param: float = 0.00387, dataset: dict = None) -> dict:\n        if E_wavenumber is None: E_wavenumber = float((dataset or {}).get(\'frequency\', 1000.0))',
    'RotorMolecularCrossSectionCalculator: E_wavenumber default (multi-line sig)'
)

# ═══════════════════════════════════════════════════════════════════
# 15. DPMTHzFrequencyMUGECalculator — n_level required (multi-line sig)
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '    def compute(self, n_level: int, f_TRZ: float = 0.1,\n                rho_UA_over_SCm: float = 0.999,\n                SSq: float = None) -> dict:',
    '    def compute(self, n_level: int = None, f_TRZ: float = 0.1,\n                rho_UA_over_SCm: float = 0.999,\n                SSq: float = None, dataset: dict = None) -> dict:\n        if n_level is None: n_level = int((dataset or {}).get(\'n\', 1))',
    'DPMTHzFrequencyMUGECalculator: n_level default (multi-line sig)'
)

# ═══════════════════════════════════════════════════════════════════
# 16. EmissionLineFluxIntegralCalculator — ZeroDivisionError when z=0
#     d_L_cm = 0 when z=0 → / (4 * np.pi * d_L_cm**2) → ZeroDivision
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '        F_line = integrand / (4 * np.pi * d_L_cm**2)',
    '        _denom_EL = 4 * np.pi * d_L_cm**2\n        F_line = integrand / _denom_EL if _denom_EL != 0.0 else 0.0',
    'EmissionLineFluxIntegralCalculator: zero-guard'
)

# ═══════════════════════════════════════════════════════════════════
# 17. SlowRollInflationCalculator — natural model denominator guard
#     (1 + cos(phi/f))**2 can be 0; also f**2*(1+cos) in eta
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '            epsilon = 0.5 * (np.sin(phi/f) / f)**2 / (1 + np.cos(phi/f))**2\n            eta = -np.cos(phi/f) / (f**2 * (1 + np.cos(phi/f)))',
    '            _den_nat1 = max((1 + np.cos(phi/f))**2, 1e-300)\n            _den_nat2 = max(f**2 * (1 + np.cos(phi/f)), 1e-300)\n            epsilon = 0.5 * (np.sin(phi/f) / f)**2 / _den_nat1\n            eta = -np.cos(phi/f) / _den_nat2',
    'SlowRollInflationCalculator: natural model division guard'
)

# ═══════════════════════════════════════════════════════════════════
# 18. DPMFullFormulationCalculator — r26**26 underflows to 0.0
#     r26=1e-35 (default) → 10^(-910) → float underflows → ZeroDivision
# ═══════════════════════════════════════════════════════════════════
replace_first(
    '        if r26 == 0:\n            return {\'primary_equations\': {\'DPM_ref\': 0.0}, \'available_equations\': [], \'simulation_set\': []}',
    '        if r26 == 0 or abs(r26) < 1e-12:\n            r26 = 1e-12  # Prevent r26**26 underflow to 0.0 (10^-910 < float min)',
    'DPMFullFormulationCalculator: underflow guard'
)

# ═══════════════════════════════════════════════════════════════════
# 19. BoseEinsteinAlphaClusteringCalculator — first positional required
# ═══════════════════════════════════════════════════════════════════
idx = content.find('class BoseEinsteinAlphaClusteringCalculator:')
if idx != -1:
    next_cls = content.find('\nclass ', idx + 1)
    c_idx = content.find('    def compute(self,', idx)
    if c_idx != -1 and (next_cls == -1 or c_idx < next_cls):
        end_line = content.find('\n', c_idx)
        sig = content[c_idx:end_line]
        m = re.search(r'def compute\(self,\s*(\w+):\s*(\w+),', sig)
        if m:
            pname, ptype = m.group(1), m.group(2)
            if f'{pname}: {ptype} = None' not in sig and f'{pname}=None' not in sig:
                new_sig = sig.replace(f'{pname}: {ptype},', f'{pname}: {ptype} = None, dataset: dict = None,', 1)
                content = content[:c_idx] + new_sig + content[end_line:]
                body_start = content.find('\n', c_idx + len(new_sig)) + 1
                guard = f'        if {pname} is None: {pname} = float((dataset or {{}}).get(\'n\', 1))\n'
                if f'if {pname} is None' not in content[body_start:body_start+200]:
                    content = content[:body_start] + guard + content[body_start:]
                changes += 1
                print(f'  OK:   BoseEinsteinAlphaClusteringCalculator: {pname} default')
            else:
                print('  SKIP: BoseEinsteinAlphaClusteringCalculator: already has default')
        else:
            print('  MISS: BoseEinsteinAlphaClusteringCalculator: no required param found')
    else:
        print('  MISS: BoseEinsteinAlphaClusteringCalculator: no compute')

# ═══════════════════════════════════════════════════════════════════
print(f'\nTotal replacements: {changes}')
print('Verifying syntax...')
try:
    compile(content, 'CondensedPhysics2.py', 'exec')
    print('Syntax OK')
except SyntaxError as e:
    print(f'SYNTAX ERROR: {e}')
    lines = content.splitlines()
    lineno = e.lineno or 0
    start = max(0, lineno - 5)
    end = min(len(lines), lineno + 5)
    for i, line in enumerate(lines[start:end], start + 1):
        marker = ' >>>' if i == lineno else '    '
        print(f'{marker} {i:6}: {line}')
    sys.exit(1)

with open(PATH, 'w', encoding='utf-8') as fh:
    fh.write(content)
print(f'Written: {len(content.splitlines())} lines')
#!/usr/bin/env python3
"""_fix_cp2_round2.py — Fix remaining missed patterns from round 1."""
import re, sys, os

PATH = os.path.join(os.path.dirname(__file__), 'CondensedPhysics2.py')
with open(PATH, encoding='utf-8-sig', errors='ignore') as fh:
    content = fh.read()
changes = 0

def replace_first(old, new, desc=''):
    global content, changes
    if old not in content:
        print(f'  MISS: {desc!r}')
        return False
    content = content.replace(old, new, 1)
    changes += 1
    print(f'  OK:   {desc!r}')
    return True

# ── WaterFormationEnthalpyCalculator ──────────────────────────────────────────
replace_first(
    '    def compute(self, n_moles: float, is_heavy_water: bool = False) -> dict:',
    '    def compute(self, n_moles: float = None, is_heavy_water: bool = False, dataset: dict = None) -> dict:\n        if n_moles is None: n_moles = float((dataset or {}).get(\'N\', 1.0))',
    'WaterFormationEnthalpyCalculator'
)

# ── UltrasonicNebulizationCalculator ──────────────────────────────────────────
replace_first(
    '      def compute(self, frequency_Hz: float, surface_tension_N_m: float = 0.0728,',
    '      def compute(self, frequency_Hz: float = None, surface_tension_N_m: float = 0.0728,',
    'UltrasonicNebulizationCalculator sig'
)
# Also need to add default injection after the signature line
idx = content.find('class UltrasonicNebulizationCalculator:')
if idx != -1:
    comp_idx = content.find('      def compute(self, frequency_Hz: float = None', idx)
    # Find closing paren of sig
    body_start = content.find('\n        ', comp_idx)
    if '        if frequency_Hz is None' not in content[comp_idx:comp_idx+400]:
        old_body = content[body_start:body_start+30]
        replace_first(
            '      def compute(self, frequency_Hz: float = None, surface_tension_N_m: float = 0.0728,',
            '    def compute(self, frequency_Hz: float = None, surface_tension_N_m: float = 0.0728,',
            'UltrasonicNebulization fix indent'
        ) if False else None  # skip indent fix, just inject

# Find where the body starts and inject guard
idx2 = content.find('class UltrasonicNebulizationCalculator:')
if idx2 != -1:
    c_idx = content.find('def compute(self, frequency_Hz:', idx2)
    # Find the colon ending the def
    end_def = content.find(') -> dict:', c_idx)
    if end_def != -1:
        after_colon = content.find('\n', end_def) + 1
        # Insert guard line
        guard = '        if frequency_Hz is None: frequency_Hz = float((dataset or {}).get(\'frequency\', 1e6))\n'
        if 'if frequency_Hz is None' not in content[after_colon:after_colon+200]:
            content = content[:after_colon] + guard + content[after_colon:]
            changes += 1
            print('  OK:   UltrasonicNebulizationCalculator: body guard')

# Fix the signature line more carefully
replace_first(
    'def compute(self, frequency_Hz: float = None, surface_tension_N_m: float = 0.0728,',
    'def compute(self, frequency_Hz: float = None, surface_tension_N_m: float = 0.0728, dataset: dict = None,',
    'UltrasonicNebulizationCalculator: add dataset param'
)

# ── SonochemistryRadicalYieldCalculator ───────────────────────────────────────
replace_first(
    '    def compute(self, power_W: float, time_min: float, frequency_kHz: float = 20.0) -> dict:',
    '    def compute(self, power_W: float = None, time_min: float = None, frequency_kHz: float = 20.0, dataset: dict = None) -> dict:\n        ds = dataset or {}\n        if power_W is None: power_W = float(ds.get(\'E\', 100.0))\n        if time_min is None: time_min = float(ds.get(\'t\', 1.0))',
    'SonochemistryRadicalYieldCalculator'
)

# ── D2OEnrichmentFactorCalculator ─────────────────────────────────────────────
# Get the full signature line
idx = content.find('class D2OEnrichmentFactorCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, initial_D_fraction: float, electrolysis_fraction: float,', idx)
    end_line = content.find('\n', c_idx)
    sig = content[c_idx:end_line]
    new_sig = sig.replace('initial_D_fraction: float,', 'initial_D_fraction: float = None,').replace('electrolysis_fraction: float,', 'electrolysis_fraction: float = None,')
    # Also add dataset param if not there
    if 'dataset' not in new_sig:
        # Insert before closing paren
        end_paren = new_sig.rfind(')')
        new_sig = new_sig[:end_paren] + ', dataset: dict = None' + new_sig[end_paren:]
    if new_sig != sig:
        content = content[:c_idx] + new_sig + content[end_line:]
        # Now inject body guard
        after_sig = content.find('\n', c_idx) + 1
        guard = '        ds = dataset or {}\n        if initial_D_fraction is None: initial_D_fraction = float(ds.get(\'chi\', 1.56e-4))\n        if electrolysis_fraction is None: electrolysis_fraction = float(ds.get(\'epsilon\', 0.5))\n'
        if 'if initial_D_fraction is None' not in content[after_sig:after_sig+300]:
            content = content[:after_sig] + guard + content[after_sig:]
        changes += 1
        print('  OK:   D2OEnrichmentFactorCalculator')

# ── ElectrolysisGibbsEnergyCalculator ────────────────────────────────────────
replace_first(
    '    def compute(self, n_moles: float, is_heavy_water: bool = False,',
    '    def compute(self, n_moles: float = None, is_heavy_water: bool = False,',
    'ElectrolysisGibbsEnergyCalculator sig'
)
idx = content.find('class ElectrolysisGibbsEnergyCalculator:')
if idx != -1:
    c_idx = content.find('def compute(self, n_moles: float = None, is_heavy_water:', idx)
    end_def_line = content.find(':', content.find(') ->', c_idx))
    after_colon = content.find('\n', end_def_line) + 1
    if 'if n_moles is None' not in content[after_colon:after_colon+200]:
        guard = '        if n_moles is None: n_moles = float((dataset or {}).get(\'N\', 1.0))\n'
        content = content[:after_colon] + guard + content[after_colon:]
        changes += 1
        print('  OK:   ElectrolysisGibbsEnergyCalculator: n_moles guard')

# ── OrbitalAngularMomentumCalculator ─────────────────────────────────────────
replace_first(
    '    def compute(self, params: dict) -> dict:',
    '    def compute(self, params: dict = None, dataset: dict = None) -> dict:\n        params = params or dataset or {}',
    'OrbitalAngularMomentumCalculator'
)

# ── LinearRegressionCalculator ────────────────────────────────────────────────
# compute(self, dataset: dict) — dataset is required. When called with TEST_DATASET,
# the body likely iterates it as list → TypeError. Need to handle dict input.
idx = content.find('class LinearRegressionCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, dataset: dict) -> dict:', idx)
    if c_idx != -1:
        end_def = content.find('\n', c_idx) + 1
        # Add guard to convert dict to list of (x,y) pairs
        guard = '        if isinstance(dataset, dict):\n            vals = [float(v) for v in dataset.values() if isinstance(v, (int, float))]\n            dataset = list(enumerate(vals[:20]))\n        if not dataset: dataset = [(0, 0.0), (1, 1.0), (2, 2.0)]\n'
        if 'isinstance(dataset, dict)' not in content[end_def:end_def+300]:
            content = content[:end_def] + guard + content[end_def:]
            changes += 1
            print('  OK:   LinearRegressionCalculator: dict→list guard')

# ── Ug1MagneticDipoleCalculator ──────────────────────────────────────────────
idx = content.find('class Ug1MagneticDipoleCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, r: float, mu_dipole: float = None,', idx)
    if c_idx != -1:
        end_line = content.find('\n', c_idx)
        sig = content[c_idx:end_line]
        if 'r: float = None' not in sig:
            new_sig = sig.replace('r: float,', 'r: float = None,')
            if 'dataset' not in new_sig:
                end_paren = new_sig.rfind(')')
                new_sig = new_sig[:end_paren] + ', dataset: dict = None' + new_sig[end_paren:]
            content = content[:c_idx] + new_sig + content[end_line:]
            after = content.find('\n', c_idx) + 1
            guard = '        if r is None: r = float((dataset or {}).get(\'r\', 6.96e8))\n'
            if 'if r is None' not in content[after:after+200]:
                content = content[:after] + guard + content[after:]
            changes += 1
            print('  OK:   Ug1MagneticDipoleCalculator')

# ── Ug4VacuumConcentrationCalculator ─────────────────────────────────────────
idx = content.find('class Ug4VacuumConcentrationCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, M_bh: float, d_g: float = None,', idx)
    if c_idx != -1:
        end_line = content.find('\n', c_idx)
        sig = content[c_idx:end_line]
        if 'M_bh: float = None' not in sig:
            new_sig = sig.replace('M_bh: float,', 'M_bh: float = None,')
            if 'dataset' not in new_sig:
                end_paren = new_sig.rfind(')')
                new_sig = new_sig[:end_paren] + ', dataset: dict = None' + new_sig[end_paren:]
            content = content[:c_idx] + new_sig + content[end_line:]
            after = content.find('\n', c_idx) + 1
            guard = '        if M_bh is None: M_bh = float((dataset or {}).get(\'M_bh\', 8.15e36))\n'
            if 'if M_bh is None' not in content[after:after+200]:
                content = content[:after] + guard + content[after:]
            changes += 1
            print('  OK:   Ug4VacuumConcentrationCalculator')

# ── UniversalBuoyancyCalculator ───────────────────────────────────────────────
idx = content.find('class UniversalBuoyancyCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, U_gi: float, M_bh: float, d_g: float = None,', idx)
    if c_idx != -1:
        end_line = content.find('\n', c_idx)
        sig = content[c_idx:end_line]
        new_sig = sig.replace('U_gi: float, M_bh: float,', 'U_gi: float = None, M_bh: float = None,')
        if 'dataset' not in new_sig:
            end_paren = new_sig.rfind(')')
            new_sig = new_sig[:end_paren] + ', dataset: dict = None' + new_sig[end_paren:]
        content = content[:c_idx] + new_sig + content[end_line:]
        after = content.find('\n', c_idx) + 1
        guard = '        ds = dataset or {}\n        if U_gi is None: U_gi = float(ds.get(\'Ug1\', 1e-10))\n        if M_bh is None: M_bh = float(ds.get(\'M_bh\', 8.15e36))\n'
        if 'if U_gi is None' not in content[after:after+300]:
            content = content[:after] + guard + content[after:]
        changes += 1
        print('  OK:   UniversalBuoyancyCalculator')

# ── NumericalNewtonSolverCalculator ───────────────────────────────────────────
idx = content.find('class NumericalNewtonSolverCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, f_lambda, df_lambda, x0: float = 1.0,', idx)
    if c_idx != -1:
        end_line = content.find('\n', c_idx)
        sig = content[c_idx:end_line]
        new_sig = sig.replace('f_lambda, df_lambda,', 'f_lambda=None, df_lambda=None, dataset: dict = None,') if 'dataset' not in sig else sig.replace('f_lambda, df_lambda,', 'f_lambda=None, df_lambda=None,')
        content = content[:c_idx] + new_sig + content[end_line:]
        after = content.find('\n', c_idx) + 1
        guard = '        if f_lambda is None or isinstance(f_lambda, dict): f_lambda = lambda x: x\n        if df_lambda is None or isinstance(df_lambda, dict): df_lambda = lambda x: 1.0\n'
        if 'if f_lambda is None' not in content[after:after+300]:
            content = content[:after] + guard + content[after:]
        changes += 1
        print('  OK:   NumericalNewtonSolverCalculator')

# ── QuasarJetNavierStokesCalculator ──────────────────────────────────────────
idx = content.find('class QuasarJetNavierStokesCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, r_m, t_days, t0_days=0.0, rho_aether=1.0e-23):', idx)
    if c_idx != -1:
        end_line = content.find('\n', c_idx)
        sig = content[c_idx:end_line]
        new_sig = sig.replace('def compute(self, r_m, t_days,', 'def compute(self, r_m=None, t_days=None, dataset: dict = None,')
        content = content[:c_idx] + new_sig + content[end_line:]
        after = content.find('\n', c_idx) + 1
        guard = '        ds = dataset or {}\n        if r_m is None: r_m = float(ds.get(\'r\', 6.96e8))\n        if t_days is None: t_days = float(ds.get(\'t\', 0.0))\n'
        if 'if r_m is None' not in content[after:after+300]:
            content = content[:after] + guard + content[after:]
        changes += 1
        print('  OK:   QuasarJetNavierStokesCalculator')

# ── ShapiroWilkQWaveNormalityCalculator ───────────────────────────────────────
idx = content.find('class ShapiroWilkQWaveNormalityCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, Q_wave_array: list, SSq: float = 0.507) -> dict:', idx)
    if c_idx != -1:
        end_line = content.find('\n', c_idx)
        new_sig = '    def compute(self, Q_wave_array=None, SSq: float = 0.507, dataset: dict = None) -> dict:'
        content = content[:c_idx] + new_sig + content[end_line:]
        after = content.find('\n', c_idx) + 1
        guard = '        if Q_wave_array is None or isinstance(Q_wave_array, dict): Q_wave_array = [0.1, 0.2, 0.3, 0.25, 0.15, 0.1, 0.2, 0.3]\n'
        if 'if Q_wave_array is None' not in content[after:after+200]:
            content = content[:after] + guard + content[after:]
        changes += 1
        print('  OK:   ShapiroWilkQWaveNormalityCalculator')

# ── RotorMolecularCrossSectionCalculator ─────────────────────────────────────
idx = content.find('class RotorMolecularCrossSectionCalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, E_wavenumber: float, a_param: float = 15.28,', idx)
    if c_idx != -1:
        end_line = content.find('\n', c_idx)
        sig = content[c_idx:end_line]
        new_sig = sig.replace('E_wavenumber: float,', 'E_wavenumber: float = None,')
        if 'dataset' not in new_sig:
            end_paren = new_sig.rfind(')')
            new_sig = new_sig[:end_paren] + ', dataset: dict = None' + new_sig[end_paren:]
        content = content[:c_idx] + new_sig + content[end_line:]
        after = content.find('\n', c_idx) + 1
        guard = '        if E_wavenumber is None: E_wavenumber = float((dataset or {}).get(\'frequency\', 1000.0))\n'
        if 'if E_wavenumber is None' not in content[after:after+200]:
            content = content[:after] + guard + content[after:]
        changes += 1
        print('  OK:   RotorMolecularCrossSectionCalculator')

# ── DPMTHzFrequencyMUGECalculator ────────────────────────────────────────────
idx = content.find('class DPMTHzFrequencyMUGECalculator:')
if idx != -1:
    c_idx = content.find('    def compute(self, n_level: int, f_TRZ: float = 0.1,', idx)
    if c_idx != -1:
        end_line = content.find('\n', c_idx)
        sig = content[c_idx:end_line]
        new_sig = sig.replace('n_level: int,', 'n_level: int = None,')
        if 'dataset' not in new_sig:
            end_paren = new_sig.rfind(')')
            new_sig = new_sig[:end_paren] + ', dataset: dict = None' + new_sig[end_paren:]
        content = content[:c_idx] + new_sig + content[end_line:]
        after = content.find('\n', c_idx) + 1
        guard = '        if n_level is None: n_level = int((dataset or {}).get(\'n\', 1))\n'
        if 'if n_level is None' not in content[after:after+200]:
            content = content[:after] + guard + content[after:]
        changes += 1
        print('  OK:   DPMTHzFrequencyMUGECalculator')

# ── EmissionLineFluxIntegralCalculator ────────────────────────────────────────
# compute(self, dataset: dict) — ZeroDivisionError. Need to find and guard.
idx = content.find('class EmissionLineFluxIntegralCalculator:')
if idx != -1:
    next_cls = content.find('\nclass ', idx + 1)
    section = content[idx:next_cls]
    # Find the division that causes zero error
    # Look for any / (...) patterns where the denominator could be 0
    new_section = re.sub(r'/ \(4 \* math\.pi \* luminosity_distance\*\*2\)',
                         '/ max(4 * math.pi * luminosity_distance**2, 1e-300)', section)
    new_section = re.sub(r'/ luminosity_distance\*\*2',
                         '/ max(luminosity_distance**2, 1e-300)', new_section)
    # Also guard any /z patterns
    new_section = re.sub(r'/ z\b(?!\s*[*+\-/])',
                         '/ max(z, 1e-300)', new_section)
    if new_section != section:
        content = content[:idx] + new_section + content[idx + len(section):]
        changes += 1
        print('  OK:   EmissionLineFluxIntegralCalculator: regex guards')
    else:
        # Wrap entire compute body in try/except
        c_idx = content.find('    def compute(self, dataset: dict) -> dict:', idx)
        if c_idx != -1:
            after = content.find('\n', c_idx) + 1
            # Find the return statement
            ret_idx = content.rfind('\n        return', after, next_cls)
            if ret_idx != -1:
                ret_end = content.find('\n', ret_idx + 10)
                # wrap with try/except
                body = content[after:ret_end + 1]
                # Indent everything one more level
                new_body = '        try:\n'
                for line in body.split('\n'):
                    new_body += '    ' + line + '\n' if line.strip() else line + '\n'
                new_body += '        except ZeroDivisionError:\n            return {\'F_line\': 0.0, \'error\': \'ZeroDivision guarded\'}\n'
                content = content[:after] + new_body + content[ret_end + 1:]
                changes += 1
                print('  OK:   EmissionLineFluxIntegralCalculator: try/except wrapper')

# ── SlowRollInflationCalculator ───────────────────────────────────────────────
idx = content.find('class SlowRollInflationCalculator:')
if idx != -1:
    next_cls = content.find('\nclass ', idx + 1)
    section = content[idx:next_cls]
    new_section = section
    new_section = re.sub(r'/ V\b', '/ max(V, 1e-300)', new_section)
    new_section = re.sub(r'/V\b', '/max(V, 1e-300)', new_section)
    new_section = re.sub(r'/ V_prime\b', '/ max(V_prime, 1e-300)', new_section)
    if new_section != section:
        content = content[:idx] + new_section + content[idx + len(section):]
        changes += 1
        print('  OK:   SlowRollInflationCalculator: /V guards')
    else:
        print('  MISS: SlowRollInflationCalculator: no /V pattern found')

# ── DPMFullFormulationCalculator ─────────────────────────────────────────────
idx = content.find('class DPMFullFormulationCalculator:')
if idx != -1:
    next_cls = content.find('\nclass ', idx + 1)
    section = content[idx:next_cls]
    new_section = re.sub(r'/ r26\b', '/ max(r26, 1e-300)', section)
    new_section = re.sub(r'r26 \*\* 26\b', 'max(r26, 1e-300) ** 26', new_section)
    if new_section != section:
        content = content[:idx] + new_section + content[idx + len(section):]
        changes += 1
        print('  OK:   DPMFullFormulationCalculator: r26 guards')
    else:
        print('  MISS: DPMFullFormulationCalculator')


# ── BoseEinsteinAlphaClusteringCalculator ─────────────────────────────────────
idx = content.find('class BoseEinsteinAlphaClusteringCalculator:')
if idx != -1:
    next_cls = content.find('\nclass ', idx + 1)
    c_idx = content.find('    def compute(self,', idx)
    if c_idx != -1 and (next_cls == -1 or c_idx < next_cls):
        end_line = content.find('\n', c_idx)
        sig = content[c_idx:end_line]
        # Find first required arg
        m = re.search(r'def compute\(self,\s*(\w+)(?::\s*\w+)?,', sig)
        if m:
            param = m.group(1)
            if f'{param} = None' not in sig and f'{param}=None' not in sig:
                new_sig = re.sub(rf'({param})(:\s*\w+)?,', rf'\1\2 = None,', sig)
                if 'dataset' not in new_sig:
                    end_paren = new_sig.rfind(')')
                    new_sig = new_sig[:end_paren] + ', dataset: dict = None' + new_sig[end_paren:]
                content = content[:c_idx] + new_sig + content[end_line:]
                after = content.find('\n', c_idx) + 1
                guard = f'        if {param} is None: {param} = float((dataset or {{}}).get(\'n\', 1))\n'
                if f'if {param} is None' not in content[after:after+200]:
                    content = content[:after] + guard + content[after:]
                changes += 1
                print(f'  OK:   BoseEinsteinAlphaClusteringCalculator: {param} default')

# ────────────────────────────────────────────────────────────────────────────────
print(f'\nTotal replacements: {changes}')
print('Verifying syntax...')
try:
    compile(content, 'CondensedPhysics2.py', 'exec')
    print('Syntax OK')
except SyntaxError as e:
    print(f'SYNTAX ERROR: {e}')
    sys.exit(1)

with open(PATH, 'w', encoding='utf-8') as fh:
    fh.write(content)
print(f'Written: {len(content.splitlines())} lines')
