"""Append Session 125 Phase C + Session 126 calculators to CondensedPhysics2.py."""
import os

APPEND_TEXT = '''

class UQFFBuoyancyCNBCalculator:
    """PAPER_480 — 6-system UQFFBuoyancyCNBModule CNB neutrino coupling: CentaurusA + 5 prior + Sweet/Kozima LENR."""
    def calculate(self, system='CentaurusA'):
        k_neutrino = 1e-10
        sigma_CNB  = 1e-49     # m^2
        n_CNB      = 3.36e8    # m^-3
        E_CNB      = 2.69e-23  # J
        F_nu   = k_neutrino * sigma_CNB * n_CNB * E_CNB
        F_Sweet = 1e-25 * 7.09e-36
        F_Koz   = 1e-18 * 1e-4
        return {
            'primary_equations': {
                'F_neutrino_CNB': F_nu,
                'F_Sweet_vacuum': F_Sweet,
                'F_Kozima_TNCF':  F_Koz,
                'CentaurusA_M': 1.094e38, 'CentaurusA_r': 6.17e17,
            },
            'available_equations': [f'F_nu={F_nu:.3e} N', f'F_Sweet={F_Sweet:.3e} N', f'F_Koz={F_Koz:.3e} N'],
            'simulation_set': [{'system': system, 'F_nu': F_nu}],
        }


# --- SESSION 125 CP2 Registry ---
SOURCE_SESSION125_CP2 = {
    'UQFFBuoyancyAstroCalculator': UQFFBuoyancyAstroCalculator(),
    'UQFFBuoyancyCNBCalculator':   UQFFBuoyancyCNBCalculator(),
}


# --- SESSION 126: Individual + Multi-System UQFF Module Calculators ---
class IndividualSystemUQFF18Calculator:
    """PAPER_481 — Oct 2025 batch: 16 new UQFF modules + Abell2256 + CentaurusA = 18 total systems."""
    SYSTEMS = {
        'ASASSN14li':    {'M': 1.989e37, 'r': 3.09e18, 'omega0': 1e-12, 'L_X': 1e37},
        'CrabNebula':    {'M': 1e31,     'r': 4.73e16, 'omega0': 1e-12, 'L_X': 1e27},
        'ElGordo':       {'M': 4.97e45,  'r': 3.09e22, 'omega0': 1e-15, 'L_X': 2e38},
        'ESO137':        {'M': 2e41,     'r': 6.17e21, 'omega0': 1e-15, 'L_X': 1e34},
        'IC2163':        {'M': 1.989e40, 'r': 3.09e20, 'omega0': 1e-12, 'L_X': 1e37},
        'J1610':         {'M': 1.73e40,  'r': 9.63e20, 'omega0': 1e-15, 'L_X': 1e39},
        'JupiterAurorae':{'M': 1.898e27, 'r': 7.1492e7,'omega0': 1e-12, 'L_X': 1e26},
        'LagoonNebula':  {'M': 1e36,     'r': 2.36e17, 'omega0': 1e-12, 'L_X': 1e32},
        'M87Jet':        {'M': 1.29e40,  'r': 4.63e19, 'omega0': 1e-15, 'L_X': 1e34},
        'NGC1365':       {'M': 7.17e41,  'r': 9.46e20, 'omega0': 1e-15, 'L_X': 1e36},
        'NGC2207':       {'M': 3.978e40, 'r': 4.40e20, 'omega0': 1e-12, 'L_X': 1e37},
        'RAquarii':      {'M': 3.978e30, 'r': 2.18e15, 'omega0': 1e-12, 'L_X': 1e32},
        'SgrAStar':      {'M': 8.56e36,  'r': 6.17e18, 'omega0': 1e-15, 'L_X': 1e36},
        'SPTCLJ2215':    {'M': 1.46e45,  'r': 3.09e22, 'omega0': 1e-15, 'L_X': 2e38},
        'StephanQuintet':{'M': 2e39,     'r': 3.09e22, 'omega0': 1e-15, 'L_X': 1e38},
        'VelaPulsar':    {'M': 2.8e30,   'r': 1.7e17,  'omega0': 1e-12, 'L_X': 1e27},
        'Abell2256':     {'M': 1.23e45,  'r': 3.93e22, 'omega0': 1e-15, 'L_X': 3.7e37},
        'CentaurusA':    {'M': 1.094e38, 'r': 6.17e17, 'omega0': 1e-12, 'L_X': 1e35},
    }

    def calculate(self, system='SgrAStar'):
        p = self.SYSTEMS.get(system, self.SYSTEMS['SgrAStar'])
        k_LENR, omega_LENR = 1e-10, 7.854e12
        E_cm = 3.0264e-8  # 189 GeV LEP reference
        F_LENR = k_LENR * (omega_LENR / p['omega0'])**2
        g_base = -6.6743e-11 * p['M'] / p['r']**2
        F_rel  = k_LENR * (1.24e24 / E_cm)**2
        return {
            'primary_equations': {
                'F_LENR': F_LENR, 'g_base': g_base, 'F_rel': F_rel,
                'system': system, 'params': p,
                'E_cm_LEP': E_cm,
            },
            'available_equations': [f'F_LENR={F_LENR:.3e} N', f'g={g_base:.3e} m/s2', f'F_rel={F_rel:.3e} N'],
            'simulation_set': list(self.SYSTEMS.keys()),
        }


class MultiSystemUQFFCompilerCalculator:
    """PAPER_483 — multi-system UQFF compilers: 4/5/7/8 AstroSystems plus computeGasNebulaIntegration and computeTriadicSolution."""
    REGISTRIES = {
        '4AstroSystems':  ['NGC685', 'NGC3507', 'NGC3511', 'AT2024tvd'],
        '5NebulaSystems': ['NGC3596', 'NGC1961', 'NGC5335', 'NGC2014', 'Carina'],
        '7NebulaSystems': ['NGC685', 'NGC3507', 'NGC3511', 'Carina', 'NGC1961', 'NGC5335', 'NGC2014'],
        '8AstroSystems':  ['NGC4826', 'NGC1805', 'NGC6307', 'NGC7027', 'NGC685', 'NGC3507', 'NGC3511', 'AT2024tvd'],
        'BuoyancyMulti':  ['M74', 'M16', 'M84', 'CentaurusA'],
    }

    def calculate(self, registry='8AstroSystems'):
        systems = self.REGISTRIES.get(registry, self.REGISTRIES['8AstroSystems'])
        # Source160 universal validation result for all template-mass systems
        F_UBi_i_template = -8.32e217    # N (real component)
        F_UBi_i_imag     = -6.75e160    # N (imaginary component)
        compressed_integrand = 6.16e39  # N (template mass)
        DPM_resonance = 1.76e17
        buoyancy_Ub1 = 6e-19
        superconductive_Ui = 1.7e-43
        return {
            'primary_equations': {
                'F_UBi_i_template_real': F_UBi_i_template,
                'F_UBi_i_template_imag': F_UBi_i_imag,
                'compressed_integrand':  compressed_integrand,
                'DPM_resonance':         DPM_resonance,
                'buoyancy_Ub1':          buoyancy_Ub1,
                'superconductive_Ui':    superconductive_Ui,
                'registry':              registry,
                'n_systems':             len(systems),
                'systems':               systems,
            },
            'available_equations': [
                f'F={F_UBi_i_template:.3e} + i{F_UBi_i_imag:.3e} N (Source160 validation)',
                'computeGasNebulaIntegration(system) => F_gas = rho_gas * g_UQFF * V_nebula',
                'computeTriadicSolution(system, t) => (F_compressed, F_resonant, F_buoyancy)',
                'computeMasterEquations(system) => dispatch to module-specific solver',
            ],
            'simulation_set': [{'registry': registry, 'system': s} for s in systems],
        }


class HydrogenResonanceUQFFCalculator:
    """PAPER_482 — HydrogenResonanceUQFFModule PTOE nuclear resonance: H_res(Z,A,t) for Z=1 to 118."""
    def calculate(self, Z=1, A=1, t=1e-15):
        import math
        # Constants from HydrogenResonanceUQFFModule.cpp
        k_A    = 0.4604
        k_dp   = 1.325e-6
        f_dp   = 1e15         # Hz dipole frequency
        E_bind = 7.8e6 * 1.602e-19  # J (hydrogen placeholder, 7.8 MeV)
        h      = 6.626e-34    # J*s
        x2_per_ZA = -1.35e172  # scaling

        # Compute components
        A_H    = 1  # reference mass number
        A_res  = k_A * Z * (A / A_H)  # amplitude resonance
        f_res  = (E_bind / h) * (A_H / A)  # resonance frequency
        sin_t  = math.sin(2.0 * math.pi * f_res * t)
        U_dp   = k_dp * A / (f_dp**2)
        SC_m   = 1.0  # superconductive parameter (normalized)
        k_nuc  = (A - Z) / max(Z, 1)  # neutron-to-proton ratio
        # Shell correction estimate: 0.1 * (Z_magic + N_magic) for magic numbers
        S_shell = 0.1 * (2 + 2)

        H_res_integrand = A_res * sin_t + U_dp * SC_m * k_nuc + S_shell
        x2 = x2_per_ZA * (Z + A)
        H_res = H_res_integrand * x2

        # Nuclear gravity analog
        G      = 6.6743e-11
        m_nuc  = A * 1.6735e-27  # kg
        rho_nuc = 1e17            # kg/m^3 nuclear density
        r_nuc  = (3 * m_nuc / (4 * math.pi * rho_nuc))**(1/3)
        g_nuc  = -G * m_nuc * rho_nuc / r_nuc

        return {
            'primary_equations': {
                'H_res':           H_res,
                'H_res_integrand': H_res_integrand,
                'A_res':           A_res,
                'f_res_Hz':        f_res,
                'U_dp':            U_dp,
                'k_nuc':           k_nuc,
                'S_shell':         S_shell,
                'g_nuc':           g_nuc,
                'r_nuc_m':         r_nuc,
                'Z': Z, 'A': A, 't': t,
            },
            'available_equations': [
                f'H_res = {H_res:.3e}',
                f'A_res = k_A*Z*(A/A_H) = {A_res:.4f}',
                f'f_res = (E_bind/h)*(A_H/A) = {f_res:.3e} Hz',
                f'g_nuc = -G*M_nuc*rho_nuc/r_nuc = {g_nuc:.3e} m/s2',
                'g_nuc full: -G*M_nuc*rho_nuc/r_nuc - k_B*T*rho/(m_e*c^2) + kappa*c^4/(G*r^2)',
                'Scale: Z=1 (Protium, r=1e-15m) to Z=92 (Uranium) to galactic r=1e22m (57 decades)',
            ],
            'simulation_set': [
                {'Z': z, 'A': round(z * 2.2)} for z in [1, 2, 6, 8, 14, 26, 50, 79, 82, 92, 118]
            ],
        }


# --- SESSION 126 CP2 Registry ---
SOURCE_SESSION126_CP2 = {
    'IndividualSystemUQFF18Calculator':      IndividualSystemUQFF18Calculator(),
    'MultiSystemUQFFCompilerCalculator':     MultiSystemUQFFCompilerCalculator(),
    'HydrogenResonanceUQFFCalculator':       HydrogenResonanceUQFFCalculator(),
}
'''

cp2_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CondensedPhysics2.py')

with open(cp2_path, 'a', encoding='utf-8') as f:
    f.write(APPEND_TEXT)

with open(cp2_path, encoding='utf-8') as f:
    lines = f.readlines()

print(f"CondensedPhysics2.py now has {len(lines)} lines")
print("Last 5 lines:")
for line in lines[-5:]:
    print(repr(line))
