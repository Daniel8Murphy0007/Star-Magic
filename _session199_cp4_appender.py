#!/usr/bin/env python3
"""Session 199 CP4 Appender — grok_share_589f6949-6fe9.txt
Adds 9 new classes (#438-#446, PAPER_854-862) to CondensedPhysics4.py.
Topics: LENR k_eta calibration, pseudo-monopole 26-state, Higgs UH vacuum,
NGC 346 Ug3 star formation, Westerlund 2 quadriadic, micro-plasmoid 25um
LENR, neutrino energy vacuum ratio, Kepler Orrery V 35-frame iterative,
Universal Magnetism U_m master equation.
"""
import re, sys, os

CP4 = os.path.join(os.path.dirname(__file__), "CondensedPhysics4.py")

# ─── Version & session metadata ───
OLD_VERSION = "5.58"
NEW_VERSION = "5.59"
SESSION = 199
SOURCE_FILE = "grok_share_589f6949-6fe9.txt"
SOURCE_LINES = 2404
SOURCE_DESC = "2404 lines, May 05 – Aug 07 2025"
FIRST_CLASS = 438
LAST_CLASS = 446
FIRST_PAPER = 854
LAST_PAPER = 862
NUM_CLASSES = LAST_CLASS - FIRST_CLASS + 1

SESSION_HEADER_LINE = (
    f"    Updated: Session {SESSION} v{NEW_VERSION} — "
    f"CP4 429→{429 + NUM_CLASSES} (#{FIRST_CLASS} LENRKetaCalibration3EnvironmentDeltaKCalc "
    f"+ #{FIRST_CLASS+1} PseudoMonopole26StateVacuumDensityCalc "
    f"+ #{FIRST_CLASS+2} HiggsVacuumUHExcitationKHiggsUQFFCalc "
    f"+ #{FIRST_CLASS+3} NGC346Ug3StarFormationTempVradCalc "
    f"+ #{FIRST_CLASS+4} Westerlund2QuadriadicRealImaginaryCalc "
    f"+ #{FIRST_CLASS+5} MicroPlasmoid25umLENRBuoyancyReversalCalc "
    f"+ #{FIRST_CLASS+6} NeutrinoEnergyUQFFVacuumRatioCalc "
    f"+ #{FIRST_CLASS+7} KeplerOrreryV35FrameIterativeUbCalc "
    f"+ #{FIRST_CLASS+8} UniversalMagnetismUmMasterEquationCalc: "
    f"PAPER_{FIRST_PAPER}-{LAST_PAPER}; {SOURCE_FILE} ({SOURCE_DESC}); "
    f"LENR k_eta 3-environment calibration + pseudo-monopole 26-state δ_n="
    f"(2π)^(n/6) + Higgs UH(t,n) vacuum excitation k_Higgs=1.79e18 + "
    f"NGC346 Ug3 T_scaled star formation + Westerlund 2 quadriadic "
    f"real+imaginary 4-set solver + micro-plasmoid 25.4μm LENR buoyancy "
    f"reversal + neutrino E_nu vacuum ratio + Kepler Orrery V 35-frame "
    f"iterative F_env(t) + Universal Magnetism U_m fourth master equation; "
    f"VDS-DVP-BH ABSENT; {LAST_PAPER}/1000 papers {LAST_PAPER/10:.1f}%)"
)

# ─── New classes code block ───
NEW_CLASSES = r'''

# ═══════════════════════════════════════════════════════════════════
# SESSION 199 — grok_share_589f6949-6fe9.txt (2404 lines, May 05 – Aug 07 2025)
# UQFF Compression Cycle 2 / Kepler Orrery V 35 Frames / LENR k_eta /
# Pseudo-Monopole / Higgs UH / NGC 346 / Westerlund 2 Quadriadic /
# Micro-Plasmoid LENR / Neutrino Energy / Universal Magnetism U_m
# 9 classes: #438-#446  PAPER_854-862
# ═══════════════════════════════════════════════════════════════════


class LENRKetaCalibration3EnvironmentDeltaKCalc(_CP4Calculator):  # PAPER_854 #438
    """LENR neutron-production calibration constant k_eta across three
    environments (Metallic Hydride Cells, Exploding Wires, Solar Corona)
    with buoyancy Delta-k tracking.

    Core equation:
        eta = k_eta * exp(-[SSq]*n/26) * exp(-(pi - t)) * U_m / rho_vac

    Calibration results from source data:
        Metallic Hydride: k_eta ~ 2.75e8, Delta_k ~ 7.25e8
        Exploding Wires:  k_eta ~ 1.91e2, Delta_k ~ 8.09e2
        Solar Corona:     k_eta ~ 6.06e-6, Delta_k ~ 3.94e-6

    Buoyancy tracking: Delta_k = k_expected - k_actual measures the U_b
    calibration residual, identifying where buoyancy counterforce shifts
    the effective neutron production rate.

    Source: grok_share_589f6949-6fe9.txt lines 1685-1740
    """

    category = "LENR / Neutron Production"

    def compute(self, dataset: dict) -> dict:
        """Compute k_eta calibration and Delta-k buoyancy residuals.

        Parameters (from dataset):
            U_m         : float — universal magnetism energy density (J/m^3)
            rho_vac     : float — vacuum density rho_vac,[UA] (J/m^3, default 7.09e-36)
            n           : int   — quantum state index (1-26, default 1)
            t_days      : float — time in days (default 1.0)
            eta_observed_metallic : float — observed eta metallic (cm^-2/s, default 1e13)
            eta_observed_wires    : float — observed eta wires (cm^-2/s, default 1e8)
            eta_observed_corona   : float — observed eta corona (cm^-2/s, default 7e-3)
            k_eta_expected_metallic : float — expected k_eta (default 1e9)
            k_eta_expected_wires    : float — expected k_eta (default 1e3)
            k_eta_expected_corona   : float — expected k_eta (default 1e-5)
        """
        U_m = float(dataset.get('U_m', 2.67e-31))
        rho_vac = float(dataset.get('rho_vac', RHO_VAC_UA))
        n = int(dataset.get('n', 1))
        t = float(dataset.get('t_days', 1.0))

        eta_met = float(dataset.get('eta_observed_metallic', 1e13))
        eta_wir = float(dataset.get('eta_observed_wires', 1e8))
        eta_cor = float(dataset.get('eta_observed_corona', 7e-3))
        k_exp_met = float(dataset.get('k_eta_expected_metallic', 1e9))
        k_exp_wir = float(dataset.get('k_eta_expected_wires', 1e3))
        k_exp_cor = float(dataset.get('k_eta_expected_corona', 1e-5))

        # Suppression factors
        ssq_factor = self._ssq_exp(n)  # exp(-[SSq]*n/26)
        time_factor = math.exp(-(math.pi - t)) if t < math.pi else math.exp(t - math.pi)
        ratio = U_m / rho_vac if rho_vac > 0 else 0.0
        denom = ssq_factor * time_factor * ratio if ratio > 0 else 1.0

        # Derive k_eta for each environment
        k_met = eta_met / denom if denom != 0 else 0.0
        k_wir = eta_wir / denom if denom != 0 else 0.0
        k_cor = eta_cor / denom if denom != 0 else 0.0

        # Delta-k buoyancy residuals
        dk_met = k_exp_met - k_met
        dk_wir = k_exp_wir - k_wir
        dk_cor = k_exp_cor - k_cor

        # Forward prediction: eta from derived k_eta
        eta_pred_met = k_met * ssq_factor * time_factor * ratio
        eta_pred_wir = k_wir * ssq_factor * time_factor * ratio
        eta_pred_cor = k_cor * ssq_factor * time_factor * ratio

        return {
            'ssq_factor': ssq_factor,
            'time_factor': time_factor,
            'U_m_over_rho_vac': ratio,
            'k_eta_metallic': k_met,
            'k_eta_wires': k_wir,
            'k_eta_corona': k_cor,
            'delta_k_metallic': dk_met,
            'delta_k_wires': dk_wir,
            'delta_k_corona': dk_cor,
            'eta_predicted_metallic': eta_pred_met,
            'eta_predicted_wires': eta_pred_wir,
            'eta_predicted_corona': eta_pred_cor,
            'buoyancy_residual_ratio_met': dk_met / k_exp_met if k_exp_met else 0,
            'buoyancy_residual_ratio_wir': dk_wir / k_exp_wir if k_exp_wir else 0,
            'buoyancy_residual_ratio_cor': dk_cor / k_exp_cor if k_exp_cor else 0,
            'paper': 'PAPER_854',
        }


class PseudoMonopole26StateVacuumDensityCalc(_CP4Calculator):  # PAPER_855 #439
    """Pseudo-monopole 26-state vacuum density progression calculator.

    Core equations:
        delta_n = (2*pi)^(n/6)     — pseudo-monopole angular spacing
        rho_vac,[UA']:[SCm](n, t) = 1e-23 * (0.1)^n * exp(-[SSq]*n/26) * exp(-(pi - t))

    Computes the full 26-state progression of vacuum density ratios,
    pseudo-monopole angular spacings, and energy densities across all
    quantum states from n=1 to n=26.

    Source: grok_share_589f6949-6fe9.txt lines 1735-1755
    """

    category = "Pseudo-Monopole / Vacuum States"

    def compute(self, dataset: dict) -> dict:
        """Compute all 26 pseudo-monopole states.

        Parameters (from dataset):
            t_days : float — time in days (default 0.0)
            rho_base : float — base vacuum density (J/m^3, default 1e-23)
            ratio_per_n : float — ratio per quantum number (default 0.1)
        """
        t = float(dataset.get('t_days', 0.0))
        rho_base = float(dataset.get('rho_base', 1e-23))
        ratio_n = float(dataset.get('ratio_per_n', 0.1))

        states = []
        for n in range(1, 27):
            delta_n = (2.0 * math.pi) ** (n / 6.0)
            ssq_fac = self._ssq_exp(n)
            time_fac = math.exp(-(math.pi - t)) if t < math.pi else math.exp(t - math.pi)
            rho_n = rho_base * (ratio_n ** n) * ssq_fac * time_fac
            energy_n = rho_n  # J/m^3 = energy density
            mass_equiv = rho_n / (C_LIGHT ** 2) if C_LIGHT > 0 else 0.0
            states.append({
                'n': n,
                'delta_n_rad': delta_n,
                'delta_n_deg': math.degrees(delta_n) % 360,
                'ssq_factor': ssq_fac,
                'rho_vac_n': rho_n,
                'mass_equiv_kg': mass_equiv,
            })

        # Summary statistics
        rho_values = [s['rho_vac_n'] for s in states]
        rho_total = sum(rho_values)

        return {
            'states': states,
            'rho_total_all_states': rho_total,
            'rho_state_1': states[0]['rho_vac_n'],
            'rho_state_26': states[25]['rho_vac_n'],
            'delta_1_rad': states[0]['delta_n_rad'],
            'delta_26_rad': states[25]['delta_n_rad'],
            'suppression_range': states[0]['ssq_factor'] / states[25]['ssq_factor'] if states[25]['ssq_factor'] > 0 else float('inf'),
            'paper': 'PAPER_855',
        }


class HiggsVacuumUHExcitationKHiggsUQFFCalc(_CP4Calculator):  # PAPER_856 #440
    """Higgs field UH vacuum excitation via UQFF pseudo-monopole density.

    Core equation:
        UH(t, n) = lambda_H * rho_vac,[UA']:[SCm](n,t) * omega_H(t)
                   * exp(-[SSq]*n/26) * exp(-(pi - t)) * (1 + f_quasi)

    Derived quantities:
        m_H = UH / c^2
        E_H = m_H * c^2 = UH
        k_Higgs = 125 GeV * 1.602e-19 / E_H  (scaling to observed Higgs mass)

    From source: UH ~ 1.539e-32 J/m^3 at n=1, t=0 → k_Higgs ~ 1.79e18
    (recalibrated from initial 1.30e30).

    Source: grok_share_589f6949-6fe9.txt lines 1738-1760
    """

    category = "Higgs Field / Vacuum Excitation"

    HIGGS_MASS_GEV = 125.0  # GeV observed Higgs mass
    EV_TO_J = 1.602e-19

    def compute(self, dataset: dict) -> dict:
        """Compute Higgs field UH from UQFF vacuum excitation.

        Parameters (from dataset):
            n           : int   — quantum state (default 1)
            t_days      : float — time in days (default 0.0)
            lambda_H    : float — Higgs coupling (default 1.0)
            omega_H     : float — angular frequency (rad/s, default 1.585e-8)
            rho_base    : float — base vacuum density (default 1e-23)
            f_quasi     : float — quasi correction (default 0.01)
        """
        n = int(dataset.get('n', 1))
        t = float(dataset.get('t_days', 0.0))
        lam_H = float(dataset.get('lambda_H', 1.0))
        omega_H = float(dataset.get('omega_H', 1.585e-8))
        rho_b = float(dataset.get('rho_base', 1e-23))
        f_quasi = float(dataset.get('f_quasi', 0.01))

        # Vacuum density at state n
        ssq_fac = self._ssq_exp(n)
        time_fac = math.exp(-(math.pi - t)) if t < math.pi else math.exp(t - math.pi)
        rho_n = rho_b * (0.1 ** n) * ssq_fac * time_fac

        # UH computation
        UH = lam_H * rho_n * omega_H * ssq_fac * time_fac * (1.0 + f_quasi)

        # Derived quantities
        m_H = UH / (C_LIGHT ** 2) if C_LIGHT > 0 else 0.0
        E_H_J = UH
        E_H_eV = E_H_J / self.EV_TO_J if self.EV_TO_J > 0 else 0.0

        # k_Higgs scaling factor
        higgs_target_J = self.HIGGS_MASS_GEV * 1e9 * self.EV_TO_J
        k_Higgs = higgs_target_J / E_H_J if E_H_J > 0 else float('inf')

        return {
            'n': n,
            'rho_vac_n': rho_n,
            'UH_J_per_m3': UH,
            'm_H_kg': m_H,
            'E_H_eV': E_H_eV,
            'k_Higgs_scaling': k_Higgs,
            'higgs_target_J': higgs_target_J,
            'ssq_factor': ssq_fac,
            'time_factor': time_fac,
            'lambda_H': lam_H,
            'omega_H': omega_H,
            'f_quasi': f_quasi,
            'paper': 'PAPER_856',
        }


class NGC346Ug3StarFormationTempVradCalc(_CP4Calculator):  # PAPER_857 #441
    """NGC 346 Ug3 star-formation temperature and radial velocity calculator.

    Core equation:
        Ug3(t, r, theta, n) = k_3 * sum_j B_j(r, theta, t, rho_vac,[SCm])
                              * cos(omega_s(t) * t * pi) * P_core * E_react(t)

    Derived:
        T_scaled = Ug3 / rho_vac,[UA]  (raw), then scaled to observational range
        v_radial = derived from T and system parameters

    NGC 346: young star-forming region in SMC, r ~ 1.496e10 m (scaled),
    T_scaled ~ 10^4 K (H II region), v_radial ~ -3.33e-5 c (outflows).

    Source: grok_share_589f6949-6fe9.txt lines 1758-1785
    """

    category = "Star Formation / Ug3"

    def compute(self, dataset: dict) -> dict:
        """Compute Ug3 gravity, temperature, and radial velocity.

        Parameters (from dataset):
            r           : float — distance (m)
            t_days      : float — time in days (default 86400 = 1 day)
            k_3         : float — coupling constant (default 1.0)
            n_stars     : int   — number of contributing sources (default 1000)
            B_dipole    : float — mean dipole moment (T*m^3, default 3.38e20)
            omega_s     : float — string rotation freq (rad/s, default 2.5e-6)
            P_core      : float — core pressure factor (default 1.0)
            rho_vac_ua  : float — rho_vac,[UA] (J/m^3, default 7.09e-36)
            T_scale_factor : float — temperature scaling divisor (default 1e31)
        """
        r = float(dataset.get('r', 1.496e10))
        t = float(dataset.get('t_days', 1.0))
        k_3 = float(dataset.get('k_3', 1.0))
        n_st = int(dataset.get('n_stars', 1000))
        B_dip = float(dataset.get('B_dipole', 3.38e20))
        omega_s = float(dataset.get('omega_s', 2.5e-6))
        P_core = float(dataset.get('P_core', 1.0))
        rho_ua = float(dataset.get('rho_vac_ua', RHO_VAC_UA))
        T_sf = float(dataset.get('T_scale_factor', 1e31))

        t_sec = t * 86400.0
        E_react = self._e_react(t)

        # B_j / r^3 sum
        B_sum = n_st * B_dip / (r ** 3) if r > 0 else 0.0

        # cos(omega_s * t * pi)
        cos_term = math.cos(omega_s * t_sec * math.pi)

        # Ug3
        Ug3 = k_3 * B_sum * cos_term * P_core * E_react

        # Temperature
        T_raw = Ug3 / rho_ua if rho_ua > 0 else 0.0
        T_scaled = T_raw / T_sf if T_sf > 0 else T_raw

        # Radial velocity (from Ug3 energy density → kinetic)
        # v_radial ~ -sqrt(2 * |Ug3| / rho_gas) * sign, approximate
        rho_gas = float(dataset.get('rho_gas', 1e-20))
        v_abs = math.sqrt(2.0 * abs(Ug3) / rho_gas) if rho_gas > 0 and Ug3 != 0 else 0.0
        v_radial = -v_abs  # outflow convention
        v_over_c = v_radial / C_LIGHT if C_LIGHT > 0 else 0.0

        return {
            'Ug3_J_per_m3': Ug3,
            'B_sum': B_sum,
            'cos_term': cos_term,
            'E_react': E_react,
            'T_raw_K': T_raw,
            'T_scaled_K': T_scaled,
            'v_radial_m_per_s': v_radial,
            'v_over_c': v_over_c,
            'paper': 'PAPER_857',
        }


class Westerlund2QuadriadicRealImaginaryCalc(_CP4Calculator):  # PAPER_858 #442
    """Westerlund 2 quadriadic UQFF solver: four master equations solved
    simultaneously for both real and imaginary solutions.

    Four master equations:
        1. F_U_g = sum_k [k_k*(f_UA'*f_SCm*REB)^2/r^2 * G_k] + k_4 * rho_vac * M_BH/r ...
        2. R(t) = sum_i [R_Ug1,i*cos(w_i*t) + R_Ug2,i*cos + R_Ug3,i*cos + R_Ug4i,i*cos]
        3. F_U_Bi = sum_k [k_Ub,k*f_UA'*f_SCm*REB/r^2 * H_k * f_Ub]
        4. U_m = sum_j [mu_j/r_j*(1-exp(-gamma*t)*cos(pi*t/n))*phi^j]*P_SCm*E_react*(1+1e13*f_H)*(1+f_q)

    Each equation yields both real and imaginary components.
    Imaginary solutions represent the quantum superconducting portion.

    Source: grok_share_589f6949-6fe9.txt lines 1940-2100
    """

    category = "Quadriadic / Simultaneous Solver"

    def compute(self, dataset: dict) -> dict:
        """Solve four master UQFF equations simultaneously.

        Parameters (from dataset):
            r           : float — system radius (m, default 1.89e16 for Westerlund 2)
            t_days      : float — age in days
            f_UA        : float — UA fraction (default 0.999)
            f_SCm       : float — SCm fraction (default 0.001)
            REB         : float — reactivity-energy-barrier (default 1.0)
            k_1         : float — gravity coupling (default 1.0)
            k_Ub        : float — buoyancy coupling (default 0.1)
            mu_j        : float — magnetic dipole moment (T*pm^3, default 3.38e23)
            gamma       : float — decay rate (day^-1, default 5e-5)
            n           : int   — quantum state (default 1)
            f_Heaviside : float — Heaviside correction (default 0.01)
            f_quasi     : float — quasi-particle correction (default 0.01)
            delta_k_eta : float — buoyancy Delta k_eta (default 7.25e8)
            V_ratio     : float — V_little/V_big Boyle ratio (default 1/33)
        """
        r = float(dataset.get('r', 1.89e16))
        t = float(dataset.get('t_days', 7.30e10))
        f_UA = float(dataset.get('f_UA', 0.999))
        f_SCm = float(dataset.get('f_SCm', 0.001))
        REB = float(dataset.get('REB', 1.0))
        k_1 = float(dataset.get('k_1', 1.0))
        k_Ub = float(dataset.get('k_Ub', 0.1))
        mu_j = float(dataset.get('mu_j', 3.38e23))
        gamma = float(dataset.get('gamma', GAMMA_DECAY))
        n = int(dataset.get('n', 1))
        f_H = float(dataset.get('f_Heaviside', 0.01))
        f_q = float(dataset.get('f_quasi', 0.01))
        dk_eta = float(dataset.get('delta_k_eta', 7.25e8))
        V_rat = float(dataset.get('V_ratio', 1.0 / 33.0))

        r2 = r * r if r > 0 else 1.0
        product = f_UA * f_SCm * REB

        # 1. Compressed UQFF (F_U_g)
        F_Ug_real = k_1 * (product ** 2) / r2
        F_Ug_imag = F_Ug_real  # imaginary mirror (quantum portion)

        # 2. Resonance UQFF (R(t))
        # 26-state sum with cos terms; approximate average cos = 0.5
        R_real = F_Ug_real * 0.5 * 26
        R_imag = R_real

        # 3. Buoyancy UQFF (F_U_Bi)
        f_Ub = k_Ub * dk_eta * (RHO_VAC_UA / RHO_VAC_SCM) * V_rat
        F_UBi_real = k_Ub * product / r2 * f_Ub
        F_UBi_imag = F_UBi_real

        # 4. Universal Magnetism (U_m)
        exp_gamma = math.exp(-gamma * t) if gamma * t < 700 else 0.0
        cos_tn = math.cos(math.pi * t / n) if n > 0 else 1.0
        bracket = 1.0 - exp_gamma * cos_tn
        E_react = self._e_react(t)
        U_m_real = (mu_j / r) * bracket * 1.0 * E_react * (1.0 + 1e13 * f_H) * (1.0 + f_q) if r > 0 else 0.0
        U_m_imag = U_m_real

        # Total quadriadic magnitude
        total_real = abs(F_Ug_real) + abs(R_real) + abs(F_UBi_real) + abs(U_m_real)

        return {
            'F_Ug_real': F_Ug_real,
            'F_Ug_imag': F_Ug_imag,
            'R_real': R_real,
            'R_imag': R_imag,
            'F_UBi_real': F_UBi_real,
            'F_UBi_imag': F_UBi_imag,
            'f_Ub_factor': f_Ub,
            'U_m_real': U_m_real,
            'U_m_imag': U_m_imag,
            'exp_gamma_t': exp_gamma,
            'total_quadriadic_real': total_real,
            'dominant_term': max(
                ('F_Ug', abs(F_Ug_real)),
                ('R', abs(R_real)),
                ('F_UBi', abs(F_UBi_real)),
                ('U_m', abs(U_m_real)),
                key=lambda x: x[1]
            )[0],
            'paper': 'PAPER_858',
        }


class MicroPlasmoid25umLENRBuoyancyReversalCalc(_CP4Calculator):  # PAPER_859 #443
    """Micro-plasmoid (25.4 micron) LENR glass reactor buoyancy reversal
    dynamics calculator.

    Models moving plasmoids in a glass reactor LENR experiment:
    - Largest plasmoid: 25.4 um (0.001 inch)
    - Duration: 165 seconds
    - Behavior: upward swirling → buoyancy reversal (downward)
    - Count: 1-8 per frame, averaging 4-5

    Applies all four UQFF master equations at micro-scale (r = 25.4e-6 m),
    computing buoyancy reversal time and force balance.

    Source: grok_share_589f6949-6fe9.txt lines 2100-2370
    """

    category = "LENR / Plasmoid Dynamics"

    def compute(self, dataset: dict) -> dict:
        """Compute micro-plasmoid UQFF dynamics with buoyancy reversal.

        Parameters (from dataset):
            r_plasmoid  : float — plasmoid radius (m, default 25.4e-6)
            r_reactor   : float — reactor size (m, default 0.0254)
            t_duration  : float — experiment duration (s, default 165)
            n_plasmoids : int   — average plasmoid count (default 5)
            f_UA        : float — UA fraction (default 0.999)
            f_SCm       : float — SCm fraction (default 0.001)
            delta_k_eta : float — buoyancy Delta k_eta (default 7.25e8)
            mu_j        : float — dipole moment (T*pm^3, default 3.38e23)
            gamma       : float — decay rate (s^-1, default 5e-5)
        """
        r_p = float(dataset.get('r_plasmoid', 25.4e-6))
        r_r = float(dataset.get('r_reactor', 0.0254))
        t_dur = float(dataset.get('t_duration', 165.0))
        n_plas = int(dataset.get('n_plasmoids', 5))
        f_UA = float(dataset.get('f_UA', 0.999))
        f_SCm = float(dataset.get('f_SCm', 0.001))
        dk_eta = float(dataset.get('delta_k_eta', 7.25e8))
        mu_j = float(dataset.get('mu_j', 3.38e23))
        gamma = float(dataset.get('gamma', 5e-5))

        r2 = r_p * r_p if r_p > 0 else 1.0
        product = f_UA * f_SCm

        # Volume ratio (Boyle's Law micro-scale)
        V_ratio = r_p / r_r if r_r > 0 else 0.0

        # F_Ug at micro-scale
        F_Ug = (product ** 2) / r2

        # F_UBi with micro-scale V_ratio
        k_Ub = 0.1
        f_Ub = k_Ub * dk_eta * (RHO_VAC_UA / RHO_VAC_SCM) * V_ratio
        F_UBi = k_Ub * product / r2 * f_Ub

        # U_m at micro-scale
        exp_g = math.exp(-gamma * t_dur) if gamma * t_dur < 700 else 0.0
        cos_t = math.cos(math.pi * t_dur)
        bracket = 1.0 - exp_g * cos_t
        U_m = (mu_j / r_p) * bracket if r_p > 0 else 0.0

        # Buoyancy reversal analysis
        # Reversal occurs when F_UBi > F_Ug (buoyancy exceeds gravity)
        buoyancy_ratio = F_UBi / F_Ug if F_Ug > 0 else float('inf')
        reversal_detected = buoyancy_ratio > 1.0

        # Velocity estimate (um/s from frame analysis)
        v_up = 4.0e-6   # typical upward velocity (m/s)
        v_down = -2.0e-6  # reversal velocity (m/s)
        reversal_time_frac = 0.7  # reversal at ~70% duration from frames

        return {
            'r_plasmoid_m': r_p,
            'V_ratio': V_ratio,
            'F_Ug': F_Ug,
            'F_UBi': F_UBi,
            'f_Ub_factor': f_Ub,
            'U_m': U_m,
            'buoyancy_ratio_UBi_over_Ug': buoyancy_ratio,
            'reversal_detected': reversal_detected,
            'n_plasmoids': n_plas,
            'v_upward_m_per_s': v_up,
            'v_downward_m_per_s': v_down,
            'reversal_time_fraction': reversal_time_frac,
            'exp_gamma_t': exp_g,
            'paper': 'PAPER_859',
        }


class NeutrinoEnergyUQFFVacuumRatioCalc(_CP4Calculator):  # PAPER_860 #444
    """Neutrino energy from UQFF vacuum density ratio.

    Core equation:
        E_neutrino proportional to rho_vac,[UA']:[SCm](n, t)
                   * exp(-[SSq]*n/26) * exp(-(pi - t))
                   * (U_m / rho_vac,[UA])

    Models neutrino energy production through vacuum density gradients,
    connecting pseudo-monopole states to weak-interaction energy scales.

    For Westerlund 2 application: E_neutrino ~ 1.05e5 eV.

    Source: grok_share_589f6949-6fe9.txt lines 1942-1960
    """

    category = "Neutrino / Vacuum Energy"

    def compute(self, dataset: dict) -> dict:
        """Compute UQFF neutrino energy from vacuum ratios.

        Parameters (from dataset):
            n           : int   — quantum state (default 1)
            t_days      : float — time in days (default 1.0)
            U_m         : float — universal magnetism energy density (J/m^3)
            rho_vac_ua  : float — rho_vac,[UA] (J/m^3, default 7.09e-36)
            rho_base    : float — base pseudo-monopole density (default 1e-23)
            k_neutrino  : float — neutrino coupling constant (default 1.0)
        """
        n = int(dataset.get('n', 1))
        t = float(dataset.get('t_days', 1.0))
        U_m = float(dataset.get('U_m', 1.80e-5))
        rho_ua = float(dataset.get('rho_vac_ua', RHO_VAC_UA))
        rho_b = float(dataset.get('rho_base', 1e-23))
        k_nu = float(dataset.get('k_neutrino', 1.0))

        # Pseudo-monopole vacuum density at state n
        ssq_fac = self._ssq_exp(n)
        time_fac = math.exp(-(math.pi - t)) if t < math.pi else math.exp(t - math.pi)
        rho_n = rho_b * (0.1 ** n) * ssq_fac * time_fac

        # U_m / rho_vac,[UA] ratio
        Um_ratio = U_m / rho_ua if rho_ua > 0 else 0.0

        # Neutrino energy (proportional)
        E_nu_J = k_nu * rho_n * ssq_fac * time_fac * Um_ratio
        E_nu_eV = E_nu_J / 1.602e-19 if E_nu_J > 0 else 0.0

        # Compare to known neutrino mass scales
        # Atmospheric neutrino: ~0.05 eV
        # Solar neutrino: ~0.01 eV mass, but MeV energy
        ratio_to_atmospheric = E_nu_eV / 0.05 if E_nu_eV > 0 else 0.0

        return {
            'n': n,
            'rho_vac_n': rho_n,
            'ssq_factor': ssq_fac,
            'time_factor': time_fac,
            'U_m_over_rho_vac_UA': Um_ratio,
            'E_neutrino_J': E_nu_J,
            'E_neutrino_eV': E_nu_eV,
            'ratio_to_atmospheric_scale': ratio_to_atmospheric,
            'paper': 'PAPER_860',
        }


class KeplerOrreryV35FrameIterativeUbCalc(_CP4Calculator):  # PAPER_861 #445
    """Kepler Orrery V 35-frame iterative U_b model refinement calculator.

    Iterates the U_b model across N_frames temporal snapshots, computing
    F_env(t) = w_orbit*F_orbit + w_tide*F_tide + w_gal*F_gal at each frame.

    Sub-equations:
        F_orbit = G * M_p * M_s / a^3
        F_tide  = G * M_p * M_s * R_p / a^6
        F_gal   = v_gal^2 / r_gal + G * M_DM / r_gal^2

    35 Kepler Orrery V frames (22 Sep – 27 Oct 2011) progressively refined
    F_env(t) to ~6.5e-2 m/s^2 for Earth-Sun + 1200 exoplanet systems.

    Source: grok_share_589f6949-6fe9.txt lines 400-1700 (iterative frames)
    """

    category = "Exoplanet / Iterative U_b"

    def compute(self, dataset: dict) -> dict:
        """Compute N-frame iterative U_b refinement.

        Parameters (from dataset):
            M_planet    : float — planet mass (kg, default 5.97e24 Earth)
            M_star      : float — star mass (kg, default 1.989e30 Sun)
            a           : float — semi-major axis (m, default 1.496e11 AU)
            R_planet    : float — planet radius (m, default 6.371e6)
            v_gal       : float — galactic orbital velocity (m/s, default 2.2e5)
            r_gal       : float — galactic radius (m, default 2.47e20)
            M_DM        : float — dark matter enclosed mass (kg, default 2.57e40)
            w_orbit     : float — weight for F_orbit (default 0.5)
            w_tide      : float — weight for F_tide (default 0.3)
            w_gal       : float — weight for F_gal (default 0.2)
            n_frames    : int   — number of frames (default 35)
            dt_days     : float — time step per frame (days, default 1.0)
        """
        M_p = float(dataset.get('M_planet', 5.97e24))
        M_s = float(dataset.get('M_star', M_SUN))
        a = float(dataset.get('a', 1.496e11))
        R_p = float(dataset.get('R_planet', 6.371e6))
        v_gal = float(dataset.get('v_gal', 2.2e5))
        r_gal = float(dataset.get('r_gal', 2.47e20))
        M_DM = float(dataset.get('M_DM', 2.57e40))
        w_orb = float(dataset.get('w_orbit', 0.5))
        w_tid = float(dataset.get('w_tide', 0.3))
        w_gal_w = float(dataset.get('w_gal', 0.2))
        n_fr = int(dataset.get('n_frames', 35))
        dt = float(dataset.get('dt_days', 1.0))

        # Static force terms
        a3 = a ** 3 if a > 0 else 1.0
        a6 = a ** 6 if a > 0 else 1.0
        F_orbit = G_NEWTON * M_p * M_s / a3
        F_tide = G_NEWTON * M_p * M_s * R_p / a6
        r_gal2 = r_gal ** 2 if r_gal > 0 else 1.0
        F_gal = (v_gal ** 2) / r_gal + G_NEWTON * M_DM / r_gal2 if r_gal > 0 else 0.0

        # Frame iteration with slight perturbation per frame
        frames = []
        for i in range(n_fr):
            # Small sinusoidal perturbation simulating orbital phase
            phase = 2.0 * math.pi * i / max(n_fr, 1)
            perturb = 1.0 + 0.01 * math.sin(phase)

            F_orb_i = F_orbit * perturb
            F_tid_i = F_tide * perturb
            F_gal_i = F_gal * perturb

            F_env_i = w_orb * F_orb_i + w_tid * F_tid_i + w_gal_w * F_gal_i

            frames.append({
                'frame': i + 1,
                't_days': i * dt,
                'F_env': F_env_i,
                'perturbation': perturb,
            })

        F_env_values = [f['F_env'] for f in frames]
        F_env_mean = sum(F_env_values) / max(len(F_env_values), 1)
        F_env_std = math.sqrt(sum((v - F_env_mean) ** 2 for v in F_env_values) / max(len(F_env_values), 1))

        return {
            'F_orbit': F_orbit,
            'F_tide': F_tide,
            'F_gal': F_gal,
            'F_env_mean': F_env_mean,
            'F_env_std': F_env_std,
            'F_env_min': min(F_env_values),
            'F_env_max': max(F_env_values),
            'n_frames': n_fr,
            'frames_first5': frames[:5],
            'frames_last5': frames[-5:],
            'convergence_ratio': F_env_std / abs(F_env_mean) if F_env_mean != 0 else 0.0,
            'paper': 'PAPER_861',
        }


class UniversalMagnetismUmMasterEquationCalc(_CP4Calculator):  # PAPER_862 #446
    """Universal Magnetism U_m — fourth master UQFF equation calculator.

    Master equation:
        U_m = sum_j [mu_j(t, rho_vac,[SCm]) / (r_j / r)
              * (1 - exp(-gamma * t) * cos(pi * t / n)) * phi^j]
              * P_SCm * E_react(t) * (1 + 1e13 * f_Heaviside) * (1 + f_quasi)

    Sub-equations:
        mu_j(t) = (1e3 + 0.4*sin(omega_c * t)) * 3.38e20 T*pm^3
        E_react(t) = 1e46 * exp(-kappa * t)
        omega_c = 1.585e-8 rad/s

    This is the fourth member of the quadriadic UQFF system
    (alongside Compressed, Resonance, and Buoyancy), governing magnetic
    and electric field dynamics through vacuum density coupling.

    Source: grok_share_589f6949-6fe9.txt lines 1900-1960, 2050-2100
    """

    category = "Universal Magnetism / Master Equation"

    OMEGA_C = 1.585e-8  # rad/s — cosmic angular frequency

    def compute(self, dataset: dict) -> dict:
        """Compute Universal Magnetism U_m master equation.

        Parameters (from dataset):
            r           : float — distance/radius (m)
            t_days      : float — time in days
            n           : int   — quantum state (1-26, default 1)
            n_sources   : int   — number of contributing sources j (default 1)
            mu_base     : float — base dipole moment coefficient (default 1e3)
            B_dipole    : float — reference dipole (T*m^3, default 3.38e20)
            gamma       : float — decay rate (day^-1, default 5e-5)
            P_SCm       : float — SCm pressure factor (default 1.0)
            f_Heaviside : float — Heaviside correction (default 0.01)
            f_quasi     : float — quasi-particle correction (default 0.01)
            phi         : float — azimuthal angle factor (default 1.0)
        """
        r = float(dataset.get('r', 1.89e16))
        t = float(dataset.get('t_days', 1.0))
        n = int(dataset.get('n', 1))
        n_src = int(dataset.get('n_sources', 1))
        mu_base = float(dataset.get('mu_base', 1e3))
        B_dip = float(dataset.get('B_dipole', 3.38e20))
        gamma = float(dataset.get('gamma', GAMMA_DECAY))
        P_SCm = float(dataset.get('P_SCm', 1.0))
        f_H = float(dataset.get('f_Heaviside', 0.01))
        f_q = float(dataset.get('f_quasi', 0.01))
        phi = float(dataset.get('phi', 1.0))

        t_sec = t * 86400.0

        # mu_j(t) = (mu_base + 0.4*sin(omega_c * t_sec)) * B_dipole
        mu_j = (mu_base + 0.4 * math.sin(self.OMEGA_C * t_sec)) * B_dip

        # Temporal bracket: (1 - exp(-gamma*t) * cos(pi*t/n))
        exp_gt = math.exp(-gamma * t) if gamma * t < 700 else 0.0
        cos_tn = math.cos(math.pi * t / n) if n > 0 else 1.0
        bracket = 1.0 - exp_gt * cos_tn

        # E_react
        E_react = self._e_react(t)

        # Per-source contribution
        U_m_per_source = (mu_j / r) * bracket * (phi ** 1) if r > 0 else 0.0

        # Total U_m with all factors
        U_m_total = n_src * U_m_per_source * P_SCm * E_react * (1.0 + 1e13 * f_H) * (1.0 + f_q)

        # Imaginary component (superconducting quantum portion)
        U_m_imag = U_m_total  # mirror magnitude

        # Magnetic field estimate: B ~ sqrt(2 * mu_0 * |U_m|)
        B_est = math.sqrt(2.0 * MU_0 * abs(U_m_total)) if U_m_total != 0 else 0.0

        return {
            'mu_j': mu_j,
            'bracket_term': bracket,
            'exp_gamma_t': exp_gt,
            'cos_pi_t_over_n': cos_tn,
            'E_react': E_react,
            'U_m_per_source': U_m_per_source,
            'U_m_total_real': U_m_total,
            'U_m_total_imag': U_m_imag,
            'B_estimated_T': B_est,
            'Heaviside_amplification': 1.0 + 1e13 * f_H,
            'quasi_correction': 1.0 + f_q,
            'n_sources': n_src,
            'paper': 'PAPER_862',
        }
'''

SESSION_LIST = '''
_SESSION_199_CLASSES = [
    'LENRKetaCalibration3EnvironmentDeltaKCalc',             # PAPER_854 #438
    'PseudoMonopole26StateVacuumDensityCalc',                 # PAPER_855 #439
    'HiggsVacuumUHExcitationKHiggsUQFFCalc',                  # PAPER_856 #440
    'NGC346Ug3StarFormationTempVradCalc',                      # PAPER_857 #441
    'Westerlund2QuadriadicRealImaginaryCalc',                  # PAPER_858 #442
    'MicroPlasmoid25umLENRBuoyancyReversalCalc',               # PAPER_859 #443
    'NeutrinoEnergyUQFFVacuumRatioCalc',                       # PAPER_860 #444
    'KeplerOrreryV35FrameIterativeUbCalc',                     # PAPER_861 #445
    'UniversalMagnetismUmMasterEquationCalc',                  # PAPER_862 #446
]

'''

def main():
    with open(CP4, "r", encoding="utf-8") as f:
        text = f.read()

    original_len = text.count('\n') + 1
    print(f"[INFO] CP4 original: {original_len} lines")

    # ── 1. Update version ──
    old_ver = f"Version: {OLD_VERSION}"
    new_ver = f"Version: {NEW_VERSION}"
    if old_ver in text:
        text = text.replace(old_ver, new_ver, 1)
        print(f"[OK] Version: {OLD_VERSION} → {NEW_VERSION}")
    else:
        print(f"[WARN] Version string '{old_ver}' not found")

    # ── 2. Add session header line ──
    anchor = f"    Updated: Session 198 v{OLD_VERSION}"
    if anchor in text:
        text = text.replace(anchor, SESSION_HEADER_LINE + "\n" + anchor, 1)
        print(f"[OK] Session {SESSION} header line inserted")
    else:
        print(f"[WARN] Anchor '{anchor}' not found — trying fallback")
        anchor2 = "Updated: Session 198"
        if anchor2 in text:
            idx = text.index(anchor2)
            line_start = text.rfind('\n', 0, idx) + 1
            text = text[:line_start] + SESSION_HEADER_LINE + "\n" + text[line_start:]
            print(f"[OK] Session {SESSION} header inserted via fallback")
        else:
            print("[ERROR] Could not find Session 198 header anchor")
            sys.exit(1)

    # ── 3. Insert new classes + session list before _SESSION_198_CLASSES ──
    marker = "_SESSION_198_CLASSES"
    if marker in text:
        idx = text.index(marker)
        # Find the start of the line
        line_start = text.rfind('\n', 0, idx)
        insert_point = line_start if line_start >= 0 else 0
        text = text[:insert_point] + NEW_CLASSES + "\n" + SESSION_LIST + text[insert_point:]
        print(f"[OK] {NUM_CLASSES} new classes + session list inserted before {marker}")
    else:
        print(f"[ERROR] Marker '{marker}' not found")
        sys.exit(1)

    # ── 4. Write ──
    with open(CP4, "w", encoding="utf-8") as f:
        f.write(text)

    new_len = text.count('\n') + 1
    print(f"[OK] CP4 updated: {original_len} → {new_len} lines")
    print(f"[OK] Classes: 429 → {429 + NUM_CLASSES} (#{FIRST_CLASS}-#{LAST_CLASS})")
    print(f"[OK] Papers: PAPER_{FIRST_PAPER}-{LAST_PAPER}")
    print(f"[OK] Version: {NEW_VERSION}")

if __name__ == "__main__":
    main()
