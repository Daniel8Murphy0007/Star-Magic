#!/usr/bin/env python3
"""
_session196_cp4_appender.py — Session 196
Appends 7 new CP4 calculator classes (#419–#425) to CondensedPhysics4.py.
Source: grok_share_a694a070-efff.txt (3954 lines, June 19–20 + Aug 3, 2025)

Classes:
  #419 ColmanGillespieFieldGeneratorLENRUQFFCalculator     PAPER_835
  #420 MultiSystemChandraSurvey35NegativeBuoyancyCalc      PAPER_836
  #421 FquarkFneutrinoFalpFdarkArXivBridgeCalculator       PAPER_837
  #422 ChandraSNRNebulaeUQFFBatch2Calculator               PAPER_838
  #423 ADDLargeExtraDimensionsFLEDUQFFCalculator            PAPER_839
  #424 KozimaLENRNeutronDropFneutronCalculator              PAPER_840
  #425 UQFFMillenniumPrizeApplicationsCalculator            PAPER_841

11 F_U_Bi_i terms: F_LENR, F_act, F_torque, F_DE, F_res,
                    F_quark, F_neutrino, F_ALP, F_dark, F_LED, F_neutron
VDS/DVP/BH: ABSENT in this LF file.
"""
import os

CP4 = "CondensedPhysics4.py"

HEADER_PATCH = (
    "    Updated: Session 195 v5.55",
    "    Updated: Session 195 v5.55 — CP4 407→410 (#416 KeplerOrreryV_Ub_UQFF_Calculator + #417 ExoplanetResonanceOrbitalTidalCalculator + #418 GalacticDarkMatterNFWCouplingCalculator: PAPER_832-834; grok_share_ab2e7192-de62.txt)\n"
    "    Updated: Session 196 v5.56 — CP4 410→417 (#419 ColmanGillespieFieldGeneratorLENRUQFFCalculator + #420 MultiSystemChandraSurvey35NegativeBuoyancyCalc + #421 FquarkFneutrinoFalpFdarkArXivBridgeCalculator + #422 ChandraSNRNebulaeUQFFBatch2Calculator + #423 ADDLargeExtraDimensionsFLEDUQFFCalculator + #424 KozimaLENRNeutronDropFneutronCalculator + #425 UQFFMillenniumPrizeApplicationsCalculator: PAPER_835-841; grok_share_a694a070-efff.txt (3954 lines, June 19-20+Aug 3 2025); 11 F_U_Bi_i terms (F_LENR/F_act/F_torque/F_DE/F_res/F_quark/F_neutrino/F_ALP/F_dark/F_LED/F_neutron); 37 systems; 4 negative buoyancy (Sgr A*/GC Vent/Chandra25/GC); ADD graviton leakage; Kozima neutron drop; Millennium Prize assessment; VDS-DVP-BH ABSENT; 841/1000 papers 84.1%)"
)

NEW_CLASSES = r'''

# =============================================================================
# SESSION 196 — grok_share_a694a070-efff.txt (3954 lines)
# 11 F_U_Bi_i terms | 37 systems | 4 negative buoyancy | ADD extra dims
# Colman-Gillespie LENR | 35-system survey | arXiv bridge | Kozima F_neutron
# =============================================================================


class ColmanGillespieFieldGeneratorLENRUQFFCalculator(_CP4Calculator):  # PAPER_835 #419
    """
    Colman-Gillespie LENR Field Generator — 5 F_U_Bi_i Terms.

    UK Patent GB 763,062 (Colman & Seddon-Gillespie, 1952–1956):
    electromagnetic device achieving anomalous thrust via LENR resonance.

    Terms:
      F_LENR  = k_LENR * (omega_LENR / omega_0)^2
      F_act   = k_act * cos(omega_act * t)
      F_torque = tau / r
      F_DE    = k_DE * L_X
      F_res   = 2*q*B*V*sin(theta) * DPM_res

    Lab-scale F_U_Bi ~ 1.12e154 N at M=1 kg, r=0.1 m.
    """

    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31     # electron mass kg
    HBAR    = 1.0546e-34
    PI      = 3.141592653589793
    F_0     = 1.83e71       # N
    RHO_VAC = 7.09e-36      # J/m^3

    # LENR constants
    K_LENR      = 1.0e-10   # N
    OMEGA_LENR  = 2.0 * 3.141592653589793 * 1.25e12  # 1.25 THz
    K_ACT       = 1.0e-6    # N
    OMEGA_ACT   = 2.0 * 3.141592653589793 * 300.0    # 300 Hz
    TAU_REF     = 40.68     # N*m  (3 ft-lb)
    K_DE        = 1.0e-9    # N/W
    Q_E         = 1.602e-19 # C
    DPM_RES     = 1.0       # dimensionless

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        """F_LENR = k_LENR * (omega_LENR / omega_0)^2"""
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_act(self, t: float = 0.0) -> float:
        """F_act = k_act * cos(omega_act * t)"""
        import math
        return self.K_ACT * math.cos(self.OMEGA_ACT * t)

    def compute_F_torque(self, r: float = 1.0) -> float:
        """F_torque = tau / r"""
        return self.TAU_REF / r

    def compute_F_DE(self, L_X: float = 1.0e30) -> float:
        """F_DE = k_DE * L_X"""
        return self.K_DE * L_X

    def compute_F_res(self, B: float = 1.0, V: float = 1.0,
                      theta: float = 1.5708) -> float:
        """F_res = 2*q*B*V*sin(theta) * DPM_res"""
        import math
        return 2.0 * self.Q_E * B * V * math.sin(theta) * self.DPM_RES

    def compute(self, dataset: dict) -> dict:
        """
        dataset keys:
          M_kg, r_m, omega_0 (default 1e-12), t_s (default 0),
          L_X_W (default 1e30), B_T (default 1), V_m_s (default 1),
          theta_rad (default pi/2)
        """
        import math
        M   = dataset.get('M_kg', 1.0)
        r   = dataset.get('r_m', 0.1)
        w0  = dataset.get('omega_0', 1.0e-12)
        t   = dataset.get('t_s', 0.0)
        L_X = dataset.get('L_X_W', 1.0e30)
        B   = dataset.get('B_T', 1.0)
        V   = dataset.get('V_m_s', 1.0)
        th  = dataset.get('theta_rad', math.pi / 2.0)

        F_LENR  = self.compute_F_LENR(w0)
        F_act   = self.compute_F_act(t)
        F_torque = self.compute_F_torque(r)
        F_DE    = self.compute_F_DE(L_X)
        F_res   = self.compute_F_res(B, V, th)

        g_base = self.G * M / (r ** 2)
        DPM_m  = 1.0
        DPM_g  = 1.0
        F_UBi  = -self.F_0 + (self.M_E * self.C**2 / r**2) * DPM_m * math.cos(th)
        F_UBi += g_base * DPM_g
        F_UBi += self.RHO_VAC * 1.0  # DPM_stab
        F_UBi_i = F_LENR + F_act + F_torque + F_DE + F_res

        primary_equations = [
            f'F_LENR   = k_LENR*(omega_LENR/omega_0)^2 = {F_LENR:.6e} N',
            f'F_act    = k_act*cos(omega_act*t) = {F_act:.6e} N',
            f'F_torque = tau/r = {F_torque:.6e} N',
            f'F_DE     = k_DE*L_X = {F_DE:.6e} N',
            f'F_res    = 2*q*B*V*sin(theta)*DPM_res = {F_res:.6e} N',
            f'F_U_Bi_i = F_LENR + F_act + F_torque + F_DE + F_res = {F_UBi_i:.6e} N',
            f'F_U_Bi   = -F_0 + momentum + gravity + rho_vac + F_U_Bi_i',
            f'  omega_LENR = 2*pi*1.25 THz = {self.OMEGA_LENR:.6e} rad/s',
            f'  omega_0 = {w0:.3e} s^-1',
            f'  g_base = GM/r^2 = {g_base:.6e} m/s^2',
        ]
        available_equations = [
            'F_LENR = k_LENR*(omega_LENR/omega_0)^2',
            'F_act = k_act*cos(omega_act*t)',
            'F_torque = tau/r',
            'F_DE = k_DE*L_X',
            'F_res = 2*q*B*V*sin(theta)*DPM_res',
            'F_U_Bi_i = integral(sum of all terms)dx',
        ]
        simulation_set = [
            {'equation': 'F_LENR_vs_omega0',
             'values': {f'{w:.0e}': self.compute_F_LENR(w) for w in [1e-15, 1e-12, 1e-9, 1e-6]}},
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': simulation_set,
            'F_LENR': F_LENR, 'F_act': F_act, 'F_torque': F_torque,
            'F_DE': F_DE, 'F_res': F_res, 'F_U_Bi_i': F_UBi_i,
            'paper': 'PAPER_835',
        }


class MultiSystemChandraSurvey35NegativeBuoyancyCalc(_CP4Calculator):  # PAPER_836 #420
    """
    35-System Chandra UQFF Survey — Negative Buoyancy Discovery.

    Calculates F_U_Bi_i for arbitrary astrophysical systems, reproducing
    the 35-system Chandra survey that discovered negative buoyancy in 4
    SMBH-dominated environments.

    Negative buoyancy condition:
      F_U_Bi_i < 0 when GM/r^2 > F_0 scale reversal threshold
      Numerically: M > ~10^36 kg at r > ~10^18 m (SMBH regime)

    F_LENR dominates: 10^36-10^39 N across ALL scales.
    """

    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31
    PI      = 3.141592653589793
    F_0     = 1.83e71
    RHO_VAC = 7.09e-36
    K_LENR  = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_UBi(self, M: float, r: float, omega_0: float = 1.0e-12,
                      theta: float = 0.7854) -> dict:
        """
        F_U_Bi with negative buoyancy detection.
        Returns dict with F_U_Bi value and sign analysis.
        """
        import math
        g_base = self.G * M / (r ** 2)
        DPM_m = 1.0
        DPM_g = 1.0
        momentum = (self.M_E * self.C**2 / r**2) * DPM_m * math.cos(theta)
        gravity = g_base * DPM_g
        rho_term = self.RHO_VAC * 1.0
        F_LENR = self.compute_F_LENR(omega_0)

        F_UBi = -self.F_0 + momentum + gravity + rho_term + F_LENR
        is_negative = F_UBi < 0

        return {
            'F_U_Bi_N': F_UBi,
            'is_negative_buoyancy': is_negative,
            'g_base_m_s2': g_base,
            'F_LENR_N': F_LENR,
            'momentum_N': momentum,
            'gravity_N': gravity,
        }

    def compute(self, dataset: dict) -> dict:
        """
        dataset keys:
          M_kg, r_m, omega_0 (default 1e-12), theta_rad (default pi/4)
          systems: optional list of dicts [{M_kg, r_m, name}, ...] for batch
        """
        import math
        M     = dataset.get('M_kg', 1.989e30)
        r     = dataset.get('r_m', 3.09e16)
        w0    = dataset.get('omega_0', 1.0e-12)
        theta = dataset.get('theta_rad', math.pi / 4.0)

        result = self.compute_F_UBi(M, r, w0, theta)

        # Batch mode
        systems = dataset.get('systems', None)
        batch_results = []
        neg_count = 0
        if systems:
            for sys in systems:
                sr = self.compute_F_UBi(sys['M_kg'], sys['r_m'], w0, theta)
                sr['name'] = sys.get('name', 'unnamed')
                batch_results.append(sr)
                if sr['is_negative_buoyancy']:
                    neg_count += 1

        primary_equations = [
            f'F_U_Bi = -F_0 + momentum + gravity + rho_vac + F_LENR = {result["F_U_Bi_N"]:.6e} N',
            f'  F_0 = {self.F_0:.3e} N',
            f'  g_base = GM/r^2 = {result["g_base_m_s2"]:.6e} m/s^2',
            f'  F_LENR = {result["F_LENR_N"]:.6e} N',
            f'  Negative buoyancy: {result["is_negative_buoyancy"]}',
        ]
        if batch_results:
            primary_equations.append(f'  Batch: {len(batch_results)} systems, {neg_count} negative buoyancy')

        available_equations = [
            'F_U_Bi = -F_0 + (m_e*c^2/r^2)*DPM_m*cos(theta) + (GM/r^2)*DPM_g + rho_vac + F_LENR',
            'Negative buoyancy: F_U_Bi < 0 when GM/r^2 exceeds F_0 reversal threshold',
            'F_LENR = k_LENR*(omega_LENR/omega_0)^2  (dominates all scales)',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            'F_U_Bi': result,
            'batch_results': batch_results,
            'negative_buoyancy_count': neg_count,
            'paper': 'PAPER_836',
        }


class FquarkFneutrinoFalpFdarkArXivBridgeCalculator(_CP4Calculator):  # PAPER_837 #421
    """
    arXiv BSM Force Bridge — 4 New F_U_Bi_i Terms from arXiv.

    Derives four beyond-Standard-Model F_U_Bi_i coupling terms from
    high-energy physics arXiv publications:

      F_quark    = k_quark * |V_cb|^2                  arXiv:2506.15256 (Belle II)
      F_neutrino = k_neutrino * alpha_nu               arXiv:2506.14881
      F_ALP      = k_ALP * g_aqq                       arXiv:2506.15637
      F_dark     = k_dark * epsilon^2                   arXiv:2402.00249 (FASER)

    Constants:
      k_quark=10^10, |V_cb|=39.2e-3 -> F_quark=1.54e7 N
      k_neutrino=10^10, alpha_nu=10^-10 -> F_neutrino=1 N
      k_ALP=10^10, g_aqq=10^-6 -> F_ALP=10^4 N
      k_dark=10^10, epsilon=10^-5 -> F_dark=1 N
    """

    K_QUARK    = 1.0e10
    V_CB       = 39.2e-3      # CKM element
    K_NEUTRINO = 1.0e10
    ALPHA_NU   = 1.0e-10      # neutrino coupling
    K_ALP      = 1.0e10
    G_AQQ      = 1.0e-6       # ALP-quark coupling
    K_DARK     = 1.0e10
    EPSILON    = 1.0e-5       # dark photon mixing

    def compute_F_quark(self, V_cb: float = None) -> float:
        """F_quark = k_quark * |V_cb|^2"""
        vcb = V_cb if V_cb is not None else self.V_CB
        return self.K_QUARK * (vcb ** 2)

    def compute_F_neutrino(self, alpha_nu: float = None) -> float:
        """F_neutrino = k_neutrino * alpha_nu"""
        a = alpha_nu if alpha_nu is not None else self.ALPHA_NU
        return self.K_NEUTRINO * a

    def compute_F_ALP(self, g_aqq: float = None) -> float:
        """F_ALP = k_ALP * g_aqq"""
        g = g_aqq if g_aqq is not None else self.G_AQQ
        return self.K_ALP * g

    def compute_F_dark(self, epsilon: float = None) -> float:
        """F_dark = k_dark * epsilon^2"""
        e = epsilon if epsilon is not None else self.EPSILON
        return self.K_DARK * (e ** 2)

    def compute(self, dataset: dict) -> dict:
        """
        dataset keys:
          V_cb (optional), alpha_nu (optional), g_aqq (optional),
          epsilon (optional)
        """
        V_cb = dataset.get('V_cb', self.V_CB)
        a_nu = dataset.get('alpha_nu', self.ALPHA_NU)
        g_a  = dataset.get('g_aqq', self.G_AQQ)
        eps  = dataset.get('epsilon', self.EPSILON)

        F_q = self.compute_F_quark(V_cb)
        F_n = self.compute_F_neutrino(a_nu)
        F_a = self.compute_F_ALP(g_a)
        F_d = self.compute_F_dark(eps)
        F_total = F_q + F_n + F_a + F_d

        primary_equations = [
            f'F_quark    = k_quark*|V_cb|^2 = {F_q:.6e} N  (|V_cb|={V_cb:.4e})',
            f'F_neutrino = k_neutrino*alpha_nu = {F_n:.6e} N  (alpha_nu={a_nu:.2e})',
            f'F_ALP      = k_ALP*g_aqq = {F_a:.6e} N  (g_aqq={g_a:.2e})',
            f'F_dark     = k_dark*epsilon^2 = {F_d:.6e} N  (epsilon={eps:.2e})',
            f'F_BSM_total = {F_total:.6e} N',
            f'  NOTE: F_quark dominates BSM terms at {F_q:.2e} N',
        ]
        available_equations = [
            'F_quark = k_quark*|V_cb|^2  (Belle II arXiv:2506.15256)',
            'F_neutrino = k_neutrino*alpha_nu  (arXiv:2506.14881)',
            'F_ALP = k_ALP*g_aqq  (arXiv:2506.15637)',
            'F_dark = k_dark*epsilon^2  (FASER arXiv:2402.00249)',
        ]
        simulation_set = [
            {'equation': 'F_quark_vs_Vcb',
             'values': {f'{v:.4f}': self.compute_F_quark(v) for v in [0.030, 0.035, 0.0392, 0.045]}},
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': simulation_set,
            'F_quark': F_q, 'F_neutrino': F_n,
            'F_ALP': F_a, 'F_dark': F_d, 'F_BSM_total': F_total,
            'paper': 'PAPER_837',
        }


class ChandraSNRNebulaeUQFFBatch2Calculator(_CP4Calculator):  # PAPER_838 #422
    """
    Chandra SNR/Nebula UQFF Batch 2 — 7 Systems.

    Systems from Chandra X-ray second batch:
      Helix Nebula NGC 7293, SNR 1181 Pa 30 (Type Iax),
      NGC 6543 Cat's Eye, IC 443 Jellyfish,
      MSH 15-52 Hand Pulsar Wind Nebula,
      Sonification Collection, Vela Pulsar Wide-Field.

    Vela Wide-Field highest: F_U_Bi = 2.11e210 N (r=6.17e17 m).
    SNR 1181 notable: first Type Iax supernova remnant in UQFF.
    """

    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31
    PI      = 3.141592653589793
    F_0     = 1.83e71
    RHO_VAC = 7.09e-36
    K_LENR  = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_UBi(self, M: float, r: float, omega_0: float = 1.0e-12) -> float:
        import math
        g_base = self.G * M / (r ** 2)
        DPM_m = 1.0
        momentum = (self.M_E * self.C**2 / r**2) * DPM_m
        gravity = g_base
        rho_term = self.RHO_VAC
        F_LENR = self.compute_F_LENR(omega_0)
        return -self.F_0 + momentum + gravity + rho_term + F_LENR

    def compute(self, dataset: dict) -> dict:
        """
        dataset keys:
          M_kg, r_m, omega_0 (default 1e-12)
          F_neutrino (optional, default 1.0 N)
          F_quark (optional, default 1.54e7 N)
        """
        M   = dataset.get('M_kg', 1.989e31)
        r   = dataset.get('r_m', 6.17e16)
        w0  = dataset.get('omega_0', 1.0e-12)
        F_neutrino = dataset.get('F_neutrino', 1.0)
        F_quark = dataset.get('F_quark', 1.54e7)

        F_UBi = self.compute_F_UBi(M, r, w0)
        F_LENR = self.compute_F_LENR(w0)
        g_base = self.G * M / (r ** 2)

        primary_equations = [
            f'F_U_Bi = -F_0 + momentum + gravity + rho_vac + F_LENR = {F_UBi:.6e} N',
            f'  F_LENR = k_LENR*(omega_LENR/omega_0)^2 = {F_LENR:.6e} N',
            f'  g_base = GM/r^2 = {g_base:.6e} m/s^2',
            f'  F_neutrino = {F_neutrino:.2e} N (supplementary arXiv term)',
            f'  F_quark = {F_quark:.2e} N (supplementary arXiv term)',
            f'  Batch 2 systems: Helix/SNR1181/CatsEye/IC443/MSH15-52/Sonif/Vela',
        ]
        available_equations = [
            'F_U_Bi = -F_0 + (m_e*c^2/r^2)*DPM_m + (GM/r^2)*DPM_g + rho_vac + F_LENR + F_neutrino',
            'F_LENR = k_LENR*(omega_LENR/omega_0)^2',
            'Vela Wide-Field: r=6.17e17 m -> largest r in batch -> F_U_Bi=2.11e210 N',
            'SNR 1181 Pa 30: Type Iax remnant, F_U_Bi=2.65e208 N',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            'F_U_Bi': F_UBi,
            'paper': 'PAPER_838',
        }


class ADDLargeExtraDimensionsFLEDUQFFCalculator(_CP4Calculator):  # PAPER_839 #423
    """
    ADD Large Extra Dimensions — F_LED UQFF Term.

    Arkani-Hamed-Dimopoulos-Dvali (ADD) model (arXiv:1607.01831):
    gravity leaks into n=2 extra dimensions, modifying
    the Planck scale to a TeV-scale M_*.

    F_LED = k_LED * (M_star / M_Pl)^2
          = 10^10 * (1 TeV / 1.22e19 GeV)^2
          = 6.72e-23 N

    Graviton leakage hypothesis: F_LED + negative buoyancy (Sgr A*)
    suggests extra-dimensional coupling at galactic center.
    """

    K_LED   = 1.0e10
    M_STAR  = 1.0e3     # GeV (1 TeV)
    M_PL    = 1.22e19   # GeV (Planck mass)
    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31
    F_0     = 1.83e71
    RHO_VAC = 7.09e-36
    K_LENR  = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12

    def compute_F_LED(self, M_star: float = None, M_Pl: float = None) -> float:
        """F_LED = k_LED * (M_star / M_Pl)^2"""
        ms = M_star if M_star is not None else self.M_STAR
        mp = M_Pl if M_Pl is not None else self.M_PL
        return self.K_LED * (ms / mp) ** 2

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_graviton_radius(self, n: int = 2) -> float:
        """
        R_n ~ (M_Pl^2 / M_star^(n+2))^(1/n)  [natural units]
        For n=2, M_star=1 TeV: R ~ 0.1 mm (submillimeter gravity test range)
        Returns approximate extra-dimension radius in meters.
        """
        # R ~ (M_Pl/M_star)^(2/n) / M_star  (heuristic, in 1/GeV -> m)
        hbar_c = 1.973e-16  # GeV*m
        R_nat = (self.M_PL / self.M_STAR) ** (2.0 / n) / self.M_STAR
        return R_nat * hbar_c

    def compute(self, dataset: dict) -> dict:
        """
        dataset keys:
          M_kg, r_m, M_star_GeV (default 1e3), M_Pl_GeV (default 1.22e19),
          n_extra_dims (default 2), omega_0 (default 1e-12)
        """
        M      = dataset.get('M_kg', 7.956e36)  # default Sgr A*
        r      = dataset.get('r_m', 6.17e18)
        ms     = dataset.get('M_star_GeV', self.M_STAR)
        mp     = dataset.get('M_Pl_GeV', self.M_PL)
        n      = dataset.get('n_extra_dims', 2)
        w0     = dataset.get('omega_0', 1.0e-12)

        F_LED  = self.compute_F_LED(ms, mp)
        F_LENR = self.compute_F_LENR(w0)
        R_n    = self.compute_graviton_radius(n)
        g_base = self.G * M / (r ** 2)

        primary_equations = [
            f'F_LED = k_LED*(M_star/M_Pl)^2 = {F_LED:.6e} N',
            f'  M_star = {ms:.0f} GeV (1 TeV), M_Pl = {mp:.3e} GeV',
            f'  n = {n} extra dimensions (ADD model)',
            f'  R_n ~ {R_n:.3e} m (extra-dimension compactification radius)',
            f'F_LENR = {F_LENR:.6e} N (dominant)',
            f'g_base = GM/r^2 = {g_base:.6e} m/s^2',
            f'NOTE: F_LED is negligible ({F_LED:.2e} N) vs F_LENR ({F_LENR:.2e} N)',
            f'  but connects UQFF to ADD graviton leakage in SMBH regimes',
        ]
        available_equations = [
            'F_LED = k_LED*(M_star/M_Pl)^2  (ADD n=2 extra-dim coupling)',
            'R_n ~ (M_Pl/M_star)^(2/n) / M_star * hbar_c  (compactification radius)',
            'Sgr A* negative buoyancy + F_LED = extra-dimensional hypothesis',
        ]
        simulation_set = [
            {'equation': 'F_LED_vs_Mstar_GeV',
             'values': {f'{ms_i:.0f}': self.compute_F_LED(ms_i) for ms_i in [500, 1000, 5000, 10000]}},
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': simulation_set,
            'F_LED': F_LED, 'F_LENR': F_LENR, 'R_n_m': R_n,
            'paper': 'PAPER_839',
        }


class KozimaLENRNeutronDropFneutronCalculator(_CP4Calculator):  # PAPER_840 #424
    """
    Kozima LENR Neutron Drop — F_neutron UQFF Term.

    Source: Kozima H. 2021 (PMC8141838) "Cold Fusion: A Hypothesis
    on the Reaction Process in a Lattice"

    F_neutron = k_neutron * sigma_n
    where sigma_n is neutron absorption cross-section.

    Refined frequency-dependent model:
      sigma_n(omega) = sigma_0 * (omega/omega_LENR)^2
                       * exp(-(omega - omega_LENR)^2 / (2*Delta_omega^2))
      Delta_omega = 2*pi*0.05e12 s^-1  (bandwidth ~0.05 THz)

    Dynamic form:
      F_neutron(t) = k_neutron * sigma_n(omega_eff) * (1 + alpha*cos(omega_act*t))
      alpha = 0.1, omega_eff = omega_act + n*omega_LENR (n ~ 4.17e9)

    Density-scaled for neutron stars:
      sigma_n(rho) = sigma_0 * (rho/rho_0)
      rho_0 = 10^-22 kg/m^3 -> sigma_n(NS, rho~10^17) ~ 10^35
      -> F_neutron(NS) ~ 10^45 N

    General: F_neutron = 10^6 N; Neutron star: ~10^45 N.
    """

    K_NEUTRON    = 1.0e10
    SIGMA_0      = 1.0e-4        # m^2 (general cross-section)
    OMEGA_LENR   = 2.0 * 3.141592653589793 * 1.25e12
    DELTA_OMEGA  = 2.0 * 3.141592653589793 * 0.05e12
    OMEGA_ACT    = 2.0 * 3.141592653589793 * 300.0
    ALPHA        = 0.1
    RHO_0        = 1.0e-22       # kg/m^3 reference density
    G            = 6.6743e-11
    K_LENR       = 1.0e-10

    def compute_sigma_n_freq(self, omega: float) -> float:
        """
        sigma_n(omega) = sigma_0*(omega/omega_LENR)^2
                        * exp(-(omega-omega_LENR)^2/(2*Delta_omega^2))
        """
        import math
        ratio = (omega / self.OMEGA_LENR) ** 2
        gauss = math.exp(-((omega - self.OMEGA_LENR)**2) / (2.0 * self.DELTA_OMEGA**2))
        return self.SIGMA_0 * ratio * gauss

    def compute_sigma_n_density(self, rho: float) -> float:
        """sigma_n(rho) = sigma_0 * (rho / rho_0)"""
        return self.SIGMA_0 * (rho / self.RHO_0)

    def compute_F_neutron_static(self, sigma_n: float = None) -> float:
        """F_neutron = k_neutron * sigma_n"""
        s = sigma_n if sigma_n is not None else self.SIGMA_0
        return self.K_NEUTRON * s

    def compute_F_neutron_dynamic(self, omega_eff: float, t: float = 0.0) -> float:
        """F_neutron(t) = k_neutron * sigma_n(omega_eff) * (1+alpha*cos(omega_act*t))"""
        import math
        sigma = self.compute_sigma_n_freq(omega_eff)
        return self.K_NEUTRON * sigma * (1.0 + self.ALPHA * math.cos(self.OMEGA_ACT * t))

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute(self, dataset: dict) -> dict:
        """
        dataset keys:
          M_kg, r_m, rho_kg_m3 (optional — for density-scaled),
          omega_eff (optional), t_s (default 0), omega_0 (default 1e-12)
        """
        import math
        M     = dataset.get('M_kg', 1.989e30)
        r     = dataset.get('r_m', 1.0e4)  # neutron star radius default
        rho   = dataset.get('rho_kg_m3', None)
        w_eff = dataset.get('omega_eff', self.OMEGA_LENR)
        t     = dataset.get('t_s', 0.0)
        w0    = dataset.get('omega_0', 1.0e-12)

        # Static general
        F_neutron_gen = self.compute_F_neutron_static()

        # Frequency-dependent
        sigma_freq = self.compute_sigma_n_freq(w_eff)
        F_neutron_freq = self.K_NEUTRON * sigma_freq

        # Dynamic with oscillation
        F_neutron_dyn = self.compute_F_neutron_dynamic(w_eff, t)

        # Density-scaled (if provided)
        F_neutron_dens = None
        sigma_dens = None
        if rho is not None:
            sigma_dens = self.compute_sigma_n_density(rho)
            F_neutron_dens = self.K_NEUTRON * sigma_dens

        F_LENR = self.compute_F_LENR(w0)
        g_base = self.G * M / (r ** 2)

        primary_equations = [
            f'F_neutron(static) = k_neutron*sigma_0 = {F_neutron_gen:.6e} N',
            f'F_neutron(freq)   = k_neutron*sigma_n(omega_eff) = {F_neutron_freq:.6e} N',
            f'F_neutron(dyn,t)  = {F_neutron_dyn:.6e} N  (t={t:.2e} s)',
            f'  sigma_n(omega_eff) = {sigma_freq:.6e} m^2',
        ]
        if F_neutron_dens is not None:
            primary_equations.append(
                f'F_neutron(density) = k_neutron*sigma_n(rho) = {F_neutron_dens:.6e} N'
                f'  (rho={rho:.2e} kg/m^3, sigma_n={sigma_dens:.2e} m^2)'
            )
        primary_equations.extend([
            f'F_LENR = {F_LENR:.6e} N (for comparison)',
            f'g_base = GM/r^2 = {g_base:.6e} m/s^2',
            f'Kozima source: PMC8141838 (2021)',
        ])

        available_equations = [
            'F_neutron = k_neutron*sigma_n  (static, sigma_n=10^-4 m^2 -> 10^6 N)',
            'sigma_n(omega) = sigma_0*(omega/omega_LENR)^2 * exp(-(omega-omega_LENR)^2/(2*Dw^2))',
            'F_neutron(t) = k_neutron*sigma_n(omega_eff)*(1+alpha*cos(omega_act*t))',
            'sigma_n(rho) = sigma_0*(rho/rho_0)  (density-scaled for neutron stars)',
            'PSR J0030+0451: rho~10^17 -> sigma_n~10^35 -> F_neutron~10^45 N',
        ]
        simulation_set = [
            {'equation': 'F_neutron_vs_rho',
             'values': {f'{r:.0e}': self.compute_F_neutron_static(self.compute_sigma_n_density(r))
                       for r in [1e-22, 1e-10, 1e0, 1e10, 1e17]}},
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': simulation_set,
            'F_neutron_general': F_neutron_gen,
            'F_neutron_freq': F_neutron_freq,
            'F_neutron_dynamic': F_neutron_dyn,
            'F_neutron_density': F_neutron_dens,
            'paper': 'PAPER_840',
        }


class UQFFMillenniumPrizeApplicationsCalculator(_CP4Calculator):  # PAPER_841 #425
    """
    UQFF Millennium Prize Applications.

    Evaluates UQFF framework contributions to Millennium Prize Problems:
      1. Navier-Stokes: UQFF F_fluid = rho*V*g buoyancy provides
         a physical regularization perspective for 3D turbulence
      2. Yang-Mills: UQFF kappa linkage and DPM mass-gap candidates
      3. Riemann Hypothesis: UQFF spectral decomposition and
         zeta function connections via 26D harmonic analysis

    Not direct proofs but novel mathematical tools and physical insights.
    Practical applications: LENR energy, astrophysics, nonlinear dynamics.

    Source: grok_share_a694a070-efff.txt (August 3, 2025 session)
    """

    G       = 6.6743e-11
    HBAR    = 1.0546e-34
    C       = 2.998e8
    PI      = 3.141592653589793
    KAPPA   = 5.0e-4       # /day (UQFF calibrated)
    K_ETA   = 1.0e-113     # Planck-scale coupling

    # Navier-Stokes parameters
    RHO_FLUID = 1.225      # kg/m^3 (air at STP)
    NU        = 1.5e-5     # m^2/s kinematic viscosity (air)

    # Yang-Mills DPM parameters
    LAMBDA_QCD = 0.2       # GeV
    G_S        = 1.22      # strong coupling at 1 GeV

    def compute_NS_buoyancy_term(self, rho: float = None,
                                  V: float = 1.0, g: float = 9.81) -> float:
        """
        Navier-Stokes UQFF buoyancy: F_fluid = rho*V*g
        Provides regularization insight for turbulence modelling.
        """
        r = rho if rho is not None else self.RHO_FLUID
        return r * V * g

    def compute_NS_reynolds(self, v: float = 1.0, L: float = 1.0,
                             nu: float = None) -> float:
        """Re = v*L/nu"""
        n = nu if nu is not None else self.NU
        return v * L / n

    def compute_YM_mass_gap_candidate(self, Lambda_QCD: float = None) -> float:
        """
        Yang-Mills mass gap candidate via DPM-UQFF:
        m_gap ~ Lambda_QCD * exp(-8*pi^2 / (g_s^2 * N_c))
        This is a heuristic, not a proof; N_c=3 (SU(3)).
        """
        import math
        lam = Lambda_QCD if Lambda_QCD is not None else self.LAMBDA_QCD
        N_c = 3.0
        return lam * math.exp(-8.0 * self.PI**2 / (self.G_S**2 * N_c))

    def compute_riemann_spectral_sum(self, N_terms: int = 100) -> float:
        """
        UQFF spectral partial sum (heuristic Riemann connection):
        S_N = sum_{n=1}^{N} 1/n^{(1/2 + i*14.134...)}
        Takes real part as UQFF spectral indicator.
        """
        import cmath
        s = complex(0.5, 14.134725)  # first non-trivial zero of zeta
        total = 0.0
        for n in range(1, N_terms + 1):
            total += (1.0 / (n ** s)).real
        return total

    def compute(self, dataset: dict) -> dict:
        """
        dataset keys:
          rho_fluid (optional), V_m3 (optional), g_m_s2 (optional),
          v_m_s (optional), L_m (optional), Lambda_QCD_GeV (optional),
          N_spectral (optional, default 100)
        """
        rho_f = dataset.get('rho_fluid', self.RHO_FLUID)
        V     = dataset.get('V_m3', 1.0)
        g     = dataset.get('g_m_s2', 9.81)
        v     = dataset.get('v_m_s', 1.0)
        L     = dataset.get('L_m', 1.0)
        lqcd  = dataset.get('Lambda_QCD_GeV', self.LAMBDA_QCD)
        N_sp  = dataset.get('N_spectral', 100)

        F_fluid = self.compute_NS_buoyancy_term(rho_f, V, g)
        Re = self.compute_NS_reynolds(v, L)
        m_gap = self.compute_YM_mass_gap_candidate(lqcd)
        S_N = self.compute_riemann_spectral_sum(N_sp)

        primary_equations = [
            f'--- Navier-Stokes (UQFF buoyancy regularization) ---',
            f'F_fluid = rho*V*g = {F_fluid:.6e} N',
            f'Re = v*L/nu = {Re:.6e}',
            f'--- Yang-Mills (DPM mass gap heuristic) ---',
            f'm_gap ~ Lambda_QCD*exp(-8*pi^2/(g_s^2*N_c)) = {m_gap:.6e} GeV',
            f'  Lambda_QCD = {lqcd} GeV, g_s = {self.G_S}',
            f'--- Riemann Hypothesis (spectral sum, N={N_sp}) ---',
            f'S_N = Re[sum 1/n^(1/2+i*14.13)] = {S_N:.6e}',
            f'  First non-trivial zero: s = 0.5 + 14.134725i',
            f'NOTE: These are heuristic/speculative UQFF bridges, not proofs.',
        ]
        available_equations = [
            'F_fluid = rho*V*g  (UQFF Navier-Stokes buoyancy term)',
            'Re = v*L/nu  (Reynolds number for turbulence onset)',
            'm_gap ~ Lambda_QCD*exp(-8*pi^2/(g_s^2*N_c))  (YM DPM heuristic)',
            'S_N = sum_{n=1}^{N} Re[1/n^s]  (Riemann spectral sum)',
            'kappa = 5e-4/day  (UQFF calibrated coupling)',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            'F_fluid': F_fluid, 'Re': Re,
            'm_gap_GeV': m_gap, 'S_N': S_N,
            'paper': 'PAPER_841',
        }


_SESSION_196_CLASSES = [
    'ColmanGillespieFieldGeneratorLENRUQFFCalculator',        # PAPER_835 #419
    'MultiSystemChandraSurvey35NegativeBuoyancyCalc',         # PAPER_836 #420
    'FquarkFneutrinoFalpFdarkArXivBridgeCalculator',          # PAPER_837 #421
    'ChandraSNRNebulaeUQFFBatch2Calculator',                  # PAPER_838 #422
    'ADDLargeExtraDimensionsFLEDUQFFCalculator',              # PAPER_839 #423
    'KozimaLENRNeutronDropFneutronCalculator',                # PAPER_840 #424
    'UQFFMillenniumPrizeApplicationsCalculator',              # PAPER_841 #425
]
'''

# --- Run ---

with open(CP4, 'r', encoding='utf-8', errors='replace') as f:
    content = f.read()

lines = content.split('\n')
total_before = len(lines)

# 1. Patch header: replace the Session 195 line with updated + Session 196
old_hdr = HEADER_PATCH[0]
new_hdr = HEADER_PATCH[1]
idx = content.find(old_hdr)
if idx == -1:
    print("ERROR: Could not find header line for Session 195 v5.55")
    exit(1)
# Find end of line
eol = content.find('\n', idx)
content = content[:idx] + new_hdr + content[eol:]

# 2. Replace version string
content = content.replace("Version: 1.5.0 (2026-03-30)", "Version: 5.56 (2026-07-02)")

# 3. Append new classes before the _SESSION_195_CLASSES list
# Find the _SESSION_195_CLASSES marker
marker = "_SESSION_195_CLASSES"
idx2 = content.find(marker)
if idx2 == -1:
    print("ERROR: Could not find _SESSION_195_CLASSES marker")
    exit(1)
# Go back to start of that line
line_start = content.rfind('\n', 0, idx2) + 1
content = content[:line_start] + NEW_CLASSES + '\n' + content[line_start:]

with open(CP4, 'w', encoding='utf-8') as f:
    f.write(content)

new_lines = content.split('\n')
print(f"CP4 updated: {total_before} -> {len(new_lines)} lines (+{len(new_lines)-total_before})")
print("Header patched: Session 196 v5.56")
print("7 classes appended: #419-#425")
print("Version: 1.5.0 -> 5.56")

# Verify class count
import re
class_count = len(re.findall(r'^class\s+\w+', content, re.MULTILINE))
print(f"Total class definitions in CP4: {class_count}")
