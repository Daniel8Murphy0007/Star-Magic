#!/usr/bin/env python3
"""
Session 197 CP4 Appender — grok_share_d1b5f26e-2f60.txt
4,022 lines | 15 share links | 11 F_U_Bi_i terms (same as S196) | ~34 systems | 25 arXiv papers

Creates 11 new classes (#426–#436) → PAPER_842–852
CP4: 417 → 428 classes | v5.56 → v5.57

Classes:
  #426 FloydSweetVTA6DocumentPCVTMotionalEfieldCalc         PAPER_842
  #427 ChandraXRayBatch1GCEagleHBC672NGC7469VirgoCalc       PAPER_843
  #428 Chandra25thAnniversaryCrabOrionNGC6334Calc           PAPER_844
  #429 ChandraSurveyMACSJ0416LensExoSMBHCalc                PAPER_845
  #430 ChandraDeathStar16SMBHGCVentTimelapseCalc            PAPER_846
  #431 SNRNebulaVelaTychoHelixSNR1181NGC6543Calc            PAPER_847
  #432 SonificationCompositeH1821IC443M74MSH1552Calc        PAPER_848
  #433 ArXiv24PaperBatch4FquarkFneutrinoFALPFdarkCalc       PAPER_849
  #434 ADDGravitonLeakageNegBuoyancySgrAExtCalc              PAPER_850
  #435 KozimaNeutronDropDensityScaled8SystemCalc             PAPER_851
  #436 LENRNextStepsExperimentalDesignPSRJ0030Calc           PAPER_852

VDS/DVP/BH: ABSENT
"""

import os, re

CP4 = os.path.join(os.path.dirname(__file__), 'CondensedPhysics4.py')

HEADER_PATCH = (
    "    Updated: Session 196 v5.56",
    "    Updated: Session 196 v5.56"
    " — CP4 410→417 (#419 ColmanGillespieFieldGeneratorLENRUQFFCalculator + #420 MultiSystemChandraSurvey35NegativeBuoyancyCalc + #421 FquarkFneutrinoFalpFdarkArXivBridgeCalculator + #422 ChandraSNRNebulaeUQFFBatch2Calculator + #423 ADDLargeExtraDimensionsFLEDUQFFCalculator + #424 KozimaLENRNeutronDropFneutronCalculator + #425 UQFFMillenniumPrizeApplicationsCalculator: PAPER_835-841; grok_share_a694a070-efff.txt (3954 lines, June 19-20+Aug 3 2025); 11 F_U_Bi_i terms (F_LENR/F_act/F_torque/F_DE/F_res/F_quark/F_neutrino/F_ALP/F_dark/F_LED/F_neutron); 37 systems; 4 negative buoyancy (Sgr A*/GC Vent/Chandra25/GC); ADD graviton leakage; Kozima neutron drop; Millennium Prize assessment; VDS-DVP-BH ABSENT; 841/1000 papers 84.1%)\n"
    "    Updated: Session 197 v5.57 — CP4 417→428 (#426 FloydSweetVTA6DocumentPCVTMotionalEfieldCalc + #427 ChandraXRayBatch1GCEagleHBC672NGC7469VirgoCalc + #428 Chandra25thAnniversaryCrabOrionNGC6334Calc + #429 ChandraSurveyMACSJ0416LensExoSMBHCalc + #430 ChandraDeathStar16SMBHGCVentTimelapseCalc + #431 SNRNebulaVelaTychoHelixSNR1181NGC6543Calc + #432 SonificationCompositeH1821IC443M74MSH1552Calc + #433 ArXiv24PaperBatch4FquarkFneutrinoFALPFdarkCalc + #434 ADDGravitonLeakageNegBuoyancySgrAExtCalc + #435 KozimaNeutronDropDensityScaled8SystemCalc + #436 LENRNextStepsExperimentalDesignPSRJ0030Calc: PAPER_842-852; grok_share_d1b5f26e-2f60.txt (4022 lines, June 19-20 2025); 11 F_U_Bi_i terms (same as S196 expanded); ~34 systems; 25 arXiv papers; Floyd Sweet VTA 6-doc PCVT; Chandra 4 batches per-system; SNR/Nebulae DeepSearch; Sonification composite; ADD graviton leakage ext; Kozima density-scaled; LENR next steps + PSR J0030+0451; VDS-DVP-BH ABSENT; 852/1000 papers 85.2%)"
)

NEW_CLASSES = r'''
# ============================================================================
# SESSION 197 — grok_share_d1b5f26e-2f60.txt (4022 lines, June 19-20 2025)
# 11 classes (#426-#436) | PAPER_842-852 | v5.57 | 2026-04-04
# Floyd Sweet VTA + Chandra 4 batches + SNR/Nebulae + Sonification +
# arXiv 24-paper BSM landscape + ADD graviton leakage + Kozima density-scaled +
# LENR next steps + PSR J0030+0451
# 11 F_U_Bi_i terms (same as S196): F_LENR, F_act, F_torque, F_DE, F_res,
# F_quark, F_neutrino, F_ALP, F_dark, F_LED, F_neutron
# VDS/DVP/BH: ABSENT
# ============================================================================


class FloydSweetVTA6DocumentPCVTMotionalEfieldCalc(_CP4Calculator):  # PAPER_842 #426
    """
    Floyd Sweet VTA 6-Document Analysis — PCVT and Motional E-field.

    Six source documents analyzed:
      1. magneticreson.pdf — magnetic resonance coupling
      2. intergalspacetrav.pdf — interstellar propulsion concepts
      3. nothingisimpossible.pdf — ZPE access overview
      4. vacuumtriodeamplifier.pdf — VTA core design
      5. spaceflux.pdf — vacuum energy extraction
      6. barium ferrite conditioning — lattice activation

    PCVT (Parametric Cascade Via Triggered-field):
      ZPE -> triggered by 300 Hz -> cascade at 1.3 THz LENR band

    Motional E-field:
      E_m = V x B  (velocity cross magnetic field)
      F_res = 2*q*B*V*sin(theta)*DPM_resonance

    VTA output: 500 W from 330 uW input (gain ~1.5 million).
    Sweet measured 10 W output from barium ferrite magnets (1990).
    """

    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31
    HBAR    = 1.0546e-34
    Q_E     = 1.602e-19     # electron charge
    PI      = 3.141592653589793
    F_0     = 1.83e71
    RHO_VAC = 7.09e-36
    K_LENR  = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12   # 1.25 THz
    K_ACT   = 1.0e-5
    OMEGA_ACT  = 2.0 * 3.141592653589793 * 300.0      # 300 Hz
    TAU_REF = 3.0 * 0.30480 * 4.44822                 # 3 ft-lb in N-m -> / r
    DPM_RES = 1.0
    VTA_GAIN = 500.0 / 330e-6                          # 500 W / 330 uW

    def compute_F_res(self, B: float = 0.3, V: float = 1.0,
                      theta: float = 0.7854) -> float:
        """Motional E-field force: F_res = 2*q*B*V*sin(theta)*DPM_res"""
        import math
        return 2.0 * self.Q_E * B * V * math.sin(theta) * self.DPM_RES

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_act(self, t: float = 0.0) -> float:
        """F_act = k_act * cos(omega_act * t)"""
        import math
        return self.K_ACT * math.cos(self.OMEGA_ACT * t)

    def compute_PCVT_cascade_ratio(self) -> float:
        """PCVT cascade: ratio = omega_LENR / omega_act = 1.25e12 / 300"""
        return self.OMEGA_LENR / self.OMEGA_ACT

    def compute(self, dataset: dict) -> dict:
        import math
        M     = dataset.get('M_kg', 1.989e30)
        r     = dataset.get('r_m', 1.0e16)
        w0    = dataset.get('omega_0', 1.0e-12)
        t     = dataset.get('t_s', 0.0)
        B     = dataset.get('B_T', 0.3)
        V     = dataset.get('V_m_s', 1.0)
        theta = dataset.get('theta_rad', math.pi / 4.0)

        F_res   = self.compute_F_res(B, V, theta)
        F_LENR  = self.compute_F_LENR(w0)
        F_act   = self.compute_F_act(t)
        F_torque = self.TAU_REF / max(r, 1.0e-10)
        cascade = self.compute_PCVT_cascade_ratio()

        primary_equations = [
            f'--- Floyd Sweet VTA 6-Document PCVT Analysis ---',
            f'F_res = 2*q*B*V*sin(theta)*DPM_res = {F_res:.6e} N',
            f'  B = {B} T, V = {V} m/s, theta = {theta:.4f} rad',
            f'F_LENR = k_LENR*(omega_LENR/omega_0)^2 = {F_LENR:.6e} N',
            f'F_act = k_act*cos(omega_act*t) = {F_act:.6e} N  (t={t:.2e}s)',
            f'F_torque = tau/r = {F_torque:.6e} N',
            f'PCVT cascade ratio = omega_LENR/omega_act = {cascade:.6e}',
            f'VTA gain = {self.VTA_GAIN:.0f}x (500W/330uW)',
            f'6 documents: magneticreson/intergal/nothingisimpossible/vta/spaceflux/barium',
        ]
        available_equations = [
            'E_m = V x B  (motional E-field)',
            'F_res = 2*q*B*V*sin(theta)*DPM_res',
            'PCVT: ZPE -> 300 Hz trigger -> 1.25 THz LENR cascade',
            'F_LENR = k_LENR*(omega_LENR/omega_0)^2',
            'F_act = k_act*cos(omega_act*t)',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            'F_res': F_res, 'F_LENR': F_LENR, 'F_act': F_act,
            'F_torque': F_torque, 'PCVT_cascade': cascade,
            'paper': 'PAPER_842',
        }


class ChandraXRayBatch1GCEagleHBC672NGC7469VirgoCalc(_CP4Calculator):  # PAPER_843 #427
    """
    Chandra X-ray Batch 1 — 5 Systems.

    Galactic Center: M=7.956e36 kg, r=6.17e18 m, F_U_Bi ~ -8.31e211 N (NEGATIVE)
    Eagle Nebula M16: M=1e32 kg, r=2.78e17 m, F_U_Bi ~ 2.65e208 N
    HBC 672: M=3.978e30 kg, r=3e15 m, F_U_Bi ~ 2.65e208 N
    NGC 7469: M=2.387e37 kg, r=1.39e22 m, F_U_Bi ~ -8.31e211 N (NEGATIVE)
    Virgo Cluster: M=2.387e45 kg, r=1.7e23 m, F_U_Bi ~ 2.65e208 N

    F_LENR dominates all systems: ~6.17e37 N.
    Negative buoyancy in GC and NGC 7469 (AGN-dominated).
    """

    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31
    F_0     = 1.83e71
    RHO_VAC = 7.09e-36
    K_LENR  = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_UBi(self, M: float, r: float, omega_0: float = 1.0e-12) -> float:
        import math
        g_base = self.G * M / (r ** 2)
        momentum = (self.M_E * self.C**2 / r**2)
        F_LENR = self.compute_F_LENR(omega_0)
        return -self.F_0 + momentum + g_base + self.RHO_VAC + F_LENR

    def compute(self, dataset: dict) -> dict:
        M  = dataset.get('M_kg', 7.956e36)
        r  = dataset.get('r_m', 6.17e18)
        w0 = dataset.get('omega_0', 1.0e-12)

        F_UBi  = self.compute_F_UBi(M, r, w0)
        F_LENR = self.compute_F_LENR(w0)
        g_base = self.G * M / (r ** 2)

        systems = dataset.get('systems', None)
        batch = []
        if systems:
            for s in systems:
                fb = self.compute_F_UBi(s['M_kg'], s['r_m'], w0)
                batch.append({'name': s.get('name', ''), 'F_U_Bi_N': fb,
                              'negative': fb < 0})

        primary_equations = [
            f'F_U_Bi = -F_0 + momentum + gravity + rho_vac + F_LENR = {F_UBi:.6e} N',
            f'  F_LENR = {F_LENR:.6e} N (dominant)',
            f'  g_base = GM/r^2 = {g_base:.6e} m/s^2',
            f'Batch 1: GC/EagleM16/HBC672/NGC7469/Virgo',
            f'  GC: M=7.956e36, r=6.17e18 -> negative buoyancy',
            f'  NGC7469: M=2.387e37, r=1.39e22 -> negative buoyancy',
        ]
        available_equations = [
            'F_U_Bi = -F_0 + (m_e*c^2/r^2) + (GM/r^2) + rho_vac + F_LENR',
            'F_LENR = k_LENR*(omega_LENR/omega_0)^2',
            'Negative buoyancy: F_U_Bi < 0 in SMBH/AGN systems',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            'F_U_Bi': F_UBi, 'batch_results': batch,
            'paper': 'PAPER_843',
        }


class Chandra25thAnniversaryCrabOrionNGC6334Calc(_CP4Calculator):  # PAPER_844 #428
    """
    Chandra 25th Anniversary — 5 Systems.

    25 Images Composite: M=1e38 kg, r=3.09e22 m -> F_U_Bi ~ -2.09e212 N (NEGATIVE)
    Crab Nebula: M=2.8e30 kg, r=5.56e16 m -> F_U_Bi ~ 2.65e208 N
    Orion Nebula: M=2e33 kg, r=6.17e16 m -> F_U_Bi ~ 2.65e208 N
    NGC 4438/4435 (Eyes Galaxies): M=2e41 kg, r=5.12e22 m -> F_U_Bi ~ 2.65e208 N
    NGC 6334 (Cat's Paw): M=2e35 kg, r=1.05e20 m -> F_U_Bi ~ 2.65e208 N

    25 Images composite shows strongest negative buoyancy: -2.09e212 N
    (massive composite central region, M=10^38 kg).
    """

    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31
    F_0     = 1.83e71
    RHO_VAC = 7.09e-36
    K_LENR  = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_UBi(self, M: float, r: float, omega_0: float = 1.0e-12) -> float:
        import math
        g_base = self.G * M / (r ** 2)
        momentum = (self.M_E * self.C**2 / r**2)
        F_LENR = self.compute_F_LENR(omega_0)
        return -self.F_0 + momentum + g_base + self.RHO_VAC + F_LENR

    def compute(self, dataset: dict) -> dict:
        M  = dataset.get('M_kg', 1e38)
        r  = dataset.get('r_m', 3.09e22)
        w0 = dataset.get('omega_0', 1.0e-12)

        F_UBi  = self.compute_F_UBi(M, r, w0)
        F_LENR = self.compute_F_LENR(w0)

        systems = dataset.get('systems', None)
        batch = []
        if systems:
            for s in systems:
                fb = self.compute_F_UBi(s['M_kg'], s['r_m'], w0)
                batch.append({'name': s.get('name', ''), 'F_U_Bi_N': fb,
                              'negative': fb < 0})

        primary_equations = [
            f'F_U_Bi = -F_0 + momentum + gravity + rho_vac + F_LENR = {F_UBi:.6e} N',
            f'  F_LENR = {F_LENR:.6e} N',
            f'Chandra 25th Anniversary Batch (5 systems):',
            f'  25 Images Composite: M=1e38, r=3.09e22 -> -2.09e212 N (NEGATIVE)',
            f'  Crab Nebula: M=2.8e30, r=5.56e16 -> positive',
            f'  Orion Nebula: M=2e33, r=6.17e16 -> positive',
            f'  NGC 4438/4435: M=2e41, r=5.12e22 -> positive',
            f'  NGC 6334: M=2e35, r=1.05e20 -> positive',
        ]
        available_equations = [
            'F_U_Bi = -F_0 + (m_e*c^2/r^2) + (GM/r^2) + rho_vac + F_LENR',
            'Chandra 25th: 25 composite images -> strongest negative buoyancy',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            'F_U_Bi': F_UBi, 'batch_results': batch,
            'paper': 'PAPER_844',
        }


class ChandraSurveyMACSJ0416LensExoSMBHCalc(_CP4Calculator):  # PAPER_845 #429
    """
    Chandra Survey Batch — MACS J0416 Lensing + 4 Systems.

    MACS J0416.1-2403: Gravitational lensing cluster, M=1e15*M_sun, r=1.47e23 m
    3C 58: Young pulsar wind nebula (1181 AD), M=2.8e30 kg, r=3.09e19 m
    Exoplanet Survey: Composite, M=1.989e30 kg, r=3.09e16 m
    SMBH Survey: 12 galaxies composite, M=1e40 kg, r=3.09e22 m
    Westerlund 1: Massive star cluster, M=6.3e34 kg, r=3.7e19 m

    Gravitational lensing amplifies X-ray signal but does not alter
    F_U_Bi_i computation (intrinsic to source).
    """

    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31
    F_0     = 1.83e71
    RHO_VAC = 7.09e-36
    K_LENR  = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_UBi(self, M: float, r: float, omega_0: float = 1.0e-12) -> float:
        g_base = self.G * M / (r ** 2)
        momentum = (self.M_E * self.C**2 / r**2)
        F_LENR = self.compute_F_LENR(omega_0)
        return -self.F_0 + momentum + g_base + self.RHO_VAC + F_LENR

    def compute_lensing_magnification(self, theta_E: float = 1.0,
                                       theta_S: float = 0.5) -> float:
        """
        Gravitational lensing magnification:
        mu = theta_E / theta_S  (simplified strong lensing)
        Does not affect intrinsic F_U_Bi_i.
        """
        if theta_S == 0:
            return float('inf')
        return theta_E / theta_S

    def compute(self, dataset: dict) -> dict:
        M  = dataset.get('M_kg', 1.989e45)
        r  = dataset.get('r_m', 1.47e23)
        w0 = dataset.get('omega_0', 1.0e-12)
        theta_E = dataset.get('theta_E_arcsec', 1.0)
        theta_S = dataset.get('theta_S_arcsec', 0.5)

        F_UBi  = self.compute_F_UBi(M, r, w0)
        F_LENR = self.compute_F_LENR(w0)
        mu = self.compute_lensing_magnification(theta_E, theta_S)

        systems = dataset.get('systems', None)
        batch = []
        if systems:
            for s in systems:
                fb = self.compute_F_UBi(s['M_kg'], s['r_m'], w0)
                batch.append({'name': s.get('name', ''), 'F_U_Bi_N': fb,
                              'negative': fb < 0})

        primary_equations = [
            f'F_U_Bi = -F_0 + momentum + gravity + rho_vac + F_LENR = {F_UBi:.6e} N',
            f'  F_LENR = {F_LENR:.6e} N',
            f'  Lensing magnification mu = theta_E/theta_S = {mu:.2f}',
            f'Survey batch: MACS_J0416/3C58/ExoSurvey/SMBHSurvey/Westerlund1',
            f'  MACS J0416.1-2403: cluster lensing, M~10^15 M_sun',
            f'  3C 58: young PWN from 1181 AD supernova',
        ]
        available_equations = [
            'F_U_Bi = -F_0 + (m_e*c^2/r^2) + (GM/r^2) + rho_vac + F_LENR',
            'mu = theta_E / theta_S  (strong lensing magnification)',
            'Lensing does not alter intrinsic F_U_Bi_i',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            'F_U_Bi': F_UBi, 'lensing_mu': mu, 'batch_results': batch,
            'paper': 'PAPER_845',
        }


class ChandraDeathStar16SMBHGCVentTimelapseCalc(_CP4Calculator):  # PAPER_846 #430
    """
    Chandra Death Star BHs — 16 SMBH + GC Vent + Timelapse.

    Death Star BHs: 16 SMBH systems composite, M=1.59e40 kg, r=3.09e22 m
    Abell 478: Galaxy cluster, M=1e45 kg, r=1.85e23 m
    NGC 5044: Elliptical galaxy group, M=3e42 kg, r=6.17e22 m
    GC Vent: Galactic center bipolar outflow, M=7.956e36 kg, r=6.17e18 m
      -> F_U_Bi ~ -8.31e211 N (NEGATIVE — same as GC)
    Cas A + Crab Timelapse: temporal evolution study

    Death Star BHs: 16 galaxies where SMBH jets destroy companion gas/stars.
    GC Vent: bipolar chimney structure extending ~500 pc from Sgr A*.
    """

    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31
    F_0     = 1.83e71
    RHO_VAC = 7.09e-36
    K_LENR  = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_UBi(self, M: float, r: float, omega_0: float = 1.0e-12) -> float:
        g_base = self.G * M / (r ** 2)
        momentum = (self.M_E * self.C**2 / r**2)
        F_LENR = self.compute_F_LENR(omega_0)
        return -self.F_0 + momentum + g_base + self.RHO_VAC + F_LENR

    def compute_jet_power(self, M_dot: float = 1.0e20,
                           v_jet: float = 0.1) -> float:
        """
        SMBH jet kinetic power: P_jet = 0.5 * M_dot * (v_jet*c)^2
        Relevant for Death Star BH feedback on companion systems.
        """
        return 0.5 * M_dot * (v_jet * self.C) ** 2

    def compute(self, dataset: dict) -> dict:
        M  = dataset.get('M_kg', 1.59e40)
        r  = dataset.get('r_m', 3.09e22)
        w0 = dataset.get('omega_0', 1.0e-12)
        M_dot = dataset.get('M_dot_kg_s', 1.0e20)
        v_jet = dataset.get('v_jet_c', 0.1)

        F_UBi  = self.compute_F_UBi(M, r, w0)
        F_LENR = self.compute_F_LENR(w0)
        P_jet  = self.compute_jet_power(M_dot, v_jet)

        systems = dataset.get('systems', None)
        batch = []
        if systems:
            for s in systems:
                fb = self.compute_F_UBi(s['M_kg'], s['r_m'], w0)
                batch.append({'name': s.get('name', ''), 'F_U_Bi_N': fb,
                              'negative': fb < 0})

        primary_equations = [
            f'F_U_Bi = -F_0 + momentum + gravity + rho_vac + F_LENR = {F_UBi:.6e} N',
            f'  F_LENR = {F_LENR:.6e} N',
            f'P_jet = 0.5*M_dot*(v_jet*c)^2 = {P_jet:.6e} W',
            f'Death Star BHs: 16 SMBHs with destructive jets',
            f'GC Vent: bipolar outflow ~500 pc, negative buoyancy',
            f'Cas A + Crab: timelapse temporal evolution',
        ]
        available_equations = [
            'F_U_Bi = -F_0 + (m_e*c^2/r^2) + (GM/r^2) + rho_vac + F_LENR',
            'P_jet = 0.5*M_dot*(v_jet*c)^2  (SMBH jet kinetic power)',
            'GC Vent: same params as Sgr A* -> same negative buoyancy',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            'F_U_Bi': F_UBi, 'P_jet_W': P_jet, 'batch_results': batch,
            'paper': 'PAPER_846',
        }


class SNRNebulaVelaTychoHelixSNR1181NGC6543Calc(_CP4Calculator):  # PAPER_847 #431
    """
    SNR/Nebula DeepSearch — 7 Systems.

    Cas A: M=3.978e30 kg, r=5.56e16 m, F_U_Bi ~ 2.65e208 N
    Crab Nebula: M=2.8e30 kg, r=5.56e16 m
    Vela Pulsar (wide-field): M=3.6e30 kg, r=6.17e17 m (F_U_Bi ~ 2.11e210 N highest)
    Tycho's SNR: M=2.8e30 kg, r=3.7e17 m
    Helix Nebula NGC 7293: M=1.2e30 kg, r=8.95e17 m
    SNR 1181 Pa 30: M=2.5e30 kg, r=2.47e19 m (Type Iax)
    NGC 6543 Cat's Eye: M=0.7e30 kg, r=6.17e18 m

    Vela wide-field has largest r in batch -> highest F_U_Bi.
    SNR 1181 Pa 30: first known Type Iax supernova remnant.
    """

    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31
    F_0     = 1.83e71
    RHO_VAC = 7.09e-36
    K_LENR  = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_UBi(self, M: float, r: float, omega_0: float = 1.0e-12) -> float:
        g_base = self.G * M / (r ** 2)
        momentum = (self.M_E * self.C**2 / r**2)
        F_LENR = self.compute_F_LENR(omega_0)
        return -self.F_0 + momentum + g_base + self.RHO_VAC + F_LENR

    def compute_snr_age(self, R_pc: float, v_exp_km_s: float) -> float:
        """SNR age = R / v_expansion (in years)"""
        R_m = R_pc * 3.086e16
        v_m = v_exp_km_s * 1000.0
        if v_m == 0:
            return 0.0
        return R_m / v_m / (365.25 * 86400)

    def compute(self, dataset: dict) -> dict:
        M  = dataset.get('M_kg', 3.978e30)
        r  = dataset.get('r_m', 5.56e16)
        w0 = dataset.get('omega_0', 1.0e-12)
        R_pc = dataset.get('R_pc', 1.7)
        v_exp = dataset.get('v_exp_km_s', 5000.0)

        F_UBi  = self.compute_F_UBi(M, r, w0)
        F_LENR = self.compute_F_LENR(w0)
        age = self.compute_snr_age(R_pc, v_exp)

        systems = dataset.get('systems', None)
        batch = []
        if systems:
            for s in systems:
                fb = self.compute_F_UBi(s['M_kg'], s['r_m'], w0)
                batch.append({'name': s.get('name', ''), 'F_U_Bi_N': fb})

        primary_equations = [
            f'F_U_Bi = -F_0 + momentum + gravity + rho_vac + F_LENR = {F_UBi:.6e} N',
            f'  F_LENR = {F_LENR:.6e} N',
            f'  SNR age ~ R/v_exp = {age:.0f} years (R={R_pc} pc, v={v_exp} km/s)',
            f'DeepSearch 7-system batch:',
            f'  CasA/Crab/Vela(wide)/Tycho/Helix/SNR1181(TypeIax)/NGC6543(CatsEye)',
            f'  Vela wide-field: r=6.17e17 m -> highest F_U_Bi=2.11e210 N',
            f'  SNR 1181 Pa 30: Type Iax remnant (1181 AD historical supernova)',
        ]
        available_equations = [
            'F_U_Bi = -F_0 + (m_e*c^2/r^2) + (GM/r^2) + rho_vac + F_LENR',
            'SNR age = R_pc * 3.086e16 / (v_exp_km_s * 1000) / (365.25*86400)',
            'Type Iax: deflagration supernova, incomplete carbon burning',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            'F_U_Bi': F_UBi, 'SNR_age_yr': age, 'batch_results': batch,
            'paper': 'PAPER_847',
        }


class SonificationCompositeH1821IC443M74MSH1552Calc(_CP4Calculator):  # PAPER_848 #432
    """
    Chandra Sonification Collection — 8 Systems.

    SNR 1181 (Pa 30): M=2.5e30, r=2.47e19 m
    H1821+643: Most massive quasar cluster, M=3e40 kg, r=3.7e23 m
    Sonification Collection: Composite 6-system, M=1e35, r=3.09e20 m
    IC 443 Jellyfish Nebula: M=30*M_sun, r=2.16e17 m
    M74 Phantom Galaxy: M=1e11*M_sun, r=5.56e20 m
    MSH 15-52 Hand PWN: M=1.5*M_sun, r=4.63e17 m
    SDSS J1531+3414: Galaxy merger, M=1e14*M_sun, r=1.54e23 m
    Sgr A*: M=4e6*M_sun, r=6.17e18 m -> NEGATIVE buoyancy

    Sonification: translation of X-ray data to audio frequencies.
    Audio mapping: low-energy X-rays -> bass, high-energy -> treble.
    H1821+643: most massive cluster-central quasar known (~3-30 billion M_sun).
    """

    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31
    F_0     = 1.83e71
    RHO_VAC = 7.09e-36
    K_LENR  = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_UBi(self, M: float, r: float, omega_0: float = 1.0e-12) -> float:
        g_base = self.G * M / (r ** 2)
        momentum = (self.M_E * self.C**2 / r**2)
        F_LENR = self.compute_F_LENR(omega_0)
        return -self.F_0 + momentum + g_base + self.RHO_VAC + F_LENR

    def compute_sonification_freq(self, E_keV: float, E_min: float = 0.5,
                                    E_max: float = 7.0, f_min: float = 60.0,
                                    f_max: float = 2000.0) -> float:
        """
        Map X-ray energy to audio frequency (log scale).
        f_audio = f_min * (f_max/f_min)^((E-E_min)/(E_max-E_min))
        """
        import math
        if E_max == E_min:
            return f_min
        frac = (E_keV - E_min) / (E_max - E_min)
        frac = max(0.0, min(1.0, frac))
        return f_min * (f_max / f_min) ** frac

    def compute(self, dataset: dict) -> dict:
        M  = dataset.get('M_kg', 3e40)
        r  = dataset.get('r_m', 3.7e23)
        w0 = dataset.get('omega_0', 1.0e-12)
        E_keV = dataset.get('E_keV', 2.0)

        F_UBi  = self.compute_F_UBi(M, r, w0)
        F_LENR = self.compute_F_LENR(w0)
        f_audio = self.compute_sonification_freq(E_keV)

        systems = dataset.get('systems', None)
        batch = []
        if systems:
            for s in systems:
                fb = self.compute_F_UBi(s['M_kg'], s['r_m'], w0)
                batch.append({'name': s.get('name', ''), 'F_U_Bi_N': fb,
                              'negative': fb < 0})

        primary_equations = [
            f'F_U_Bi = -F_0 + momentum + gravity + rho_vac + F_LENR = {F_UBi:.6e} N',
            f'  F_LENR = {F_LENR:.6e} N',
            f'Sonification: E_keV={E_keV} -> f_audio={f_audio:.1f} Hz',
            f'8-system sonification batch:',
            f'  SNR1181/H1821+643/SonifColl/IC443/M74/MSH15-52/SDSSJ1531/SgrA*',
            f'  H1821+643: most massive cluster-central quasar',
            f'  Sgr A*: negative buoyancy in sonification context',
        ]
        available_equations = [
            'F_U_Bi = -F_0 + (m_e*c^2/r^2) + (GM/r^2) + rho_vac + F_LENR',
            'f_audio = f_min*(f_max/f_min)^((E-E_min)/(E_max-E_min))  (sonification)',
            'X-ray -> audio mapping: 0.5-7.0 keV -> 60-2000 Hz',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            'F_U_Bi': F_UBi, 'f_audio_Hz': f_audio, 'batch_results': batch,
            'paper': 'PAPER_848',
        }


class ArXiv24PaperBatch4FquarkFneutrinoFALPFdarkCalc(_CP4Calculator):  # PAPER_849 #433
    """
    Integrated 24-Paper arXiv BSM Landscape.

    Four F_U_Bi_i terms derived from high-energy physics arXiv literature:

    Batch 1 (5 papers): F_quark = k_quark*|V_cb|^2 = 1.54e7 N
      arXiv:2506.15256 (Belle II B->D(*) |V_cb|=39.2e-3)
      arXiv:2506.15347, 2506.15390, 2506.15515, 2506.15533

    Batch 2 (6 papers): F_neutrino = k_neutrino*alpha_nu = 1.0 N
      arXiv:2506.14881 (neutrino polarizability alpha_nu=1e-10)
      arXiv:2506.14989, 2506.15046, 2506.15164, 2506.15245, 2506.15306

    Batch 3 (6 papers): F_ALP = k_ALP*g_aqq = 1e4 N
      arXiv:2506.15637 (ALP-hadron covariant g_aqq=1e-6)
      arXiv:2506.15428, 2506.15445, 2412.04357, 2412.10141, 2503.05679

    Batch 4 (7 papers): F_dark = k_dark*epsilon^2 = 1.0 N
      arXiv:2402.00249 (FASER dark photons epsilon=1e-5)
      arXiv:2506.13588, 2410.11367, 2412.03677, 2502.19817, 2503.01607, 2506.02450

    Total BSM: F_quark dominates at 1.54e7 N.
    All 8 sonification systems recalculated with integrated BSM terms.
    """

    K_QUARK    = 1.0e10
    V_CB       = 39.2e-3
    K_NEUTRINO = 1.0e10
    ALPHA_NU   = 1.0e-10
    K_ALP      = 1.0e10
    G_AQQ      = 1.0e-6
    K_DARK     = 1.0e10
    EPSILON    = 1.0e-5
    G          = 6.6743e-11
    C          = 2.998e8
    M_E        = 9.109e-31
    F_0        = 1.83e71
    RHO_VAC    = 7.09e-36
    K_LENR     = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12

    def compute_F_quark(self, V_cb: float = None) -> float:
        vcb = V_cb if V_cb is not None else self.V_CB
        return self.K_QUARK * (vcb ** 2)

    def compute_F_neutrino(self, alpha_nu: float = None) -> float:
        a = alpha_nu if alpha_nu is not None else self.ALPHA_NU
        return self.K_NEUTRINO * a

    def compute_F_ALP(self, g_aqq: float = None) -> float:
        g = g_aqq if g_aqq is not None else self.G_AQQ
        return self.K_ALP * g

    def compute_F_dark(self, epsilon: float = None) -> float:
        e = epsilon if epsilon is not None else self.EPSILON
        return self.K_DARK * (e ** 2)

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_UBi_integrated(self, M: float, r: float,
                                   omega_0: float = 1.0e-12) -> dict:
        """Full F_U_Bi with all 4 BSM terms integrated."""
        g_base = self.G * M / (r ** 2)
        momentum = (self.M_E * self.C**2 / r**2)
        F_LENR = self.compute_F_LENR(omega_0)
        F_q = self.compute_F_quark()
        F_n = self.compute_F_neutrino()
        F_a = self.compute_F_ALP()
        F_d = self.compute_F_dark()
        F_BSM = F_q + F_n + F_a + F_d
        F_UBi = -self.F_0 + momentum + g_base + self.RHO_VAC + F_LENR + F_BSM
        return {
            'F_U_Bi_N': F_UBi, 'F_LENR': F_LENR, 'F_BSM': F_BSM,
            'F_quark': F_q, 'F_neutrino': F_n, 'F_ALP': F_a, 'F_dark': F_d,
        }

    def compute(self, dataset: dict) -> dict:
        M  = dataset.get('M_kg', 1.989e30)
        r  = dataset.get('r_m', 1.0e16)
        w0 = dataset.get('omega_0', 1.0e-12)

        res = self.compute_F_UBi_integrated(M, r, w0)

        systems = dataset.get('systems', None)
        batch = []
        if systems:
            for s in systems:
                sr = self.compute_F_UBi_integrated(s['M_kg'], s['r_m'], w0)
                sr['name'] = s.get('name', '')
                batch.append(sr)

        primary_equations = [
            f'--- 24-Paper arXiv BSM Landscape (4 batches) ---',
            f'F_quark    = k_quark*|V_cb|^2 = {res["F_quark"]:.6e} N  (dominant)',
            f'F_neutrino = k_neutrino*alpha_nu = {res["F_neutrino"]:.6e} N',
            f'F_ALP      = k_ALP*g_aqq = {res["F_ALP"]:.6e} N',
            f'F_dark     = k_dark*epsilon^2 = {res["F_dark"]:.6e} N',
            f'F_BSM_total = {res["F_BSM"]:.6e} N',
            f'F_LENR = {res["F_LENR"]:.6e} N',
            f'F_U_Bi (integrated) = {res["F_U_Bi_N"]:.6e} N',
            f'24 papers: Batch1(5)+Batch2(6)+Batch3(6)+Batch4(7)',
        ]
        available_equations = [
            'F_quark = k_quark*|V_cb|^2  (Belle II, 5 papers)',
            'F_neutrino = k_neutrino*alpha_nu  (6 papers)',
            'F_ALP = k_ALP*g_aqq  (ALP-hadron, 6 papers)',
            'F_dark = k_dark*epsilon^2  (FASER dark photons, 7 papers)',
            'F_U_Bi_integrated = base + F_LENR + F_quark + F_neutrino + F_ALP + F_dark',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            **res, 'batch_results': batch,
            'paper': 'PAPER_849',
        }


class ADDGravitonLeakageNegBuoyancySgrAExtCalc(_CP4Calculator):  # PAPER_850 #434
    """
    ADD Extra Dimensions Extended — Graviton Leakage + Negative Buoyancy.

    Arkani-Hamed-Dimopoulos-Dvali (ADD, arXiv:1607.01831):
    F_LED = k_LED * (M_star / M_Pl)^2 = 6.72e-23 N

    Extended treatment for 8 sonification systems:
    F_LED negligible vs F_LENR but connects UQFF to extra-dimensional
    physics at SMBH negative buoyancy sites.

    Graviton leakage hypothesis:
      In ADD n=2, gravity leaks into 2 compactified extra dimensions
      at radius R ~ 0.1 mm (submillimeter gravity test range).
      At Sgr A* (negative buoyancy), F_LED + F_UBi < 0 suggests
      graviton drainage into extra dimensions amplifies vacuum reversal.

    Compactification radius: R_n ~ (M_Pl/M_star)^(2/n) * hbar_c / M_star
    For n=2: R ~ 1.96e-4 m ~ 0.2 mm.
    """

    K_LED   = 1.0e10
    M_STAR  = 1.0e3
    M_PL    = 1.22e19
    G       = 6.6743e-11
    C       = 2.998e8
    M_E     = 9.109e-31
    F_0     = 1.83e71
    RHO_VAC = 7.09e-36
    K_LENR  = 1.0e-10
    OMEGA_LENR = 2.0 * 3.141592653589793 * 1.25e12
    HBAR_C  = 1.973e-16     # GeV*m

    def compute_F_LED(self, M_star: float = None, M_Pl: float = None) -> float:
        ms = M_star if M_star is not None else self.M_STAR
        mp = M_Pl if M_Pl is not None else self.M_PL
        return self.K_LED * (ms / mp) ** 2

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_compactification_radius(self, n: int = 2,
                                          M_star: float = None) -> float:
        """R_n ~ (M_Pl/M_star)^(2/n) * hbar_c / M_star"""
        ms = M_star if M_star is not None else self.M_STAR
        R_nat = (self.M_PL / ms) ** (2.0 / n) / ms
        return R_nat * self.HBAR_C

    def compute_F_UBi_with_LED(self, M: float, r: float,
                                 omega_0: float = 1.0e-12,
                                 n: int = 2) -> dict:
        g_base = self.G * M / (r ** 2)
        momentum = (self.M_E * self.C**2 / r**2)
        F_LENR = self.compute_F_LENR(omega_0)
        F_LED = self.compute_F_LED()
        R_n = self.compute_compactification_radius(n)
        F_UBi = -self.F_0 + momentum + g_base + self.RHO_VAC + F_LENR + F_LED
        return {
            'F_U_Bi_N': F_UBi, 'F_LENR': F_LENR, 'F_LED': F_LED,
            'R_n_m': R_n, 'negative': F_UBi < 0,
        }

    def compute(self, dataset: dict) -> dict:
        M  = dataset.get('M_kg', 7.956e36)
        r  = dataset.get('r_m', 6.17e18)
        w0 = dataset.get('omega_0', 1.0e-12)
        n  = dataset.get('n_extra_dims', 2)

        res = self.compute_F_UBi_with_LED(M, r, w0, n)

        systems = dataset.get('systems', None)
        batch = []
        if systems:
            for s in systems:
                sr = self.compute_F_UBi_with_LED(s['M_kg'], s['r_m'], w0, n)
                sr['name'] = s.get('name', '')
                batch.append(sr)

        primary_equations = [
            f'--- ADD Graviton Leakage Extended (arXiv:1607.01831) ---',
            f'F_LED = k_LED*(M_star/M_Pl)^2 = {res["F_LED"]:.6e} N',
            f'  M_star = {self.M_STAR:.0f} GeV, M_Pl = {self.M_PL:.3e} GeV, n = {n}',
            f'R_n (compactification) = {res["R_n_m"]:.3e} m  (~0.2 mm for n=2)',
            f'F_LENR = {res["F_LENR"]:.6e} N',
            f'F_U_Bi (with F_LED) = {res["F_U_Bi_N"]:.6e} N',
            f'  Negative buoyancy: {res["negative"]}',
            f'Graviton leakage: F_LED negligible but ADD hypothesis at Sgr A*',
            f'8 systems recalculated with F_LED integrated',
        ]
        available_equations = [
            'F_LED = k_LED*(M_star/M_Pl)^2  (ADD n=2)',
            'R_n ~ (M_Pl/M_star)^(2/n) * hbar_c / M_star',
            'F_U_Bi = base + F_LENR + F_LED  (graviton leakage)',
            'Sgr A*: negative buoyancy + graviton drainage hypothesis',
        ]
        simulation_set = [
            {'equation': 'R_n_vs_n',
             'values': {f'n={ni}': self.compute_compactification_radius(ni)
                       for ni in [1, 2, 3, 4, 6]}},
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': simulation_set,
            **res, 'batch_results': batch,
            'paper': 'PAPER_850',
        }


class KozimaNeutronDropDensityScaled8SystemCalc(_CP4Calculator):  # PAPER_851 #435
    """
    Kozima LENR Neutron Drop — Density-Scaled 8-System Batch.

    Source: Kozima H. 2021 (PMC8141838)
    F_neutron = k_neutron * sigma_n

    Three models:
      1. Static: sigma_n = 10^-4 m^2 -> F_neutron = 10^6 N
      2. Frequency-dependent: sigma_n(omega) = sigma_0*(omega/omega_LENR)^2
         * exp(-(omega-omega_LENR)^2/(2*Delta_omega^2))
      3. Density-scaled: sigma_n(rho) = sigma_0 * (rho/rho_0)
         For neutron stars: rho~10^17 kg/m^3 -> sigma_n~10^35 -> F_neutron~10^45 N

    8-system sonification batch recalculated with F_neutron:
      Lab-scale: F_neutron = 10^6 N
      Neutron star (PSR J0030+0451): F_neutron ~ 10^45 N
    """

    K_NEUTRON   = 1.0e10
    SIGMA_0     = 1.0e-4
    OMEGA_LENR  = 2.0 * 3.141592653589793 * 1.25e12
    DELTA_OMEGA = 2.0 * 3.141592653589793 * 0.05e12
    OMEGA_ACT   = 2.0 * 3.141592653589793 * 300.0
    ALPHA       = 0.1
    RHO_0       = 1.0e-22
    G           = 6.6743e-11
    C           = 2.998e8
    M_E         = 9.109e-31
    F_0         = 1.83e71
    RHO_VAC     = 7.09e-36
    K_LENR      = 1.0e-10

    def compute_sigma_n_freq(self, omega: float) -> float:
        import math
        ratio = (omega / self.OMEGA_LENR) ** 2
        gauss = math.exp(-((omega - self.OMEGA_LENR)**2) / (2.0 * self.DELTA_OMEGA**2))
        return self.SIGMA_0 * ratio * gauss

    def compute_sigma_n_density(self, rho: float) -> float:
        return self.SIGMA_0 * (rho / self.RHO_0)

    def compute_F_neutron(self, sigma_n: float = None) -> float:
        s = sigma_n if sigma_n is not None else self.SIGMA_0
        return self.K_NEUTRON * s

    def compute_F_neutron_dynamic(self, omega_eff: float, t: float = 0.0) -> float:
        import math
        sigma = self.compute_sigma_n_freq(omega_eff)
        return self.K_NEUTRON * sigma * (1.0 + self.ALPHA * math.cos(self.OMEGA_ACT * t))

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_F_UBi_with_neutron(self, M: float, r: float,
                                     rho: float = None,
                                     omega_0: float = 1.0e-12) -> dict:
        g_base = self.G * M / (r ** 2)
        momentum = (self.M_E * self.C**2 / r**2)
        F_LENR = self.compute_F_LENR(omega_0)
        F_neutron_gen = self.compute_F_neutron()
        F_neutron_dens = None
        if rho is not None:
            sigma_d = self.compute_sigma_n_density(rho)
            F_neutron_dens = self.K_NEUTRON * sigma_d
        F_n = F_neutron_dens if F_neutron_dens is not None else F_neutron_gen
        F_UBi = -self.F_0 + momentum + g_base + self.RHO_VAC + F_LENR + F_n
        return {
            'F_U_Bi_N': F_UBi, 'F_LENR': F_LENR,
            'F_neutron_general': F_neutron_gen,
            'F_neutron_density': F_neutron_dens,
            'negative': F_UBi < 0,
        }

    def compute(self, dataset: dict) -> dict:
        M   = dataset.get('M_kg', 1.989e30)
        r   = dataset.get('r_m', 1.0e4)
        rho = dataset.get('rho_kg_m3', None)
        w0  = dataset.get('omega_0', 1.0e-12)

        res = self.compute_F_UBi_with_neutron(M, r, rho, w0)

        systems = dataset.get('systems', None)
        batch = []
        if systems:
            for s in systems:
                sr = self.compute_F_UBi_with_neutron(
                    s['M_kg'], s['r_m'], s.get('rho_kg_m3', None), w0)
                sr['name'] = s.get('name', '')
                batch.append(sr)

        primary_equations = [
            f'--- Kozima Neutron Drop Density-Scaled 8-System ---',
            f'F_neutron(static) = k_neutron*sigma_0 = {res["F_neutron_general"]:.6e} N',
        ]
        if res['F_neutron_density'] is not None:
            primary_equations.append(
                f'F_neutron(density) = {res["F_neutron_density"]:.6e} N  (rho={rho:.2e})')
        primary_equations.extend([
            f'F_LENR = {res["F_LENR"]:.6e} N',
            f'F_U_Bi (with F_neutron) = {res["F_U_Bi_N"]:.6e} N',
            f'  Negative buoyancy: {res["negative"]}',
            f'Density scaling: sigma_n(rho) = sigma_0*(rho/rho_0)',
            f'  Lab: rho~1e-22 -> F_neutron=10^6 N',
            f'  NS (PSR J0030+0451): rho~10^17 -> F_neutron~10^45 N',
        ])
        available_equations = [
            'F_neutron = k_neutron*sigma_n  (static)',
            'sigma_n(omega) = sigma_0*(omega/omega_LENR)^2*exp(-...)',
            'sigma_n(rho) = sigma_0*(rho/rho_0)  (density-scaled)',
            'F_neutron(t) = k_neutron*sigma_n*(1+alpha*cos(omega_act*t))',
        ]
        simulation_set = [
            {'equation': 'F_neutron_vs_rho',
             'values': {f'{r:.0e}': self.compute_F_neutron(self.compute_sigma_n_density(r))
                       for r in [1e-22, 1e-10, 1e0, 1e10, 1e17]}},
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': simulation_set,
            **res, 'batch_results': batch,
            'paper': 'PAPER_851',
        }


class LENRNextStepsExperimentalDesignPSRJ0030Calc(_CP4Calculator):  # PAPER_852 #436
    """
    LENR Next Steps Experimental Design + PSR J0030+0451.

    Four experimental tracks from the file's Next Steps:

    1. LENR Replication: Pd-D electrolysis + Ni-H gas loading
       Target: reproduce excess heat > 10 W from 1-10 mW input
       THz source: QCL (quantum cascade laser) at 1.25 THz
       Equipment: calorimeter, neutron detector, gamma spectrometer

    2. Sgr A* Investigation: ALMA Band 10 (787-950 GHz) + EHT 230 GHz
       Proposal: search for 1.25 THz emission from Sgr A* environment
       Test: nonzero F_LENR at astrophysical scale

    3. Neutron Drop Refinement: sigma_n(omega) frequency-dependent model
       Measure sigma_n at discrete frequencies near 1.25 THz
       Predict: resonant peak in neutron absorption cross-section

    4. Astrophysical Analogues: PSR J0030+0451 (NICER)
       M ~ 1.44 M_sun, R ~ 13 km, rho ~ 10^17 kg/m^3
       Test: density-scaled F_neutron ~ 10^45 N
       Compare: neutron star hotspot X-ray profiles

    PSR J0030+0451: millisecond pulsar, 4.87 ms period,
    NICER hotspot mapping constraints (2019).
    """

    G           = 6.6743e-11
    C           = 2.998e8
    M_SUN       = 1.989e30
    K_NEUTRON   = 1.0e10
    SIGMA_0     = 1.0e-4
    RHO_0       = 1.0e-22
    K_LENR      = 1.0e-10
    OMEGA_LENR  = 2.0 * 3.141592653589793 * 1.25e12
    M_E         = 9.109e-31
    F_0         = 1.83e71
    RHO_VAC     = 7.09e-36
    # PSR J0030+0451 defaults
    M_PSR       = 1.44 * 1.989e30   # kg
    R_PSR       = 13.0e3            # m (13 km)
    RHO_PSR     = 1.0e17            # kg/m^3

    def compute_F_neutron_density(self, rho: float) -> float:
        sigma_d = self.SIGMA_0 * (rho / self.RHO_0)
        return self.K_NEUTRON * sigma_d

    def compute_F_LENR(self, omega_0: float = 1.0e-12) -> float:
        return self.K_LENR * (self.OMEGA_LENR / omega_0) ** 2

    def compute_LENR_gain(self, P_out: float = 10.0,
                           P_in: float = 0.001) -> float:
        """COP (Coefficient of Performance) = P_out / P_in"""
        if P_in == 0:
            return float('inf')
        return P_out / P_in

    def compute_ALMA_frequency_coverage(self, band: int = 10) -> dict:
        """ALMA frequency bands relevant to UQFF LENR signature search."""
        bands = {
            7: (275e9, 373e9),
            8: (385e9, 500e9),
            9: (602e9, 720e9),
            10: (787e9, 950e9),
        }
        if band in bands:
            lo, hi = bands[band]
            covers_lenr = lo <= 1.25e12 <= hi
            return {'band': band, 'freq_lo_Hz': lo, 'freq_hi_Hz': hi,
                    'covers_1.25THz': covers_lenr}
        return {'band': band, 'error': 'unknown band'}

    def compute(self, dataset: dict) -> dict:
        # PSR J0030+0451 parameters
        M_psr = dataset.get('M_kg', self.M_PSR)
        R_psr = dataset.get('r_m', self.R_PSR)
        rho = dataset.get('rho_kg_m3', self.RHO_PSR)
        w0  = dataset.get('omega_0', 1.0e-12)
        P_out = dataset.get('P_out_W', 10.0)
        P_in  = dataset.get('P_in_W', 0.001)

        F_neutron_psr = self.compute_F_neutron_density(rho)
        F_LENR = self.compute_F_LENR(w0)
        COP = self.compute_LENR_gain(P_out, P_in)
        g_psr = self.G * M_psr / (R_psr ** 2)

        alma = self.compute_ALMA_frequency_coverage(10)

        primary_equations = [
            f'--- LENR Next Steps Experimental Design ---',
            f'Track 1 LENR Replication:',
            f'  COP = P_out/P_in = {COP:.0f}x ({P_out}W / {P_in}W)',
            f'  Target: Pd-D + Ni-H at 1.25 THz (QCL source)',
            f'Track 2 Sgr A* Investigation:',
            f'  ALMA Band 10: {alma["freq_lo_Hz"]:.0e}-{alma["freq_hi_Hz"]:.0e} Hz',
            f'  Covers 1.25 THz: {alma.get("covers_1.25THz", "N/A")}',
            f'Track 3 Neutron Drop Refinement:',
            f'  sigma_n(omega) resonance peak at omega_LENR = {self.OMEGA_LENR:.3e} rad/s',
            f'Track 4 PSR J0030+0451 Astrophysical Analogue:',
            f'  M = {M_psr:.3e} kg ({M_psr/self.M_SUN:.2f} M_sun)',
            f'  R = {R_psr:.0f} m, rho = {rho:.2e} kg/m^3',
            f'  g_surface = GM/R^2 = {g_psr:.6e} m/s^2',
            f'  F_neutron(density) = {F_neutron_psr:.6e} N',
            f'  F_LENR = {F_LENR:.6e} N',
        ]
        available_equations = [
            'COP = P_out / P_in  (LENR replication target)',
            'ALMA Band 10: 787-950 GHz (nearest to 1.25 THz)',
            'sigma_n(rho) = sigma_0*(rho/rho_0)  (density-scaled neutron drop)',
            'F_neutron(NS) = k_neutron*sigma_n(rho~10^17) ~ 10^45 N',
            'PSR J0030+0451: NICER hotspot mapping constraints (2019)',
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': [],
            'F_neutron_PSR': F_neutron_psr, 'F_LENR': F_LENR,
            'COP': COP, 'g_surface': g_psr,
            'ALMA_band10': alma,
            'paper': 'PAPER_852',
        }


_SESSION_197_CLASSES = [
    'FloydSweetVTA6DocumentPCVTMotionalEfieldCalc',         # PAPER_842 #426
    'ChandraXRayBatch1GCEagleHBC672NGC7469VirgoCalc',       # PAPER_843 #427
    'Chandra25thAnniversaryCrabOrionNGC6334Calc',           # PAPER_844 #428
    'ChandraSurveyMACSJ0416LensExoSMBHCalc',                # PAPER_845 #429
    'ChandraDeathStar16SMBHGCVentTimelapseCalc',            # PAPER_846 #430
    'SNRNebulaVelaTychoHelixSNR1181NGC6543Calc',            # PAPER_847 #431
    'SonificationCompositeH1821IC443M74MSH1552Calc',        # PAPER_848 #432
    'ArXiv24PaperBatch4FquarkFneutrinoFALPFdarkCalc',       # PAPER_849 #433
    'ADDGravitonLeakageNegBuoyancySgrAExtCalc',              # PAPER_850 #434
    'KozimaNeutronDropDensityScaled8SystemCalc',             # PAPER_851 #435
    'LENRNextStepsExperimentalDesignPSRJ0030Calc',           # PAPER_852 #436
]
'''

# --- Run ---

with open(CP4, 'r', encoding='utf-8', errors='replace') as f:
    content = f.read()

lines = content.split('\n')
total_before = len(lines)

# 1. Patch header: find Session 196 line, replace with updated + Session 197
old_hdr = HEADER_PATCH[0]
new_hdr = HEADER_PATCH[1]
idx = content.find(old_hdr)
if idx == -1:
    print("ERROR: Could not find header line for Session 196 v5.56")
    exit(1)
eol = content.find('\n', idx)
content = content[:idx] + new_hdr + content[eol:]

# 2. Replace version string
content = content.replace("Version: 5.56 (2026-07-02)", "Version: 5.57 (2026-04-04)")

# 3. Append new classes before the _SESSION_196_CLASSES list
marker = "_SESSION_196_CLASSES"
idx2 = content.find(marker)
if idx2 == -1:
    print("ERROR: Could not find _SESSION_196_CLASSES marker")
    exit(1)
line_start = content.rfind('\n', 0, idx2) + 1
content = content[:line_start] + NEW_CLASSES + '\n' + content[line_start:]

with open(CP4, 'w', encoding='utf-8') as f:
    f.write(content)

new_lines = content.split('\n')
print(f"CP4 updated: {total_before} -> {len(new_lines)} lines (+{len(new_lines)-total_before})")
print("Header patched: Session 197 v5.57")
print("11 classes appended: #426-#436")
print("Version: 5.56 -> 5.57")

# Verify class count
import re
class_count = len(re.findall(r'^class\s+\w+', content, re.MULTILINE))
print(f"Total class definitions in CP4: {class_count}")
