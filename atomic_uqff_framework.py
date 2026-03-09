"""
Atomic UQFF Framework — Z-Dependent Scaling
============================================

Source: Grok thread (grok_share_7b0e961f_conversation.txt), Sept 21, 2025.
Atomic phase: Z=1–118 full Periodic Table DPM-derived framework.

Physics:
    p_Z = Z / (Z+1)                          — DPM proportion (H: 0.5, Og: 0.992)
    SSq_Z = 0.507 + (Z/118) × 0.1           — Z-dependent shell quotient
    fulcrum = p_Z × (1 - SSq_Z)             — buoyancy-gravity balance ratio
    shells_Z = ceil(Z/8)                     — atomic shell count
    δρ/ρ_Z → mean 5.04×10⁻⁵, std 2.89×10⁻⁵  — density perturbation (KS p=0.741 ~ normal)
    ρ_vac,UA_Z → mean 2.69×10⁻¹⁰ J/m³ (CMB-matched)

Statistical Validation (Grok thread):
    Jarque-Bera:       stat=8.78,  p=0.012   (Q_wave — rejects normality)
    Shapiro-Wilk:      stat=0.955, p=0.00055 (rejects normality)
    KS test (δρ/ρ):    stat=0.061, p=0.741   (fails to reject — near-normal)
    Anderson-Darling:  stat=1.35 > 1.092 critical

DPM Mayan Periodic Table (5 expansion cycles):
    3 stable states (solid, liquid, gas) + 1 unstable (plasma)
    4 series: Alkaline, Halogen, Lanthanide, Actinide
    Element series defined by p_Z thresholds

Author: Daniel T. Murphy
Framework: Universal Quantum Field Superconductive Framework (UQFF)
"""

import math
import numpy as np
from typing import Dict, Optional


class AtomicUQFFFramework:
    """
    Z-dependent UQFF scaling for full Periodic Table (Z=1–118).

    Calibrated constants:
        SSq_base = 0.507 (range 0.507–0.607 over Z=1–118)
        δρ/ρ mean = 5.04×10⁻⁵  (KS p=0.741, near-normal)
        Q_wave mean = 3.97e4 J/m³  (std 6.33e4)
        F_U_Bi_i mean = -6.06×10²¹⁷ N (47-system bootstrap, 3% std)
        ρ_vac,UA_Z (CMB) = 2.69×10⁻¹⁰ J/m³ (vs. baseline 7.09e-36)
    """

    # Calibrated statistical constants from 47-system thread bootstrap
    DELTA_RHO_MEAN = 5.04e-5       # Z-dependent density perturbation mean
    DELTA_RHO_STD  = 2.89e-5       # Z-dependent density perturbation std
    Q_WAVE_MEAN    = 3.97e4        # J/m³ — DPM resonance + phase mean
    Q_WAVE_STD     = 6.33e4        # J/m³ — Q_wave standard deviation
    F_UBI_MEAN     = -6.06e217     # N — 47-system log bootstrap mean
    F_UBI_LOG_STD  = 0.030         # Log bootstrap std (3%)
    RHO_VAC_UA_CMB = 2.69e-10      # J/m³ — CMB-matched rho_vac_UA_Z
    RHO_VAC_UA_BASE = 7.09e-36     # J/m³ — baseline rho_vac_UA (non-Z)

    # DPM p_Z thresholds for Mayan Periodic Table states
    P_Z_SOLID_MAX    = 0.60
    P_Z_LIQUID_MAX   = 0.75
    P_Z_GAS_MAX      = 0.90
    # above 0.90 → plasma (unstable)

    def p_Z(self, Z: int) -> float:
        """
        DPM proportion for element Z.

        p_Z = Z / (Z + 1)

        Range: 0.5 (H, Z=1) → 0.992 (Og, Z=118)
        """
        return Z / (Z + 1.0)

    def SSq_Z(self, Z: int) -> float:
        """
        Z-dependent UQFF shell quotient [SSq]_Z.

        SSq_Z = 0.507 + (Z/118) × 0.1

        Range: 0.508 (H) → 0.607 (Og)
        SSq_Z base (0.507) matches [SSq]=0.57 calibrated at NS merger Ye~0.1
        """
        return 0.507 + (Z / 118.0) * 0.1

    def fulcrum_ratio(self, Z: int) -> float:
        """
        Buoyancy-gravity fulcrum ratio.

        F_buoy / F_grav = p_Z × (1 - SSq_Z)

        Range: H: 0.246 → Og: 0.390
        Stable atoms have fulcrum < 0.5 (gravity-dominated; buoyancy ~25–39%)
        """
        return self.p_Z(Z) * (1.0 - self.SSq_Z(Z))

    def shells_Z(self, Z: int) -> int:
        """
        Atomic shell count for element Z.

        shells_Z = ceil(Z / 8)

        Examples: H(1)→1, Ne(10)→2, Ar(18)→3, Kr(36)→5, Rn(86)→11, Og(118)→15
        """
        return math.ceil(Z / 8.0)

    def delta_rho_over_rho_Z(self, Z: int) -> float:
        """
        Z-dependent vacuum density perturbation.

        δρ/ρ_Z ≈ DELTA_RHO_MEAN + (Z/118) × DELTA_RHO_STD

        Statistical basis: Z-independent mean 5.04×10⁻⁵ with Z-gradient.
        KS test p=0.741 (δρ/ρ sample ~ normal) validated in Grok thread.
        """
        return self.DELTA_RHO_MEAN + (Z / 118.0) * self.DELTA_RHO_STD

    def rho_vac_UA_Z(self, Z: int) -> float:
        """
        Z-dependent CMB-matched aether vacuum energy density.

        Scaled from CMB-matched value 2.69×10⁻¹⁰ J/m³ using p_Z:
        ρ_vac,UA_Z ≈ RHO_VAC_UA_CMB × p_Z

        CMB-matched vs baseline 7.09e-36 J/m³ divergence:
        rho_vac_UA_Z / rho_vac_UA_base ~ 10^{26} (large due to CMB rescaling)
        """
        return self.RHO_VAC_UA_CMB * self.p_Z(Z)

    def buoyancy_gravity_balance(self, Z: int, F_grav: float) -> dict:
        """
        Full buoyancy-gravity balance for element Z.

        F_buoy = F_grav × fulcrum_ratio(Z)

        Parameters:
            Z:      Atomic number (1–118)
            F_grav: Gravitational force on nucleus (N)

        Returns:
            dict with Z, p_Z, SSq_Z, fulcrum_ratio, F_buoy, F_grav
        """
        ratio = self.fulcrum_ratio(Z)
        return {
            'Z': Z,
            'element_symbol': self._element_symbol(Z),
            'p_Z': self.p_Z(Z),
            'SSq_Z': self.SSq_Z(Z),
            'shells_Z': self.shells_Z(Z),
            'fulcrum_ratio': ratio,
            'F_buoy': F_grav * ratio,
            'F_grav': F_grav,
            'stability': 'gravity-dominated' if ratio < 0.5 else 'buoyancy-dominated',
        }

    def r_process_threshold(self, n: int = 13) -> tuple:
        """
        Check if UQFF suppression allows r-process A>140 nucleosynthesis.

        exp(-[SSq]·n/26) < 0.9 → A>140 r-process nuclei accessible
        [SSq]=0.57 calibrated from Ye~0.1 NS merger disk outflows.

        Parameters:
            n: Shell index (default 13 = half of 26 UQFF levels)

        Returns:
            (threshold_met: bool, suppression_factor: float)
        """
        SSq = 0.57  # calibrated from NS merger Ye~0.1
        suppression = np.exp(-SSq * n / 26.0)
        return suppression < 0.9, suppression

    def mayan_state(self, Z: int) -> str:
        """
        DPM Mayan Periodic Table state classification.

        5-cycle DPM derivation: 3 stable (solid, liquid, gas) + 1 unstable (plasma)
        4 series: Alkaline, Halogen, Lanthanide, Actinide (by p_Z thresholds)

        Parameters:
            Z: Atomic number

        Returns:
            State string (solid/liquid/gas/plasma)
        """
        p = self.p_Z(Z)
        if p < self.P_Z_SOLID_MAX:
            return 'solid'
        elif p < self.P_Z_LIQUID_MAX:
            return 'liquid'
        elif p < self.P_Z_GAS_MAX:
            return 'gas'
        else:
            return 'plasma (unstable)'

    def mayan_series(self, Z: int) -> str:
        """
        DPM-derived Periodic Table series classification.

        Based on electron shell filling (approximation; full DPM uses p_Z cycles).
        """
        if Z in (1, 2):
            return 'Period 1 / Light'
        elif Z <= 10:
            s_count = Z - 2
            if Z in (3, 4):
                return 'Alkaline'
            elif Z in (9, 10):
                return 'Halogen/Noble'
        elif 57 <= Z <= 71:
            return 'Lanthanide'
        elif 89 <= Z <= 103:
            return 'Actinide'
        elif Z in range(1, 30, 2):
            return 'Alkaline'
        return 'Main Group'

    def scan_periodic_table(self) -> list:
        """
        Full Z=1–118 scan with all UQFF Z-dependent quantities.

        Returns:
            list of dicts for Z=1 to 118
        """
        results = []
        for Z in range(1, 119):
            entry = {
                'Z': Z,
                'symbol': self._element_symbol(Z),
                'p_Z': round(self.p_Z(Z), 6),
                'SSq_Z': round(self.SSq_Z(Z), 6),
                'fulcrum': round(self.fulcrum_ratio(Z), 6),
                'shells_Z': self.shells_Z(Z),
                'delta_rho': round(self.delta_rho_over_rho_Z(Z), 3e-7),
                'rho_vac_UA_Z': self.rho_vac_UA_Z(Z),
                'state': self.mayan_state(Z),
            }
            results.append(entry)
        return results

    @staticmethod
    def _element_symbol(Z: int) -> str:
        """Basic Z→symbol map for Z=1–118."""
        symbols = [
            '', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
            'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
            'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
            'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
            'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
            'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
            'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
        ]
        if 1 <= Z <= 118:
            return symbols[Z]
        return f'Z{Z}'


# =============================================================================
# STANDALONE USAGE
# =============================================================================

if __name__ == '__main__':
    fw = AtomicUQFFFramework()
    print("=== Atomic UQFF Framework — Z-Dependent Scaling ===\n")

    # Key elements
    for Z in [1, 6, 8, 26, 56, 92, 118]:
        symbol = fw._element_symbol(Z)
        p = fw.p_Z(Z)
        SSq = fw.SSq_Z(Z)
        fulcrum = fw.fulcrum_ratio(Z)
        state = fw.mayan_state(Z)
        print(f"  Z={Z:3d} ({symbol:2s}): p_Z={p:.4f}, SSq_Z={SSq:.4f}, "
              f"fulcrum={fulcrum:.4f}, state={state}")

    print("\nr-process threshold check (n=13, [SSq]=0.57):")
    met, suppression = fw.r_process_threshold()
    print(f"  suppression = {suppression:.4f}, A>140 accessible: {met}")

    print(f"\nStatistical constants from 47-system bootstrap:")
    print(f"  F_U_Bi_i mean = {fw.F_UBI_MEAN:.2e} N (log std {fw.F_UBI_LOG_STD})")
    print(f"  Q_wave mean   = {fw.Q_WAVE_MEAN:.2e} J/m³")
    print(f"  δρ/ρ mean     = {fw.DELTA_RHO_MEAN:.2e}")
    print(f"  ρ_vac_UA_CMB  = {fw.RHO_VAC_UA_CMB:.2e} J/m³")
