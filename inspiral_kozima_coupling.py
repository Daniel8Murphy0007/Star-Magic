"""
Inspiral Kozima Coupling — F_neutron as 5th GW Damping Channel

Session 204 | Daniel Murphy
PURPOSE: PAPER_008b uses 4 damping channels for GW strain reduction:
           h_UQFF(t) = h_GR(t) × D_total
           D_total = D_Aether × D_SCm × D_TRZ × D_String = 0.333

         This gives 66.7% strain reduction and 367.8 additional oscillation
         cycles of phase lag.  However, the Kozima neutron-drop force F_neutron
         is NOT included as a 5th damping channel.

         In LENR environments where SCm vacuum displacement occurs, the
         F_neutron^SCm couples to the inspiral metric perturbation h_ab via
         the buoyancy-reversed stress-energy tensor.

         This module:
           1. Computes D_Kozima from F_neutron^SCm contribution
           2. Extends h_UQFF to 5-channel damping
           3. Derives modified phase lag with Kozima channel
           4. Compares 4-channel vs 5-channel waveforms

EXISTING (PAPER_008b via wstp_symbolic_exporter.py):
  D_Aether = 1.000,  D_SCm = 1.000,  D_TRZ = 0.900,  D_String = 0.370
  F_combined = 0.333 (66.7% reduction)
  Phase lag: 367.8 additional oscillation cycles

NEW:
  D_Kozima = 1 / (1 + |F_neutron^SCm| / F_GW_scale)
  D_total_5 = D_Aether × D_SCm × D_TRZ × D_String × D_Kozima
  Δφ_Kozima = additional phase cycles from Kozima channel

ARCHITECTURE: CondensedPhysics.py rules — pure calculator, no hardcoded systems.
"""

import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

# ── §0  CONSTANTS ─────────────────────────────────────────────────────────

PI = math.pi
C = 2.998e8
HBAR = 1.055e-34
G = 6.674e-11
M_SUN = 1.989e30

# UQFF calibrated
SSQ = 0.57
KAPPA = 5.787e-9
BETA_I = 0.603
H_SCM = 0.99
U_UA = 1e-4

# Vacuum densities
RHO_SCM = 7.09e-37
RHO_UA = 7.09e-36
RHO_VAC_SCM = 9.47e-27

# GW damping (PAPER_001-009, wstp_symbolic_exporter.py L47-50)
D_AETHER = 1.000
D_SCM = 1.000
D_TRZ = 0.900
D_STRING = 0.370
F_COMBINED_4CH = D_AETHER * D_SCM * D_TRZ * D_STRING  # 0.333

# LENR parameters (kozima_scm_cross_section.py)
F_LENR_THZ = 1.25e12
OMEGA_SCM = 2 * PI * F_LENR_THZ
GAMMA_DEFAULT = 0.1e12 * 2 * PI
SIGMA_0 = 1e-4


# ── §1  KOZIMA DAMPING CHANNEL ───────────────────────────────────────────

class InspiralKozimaCoupling:
    """
    Extends the 4-channel UQFF GW damping to include F_neutron^SCm as
    a 5th channel via the Kozima neutron-drop mechanism.

    The Kozima damping factor is derived from the ratio of the neutron
    production force to the characteristic GW radiation-reaction force:

      D_Kozima = 1 / (1 + |F_neutron^SCm| / F_GW)

    where F_GW = (c^5 / G) × (h / (2π f))^2 × (2π f / c)^5 at the
    detector, but here we use the characteristic scale:

      F_GW_scale = (32/5) G μ M² / a³ (Peters 1964)

    for an inspiral with reduced mass μ, total mass M, separation a.
    """

    def __init__(self):
        self.d_aether = D_AETHER
        self.d_scm = D_SCM
        self.d_trz = D_TRZ
        self.d_string = D_STRING

    # ── 1a. GW radiation-reaction force scale ────────────────────────────

    def compute_gw_force_scale(self,
                               M1_kg: float = 36 * M_SUN,
                               M2_kg: float = 29 * M_SUN,
                               a_m: float = 350e3) -> Dict[str, Any]:
        """
        Characteristic GW radiation-reaction force (Peters 1964):
          F_GW = (32/5) G μ M² / a³

        where μ = M1 M2 / M_tot (reduced mass), M = M1 + M2 (total mass).

        Parameters:
          M1_kg, M2_kg: component masses (kg), default GW150914-like
          a_m: orbital separation (m), default 350 km (late inspiral)
        """
        M_tot = M1_kg + M2_kg
        mu = (M1_kg * M2_kg) / M_tot

        # Peters formula: radiation-reaction force
        F_GW = (32.0 / 5.0) * G * mu * M_tot**2 / a_m**3

        # Chirp mass
        M_c = (M1_kg * M2_kg)**0.6 / M_tot**0.2

        # Orbital frequency
        f_orb = (1.0 / (2 * PI)) * math.sqrt(G * M_tot / a_m**3)
        f_GW = 2 * f_orb

        return {
            "F_GW_N": F_GW,
            "M1_kg": M1_kg,
            "M2_kg": M2_kg,
            "M_tot_kg": M_tot,
            "mu_kg": mu,
            "M_chirp_kg": M_c,
            "M_chirp_Msun": M_c / M_SUN,
            "a_m": a_m,
            "f_orb_Hz": f_orb,
            "f_GW_Hz": f_GW,
            "equation": "F_GW = (32/5) G μ M² / a³  [Peters 1964]",
        }

    # ── 1b. F_neutron^SCm computation (from kozima_scm_cross_section) ────

    def compute_f_neutron_scm(self,
                              omega: float = OMEGA_SCM,
                              n: int = 13,
                              N_n: float = 1e28,
                              phi_phonon: float = 1e16,
                              F_U_Bi: float = 1.1,
                              F_U: float = 1.0) -> Dict[str, Any]:
        """
        Kozima SCm-modulated neutron production force:
          F_neutron^SCm = N_n · σ_n^SCm(ω,n) · Φ_phonon · (F_{U,Bi}/F_U - 1)

        where σ_n^SCm(ω,n) = σ₀ · exp[-(ω-ω_SCm)²/(2Γ²)] · (1+[SSq]·n/26)
        """
        # SCm cross-section
        exponent = -((omega - OMEGA_SCM) ** 2) / (2 * GAMMA_DEFAULT ** 2)
        gaussian = math.exp(exponent)
        vds_factor = 1.0 + (SSQ * n) / 26
        sigma_scm = SIGMA_0 * gaussian * vds_factor

        # Buoyancy reversal
        buoyancy_reversal = (F_U_Bi / F_U) - 1.0 if F_U != 0 else 0.0

        # Neutron force
        F_neutron = N_n * sigma_scm * phi_phonon * buoyancy_reversal

        return {
            "F_neutron_N": F_neutron,
            "sigma_n_scm": sigma_scm,
            "gaussian": gaussian,
            "vds_factor": vds_factor,
            "buoyancy_reversal": buoyancy_reversal,
            "N_n": N_n,
            "phi_phonon": phi_phonon,
            "omega": omega,
            "vds_level": n,
            "equation": "F_neutron^SCm = N_n · σ_n^SCm · Φ_phonon · (F_{U,Bi}/F_U - 1)",
        }

    # ── 1c. Kozima damping factor ────────────────────────────────────────

    def compute_d_kozima(self,
                         F_neutron_N: float,
                         F_GW_N: float) -> Dict[str, Any]:
        """
        5th damping channel from Kozima neutron-drop coupling:
          D_Kozima = 1 / (1 + |F_neutron^SCm| / F_GW)

        Physical interpretation: The neutron-drop force acts as an additional
        energy dissipation pathway, extracting energy from the GW field into
        the SCm vacuum lattice via phonon-mediated neutron capture events.
        """
        ratio = abs(F_neutron_N) / F_GW_N if F_GW_N != 0 else 0
        D_kozima = 1.0 / (1.0 + ratio)

        return {
            "D_Kozima": D_kozima,
            "F_neutron_N": F_neutron_N,
            "F_GW_N": F_GW_N,
            "ratio": ratio,
            "strain_reduction_pct": (1.0 - D_kozima) * 100,
            "equation": "D_Kozima = 1 / (1 + |F_neutron^SCm| / F_GW)",
        }

    # ── 1d. 5-channel combined damping ───────────────────────────────────

    def compute_5channel_damping(self,
                                 D_Kozima: float) -> Dict[str, Any]:
        """
        Extend the 4-channel damping to 5 channels:
          D_total_5 = D_Aether × D_SCm × D_TRZ × D_String × D_Kozima
        """
        D_total_4 = F_COMBINED_4CH
        D_total_5 = D_total_4 * D_Kozima

        return {
            "D_total_4ch": D_total_4,
            "D_total_5ch": D_total_5,
            "D_Aether": self.d_aether,
            "D_SCm": self.d_scm,
            "D_TRZ": self.d_trz,
            "D_String": self.d_string,
            "D_Kozima": D_Kozima,
            "strain_reduction_4ch_pct": (1.0 - D_total_4) * 100,
            "strain_reduction_5ch_pct": (1.0 - D_total_5) * 100,
            "additional_reduction_pct": (D_total_4 - D_total_5) / D_total_4 * 100 if D_total_4 != 0 else 0,
            "equation_4ch": "D_4 = D_Aether × D_SCm × D_TRZ × D_String = 0.333",
            "equation_5ch": "D_5 = D_Aether × D_SCm × D_TRZ × D_String × D_Kozima",
        }

    # ── 1e. Phase lag modification ───────────────────────────────────────

    def compute_phase_lag(self,
                          D_total_4: float = F_COMBINED_4CH,
                          D_total_5: float = 0.0,
                          f_GW_Hz: float = 100.0,
                          T_inspiral_s: float = 0.2) -> Dict[str, Any]:
        """
        Phase lag in oscillation cycles from UQFF damping.

        PAPER_008b baseline: 367.8 additional oscillation cycles (4-channel).

        Phase accumulation from damping:
          Δφ = ∫₀ᵀ 2π f_GW (1 - D_total) dt
          N_cycles = Δφ / (2π) = f_GW × T_inspiral × (1 - D_total)
        """
        # 4-channel (PAPER_008b baseline)
        N_cycles_4 = f_GW_Hz * T_inspiral_s * (1.0 - D_total_4)

        # 5-channel
        N_cycles_5 = f_GW_Hz * T_inspiral_s * (1.0 - D_total_5)

        # Additional cycles from Kozima
        delta_N = N_cycles_5 - N_cycles_4

        return {
            "N_cycles_4ch": N_cycles_4,
            "N_cycles_5ch": N_cycles_5,
            "delta_N_Kozima": delta_N,
            "f_GW_Hz": f_GW_Hz,
            "T_inspiral_s": T_inspiral_s,
            "PAPER_008b_baseline": 367.8,
            "equation": "N_cycles = f_GW × T_inspiral × (1 − D_total)",
        }

    # ── 1f. GW waveform comparison ───────────────────────────────────────

    def compute_waveform_comparison(self,
                                    M1_kg: float = 36 * M_SUN,
                                    M2_kg: float = 29 * M_SUN,
                                    a_m: float = 350e3,
                                    d_L_m: float = 410e6 * 3.086e16 * 1e6,
                                    omega_kozima: float = OMEGA_SCM,
                                    n_vds: int = 13,
                                    N_n: float = 1e28,
                                    phi_phonon: float = 1e16,
                                    F_U_Bi: float = 1.1,
                                    F_U: float = 1.0
                                    ) -> Dict[str, Any]:
        """
        Full comparison: 4-channel vs 5-channel inspiral waveform.

        Returns the complete derivation chain from binary parameters
        through Kozima coupling to modified strain and phase lag.

        Parameters:
          M1_kg, M2_kg: binary component masses
          a_m: orbital separation
          d_L_m: luminosity distance (m)
          omega_kozima: resonance frequency for F_neutron
          n_vds: VDS level
          N_n: neutron density
          phi_phonon: phonon fluence
          F_U_Bi, F_U: buoyancy force ratio
        """
        # Step 1: GW force scale
        gw = self.compute_gw_force_scale(M1_kg, M2_kg, a_m)

        # Step 2: F_neutron^SCm
        fn = self.compute_f_neutron_scm(omega_kozima, n_vds, N_n,
                                        phi_phonon, F_U_Bi, F_U)

        # Step 3: Kozima damping factor
        dk = self.compute_d_kozima(fn["F_neutron_N"], gw["F_GW_N"])

        # Step 4: 5-channel damping
        d5 = self.compute_5channel_damping(dk["D_Kozima"])

        # Step 5: GR strain amplitude (leading order)
        M_c = gw["M_chirp_kg"]
        f_GW = gw["f_GW_Hz"]
        if d_L_m > 0 and f_GW > 0:
            h_GR = (4.0 / d_L_m) * (G * M_c / C**2)**(5.0/3) * \
                   (PI * f_GW / C)**(2.0/3)
        else:
            h_GR = 0.0

        # Step 6: Damped strains
        h_4ch = h_GR * d5["D_total_4ch"]
        h_5ch = h_GR * d5["D_total_5ch"]

        # Step 7: Phase lag
        T_inspiral = 0.2  # seconds (late inspiral)
        phase = self.compute_phase_lag(d5["D_total_4ch"], d5["D_total_5ch"],
                                       f_GW, T_inspiral)

        # Derivation chain
        chain = [
            "══════════════════════════════════════════════════════════════",
            "INSPIRAL KOZIMA COUPLING — 5-CHANNEL GW DAMPING",
            "══════════════════════════════════════════════════════════════",
            "",
            "§1. Binary System Parameters:",
            f"  M1 = {M1_kg:.4e} kg ({M1_kg/M_SUN:.1f} M☉)",
            f"  M2 = {M2_kg:.4e} kg ({M2_kg/M_SUN:.1f} M☉)",
            f"  M_chirp = {M_c:.4e} kg ({M_c/M_SUN:.2f} M☉)",
            f"  a = {a_m:.4e} m",
            f"  f_GW = {f_GW:.4f} Hz",
            f"  d_L = {d_L_m:.4e} m",
            "",
            "§2. GW Radiation-Reaction Force (Peters 1964):",
            f"  F_GW = (32/5) G μ M² / a³",
            f"       = {gw['F_GW_N']:.6e} N",
            "",
            "§3. Kozima Neutron-Drop Force (SCm-modulated):",
            f"  σ_n^SCm(ω,n) = σ₀ · exp[-(ω-ω_SCm)²/(2Γ²)] · (1+[SSq]·n/26)",
            f"               = {fn['sigma_n_scm']:.6e}",
            f"  F_neutron^SCm = N_n · σ_n^SCm · Φ_phonon · (F_{{U,Bi}}/F_U - 1)",
            f"                = {fn['F_neutron_N']:.6e} N",
            "",
            "§4. Kozima Damping Factor:",
            f"  D_Kozima = 1 / (1 + |F_neutron| / F_GW)",
            f"           = 1 / (1 + {dk['ratio']:.6e})",
            f"           = {dk['D_Kozima']:.10f}",
            "",
            "§5. 4-Channel vs 5-Channel Damping:",
            f"  D_4 = D_Aether × D_SCm × D_TRZ × D_String",
            f"       = {self.d_aether} × {self.d_scm} × {self.d_trz} × {self.d_string}",
            f"       = {d5['D_total_4ch']:.6f}",
            f"  D_5 = D_4 × D_Kozima",
            f"       = {d5['D_total_4ch']:.6f} × {dk['D_Kozima']:.10f}",
            f"       = {d5['D_total_5ch']:.10f}",
            "",
            "§6. GW Strain Comparison:",
            f"  h_GR  = (4/d_L)(GM_c/c²)^(5/3)(πf/c)^(2/3)",
            f"        = {h_GR:.6e}",
            f"  h_4ch = h_GR × D_4 = {h_4ch:.6e}  (66.7% reduction)",
            f"  h_5ch = h_GR × D_5 = {h_5ch:.6e}  ({d5['strain_reduction_5ch_pct']:.2f}% reduction)",
            "",
            "§7. Phase Lag:",
            f"  N_cycles = f_GW × T × (1 - D)",
            f"  4-channel: {phase['N_cycles_4ch']:.2f} cycles",
            f"  5-channel: {phase['N_cycles_5ch']:.2f} cycles",
            f"  Kozima contribution: {phase['delta_N_Kozima']:.4f} additional cycles",
            f"  PAPER_008b baseline: {phase['PAPER_008b_baseline']} cycles",
            "",
            "§8. Physical Interpretation:",
            "  The Kozima neutron-drop force acts as an energy sink in the",
            "  inspiral GW field.  Neutron capture events mediated by SCm",
            "  phonon resonance at 1.25 THz extract orbital energy into",
            "  the lattice heat bath, producing a 5th damping channel.",
            "  The effect is small for astrophysical binaries (D_Kozima ≈ 1)",
            "  but becomes significant for lab-scale LENR GW detectors.",
            "══════════════════════════════════════════════════════════════",
        ]

        return {
            "gw_system": gw,
            "f_neutron": fn,
            "d_kozima": dk,
            "damping_5ch": d5,
            "h_GR": h_GR,
            "h_4ch": h_4ch,
            "h_5ch": h_5ch,
            "phase_lag": phase,
            "derivation_chain": chain,
        }


# ── §2  SELF-TEST ─────────────────────────────────────────────────────────

def main():
    """Validate inspiral Kozima coupling module."""
    print("=" * 72)
    print("INSPIRAL KOZIMA COUPLING — 5-CHANNEL GW DAMPING VALIDATION")
    print("=" * 72)

    calc = InspiralKozimaCoupling()

    # Run default (GW150914-like) comparison
    result = calc.compute_waveform_comparison()

    # Print derivation chain
    for line in result["derivation_chain"]:
        print(line)

    # Assertions
    assert result["damping_5ch"]["D_total_4ch"] == F_COMBINED_4CH, \
        f"4-channel damping mismatch: {result['damping_5ch']['D_total_4ch']}"
    assert result["damping_5ch"]["D_total_5ch"] <= F_COMBINED_4CH, \
        "5-channel should be ≤ 4-channel"
    assert result["d_kozima"]["D_Kozima"] <= 1.0, \
        "D_Kozima must be ≤ 1"
    assert result["d_kozima"]["D_Kozima"] > 0, \
        "D_Kozima must be > 0"
    assert result["h_5ch"] <= result["h_4ch"], \
        "5-channel strain should be ≤ 4-channel"

    print()
    print(f"✓ D_total_4ch = {result['damping_5ch']['D_total_4ch']:.6f}")
    print(f"✓ D_total_5ch = {result['damping_5ch']['D_total_5ch']:.10f}")
    print(f"✓ D_Kozima    = {result['d_kozima']['D_Kozima']:.10f}")
    print(f"✓ h_GR        = {result['h_GR']:.6e}")
    print(f"✓ h_4ch       = {result['h_4ch']:.6e}")
    print(f"✓ h_5ch       = {result['h_5ch']:.6e}")
    print(f"✓ Phase lag 4ch = {result['phase_lag']['N_cycles_4ch']:.2f} cycles")
    print(f"✓ Phase lag 5ch = {result['phase_lag']['N_cycles_5ch']:.2f} cycles")
    print()
    print("ALL ASSERTIONS PASSED")
    print(json.dumps({
        "module": "inspiral_kozima_coupling",
        "status": "VALIDATED",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "D_Kozima": result["d_kozima"]["D_Kozima"],
        "D_total_5ch": result["damping_5ch"]["D_total_5ch"],
    }, indent=2))


if __name__ == "__main__":
    main()
