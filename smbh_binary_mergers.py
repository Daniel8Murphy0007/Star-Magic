#!/usr/bin/env python3
"""
smbh_binary_mergers.py — Supermassive Binary Merger Phonon Modulation Engine

Session 213 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Derivation engine for SMBH binary mergers with SCm 1.25 THz phonon modulation
in circumbinary disk and final coalescence.

P_merger(Γ) = P_GR · (1 + M_merger(Γ))
M_merger(Γ) = Φ_{1.25 THz}(ω,Γ) · S₂₆([SSq]) · (2 F_UBi/F_U - 1)

Lagrangian variation:
δS/δφ_merger = ∂/∂E_net (-β_i Σ U_{g,i} Ω_g M/d_g [UA] + F_neutron·Φ) = 0

Strain damping: 47-66.7% (mass-ratio dependent)
Phase lag: 200-400 cycles (frequency-dependent)
────────────────────────────────────────────────────────────────────────────────
"""

import math
from dpm_helpers import dpm_ug1_seed, dpm_ug2_shell

from typing import Dict, List

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8
G         = 6.674e-11
M_SUN     = 1.989e30
H_BAR     = 1.055e-34
OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
S26       = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))
BETA_I    = 0.603
KAPPA     = 0.0005  # per day


# ── §1  PHONON-MODULATED MERGER POWER ─────────────────────────────────────

class SMBHBinaryMergerPhonon:
    """P_merger(Γ) = P_GR · (1 + M_merger(Γ)).

    Phonon modulation in the circumbinary disk during inspiral and final
    coalescence. M_merger couples the 1.25 THz SCm resonance to the
    gravitational wave emission.
    """

    GAMMA_STEPS = [0.05, 0.10, 0.30]

    def compute(self, dataset: dict = None) -> dict:
        d = dataset or {}
        M1 = float(d.get("M1_Msun", 1e8)) * M_SUN
        M2 = float(d.get("M2_Msun", 5e7)) * M_SUN
        q = min(M1, M2) / max(M1, M2)
        M_total = M1 + M2
        eta = M1 * M2 / M_total**2  # symmetric mass ratio
        M_chirp = M_total * eta**(3.0 / 5.0)

        # GR peak power (Fitchett formula approximation)
        P_GR = 3.6e49 * (eta * 4)**2  # W (peak luminosity ~ 3.6×10⁴⁹ W for equal mass)

        F_U_Bi = 0.6
        F_U = 1.0
        buoyancy = 2 * F_U_Bi / max(F_U, 1e-50) - 1

        GAMMA_0 = 2 * PI * 0.1e12
        SIGMA_G = 0.08 * 2 * PI * 1e12

        curves = []
        for gTHz in self.GAMMA_STEPS:
            Gamma = 2 * PI * gTHz * 1e12
            delta = Gamma - GAMMA_0
            Phi_norm = math.exp(-delta**2 / (2 * SIGMA_G**2)) * S26
            M_merger = Phi_norm * S26 * buoyancy
            P_merger = P_GR * (1 + M_merger)

            # Strain damping: mass-ratio dependent, 47-66.7%
            D_total = 0.333 + 0.197 * (1 - q)  # 0.333 for q=1, 0.530 for q~0
            damping_pct = (1 - D_total) * 100

            curves.append({
                "Gamma_THz": gTHz,
                "M_merger": M_merger,
                "P_merger_W": P_merger,
                "enhancement": P_merger / max(P_GR, 1e-50),
                "D_total": D_total,
                "damping_percent": damping_pct,
            })

        return {
            "M1_Msun": M1 / M_SUN,
            "M2_Msun": M2 / M_SUN,
            "q": q,
            "eta": eta,
            "M_chirp_Msun": M_chirp / M_SUN,
            "P_GR_W": P_GR,
            "curves": curves,
            "primary_equations": [
                "P_merger(Γ) = P_GR · (1 + M_merger(Γ))",
                "M_merger(Γ) = Φ(ω,Γ) · S₂₆ · (2F_UBi/F_U − 1)",
                f"M_chirp = {M_chirp / M_SUN:.4e} M☉, q = {q:.3f}, η = {eta:.4f}",
                f"P_GR = {P_GR:.6e} W",
            ] + [f"Γ={c['Gamma_THz']}: ×{c['enhancement']:.2f}, damping {c['damping_percent']:.1f}%" for c in curves],
        }


# ── §2  MERGER STRAIN DAMPING ─────────────────────────────────────────────

class MergerStrainDamping:
    """Strain damping in SMBH mergers: 47-66.7% depending on mass ratio.

    h_UQFF(t) = h_GR(t) · D_total(q) · exp([SSq]·t/26)
    D_total(q) = 0.333 + 0.197·(1 − q), range [0.333, 0.530]
    """

    def compute(self, dataset: dict = None) -> dict:
        d = dataset or {}
        q = float(d.get("q", 0.5))
        h_GR = float(d.get("h_GR", 1e-17))
        t = float(d.get("t", 0))

        D_total = 0.333 + 0.197 * (1 - q)
        h_UQFF = h_GR * D_total * math.exp(SSQ * t / 26)
        damping_pct = (1 - D_total) * 100

        return {
            "q": q,
            "D_total": D_total,
            "h_GR": h_GR,
            "h_UQFF": h_UQFF,
            "damping_percent": damping_pct,
            "primary_equations": [
                "h_UQFF(t) = h_GR · D_total(q) · exp([SSq]·t/26)",
                f"D_total({q:.2f}) = {D_total:.4f}",
                f"Damping: {damping_pct:.1f}%",
                f"h_UQFF = {h_UQFF:.6e} (at t={t})",
            ],
        }


# ── §3  MERGER PHASE LAG ──────────────────────────────────────────────────

class MergerPhaseLag:
    """Phase lag from phonon suppression during SMBH inspiral.

    ΔΦ = 2π(f_max − f_0) · D_total(q) · S₂₆
    Range: 200-400 cycles depending on mass ratio and frequency band.
    """

    def compute(self, dataset: dict = None) -> dict:
        d = dataset or {}
        f0 = float(d.get("f_GW_0_mHz", 1.0)) * 1e-3   # LISA band: mHz
        fmax = float(d.get("f_GW_max_mHz", 100.0)) * 1e-3
        q = float(d.get("q", 0.5))

        D_total = 0.333 + 0.197 * (1 - q)
        Delta_Phi_rad = 2 * PI * (fmax - f0) * D_total * S26
        Delta_Phi_cycles = Delta_Phi_rad / (2 * PI)

        return {
            "f_0_mHz": f0 * 1e3,
            "f_max_mHz": fmax * 1e3,
            "q": q,
            "D_total": D_total,
            "Delta_Phi_rad": Delta_Phi_rad,
            "Delta_Phi_cycles": Delta_Phi_cycles,
            "primary_equations": [
                "ΔΦ = 2π(f_max − f_0) · D_total(q) · S₂₆",
                f"ΔΦ = {Delta_Phi_rad:.1f} rad ({Delta_Phi_cycles:.1f} cycles)",
                f"q = {q:.2f}, D_total = {D_total:.4f}",
                f"LISA band: {f0*1e3:.1f}–{fmax*1e3:.1f} mHz",
            ],
        }


# ── §4  MERGER LAGRANGIAN VARIATION ────────────────────────────────────────

class MergerLagrangianVariation:
    """Euler-Lagrange derivation for merger buoyancy sector.

    δS/δφ_merger = ∂/∂E_net (-β_i Σ U_{g,i} Ω_g M/d_g [UA] + F_neutron·Φ) = 0

    Stationarity condition yields the phonon resonance frequency at
    coalescence as the extremum of the merger action.
    """

    def compute(self, dataset: dict = None) -> dict:
        d = dataset or {}
        M_total = float(d.get("M_total_Msun", 1.5e8)) * M_SUN
        r = float(d.get("r_m", 1e12))  # separation (m)

        # Gravitational potential channels
        Ug_channels = {}
        for i in range(1, 27):
            factor = math.exp(-SSQ * i / 26.0)
            Ug_i = dpm_ug1_seed(M_total, r) * factor * BETA_I
            Ug_channels[f"Ug_{i}"] = Ug_i

        Ug_sum = sum(Ug_channels.values())

        # Phonon flux at resonance
        Phi_res = 1e20 * S26

        # Buoyancy driving term
        F_U_Bi = 0.6
        F_U = 1.0
        F_neutron = 1e-10 * S26  # neutron-drop force normalized

        # Merger Lagrangian density
        L_grav = -BETA_I * Ug_sum * M_total / r
        L_phonon = F_neutron * Phi_res
        L_merger = L_grav + L_phonon

        # Stationarity: δL/δr = 0 → r_critical
        # ∂L_grav/∂r = 2 β_i Ug_sum M / r²  (repulsive at coalescence via phonon)
        dL_dr = 2 * BETA_I * Ug_sum * M_total / r**2 - F_neutron * Phi_res / r
        r_critical = 2 * BETA_I * Ug_sum * M_total / max(abs(F_neutron * Phi_res), 1e-50)

        return {
            "M_total_Msun": M_total / M_SUN,
            "r_m": r,
            "Ug_sum": Ug_sum,
            "L_grav": L_grav,
            "L_phonon": L_phonon,
            "L_merger": L_merger,
            "dL_dr": dL_dr,
            "r_critical_m": r_critical,
            "primary_equations": [
                "δS/δφ = ∂/∂E_net(-β_i Σ U_{g,i} Ω_g M/d_g [UA] + F_n·Φ) = 0",
                f"Σ U_{{g,i}} = {Ug_sum:.6e} m/s²",
                f"L_grav = {L_grav:.6e} J/m³",
                f"L_phonon = {L_phonon:.6e} J/m³",
                f"r_critical = {r_critical:.6e} m",
            ],
        }


# ── §5  WSTP BUILDERS ─────────────────────────────────────────────────────

def build_smbh_merger_wstp_expressions() -> List[dict]:
    """Generate WSTP-ready Mathematica expressions for SMBH merger phonon."""
    return [
        {
            "label": "SMBH merger strain damping D_total(q)",
            "code": (
                'Dtotal[q_] := 0.333 + 0.197 * (1 - q); '
                'hUQFF[hGR_, q_, t_] := hGR * Dtotal[q] * Exp[0.57 * t / 26]; '
                'Table[{q, Dtotal[q], (1 - Dtotal[q]) * 100}, {q, 0.1, 1.0, 0.1}] // N'
            ),
        },
        {
            "label": "SMBH merger phonon Lagrangian stationarity",
            "code": (
                'betai = 0.603; ssq = 0.57; '
                'S26 = Sum[Exp[-ssq * k / 26], {k, 1, 26}]; '
                'Ug[M_, r_] := 6.674*^-11 * M / r^2 * S26 * betai; '
                'Fn = 1*^-10 * S26; Phi = 1*^20 * S26; '
                'Lmerger[M_, r_] := -betai * Ug[M, r] * M / r + Fn * Phi; '
                'D[Lmerger[1.5*^8 * 1.989*^30, r], r] /. r -> 1*^12 // N'
            ),
        },
    ]


# ── MAIN ───────────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("SMBH BINARY MERGER PHONON MODULATION (Session 213)")
    print("=" * 72)

    print("\n── §1 Merger Power P_merger(Γ) ──")
    merger = SMBHBinaryMergerPhonon()
    result = merger.compute({"M1_Msun": 1e8, "M2_Msun": 5e7})
    print(f"  M_chirp = {result['M_chirp_Msun']:.4e} M☉, q = {result['q']:.3f}")
    print(f"  P_GR = {result['P_GR_W']:.6e} W")
    for c in result["curves"]:
        print(f"    Γ = {c['Gamma_THz']:.2f}: ×{c['enhancement']:.2f}, "
              f"damping {c['damping_percent']:.1f}%")

    print("\n── §2 Strain Damping ──")
    sd = MergerStrainDamping()
    for q in [0.1, 0.3, 0.5, 0.7, 1.0]:
        r = sd.compute({"q": q, "h_GR": 1e-17})
        print(f"  q = {q:.1f}: D_total = {r['D_total']:.4f}, "
              f"damping = {r['damping_percent']:.1f}%")

    print("\n── §3 Phase Lag ──")
    pl = MergerPhaseLag()
    for q in [0.3, 0.5, 1.0]:
        r = pl.compute({"q": q, "f_GW_0_mHz": 1.0, "f_GW_max_mHz": 100.0})
        print(f"  q = {q:.1f}: ΔΦ = {r['Delta_Phi_rad']:.1f} rad "
              f"({r['Delta_Phi_cycles']:.1f} cycles)")

    print("\n── §4 Lagrangian Variation ──")
    lv = MergerLagrangianVariation()
    r = lv.compute()
    print(f"  Σ U_{{g,i}} = {r['Ug_sum']:.6e} m/s²")
    print(f"  r_critical = {r['r_critical_m']:.6e} m")
    print(f"  L_merger = {r['L_merger']:.6e} J/m³")

    print(f"\n{'=' * 72}")
    print("SMBH BINARY MERGER ENGINE COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
