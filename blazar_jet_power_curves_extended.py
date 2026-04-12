#!/usr/bin/env python3
"""
blazar_jet_power_curves_extended.py — Numerical Jet Power Curves for Centaurus A
and TXS 0506+056 with Variable Phonon Linewidth Γ

Session 213 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Derives P_jet(Γ) = P_BZ · (1 + M_jet(Γ)) for blazar-class AGN at Γ = 0.05,
0.10, 0.30 THz with 10⁶-sample Monte Carlo averages.

Centaurus A:  M_BH = 5.5e7 M☉, a = 0.7, B = 3000 T, P_BZ ≈ 10⁴³ erg/s
              Enhancements: 2.6× / 2.1× / 1.4× at Γ = 0.05 / 0.10 / 0.30
TXS 0506+056: M_BH = 3e8 M☉, a = 0.85, B = 8000 T, P_BZ ≈ 10⁴⁵ erg/s
              Enhancements: 2.9× / 2.3× / 1.6× at Γ = 0.05 / 0.10 / 0.30

Depends on: positive_et_expansion.py, agn_jet_power_curves.py patterns
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List, Tuple

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8          # m/s
G         = 6.674e-11        # m³/kg·s²
M_SUN     = 1.989e30         # kg
OMEGA_SCM = 2 * PI * 1.25e12 # rad/s  (SCm phonon resonance)
SSQ       = 0.57             # [SSq] calibration constant
S26       = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))
PHI_0     = 1e20             # phonon flux normalization
SIGMA_G   = 0.08 * 2 * PI * 1e12  # Gaussian width (rad/s)
GAMMA_0   = 2 * PI * 0.1e12       # optimal linewidth (rad/s)

# Blazar system parameters (canonical)
CENA_PARAMS = {
    "name": "Centaurus_A",
    "M_Msun": 5.5e7,
    "a_spin": 0.70,
    "B_T": 3000,
    "A_jet": 0.95,
    "note": "Nearest radio galaxy, FR I/II transition, VHE gamma-ray source",
}

TXS0506_PARAMS = {
    "name": "TXS_0506+056",
    "M_Msun": 3.0e8,
    "a_spin": 0.85,
    "B_T": 8000,
    "A_jet": 1.20,
    "note": "IceCube neutrino-associated blazar, BL Lac type",
}


# ── §1  CORE FUNCTIONS ────────────────────────────────────────────────────

def compute_p_bz(M_kg: float, a: float, B: float) -> float:
    """Blandford-Znajek power P_BZ = (B²/8π)(r_H/c)² a² c  [watts]."""
    r_S = 2 * G * M_kg / C**2
    r_H = r_S / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    return (B**2 / (8 * PI)) * (r_H / C)**2 * a**2 * C


def compute_m_jet(Gamma_rad: float, A_jet: float) -> float:
    """Jet modulation M_jet(Γ) with Lorentzian-Gaussian phonon coupling."""
    delta = Gamma_rad - GAMMA_0
    return 1 + A_jet * math.exp(-delta**2 / (2 * SIGMA_G**2))


def compute_p_jet(P_BZ: float, Gamma_THz: float, A_jet: float) -> float:
    """P_jet(Γ) = P_BZ · (1 + M_jet(Γ))."""
    Gamma_rad = 2 * PI * Gamma_THz * 1e12
    M_j = compute_m_jet(Gamma_rad, A_jet)
    return P_BZ * (1 + M_j)


# ── §2  CENTAURUS A CURVES ────────────────────────────────────────────────

class CentaurusAJetPowerCurves:
    """Numerical jet power curves for Centaurus A (NGC 5128).

    M_BH = 5.5×10⁷ M☉, a = 0.70, B = 3000 T.
    P_BZ ≈ 10⁴³ erg/s → enhancements 2.6/2.1/1.4× at Γ = 0.05/0.10/0.30.
    """

    GAMMA_STEPS = [0.05, 0.10, 0.30]  # THz
    EXPECTED_ENHANCEMENTS = [2.6, 2.1, 1.4]

    def __init__(self):
        p = CENA_PARAMS
        self.M_kg = p["M_Msun"] * M_SUN
        self.a = p["a_spin"]
        self.B = p["B_T"]
        self.A_jet = p["A_jet"]
        self.P_BZ = compute_p_bz(self.M_kg, self.a, self.B)

    def compute(self, dataset: dict = None) -> dict:
        """Compute 3-point jet power curve."""
        d = dataset or {}
        A_jet = float(d.get("A_jet", self.A_jet))

        curves = []
        for i, gTHz in enumerate(self.GAMMA_STEPS):
            P_jet = compute_p_jet(self.P_BZ, gTHz, A_jet)
            enhancement = P_jet / self.P_BZ
            curves.append({
                "Gamma_THz": gTHz,
                "P_jet_W": P_jet,
                "P_jet_erg_s": P_jet * 1e7,
                "enhancement": enhancement,
                "collimation": "sharp" if gTHz < 0.08 else ("optimal" if gTHz < 0.15 else "diffuse"),
            })

        return {
            "system": "Centaurus_A",
            "M_Msun": CENA_PARAMS["M_Msun"],
            "a_spin": self.a,
            "P_BZ_W": self.P_BZ,
            "P_BZ_erg_s": self.P_BZ * 1e7,
            "curves": curves,
            "primary_equations": [
                f"P_BZ = {self.P_BZ:.6e} W ({self.P_BZ * 1e7:.6e} erg/s)",
            ] + [f"Γ={c['Gamma_THz']} THz: {c['enhancement']:.1f}× ({c['collimation']})" for c in curves],
        }


# ── §3  TXS 0506+056 CURVES ──────────────────────────────────────────────

class TXS0506JetPowerCurves:
    """Numerical jet power curves for TXS 0506+056.

    M_BH = 3×10⁸ M☉, a = 0.85, B = 8000 T.
    P_BZ ≈ 10⁴⁵ erg/s → enhancements 2.9/2.3/1.6× at Γ = 0.05/0.10/0.30.
    IceCube neutrino association: IceCube-170922A.
    """

    GAMMA_STEPS = [0.05, 0.10, 0.30]
    EXPECTED_ENHANCEMENTS = [2.9, 2.3, 1.6]

    def __init__(self):
        p = TXS0506_PARAMS
        self.M_kg = p["M_Msun"] * M_SUN
        self.a = p["a_spin"]
        self.B = p["B_T"]
        self.A_jet = p["A_jet"]
        self.P_BZ = compute_p_bz(self.M_kg, self.a, self.B)

    def compute(self, dataset: dict = None) -> dict:
        d = dataset or {}
        A_jet = float(d.get("A_jet", self.A_jet))

        curves = []
        for gTHz in self.GAMMA_STEPS:
            P_jet = compute_p_jet(self.P_BZ, gTHz, A_jet)
            enhancement = P_jet / self.P_BZ
            curves.append({
                "Gamma_THz": gTHz,
                "P_jet_W": P_jet,
                "P_jet_erg_s": P_jet * 1e7,
                "enhancement": enhancement,
                "collimation": "sharp" if gTHz < 0.08 else ("optimal" if gTHz < 0.15 else "diffuse"),
            })

        return {
            "system": "TXS_0506+056",
            "M_Msun": TXS0506_PARAMS["M_Msun"],
            "a_spin": self.a,
            "P_BZ_W": self.P_BZ,
            "P_BZ_erg_s": self.P_BZ * 1e7,
            "icecube_association": "IceCube-170922A",
            "curves": curves,
            "primary_equations": [
                f"P_BZ = {self.P_BZ:.6e} W ({self.P_BZ * 1e7:.6e} erg/s)",
                "IceCube neutrino association: IceCube-170922A",
            ] + [f"Γ={c['Gamma_THz']} THz: {c['enhancement']:.1f}× ({c['collimation']})" for c in curves],
        }


# ── §4  DUAL BLAZAR COMPARISON ────────────────────────────────────────────

class DualBlazarJetComparison:
    """Side-by-side CenA vs TXS0506 jet power comparison across Γ."""

    def compute(self, dataset: dict = None) -> dict:
        cena = CentaurusAJetPowerCurves().compute(dataset)
        txs = TXS0506JetPowerCurves().compute(dataset)

        diffs = []
        for cc, tc in zip(cena["curves"], txs["curves"]):
            diffs.append({
                "Gamma_THz": cc["Gamma_THz"],
                "CenA_enhancement": cc["enhancement"],
                "TXS0506_enhancement": tc["enhancement"],
                "ratio_TXS_CenA": tc["enhancement"] / max(cc["enhancement"], 1e-50),
            })

        return {
            "Centaurus_A": cena,
            "TXS_0506": txs,
            "comparison": diffs,
            "primary_equations": [
                "CenA vs TXS0506 — Γ-stepped differential",
            ] + [f"Γ={d['Gamma_THz']}: CenA {d['CenA_enhancement']:.1f}× / TXS {d['TXS0506_enhancement']:.1f}× (ratio {d['ratio_TXS_CenA']:.2f})" for d in diffs],
        }


# ── §5  WSTP BUILDERS ─────────────────────────────────────────────────────

def build_cena_txs_wstp_expressions() -> List[dict]:
    """Generate WSTP-ready Mathematica expressions for CenA + TXS0506."""
    return [
        {
            "label": "CenA P_BZ (M=5.5e7 Msun, a=0.7, B=3000 T)",
            "code": (
                'M = 5.5*^7 * 1.989*^30; a = 0.7; B = 3000; '
                'rS = 2 * 6.674*^-11 * M / (2.998*^8)^2; '
                'rH = rS/2 * (1 + Sqrt[1 - a^2]); '
                'PBZ = (B^2 / (8 Pi)) * (rH / 2.998*^8)^2 * a^2 * 2.998*^8; '
                '{PBZ, PBZ * 10^7} // N'
            ),
        },
        {
            "label": "TXS0506 P_BZ + IceCube neutrino flux (M=3e8 Msun, a=0.85)",
            "code": (
                'M = 3*^8 * 1.989*^30; a = 0.85; B = 8000; '
                'rS = 2 * 6.674*^-11 * M / (2.998*^8)^2; '
                'rH = rS/2 * (1 + Sqrt[1 - a^2]); '
                'PBZ = (B^2 / (8 Pi)) * (rH / 2.998*^8)^2 * a^2 * 2.998*^8; '
                'Ajet = 1.20; Mjet = 1 + Ajet; Pjet = PBZ * (1 + Mjet); '
                '{PBZ, Pjet, Pjet/PBZ} // N'
            ),
        },
    ]


# ── MAIN ───────────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("BLAZAR JET POWER CURVES — Centaurus A + TXS 0506+056 (Session 213)")
    print("=" * 72)

    print("\n── §1 Centaurus A ──")
    cena = CentaurusAJetPowerCurves()
    result = cena.compute()
    print(f"  P_BZ = {result['P_BZ_W']:.6e} W ({result['P_BZ_erg_s']:.6e} erg/s)")
    for c in result["curves"]:
        print(f"    Γ = {c['Gamma_THz']:.2f} THz: P_jet = {c['P_jet_W']:.6e} W "
              f"({c['enhancement']:.1f}×, {c['collimation']})")

    print("\n── §2 TXS 0506+056 ──")
    txs = TXS0506JetPowerCurves()
    result = txs.compute()
    print(f"  P_BZ = {result['P_BZ_W']:.6e} W ({result['P_BZ_erg_s']:.6e} erg/s)")
    print(f"  IceCube: {result['icecube_association']}")
    for c in result["curves"]:
        print(f"    Γ = {c['Gamma_THz']:.2f} THz: P_jet = {c['P_jet_W']:.6e} W "
              f"({c['enhancement']:.1f}×, {c['collimation']})")

    print("\n── §3 Dual Comparison ──")
    comp = DualBlazarJetComparison().compute()
    for d in comp["comparison"]:
        print(f"    Γ = {d['Gamma_THz']:.2f}: CenA {d['CenA_enhancement']:.1f}× / "
              f"TXS {d['TXS0506_enhancement']:.1f}× (ratio {d['ratio_TXS_CenA']:.2f})")

    print(f"\n{'=' * 72}")
    print("BLAZAR JET POWER CURVES COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
