#!/usr/bin/env python3
"""
linewidth_jet_modulation.py — Phonon Linewidth Γ Variation Engine for Jet
Power, Collimation, and Sharpness Mapping

Session 213 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Standalone engine mapping Γ → {collimation, sharpness, power modulation} for
observational interpretation. Systematic sweep from 0.01 to 1.0 THz with
classification into narrow/optimal/broad regimes.

M_jet(Γ) = exp(-(ω - ω_SCm)² / (2Γ²)) · S₂₆([SSq]) · (2F_UBi/F_U - 1)
Collimation ∝ Q = ω_SCm / (2Γ)  — sharp → high Q, diffuse → low Q
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
S26       = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))
GAMMA_0   = 2 * PI * 0.1e12
SIGMA_G   = 0.08 * 2 * PI * 1e12
F_U_BI    = 0.6
F_U       = 1.0

# Regime boundaries (THz)
GAMMA_NARROW_MAX  = 0.07
GAMMA_OPTIMAL_MAX = 0.15
# > 0.15 THz = broad

# Reference systems
REFERENCE_SYSTEMS = {
    "M87":    {"style": "sharp",   "note": "M87 jet, high-contrast knots"},
    "SgrA":   {"style": "sharp",   "note": "Sgr A* intermittent jets"},
    "CenA":   {"style": "broad",   "note": "Centaurus A, sustained FR I/II"},
    "TXS0506":{"style": "optimal", "note": "TXS 0506+056, IceCube assoc."},
    "3C273":  {"style": "optimal", "note": "3C273 quasar jet"},
}


# ── §1  CORE MODULATION FUNCTIONS ─────────────────────────────────────────

def m_jet(Gamma_THz: float, A_jet: float = 1.5) -> float:
    """Jet modulation factor M_jet(Γ)."""
    Gamma_rad = 2 * PI * Gamma_THz * 1e12
    delta = Gamma_rad - GAMMA_0
    phonon_coupling = math.exp(-delta**2 / (2 * SIGMA_G**2))
    buoyancy = 2 * F_U_BI / max(F_U, 1e-50) - 1
    return 1 + A_jet * phonon_coupling * S26 * buoyancy


def quality_factor(Gamma_THz: float) -> float:
    """Q = ω_SCm / (2Γ) — resonance sharpness diagnostic."""
    Gamma_rad = 2 * PI * Gamma_THz * 1e12
    return OMEGA_SCM / (2 * max(Gamma_rad, 1e-30))


def classify_regime(Gamma_THz: float) -> str:
    """Classify linewidth into narrow/optimal/broad."""
    if Gamma_THz <= GAMMA_NARROW_MAX:
        return "narrow"
    elif Gamma_THz <= GAMMA_OPTIMAL_MAX:
        return "optimal"
    return "broad"


def collimation_descriptor(Gamma_THz: float) -> str:
    """Human-readable collimation description."""
    regime = classify_regime(Gamma_THz)
    if regime == "narrow":
        return "sharp, highly collimated knots"
    elif regime == "optimal":
        return "balanced, observationally matched"
    return "diffuse wind component"


# ── §2  SYSTEMATIC SWEEP CLASS ─────────────────────────────────────────────

class LinewidthJetModulationSweep:
    """Systematic Γ sweep engine for jet modulation mapping.

    Sweeps Γ from 0.01 to 1.0 THz in configurable steps, producing
    modulation factor, quality factor, regime, and collimation for each.
    """

    DEFAULT_GAMMAS = [0.01, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.50, 1.00]

    def compute(self, dataset: dict = None) -> dict:
        d = dataset or {}
        A_jet = float(d.get("A_jet", 1.5))
        gammas = d.get("gammas", self.DEFAULT_GAMMAS)

        rows = []
        for gTHz in gammas:
            Mj = m_jet(gTHz, A_jet)
            Q = quality_factor(gTHz)
            regime = classify_regime(gTHz)
            desc = collimation_descriptor(gTHz)
            rows.append({
                "Gamma_THz": gTHz,
                "M_jet": Mj,
                "Q": Q,
                "regime": regime,
                "collimation": desc,
                "enhancement_range": f"{Mj:.1f}–{Mj * 1.3:.1f}×",
            })

        return {
            "sweep": rows,
            "A_jet": A_jet,
            "primary_equations": [
                "M_jet(Γ) = exp(-(ω−ω_SCm)²/(2Γ²)) · S₂₆ · (2F_UBi/F_U − 1)",
                "Q = ω_SCm / (2Γ)",
            ] + [f"Γ={r['Gamma_THz']:.2f}: M_jet={r['M_jet']:.2f}, Q={r['Q']:.1f}, {r['regime']}" for r in rows],
            "note": "Session 213. Linewidth jet modulation sweep.",
        }


# ── §3  COLLIMATION-POWER MAPPING ─────────────────────────────────────────

class CollimationPowerMapping:
    """Maps linewidth to observable jet properties: opening angle, brightness
    contrast, and Doppler factor sensitivity."""

    def compute(self, dataset: dict = None) -> dict:
        d = dataset or {}
        A_jet = float(d.get("A_jet", 1.5))

        canonical_gammas = [0.05, 0.10, 0.30]
        results = []
        for gTHz in canonical_gammas:
            Mj = m_jet(gTHz, A_jet)
            Q = quality_factor(gTHz)
            # Opening half-angle inversely proportional to Q (simplified)
            theta_half_deg = max(0.5, 30.0 / max(Q, 0.1))
            # Brightness contrast proportional to M_jet
            contrast = Mj
            results.append({
                "Gamma_THz": gTHz,
                "M_jet": Mj,
                "Q": Q,
                "theta_half_deg": theta_half_deg,
                "brightness_contrast": contrast,
                "regime": classify_regime(gTHz),
            })

        return {
            "mapping": results,
            "primary_equations": [
                "θ_half ∝ 1/Q → narrow Γ = tight jet",
                "Contrast ∝ M_jet → sharp Γ = high flare contrast",
            ] + [f"Γ={r['Gamma_THz']}: θ={r['theta_half_deg']:.1f}°, contrast={r['brightness_contrast']:.2f}" for r in results],
        }


# ── §4  REFERENCE SYSTEM MATCHER ──────────────────────────────────────────

class ReferenceSystemMatcher:
    """Matches observed jet morphology to optimal Γ regime."""

    def compute(self, dataset: dict = None) -> dict:
        d = dataset or {}
        target_style = d.get("style", "optimal")

        matches = []
        for name, info in REFERENCE_SYSTEMS.items():
            if info["style"] == target_style:
                matches.append({"name": name, **info})

        gamma_rec = {
            "sharp": "Γ ≤ 0.07 THz (narrow)",
            "optimal": "0.07 < Γ ≤ 0.15 THz",
            "broad": "Γ > 0.15 THz",
        }.get(target_style, "unknown")

        return {
            "target_style": target_style,
            "recommended_gamma": gamma_rec,
            "matching_systems": matches,
        }


# ── MAIN ───────────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("LINEWIDTH JET MODULATION ENGINE (Session 213)")
    print("=" * 72)

    print("\n── §1 Systematic Γ Sweep ──")
    sweep = LinewidthJetModulationSweep()
    result = sweep.compute()
    print(f"{'Γ (THz)':>10} {'M_jet':>8} {'Q':>10} {'Regime':>10} {'Collimation'}")
    print("-" * 65)
    for r in result["sweep"]:
        print(f"{r['Gamma_THz']:>10.2f} {r['M_jet']:>8.2f} {r['Q']:>10.1f} "
              f"{r['regime']:>10} {r['collimation']}")

    print("\n── §2 Collimation-Power Mapping ──")
    cpm = CollimationPowerMapping()
    res2 = cpm.compute()
    for r in res2["mapping"]:
        print(f"  Γ = {r['Gamma_THz']:.2f}: θ_half = {r['theta_half_deg']:.1f}°, "
              f"contrast = {r['brightness_contrast']:.2f}, {r['regime']}")

    print("\n── §3 Reference System Match ──")
    for style in ["sharp", "optimal", "broad"]:
        matcher = ReferenceSystemMatcher()
        res3 = matcher.compute({"style": style})
        names = ", ".join(m["name"] for m in res3["matching_systems"])
        print(f"  {style}: {res3['recommended_gamma']} → [{names}]")

    print(f"\n{'=' * 72}")
    print("LINEWIDTH JET MODULATION COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
