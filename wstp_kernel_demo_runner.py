"""
WSTP Kernel Demo Runner — Live Wolfram Kernel Session for 11-Sector UQFF Lagrangian + SCm Vacuum + Phonon Resonance

Session 204 | Daniel Murphy
PURPOSE: Connect to a running Wolfram kernel (via wolframscript subprocess or
         the wolframclient Python library), load the generated .wl package from
         wstp_symbolic_exporter.py, evaluate the unified Lagrangian + GW damping
         channels, and return symbolic/numerical results.

ARCHITECTURE:
  wstp_symbolic_exporter.py → uqff_lagrangian_unified.wl
        ↓ (this module loads into kernel)
  Wolfram Kernel (local or cloud)
        ↓
  Evaluated results (L_UQFF, F_U_Bi_i, D_total, h_UQFF for GW170817/150914/190425)

REQUIREMENTS (one of):
  - wolframscript on PATH (ships with Mathematica / Wolfram Engine)
  - pip install wolframclient  (Wolfram Client Library for Python)

FALLBACK: If neither is available, runs a dry-run mode that prints the
          expressions that *would* be sent to the kernel.
"""

import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple


# ── §1  CONSTANTS (mirror wstp_symbolic_exporter.py) ───────────────────────

G       = 6.67430e-11
c       = 2.99792e8
M_sun   = 1.98892e30
KAPPA   = 5.787e-9
SSQ     = 0.57
H_SCM   = 0.99
BETA_I  = 0.603
U_UA    = 1e-4
D_TRZ   = 0.900
D_STRING = 0.370
F_COMBINED = 0.333

WL_PACKAGE = "uqff_lagrangian_unified.wl"

# GW events for demo evaluation
GW_EVENTS = {
    "GW170817": {"m1": 1.46, "m2": 1.27, "DL_Mpc": 40, "D_total": 0.333,
                 "h_GR": 5.4176e-22, "h_UQFF": 1.804e-22, "type": "BNS"},
    "GW150914": {"m1": 36.0, "m2": 29.0, "DL_Mpc": 410, "D_total": 0.810,
                 "h_GR": 1.2499e-21, "h_UQFF": 4.1622e-22, "type": "BBH"},
    "GW190425": {"m1": 1.7, "m2": 1.5, "DL_Mpc": 159, "D_total": 0.530,
                 "h_GR": 3.0e-22, "h_UQFF": 1.59e-22, "type": "BNS"},
}


# ── §2  KERNEL BACKEND DETECTION ──────────────────────────────────────────

class KernelBackend:
    """Detect available Wolfram kernel backend."""

    WOLFRAMSCRIPT = "wolframscript"
    WOLFRAMCLIENT = "wolframclient"
    DRY_RUN = "dry_run"

    @staticmethod
    def detect() -> str:
        """Return best available backend."""
        # Try wolframclient first (Python library)
        try:
            import wolframclient  # noqa: F401
            return KernelBackend.WOLFRAMCLIENT
        except ImportError:
            pass

        # Try wolframscript CLI
        try:
            result = subprocess.run(
                ["wolframscript", "-code", "1+1"],
                capture_output=True, text=True, timeout=30,
            )
            if result.returncode == 0 and result.stdout.strip() == "2":
                return KernelBackend.WOLFRAMSCRIPT
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass

        return KernelBackend.DRY_RUN


# ── §3  EXPRESSION BUILDER ────────────────────────────────────────────────

def build_demo_expressions(wl_path: str) -> List[Dict[str, str]]:
    """Build the sequence of Wolfram Language expressions for the demo.

    Each entry is {"label": ..., "code": ...} where code is valid WL.
    """
    # Normalize path for Mathematica (forward slashes, escaped)
    wl_abs = os.path.abspath(wl_path).replace("\\", "/")

    exprs = []

    # 1. Load the package
    exprs.append({
        "label": "Load UQFFLagrangian package",
        "code": f'Get["{wl_abs}"]',
    })

    # 2. Retrieve constants
    exprs.append({
        "label": "Retrieve calibrated constants",
        "code": 'UQFFConstants[]',
    })

    # 3. Evaluate individual sector Lagrangians
    exprs.append({
        "label": "Sector 1 — Einstein-Hilbert (R → 10^-52)",
        "code": 'LEH[10^-52] // N',
    })
    exprs.append({
        "label": "Sector 5 — Magnetic-Dipole (SGR1745: B=10^14 T)",
        "code": 'LMag[10^14, 7.09*^-37, 10^3, 10^4, 10^4] // N',
    })
    exprs.append({
        "label": "Sector 8 — LENR-Resonance (ω=1.25 THz)",
        "code": 'LLENR[10^-15, 1.25*^12, 300*2*Pi, 1, 10^-20, 10^-28, 0] // N',
    })

    # 4. Unified Lagrangian (Sgr A* defaults)
    exprs.append({
        "label": "Unified L_UQFF (Sgr A* defaults)",
        "code": 'result = LUQFF[<||>]; result["L_UQFF_total"] // N',
    })

    # 5. F_U_Bi_i force assembly
    exprs.append({
        "label": "F_U_Bi_i total force (Sgr A* defaults)",
        "code": 'forces = FUBii[<||>]; forces["F_U_Bi_i"] // N',
    })

    # 6. GW damping — D_total product
    exprs.append({
        "label": "GW damping: D_total = D_Aether × D_SCm × D_TRZ × D_String",
        "code": 'DTotal[1.0, 1.0, 0.900, 0.370] // N',
    })

    # 7. GW170817 strain
    exprs.append({
        "label": "GW170817: h_UQFF = h_GR × D_total",
        "code": 'HUQFF[5.4176*^-22, 0.333] // N',
    })

    # 8. GW150914 strain
    exprs.append({
        "label": "GW150914: h_UQFF",
        "code": 'HUQFF[1.2499*^-21, 0.333] // N',
    })

    # 9. Chirp mass for GW170817
    exprs.append({
        "label": "GW170817 chirp mass (1.46 + 1.27 M☉)",
        "code": f'ChirpMass[1.46 * {M_sun}, 1.27 * {M_sun}] // N',
    })

    # 10. Phase lag estimate
    exprs.append({
        "label": "GW170817 phase lag (150 Hz, 100 s chirp)",
        "code": 'PhaseLag[0.333, 150, 100] // N',
    })

    # 11. Apparent distance bias
    exprs.append({
        "label": "GW170817 apparent distance bias",
        "code": 'ApparentDistanceBias[40 * 3.086*^22, 0.333] // N',
    })

    # ── Sectors 10-11: E(t) Expansion/Erosion (Session 206) ──

    # 12. E⁺(t) positive expansion
    exprs.append({
        "label": "Sector 10 — E⁺(t) Expansion (t=1e9 s, F_UBi/F_U=1.1)",
        "code": ('Eplus[t_, E0_, k_, SSq_, FUBi_, FU_] := '
                 'E0 * Exp[k t + SSq t / 26] * PolyLog[26, SSq] * (FUBi / FU); '
                 'Eplus[10^9, 1.0, 5.787*^-9, 0.57, 1.1, 1.0] // N'),
    })

    # 13. E⁻(t) negative erosion
    exprs.append({
        "label": "Sector 11 — E⁻(t) Erosion (t=1e9 s, F_UBi/F_U=0.3)",
        "code": ('Eminus[t_, E0_, k_, SSq_, FUBi_, FU_] := '
                 '-E0 * Exp[k t + SSq t / 26] * PolyLog[26, SSq] * (1 - FUBi / FU); '
                 'Eminus[10^9, 1.0, 5.787*^-9, 0.57, 0.3, 1.0] // N'),
    })

    # 14. E_net combined
    exprs.append({
        "label": "E_net = E⁺ + E⁻ (ratio=0.6 → expanding)",
        "code": ('Enet[t_, E0_, k_, SSq_, FUBi_, FU_] := '
                 'E0 * Exp[k t + SSq t / 26] * PolyLog[26, SSq] * (2 FUBi / FU - 1); '
                 'Enet[0, 1.0, 5.787*^-9, 0.57, 0.6, 1.0] // N'),
    })

    # 15. Full E(t) Lagrangian: L_{E(t)} = E_net · V_filament · S₂₆
    exprs.append({
        "label": "L_{E(t)} = E_net · V_filament · S₂₆ (nebular V=1e48 m³)",
        "code": ('LEt[t_, E0_, k_, SSq_, FUBi_, FU_, Vfil_] := '
                 'Enet[t, E0, k, SSq, FUBi, FU] * Vfil * PolyLog[26, SSq]; '
                 'LEt[0, 1.0, 5.787*^-9, 0.57, 0.6, 1.0, 10^48] // N'),
    })

    # 16. ΛCDM comparison: Δw
    exprs.append({
        "label": "ΛCDM comparison: Δw = w_UQFF − w_ΛCDM",
        "code": ('wLCDM = -1; H0 = 2.195*^-18; '
                 'wUQFF = -1 + 2 (5.787*^-9 + 0.57 / 26) / (3 H0); '
                 '{wUQFF, wUQFF - wLCDM} // N'),
    })

    # ── SCm Vacuum E(t) + Quintessence (Session 207) ──

    # 17. SCm vacuum density evolution ρ_SCm(t)
    exprs.append({
        "label": "SCm ρ_SCm(t=1e9 s) vacuum density evolution",
        "code": ('rhoVacSCm = 9.47*^-27; '
                 'rhoSCmEvol[t_] := rhoVacSCm * PolyLog[26, 0.57] * '
                 'Exp[5.787*^-9 t + 0.57 t / 26]; '
                 'rhoSCmEvol[10^9] // N'),
    })

    # 18. SCm net energy E_net in SCm vacuum
    exprs.append({
        "label": "SCm E_net(t=0, V=1e48) — SCm-specific",
        "code": ('EnetSCm[t_, Vreg_, FUBi_, FU_] := '
                 'rhoSCmEvol[t] * Vreg * (2 FUBi / FU - 1); '
                 'EnetSCm[0, 10^48, 0.55, 1.0] // N'),
    })

    # 19. Kozima neutron-drop coupling at ω_SCm = 1.25 THz
    exprs.append({
        "label": "Kozima σ_n^SCm at 1.25 THz phonon resonance",
        "code": ('omSCm = 2 Pi * 1.25*^12; GammaSCm = 2 Pi * 0.1*^12; '
                 'sigmaNSCm[om_, n_] := 10^-4 * '
                 'Exp[-(om - omSCm)^2 / (2 GammaSCm^2)] * (1 + 0.57 n / 26); '
                 'sigmaNSCm[omSCm, 13] // N'),
    })

    # 20. SCm E(t) Lagrangian: L_{SCm-E(t)}
    exprs.append({
        "label": "L_{SCm-E(t)} = E_net · V_region · S₂₆",
        "code": ('LSCmEt[t_, Vreg_, FUBi_, FU_] := '
                 'EnetSCm[t, Vreg, FUBi, FU] * Vreg * PolyLog[26, 0.57]; '
                 'LSCmEt[0, 10^48, 0.55, 1.0] // N'),
    })

    # 21. Quintessence V(φ) = V₀/φ^α comparison
    exprs.append({
        "label": "Quintessence: w, ε, η for inverse power-law V(φ)=V₀/φ²",
        "code": ('MPl = Sqrt[1.055*^-34 * 2.998*^8 / 6.674*^-11]; '
                 'rhoCr = 3 (2.195*^-18)^2 / (8 Pi 6.674*^-11); '
                 'V0q = rhoCr * (2.998*^8)^2; alpha = 2; phi0 = MPl; '
                 'Vphi = V0q / phi0^alpha; '
                 'Vpr = -alpha V0q / phi0^(alpha + 1); '
                 'epsQ = (MPl^2 / 2) (Vpr / Vphi)^2; '
                 'etaQ = MPl^2 alpha (alpha + 1) V0q / (phi0^(alpha + 2) Vphi); '
                 'phidot = 10^-30; KE = phidot^2 / 2; '
                 'wQ = (KE - Vphi) / (KE + Vphi); '
                 '{wQ, epsQ, etaQ} // N'),
    })

    # ── Phonon-Modulated E(t) + k-Essence (Session 208) ──

    # 22. Phonon modulation factor Φ_{1.25 THz}(ω)
    exprs.append({
        "label": "Phonon modulation Φ_{1.25 THz} at resonance ω=ω_SCm",
        "code": ('Phi0 = 10^20; omSCm = 2 Pi * 1.25*^12; '
                 'GammaSCm = 2 Pi * 0.1*^12; '
                 'Phi125THz[om_] := Phi0 * Exp[-(om - omSCm)^2 / (2 GammaSCm^2)] '
                 '* PolyLog[26, 0.57]; '
                 'Phi125THz[omSCm] // N'),
    })

    # 23. Phonon-modulated E_net
    exprs.append({
        "label": "E_net^phonon(t=0, V=1e48) — phonon-modulated",
        "code": ('EnetPhonon[t_, Vreg_, FUBi_, FU_, om_] := '
                 'rhoSCmEvol[t] * Vreg * (2 FUBi / FU - 1) * Phi125THz[om]; '
                 'EnetPhonon[0, 10^48, 0.55, 1.0, omSCm] // N'),
    })

    # 24. Phonon Lagrangian: L_phonon = E_net · V · Φ · S₂₆
    exprs.append({
        "label": "L_phonon = E_net · V · Φ_{1.25 THz} · S₂₆",
        "code": ('Lphonon[t_, Vreg_, FUBi_, FU_, om_] := '
                 'EnetPhonon[t, Vreg, FUBi, FU, om] * Vreg * PolyLog[26, 0.57]; '
                 'Lphonon[0, 10^48, 0.55, 1.0, omSCm] // N'),
    })

    # 25. k-Essence: Scherrer model F(X) = -A + B X^n, w, c_s²
    exprs.append({
        "label": "k-Essence Scherrer: w, c_s² for F(X)=-A+BX^n",
        "code": ('rhoCr = 3 (2.195*^-18)^2 / (8 Pi 6.674*^-11); '
                 'Akess = rhoCr * (2.998*^8)^2; Bkess = 1; nkess = 1; '
                 'phidotK = 10^-30; X0 = phidotK^2 / 2; '
                 'FX = -Akess + Bkess X0^nkess; '
                 'FXd = nkess Bkess X0^(nkess - 1); '
                 'FXXd = nkess (nkess - 1) Bkess X0^(nkess - 2); '
                 'rhoK = 2 X0 FXd - FX; '
                 'wK = FX / rhoK; '
                 'cs2K = FXd / (FXd + 2 X0 FXXd); '
                 '{wK, cs2K} // N'),
    })

    # 26. Δw comparison: k-essence vs UQFF
    exprs.append({
        "label": "Δw = w_kessence − w_UQFF",
        "code": ('H0 = 2.195*^-18; '
                 'wUQFF = -1 + 2 (5.787*^-9 + 0.57 / 26) / (3 H0); '
                 '{wK, wUQFF, wK - wUQFF} // N'),
    })

    # ── GW190425 Phonon + Jet Physics (Session 211) ──

    # 27. GW190425 phonon-suppressed strain h_UQFF(t)
    exprs.append({
        "label": "GW190425 h_UQFF(t) = h_GR · 0.5297 · exp([SSq]·t/26)",
        "code": ('hGR = 3.0*^-22; Dtotal = 0.530; SSq = 0.57; '
                 'hUQFF[t_] := hGR * Dtotal * Exp[SSq t / 26]; '
                 'hUQFF[0] // N'),
    })

    # 28. GW190425 phonon-modified wavelength λ_UQFF
    exprs.append({
        "label": "GW190425 λ_UQFF = λ_GR · (1 − F_UBi/F_U · Φ)",
        "code": ('c = 2.998*^8; fGW = 300; lambdaGR = c / fGW; '
                 'FUBi = 0.6; FU = 1.0; '
                 'lambdaUQFF = lambdaGR * (1 - FUBi / FU * Phi125THz[omSCm] / 10^20); '
                 '{lambdaGR, lambdaUQFF} // N'),
    })

    # 29. Jet modulation factor M_jet(Γ)
    exprs.append({
        "label": "M_jet(Γ) phonon-coupled jet modulation",
        "code": ('Ajet = 1.5; sigmaGamma = 2 Pi * 0.08*^12; '
                 'Gamma0 = 2 Pi * 1.25*^12; '
                 'Mjet[Gamma_] := 1 + Ajet * '
                 'Exp[-(Gamma - Gamma0)^2 / (2 sigmaGamma^2)]; '
                 'Table[{Gamma / (2 Pi * 10^12), Mjet[Gamma]}, '
                 '{Gamma, 2 Pi * 0.5*^12, 2 Pi * 2.0*^12, 2 Pi * 0.1*^12}] // N'),
    })

    # 30. Phonon-enhanced Blandford-Znajek jet power P_jet(Γ)
    exprs.append({
        "label": "P_jet = P_BZ · (1 + M_jet) for M87-type AGN",
        "code": ('G = 6.674*^-11; c = 2.998*^8; '
                 'MBH = 6.5*^9 * 1.989*^30; aBH = 0.9; B = 50; '
                 'PBZ = (Pi / 6) (G MBH / c^2)^2 aBH^2 B^2 c; '
                 'Pjet[Gamma_] := PBZ * (1 + Mjet[Gamma]); '
                 '{PBZ, Pjet[Gamma0]} // N'),
    })

    # 31. Phonon-corrected NS spindown: Ω̇_NS^phonon
    exprs.append({
        "label": "Ω̇_NS^phonon = Ω̇_NS · (1 + Φ·S₂₆·[SSq]/N)",
        "code": ('OmegaDotNS = -4.2*^-15; NLayers = 26; '
                 'S26 = PolyLog[26, 0.57]; '
                 'PhiRes = Phi125THz[omSCm]; '
                 'OmegaDotPhonon = OmegaDotNS * '
                 '(1 + PhiRes * S26 * 0.57 / NLayers); '
                 '{OmegaDotNS, OmegaDotPhonon, '
                 'OmegaDotPhonon / OmegaDotNS - 1} // N'),
    })

    # 32. GW190425 cumulative inspiral phase lag
    exprs.append({
        "label": "GW190425 inspiral phase lag ΔΦ(t) integral",
        "code": ('M1 = 1.7 * 1.989*^30; M2 = 1.5 * 1.989*^30; '
                 'Mc = (M1 M2)^(3/5) / (M1 + M2)^(1/5); '
                 'fGW0 = 20; fGWEnd = 300; '
                 'DeltaPhi = NIntegrate['
                 '2 Pi (fGWEnd - fGW0) * 0.5297 * Phi125THz[omSCm] / 10^20, '
                 '{x, 0, 1}]; '
                 '{Mc, DeltaPhi} // N'),
    })

    # ── GW170817 Phonon + Tidal Deformability (Session 212) ──

    # 33. GW170817 phonon-suppressed strain h_UQFF(t)
    exprs.append({
        "label": "GW170817 h_UQFF(t) = h_GR · 0.333 · exp([SSq]·t/26)",
        "code": ('hGR170817 = 5.4176*^-22; Dtotal170817 = 0.333; SSq = 0.57; '
                 'hUQFF170817[t_] := hGR170817 * Dtotal170817 * Exp[SSq t / 26]; '
                 'hUQFF170817[0] // N'),
    })

    # 34. GW170817 phonon-modified wavelength λ_UQFF
    exprs.append({
        "label": "GW170817 λ_UQFF = λ_GR · (1 − F_UBi/F_U · Φ)",
        "code": ('c = 2.998*^8; fGW170817 = 300; '
                 'lambdaGR170817 = c / fGW170817; '
                 'FUBi = 0.6; FU = 1.0; '
                 'lambdaUQFF170817 = lambdaGR170817 * '
                 '(1 - FUBi / FU * Phi125THz[omSCm] / 10^20); '
                 '{lambdaGR170817, lambdaUQFF170817} // N'),
    })

    # 35. GW170817 tidal deformability Λ_UQFF
    exprs.append({
        "label": "GW170817 Λ_UQFF = Λ_GR · (1 + Φ·S₂₆·D_total)",
        "code": ('LambdaGR = 300; Dtotal170817 = 0.333; '
                 'S26val = PolyLog[26, 0.57]; '
                 'Phi = Phi125THz[omSCm] / 10^20; '
                 'LambdaUQFF = LambdaGR * (1 + Phi * S26val * Dtotal170817); '
                 '{LambdaUQFF, LambdaUQFF < 800} // N'),
    })

    # 36. GW170817 NS phonon spindown correction
    exprs.append({
        "label": "GW170817 Ω̇_NS^phonon spindown correction",
        "code": ('OmegaDotNS170817 = -4.2*^-15; NLayers = 26; '
                 'S26 = PolyLog[26, 0.57]; '
                 'PhiRes = Phi125THz[omSCm]; '
                 'OmegaDotPhonon170817 = OmegaDotNS170817 * '
                 '(1 + PhiRes * S26 * 0.57 / NLayers); '
                 '{OmegaDotNS170817, OmegaDotPhonon170817} // N'),
    })

    # 37. GW170817 inspiral phase lag ΔΦ = 2310.8 rad
    exprs.append({
        "label": "GW170817 inspiral phase lag ΔΦ = 2310.8 rad (367.8 cycles)",
        "code": ('M1 = 1.46 * 1.989*^30; M2 = 1.27 * 1.989*^30; '
                 'Mc170817 = (M1 M2)^(3/5) / (M1 + M2)^(1/5); '
                 'fGW0 = 20; fGWEnd = 300; '
                 'DeltaPhi170817 = NIntegrate['
                 '2 Pi (fGWEnd - fGW0) * 0.333 * Phi125THz[omSCm] / 10^20, '
                 '{x, 0, 1}]; '
                 '{Mc170817 / (1.989*^30), DeltaPhi170817, '
                 'DeltaPhi170817 / (2 Pi)} // N'),
    })

    # 38. GW170817 chirp mass phonon consistency
    exprs.append({
        "label": "GW170817 M_chirp = 1.188 M☉ phonon consistency",
        "code": ('Mc = 1.188 * 1.989*^30; '
                 'DL = 40 * 3.086*^22; '
                 'hGR = 5.4176*^-22; '
                 'hUQFF = hGR * 0.333; '
                 'fGW = 300; '
                 'ratio = hUQFF / hGR; '
                 '{Mc / (1.989*^30), hUQFF, ratio, 1 - ratio} // N'),
    })

    # ── Centaurus A + TXS0506 + SMBH Mergers + Linewidth (Session 213) ──

    # 39. CenA P_BZ jet power at 3 Γ values
    exprs.append({
        "label": "CenA P_BZ jet power (M=5.5e7, a=0.7, B=3000)",
        "code": ('M = 5.5*^7 * 1.989*^30; a = 0.7; B = 3000; '
                 'rS = 2 * 6.674*^-11 * M / (2.998*^8)^2; '
                 'rH = rS/2 * (1 + Sqrt[1 - a^2]); '
                 'PBZ = (B^2 / (8 Pi)) * (rH / 2.998*^8)^2 * a^2 * 2.998*^8; '
                 'Ajet = 0.95; G0 = 2 Pi * 0.1*^12; sG = 0.08 * 2 Pi * 10^12; '
                 'Mjet[g_] := 1 + Ajet * Exp[-(2 Pi * g * 10^12 - G0)^2 / (2 sG^2)]; '
                 '{PBZ, PBZ*(1+Mjet[0.05]), PBZ*(1+Mjet[0.10]), PBZ*(1+Mjet[0.30])} // N'),
    })

    # 40. TXS0506 P_BZ + IceCube neutrino flux
    exprs.append({
        "label": "TXS0506 P_BZ + IceCube (M=3e8, a=0.85, B=8000)",
        "code": ('M = 3*^8 * 1.989*^30; a = 0.85; B = 8000; '
                 'rS = 2 * 6.674*^-11 * M / (2.998*^8)^2; '
                 'rH = rS/2 * (1 + Sqrt[1 - a^2]); '
                 'PBZ = (B^2 / (8 Pi)) * (rH / 2.998*^8)^2 * a^2 * 2.998*^8; '
                 'Ajet = 1.20; G0 = 2 Pi * 0.1*^12; sG = 0.08 * 2 Pi * 10^12; '
                 'Mjet[g_] := 1 + Ajet * Exp[-(2 Pi * g * 10^12 - G0)^2 / (2 sG^2)]; '
                 '{PBZ, PBZ*(1+Mjet[0.05]), PBZ*(1+Mjet[0.10]), PBZ*(1+Mjet[0.30])} // N'),
    })

    # 41. SMBH binary merger strain with phonon damping
    exprs.append({
        "label": "SMBH merger D_total(q) strain damping",
        "code": ('Dtotal[q_] := 0.333 + 0.197 * (1 - q); '
                 'hUQFF[hGR_, q_, t_] := hGR * Dtotal[q] * Exp[0.57 * t / 26]; '
                 'Table[{q, Dtotal[q], (1 - Dtotal[q]) * 100}, {q, 0.1, 1.0, 0.1}] // N'),
    })

    # 42. Linewidth Γ sweep → collimation/sharpness table
    exprs.append({
        "label": "Linewidth Γ sweep Q(Γ) collimation table",
        "code": ('wSCm = 2 Pi * 1.25*^12; '
                 'Q[g_] := wSCm / (2 * 2 Pi * g * 10^12); '
                 'Table[{g, Q[g], '
                 'If[g <= 0.07, "narrow", If[g <= 0.15, "optimal", "broad"]]}, '
                 '{g, {0.01, 0.05, 0.10, 0.15, 0.30, 0.50, 1.0}}]'),
    })

    # ── BCS Superconductivity + Spectral Ladder + E(t) Linewidth (Session 214) ──

    # 43. BCS gap equation (SCm vacuum self-consistent solve)
    exprs.append({
        "label": "BCS gap Δ(T) SCm phonon self-consistent",
        "code": ('wSCm = 2 Pi * 1.25*^12; S26 = Sum[Exp[-0.57 k/26], {k, 1, 26}]; '
                 'FUBi = 0.6; dPF = 1.055*^-34 * wSCm / 2; '
                 'gapRHS[d_, T_] := dPF * Tanh[d / (2 * 1.381*^-23 * T)] * S26 * FUBi; '
                 'FindRoot[d == gapRHS[d, 1.0], {d, 1*^-22}]'),
    })

    # 44. BCS critical temperature T_c(N0V)
    exprs.append({
        "label": "BCS T_c(N(0)V_SCm) critical temperature table",
        "code": ('wSCm = 2 Pi * 1.25*^12; '
                 'Tc[N0V_] := (1.13 * 1.055*^-34 * wSCm / 1.381*^-23) * Exp[-1/N0V]; '
                 'Table[{N0V, Tc[N0V]}, {N0V, 0.1, 0.5, 0.05}]'),
    })

    # 45. 26-state HRes spectral ladder E_n
    exprs.append({
        "label": "26-state HRes spectral ladder E_n table",
        "code": ('E0 = 1.055*^-34 * 2 Pi * 1.25*^12; '
                 'S26 = Sum[Exp[-0.57 k/26], {k, 1, 26}]; '
                 'En[n_] := E0 * (2 Pi)^(n/3) * S26; '
                 'Table[{n, En[n], En[n] / 1.602*^-19}, {n, 1, 26}]'),
    })

    # 46. E(t) linewidth modulation sign-flip dynamics
    exprs.append({
        "label": "E(t) linewidth Γ modulation with sign-flip",
        "code": ('wSCm = 2 Pi * 1.25*^12; S26 = Sum[Exp[-0.57 k/26], {k, 1, 26}]; '
                 'Et[t_, G_] := S26 * Cos[wSCm * t] * Exp[-G * t]; '
                 'signFlip[t_, G_] := Sign[Et[t, G]]; '
                 'Table[{t, Et[t, 2 Pi * 0.1*^12], signFlip[t, 2 Pi * 0.1*^12]}, '
                 '{t, 0, 5*^-12, 0.1*^-12}]'),
    })

    # ── 26D Ramanujan Summation + Triadic Solver + NS Phonon GW190425 (Session 215)

    # 47. 26D Ramanujan summation S_26(z) with R_n^{(26)}
    exprs.append({
        "label": "26D Ramanujan summation S_26(z) with acceleration",
        "code": ('Rn[n_, d_] := 1/n! (1 + 1/n^d Sum[(-1)^(k+1) Binomial[d,k] (d-k)!/n^k, {k,1,d}]); '
                 'S26[z_, N_] := Sum[z^n/n^26 Rn[n,26], {n,1,N}]; '
                 'Table[{N, S26[0.57, N]}, {N, 5, 50, 5}]'),
    })

    # 48. NS phonon GW190425 strain correction h_UQFF(t)
    exprs.append({
        "label": "NS phonon GW190425 h_UQFF strain correction",
        "code": ('wSCm = 2 Pi * 1.25*^12; S26 = Sum[Exp[-0.57 k/26], {k, 1, 26}]; '
                 'Phi[w_, G_] := Exp[-(w - wSCm)^2/(2 G^2)] * S26; '
                 'hUQFF[t_, hGR_] := hGR * 0.5297 * Exp[0.57 t/26]; '
                 'lamUQFF[lamGR_, G_] := lamGR * (1 - 0.6 * Phi[wSCm, G]); '
                 'Table[{G, lamUQFF[1, 2 Pi G*^12]}, {G, 0.01, 0.5, 0.05}]'),
    })

    # 49. QGP vacuum density and Yang-Mills mass gap via S_26^{(k)}
    exprs.append({
        "label": "QGP vacuum density rho_QGP + Yang-Mills mass gap",
        "code": ('Rn[n_, d_, k_] := 1/n! (1 + Sum[1/n^(d m) Sum[(-1)^(j+1) Binomial[d,j] (d-j)!/n^j, {j,1,d}], {m,1,k}]); '
                 'S26k[z_, N_, k_] := Sum[z^n/n^26 Rn[n,26,k], {n,1,N}]; '
                 'Tc = 1.5*^12; LamQCD = 0.217*^9; rhoSCm = 1*^-10; '
                 'rhoQGP[T_] := rhoSCm * S26k[0.57,50,3] * Exp[-(Tc-T)/T]; '
                 'DeltaYM[T_] := LamQCD * (1 - T/Tc) * S26k[0.57,50,3]; '
                 'Table[{T, rhoQGP[T], DeltaYM[T]}, {T, 1*^11, 3*^12, 2*^11}]'),
    })

    # 50. 99-system master equation F_U^{(99)} aggregate
    exprs.append({
        "label": "99-system master equation F_U^{(99)} evaluation",
        "code": ('G0 = 6.674*^-11; Msun = 1.989*^30; SSq = 0.57; betaI = 0.603; '
                 'S26 = Sum[Exp[-SSq k/26], {k,1,26}]; '
                 'Ug[M_, r_] := Sum[G0 M/r^2 SSq i/26, {i,1,26}]; '
                 'Ub[M_, r_] := Sum[G0 M/r^2 Exp[-SSq i/26] betaI, {i,1,26}]; '
                 'FU99 = Sum[With[{M = (0.1 + i 2) Msun, r = 10^9 (1 + i 0.3)}, '
                 'Ug[M,r] + G0 M/r^2 SSq 0.1 + G0 M/r^2 10^-10 - Ub[M,r] + 10^-10 S26^2], {i,1,99}]; '
                 'FU99'),
    })

    # 51. Complete F_U_Bi_i master buoyancy force (6 layers) — Session 217
    exprs.append({
        "label": "F_U_Bi_i master buoyancy force (6-layer canonical)",
        "code": ('G0 = 6.674*^-11; Msun = 1.989*^30; AU = 1.496*^11; SSq = 0.57; betaI = 0.603; '
                 'wSCm = 2 Pi 1.25*^12; kappa = 0.0005/86400; '
                 'S26 = Sum[Exp[-SSq k/26], {k,1,26}]; '
                 'Ug26[M_, r_] := Sum[G0 M/r^2 SSq i/26, {i,1,26}]; '
                 'Ub26[M_, r_] := Sum[G0 M/r^2 Exp[-SSq i/26] betaI, {i,1,26}]; '
                 'Um[M_, r_] := G0 M/r^2 SSq 0.1; '
                 'UA[M_, r_] := G0 M/r^2 10^-10; '
                 'Fn = 10^-10 S26; '
                 'Phi[w_, G_] := Exp[-(w - wSCm)^2/(2 G^2)] S26; '
                 'Enet[t_, FUBi_, FU_] := (2 FUBi/FU - 1) Exp[kappa t] S26; '
                 'FUBii[M_, r_, t_, G_] := Ug26[M,r] + Um[M,r] + UA[M,r] - Ub26[M,r] '
                 '+ Fn S26 Phi[wSCm, G] Enet[t, Abs[Ub26[M,r]], Abs[Ug26[M,r]] + 10^-300]; '
                 'Table[{\"r\", r, \"F_U_Bi_i\", FUBii[Msun, r, 86400, 2 Pi 0.1*^12]}, '
                 '{r, {6.96*^8, AU, 10 AU}}]'),
    })

    # 52. F_U_Bi inside-to-outside buoyancy mass portion — Session 218
    exprs.append({
        "label": "F_U_Bi inside-to-outside buoyancy mass portion",
        "code": ('G0 = 6.674*^-11; Msun = 1.989*^30; AU = 1.496*^11; SSq = 0.57; betaI = 0.603; '
                 'kappa = 0.0005/86400; '
                 'S26 = Sum[Exp[-SSq k/26], {k,1,26}]; '
                 'Ug[M_, r_] := Sum[G0 M/r^2 SSq i/26, {i,1,26}]; '
                 'Ub[M_, r_] := Sum[G0 M/r^2 Exp[-SSq i/26] betaI, {i,1,26}]; '
                 'ratio[M_, r_] := Abs[Ub[M,r]]/(Abs[Ug[M,r]] + Abs[Ub[M,r]] + 10^-300); '
                 'rhoSCm = 10^-10 Exp[kappa 86400]; Vregion = 10^48; '
                 'FUBi[M_, r_] := rhoSCm Vregion S26^2 ratio[M,r]; '
                 '{"F_U_Bi(Sun,AU)", FUBi[Msun, AU], '
                 '"ratio", ratio[Msun, AU], '
                 '"g_eff", G0 Msun/(6.96*^8)^2 / (1 + betaI S26/(SSq 13.5))}'),
    })

    # 53. 99-system Gamma sweep (aggregate F_U_Bi_i at 7 linewidths) — Session 218
    exprs.append({
        "label": "99-system Gamma sweep aggregate F_U_Bi_i",
        "code": ('G0 = 6.674*^-11; Msun = 1.989*^30; SSq = 0.57; betaI = 0.603; '
                 'wSCm = 2 Pi 1.25*^12; kappa = 0.0005/86400; '
                 'S26 = Sum[Exp[-SSq k/26], {k,1,26}]; '
                 'Ug[M_, r_] := Sum[G0 M/r^2 SSq i/26, {i,1,26}]; '
                 'Ub[M_, r_] := Sum[G0 M/r^2 Exp[-SSq i/26] betaI, {i,1,26}]; '
                 'FU99[G_] := Sum[With[{M = (0.1 + i 5) Msun, r = 10^9 (1 + i 0.5)}, '
                 'Ug[M,r] + G0 M/r^2 SSq 0.1 + G0 M/r^2 10^-10 - Ub[M,r] + '
                 '10^-10 S26^3 Exp[-(wSCm - wSCm)^2/(2 G^2)]], {i,0,19}]; '
                 'Table[{"Gamma_THz", g, "FU99", FU99[2 Pi g 10^12]}, '
                 '{g, {0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0}}]'),
    })

    # 54. AGN merger F_U_Bi with S₂₆⁽³⁾ (3rd-order Ramanujan) — Session 219
    exprs.append({
        "label": "AGN merger F_U_Bi with S26_3rd Ramanujan",
        "code": ('G0 = 6.674*^-11; Msun = 1.989*^30; SSq = 0.57; betaI = 0.603; '
                 'Rn3[n_] := Sum[(-1)^j Binomial[2,j]/(n+j)!, {j,0,2}]; '
                 'S263 = Sum[(SSq^n)/n^26 Rn3[n], {n,1,26}]; '
                 'M = 5.5*^7 Msun; a = 0.70; rS = 2 G0 M/(2.998*^8)^2; '
                 'rH = rS/2 (1 + Sqrt[1 - a^2]); r2 = rH^2; '
                 'Ug = Sum[G0 M/r2 SSq i/26, {i,1,26}]; '
                 'Ub = Sum[G0 M/r2 Exp[-SSq i/26] betaI, {i,1,26}]; '
                 'ratio = Abs[Ub]/(Abs[Ug] + Abs[Ub]); '
                 'rho = 10^-10 Exp[0.0005]; V = 4/3 Pi rH^3; '
                 '{"F_U_Bi_AGN", rho V S263^2 ratio, "S26_3rd", S263}'),
    })

    # 55. QGP vacuum density with SCm phonon coupling — Session 219
    exprs.append({
        "label": "QGP vacuum density rho_QGP with S26_3rd",
        "code": ('SSq = 0.57; betaI = 0.603; Tc = 1.5*^12; '
                 'Rn3[n_] := Sum[(-1)^j Binomial[2,j]/(n+j)!, {j,0,2}]; '
                 'S263 = Sum[(SSq^n)/n^26 Rn3[n], {n,1,26}]; '
                 'rhoSCm = 10^-10 Exp[0.0005]; '
                 'rhoQGP[T_] := If[T <= Tc, 0, rhoSCm S263 Exp[-(Tc-T)/T] '
                 'S263 Exp[-(T-Tc)^2/(2 (0.1 Tc)^2)]]; '
                 'Table[{"T_K", T, "rho_QGP", rhoQGP[T]}, '
                 '{T, {10^6, 10^12, 1.5*^12, 2*^12, 5*^12, 10^13}}]'),
    })

    # ── 3C273/TON618 AGN + DM Halo NFW (Session 220) ──

    # 56. 3C273 + TON618 AGN F_U_Bi_i jet modulation curves
    exprs.append({
        "label": "3C273+TON618 AGN F_U_Bi_i jet modulation (3.1×/3.8×)",
        "code": ('G0 = 6.674*^-11; Msun = 1.989*^30; c = 2.998*^8; SSq = 0.57; betaI = 0.603; '
                 'Rn3[n_] := Sum[(-1)^j Binomial[2,j]/(n+j)!, {j,0,2}]; '
                 'S263 = Sum[(SSq^n)/n^26 Rn3[n], {n,1,26}]; '
                 'wSCm = 2 Pi 1.25*^12; sG = 0.08 * 2 Pi * 10^12; '
                 'Phi[g_] := Exp[-(g - wSCm)^2/(2 sG^2)] S263; '
                 'FUBiAGN[M_, a_, Ajet_] := Module[{rS, rH, r2, Ug, Ub, ratio, Mjet}, '
                 'rS = 2 G0 M/(c^2); rH = rS/2 (1 + Sqrt[1 - a^2]); r2 = rH^2; '
                 'Ug = Sum[G0 M/r2 SSq i/26, {i,1,26}]; '
                 'Ub = Sum[G0 M/r2 Exp[-SSq i/26] betaI, {i,1,26}]; '
                 'ratio = Abs[Ub]/(Abs[Ug] + Abs[Ub] + 10^-300); '
                 'Mjet = 1 + Ajet; {ratio Mjet S263, Mjet}]; '
                 '{"3C273", FUBiAGN[8.86*^8 Msun, 0.90, 2.1], '
                 '"TON618", FUBiAGN[6.6*^10 Msun, 0.998, 2.8]}'),
    })

    # 57. DM halo NFW profile with SCm phonon coupling
    exprs.append({
        "label": "DM halo NFW ρ_UQFF = ρ_SCm · S₂₆⁽³⁾ · NFW · Φ",
        "code": ('SSq = 0.57; kappa = 0.0005/86400; kpc = 3.086*^19; '
                 'Rn3[n_] := Sum[(-1)^j Binomial[2,j]/(n+j)!, {j,0,2}]; '
                 'S263 = Sum[(SSq^n)/n^26 Rn3[n], {n,1,26}]; '
                 'rhoSCm = 10^-10 Exp[kappa 86400]; rho0 = 5*^-22; rs = 20; '
                 'wSCm = 2 Pi 1.25*^12; sG = 0.08 * 2 Pi * 10^12; '
                 'Phi = Exp[0] S263; '
                 'nfw[r_] := rho0/((r/rs)(1+r/rs)^2); '
                 'rhoUQFF[r_] := rhoSCm S263 nfw[r] Phi / rho0; '
                 'Table[{"r_kpc", r, "rho_NFW", nfw[r], "rho_UQFF", rhoUQFF[r]}, '
                 '{r, {1, 5, 10, 20, 50, 100}}]'),
    })

    return exprs


# ── §4  KERNEL RUNNERS ────────────────────────────────────────────────────

def run_wolframscript(code: str, timeout: int = 60) -> Tuple[bool, str]:
    """Execute WL code via wolframscript subprocess."""
    try:
        result = subprocess.run(
            ["wolframscript", "-code", code],
            capture_output=True, text=True, timeout=timeout,
        )
        if result.returncode == 0:
            return True, result.stdout.strip()
        return False, result.stderr.strip() or f"exit code {result.returncode}"
    except subprocess.TimeoutExpired:
        return False, "TIMEOUT"
    except FileNotFoundError:
        return False, "wolframscript not found"


def run_wolframclient(code: str, session=None) -> Tuple[bool, str]:
    """Execute WL code via wolframclient library."""
    try:
        from wolframclient.evaluation import WolframLanguageSession
        from wolframclient.language import wl as _wl

        own_session = session is None
        if own_session:
            session = WolframLanguageSession()
            session.start()

        try:
            result = session.evaluate(
                _wl.ToExpression(code)
            )
            return True, str(result)
        finally:
            if own_session:
                session.terminate()
    except ImportError:
        return False, "wolframclient not installed"
    except Exception as e:
        return False, str(e)


def run_dry(code: str) -> Tuple[bool, str]:
    """Dry-run: echo the expression without evaluation."""
    return True, f"[DRY-RUN] {code}"


# ── §5  DEMO ENGINE ──────────────────────────────────────────────────────

@dataclass
class DemoResult:
    """Result of a single demo expression evaluation."""
    label: str
    code: str
    success: bool
    output: str
    elapsed_ms: float


class WSTPKernelDemoRunner:
    """Execute live WSTP kernel demos against the generated .wl package.

    Supports three backends: wolframscript (subprocess), wolframclient
    (Python library), or dry-run (prints expressions without kernel).
    """

    def __init__(self, wl_path: str = WL_PACKAGE, backend: Optional[str] = None):
        self.wl_path = wl_path
        self.backend = backend or KernelBackend.detect()
        self._wc_session = None
        self.results: List[DemoResult] = []

    def _ensure_package(self):
        """Generate .wl package if missing."""
        if not os.path.exists(self.wl_path):
            try:
                from wstp_symbolic_exporter import WSTPSymbolicExporter
                exporter = WSTPSymbolicExporter(self.wl_path)
                exporter.export()
                print(f"  [auto-generated] {self.wl_path}")
            except ImportError:
                raise FileNotFoundError(
                    f"{self.wl_path} not found and wstp_symbolic_exporter.py unavailable"
                )

    def _eval(self, code: str) -> Tuple[bool, str]:
        """Route to the active backend."""
        if self.backend == KernelBackend.WOLFRAMCLIENT:
            return run_wolframclient(code, self._wc_session)
        elif self.backend == KernelBackend.WOLFRAMSCRIPT:
            return run_wolframscript(code)
        else:
            return run_dry(code)

    def run_demo(self) -> Dict[str, Any]:
        """Execute the full demo sequence and return results."""
        self._ensure_package()
        self.results.clear()

        expressions = build_demo_expressions(self.wl_path)

        # Open persistent session for wolframclient
        if self.backend == KernelBackend.WOLFRAMCLIENT:
            try:
                from wolframclient.evaluation import WolframLanguageSession
                self._wc_session = WolframLanguageSession()
                self._wc_session.start()
            except Exception as e:
                print(f"  [WARN] wolframclient session failed, falling back to dry-run: {e}")
                self.backend = KernelBackend.DRY_RUN

        t0_all = time.perf_counter()

        for expr in expressions:
            t0 = time.perf_counter()
            success, output = self._eval(expr["code"])
            elapsed = (time.perf_counter() - t0) * 1000

            self.results.append(DemoResult(
                label=expr["label"],
                code=expr["code"],
                success=success,
                output=output,
                elapsed_ms=round(elapsed, 2),
            ))

        total_elapsed = time.perf_counter() - t0_all

        # Close session
        if self._wc_session is not None:
            try:
                self._wc_session.terminate()
            except Exception:
                pass
            self._wc_session = None

        # GW validation cross-check (Python side)
        gw_validation = self._validate_gw_events()

        return {
            "backend": self.backend,
            "wl_package": os.path.abspath(self.wl_path),
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "total_elapsed_ms": round(total_elapsed * 1000, 2),
            "n_expressions": len(self.results),
            "n_success": sum(1 for r in self.results if r.success),
            "n_failed": sum(1 for r in self.results if not r.success),
            "results": [
                {
                    "label": r.label,
                    "code": r.code,
                    "success": r.success,
                    "output": r.output,
                    "elapsed_ms": r.elapsed_ms,
                }
                for r in self.results
            ],
            "gw_validation": gw_validation,
        }

    def _validate_gw_events(self) -> List[Dict[str, Any]]:
        """Python-side GW cross-validation against PAPER_001-009 values."""
        validations = []
        for name, ev in GW_EVENTS.items():
            d_total = ev["D_total"]
            h_predicted = ev["h_GR"] * d_total
            h_observed = ev["h_UQFF"]
            residual = abs(h_predicted - h_observed) / max(h_observed, 1e-50)
            validations.append({
                "event": name,
                "type": ev["type"],
                "D_total": d_total,
                "h_GR": ev["h_GR"],
                "h_UQFF_predicted": h_predicted,
                "h_UQFF_observed": h_observed,
                "residual": round(residual, 6),
                "pass": residual < 0.15,
            })
        return validations

    def export_results(self, report: Dict[str, Any],
                       filepath: str = "wstp_kernel_demo_results.json") -> str:
        """Export demo results to JSON."""
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2, default=str)
        return filepath


# ── §6  MAIN ──────────────────────────────────────────────────────────────

def main():
    """Run live WSTP kernel demo on the 11-sector UQFF Lagrangian."""
    print("=" * 72)
    print("WSTP Kernel Demo Runner — 9-Sector UQFF Lagrangian Live Session")
    print("=" * 72)

    backend = KernelBackend.detect()
    print(f"  Backend detected: {backend}")

    if backend == KernelBackend.DRY_RUN:
        print("  [INFO] No Wolfram kernel available — running in dry-run mode.")
        print("         Install wolframscript or 'pip install wolframclient' for live evaluation.")

    runner = WSTPKernelDemoRunner(backend=backend)

    print(f"\n── Running demo ({runner.wl_path}) ──\n")
    report = runner.run_demo()

    # Print results
    for r in report["results"]:
        status = "OK" if r["success"] else "FAIL"
        print(f"  [{status:4s}] {r['label']}")
        print(f"         Code:   {r['code'][:80]}{'...' if len(r['code']) > 80 else ''}")
        print(f"         Output: {r['output'][:100]}{'...' if len(r['output']) > 100 else ''}")
        print(f"         Time:   {r['elapsed_ms']:.2f} ms")
        print()

    # GW validation
    print("── GW Event Validation (Python cross-check) ──\n")
    for ev in report["gw_validation"]:
        status = "PASS" if ev["pass"] else "FAIL"
        print(f"  {ev['event']:10s} ({ev['type']}):  D={ev['D_total']:.3f}  "
              f"h_pred={ev['h_UQFF_predicted']:.4e}  "
              f"h_obs={ev['h_UQFF_observed']:.4e}  "
              f"resid={ev['residual']:.4f}  [{status}]")

    # Summary
    print(f"\n── Summary ──")
    print(f"  Backend:     {report['backend']}")
    print(f"  Expressions: {report['n_success']}/{report['n_expressions']} succeeded")
    print(f"  Total time:  {report['total_elapsed_ms']:.2f} ms")

    # Export
    out_path = runner.export_results(report)
    print(f"\n[OK] Results exported: {out_path}")

    print(f"\n{'=' * 72}")
    print(f"DEMO COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
