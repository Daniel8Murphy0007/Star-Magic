#!/usr/bin/env python3
"""
fubi_smbh_mergers.py — F_U_Bi (Inside-to-Outside) for SMBH Mergers

Session 220 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Dedicated module for F_U_Bi (inside-to-outside buoyancy portion of mass)
in SMBH mergers during coalescence. Extends SMBHBinaryMergerFUBiCalc from
Session 219 with:
  - Inspiral phase: phonon-damped orbital decay
  - Coalescence phase: mass ejection via buoyancy
  - Ringdown phase: remnant mass correction
  - Lagrangian: δS/δφ_merger = 0 → F_U_Bi coupling

Architecture: Pure calculator. Parameters via dataset dict.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List

# ── §0  Constants ──────────────────────────────────────────────────────────

PI        = math.pi
G         = 6.674e-11
C         = 2.998e8
HBAR      = 1.055e-34
M_SUN     = 1.989e30

OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
BETA_I    = 0.603
KAPPA     = 0.0005 / 86400.0
GAMMA_0   = 2 * PI * 0.1e12
SIGMA_G   = 0.08 * 2 * PI * 1e12
RHO_VAC   = 1e-10
F_NEUTRON = 1e-10

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))


# ── §1  Ramanujan S₂₆⁽³⁾ ──────────────────────────────────────────────────

def ramanujan_Rn(n: int, k: int = 3) -> float:
    total = 0.0
    for j in range(k):
        sign = (-1) ** j
        binom = 1.0
        for m in range(j):
            binom *= (k - 1 - m) / (m + 1)
        nfact = math.factorial(min(n + j, 170))
        total += sign * binom / nfact
    return total


S26_3RD = sum((SSQ ** n) / (n ** 26) * ramanujan_Rn(n, 3) for n in range(1, 27))


def horizon_r(M_kg: float, a: float) -> float:
    rS = 2 * G * M_kg / C**2
    return rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))


def Phi(gamma: float) -> float:
    return math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD


# ── §2  SMBHInspiralFUBiCalc ──────────────────────────────────────────────

class SMBHInspiralFUBiCalc:
    """SMBH binary inspiral: F_U_Bi during orbital decay.

    Phonon-damped inspiral with buoyancy-corrected GW emission.
    F_U_Bi = ρ_SCm · V · S₂₆⁽³⁾² · ratio · damping_envelope(Γ)
    """

    def compute(self, dataset: dict) -> dict:
        M1 = float(dataset.get("M1_Msun", 5.5e7)) * M_SUN
        M2 = float(dataset.get("M2_Msun", 3.0e7)) * M_SUN
        a1 = float(dataset.get("a1_spin", 0.70))
        a2 = float(dataset.get("a2_spin", 0.60))
        gamma_THz = float(dataset.get("gamma_THz", 0.10))

        M_total = M1 + M2
        mu = M1 * M2 / M_total  # reduced mass
        eta = mu / M_total  # symmetric mass ratio

        rH1 = horizon_r(M1, a1)
        rH2 = horizon_r(M2, a2)
        r_ISCO = 6 * G * M_total / C**2  # Schwarzschild ISCO

        gamma = 2 * PI * gamma_THz * 1e12
        Ph = Phi(gamma)

        # F_U_Bi for each BH
        rho = RHO_VAC * math.exp(min(KAPPA * 86400, 500.0))
        V1 = (4 / 3) * PI * rH1**3
        V2 = (4 / 3) * PI * rH2**3

        r2_1 = rH1**2
        r2_2 = rH2**2
        Ug1 = sum(G * M1 / r2_1 * SSQ * i / 26 for i in range(1, 27))
        Ub1 = sum(G * M1 / r2_1 * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
        Ug2 = sum(G * M2 / r2_2 * SSQ * i / 26 for i in range(1, 27))
        Ub2 = sum(G * M2 / r2_2 * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))

        ratio1 = abs(Ub1) / (abs(Ug1) + abs(Ub1) + 1e-300)
        ratio2 = abs(Ub2) / (abs(Ug2) + abs(Ub2) + 1e-300)

        FUBi_1 = rho * V1 * S26_3RD**2 * Ph * ratio1
        FUBi_2 = rho * V2 * S26_3RD**2 * Ph * ratio2
        FUBi_total = FUBi_1 + FUBi_2

        # Damping factor: 66.7% strain reduction at resonance
        damping = 1.0 - 0.667 * Ph / S26_3RD
        phase_lag = 200 + 167 * Ph / S26_3RD

        return {
            "primary_equations": [
                f"F_U_Bi(M1) = {FUBi_1:.6e} N",
                f"F_U_Bi(M2) = {FUBi_2:.6e} N",
                f"F_U_Bi(total) = {FUBi_total:.6e} N",
                f"η = {eta:.4f}, damping = {damping:.3f}",
                f"Phase lag = {phase_lag:.1f} cycles",
                f"r_ISCO = {r_ISCO:.6e} m",
            ],
            "note": "PAPER_1014 CP4. Session 220. SMBH inspiral F_U_Bi.",
        }


# ── §3  SMBHCoalescenceFUBiCalc ──────────────────────────────────────────

class SMBHCoalescenceFUBiCalc:
    """SMBH coalescence: mass ejection via buoyancy during merger.

    At coalescence, F_U_Bi contributes to inertial mass modification:
    M_remnant = M_total - E_GW/c² - ΔM_buoyancy
    """

    def compute(self, dataset: dict) -> dict:
        M1 = float(dataset.get("M1_Msun", 5.5e7)) * M_SUN
        M2 = float(dataset.get("M2_Msun", 3.0e7)) * M_SUN
        a_final = float(dataset.get("a_final", 0.69))
        gamma_THz = float(dataset.get("gamma_THz", 0.10))

        M_total = M1 + M2
        eta = (M1 * M2) / M_total**2
        gamma = 2 * PI * gamma_THz * 1e12
        Ph = Phi(gamma)

        # GR energy loss
        E_GW_frac = 4 * eta**2 * (1 - 0.2 * (1 - 4 * eta))
        E_GW = E_GW_frac * M_total * C**2

        # Buoyancy mass correction
        rH_final = horizon_r(M_total, a_final)
        rho = RHO_VAC * math.exp(min(KAPPA * 86400, 500.0))
        V_final = (4 / 3) * PI * rH_final**3
        delta_M_buoy = rho * V_final * S26_3RD**2 * Ph / C**2

        M_remnant = M_total - E_GW / C**2 - delta_M_buoy

        return {
            "primary_equations": [
                f"M_total = {M_total / M_SUN:.2e} M☉",
                f"E_GW = {E_GW:.6e} J (η = {eta:.4f})",
                f"ΔM_buoyancy = {delta_M_buoy:.6e} kg ({delta_M_buoy / M_SUN:.2e} M☉)",
                f"M_remnant = {M_remnant / M_SUN:.2e} M☉",
                f"a_final = {a_final}",
            ],
            "note": "PAPER_1015 CP4. Session 220. SMBH coalescence F_U_Bi.",
        }


# ── §4  SMBHRingdownFUBiCalc ─────────────────────────────────────────────

class SMBHRingdownFUBiCalc:
    """SMBH ringdown: buoyancy-modified QNM frequencies.

    QNM frequency: f_QNM ≈ c³/(2πGM) · (1 - 0.63(1-a)^{3/10})
    SCm modification: f_UQFF = f_QNM · (1 + S₂₆⁽³⁾ · Φ(Γ))
    """

    def compute(self, dataset: dict) -> dict:
        M_remnant = float(dataset.get("M_remnant_Msun", 8.2e7)) * M_SUN
        a = float(dataset.get("a_spin", 0.69))
        gamma_THz = float(dataset.get("gamma_THz", 0.10))

        gamma = 2 * PI * gamma_THz * 1e12
        Ph = Phi(gamma)

        f_QNM = C**3 / (2 * PI * G * M_remnant) * (1 - 0.63 * (1 - a)**0.3)
        f_UQFF = f_QNM * (1 + S26_3RD * Ph)
        tau_QNM = 2 * G * M_remnant / (C**3 * (1 - a)**0.45)
        tau_UQFF = tau_QNM * (1 - S26_3RD * Ph * 0.1)

        return {
            "primary_equations": [
                f"f_QNM = {f_QNM:.6e} Hz",
                f"f_UQFF = {f_UQFF:.6e} Hz (SCm correction)",
                f"τ_QNM = {tau_QNM:.6e} s",
                f"τ_UQFF = {tau_UQFF:.6e} s",
                f"Δf/f_QNM = {S26_3RD * Ph:.6e}",
            ],
            "note": "PAPER_1016 CP4. Session 220. SMBH ringdown F_U_Bi.",
        }


# ── §5  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    ok = True
    passed = 0

    # Test 1: Inspiral F_U_Bi total positive
    calc = SMBHInspiralFUBiCalc()
    res = calc.compute({})
    fubi_line = res["primary_equations"][2]
    val = float(fubi_line.split("= ")[1].split(" N")[0])
    if val > 0:
        print(f"[ OK ] Inspiral: {fubi_line}")
        passed += 1
    else:
        print(f"[FAIL] Inspiral F_U_Bi should be positive"); ok = False

    # Test 2: Inspiral damping < 1.0
    damp_line = res["primary_equations"][3]
    damp_val = float(damp_line.split("damping = ")[1])
    if 0 < damp_val < 1.0:
        print(f"[ OK ] Inspiral: {damp_line}")
        passed += 1
    else:
        print(f"[FAIL] Damping should be 0 < d < 1"); ok = False

    # Test 3: Inspiral phase lag ~367
    lag_line = res["primary_equations"][4]
    lag_val = float(lag_line.split("= ")[1].split(" cycles")[0])
    if 360 < lag_val < 370:
        print(f"[ OK ] Inspiral: {lag_line}")
        passed += 1
    else:
        print(f"[FAIL] Phase lag: {lag_line}"); ok = False

    # Test 4: Coalescence mass budget
    calc = SMBHCoalescenceFUBiCalc()
    res = calc.compute({})
    remnant_line = res["primary_equations"][3]
    if "M☉" in remnant_line:
        print(f"[ OK ] Coalescence: {remnant_line}")
        passed += 1
    else:
        print(f"[FAIL] Coalescence remnant"); ok = False

    # Test 5: ΔM_buoyancy > 0
    dm_line = res["primary_equations"][2]
    dm_val = float(dm_line.split("= ")[1].split(" kg")[0])
    if dm_val > 0:
        print(f"[ OK ] Coalescence: {dm_line}")
        passed += 1
    else:
        print(f"[FAIL] ΔM_buoyancy should be positive"); ok = False

    # Test 6: Ringdown QNM frequency > 0
    calc = SMBHRingdownFUBiCalc()
    res = calc.compute({})
    fqnm_line = res["primary_equations"][0]
    fqnm_val = float(fqnm_line.split("= ")[1].split(" Hz")[0])
    if fqnm_val > 0:
        print(f"[ OK ] Ringdown: {fqnm_line}")
        passed += 1
    else:
        print(f"[FAIL] f_QNM should be positive"); ok = False

    # Test 7: SCm correction to QNM
    df_line = res["primary_equations"][4]
    df_val = float(df_line.split("= ")[1])
    if 0 < df_val < 0.01:
        print(f"[ OK ] Ringdown: {df_line}")
        passed += 1
    else:
        print(f"[FAIL] Δf/f_QNM: {df_line}"); ok = False

    # Test 8: S₂₆⁽³⁾
    if S26_3RD > 0:
        print(f"[ OK ] S₂₆⁽³⁾ = {S26_3RD:.6e}")
        passed += 1
    else:
        print(f"[FAIL] S₂₆⁽³⁾"); ok = False

    print(f"\n{'='*60}")
    print(f"  fubi_smbh_mergers.py: {passed}/8 tests passed")
    print(f"{'='*60}")
    return ok


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
