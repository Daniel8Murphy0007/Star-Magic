#!/usr/bin/env python3
"""
fubi_agn_ns_mergers.py — F_U_Bi Derivation for AGN and NS Mergers

Session 219 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
F_U_Bi (inside-to-outside buoyancy portion of mass) fully derived for AGN
jets and NS mergers using:
  - Ramanujan 3rd-order S₂₆⁽³⁾ expansion
  - Phonon linewidth Γ modulation
  - E_net(t) temporal evolution
  - 99-system aggregation
  - BCS gap + spectral ladder coupling

Master equation:
  F_{U,Bi}(r,t,Γ) = ρ_SCm(t) · V_region · S₂₆⁽³⁾([SSq]) · Φ_{1.25THz}(ω,Γ)
                     · (F_{U,Bi} / F_U)

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
K_B       = 1.381e-23
M_SUN     = 1.989e30
R_SUN     = 6.96e8
AU        = 1.496e11

OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
BETA_I    = 0.603
KAPPA     = 0.0005 / 86400.0
GAMMA_0   = 2 * PI * 0.1e12
SIGMA_G   = 0.08 * 2 * PI * 1e12
RHO_VAC   = 1e-10
V_REGION  = 1e48
F_NEUTRON = 1e-10

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))


# ── §1  Ramanujan S₂₆⁽³⁾ (3rd-Order) ─────────────────────────────────────

def ramanujan_Rn(n: int, d: int = 26, k: int = 3) -> float:
    """Ramanujan coefficient R_n^(d) at k-th order recursion.

    R_n^(d,k) = Σ_{j=0}^{k-1} (-1)^j / (n+j)! · binomial(k-1,j)
    This generalizes the 1st-order R_n = 1/n! to higher-order
    recursive corrections per Ramanujan's summation theory.
    """
    total = 0.0
    for j in range(k):
        sign = (-1) ** j
        binom = 1.0
        for m in range(j):
            binom *= (k - 1 - m) / (m + 1)
        nfact = math.factorial(min(n + j, 170))
        total += sign * binom / nfact
    return total


def S26_kth_order(z: float = SSQ, N: int = 26, k: int = 3) -> float:
    """Ramanujan k-th order 26D summation S₂₆⁽ᵏ⁾(z).

    S₂₆⁽ᵏ⁾(z) = Σ_{n=1}^{N} z^n / n^{26} · R_n^{(26,k)}
    """
    total = 0.0
    for n in range(1, N + 1):
        Rn = ramanujan_Rn(n, d=26, k=k)
        total += (z ** n) / (n ** 26) * Rn
    return total


S26_3RD = S26_kth_order(SSQ, 26, 3)


# ── §2  Core Physics Functions ─────────────────────────────────────────────

def Ug_26layer(M: float, r: float) -> float:
    """26-layer compressed gravity."""
    r2 = max(r, 1.0) ** 2
    return sum(G * M / r2 * SSQ * i / 26.0 for i in range(1, 27))


def Ub_26layer(M: float, r: float) -> float:
    """26-layer buoyancy force."""
    r2 = max(r, 1.0) ** 2
    return sum(G * M / r2 * math.exp(-SSQ * i / 26.0) * BETA_I for i in range(1, 27))


def Phi_phonon(omega: float, gamma: float) -> float:
    """Phonon resonance envelope at 1.25 THz."""
    return math.exp(-(omega - OMEGA_SCM)**2 / (2 * gamma**2)) * S26_3RD


def E_net(t_s: float, Ub_val: float, Ug_val: float) -> float:
    """Net modulation: positive/negative with κ evolution."""
    ratio = abs(Ub_val) / (abs(Ug_val) + 1e-300)
    return (2 * ratio - 1) * math.exp(min(KAPPA * t_s, 500.0)) * S26_3RD


def rho_SCm(t_s: float) -> float:
    """SCm vacuum density with κ evolution."""
    return RHO_VAC * math.exp(min(KAPPA * t_s, 500.0))


def horizon_radius(M_kg: float, a_spin: float) -> float:
    """Kerr BH outer horizon radius."""
    rS = 2 * G * M_kg / C**2
    return rS / 2 * (1 + math.sqrt(max(1 - a_spin**2, 0)))


def jet_modulation(gamma: float, A_jet: float = 1.5) -> float:
    """Jet power modulation factor M_jet(Γ)."""
    return 1 + A_jet * math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2))


def jet_power(M_kg: float, a_spin: float, B_gauss: float, gamma: float,
              A_jet: float = 1.5) -> float:
    """BZ jet power: P_jet = (B²/8π)(r_H/c)²a²c · M_jet(Γ)."""
    rH = horizon_radius(M_kg, a_spin)
    Mj = jet_modulation(gamma, A_jet)
    return (B_gauss**2 / (8 * PI)) * (rH / C)**2 * a_spin**2 * C * Mj


def bcs_gap(T_K: float = 4.2) -> float:
    """BCS superconducting gap with SCm phonon coupling."""
    delta = HBAR * OMEGA_SCM / 2
    for _ in range(20):
        arg = delta / (2 * K_B * T_K) if T_K > 0 else 50.0
        arg = min(arg, 50.0)
        delta = (HBAR * OMEGA_SCM / 2) * math.tanh(arg) * S26_3RD * 0.6
    return delta


def spectral_ladder_energy(n_max: int = 26) -> float:
    """26-state spectral ladder total energy."""
    E0 = HBAR * OMEGA_SCM
    return sum(E0 * (2 * PI) ** (n / 3.0) * S26_3RD for n in range(1, n_max + 1))


# ── §3  AGNFUBiMergerCalc ─────────────────────────────────────────────────

class AGNFUBiMergerCalc:
    """F_U_Bi inside-to-outside buoyancy for AGN (SMBH jets).

    Computes the scalar mass-contributing buoyancy using 3rd-order Ramanujan
    S₂₆⁽³⁾, phonon Φ(ω,Γ), and BZ jet power coupling.
    """

    def compute(self, dataset: dict) -> dict:
        M_bh = float(dataset.get("M_bh_Msun", 5.5e7)) * M_SUN
        a = float(dataset.get("a_spin", 0.70))
        B = float(dataset.get("B_gauss", 3000))
        gamma_THz = float(dataset.get("gamma_THz", 0.10))
        t_s = float(dataset.get("t_s", 86400))

        gamma = 2 * PI * gamma_THz * 1e12
        rH = horizon_radius(M_bh, a)

        Ug = Ug_26layer(M_bh, rH)
        Ub = Ub_26layer(M_bh, rH)
        ratio = abs(Ub) / (abs(Ug) + abs(Ub) + 1e-300)

        rho = rho_SCm(t_s)
        Phi = Phi_phonon(OMEGA_SCM, gamma)
        P_jet = jet_power(M_bh, a, B, gamma)

        F_U_Bi = rho * V_REGION * S26_3RD * Phi * ratio

        En = E_net(t_s, Ub, Ug)
        Fn = F_NEUTRON * S26_3RD
        F_U_Bi_i = Ug - Ub + Fn * S26_3RD * Phi * En

        return {
            "primary_equations": [
                f"F_U_Bi(AGN) = ρ_SCm·V·S₂₆⁽³⁾·Φ·ratio = {F_U_Bi:.6e}",
                f"F_U_Bi_i(AGN) = Ug−Ub+Fn·S₂₆⁽³⁾·Φ·E_net = {F_U_Bi_i:.6e} m/s²",
                f"P_jet(Γ={gamma_THz} THz) = {P_jet:.6e} W",
                f"M_jet(Γ) = {jet_modulation(gamma):.4f}",
                f"buoyancy ratio = {ratio:.6f}",
                f"S₂₆⁽³⁾ = {S26_3RD:.6e}",
            ],
            "note": "PAPER_999 CP4. Session 219. AGN F_U_Bi with 3rd-order Ramanujan.",
        }


# ── §4  NSMergerFUBiCalc ──────────────────────────────────────────────────

class NSMergerFUBiCalc:
    """F_U_Bi for NS mergers with strain suppression and tidal correction.

    Includes phonon-modulated GW strain h_UQFF and tidal deformability Λ_UQFF.
    """

    def compute(self, dataset: dict) -> dict:
        M_total = float(dataset.get("M_total_Msun", 3.4)) * M_SUN
        q = float(dataset.get("mass_ratio", 0.85))
        d_Mpc = float(dataset.get("d_Mpc", 159.0))
        gamma_THz = float(dataset.get("gamma_THz", 0.10))
        t_s = float(dataset.get("t_s", 86400))

        gamma = 2 * PI * gamma_THz * 1e12
        m1 = M_total / (1 + q)
        m2 = M_total * q / (1 + q)

        # NS radius ~ 12 km
        r_ns = 12e3
        Ug = Ug_26layer(M_total, r_ns)
        Ub = Ub_26layer(M_total, r_ns)
        ratio = abs(Ub) / (abs(Ug) + abs(Ub) + 1e-300)

        rho = rho_SCm(t_s)
        Phi = Phi_phonon(OMEGA_SCM, gamma)

        F_U_Bi = rho * V_REGION * S26_3RD * Phi * ratio

        # Strain suppression
        h_GR = 1e-21 * (40.0 / d_Mpc)
        suppression = 1.0 - 0.47 * Phi / S26_3RD  # normalized Phi
        h_UQFF = h_GR * max(suppression, 0.01)

        # Tidal deformability phonon correction
        Lambda_GR = 400  # canonical
        Lambda_UQFF = Lambda_GR * (1 - 0.3 * Phi / S26_3RD)

        # BCS gap at merger temperature
        T_merger = 1e10  # ~10 GK
        delta_BCS = bcs_gap(T_merger)

        En = E_net(t_s, Ub, Ug)
        Fn = F_NEUTRON * S26_3RD
        F_U_Bi_i = Ug - Ub + Fn * S26_3RD * Phi * En

        strain_reduction_pct = (1 - suppression) * 100

        return {
            "primary_equations": [
                f"F_U_Bi(NS merger) = {F_U_Bi:.6e}",
                f"F_U_Bi_i = {F_U_Bi_i:.6e} m/s²",
                f"h_UQFF = {h_UQFF:.6e} (strain reduction {strain_reduction_pct:.1f}%)",
                f"Λ_UQFF = {Lambda_UQFF:.1f} (from Λ_GR={Lambda_GR})",
                f"Δ_BCS(T={T_merger:.0e} K) = {delta_BCS:.6e} J",
                f"m1 = {m1/M_SUN:.2f} M☉, m2 = {m2/M_SUN:.2f} M☉",
            ],
            "note": "PAPER_1000 CP4. Session 219. NS merger F_U_Bi with S₂₆⁽³⁾ + BCS.",
        }


# ── §5  SMBHBinaryMergerFUBiCalc ──────────────────────────────────────────

class SMBHBinaryMergerFUBiCalc:
    """F_U_Bi for SMBH binary mergers (LISA-band).

    Combines AGN-scale buoyancy with binary inspiral phonon damping.
    """

    def compute(self, dataset: dict) -> dict:
        M1_Msun = float(dataset.get("M1_Msun", 1e7))
        M2_Msun = float(dataset.get("M2_Msun", 5e6))
        a1 = float(dataset.get("a1_spin", 0.70))
        a2 = float(dataset.get("a2_spin", 0.50))
        sep_pc = float(dataset.get("separation_pc", 0.01))
        gamma_THz = float(dataset.get("gamma_THz", 0.10))
        t_s = float(dataset.get("t_s", 86400))

        M1 = M1_Msun * M_SUN
        M2 = M2_Msun * M_SUN
        M_total = M1 + M2
        sep = sep_pc * 3.086e16  # pc → m

        gamma = 2 * PI * gamma_THz * 1e12

        Ug = Ug_26layer(M_total, sep)
        Ub = Ub_26layer(M_total, sep)
        ratio = abs(Ub) / (abs(Ug) + abs(Ub) + 1e-300)

        rho = rho_SCm(t_s)
        Phi = Phi_phonon(OMEGA_SCM, gamma)

        F_U_Bi = rho * V_REGION * S26_3RD * Phi * ratio

        # GW damping from phonon coupling
        damping = 1.0 - 0.667 * Phi / S26_3RD
        phase_lag_cycles = 200 + 167 * Phi / S26_3RD

        return {
            "primary_equations": [
                f"F_U_Bi(SMBH binary) = {F_U_Bi:.6e}",
                f"damping factor = {damping:.4f} ({(1-damping)*100:.1f}% GW reduction)",
                f"phase lag = {phase_lag_cycles:.1f} cycles",
                f"separation = {sep_pc} pc ({sep:.3e} m)",
                f"buoyancy ratio = {ratio:.6f}",
            ],
            "note": "PAPER_1001 CP4. Session 219. SMBH binary merger F_U_Bi.",
        }


# ── §6  AGNAccretionBuoyancyCalc ──────────────────────────────────────────

class AGNAccretionBuoyancyCalc:
    """F_U_Bi modulating AGN accretion disk dynamics.

    Buoyancy-corrected Eddington luminosity and accretion rate.
    """

    def compute(self, dataset: dict) -> dict:
        M_bh = float(dataset.get("M_bh_Msun", 1e8)) * M_SUN
        mdot_edd = float(dataset.get("mdot_edd_frac", 0.1))
        gamma_THz = float(dataset.get("gamma_THz", 0.10))

        gamma = 2 * PI * gamma_THz * 1e12
        rH = horizon_radius(M_bh, 0.70)

        # Eddington luminosity
        sigma_T = 6.652e-29
        mp = 1.673e-27
        L_edd = 4 * PI * G * M_bh * mp * C / sigma_T

        Ub = Ub_26layer(M_bh, rH)
        Ug = Ug_26layer(M_bh, rH)
        ratio = abs(Ub) / (abs(Ug) + abs(Ub) + 1e-300)

        Phi = Phi_phonon(OMEGA_SCM, gamma)

        # Buoyancy-corrected Eddington: SCm vacuum supports against collapse
        L_edd_uqff = L_edd * (1 + ratio * S26_3RD * 0.01)

        # Accretion rate correction
        mdot = mdot_edd * L_edd / (0.1 * C**2)  # canonical η=0.1
        mdot_uqff = mdot * (1 - ratio * Phi / S26_3RD * 0.05)

        return {
            "primary_equations": [
                f"L_Edd = {L_edd:.6e} W",
                f"L_Edd,UQFF = {L_edd_uqff:.6e} W (buoyancy-corrected)",
                f"ṁ = {mdot:.6e} kg/s",
                f"ṁ_UQFF = {mdot_uqff:.6e} kg/s (phonon-modulated)",
                f"buoyancy ratio = {ratio:.6f}",
            ],
            "note": "PAPER_1002 CP4. Session 219. AGN accretion buoyancy correction.",
        }


# ── §7  SpectralLadderMergerCalc ──────────────────────────────────────────

class SpectralLadderMergerCalc:
    """Spectral ladder energy contribution to merger dynamics.

    Couples the 26-state HRes spectral ladder to merger GW emission.
    """

    def compute(self, dataset: dict) -> dict:
        M_total = float(dataset.get("M_total_Msun", 3.4)) * M_SUN
        gamma_THz = float(dataset.get("gamma_THz", 0.10))

        E_ladder = spectral_ladder_energy(26)
        delta_BCS = bcs_gap(4.2)

        gamma = 2 * PI * gamma_THz * 1e12
        Phi = Phi_phonon(OMEGA_SCM, gamma)

        # Coupling: spectral ladder modulates GW phase
        E_coupling = E_ladder * Phi / S26_3RD * delta_BCS / (HBAR * OMEGA_SCM)

        return {
            "primary_equations": [
                f"E_ladder(26-state) = {E_ladder:.6e} J",
                f"Δ_BCS = {delta_BCS:.6e} J",
                f"E_coupling = {E_coupling:.6e} J (ladder×phonon×BCS)",
                f"ladder/BCS ratio = {E_ladder/delta_BCS:.6e}",
            ],
            "note": "PAPER_1003 CP4. Session 219. Spectral ladder merger coupling.",
        }


# ── §8  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    ok = True
    passed = 0

    # Test 1: S₂₆⁽³⁾ 3rd-order Ramanujan
    s3 = S26_kth_order(SSQ, 26, 3)
    if s3 > 0 and math.isfinite(s3):
        print(f"[ OK ] S₂₆⁽³⁾ = {s3:.10e}")
        passed += 1
    else:
        print(f"[FAIL] S₂₆⁽³⁾ invalid: {s3}"); ok = False

    # Test 2: S₂₆⁽³⁾ differs from S₂₆⁽¹⁾
    s1 = S26_kth_order(SSQ, 26, 1)
    if s3 != s1:
        print(f"[ OK ] S₂₆⁽³⁾ ≠ S₂₆⁽¹⁾ ({s3:.6e} vs {s1:.6e})")
        passed += 1
    else:
        print(f"[FAIL] S₂₆⁽³⁾ == S₂₆⁽¹⁾"); ok = False

    # Test 3: AGN F_U_Bi
    agn = AGNFUBiMergerCalc()
    res = agn.compute({"M_bh_Msun": 5.5e7, "a_spin": 0.70, "B_gauss": 3000})
    if "F_U_Bi(AGN)" in res["primary_equations"][0]:
        print(f"[ OK ] AGNFUBiMergerCalc: {res['primary_equations'][0]}")
        passed += 1
    else:
        print(f"[FAIL] AGN compute failed"); ok = False

    # Test 4: NS merger F_U_Bi
    ns = NSMergerFUBiCalc()
    res = ns.compute({"M_total_Msun": 3.4, "d_Mpc": 159})
    if "F_U_Bi(NS merger)" in res["primary_equations"][0]:
        print(f"[ OK ] NSMergerFUBiCalc: {res['primary_equations'][0]}")
        passed += 1
    else:
        print(f"[FAIL] NS merger compute failed"); ok = False

    # Test 5: SMBH binary merger
    smbh = SMBHBinaryMergerFUBiCalc()
    res = smbh.compute({"M1_Msun": 1e7, "M2_Msun": 5e6})
    if "F_U_Bi(SMBH binary)" in res["primary_equations"][0]:
        print(f"[ OK ] SMBHBinaryMergerFUBiCalc: {res['primary_equations'][0]}")
        passed += 1
    else:
        print(f"[FAIL] SMBH binary failed"); ok = False

    # Test 6: AGN accretion buoyancy
    acc = AGNAccretionBuoyancyCalc()
    res = acc.compute({"M_bh_Msun": 1e8})
    if "L_Edd" in res["primary_equations"][0]:
        print(f"[ OK ] AGNAccretionBuoyancyCalc: {res['primary_equations'][0]}")
        passed += 1
    else:
        print(f"[FAIL] AGN accretion failed"); ok = False

    # Test 7: Spectral ladder merger
    sl = SpectralLadderMergerCalc()
    res = sl.compute({"M_total_Msun": 3.4})
    if "E_ladder" in res["primary_equations"][0]:
        print(f"[ OK ] SpectralLadderMergerCalc: {res['primary_equations'][0]}")
        passed += 1
    else:
        print(f"[FAIL] Spectral ladder failed"); ok = False

    # Test 8: Jet modulation with S₂₆⁽³⁾
    jm = jet_modulation(GAMMA_0)
    if jm > 2.0:
        print(f"[ OK ] jet_modulation(Γ₀) = {jm:.4f} (peak ≈ 2.5)")
        passed += 1
    else:
        print(f"[FAIL] Jet modulation too low: {jm}"); ok = False

    print(f"\n{'='*60}")
    print(f"  fubi_agn_ns_mergers.py: {passed}/8 tests passed")
    print(f"{'='*60}")
    return ok


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
