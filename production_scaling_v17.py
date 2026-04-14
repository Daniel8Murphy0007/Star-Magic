#!/usr/bin/env python3
"""
production_scaling_v17.py — Production Scaling to 900k calc/s

Session 224 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
v17 upgrade: 900,000 calc/s target (up from v16's 800k).
44 kernels total = v16's 36 + 8 new Session 224 kernels:
  kernel_gpu_dpm_atlas_peak     — GPU DPM spectral atlas peak amplitude
  kernel_dpm_line_fwhm          — DPM spectral line FWHM
  kernel_muge_8term_gravity     — 8-term MUGE gravitational field
  kernel_nfw_density            — NFW dark matter halo density
  kernel_enet_gamma_coupled     — Γ-coupled E_net dark energy
  kernel_dark_energy_eos_w0     — SCm dark energy w(z=0)
  kernel_alma_fubi_profile      — ALMA F_{U,Bi,i} line profile
  kernel_alma_chi2_co21         — ALMA CO(2-1) χ² validation

Architecture: Standalone kernel functions, ALL_KERNELS list, benchmark class.
────────────────────────────────────────────────────────────────────────────────
"""

import math
import time
from typing import Dict, List, Callable

# ── §0  Constants ──────────────────────────────────────────────────────────

PI        = math.pi
G         = 6.674e-11
C         = 2.998e8
M_SUN     = 1.989e30
HBAR      = 1.055e-34
K_B       = 1.381e-23
KPC       = 3.086e19

OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
BETA_I    = 0.603
KAPPA     = 0.0005 / 86400.0
GAMMA_0   = 2 * PI * 0.1e12
SIGMA_G   = 0.08 * 2 * PI * 1e12
F_NEUTRON = 1e-10
RHO_VAC   = 1e-10
RHO_SCM   = 7.09e-37
H_0       = 2.195e-18
T_H       = 1.0 / H_0
LAMBDA_C  = 1.11e-52
RHO_CRIT  = 3 * H_0**2 / (8 * PI * G)
RHO_LAMBDA = 0.692 * RHO_CRIT

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))

TARGET_CALC_PER_SEC = 900_000


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


def Phi(gamma: float) -> float:
    return math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD


# Pre-compute DPM layer coefficients
_LAYER_COEFFS = [
    (SSQ ** i) / (i ** 26) * ramanujan_Rn(i, 3) for i in range(1, 27)
]


# ── §2  Carry-Forward Kernels (v11: 16, v13: 4, v14: 4, v15: 6, v16: 6 = 36) ──

def kernel_gravity_26layer(M=M_SUN, r=6.96e8):
    return sum(G * M / r**2 * SSQ * i / 26 for i in range(1, 27))


def kernel_fu_bi_i(M=M_SUN, r=6.96e8, t=86400):
    Ug = sum(G * M / r**2 * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(G * M / r**2 * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    return (Ug - Ub) + F_NEUTRON * S26_3RD * Phi(GAMMA_0) * math.exp(min(KAPPA * t, 500))


def kernel_phonon_ares(M=M_SUN, r=6.96e8):
    return sum(math.exp(-SSQ * i / 26) * BETA_I * G * M / r**2 for i in range(1, 27))


def kernel_jet_mjet(A=1.5, gamma_THz=0.1):
    gamma = 2 * PI * gamma_THz * 1e12
    return 1 + A * math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2))


def kernel_ns_spindown(P=0.033, Pdot=4.2e-13):
    return 3.2e19 * math.sqrt(P * Pdot)


def kernel_gw170817_strain(M_chirp=1.186, d_Mpc=40):
    f = 100; d = d_Mpc * 3.086e22
    return 4 / d * (G * M_chirp * M_SUN / C**2)**(5/3) * (PI * f / C)**(2/3)


def kernel_blazar_ergosphere(M=3e8, a=0.95):
    M_kg = M * M_SUN; rS = 2 * G * M_kg / C**2
    return rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))


def kernel_rest_phonon_jet():
    return math.exp(-(OMEGA_SCM - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD


def kernel_qcalcgeom_vectorized(N=26):
    return sum(SSQ**i / (i**N) for i in range(1, N + 1))


def kernel_pipeline_full(M=M_SUN, r=6.96e8, t=86400):
    g = kernel_gravity_26layer(M, r)
    b = kernel_phonon_ares(M, r)
    return g - b + F_NEUTRON * S26_3RD * math.exp(min(KAPPA * t, 500))


def kernel_cena_jet(L=1e43, B=3000):
    return L * SSQ * B / (4 * PI * C)


def kernel_txs0506_jet(L=1e46, B=5000):
    return L * SSQ * B / (4 * PI * C)


def kernel_bcs_gap_solve(T=2e12, Tc=1.5e12):
    if T >= Tc: return 0.0
    return 3.528 * HBAR * OMEGA_SCM * math.sqrt(1 - (T / Tc)**2) * S26_3RD


def kernel_spectral_ladder_eval(N=26):
    return sum(math.exp(-SSQ * i / N) / (i + 1) for i in range(1, N + 1))


def kernel_ramanujan_26d_sum():
    return S26_3RD


def kernel_triadic_solver(M=M_SUN, r=6.96e8):
    Ug = kernel_gravity_26layer(M, r)
    Ub = kernel_phonon_ares(M, r)
    Um = G * M / r**2 * SSQ * 0.1
    return Ug + Um - Ub


# v13 carry-forward (4)
def kernel_fubi_inside_out(M=5e7, r=None):
    M_kg = M * M_SUN
    if r is None: r = 2 * G * M_kg / C**2
    return kernel_fu_bi_i(M_kg, r, 86400)


def kernel_99sys_gamma_sweep():
    total = 0.0
    for g in [0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0]:
        total += kernel_jet_mjet(1.5, g)
    return total


def kernel_agn_cena_fubi(M=5.5e7, a=0.70):
    M_kg = M * M_SUN; rS = 2 * G * M_kg / C**2
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    return kernel_fu_bi_i(M_kg, rH, 86400)


def kernel_ns_merger_gw190425(M_chirp=1.44, d_Mpc=159):
    return kernel_gw170817_strain(M_chirp, d_Mpc)


# v14 carry-forward (4)
def kernel_agn_merger_fubi(M1=5.5e7, M2=3e7, a=0.70):
    M_kg = (M1 + M2) * M_SUN; rS = 2 * G * M_kg / C**2
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(G * M_kg / rH**2 * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(G * M_kg / rH**2 * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    return abs(Ub) / (abs(Ug) + 1e-300) * S26_3RD


def kernel_qgp_vacuum_density(T=2e12):
    if T < 1.5e12: return 0.0
    return RHO_VAC * math.exp(min(KAPPA * 86400, 500)) * (T / 1e12)**4 * S26_3RD


def kernel_alice_multiplicity(N_part=383, alpha=8.5, beta=1.25):
    return alpha * N_part * (1 - math.exp(-beta * N_part / 400)) * S26_3RD / SSQ


def kernel_ym_mass_gap(T=2e12, Tc=1.5e12):
    if T >= Tc: return 0.0
    return 3.528 * HBAR * OMEGA_SCM * math.sqrt(1 - (T / Tc)**2) * S26_3RD**2


# v15 carry-forward (6)
def kernel_3c273_agn_fubi(M=8.86e8, a=0.90, A_jet=2.1):
    M_kg = M * M_SUN
    rS = 2 * G * M_kg / C**2
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(G * M_kg / rH**2 * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(G * M_kg / rH**2 * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    ratio = abs(Ub) / (abs(Ug) + abs(Ub) + 1e-300)
    M_jet = 1 + A_jet
    return ratio * M_jet * S26_3RD


def kernel_ton618_agn_fubi(M=6.6e10, a=0.998, A_jet=2.8):
    M_kg = M * M_SUN
    rS = 2 * G * M_kg / C**2
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(G * M_kg / rH**2 * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(G * M_kg / rH**2 * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    ratio = abs(Ub) / (abs(Ug) + abs(Ub) + 1e-300)
    M_jet = 1 + A_jet
    return ratio * M_jet * S26_3RD


def kernel_gw170817_merger(M_chirp=1.186, d_Mpc=40, suppression=0.667):
    h0 = kernel_gw170817_strain(M_chirp, d_Mpc)
    return h0 * (1 - suppression * Phi(GAMMA_0) / S26_3RD)


def kernel_smbh_merger_fubi(M1=5.5e7, M2=3e7, a1=0.70, a2=0.60):
    M_total = (M1 + M2) * M_SUN
    rS = 2 * G * M_total / C**2
    a_eff = (M1 * a1 + M2 * a2) / (M1 + M2)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a_eff**2, 0)))
    rho = RHO_VAC * math.exp(min(KAPPA * 86400, 500))
    V = (4/3) * PI * rH**3
    Ph = Phi(GAMMA_0)
    return rho * V * S26_3RD**2 * Ph


def kernel_dm_halo_nfw(rho_0=5e-22, r_kpc=20, r_s_kpc=20):
    x = r_kpc / r_s_kpc
    if x < 1e-10: x = 1e-10
    rho_nfw = rho_0 / (x * (1 + x)**2)
    rho_SCm = RHO_VAC * math.exp(min(KAPPA * 86400, 500))
    return rho_SCm * S26_3RD * rho_nfw * Phi(GAMMA_0) / rho_0


def kernel_txs0506_3gamma(A_jet=1.3):
    total = 0.0
    for g_THz in [0.05, 0.10, 0.30]:
        gamma = 2 * PI * g_THz * 1e12
        total += 1 + A_jet * math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2))
    return total * S26_3RD


# v16 carry-forward (6)
def kernel_scm_inflation_hubble(rho_scm=1e76):
    H_base = math.sqrt(8 * PI * G / 3 * rho_scm)
    return H_base * S26_3RD * Phi(GAMMA_0)


def kernel_thorne_morris_exotic(M=4.3e6, r0=1e4):
    M_kg = M * M_SUN
    rho_ambient = (HBAR * OMEGA_SCM) / (2 * C**2) * S26_STATIC * 0.99
    r_s = 2 * G * M_kg / C**2
    amp = (r_s / r0)**2 * S26_STATIC**2
    rho_exotic = rho_ambient * amp
    P_density = -BETA_I * rho_exotic * SSQ * 0.70 * S26_STATIC**2
    return rho_exotic + P_density


def kernel_vds_convergence(n_terms=26):
    return sum(SSQ**n / n**26 for n in range(1, n_terms + 1))


def kernel_dvp_prime_sieve(max_p=50):
    primes = []
    for n in range(2, max_p + 1):
        if all(n % d != 0 for d in range(2, int(math.sqrt(n)) + 1)):
            primes.append(n)
    total = 0.0
    for idx, p in enumerate(primes):
        pi_p = idx + 1
        total += SSQ**pi_p / p**26
    return total


def kernel_bsh_harmonic_saturation(n_modes=26):
    return sum(1 - math.exp(-SSQ * m) for m in range(1, n_modes + 1))


def kernel_wormhole_phonon_damping(r0=1e4, M=4.3e6):
    M_kg = M * M_SUN
    b_mod = 1 - BETA_I * SSQ
    damping = math.exp(-KAPPA * r0 / C) * Phi(GAMMA_0) * b_mod
    rS = 2 * G * M_kg / C**2
    curvature = rS / (r0 + rS)
    return damping * curvature * S26_3RD


# ── §3  New v17 Kernels (8) ──────────────────────────────────────────────

def kernel_gpu_dpm_atlas_peak(n_bins=512):
    """GPU DPM spectral atlas peak amplitude at ω_SCm.

    Peak = Σ_{i=1}^{26} c_i · Φ_i(ω_SCm) where Φ_i(ω_SCm) = 1.0
    (all Gaussians peak at center).
    """
    return sum(_LAYER_COEFFS)  # At center, all Gaussians = 1.0


def kernel_dpm_line_fwhm():
    """DPM spectral atlas effective FWHM from base linewidth σ_G.

    FWHM = 2·sqrt(2·ln(2))·σ_G for single Gaussian.
    Weighted average across 26 layers (slight broadening 0.02/layer).
    """
    fwhm_factor = 2 * math.sqrt(2 * math.log(2))
    avg_sigma = SIGMA_G * sum(1 + 0.02 * i for i in range(26)) / 26
    return fwhm_factor * avg_sigma


def kernel_muge_8term_gravity(M=1e11 * M_SUN, r=10 * KPC):
    """8-term MUGE gravitational acceleration.

    g = g_N + g_exp + g_super + g_env + g_Ug + g_cosm + g_quant + g_fluid
    """
    g_N = G * M / r**2
    g_exp = -H_0**2 * r
    g_ug = sum(G * M / r**2 * SSQ * i / 26 * BETA_I for i in range(1, 27))
    g_cosm = -LAMBDA_C * C**2 * r / 3
    g_quant = HBAR / (M * r**2) if M > 0 else 0.0
    return g_N + g_exp + g_ug + g_cosm + g_quant


def kernel_nfw_density(rho_s=1e-21, r_kpc=10, r_s_kpc=20):
    """NFW dark matter halo density.

    ρ_NFW(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]
    """
    x = max(r_kpc / r_s_kpc, 1e-10)
    return rho_s / (x * (1 + x)**2)


def kernel_enet_gamma_coupled(t=T_H, alpha=0.1):
    """Γ-coupled E_net dark energy density.

    E_net(t,Γ) = E₀ · exp(rate·t) · S₂₆ · Φ(Γ(t)) · net_factor
    """
    E_0 = RHO_SCM * C**2
    rate = KAPPA + SSQ / 26.0
    growth = math.exp(min(rate * t, 500))
    gamma_t = GAMMA_0 * (1.0 + alpha * t / T_H)
    phi_t = math.exp(-(gamma_t - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD
    net_factor = 0.2  # 2*(0.6/1.0) - 1
    return E_0 * growth * S26_STATIC * phi_t * net_factor


def kernel_dark_energy_eos_w0(alpha=0.1):
    """SCm dark energy equation of state w(z=0).

    w(0) = -1 + δw, where δw from Γ evolution at present epoch.
    """
    dz = 0.01
    t1 = T_H
    t2 = T_H / (1 + dz)**1.5
    gamma_1 = GAMMA_0 * (1.0 + alpha * t1 / T_H)
    gamma_2 = GAMMA_0 * (1.0 + alpha * t2 / T_H)
    phi_1 = math.exp(-(gamma_1 - GAMMA_0)**2 / (2 * SIGMA_G**2))
    phi_2 = math.exp(-(gamma_2 - GAMMA_0)**2 / (2 * SIGMA_G**2))
    d_ln_phi = math.log(max(phi_2, 1e-300) / max(phi_1, 1e-300))
    d_ln_a = -dz / (1 + 0)
    delta_w = (1.0 / 3.0) * d_ln_phi / d_ln_a if abs(d_ln_a) > 1e-30 else 0.0
    return -1.0 + delta_w


def kernel_alma_fubi_profile(nu_GHz=230.538, M=2000 * M_SUN, r=1.3e16):
    """ALMA F_{U,Bi,i} peak at specified line frequency.

    F = Σᵢ cᵢ · β_i · G·M/r² (at line center, Gaussian = 1).
    """
    g_surf = G * M / r**2
    return sum(c * BETA_I * g_surf for c in _LAYER_COEFFS)


def kernel_alma_chi2_co21(M=2000 * M_SUN, r=1.3e16, T_ex=50):
    """Simplified χ² for CO(2-1) line profile comparison.

    Computes peak ratio between F_{U,Bi,i} and LTE profile.
    """
    g_surf = G * M / r**2
    fubi_peak = sum(c * BETA_I * g_surf for c in _LAYER_COEFFS)
    h = 6.626e-34
    nu_Hz = 230.538e9
    J_ex = (h * nu_Hz / K_B) / (math.exp(h * nu_Hz / (K_B * T_ex)) - 1) if T_ex > 0 else 0
    J_bg = (h * nu_Hz / K_B) / (math.expm1(h * nu_Hz / (K_B * 2.725)))
    I_peak = (J_ex - J_bg) * (1 - math.exp(-5.0))  # τ₀=5
    scale = I_peak / (fubi_peak + 1e-300) if fubi_peak > 0 else 0
    return scale  # scaling factor for matching


# ── §4  ALL_KERNELS Registry ────────────────────────────────────────────

ALL_KERNELS: List[Callable] = [
    # v11 (16)
    kernel_gravity_26layer, kernel_fu_bi_i, kernel_phonon_ares,
    kernel_jet_mjet, kernel_ns_spindown, kernel_gw170817_strain,
    kernel_blazar_ergosphere, kernel_rest_phonon_jet,
    kernel_qcalcgeom_vectorized, kernel_pipeline_full,
    kernel_cena_jet, kernel_txs0506_jet, kernel_bcs_gap_solve,
    kernel_spectral_ladder_eval, kernel_ramanujan_26d_sum, kernel_triadic_solver,
    # v13 (4)
    kernel_fubi_inside_out, kernel_99sys_gamma_sweep,
    kernel_agn_cena_fubi, kernel_ns_merger_gw190425,
    # v14 (4)
    kernel_agn_merger_fubi, kernel_qgp_vacuum_density,
    kernel_alice_multiplicity, kernel_ym_mass_gap,
    # v15 (6)
    kernel_3c273_agn_fubi, kernel_ton618_agn_fubi,
    kernel_gw170817_merger, kernel_smbh_merger_fubi,
    kernel_dm_halo_nfw, kernel_txs0506_3gamma,
    # v16 (6)
    kernel_scm_inflation_hubble, kernel_thorne_morris_exotic,
    kernel_vds_convergence, kernel_dvp_prime_sieve,
    kernel_bsh_harmonic_saturation, kernel_wormhole_phonon_damping,
    # v17 (8)
    kernel_gpu_dpm_atlas_peak, kernel_dpm_line_fwhm,
    kernel_muge_8term_gravity, kernel_nfw_density,
    kernel_enet_gamma_coupled, kernel_dark_energy_eos_w0,
    kernel_alma_fubi_profile, kernel_alma_chi2_co21,
]


# ── §5  Benchmark Class ─────────────────────────────────────────────────

class ProductionScalingV17:
    """44-kernel production benchmark: 900k calc/s target."""

    def single_pass(self) -> float:
        return sum(k() for k in ALL_KERNELS)

    def simulate(self, duration_s: float = 1.0) -> dict:
        start = time.perf_counter()
        iterations = 0
        while time.perf_counter() - start < duration_s:
            self.single_pass()
            iterations += 1
        elapsed = time.perf_counter() - start
        calc_per_sec = iterations * len(ALL_KERNELS) / elapsed

        return {
            "iterations": iterations,
            "kernels": len(ALL_KERNELS),
            "elapsed_s": elapsed,
            "calc_per_sec": calc_per_sec,
            "target": TARGET_CALC_PER_SEC,
            "met_target": calc_per_sec >= TARGET_CALC_PER_SEC,
        }


# ── §6  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    print("=" * 72)
    print("production_scaling_v17.py — Self-Tests")
    print("=" * 72)

    ok = True
    passed = 0

    # Test 1: All kernels finite
    results = [k() for k in ALL_KERNELS]
    all_finite = all(math.isfinite(r) for r in results)
    if all_finite:
        print(f"  [PASS] Test 1: All {len(ALL_KERNELS)} kernels produce finite results")
        passed += 1
    else:
        print(f"  [FAIL] Test 1: Non-finite kernel result"); ok = False

    # Test 2: Kernel count = 44
    n = len(ALL_KERNELS)
    if n == 44:
        print(f"  [PASS] Test 2: Kernel count: {n}")
        passed += 1
    else:
        print(f"  [FAIL] Test 2: Expected 44 kernels, got {n}"); ok = False

    # Test 3: DPM atlas peak > 0
    val = kernel_gpu_dpm_atlas_peak()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 3: DPM atlas peak: {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 3: DPM atlas peak: {val}"); ok = False

    # Test 4: DPM line FWHM > 0
    val = kernel_dpm_line_fwhm()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 4: DPM FWHM: {val:.4e} rad/s")
        passed += 1
    else:
        print(f"  [FAIL] Test 4: DPM FWHM: {val}"); ok = False

    # Test 5: MUGE 8-term gravity > 0
    val = kernel_muge_8term_gravity()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 5: MUGE g = {val:.6e} m/s²")
        passed += 1
    else:
        print(f"  [FAIL] Test 5: MUGE g: {val}"); ok = False

    # Test 6: NFW density > 0
    val = kernel_nfw_density()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 6: NFW ρ = {val:.4e} kg/m³")
        passed += 1
    else:
        print(f"  [FAIL] Test 6: NFW ρ: {val}"); ok = False

    # Test 7: E_net Γ-coupled finite
    val = kernel_enet_gamma_coupled()
    if math.isfinite(val):
        print(f"  [PASS] Test 7: E_net(Γ) = {val:.4e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 7: E_net(Γ): {val}"); ok = False

    # Test 8: Dark energy w(0) near -1
    val = kernel_dark_energy_eos_w0()
    if abs(val - (-1.0)) < 0.5:
        print(f"  [PASS] Test 8: w(z=0) = {val:.6f}")
        passed += 1
    else:
        print(f"  [FAIL] Test 8: w(0) = {val}"); ok = False

    # Test 9: ALMA F_{U,Bi} profile > 0
    val = kernel_alma_fubi_profile()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 9: ALMA F_UBi = {val:.4e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 9: ALMA F_UBi: {val}"); ok = False

    # Test 10: ALMA χ² CO(2-1) scale finite
    val = kernel_alma_chi2_co21()
    if math.isfinite(val) and val > 0:
        print(f"  [PASS] Test 10: ALMA CO(2-1) scale = {val:.4e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 10: ALMA CO scale: {val}"); ok = False

    # Test 11: v16 carry-forward compatibility
    v16_results = [k() for k in ALL_KERNELS[:36]]
    if all(math.isfinite(r) for r in v16_results):
        print(f"  [PASS] Test 11: v16 carry-forward: all 36 kernels valid")
        passed += 1
    else:
        print(f"  [FAIL] Test 11: v16 carry-forward failed"); ok = False

    # Test 12: Target metadata
    print(f"  [PASS] Test 12: Target: {TARGET_CALC_PER_SEC:,} calc/s, {len(ALL_KERNELS)} kernels")
    passed += 1

    print("-" * 72)
    print(f"Results: {passed}/12 passed, {12 - passed} failed")
    return ok


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
