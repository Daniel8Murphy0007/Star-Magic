#!/usr/bin/env python3
"""
production_scaling_v18.py — Production Scaling to 950k calc/s

Session 225 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
v18 upgrade: 950,000 calc/s target (up from v17's 900k).
48 kernels total = v17's 44 + 4 new Session 225 kernels:
  kernel_qcalcgeom_master        — QCalcGeom master equation (r_cross × compact × S26 × Phi)
  kernel_cluster_cooling_suppression — Galaxy cluster cooling-flow suppression S(Gamma)
  kernel_ramanujan_binomial_Rn   — Explicit double-sum binomial R_n^{(26,3)}
  kernel_lenr_cop_gamma          — LENR COP as f(Gamma)

Architecture: Standalone kernel functions, ALL_KERNELS list, benchmark class.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from dpm_helpers import dpm_emergent_ug1, dpm_emergent_ug2

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
SIGMA_T   = 6.652e-29
M_PROTON  = 1.673e-27
ETA_BSFG  = 1e-22

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))

TARGET_CALC_PER_SEC = 950_000


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


_LAYER_COEFFS = [
    (SSQ ** i) / (i ** 26) * ramanujan_Rn(i, 3) for i in range(1, 27)
]


# ── §2  Carry-Forward Kernels (v11: 16, v13: 4, v14: 4, v15: 6, v16: 6, v17: 8 = 44) ──

def kernel_gravity_26layer(M=M_SUN, r=6.96e8):
    return sum(dpm_emergent_ug1(M, r) * SSQ * i / 26 for i in range(1, 27))

def kernel_fu_bi_i(M=M_SUN, r=6.96e8, t=86400):
    Ug = sum(dpm_emergent_ug1(M, r) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_emergent_ug1(M, r) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    return (Ug - Ub) + F_NEUTRON * S26_3RD * Phi(GAMMA_0) * math.exp(min(KAPPA * t, 500))

def kernel_phonon_ares(M=M_SUN, r=6.96e8):
    return sum(math.exp(-SSQ * i / 26) * BETA_I * dpm_emergent_ug1(M, r) for i in range(1, 27))

def kernel_jet_mjet(A=1.5, gamma_THz=0.1):
    gamma = 2 * PI * gamma_THz * 1e12
    return 1 + A * math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2))

def kernel_ns_spindown(P=0.033, Pdot=4.2e-13):
    return 3.2e19 * math.sqrt(P * Pdot)

def kernel_gw170817_strain(M_chirp=1.186, d_Mpc=40):
    f = 100; d = d_Mpc * 3.086e22
    return 4 / d * (G * M_chirp * M_SUN / C**2)**(5/3) * (PI * f / C)**(2/3)

def kernel_blazar_ergosphere(M=3e8, a=0.95):
    M_kg = M * M_SUN; rS = 2 * dpm_emergent_ug1(M_kg, C)
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
    Um = dpm_emergent_ug1(M, r) * SSQ * 0.1
    return Ug + Um - Ub

# v13 carry-forward (4)
def kernel_fubi_inside_out(M=5e7, r=None):
    M_kg = M * M_SUN
    if r is None: r = 2 * dpm_emergent_ug1(M_kg, C)
    return kernel_fu_bi_i(M_kg, r, 86400)

def kernel_99sys_gamma_sweep():
    total = 0.0
    for g in [0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0]:
        total += kernel_jet_mjet(1.5, g)
    return total

def kernel_agn_cena_fubi(M=5.5e7, a=0.70):
    M_kg = M * M_SUN; rS = 2 * dpm_emergent_ug1(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    return kernel_fu_bi_i(M_kg, rH, 86400)

def kernel_ns_merger_gw190425(M_chirp=1.44, d_Mpc=159):
    return kernel_gw170817_strain(M_chirp, d_Mpc)

# v14 carry-forward (4)
def kernel_agn_merger_fubi(M1=5.5e7, M2=3e7, a=0.70):
    M_kg = (M1 + M2) * M_SUN; rS = 2 * dpm_emergent_ug1(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(dpm_emergent_ug1(M_kg, rH) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_emergent_ug1(M_kg, rH) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
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
    M_kg = M * M_SUN; rS = 2 * dpm_emergent_ug1(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(dpm_emergent_ug1(M_kg, rH) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_emergent_ug1(M_kg, rH) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    ratio = abs(Ub) / (abs(Ug) + abs(Ub) + 1e-300)
    return ratio * (1 + A_jet) * S26_3RD

def kernel_ton618_agn_fubi(M=6.6e10, a=0.998, A_jet=2.8):
    M_kg = M * M_SUN; rS = 2 * dpm_emergent_ug1(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(dpm_emergent_ug1(M_kg, rH) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_emergent_ug1(M_kg, rH) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    ratio = abs(Ub) / (abs(Ug) + abs(Ub) + 1e-300)
    return ratio * (1 + A_jet) * S26_3RD

def kernel_gw170817_merger(M_chirp=1.186, d_Mpc=40, suppression=0.667):
    h0 = kernel_gw170817_strain(M_chirp, d_Mpc)
    return h0 * (1 - suppression * Phi(GAMMA_0) / S26_3RD)

def kernel_smbh_merger_fubi(M1=5.5e7, M2=3e7, a1=0.70, a2=0.60):
    M_total = (M1 + M2) * M_SUN
    rS = 2 * dpm_emergent_ug1(M_total, C)
    a_eff = (M1 * a1 + M2 * a2) / (M1 + M2)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a_eff**2, 0)))
    rho = RHO_VAC * math.exp(min(KAPPA * 86400, 500))
    V = (4/3) * PI * rH**3
    return rho * V * S26_3RD**2 * Phi(GAMMA_0)

def kernel_dm_halo_nfw(rho_0=5e-22, r_kpc=20, r_s_kpc=20):
    x = max(r_kpc / r_s_kpc, 1e-10)
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
    r_s = 2 * dpm_emergent_ug1(M_kg, C)
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
    rS = 2 * dpm_emergent_ug1(M_kg, C)
    curvature = rS / (r0 + rS)
    return damping * curvature * S26_3RD

# v17 carry-forward (8)
def kernel_gpu_dpm_atlas_peak(n_bins=512):
    return sum(_LAYER_COEFFS)

def kernel_dpm_line_fwhm():
    fwhm_factor = 2 * math.sqrt(2 * math.log(2))
    avg_sigma = SIGMA_G * sum(1 + 0.02 * i for i in range(26)) / 26
    return fwhm_factor * avg_sigma

def kernel_muge_8term_gravity(M=1e11 * M_SUN, r=10 * KPC):
    g_N = dpm_emergent_ug1(M, r)
    g_exp = -H_0**2 * r
    g_ug = sum(dpm_emergent_ug1(M, r) * SSQ * i / 26 * BETA_I for i in range(1, 27))
    g_cosm = -LAMBDA_C * C**2 * r / 3
    g_quant = HBAR / (M * r**2) if M > 0 else 0.0
    return g_N + g_exp + g_ug + g_cosm + g_quant

def kernel_nfw_density(rho_s=1e-21, r_kpc=10, r_s_kpc=20):
    x = max(r_kpc / r_s_kpc, 1e-10)
    return rho_s / (x * (1 + x)**2)

def kernel_enet_gamma_coupled(t=T_H, alpha=0.1):
    E_0 = RHO_SCM * C**2
    rate = KAPPA + SSQ / 26.0
    growth = math.exp(min(rate * t, 500))
    gamma_t = GAMMA_0 * (1.0 + alpha * t / T_H)
    phi_t = math.exp(-(gamma_t - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD
    net_factor = 0.2
    return E_0 * growth * S26_STATIC * phi_t * net_factor

def kernel_dark_energy_eos_w0(alpha=0.1):
    dz = 0.01
    t1 = T_H; t2 = T_H / (1 + dz)**1.5
    gamma_1 = GAMMA_0 * (1.0 + alpha * t1 / T_H)
    gamma_2 = GAMMA_0 * (1.0 + alpha * t2 / T_H)
    phi_1 = math.exp(-(gamma_1 - GAMMA_0)**2 / (2 * SIGMA_G**2))
    phi_2 = math.exp(-(gamma_2 - GAMMA_0)**2 / (2 * SIGMA_G**2))
    d_ln_phi = math.log(max(phi_2, 1e-300) / max(phi_1, 1e-300))
    d_ln_a = -dz / (1 + 0)
    delta_w = (1.0 / 3.0) * d_ln_phi / d_ln_a if abs(d_ln_a) > 1e-30 else 0.0
    return -1.0 + delta_w

def kernel_alma_fubi_profile(nu_GHz=230.538, M=2000 * M_SUN, r=1.3e16):
    g_surf = dpm_emergent_ug1(M, r)
    return sum(c * BETA_I * g_surf for c in _LAYER_COEFFS)

def kernel_alma_chi2_co21(M=2000 * M_SUN, r=1.3e16, T_ex=50):
    g_surf = dpm_emergent_ug1(M, r)
    fubi_peak = sum(c * BETA_I * g_surf for c in _LAYER_COEFFS)
    h = 6.626e-34; nu_Hz = 230.538e9
    J_ex = (h * nu_Hz / K_B) / (math.exp(h * nu_Hz / (K_B * T_ex)) - 1) if T_ex > 0 else 0
    J_bg = (h * nu_Hz / K_B) / (math.expm1(h * nu_Hz / (K_B * 2.725)))
    I_peak = (J_ex - J_bg) * (1 - math.exp(-5.0))
    scale = I_peak / (fubi_peak + 1e-300) if fubi_peak > 0 else 0
    return scale


# ── §3  New v18 Kernels (4) ──────────────────────────────────────────────

def kernel_qcalcgeom_master(M=M_SUN, gamma=GAMMA_0):
    """QCalcGeom master equation: r_cross * (26!)^{-1/13} * S26^{(3)} * Phi(Gamma).

    Combines BSFG crossover radius, 26D compactification scale,
    Ramanujan polylogarithmic sum, and phonon fluence.
    """
    r_s = 2 * dpm_emergent_ug1(M, C)
    r_cross = math.sqrt(ETA_BSFG) * r_s
    compact_scale = math.factorial(26) ** (-1.0 / 13.0)
    phi_val = Phi(gamma)
    return r_cross * compact_scale * S26_3RD * phi_val


def kernel_cluster_cooling_suppression(T_keV=4.0, n_e=3e4, r_cool_kpc=100,
                                        M_bh=3.4e8, eta_jet=0.1):
    """Galaxy cluster cooling-flow suppression factor S(Gamma).

    S = min(Q_heat / L_cool, 1) where:
      L_cool = n_e^2 * Lambda(T) * V_cool
      Q_heat ~ eta_jet * L_Edd (Eddington-limited heating)
    """
    T_K = T_keV * 1.1605e7
    g_ff = 1.2
    Lambda_ff = 1.42e-40 * math.sqrt(T_K) * g_ff * (1 + 1.8 * 0.3)
    r_cool = r_cool_kpc * KPC
    V_cool = (4.0 / 3.0) * PI * r_cool**3
    L_cool = n_e**2 * Lambda_ff * V_cool

    M_bh_kg = M_bh * M_SUN
    L_edd = 4 * PI * G * M_bh_kg * M_PROTON * C / SIGMA_T
    Q_heat = eta_jet * L_edd

    return min(Q_heat / L_cool, 1.0) if L_cool > 0 else 0.0


def kernel_ramanujan_binomial_Rn(n=5, D=26, k=3):
    """Explicit double-sum binomial expansion R_n^{(D,k)}.

    R_n = (2*pi)^{n/6} / n! * [1 + sum_m sum_j (-1)^{j+1} C(D,j)(D-j)!/n^{j+Dm}]
    """
    prefactor = (2 * PI) ** (n / 6.0) / math.factorial(min(n, 170))
    correction = 1.0
    for m in range(1, k + 1):
        inner = 0.0
        for j in range(1, D + 1):
            sign = (-1) ** (j + 1)
            comb = math.comb(D, j)
            fact_dj = math.factorial(D - j)
            inner += sign * comb * fact_dj / (n ** j)
        correction += inner / (n ** (D * m))
    return prefactor * correction


def kernel_lenr_cop_gamma(gamma=GAMMA_0, T_cell=350.0):
    """LENR COP as function of Gamma.

    COP = P_out / P_in where:
      P_in = hbar * omega_SCm * phi_0 * Phi(Gamma) * A_cell
      P_out = E_dd * R_nd * V_active
    """
    phi_0 = 1e20; A_cell = 1e-4; V_active = 1e-6
    sigma_n = 1e-4; n_H = 6e28
    E_dd = 23.8e6 * 1.602e-19
    E_a = 0.3 * 1.602e-19

    phi_val = Phi(gamma)
    P_in = HBAR * OMEGA_SCM * phi_0 * phi_val * A_cell

    boltzmann = math.exp(-E_a / (K_B * T_cell))
    P_plasmoid = phi_val / S26_3RD if S26_3RD > 0 else 0
    R_nd = sigma_n * n_H * phi_0 * phi_val * boltzmann * P_plasmoid
    P_out = E_dd * R_nd * V_active

    return P_out / P_in if P_in > 0 else 0.0


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
    # v18 (4)
    kernel_qcalcgeom_master, kernel_cluster_cooling_suppression,
    kernel_ramanujan_binomial_Rn, kernel_lenr_cop_gamma,
]


# ── §5  Benchmark Class ─────────────────────────────────────────────────

class ProductionScalingV18:
    """48-kernel production benchmark: 950k calc/s target."""

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
    print("production_scaling_v18.py — Self-Tests")
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

    # Test 2: Kernel count = 48
    n = len(ALL_KERNELS)
    if n == 48:
        print(f"  [PASS] Test 2: Kernel count: {n}")
        passed += 1
    else:
        print(f"  [FAIL] Test 2: Expected 48 kernels, got {n}"); ok = False

    # Test 3: QCalcGeom master > 0
    val = kernel_qcalcgeom_master()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 3: QCalcGeom master = {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 3: QCalcGeom = {val}"); ok = False

    # Test 4: Cluster cooling suppression in [0, 1]
    val = kernel_cluster_cooling_suppression()
    if 0 <= val <= 1:
        print(f"  [PASS] Test 4: Cooling suppression S = {val:.6f}")
        passed += 1
    else:
        print(f"  [FAIL] Test 4: S = {val}"); ok = False

    # Test 5: Ramanujan binomial finite
    val = kernel_ramanujan_binomial_Rn()
    if math.isfinite(val):
        print(f"  [PASS] Test 5: R_5^(26,3) = {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 5: R_n = {val}"); ok = False

    # Test 6: LENR COP > 0
    val = kernel_lenr_cop_gamma()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 6: LENR COP = {val:.4e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 6: COP = {val}"); ok = False

    # Test 7: v17 carry-forward kernels still work
    v17_vals = [kernel_gpu_dpm_atlas_peak(), kernel_dpm_line_fwhm(),
                kernel_muge_8term_gravity(), kernel_nfw_density()]
    if all(math.isfinite(v) for v in v17_vals):
        print(f"  [PASS] Test 7: v17 carry-forward kernels OK")
        passed += 1
    else:
        print(f"  [FAIL] Test 7: v17 carry-forward failure"); ok = False

    # Test 8: Single pass produces finite total
    bench = ProductionScalingV18()
    total = bench.single_pass()
    if math.isfinite(total):
        print(f"  [PASS] Test 8: Single pass total = {total:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 8: Non-finite pass total"); ok = False

    # Test 9: Benchmark runs without error
    result = bench.simulate(0.5)
    if result['calc_per_sec'] > 0:
        label = "MET TARGET" if result['met_target'] else "below target"
        print(f"  [PASS] Test 9: {result['calc_per_sec']:.0f} calc/s ({label})")
        passed += 1
    else:
        print(f"  [FAIL] Test 9: Benchmark failed"); ok = False

    # Test 10: Kernel names are unique
    names = [k.__name__ for k in ALL_KERNELS]
    if len(names) == len(set(names)):
        print(f"  [PASS] Test 10: All {len(names)} kernel names unique")
        passed += 1
    else:
        print(f"  [FAIL] Test 10: Duplicate kernel names"); ok = False

    # Test 11: LENR COP peaks at Gamma_0
    cop_peak = kernel_lenr_cop_gamma(gamma=GAMMA_0)
    cop_off = kernel_lenr_cop_gamma(gamma=2 * PI * 0.50 * 1e12)
    if cop_peak > cop_off:
        print(f"  [PASS] Test 11: COP peaks at Gamma_0 ({cop_peak:.2e} > {cop_off:.2e})")
        passed += 1
    else:
        print(f"  [FAIL] Test 11: COP peak location wrong"); ok = False

    # Test 12: Cluster suppression increases with larger BH
    s1 = kernel_cluster_cooling_suppression(M_bh=3.4e8)
    s2 = kernel_cluster_cooling_suppression(M_bh=6.5e9)
    if s2 >= s1:
        print(f"  [PASS] Test 12: Suppression increases with M_BH ({s1:.4f} -> {s2:.4f})")
        passed += 1
    else:
        print(f"  [FAIL] Test 12: S not monotone in M_BH"); ok = False

    print(f"\n  production_scaling_v18.py: {passed}/12 tests passed")
    return ok


if __name__ == "__main__":
    success = _run_tests()
    raise SystemExit(0 if success else 1)
