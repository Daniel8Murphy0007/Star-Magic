#!/usr/bin/env python3
"""
production_scaling_v27.py — Production Scaling to 1.1M calc/s

Session 226 | Daniel Murphy
────────────────────────────────────────────────────────────────────────────────
v27 upgrade from v26 pattern: 1,100,000 calc/s target.
64 kernels total = v26's 56 + 8 new Session 226-B kernels:
  kernel_exp_gate_fidelity       — Exponential F_gate = exp(-Γt/T₂)·S₂₆⁽³⁾
  kernel_pimath_key_gen          — PImath K = S₂₆⁽³⁾·π^(n/26) mod 113
  kernel_muge_3d_gravity         — 3D MUGE 26-layer gravity field
  kernel_sigma_n_density         — σ_n(ω,ρ) cross section (2-arg)
  kernel_bsfg_geodesic_trace     — BSFG wormhole geodesic throughput
  kernel_triadic_master_gpu      — Triadic master gravity function
  kernel_phi_qubit_lagrangian    — φ_qubit Lagrangian sector
  kernel_appendix_template       — Whitepaper appendix template engine

Architecture: Standalone kernel functions, ALL_KERNELS list, benchmark class.
Endpoint: /api/phonon/jet for REST integration via uqff_server.js
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
eV        = 1.602e-19

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
F_UBI_RATIO = 0.6
H_SCM     = 0.99

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))

TARGET_CALC_PER_SEC = 1_100_000


# ── §1  Ramanujan S₂₆⁽³⁾ ──────────────────────────────────────────────────

def ramanujan_Rn(n: int, k: int = 3) -> float:
    prefactor = (2 * PI) ** (n / 6.0) / math.factorial(min(n, 170))
    correction = 0.0
    for m in range(1, k + 1):
        inner = 0.0
        for j in range(1, 27):
            sign = (-1) ** (j + 1)
            binom = math.comb(26, j)
            inner += sign * binom * math.factorial(26 - j) / n ** j
        correction += inner / n ** (26 * m)
    return prefactor * (1.0 + correction)


S26_3RD = sum((SSQ ** n) / (n ** 26) * ramanujan_Rn(n, 3) for n in range(1, 28))


def Phi(gamma: float) -> float:
    return math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD


_LAYER_COEFFS = [
    (SSQ ** i) / (i ** 26) * ramanujan_Rn(i, 3) for i in range(1, 27)
]


# ── §2  DPM-FOUNDATION GRAVITY helpers ───────────────────────────────────────

def _dpm_ug1(r: float, mu: float = 1e15, R_body: float = 1e7) -> float:
    if r <= 0:
        return 0.0
    return mu / (4 * PI * r ** 3) * S26_3RD


def _dpm_ug2(r: float, Z: float = 26, alpha: float = 1.0) -> float:
    if r <= 0:
        return 0.0
    return (Z * alpha / r ** 2) * SSQ * BETA_I


def _dpm_ug3(r: float, t: float = 0.0, omega: float = OMEGA_SCM) -> float:
    if r <= 0:
        return 0.0
    return (HBAR * omega / r) * math.cos(omega * t) * S26_3RD


def _dpm_ug4(r: float, rho_vac: float = RHO_VAC) -> float:
    if r <= 0:
        return 0.0
    return (KAPPA * rho_vac / r) * SSQ


def _dpm_full(r: float = 1e10, M: float = 4.3e6 * M_SUN) -> float:
    return _dpm_ug1(r) + _dpm_ug2(r) + _dpm_ug3(r) + _dpm_ug4(r)


# ── §3  v26 Kernel Functions (inherited) ───────────────────────────────────

def kernel_gravity_26layer(r: float = 1e10) -> float:
    total = 0.0
    for i in range(1, 27):
        ri = r * (1 + 0.01 * i)
        total += _LAYER_COEFFS[i - 1] * (_dpm_ug1(ri) + _dpm_ug2(ri) + _dpm_ug3(ri) + _dpm_ug4(ri))
    return total

def kernel_fu_bi_i(r: float = 1e10, M: float = 4.3e6 * M_SUN) -> float:
    Ug1 = _dpm_ug1(r)
    Ubi = BETA_I * Ug1 * (M / r)
    return Ubi * S26_3RD * F_UBI_RATIO

def kernel_phonon_ares(gamma: float = GAMMA_0) -> float:
    return Phi(gamma) * H_SCM * BETA_I

def kernel_jet_mjet(gamma: float = GAMMA_0) -> float:
    A_jet = 1.5
    dg = gamma - GAMMA_0
    return 1 + A_jet * math.exp(-dg**2 / (2 * SIGMA_G**2))

def kernel_ns_spindown() -> float:
    omega_dot = -4.2e-15
    return omega_dot * (1 + Phi(GAMMA_0) * SSQ / 26)

def kernel_gw170817_strain() -> float:
    return 5.4176e-22 * 0.333

def kernel_blazar_ergosphere(M: float = 6.5e9 * M_SUN, a: float = 0.9) -> float:
    r_ergo = 2 * G * M / C ** 2
    return r_ergo * (1 + a * SSQ * BETA_I)

def kernel_rest_phonon_jet() -> float:
    P_BZ = (PI / 6) * (G * 6.5e9 * M_SUN / C ** 2) ** 2 * 0.9 ** 2 * 50 ** 2 * C
    return P_BZ * kernel_jet_mjet()

def kernel_qcalcgeom_vectorized() -> float:
    return sum(_LAYER_COEFFS[i] * (SSQ ** i) for i in range(26))

def kernel_pipeline_full() -> float:
    return kernel_gravity_26layer() + kernel_fu_bi_i() + kernel_phonon_ares()

def kernel_cena_jet() -> float:
    M = 5.5e7 * M_SUN
    r_s = 2 * G * M / C ** 2
    P = (50 ** 2 / (8 * PI)) * (r_s / C) ** 2 * 0.7 ** 2 * C
    return P * (1 + Phi(GAMMA_0))

def kernel_txs0506_jet() -> float:
    M = 3e8 * M_SUN
    r_s = 2 * G * M / C ** 2
    P = (8000 ** 2 / (8 * PI)) * (r_s / C) ** 2 * 0.85 ** 2 * C
    return P * Phi(GAMMA_0)

def kernel_bcs_gap_solve(T: float = 0.015) -> float:
    delta_scm = HBAR * OMEGA_SCM
    ratio = delta_scm / (K_B * max(T, 1e-30))
    ratio = min(ratio, 700)
    return delta_scm * math.tanh(1.74 * math.sqrt(ratio - 1)) if ratio > 1 else delta_scm

def kernel_spectral_ladder_eval() -> float:
    return sum(c * math.exp(-i * 0.01) for i, c in enumerate(_LAYER_COEFFS))

def kernel_ramanujan_26d_sum() -> float:
    return S26_3RD

def kernel_triadic_solver() -> float:
    return _dpm_ug1(1e10) + _dpm_ug2(1e10) + _dpm_ug3(1e10) + _dpm_ug4(1e10)

def kernel_fubi_inside_out() -> float:
    total = 0.0
    for i in range(1, 27):
        r_i = 1e10 * (26 - i + 1) / 26
        total += _dpm_full(max(r_i, 1e6)) * _LAYER_COEFFS[i - 1]
    return total * F_UBI_RATIO

def kernel_99sys_gamma_sweep() -> float:
    total = 0.0
    for g in range(1, 100):
        gamma = 2 * PI * g * 0.01 * 1e12
        total += Phi(gamma)
    return total / 99

def kernel_agn_cena_fubi() -> float:
    return kernel_fu_bi_i(r=1e11, M=5.5e7 * M_SUN)

def kernel_ns_merger_gw190425() -> float:
    return 3.0e-22 * 0.530

def kernel_agn_merger_fubi() -> float:
    return kernel_fu_bi_i(r=1e12, M=1e8 * M_SUN)

def kernel_qgp_vacuum_density() -> float:
    return RHO_SCM * C ** 2 * SSQ * S26_3RD

def kernel_alice_multiplicity() -> float:
    s_NN = 5020  # GeV
    return 2.0 * math.log(s_NN) ** 2 * SSQ * BETA_I

def kernel_ym_mass_gap() -> float:
    Lambda_QCD = 0.2  # GeV
    return Lambda_QCD * (1 + SSQ * S26_3RD)

def kernel_3c273_agn_fubi() -> float:
    return kernel_fu_bi_i(r=1e13, M=8.86e8 * M_SUN)

def kernel_ton618_agn_fubi() -> float:
    return kernel_fu_bi_i(r=1e14, M=6.6e10 * M_SUN)

def kernel_gw170817_merger() -> float:
    return 5.4176e-22 * 0.333 * math.exp(SSQ * 0 / 26)

def kernel_smbh_merger_fubi() -> float:
    return kernel_fu_bi_i(r=1e15, M=1e10 * M_SUN)

def kernel_dm_halo_nfw() -> float:
    rho_s = 1e7 * M_SUN / KPC ** 3
    r_s = 20 * KPC
    r = 8.5 * KPC
    x = r / r_s
    return rho_s / (x * (1 + x) ** 2)

def kernel_txs0506_3gamma() -> float:
    total = 0.0
    for g_THz in [0.05, 0.10, 0.30]:
        gamma = 2 * PI * g_THz * 1e12
        total += Phi(gamma)
    return total / 3

def kernel_scm_inflation_hubble() -> float:
    rho_inf = RHO_CRIT * 1e60
    return math.sqrt(8 * PI * G * rho_inf / 3) * SSQ * S26_3RD

def kernel_thorne_morris_exotic() -> float:
    r0 = 1e4
    b0 = 1e4
    b_mod = 1 - BETA_I * SSQ
    return -(1 / (8 * PI * G)) * (b0 * b_mod / r0 ** 2)

def kernel_vds_convergence() -> float:
    return sum((SSQ ** n) / (n ** 26) for n in range(1, 28))

def kernel_dvp_prime_sieve() -> float:
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23]
    return sum(SSQ ** p / p ** 26 for p in primes) * S26_3RD

def kernel_bsh_harmonic_saturation() -> float:
    return sum(math.sin(n * PI / 26) / n for n in range(1, 27)) * SSQ

def kernel_wormhole_phonon_damping() -> float:
    D_aether = 1.0
    D_scm = 1.0
    D_trz = 0.900
    D_string = 0.370
    return D_aether * D_scm * D_trz * D_string * Phi(GAMMA_0)

def kernel_gpu_dpm_atlas_peak() -> float:
    return sum(_LAYER_COEFFS) * math.exp(0)

def kernel_dpm_line_fwhm() -> float:
    return SIGMA_G * 2 * math.sqrt(2 * math.log(2))

def kernel_muge_8term_gravity(r: float = 1e10) -> float:
    a_dpm = _dpm_full(r)
    a_exp = H_0 * C * r
    a_super = 1e-15 * math.exp(-r / (1e15))
    a_envelope = SSQ * BETA_I * a_dpm
    a_ug_sum = sum(_dpm_full(r * (1 + 0.01 * i)) for i in range(8))
    a_cosm = LAMBDA_C * C ** 2 * r / 3
    a_quantum = HBAR / (M_PROTON * r ** 2)
    a_fluid = ETA_BSFG * C / r
    return a_dpm + a_exp + a_super + a_envelope + a_ug_sum + a_cosm + a_quantum + a_fluid

def kernel_nfw_density() -> float:
    return kernel_dm_halo_nfw()

def kernel_enet_gamma_coupled() -> float:
    E0 = 1.0
    t = 0.0
    return E0 * math.exp(KAPPA * t + SSQ * t / 26) * kernel_vds_convergence() * (2 * F_UBI_RATIO - 1)

def kernel_dark_energy_eos_w0() -> float:
    return -1 + 2 * (KAPPA + SSQ / 26) / (3 * H_0)

def kernel_alma_fubi_profile() -> float:
    return kernel_fu_bi_i(r=1e14, M=50 * M_SUN) * Phi(GAMMA_0)

def kernel_alma_chi2_co21() -> float:
    obs = 1.5e-20
    model = kernel_alma_fubi_profile()
    sigma = obs * 0.1
    return ((obs - model) / sigma) ** 2 if sigma > 0 else 0.0

def kernel_qcalcgeom_master() -> float:
    return kernel_qcalcgeom_vectorized() * kernel_gravity_26layer()

def kernel_cluster_cooling_suppression() -> float:
    T_cluster = 5e7  # K
    n_e = 1e3  # m^-3
    Lambda = 1e-23 * math.sqrt(T_cluster) * SSQ * S26_3RD
    return n_e ** 2 * Lambda

def kernel_ramanujan_binomial_Rn() -> float:
    return sum(ramanujan_Rn(n, 3) for n in range(1, 27))

def kernel_lenr_cop_gamma() -> float:
    sigma_0 = 1e-28
    dw = 0
    sigma = sigma_0 * math.exp(-dw ** 2 / (2 * GAMMA_0 ** 2)) * S26_3RD
    P_out = sigma * 1e20  # flux
    P_in = 100  # W
    return P_out / P_in if P_in > 0 else 0.0

def kernel_phonon_gate_unitary() -> float:
    theta = PI / 4
    phi_val = Phi(GAMMA_0)
    return abs(math.exp(-1j * theta / 2 * phi_val))

def kernel_dvp_vorticity_bound() -> float:
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
    enstrophy = sum(1 / p ** 2 for p in primes) * SSQ * BETA_I
    return enstrophy * S26_3RD

def kernel_bsm_mass_bound() -> float:
    delta_scm = HBAR * OMEGA_SCM
    return delta_scm / eV * SSQ * BETA_I

def kernel_pcore_composite() -> float:
    rho_Fe = 12800
    rho_Si = 4500
    f_Fe = 0.35
    f_Si = 0.65
    return (f_Fe * rho_Fe + f_Si * rho_Si) * (1 + SSQ * BETA_I * S26_3RD)

def kernel_cooper_pair_gap(T: float = 0.015) -> float:
    delta_0 = HBAR * OMEGA_SCM * SSQ * BETA_I
    return delta_0 / eV

def kernel_ramanujan_router() -> float:
    domains = [
        kernel_gravity_26layer, kernel_phonon_ares,
        kernel_bcs_gap_solve, kernel_ns_spindown,
        kernel_dm_halo_nfw, kernel_vds_convergence,
        kernel_dvp_prime_sieve, kernel_bsh_harmonic_saturation,
    ]
    total = sum(abs(d()) for d in domains)
    return total * S26_3RD

def kernel_cfd_stam_step() -> float:
    N = 8
    v = [[math.sin(PI * i / N) * math.cos(PI * j / N) for j in range(N)] for i in range(N)]
    max_v = max(abs(v[i][j]) for i in range(N) for j in range(N))
    return max_v * SSQ * BETA_I

def kernel_vds_dvp_bsh_cross() -> float:
    return kernel_vds_convergence() * kernel_dvp_prime_sieve() * kernel_bsh_harmonic_saturation()


# ── §4  NEW v27 Kernels ────────────────────────────────────────────────────

def kernel_exp_gate_fidelity(T: float = 0.015, t_gate: float = 1e-9) -> float:
    """Exponential gate fidelity: F = exp(-Γ·t/T₂)·S₂₆⁽³⁾·(1-F_{UBi}/F_U)."""
    delta_scm = HBAR * OMEGA_SCM
    ratio = delta_scm / (K_B * max(T, 1e-30))
    ratio = min(ratio, 700)
    T2 = (HBAR / delta_scm) * math.exp(ratio) * S26_3RD * F_UBI_RATIO
    exponent = max(-GAMMA_0 * t_gate / T2, -700)
    return math.exp(exponent) * S26_3RD * (1.0 - F_UBI_RATIO)


def kernel_pimath_key_gen() -> float:
    """PImath key: K = floor(S₂₆⁽³⁾·π^(n/26)·10¹²) mod 113, summed n=1..26."""
    total = 0
    for n in range(1, 27):
        raw = S26_3RD * (PI ** (n / 26.0))
        total += int(abs(raw) * 1e12) % 113
    return float(total)


def kernel_muge_3d_gravity(r: float = 1e10) -> float:
    """3D MUGE 26-layer gravity field magnitude."""
    total = 0.0
    for i in range(1, 27):
        ri = r * (1 + 0.01 * i)
        total += _LAYER_COEFFS[i - 1] * (_dpm_ug1(ri) + _dpm_ug2(ri) + _dpm_ug3(ri) + _dpm_ug4(ri))
    return abs(total)


def kernel_sigma_n_density(omega: float = OMEGA_SCM, rho: float = 1e4) -> float:
    """σ_n(ω,ρ) — cross section with frequency AND density dependence."""
    sigma_0 = 1e-28
    dw = omega - OMEGA_SCM
    gaussian = math.exp(-dw ** 2 / (2 * GAMMA_0 ** 2))
    density_factor = 1 + SSQ * rho / (rho + 1e6)
    return sigma_0 * gaussian * S26_3RD * density_factor


def kernel_bsfg_geodesic_trace() -> float:
    """BSFG wormhole geodesic: single null geodesic step."""
    M = 4.3e6 * M_SUN
    r = 5e4
    b0 = 1e4
    r0 = 1e4
    b = b0 * (r0 / r) * (1 - BETA_I * SSQ)
    Phi_r = -G * M / (r * C ** 2) * (1 + ETA_BSFG * 4.04)
    g_rr = 1.0 / (1 - b / r) if abs(1 - b / r) > 1e-30 else 1e30
    V_eff = (1 - b / r) * (1e10 ** 2 / (b ** 2 + r ** 2))
    return V_eff * math.exp(2 * Phi_r)


def kernel_triadic_master_gpu() -> float:
    """Triadic master: g(r) = Ug1+Ug2+Ug3+Ug4 with 26D sum."""
    r = 1e10
    total = 0.0
    for i in range(1, 27):
        omega_i = OMEGA_SCM * (1 + 0.01 * (i - 1))
        ri = r * i
        u1 = _dpm_ug1(ri)
        u2 = _dpm_ug2(ri)
        u3 = _dpm_ug3(ri, omega=omega_i)
        u4 = _dpm_ug4(ri)
        total += _LAYER_COEFFS[i - 1] * (u1 + u2 + u3 + u4)
    return total


def kernel_phi_qubit_lagrangian() -> float:
    """φ_qubit Lagrangian sector: L = (-β_i·Ug1·M/r·U_UA + F_n·Φ)·F_gate."""
    r = 1e4
    M = 1.989e30
    mu = 1e11 * r ** 3
    Ug1 = _dpm_ug1(r, mu=mu)
    buoyancy = -BETA_I * Ug1 * (M / r) * 1e-4
    phonon = F_NEUTRON * Phi(GAMMA_0)
    F_gate = kernel_exp_gate_fidelity()
    return (buoyancy + phonon) * F_gate


def kernel_appendix_template() -> float:
    """Whitepaper appendix section count template."""
    sections = [
        "SCm_Qubit_Gate_Fidelity", "PImath_Cryptography",
        "LENR_Cross_Section", "MUGE_3D_Systems",
    ]
    return float(len(sections)) * S26_3RD


# ── §5  ALL_KERNELS Registry ────────────────────────────────────────────

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
    # v26 (8)
    kernel_phonon_gate_unitary, kernel_dvp_vorticity_bound,
    kernel_bsm_mass_bound, kernel_pcore_composite,
    kernel_cooper_pair_gap, kernel_ramanujan_router,
    kernel_cfd_stam_step, kernel_vds_dvp_bsh_cross,
    # v27 (8) — NEW
    kernel_exp_gate_fidelity, kernel_pimath_key_gen,
    kernel_muge_3d_gravity, kernel_sigma_n_density,
    kernel_bsfg_geodesic_trace, kernel_triadic_master_gpu,
    kernel_phi_qubit_lagrangian, kernel_appendix_template,
]


# ── §6  Benchmark Class ─────────────────────────────────────────────────

class ProductionScalingV27:
    """64-kernel production benchmark: 1.1M calc/s target.

    REST endpoint: /api/phonon/jet via uqff_server.js integration.
    """

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

    def compute(self, dataset: dict) -> dict:
        duration = float(dataset.get("duration_s", 1.0))
        result = self.simulate(duration)
        result["primary_equations"] = [
            f"64 kernels × {result['iterations']} iterations in {result['elapsed_s']:.3f}s",
            f"Throughput: {result['calc_per_sec']:.0f} calc/s",
            f"Target: {TARGET_CALC_PER_SEC:,} calc/s",
            f"Met target: {result['met_target']}",
        ]
        return result


# ── §7  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    print("=" * 72)
    print("production_scaling_v27.py — Self-Tests")
    print("=" * 72)

    passed = 0
    ok = True

    # Test 1: all kernels produce finite results
    results = [k() for k in ALL_KERNELS]
    if all(math.isfinite(r) for r in results):
        print(f"  [PASS] Test 1: All {len(ALL_KERNELS)} kernels produce finite results")
        passed += 1
    else:
        print(f"  [FAIL] Test 1: Non-finite kernel results")
        ok = False

    # Test 2: kernel count = 64
    n = len(ALL_KERNELS)
    if n == 64:
        print(f"  [PASS] Test 2: Kernel count = {n}")
        passed += 1
    else:
        print(f"  [FAIL] Test 2: Expected 64 kernels, got {n}")
        ok = False

    # Test 3: exp gate fidelity in (0,1)
    val = kernel_exp_gate_fidelity()
    if 0 < val < 1:
        print(f"  [PASS] Test 3: F_gate = {val:.10f}")
        passed += 1
    else:
        print(f"  [FAIL] Test 3: F_gate = {val}")
        ok = False

    # Test 4: PImath key sum > 0
    val = kernel_pimath_key_gen()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 4: PImath key sum = {val:.0f}")
        passed += 1
    else:
        print(f"  [FAIL] Test 4: PImath = {val}")
        ok = False

    # Test 5: MUGE 3D gravity finite
    val = kernel_muge_3d_gravity()
    if math.isfinite(val):
        print(f"  [PASS] Test 5: MUGE 3D = {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 5: MUGE 3D = {val}")
        ok = False

    # Test 6: σ_n(ω,ρ) > σ_n(ω) (density enhances)
    sigma_base = kernel_sigma_n_density(rho=0)
    sigma_dense = kernel_sigma_n_density(rho=1e4)
    if sigma_dense >= sigma_base:
        print(f"  [PASS] Test 6: σ_n(ρ=1e4) = {sigma_dense:.6e} >= σ_n(ρ=0) = {sigma_base:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 6: Density should enhance cross section")
        ok = False

    # Test 7: BSFG geodesic finite and > 0
    val = kernel_bsfg_geodesic_trace()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 7: V_eff = {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 7: V_eff = {val}")
        ok = False

    # Test 8: Triadic master finite
    val = kernel_triadic_master_gpu()
    if math.isfinite(val):
        print(f"  [PASS] Test 8: Triadic master = {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 8: Triadic = {val}")
        ok = False

    # Test 9: φ_qubit Lagrangian finite
    val = kernel_phi_qubit_lagrangian()
    if math.isfinite(val):
        print(f"  [PASS] Test 9: L_φ_qubit = {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 9: L_φ_qubit = {val}")
        ok = False

    # Test 10: Appendix template > 0
    val = kernel_appendix_template()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 10: Appendix sections·S₂₆ = {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 10: appendix = {val}")
        ok = False

    # Test 11: Kernel names unique
    names = [k.__name__ for k in ALL_KERNELS]
    if len(names) == len(set(names)):
        print(f"  [PASS] Test 11: All {len(names)} kernel names unique")
        passed += 1
    else:
        print(f"  [FAIL] Test 11: Duplicate kernel names")
        ok = False

    # Test 12: Single pass finite
    bench = ProductionScalingV27()
    total = bench.single_pass()
    if math.isfinite(total):
        print(f"  [PASS] Test 12: Single pass total = {total:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 12: Non-finite pass total")
        ok = False

    # Test 13: Benchmark runs
    result = bench.simulate(0.5)
    if result['calc_per_sec'] > 0:
        label = "MET TARGET" if result['met_target'] else "below target"
        print(f"  [PASS] Test 13: {result['calc_per_sec']:.0f} calc/s ({label})")
        passed += 1
    else:
        print(f"  [FAIL] Test 13: Benchmark failed")
        ok = False

    print(f"\n  production_scaling_v27.py: {passed}/13 tests passed")
    return ok


if __name__ == "__main__":
    success = _run_tests()
    raise SystemExit(0 if success else 1)
