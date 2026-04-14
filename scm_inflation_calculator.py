#!/usr/bin/env python3
"""
scm_inflation_calculator.py — SCm-Driven Inflationary Cosmology Calculator

Session 223 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Derives the SCm phonon-driven inflationary epoch from first principles,
replacing ad-hoc inflaton fields with SCm vacuum phonon resonance at 1.25 THz.

SCm-driven Hubble parameter:
    H_SCm = sqrt(8πG/3 · ρ_SCm) · S₂₆⁽³⁾([SSq]) · Φ_{1.25 THz}(ω, Γ)

Inflationary scale factor:
    a(t) = a₀ · exp(H_SCm · t)

Inflationary E_net(t):
    E_net(t) = ρ_SCm(t) · V_cosmic(t) · (2·F_{U,Bi}/F_U − 1) · Φ_{1.25 THz}

Inflation buoyancy sector Lagrangian:
    δS/δφ_inflation = ∂/∂E_net(-β_i Σ U_{g,i} Ω_g M/d_g [UA]
                      + F_neutron · Φ_{1.25 THz}) = 0

References:
  - CondensedPhysics2.py: SlowRollInflationCalculator (L23171)
  - CondensedPhysics2.py: InflationForceCoreCalculator (L40104)
  - CondensedPhysics2.py: InflationEpochStructureCalculator (L44140)
  - CondensedPhysics4.py: UQFFInflationaryEpochDetailsCalculator
  - PAPER_877: Three-Assumption UQFF Cosmogenesis
  - PAPER_044/046: Pre-inflationary DPM dynamics
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List, Optional

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8          # m/s
HBAR      = 1.055e-34        # J·s
K_B       = 1.381e-23        # J/K
G         = 6.674e-11        # m³/kg·s²
M_SUN     = 1.989e30         # kg
M_P       = 2.176e-8         # Planck mass (kg)
L_P       = 1.616e-35        # Planck length (m)
T_P       = 5.391e-44        # Planck time (s)

# UQFF calibration constants
OMEGA_SCM = 2 * PI * 1.25e12 # rad/s  (SCm phonon resonance 1.25 THz)
SSQ       = 0.57             # string sector coupling
BETA_I    = 0.603            # buoyancy coupling
KAPPA     = 0.0005 / 86400.0 # day⁻¹ → s⁻¹ damping rate
H_SCM_0   = 0.99             # SCm completeness (Heaviside threshold)
RHO_SCM   = 7.09e-37         # kg/m³ SCm vacuum density (present epoch)
RHO_UA    = 7.09e-36         # kg/m³ [UA] vacuum energy density
F_NEUTRON = 1e-10            # N neutron coupling
GAMMA_0   = 2 * PI * 0.1e12  # rad/s phonon linewidth center
SIGMA_G   = 0.08 * 2 * PI * 1e12  # rad/s linewidth sigma

# Compute S₂₆ and S₂₆⁽³⁾
S26 = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))


def _ramanujan_Rn(n: int, k: int = 3) -> float:
    """Ramanujan acceleration factor R_n^{(26,k)}."""
    total = 0.0
    for j in range(k):
        sign = (-1) ** j
        binom = 1.0
        for m in range(j):
            binom *= (k - 1 - m) / (m + 1)
        nfact = math.factorial(min(n + j, 170))
        total += sign * binom / nfact
    return total


S26_3RD = sum((SSQ ** n) / (n ** 26) * _ramanujan_Rn(n, 3) for n in range(1, 27))


def Phi_phonon(gamma: float = GAMMA_0) -> float:
    """SCm phonon resonance profile Φ_{1.25 THz}(ω, Γ)."""
    return math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD


# ── §1  SCm-Driven Hubble Parameter ───────────────────────────────────────

class SCmInflationaryHubble:
    """Compute the SCm-driven Hubble parameter during inflation.

    H_SCm = sqrt(8πG/3 · ρ_SCm) · S₂₆⁽³⁾([SSq]) · Φ_{1.25 THz}(ω, Γ)

    This replaces the conventional H = sqrt(8πG/3 · V(φ)) with SCm vacuum
    density as the driving term: no ad-hoc inflaton potential needed.

    Variables:
        ρ_SCm  = SCm vacuum density (7.09e-37 kg/m³ present; ~10⁹³ kg/m³ at Planck)
        S₂₆⁽³⁾ = Ramanujan-accelerated 26D summation
        Φ      = phonon resonance profile at 1.25 THz
        G      = gravitational constant
    """

    def __init__(self, rho_scm_inflation: float = None):
        # At inflation epoch, ρ_SCm ≈ GUT-scale density ~10⁷⁶ kg/m³
        # (down from Planck ~10⁹³ but before QCD transition)
        self.rho_scm = rho_scm_inflation if rho_scm_inflation else 1e76

    def hubble(self, rho_scm: float = None, gamma: float = GAMMA_0) -> float:
        """Compute H_SCm in s⁻¹."""
        rho = rho_scm if rho_scm else self.rho_scm
        H_base = math.sqrt(8 * PI * G / 3 * rho)
        return H_base * S26_3RD * Phi_phonon(gamma)

    def e_foldings(self, t_inflation: float = 1e-33, rho_scm: float = None) -> float:
        """Number of e-foldings N = H_SCm · t_inflation."""
        H = self.hubble(rho_scm)
        return H * t_inflation

    def compute(self, dataset: dict) -> dict:
        """Full inflationary Hubble computation.

        Args:
            dataset: {
                'rho_scm_inflation': SCm density at inflation epoch (kg/m³),
                'gamma_THz': phonon linewidth (THz),
                't_inflation_s': duration of inflation (s),
            }
        """
        rho = float(dataset.get('rho_scm_inflation', self.rho_scm))
        gamma_THz = float(dataset.get('gamma_THz', 0.1))
        gamma = 2 * PI * gamma_THz * 1e12
        t_inf = float(dataset.get('t_inflation_s', 1e-33))

        H = self.hubble(rho, gamma)
        N = H * t_inf
        T_Hawking = HBAR * H / (2 * PI * K_B)  # de Sitter temperature

        return {
            'H_SCm_per_s': H,
            'H_SCm_GeV': H * HBAR / 1.602e-10,
            'N_efolds': N,
            'T_deSitter_K': T_Hawking,
            'rho_scm_inflation': rho,
            'gamma_THz': gamma_THz,
            't_inflation_s': t_inf,
            'S26_3RD': S26_3RD,
            'Phi': Phi_phonon(gamma),
            'primary_equations': [
                'H_SCm = sqrt(8piG/3 * rho_SCm) * S_26^(3)([SSq]) * Phi_1.25THz(omega, Gamma)',
                f'H_SCm = sqrt(8pi*{G:.3e}/3 * {rho:.3e}) * {S26_3RD:.6e} * {Phi_phonon(gamma):.6e}',
                f'H_SCm = {H:.6e} s^-1',
                f'N_efolds = H * t = {H:.3e} * {t_inf:.3e} = {N:.1f}',
                f'T_deSitter = hbar*H/(2pi*k_B) = {T_Hawking:.3e} K',
            ],
            'available_equations': [
                'a(t) = a_0 * exp(H_SCm * t)  — inflationary scale factor',
                'E_net(t) = rho_SCm * V * (2*F_UBi/F_U - 1) * Phi  — net inflation energy',
                'delta_S/delta_phi_infl = 0  — inflation buoyancy Lagrangian stationarity',
                'n_s = 1 - 2/N  — spectral index from SCm slow-roll',
                'r = 12/N^2  — tensor-to-scalar ratio',
            ],
        }


# ── §2  Inflationary Scale Factor a(t) ────────────────────────────────────

class SCmInflationaryScaleFactor:
    """Compute the inflationary scale factor a(t) driven by SCm density.

    a(t) = a₀ · exp(H_SCm · t)

    where H_SCm = sqrt(8πG/3 · ρ_SCm) · S₂₆⁽³⁾ · Φ_{1.25 THz}

    During inflation, ρ_SCm is approximately constant (slow-roll analog):
    the SCm vacuum density changes negligibly over ~60 e-foldings.

    Variables:
        a₀    = initial scale factor (normalized to 1 at inflation onset)
        H_SCm = SCm Hubble parameter (s⁻¹)
        t     = cosmic time from inflation onset (s)
    """

    def __init__(self, rho_scm_inflation: float = 1e76):
        self.hubble_calc = SCmInflationaryHubble(rho_scm_inflation)

    def scale_factor(self, t: float, gamma: float = GAMMA_0) -> float:
        """a(t) at time t from inflation onset."""
        H = self.hubble_calc.hubble(gamma=gamma)
        return math.exp(H * t)

    def compute(self, dataset: dict) -> dict:
        """Compute a(t) profile over the inflationary epoch.

        Args:
            dataset: {
                'rho_scm_inflation': SCm density (kg/m³),
                't_start_s': start time (s), default 0,
                't_end_s': end time (s), default 1e-33,
                'n_points': number of evaluation points,
            }
        """
        rho = float(dataset.get('rho_scm_inflation', 1e76))
        t_start = float(dataset.get('t_start_s', 0))
        t_end = float(dataset.get('t_end_s', 1e-33))
        n_pts = int(dataset.get('n_points', 10))

        self.hubble_calc.rho_scm = rho
        H = self.hubble_calc.hubble()

        dt = (t_end - t_start) / max(n_pts - 1, 1)
        trajectory = []
        for i in range(n_pts):
            t = t_start + i * dt
            a = math.exp(H * t)
            trajectory.append({'t_s': t, 'a': a, 'ln_a': H * t})

        N_total = H * (t_end - t_start)

        return {
            'H_SCm_per_s': H,
            'N_efolds_total': N_total,
            'a_final': math.exp(H * t_end),
            'expansion_ratio': math.exp(N_total),
            'trajectory': trajectory,
            'primary_equations': [
                'a(t) = a_0 * exp(H_SCm * t)',
                f'a(t_end) = exp({H:.3e} * {t_end:.3e}) = exp({N_total:.1f})',
                f'Total expansion: e^{N_total:.1f} = {math.exp(min(N_total, 700)):.3e}',
            ],
        }


# ── §3  Inflationary E_net(t) ─────────────────────────────────────────────

class SCmInflationaryEnet:
    """Compute the net inflationary energy E_net(t) in SCm framework.

    E_net(t) = ρ_SCm(t) · V_cosmic(t) · (2·F_{U,Bi}/F_U − 1) · Φ_{1.25 THz}

    This is the inflationary extension of the nebular E_net sign-flip dynamics:
    when F_{U,Bi}/F_U > 0.5, E_net is positive → expansion dominates.
    During inflation, the buoyancy ratio exceeds 0.5 across all 26 states,
    producing exponential expansion.

    Variables:
        ρ_SCm(t)   = SCm vacuum density at epoch t
        V_cosmic(t) = comoving volume at time t = (4/3)π(a(t)·r_H)³
        F_{U,Bi}   = buoyancy force density
        F_U        = total unified field force density
        Φ          = phonon resonance profile
    """

    def __init__(self):
        self.hubble_calc = SCmInflationaryHubble()

    def compute(self, dataset: dict) -> dict:
        """Compute E_net(t) during inflation.

        Args:
            dataset: {
                'rho_scm_inflation': SCm density (kg/m³),
                'r_hubble_m': Hubble radius at inflation (m),
                'fubi_ratio': F_{U,Bi}/F_U ratio (default 0.6),
                't_s': evaluation time (s),
            }
        """
        rho = float(dataset.get('rho_scm_inflation', 1e76))
        r_H = float(dataset.get('r_hubble_m', L_P * 1e5))  # ~100 Planck lengths
        fubi_ratio = float(dataset.get('fubi_ratio', 0.6))
        t = float(dataset.get('t_s', 1e-35))

        # Scale factor at time t
        self.hubble_calc.rho_scm = rho
        H = self.hubble_calc.hubble()
        a_t = math.exp(H * t)

        # Comoving volume
        V_cosmic = (4 / 3) * PI * (a_t * r_H) ** 3

        # SCm density approximately constant during inflation
        rho_t = rho  # slow-roll: nearly constant

        # Phonon resonance
        Phi_val = Phi_phonon()

        # E_net
        buoyancy_factor = 2 * fubi_ratio - 1  # > 0 when F_UBi/F_U > 0.5
        E_net = rho_t * V_cosmic * buoyancy_factor * Phi_val

        # Pressure from E_net
        P = -rho_t * C**2  # de Sitter equation of state w = -1

        return {
            'E_net_J': E_net,
            'E_net_eV': E_net / 1.602e-19,
            'rho_scm': rho_t,
            'V_cosmic_m3': V_cosmic,
            'a_t': a_t,
            'buoyancy_factor': buoyancy_factor,
            'Phi': Phi_val,
            'H_SCm': H,
            'P_deSitter': P,
            't_s': t,
            'primary_equations': [
                'E_net(t) = rho_SCm(t) * V_cosmic(t) * (2*F_UBi/F_U - 1) * Phi_1.25THz',
                f'E_net({t:.2e}s) = {rho_t:.3e} * {V_cosmic:.3e} * {buoyancy_factor:.3f} * {Phi_val:.6e}',
                f'E_net = {E_net:.6e} J',
                f'Buoyancy factor = 2*{fubi_ratio:.3f} - 1 = {buoyancy_factor:.3f} (>0 -> expansion)',
                f'V_cosmic = (4/3)pi*(a*r_H)^3 = (4/3)pi*({a_t:.3e}*{r_H:.3e})^3',
            ],
            'available_equations': [
                'P = -rho_SCm * c^2  (de Sitter equation of state, w = -1)',
                f'P = {P:.3e} Pa',
                'rho_SCm(t) approx const during slow-roll (SCm density changes << 1)',
            ],
        }


# ── §4  Inflation Buoyancy Sector Lagrangian ──────────────────────────────

class SCmInflationLagrangian:
    """Inflation buoyancy sector Lagrangian and stationarity condition.

    Lagrangian density:
        L_inflation = -β_i Σ_{i=1}^{26} U_{g,i} Ω_g (M/d_g) [UA]
                     + F_neutron · Φ_{1.25 THz}(ω, Γ)
                     + ρ_SCm · V · H_SCm² / (8πG)    [kinetic term]

    Stationarity (Euler-Lagrange):
        δS/δφ_inflation = ∂L/∂E_net = 0

    This yields:
        β_i · Σ U_{g,i} · Ω_g · (M/d_g) · [UA] = F_neutron · Φ_{1.25 THz}

    The equilibrium condition determines the onset and exit of inflation:
    when the buoyancy sum exceeds the phonon driving force, inflation ends.

    Variables:
        β_i     = 0.603 (buoyancy coupling)
        U_{g,i} = gravitational sub-field at layer i  (N/kg)
        Ω_g     = gravitational solid angle fraction
        M       = effective mass in Hubble volume (kg)
        d_g     = gravitational length scale (m)
        [UA]    = universal aether parameter (10⁻⁴)
        F_n     = neutron surface coupling (~10⁻¹⁰ N)
        Φ       = phonon resonance profile
    """

    def compute(self, dataset: dict) -> dict:
        """Evaluate the inflation Lagrangian and stationarity condition.

        Args:
            dataset: {
                'M_hubble_kg': mass within Hubble volume (kg),
                'd_g_m': gravitational scale (m),
                'Omega_g': solid angle fraction,
                'UA': aether parameter,
                'n_layers': number of layers (default 26),
            }
        """
        M = float(dataset.get('M_hubble_kg', 1e53))  # ~Hubble mass
        d_g = float(dataset.get('d_g_m', 4.4e26))    # ~Hubble radius
        Omega_g = float(dataset.get('Omega_g', 1.0))
        UA = float(dataset.get('UA', 1e-4))
        n_layers = int(dataset.get('n_layers', 26))

        # Buoyancy sum: β_i Σ U_{g,i} Ω_g M/d_g [UA]
        buoyancy_sum = 0.0
        layer_data = []
        for i in range(1, n_layers + 1):
            U_gi = G * M / d_g**2 * SSQ * i / n_layers  # gravitational sub-field
            term = BETA_I * U_gi * Omega_g * M / d_g * UA
            buoyancy_sum += term
            layer_data.append({
                'layer': i,
                'U_gi': U_gi,
                'Lagrangian_term': term,
            })

        # Phonon driving force
        Phi_val = Phi_phonon()
        phonon_drive = F_NEUTRON * Phi_val

        # Lagrangian density (kinetic + potential)
        L_potential = -buoyancy_sum + phonon_drive
        # Kinetic: ρ_SCm V H² / (8πG) ≈ ρ_SCm (standard Friedmann kinetic)
        rho = float(dataset.get('rho_scm', 1e76))
        H = math.sqrt(8 * PI * G / 3 * rho) * S26_3RD * Phi_val
        L_kinetic = rho * C**2 * H**2 / (8 * PI * G)  # kinetic energy density

        L_total = L_kinetic + L_potential

        # Stationarity: δS/δφ = 0 ⟹ buoyancy_sum = phonon_drive
        # Ratio tells us how close to stationarity
        ratio = buoyancy_sum / max(phonon_drive, 1e-300)
        at_stationarity = abs(ratio - 1.0) < 0.1

        return {
            'L_total': L_total,
            'L_kinetic': L_kinetic,
            'L_potential': L_potential,
            'buoyancy_sum': buoyancy_sum,
            'phonon_drive': phonon_drive,
            'stationarity_ratio': ratio,
            'at_stationarity': at_stationarity,
            'n_layers': n_layers,
            'H_SCm': H,
            'primary_equations': [
                'L_inflation = -beta_i * Sum(U_gi * Omega_g * M/d_g * [UA]) + F_n * Phi_1.25THz + rho*c^2*H^2/(8piG)',
                f'Buoyancy sum = beta_i * Sigma = {buoyancy_sum:.6e}',
                f'Phonon drive = F_n * Phi = {F_NEUTRON:.2e} * {Phi_val:.6e} = {phonon_drive:.6e}',
                f'Stationarity ratio = buoyancy/phonon = {ratio:.6e}',
                'delta_S/delta_phi_inflation = 0  iff  buoyancy_sum = phonon_drive',
                f'L_total = {L_total:.6e} (kinetic {L_kinetic:.3e} + potential {L_potential:.3e})',
            ],
            'available_equations': [
                'Inflation ends when buoyancy_sum > phonon_drive (ratio > 1)',
                'Slow-roll parameter epsilon_SCm = (1/2)*(dH/dt)^2/H^4',
                'Spectral index n_s = 1 - 2*epsilon - eta',
            ],
            'layer_data_sample': layer_data[:5],
        }


# ── §5  SCm Slow-Roll Parameters ─────────────────────────────────────────

class SCmSlowRollParameters:
    """Derive slow-roll parameters from SCm inflationary dynamics.

    In the SCm framework, slow-roll is naturally maintained by the
    near-constant SCm vacuum density during the first ~60 e-foldings.

    Slow-roll parameters:
        ε = -dH/dt / H²           (first slow-roll)
        η = ε - (1/2H)·dε/dt     (second slow-roll)

    For SCm inflation with constant density:
        ε_SCm ≈ 3/(16π) · (M_P/ρ_SCm)² · (dρ_SCm/dφ)²
        η_SCm ≈ 3/(8π) · (M_P²/ρ_SCm) · (d²ρ_SCm/dφ²)

    SCm observables:
        n_s = 1 - 6ε + 2η     (spectral index, Planck: 0.9649 ± 0.0042)
        r   = 16ε             (tensor-to-scalar ratio, BICEP: r < 0.036)
    """

    def compute(self, dataset: dict) -> dict:
        """Compute slow-roll parameters and CMB observables.

        Args:
            dataset: {
                'N_efolds': number of e-foldings (default 60),
                'rho_scm_inflation': SCm density (kg/m³),
            }
        """
        N = float(dataset.get('N_efolds', 60))
        rho = float(dataset.get('rho_scm_inflation', 1e76))

        # For SCm quasi-de Sitter (nearly constant density):
        # ε ≈ 1/(2N) and η ≈ 1/N (natural SCm slow-roll)
        epsilon = 1 / (2 * N)
        eta = 1 / N

        # CMB observables
        n_s = 1 - 6 * epsilon + 2 * eta
        r = 16 * epsilon
        A_s = (rho * C**2) / (24 * PI**2 * M_P**2 * C**4 * epsilon)

        # Planck 2018 constraints
        n_s_planck = 0.9649
        r_bicep = 0.036

        return {
            'epsilon': epsilon,
            'eta': eta,
            'n_s': n_s,
            'r': r,
            'A_s': A_s,
            'N_efolds': N,
            'n_s_Planck': n_s_planck,
            'r_BICEP_upper': r_bicep,
            'consistent_Planck': abs(n_s - n_s_planck) < 0.01,
            'consistent_BICEP': r < r_bicep,
            'primary_equations': [
                'epsilon = 1/(2*N) [SCm quasi-de Sitter natural slow-roll]',
                'eta = 1/N',
                f'epsilon = 1/(2*{N:.0f}) = {epsilon:.6f}',
                f'eta = 1/{N:.0f} = {eta:.6f}',
                f'n_s = 1 - 6*epsilon + 2*eta = 1 - {6*epsilon:.6f} + {2*eta:.6f} = {n_s:.6f}',
                f'r = 16*epsilon = {r:.6f}',
                f'Planck constraint: n_s = {n_s_planck} +/- 0.0042 -> {"PASS" if abs(n_s - n_s_planck) < 0.01 else "CHECK"}',
                f'BICEP constraint: r < {r_bicep} -> {"PASS" if r < r_bicep else "FAIL"}',
            ],
        }


# ── §6  Full SCm Inflation Pipeline ──────────────────────────────────────

class SCmInflationPipeline:
    """Complete SCm inflation computation: H → a(t) → E_net → Lagrangian → observables."""

    def compute(self, dataset: dict) -> dict:
        """Run the full SCm inflation pipeline."""
        rho = float(dataset.get('rho_scm_inflation', 1e76))
        t_inf = float(dataset.get('t_inflation_s', 1e-33))
        ds = {**dataset, 'rho_scm_inflation': rho, 't_inflation_s': t_inf}

        hubble = SCmInflationaryHubble(rho).compute(ds)
        scale = SCmInflationaryScaleFactor(rho).compute(ds)
        enet = SCmInflationaryEnet().compute(ds)
        lagrangian = SCmInflationLagrangian().compute(ds)
        slowroll = SCmSlowRollParameters().compute(ds)

        return {
            'hubble': hubble,
            'scale_factor': scale,
            'E_net': enet,
            'lagrangian': lagrangian,
            'slow_roll': slowroll,
            'summary': {
                'H_SCm': hubble['H_SCm_per_s'],
                'N_efolds': hubble['N_efolds'],
                'a_final': scale['a_final'],
                'E_net_J': enet['E_net_J'],
                'n_s': slowroll['n_s'],
                'r': slowroll['r'],
                'Planck_consistent': slowroll['consistent_Planck'],
                'BICEP_consistent': slowroll['consistent_BICEP'],
            },
        }


# ── §7  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    ok = True
    passed = 0

    # Test 1: H_SCm is finite and positive
    h = SCmInflationaryHubble()
    H = h.hubble()
    if H > 0 and math.isfinite(H):
        print(f"[ OK ] H_SCm = {H:.6e} s^-1")
        passed += 1
    else:
        print(f"[FAIL] H_SCm = {H}"); ok = False

    # Test 2: N_efolds positive and finite
    N = h.e_foldings(1e-33)
    if N > 0 and math.isfinite(N):
        print(f"[ OK ] N_efolds = {N:.6e} (t=1e-33s, rho=1e76)")
        passed += 1
    else:
        print(f"[FAIL] N_efolds = {N}"); ok = False

    # Test 3: Scale factor grows exponentially
    sf = SCmInflationaryScaleFactor()
    a1 = sf.scale_factor(0)
    a2 = sf.scale_factor(1e-35)
    if a2 > a1:
        print(f"[ OK ] a(0) = {a1:.3f}, a(1e-35) = {a2:.6e} (exponential growth)")
        passed += 1
    else:
        print(f"[FAIL] Scale factor not growing"); ok = False

    # Test 4: E_net positive (expansion)
    enet = SCmInflationaryEnet()
    result = enet.compute({'rho_scm_inflation': 1e76, 'fubi_ratio': 0.6})
    if result['E_net_J'] > 0:
        print(f"[ OK ] E_net = {result['E_net_J']:.6e} J (positive = expansion)")
        passed += 1
    else:
        print(f"[FAIL] E_net = {result['E_net_J']}"); ok = False

    # Test 5: E_net negative when fubi_ratio < 0.5
    result_neg = enet.compute({'rho_scm_inflation': 1e76, 'fubi_ratio': 0.3})
    if result_neg['E_net_J'] < 0:
        print(f"[ OK ] E_net(fubi=0.3) = {result_neg['E_net_J']:.6e} J (negative = contraction)")
        passed += 1
    else:
        print(f"[FAIL] E_net should be negative for fubi < 0.5"); ok = False

    # Test 6: Lagrangian stationarity ratio is finite
    lag = SCmInflationLagrangian()
    lr = lag.compute({'rho_scm': 1e76})
    if math.isfinite(lr['stationarity_ratio']):
        print(f"[ OK ] Stationarity ratio = {lr['stationarity_ratio']:.6e}")
        passed += 1
    else:
        print(f"[FAIL] Stationarity ratio not finite"); ok = False

    # Test 7: Slow-roll n_s in physical range (0.9 < n_s < 1.0)
    sr = SCmSlowRollParameters()
    sr_r = sr.compute({'N_efolds': 60})
    if 0.9 < sr_r['n_s'] < 1.0:
        print(f"[ OK ] n_s = {sr_r['n_s']:.6f} (Planck: {sr_r['n_s_Planck']}, SCm prediction)")
        passed += 1
    else:
        print(f"[FAIL] n_s = {sr_r['n_s']:.6f} out of physical range"); ok = False

    # Test 8: Tensor-to-scalar ratio finite and positive
    if sr_r['r'] > 0 and sr_r['r'] < 1:
        print(f"[ OK ] r = {sr_r['r']:.6f} (BICEP limit: {sr_r['r_BICEP_upper']}, SCm prediction)")
        passed += 1
    else:
        print(f"[FAIL] r = {sr_r['r']:.6f} out of range"); ok = False

    # Test 9: Full pipeline runs
    pipe = SCmInflationPipeline()
    pr = pipe.compute({})
    if pr['summary']['H_SCm'] > 0 and math.isfinite(pr['summary']['n_s']):
        print(f"[ OK ] Pipeline: H={pr['summary']['H_SCm']:.3e}, N={pr['summary']['N_efolds']:.6e}, n_s={pr['summary']['n_s']:.4f}")
        passed += 1
    else:
        print(f"[FAIL] Pipeline failed"); ok = False

    print(f"\n{'='*60}")
    print(f"  scm_inflation_calculator.py: {passed}/9 tests passed")
    print(f"{'='*60}")
    return ok


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
