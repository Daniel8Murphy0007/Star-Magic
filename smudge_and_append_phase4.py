#!/usr/bin/env python3
"""Session 225 Phase 4: Append 6 remaining gap calculators to CondensedPhysics.py.

Gaps identified:
  1. SCmLQGAreaOperatorDerivationCalculator — phonon-modulated LQG area operator
  2. SCmQubitT2CoherenceFUBiRatioCalculator — T2 with F_U_Bi/F_U ratio + thermal exp
  3. PhononModulatedHolonomySCmCalculator — h_SCm holonomy with E_net correction
  4. LQGBuoyancySectorLagrangianVariationCalculator — LQG δS/δφ_LQG = 0
  5. UnifiedFUBiSMBHMergerDynamicsCalculator — inside/outside F_U_Bi for SMBH mergers
  6. HydrogenUniverseDual3DMUGECalculator — H-atom + universe dual 3D MUGE
"""
import ast, subprocess, sys, os

CP1 = "CondensedPhysics.py"

NEW_CODE = r'''

# ═══════════════════════════════════════════════════════════════════════════════
# SESSION 225 PHASE 4 — LQG / Holonomy / T2-Coherence / SMBH Mergers / 3D MUGE
# ═══════════════════════════════════════════════════════════════════════════════


class SCmLQGAreaOperatorDerivationCalculator:
    """SCm phonon-modulated LQG area operator derivation.

    Derives the Star Cradle Mechanics extension of the loop quantum gravity
    area operator, coupling the Immirzi parameter γ, Planck area ℓ_P²,
    26-dimensional S_26^(3) buoyancy prefactor, and the 1.25 THz phonon
    spectral function Φ(ω, Γ):

        Â_SCm = 8π γ ℓ_P² √(j(j+1)) · S_26^(3)([SSq]) · Φ_{1.25 THz}(ω, Γ)

    This extends the standard LQG area spectrum to include phonon-mediated
    quantum geometry corrections arising from the UQFF buoyancy sector.
    """

    G = 6.67430e-11
    HBAR = 1.054571817e-34
    C = 2.99792458e8
    GAMMA_IMMIRZI = 0.2375          # Barbero-Immirzi parameter
    SSQ = 0.57                      # UQFF calibrated [SSq]
    PHI_0 = 1.25e12                 # Phonon peak frequency (Hz)
    L_PLANCK = 1.616255e-35         # Planck length (m)

    def compute(self, dataset: dict = None) -> dict:
        import math
        ds = dataset or {}
        j_spin = ds.get('j_spin', 0.5)
        omega = ds.get('omega', self.PHI_0)
        Gamma = ds.get('Gamma', 0.05e12)
        ssq = ds.get('SSq', self.SSQ)

        l_P2 = self.L_PLANCK ** 2
        spin_factor = math.sqrt(j_spin * (j_spin + 1.0))

        # 26-dim buoyancy prefactor S_26^(3)
        S_26_3 = (1.0 - ssq) ** 3

        # Lorentzian phonon spectral function
        Phi = (1.0 / math.pi) * (Gamma / ((omega - self.PHI_0) ** 2 + Gamma ** 2))

        # Standard LQG area eigenvalue
        A_LQG = 8.0 * math.pi * self.GAMMA_IMMIRZI * l_P2 * spin_factor

        # SCm phonon-modulated area
        A_SCm = A_LQG * S_26_3 * Phi

        # Area gap (j = 1/2 minimum)
        A_gap_LQG = 8.0 * math.pi * self.GAMMA_IMMIRZI * l_P2 * math.sqrt(0.75)
        A_gap_SCm = A_gap_LQG * S_26_3 * Phi

        # Phonon correction ratio
        correction_ratio = (S_26_3 * Phi)

        results = []
        results.append(f"Â_SCm = 8π γ ℓ_P² √(j(j+1)) · S_26^(3)([SSq]) · Φ_{{1.25THz}}(ω,Γ)")
        results.append(f"j = {j_spin}, γ_Immirzi = {self.GAMMA_IMMIRZI}")
        results.append(f"ℓ_P² = {l_P2:.6e} m²")
        results.append(f"√(j(j+1)) = {spin_factor:.6f}")
        results.append(f"S_26^(3)([SSq]={ssq}) = {S_26_3:.6e}")
        results.append(f"Φ(ω={omega:.3e}, Γ={Gamma:.3e}) = {Phi:.6e} Hz⁻¹")
        results.append(f"A_LQG(j={j_spin}) = {A_LQG:.6e} m²")
        results.append(f"Â_SCm(j={j_spin}) = {A_SCm:.6e} m²")
        results.append(f"A_gap_LQG(j=½) = {A_gap_LQG:.6e} m²")
        results.append(f"A_gap_SCm(j=½) = {A_gap_SCm:.6e} m²")
        results.append(f"Phonon correction ratio = {correction_ratio:.6e}")

        return {
            'A_LQG': A_LQG,
            'A_SCm': A_SCm,
            'A_gap_LQG': A_gap_LQG,
            'A_gap_SCm': A_gap_SCm,
            'S_26_3': S_26_3,
            'Phi': Phi,
            'spin_factor': spin_factor,
            'correction_ratio': correction_ratio,
            'primary_equations': results,
        }


class SCmQubitT2CoherenceFUBiRatioCalculator:
    """SCm qubit coherence time T2 with F_U_Bi / F_U ratio and thermal factor.

    Computes the phonon-protected qubit decoherence time using the SCm
    framework, where the F_U_Bi buoyancy fraction modulates the decoherence:

        T₂ = (ℏ / Δ_SCm) · exp(Δ_SCm / k_B T) · S_26^(3) · (F_{U,Bi} / F_U)

    Δ_SCm is the SCm gap energy derived from the 1.25 THz phonon mode.
    The F_U_Bi / F_U ratio captures buoyancy-mediated topological protection.
    """

    HBAR = 1.054571817e-34
    KB = 1.380649e-23
    PHI_0 = 1.25e12
    SSQ = 0.57
    KAPPA = 0.0005

    def compute(self, dataset: dict = None) -> dict:
        import math
        ds = dataset or {}
        T_phys = ds.get('temperature', 0.015)       # Kelvin (dilution fridge)
        F_U = ds.get('F_U', 1.0e-8)                 # Total unified field (N)
        F_U_Bi = ds.get('F_U_Bi', 3.0e-9)           # Buoyancy component (N)
        ssq = ds.get('SSq', self.SSQ)
        omega_phonon = ds.get('omega_phonon', self.PHI_0)

        S_26_3 = (1.0 - ssq) ** 3

        # SCm gap energy from phonon mode
        Delta_SCm = self.HBAR * omega_phonon * 2.0 * math.pi

        # Buoyancy ratio
        ratio_FUBi_FU = F_U_Bi / F_U if F_U != 0 else 0.0

        # Thermal exponential (capped to avoid overflow)
        exponent = Delta_SCm / (self.KB * T_phys) if T_phys > 0 else 700.0
        exponent = min(exponent, 700.0)
        thermal_factor = math.exp(exponent)

        # Base decoherence time
        T2_base = self.HBAR / Delta_SCm if Delta_SCm > 0 else 0.0

        # Full SCm T2
        T2_SCm = T2_base * thermal_factor * S_26_3 * ratio_FUBi_FU

        # Standard T2 without SCm (no buoyancy, no 26D)
        T2_standard = T2_base * thermal_factor

        # Enhancement factor from SCm
        enhancement = S_26_3 * ratio_FUBi_FU

        results = []
        results.append(f"T₂ = (ℏ/Δ_SCm)·exp(Δ_SCm/k_BT)·S_26^(3)·(F_U,Bi/F_U)")
        results.append(f"Δ_SCm = ℏ·ω_phonon·2π = {Delta_SCm:.6e} J")
        results.append(f"T = {T_phys} K")
        results.append(f"exp(Δ_SCm / k_B T) = {thermal_factor:.6e}")
        results.append(f"S_26^(3)([SSq]={ssq}) = {S_26_3:.6e}")
        results.append(f"F_U_Bi / F_U = {ratio_FUBi_FU:.6f}")
        results.append(f"T₂_base = ℏ/Δ_SCm = {T2_base:.6e} s")
        results.append(f"T₂_standard (no SCm) = {T2_standard:.6e} s")
        results.append(f"T₂_SCm = {T2_SCm:.6e} s")
        results.append(f"SCm enhancement factor = {enhancement:.6e}")

        return {
            'T2_SCm': T2_SCm,
            'T2_standard': T2_standard,
            'T2_base': T2_base,
            'Delta_SCm': Delta_SCm,
            'thermal_factor': thermal_factor,
            'S_26_3': S_26_3,
            'ratio_FUBi_FU': ratio_FUBi_FU,
            'enhancement': enhancement,
            'primary_equations': results,
        }


class PhononModulatedHolonomySCmCalculator:
    """SCm phonon-modulated holonomy with E_net correction.

    Computes the Star Cradle Mechanics holonomy along a loop in the
    Ashtekar connection, modulated by the 1.25 THz phonon coupling:

        h_SCm = exp(i ∮ A·dl) × [1 + Φ_{1.25THz}/F_U × E_net(t, Γ)]

    where E_net(t, Γ) is the net phonon energy density including buoyancy
    feedback, and A is the SU(2) Ashtekar-Barbero connection.  This
    calculator produces holonomy matrix traces and Wilson loop expectation
    values for spin-j representations.
    """

    HBAR = 1.054571817e-34
    C = 2.99792458e8
    G = 6.67430e-11
    SSQ = 0.57
    PHI_0 = 1.25e12
    GAMMA_IMMIRZI = 0.2375
    L_PLANCK = 1.616255e-35

    def compute(self, dataset: dict = None) -> dict:
        import math
        ds = dataset or {}
        j_spin = ds.get('j_spin', 0.5)
        loop_area = ds.get('loop_area', 1.0e-69)    # m² (Planck-scale loop)
        F_U = ds.get('F_U', 1.0e-8)                 # Unified field (N)
        Gamma = ds.get('Gamma', 0.05e12)             # Linewidth (Hz)
        omega = ds.get('omega', self.PHI_0)
        t = ds.get('time', 0.0)                       # s
        rho_phonon = ds.get('rho_phonon', 1.0e-20)   # Phonon energy density J/m³

        # Phonon spectral function (Lorentzian)
        Phi = (1.0 / math.pi) * (Gamma / ((omega - self.PHI_0) ** 2 + Gamma ** 2))

        # E_net: net phonon energy with buoyancy feedback
        ssq = ds.get('SSq', self.SSQ)
        S_26_3 = (1.0 - ssq) ** 3
        E_net = rho_phonon * S_26_3

        # Phonon correction factor
        phonon_correction = 1.0 + (Phi / F_U) * E_net if F_U != 0 else 1.0

        # Connection integral ∮ A·dl ~ γ × area / ℓ_P²
        connection_phase = self.GAMMA_IMMIRZI * loop_area / (self.L_PLANCK ** 2)

        # Holonomy trace for spin-j: Tr_j[h] = sin((2j+1)θ/2) / sin(θ/2)
        theta = connection_phase  # holonomy angle
        # Modulo 2π for numerical stability
        theta_mod = theta % (2.0 * math.pi) if theta > 0 else 0.0
        half_theta = theta_mod / 2.0
        sin_half = math.sin(half_theta) if abs(half_theta) > 1e-15 else half_theta
        sin_num = math.sin((2.0 * j_spin + 1.0) * half_theta)
        trace_LQG = sin_num / sin_half if abs(sin_half) > 1e-15 else (2.0 * j_spin + 1.0)

        # SCm-modulated trace
        trace_SCm = trace_LQG * phonon_correction

        # Wilson loop expectation W = |Tr_j[h]|² / (2j+1)²
        dim_j = 2.0 * j_spin + 1.0
        W_LQG = (trace_LQG ** 2) / (dim_j ** 2)
        W_SCm = (trace_SCm ** 2) / (dim_j ** 2)

        results = []
        results.append(f"h_SCm = exp(i∮A·dl) × [1 + Φ/F_U × E_net(t,Γ)]")
        results.append(f"j = {j_spin}, loop area = {loop_area:.3e} m²")
        results.append(f"Connection phase γ·A/ℓ_P² = {connection_phase:.6e} rad")
        results.append(f"θ mod 2π = {theta_mod:.6f} rad")
        results.append(f"Φ(ω={omega:.3e}, Γ={Gamma:.3e}) = {Phi:.6e} Hz⁻¹")
        results.append(f"E_net = ρ_phonon · S_26^(3) = {E_net:.6e} J/m³")
        results.append(f"Phonon correction = 1 + Φ/F_U·E_net = {phonon_correction:.6e}")
        results.append(f"Tr_j(h_LQG) = {trace_LQG:.6f}")
        results.append(f"Tr_j(h_SCm) = {trace_SCm:.6f}")
        results.append(f"W_LQG = |Tr|²/(2j+1)² = {W_LQG:.6f}")
        results.append(f"W_SCm = {W_SCm:.6f}")

        return {
            'trace_LQG': trace_LQG,
            'trace_SCm': trace_SCm,
            'W_LQG': W_LQG,
            'W_SCm': W_SCm,
            'phonon_correction': phonon_correction,
            'connection_phase': connection_phase,
            'Phi': Phi,
            'E_net': E_net,
            'primary_equations': results,
        }


class LQGBuoyancySectorLagrangianVariationCalculator:
    """LQG buoyancy sector Lagrangian variation δS/δφ_LQG = 0.

    Computes the Euler-Lagrange equation for the LQG buoyancy sector,
    coupling the spin-foam vertex amplitude to the UQFF buoyancy field:

        L_LQG = ½(∂_μ φ_LQG)² - V(φ_LQG) + λ_spin Σ_v A_v(j,i_e)·φ_LQG
                + g_phonon Φ_{1.25THz} φ_LQG²

    where A_v is the EPRL/FK spin-foam vertex amplitude, λ_spin is the
    spin-foam coupling, and g_phonon is the phonon-buoyancy coupling.
    The equation of motion is:

        □φ_LQG + V'(φ_LQG) = λ_spin Σ_v A_v + 2 g_phonon Φ φ_LQG
    """

    HBAR = 1.054571817e-34
    C = 2.99792458e8
    G = 6.67430e-11
    SSQ = 0.57
    PHI_0 = 1.25e12
    GAMMA_IMMIRZI = 0.2375
    L_PLANCK = 1.616255e-35

    def compute(self, dataset: dict = None) -> dict:
        import math
        ds = dataset or {}
        phi_0 = ds.get('phi_LQG_0', 1.0e-30)        # Field amplitude (Planck units)
        m_phi = ds.get('m_phi', 1.0e-33)              # Field mass (kg)
        lambda_spin = ds.get('lambda_spin', 1.0e-10)  # Spin-foam coupling
        g_phonon = ds.get('g_phonon', 1.0e-5)         # Phonon-buoyancy coupling
        n_vertices = ds.get('n_vertices', 100)         # Spin-foam vertices
        j_avg = ds.get('j_avg', 0.5)                  # Average spin at vertices
        omega = ds.get('omega', self.PHI_0)
        Gamma = ds.get('Gamma', 0.05e12)
        ssq = ds.get('SSq', self.SSQ)

        S_26_3 = (1.0 - ssq) ** 3

        # Phonon spectral function
        Phi = (1.0 / math.pi) * (Gamma / ((omega - self.PHI_0) ** 2 + Gamma ** 2))

        # Potential V(φ) = ½ m² c² φ² / ℏ² (harmonic)
        omega_phi = m_phi * self.C ** 2 / self.HBAR
        V_phi = 0.5 * omega_phi ** 2 * phi_0 ** 2
        V_prime = omega_phi ** 2 * phi_0

        # Spin-foam vertex amplitude (simplified Ponzano-Regge/EPRL)
        # A_v ~ (2j+1) × 6j-symbol approximation
        A_v_single = (2.0 * j_avg + 1.0) * math.exp(-self.GAMMA_IMMIRZI * j_avg)
        A_v_sum = n_vertices * A_v_single

        # Source term from spin-foam coupling
        source_spinfoam = lambda_spin * A_v_sum

        # Source term from phonon-buoyancy coupling
        source_phonon = 2.0 * g_phonon * Phi * phi_0

        # Total source = λ_spin Σ A_v + 2 g_phonon Φ φ
        total_source = source_spinfoam + source_phonon

        # Equation of motion residual: □φ + V'(φ) - source = 0
        # At equilibrium (static, homogeneous): V'(φ₀) = source
        eom_residual = V_prime - total_source
        eom_balance = abs(eom_residual) / max(abs(V_prime), 1e-300)

        # Lagrangian density at φ₀ (kinetic term = 0 for static)
        L_potential = -V_phi
        L_spinfoam = lambda_spin * A_v_sum * phi_0
        L_phonon = g_phonon * Phi * phi_0 ** 2
        L_total = L_potential + L_spinfoam + L_phonon

        # Effective mass correction from phonon coupling
        m_eff_sq = omega_phi ** 2 - 2.0 * g_phonon * Phi
        m_eff = math.sqrt(abs(m_eff_sq)) * self.HBAR / self.C ** 2
        tachyonic = m_eff_sq < 0

        results = []
        results.append(f"L_LQG = ½(∂φ)² - V(φ) + λ_spin Σ_v A_v·φ + g_phonon Φ φ²")
        results.append(f"EoM: □φ + V'(φ) = λ_spin Σ_v A_v + 2g_phonon Φ φ")
        results.append(f"φ₀ = {phi_0:.6e}")
        results.append(f"ω_φ = m_φ c²/ℏ = {omega_phi:.6e} rad/s")
        results.append(f"V(φ₀) = {V_phi:.6e}")
        results.append(f"V'(φ₀) = {V_prime:.6e}")
        results.append(f"A_v(j={j_avg}) = {A_v_single:.6e}, N_v = {n_vertices}")
        results.append(f"λ_spin Σ A_v = {source_spinfoam:.6e}")
        results.append(f"2 g_phonon Φ φ₀ = {source_phonon:.6e}")
        results.append(f"EoM residual = {eom_residual:.6e}")
        results.append(f"Balance ratio |residual/V'| = {eom_balance:.6e}")
        results.append(f"L_total = {L_total:.6e}")
        results.append(f"m_eff = {m_eff:.6e} kg {'(TACHYONIC)' if tachyonic else ''}")
        results.append(f"S_26^(3) = {S_26_3:.6e}")

        return {
            'L_total': L_total,
            'V_phi': V_phi,
            'source_spinfoam': source_spinfoam,
            'source_phonon': source_phonon,
            'eom_residual': eom_residual,
            'eom_balance': eom_balance,
            'm_eff': m_eff,
            'tachyonic': tachyonic,
            'A_v_sum': A_v_sum,
            'primary_equations': results,
        }


class UnifiedFUBiSMBHMergerDynamicsCalculator:
    """Unified F_U_Bi for SMBH binary merger dynamics (inside + outside).

    Computes the complete buoyancy-mediated unified field contribution
    during SMBH binary mergers, distinguishing inside-to-outside and
    outside-to-inside dynamics:

    Inside-to-outside (inspiral → ringdown):
        F_{U,Bi}^{in→out} = Σ_{i=1}^{26} [Ug1_i + Ug2_i + Ug3_i + Ug4_i]
                            × (1 + Φ/F_U · ΔE_gw) × (M_chirp/M_total)^{5/3}

    Outside-to-inside (environment → binary):
        F_{U,Bi}^{out→in} = F_U - F_{U,Bi}^{in→out}
                            - Σ_k F_{tidal,k}(r_k, M_k)

    where ΔE_gw is the gravitational wave energy release and M_chirp is
    the chirp mass of the binary.
    """

    G = 6.67430e-11
    C = 2.99792458e8
    HBAR = 1.054571817e-34
    SOLAR_MASS = 1.989e30
    SSQ = 0.57
    PHI_0 = 1.25e12
    KAPPA = 0.0005

    def compute(self, dataset: dict = None) -> dict:
        import math
        ds = dataset or {}
        M1 = ds.get('M1', 3.6e9 * self.SOLAR_MASS)    # Primary SMBH mass (kg)
        M2 = ds.get('M2', 1.2e9 * self.SOLAR_MASS)    # Secondary SMBH mass (kg)
        r_sep = ds.get('r_separation', 1.0e13)          # Binary separation (m)
        F_U = ds.get('F_U', 1.0e-8)                     # Unified field strength (N)
        n_layers = ds.get('n_layers', 26)
        Gamma = ds.get('Gamma', 0.05e12)
        omega = ds.get('omega', self.PHI_0)
        n_tidal = ds.get('n_tidal_sources', 0)
        tidal_forces = ds.get('tidal_forces', [])        # list of F_tidal values

        M_total = M1 + M2
        eta = (M1 * M2) / (M_total ** 2)                # Symmetric mass ratio
        M_chirp = M_total * eta ** 0.6                   # Chirp mass

        # Gravitational wave energy release (Newtonian inspiral)
        E_gw = self.G * M1 * M2 / (2.0 * r_sep)

        # 26-layer gravity sum (Ug1 through Ug4 per layer)
        ssq = ds.get('SSq', self.SSQ)
        Ug_sum = 0.0
        for i in range(1, n_layers + 1):
            Q_i = 1.0 / (1.0 + self.KAPPA * i)
            Ug1_i = self.G * M_total / (r_sep ** 2) * Q_i
            Ug2_i = (1.0 - ssq) * Q_i * self.G * M_total / (r_sep ** 2) * 0.1
            Ug3_i = ssq * Q_i * self.G * M_total / (r_sep ** 2) * 0.05 * math.sin(i)
            Ug4_i = self.KAPPA * Q_i * 1.0e-15
            Ug_sum += (Ug1_i + Ug2_i + Ug3_i + Ug4_i)

        # Phonon spectral function
        Phi = (1.0 / math.pi) * (Gamma / ((omega - self.PHI_0) ** 2 + Gamma ** 2))

        # Chirp mass factor
        chirp_factor = (M_chirp / M_total) ** (5.0 / 3.0)

        # GW phonon coupling
        gw_phonon = 1.0 + (Phi / F_U) * E_gw if F_U != 0 else 1.0

        # Inside-to-outside F_U_Bi
        F_UBi_in_out = Ug_sum * gw_phonon * chirp_factor

        # Tidal forces from environment
        F_tidal_total = sum(tidal_forces[:n_tidal]) if tidal_forces else 0.0

        # Outside-to-inside F_U_Bi
        F_UBi_out_in = F_U - F_UBi_in_out - F_tidal_total

        # Net buoyancy imbalance
        Delta_F = F_UBi_in_out - F_UBi_out_in

        # Merger timescale (Peters formula simplified)
        a0 = r_sep
        t_merger = (5.0 / 256.0) * (self.C ** 5 * a0 ** 4) / (
            self.G ** 3 * M1 * M2 * M_total) if M1 * M2 * M_total > 0 else float('inf')

        results = []
        results.append(f"F_U,Bi^{{in→out}} = Σ_26 [Ug1+Ug2+Ug3+Ug4] × (1+Φ/F_U·ΔE_gw) × (M_c/M)^{{5/3}}")
        results.append(f"F_U,Bi^{{out→in}} = F_U - F_U,Bi^{{in→out}} - Σ F_tidal")
        results.append(f"M₁ = {M1:.3e} kg, M₂ = {M2:.3e} kg")
        results.append(f"M_total = {M_total:.3e} kg")
        results.append(f"M_chirp = {M_chirp:.3e} kg, η = {eta:.6f}")
        results.append(f"(M_chirp/M_total)^(5/3) = {chirp_factor:.6f}")
        results.append(f"E_gw (Newtonian) = {E_gw:.6e} J")
        results.append(f"Σ_26 Ug = {Ug_sum:.6e} N/kg")
        results.append(f"GW-phonon coupling = {gw_phonon:.6e}")
        results.append(f"F_U,Bi^{{in→out}} = {F_UBi_in_out:.6e}")
        results.append(f"F_U,Bi^{{out→in}} = {F_UBi_out_in:.6e}")
        results.append(f"ΔF(in-out) = {Delta_F:.6e}")
        results.append(f"Σ F_tidal = {F_tidal_total:.6e}")
        results.append(f"t_merger (Peters) = {t_merger:.6e} s")

        return {
            'F_UBi_in_out': F_UBi_in_out,
            'F_UBi_out_in': F_UBi_out_in,
            'Delta_F': Delta_F,
            'M_chirp': M_chirp,
            'eta': eta,
            'chirp_factor': chirp_factor,
            'E_gw': E_gw,
            'Ug_sum': Ug_sum,
            'gw_phonon': gw_phonon,
            't_merger': t_merger,
            'primary_equations': results,
        }


class HydrogenUniverseDual3DMUGECalculator:
    """Hydrogen atom + full universe dual 3D MUGE system.

    Constructs a dual 3D modified-universal-gravity-equation (MUGE) system
    comparing the hydrogen atom (quantum scale) to the observable universe
    (cosmological scale), demonstrating UQFF scale invariance:

    Hydrogen MUGE:
        g_MUGE_H(r) = g_Newton + g_expansion + g_super + g_envelope
                     + g_Ug_sum + g_cosm + g_quantum + g_fluid + g_perturb
                     (9-term compressed MUGE at Bohr radius)

    Universe MUGE:
        g_MUGE_U(r) = same 9 terms evaluated at Hubble radius R_H

    3D spatial grid with MUGE evaluated at each point for volumetric
    field comparison.
    """

    G = 6.67430e-11
    C = 2.99792458e8
    HBAR = 1.054571817e-34
    KB = 1.380649e-23
    M_PROTON = 1.67262192e-27
    M_ELECTRON = 9.1093837e-31
    A_BOHR = 5.29177210903e-11     # Bohr radius (m)
    R_HUBBLE = 4.4e26              # Hubble radius (m)
    M_UNIVERSE = 1.5e53            # Observable universe mass (kg)
    H0 = 2.2e-18                   # Hubble constant (1/s)
    SSQ = 0.57
    KAPPA = 0.0005
    LAMBDA_COSMO = 1.1056e-52      # Cosmological constant (m⁻²)

    def compute(self, dataset: dict = None) -> dict:
        import math
        ds = dataset or {}
        n_grid = ds.get('n_grid', 5)   # Grid points per axis

        def muge_9term(M, r, scale_label):
            """9-term compressed MUGE at radius r for mass M."""
            if r == 0:
                r = 1e-50
            # 1. Newtonian
            g_newton = self.G * M / (r ** 2)
            # 2. Hubble expansion
            g_expansion = self.H0 ** 2 * r
            # 3. Super (magnetic suppression)
            B_ref = 1e-6               # Reference B-field (T)
            mu_0 = 4e-7 * math.pi
            g_super = (B_ref ** 2) / (2.0 * mu_0 * M / (4.0 / 3.0 * math.pi * r ** 3)) / r if M > 0 else 0
            # 4. Envelope
            g_envelope = self.KAPPA * g_newton * math.exp(-r / (100.0 * r))
            # 5. Ug_sum (26-layer aggregate)
            Ug_sum = 0.0
            for i in range(1, 27):
                Q_i = 1.0 / (1.0 + self.KAPPA * i)
                Ug_sum += Q_i * g_newton / 26.0
            # 6. Cosmological
            g_cosm = (self.LAMBDA_COSMO * self.C ** 2 / 3.0) * r
            # 7. Quantum
            g_quantum = self.HBAR / (M * r ** 3) if M > 0 else 0
            # 8. Fluid (Navier-Stokes viscous)
            eta_visc = 1e-6            # Kinematic viscosity proxy
            g_fluid = eta_visc / (r ** 2) if r > 0 else 0
            # 9. Perturbation (dark matter)
            rho_dm = 2.25e-27          # DM density kg/m³
            g_perturb = (4.0 / 3.0) * math.pi * self.G * rho_dm * r

            g_total = (g_newton + g_expansion + g_super + g_envelope
                       + Ug_sum + g_cosm + g_quantum + g_fluid + g_perturb)
            return {
                'label': scale_label,
                'r': r,
                'M': M,
                'g_newton': g_newton,
                'g_expansion': g_expansion,
                'g_super': g_super,
                'g_envelope': g_envelope,
                'g_Ug_sum': Ug_sum,
                'g_cosm': g_cosm,
                'g_quantum': g_quantum,
                'g_fluid': g_fluid,
                'g_perturb': g_perturb,
                'g_total': g_total,
            }

        # Single-point evaluations
        H_result = muge_9term(self.M_PROTON, self.A_BOHR, "Hydrogen")
        U_result = muge_9term(self.M_UNIVERSE, self.R_HUBBLE, "Universe")

        # Scale invariance ratio
        ratio = H_result['g_total'] / U_result['g_total'] if U_result['g_total'] != 0 else float('inf')

        # 3D grid evaluations (centered on each scale)
        H_grid = []
        U_grid = []
        for ix in range(n_grid):
            for iy in range(n_grid):
                for iz in range(n_grid):
                    # Fractional offset [-0.5, +0.5] per axis
                    fx = (ix / max(n_grid - 1, 1)) - 0.5
                    fy = (iy / max(n_grid - 1, 1)) - 0.5
                    fz = (iz / max(n_grid - 1, 1)) - 0.5
                    r_frac = math.sqrt(fx ** 2 + fy ** 2 + fz ** 2)
                    r_frac = max(r_frac, 0.01)

                    r_H = self.A_BOHR * (0.5 + r_frac)
                    r_U = self.R_HUBBLE * (0.5 + r_frac)

                    H_grid.append(muge_9term(self.M_PROTON, r_H, f"H({ix},{iy},{iz})")['g_total'])
                    U_grid.append(muge_9term(self.M_UNIVERSE, r_U, f"U({ix},{iy},{iz})")['g_total'])

        H_grid_avg = sum(H_grid) / len(H_grid) if H_grid else 0
        U_grid_avg = sum(U_grid) / len(U_grid) if U_grid else 0
        grid_ratio = H_grid_avg / U_grid_avg if U_grid_avg != 0 else float('inf')

        results = []
        results.append(f"g_MUGE = g_N + g_exp + g_super + g_env + g_Ug + g_cosm + g_Q + g_fluid + g_pert")
        results.append(f"--- Hydrogen (r = a₀ = {self.A_BOHR:.3e} m, M = m_p) ---")
        for k in ['g_newton', 'g_expansion', 'g_super', 'g_envelope',
                   'g_Ug_sum', 'g_cosm', 'g_quantum', 'g_fluid', 'g_perturb']:
            results.append(f"  {k} = {H_result[k]:.6e} m/s²")
        results.append(f"  g_MUGE_H = {H_result['g_total']:.6e} m/s²")
        results.append(f"--- Universe (r = R_H = {self.R_HUBBLE:.3e} m, M = M_U) ---")
        for k in ['g_newton', 'g_expansion', 'g_super', 'g_envelope',
                   'g_Ug_sum', 'g_cosm', 'g_quantum', 'g_fluid', 'g_perturb']:
            results.append(f"  {k} = {U_result[k]:.6e} m/s²")
        results.append(f"  g_MUGE_U = {U_result['g_total']:.6e} m/s²")
        results.append(f"Scale ratio g_H / g_U = {ratio:.6e}")
        results.append(f"3D grid ({n_grid}³ = {n_grid**3} pts): avg(g_H) = {H_grid_avg:.6e}, avg(g_U) = {U_grid_avg:.6e}")
        results.append(f"3D grid ratio = {grid_ratio:.6e}")

        return {
            'H_result': H_result,
            'U_result': U_result,
            'scale_ratio': ratio,
            'H_grid_avg': H_grid_avg,
            'U_grid_avg': U_grid_avg,
            'grid_ratio': grid_ratio,
            'n_grid_total': n_grid ** 3,
            'primary_equations': results,
        }


# Registry for Session 225 Phase 4 -- LQG / Holonomy / T2 / SMBH / 3D MUGE
SESSION_225_PHASE4_LQG_HOLONOMY_SMBH_3DMUGE_CALCULATORS = {
    'SCmLQGAreaOperatorDerivationCalculator':            SCmLQGAreaOperatorDerivationCalculator(),
    'SCmQubitT2CoherenceFUBiRatioCalculator':            SCmQubitT2CoherenceFUBiRatioCalculator(),
    'PhononModulatedHolonomySCmCalculator':              PhononModulatedHolonomySCmCalculator(),
    'LQGBuoyancySectorLagrangianVariationCalculator':    LQGBuoyancySectorLagrangianVariationCalculator(),
    'UnifiedFUBiSMBHMergerDynamicsCalculator':           UnifiedFUBiSMBHMergerDynamicsCalculator(),
    'HydrogenUniverseDual3DMUGECalculator':              HydrogenUniverseDual3DMUGECalculator(),
}
'''

# ── Validate snippet AST ──
print("AST-validating new code snippet...")
try:
    ast.parse(NEW_CODE)
    print("  ✓ AST validation passed")
except SyntaxError as e:
    print(f"  ✗ SyntaxError: {e}")
    sys.exit(1)

# ── Smudge CP1 ──
print(f"Smudging {CP1} via git lfs...")
with open(CP1, "rb") as f:
    result = subprocess.run(["git", "lfs", "smudge"], stdin=f, capture_output=True)
if result.returncode != 0:
    print(f"  git lfs smudge failed: {result.stderr.decode()}")
    sys.exit(1)

data = result.stdout
if data[:2] == b"\xff\xfe":
    text = data.decode("utf-16-le")
    print("  Decoded as UTF-16-LE")
elif data[:3] == b"\xef\xbb\xbf":
    text = data[3:].decode("utf-8")
    print("  Decoded as UTF-8 with BOM")
else:
    text = data.decode("utf-8")
    print(f"  Decoded as UTF-8, {len(text)} chars")

lines_before = text.count("\n")
print(f"  Lines before: {lines_before + 1}")

# ── Append ──
text += NEW_CODE
lines_after = text.count("\n")
print(f"  Lines after: {lines_after + 1}")
print(f"  Appended {lines_after - lines_before} lines")

# ── Write back ──
with open(CP1, "w", encoding="utf-8") as f:
    f.write(text)
print(f"  ✓ Written back as UTF-8")

# ── Verify ──
with open(CP1, "r", encoding="utf-8") as f:
    final = f.read()
print(f"  Final size: {len(final)} chars, {final.count(chr(10)) + 1} lines")
# Check last 3 lines
last3 = final.strip().split("\n")[-3:]
for ln in last3:
    print(f"  | {ln}")
print("\nDone.")
