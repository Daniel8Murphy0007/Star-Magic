#!/usr/bin/env python3
"""Session 225 Phase 5: Append 3 remaining gap calculators to CondensedPhysics.py.

Gaps from workflow report analysis:
  1. SCmStringTheory26DActionCalculator — 26D string action with phonon Lagrangian + brane coupling
  2. UQFF26DGeometricFoldingOperatorCalculator — F_26(x) folding operator + folded metric
  3. VDSDVPBHUnifiedNumberSystemCalculator — unified VDS/DVP/BH derivation engine with Ramanujan acceleration
"""
import ast, subprocess, sys

CP1 = "CondensedPhysics.py"

NEW_CODE = r'''

# ═══════════════════════════════════════════════════════════════════════════════
# SESSION 225 PHASE 5 — SCm String Theory / 26D Folding / VDS-DVP-BH Unified
# ═══════════════════════════════════════════════════════════════════════════════


class SCmStringTheory26DActionCalculator:
    """SCm phonon coupling to strings and branes in 26D compactification.

    Constructs the complete SCm-String action integrating the 26D Einstein-
    Hilbert-Yang-Mills action with UQFF buoyancy and phonon Lagrangian terms:

        S_SCm-String = ∫d²⁶x √(-g) [ R - ¼F^a_μν F_a^μν
                       + ½η ρ_A v_UA² cos(πt_n)
                       + L_phonon(Φ, ω, Γ) ]

    The phonon-modulated string tension:
        T_SCm = T_0 · S_26^(3)([SSq]) · Φ_1.25THz(ω, Γ)

    The brane phonon coupling (p-brane):
        δS_brane = ∫d^{p+1}σ √(-γ_ab) · Φ · E_net(t, Γ)

    Evaluates all terms numerically for a given compactification geometry.
    """

    G = 6.67430e-11
    HBAR = 1.054571817e-34
    C = 2.99792458e8
    ALPHA_PRIME = 1e-34              # String scale α' (m²)
    SSQ = 0.57
    PHI_0 = 1.25e12                  # Phonon peak frequency (Hz)
    KAPPA = 0.0005

    def compute(self, dataset: dict = None) -> dict:
        import math
        ds = dataset or {}

        # Compactification parameters
        n_dim = ds.get('n_dim', 26)
        R_compact = ds.get('R_compact', 1.0e-35)     # Compactification radius (m)
        alpha_p = ds.get('alpha_prime', self.ALPHA_PRIME)
        ssq = ds.get('SSq', self.SSQ)
        omega = ds.get('omega', self.PHI_0)
        Gamma = ds.get('Gamma', 0.1e12)
        rho_A = ds.get('rho_A', 1.0e-26)             # Aether density (kg/m³)
        v_UA = ds.get('v_UA', 1.0e-4)                 # [UA] velocity (m/s)
        t_n = ds.get('t_n', 1.0)                       # Negative-time cycle index
        eta = ds.get('eta_coupling', 0.5)              # Buoyancy coupling
        p_brane = ds.get('p_brane', 3)                 # Brane dimension (3 = D3-brane)
        rho_phonon = ds.get('rho_phonon', 1.0e-20)    # Phonon energy density (J/m³)

        S_26_3 = (1.0 - ssq) ** 3

        # Gaussian phonon spectral function (as per report)
        omega_SCm = 2.0 * math.pi * self.PHI_0
        omega_in = 2.0 * math.pi * omega
        Gamma_rad = 2.0 * math.pi * Gamma
        Phi_gauss = math.exp(-((omega_in - omega_SCm) ** 2) / (2.0 * Gamma_rad ** 2))

        # --- String tension ---
        T_0 = 1.0 / (2.0 * math.pi * alpha_p)         # Fundamental string tension
        T_SCm = T_0 * S_26_3 * Phi_gauss

        # --- 26D action components ---
        # Volume factor of compact dimensions (torus compactification)
        V_compact = (2.0 * math.pi * R_compact) ** (n_dim - 4)

        # Ricci scalar contribution (Einstein-Hilbert in 26D → 4D)
        R_EH = 1.0 / (16.0 * math.pi * self.G) * V_compact

        # Yang-Mills term: -¼ F²
        g_YM = ds.get('g_YM', 0.5)                    # Gauge coupling
        F_sq = g_YM ** 2 / (4.0 * math.pi)            # Proxy for F² in natural units
        S_YM = -0.25 * F_sq * V_compact

        # Buoyancy term: ½ η ρ_A v_UA² cos(πt_n)
        S_buoy = 0.5 * eta * rho_A * v_UA ** 2 * math.cos(math.pi * t_n) * V_compact

        # Phonon Lagrangian: L_phonon = ½ (∂Φ)² - ½ m_phonon² Φ²
        m_phonon = self.HBAR * omega_SCm / self.C ** 2
        L_phonon = 0.5 * (Phi_gauss ** 2) * (omega_SCm ** 2 - m_phonon ** 2 * self.C ** 4 / self.HBAR ** 2)
        S_phonon = L_phonon * V_compact

        # Total action density (per unit 4-volume)
        S_total = R_EH + S_YM + S_buoy + S_phonon

        # --- Brane phonon coupling ---
        E_net = rho_phonon * S_26_3
        sigma_brane = (2.0 * math.pi * R_compact) ** (p_brane - 3) if p_brane > 3 else 1.0
        delta_S_brane = sigma_brane * Phi_gauss * E_net

        # --- String-phonon mass spectrum (Regge trajectory) ---
        n_max = ds.get('n_excitations', 5)
        mass_spectrum = []
        for n in range(n_max + 1):
            M_n_sq = n / alpha_p - 1.0 / alpha_p  # Open string: M² = (n-1)/α'
            M_n = math.sqrt(max(M_n_sq, 0.0))
            M_n_SCm = M_n * S_26_3 * Phi_gauss if M_n > 0 else 0.0
            mass_spectrum.append({'n': n, 'M_n': M_n, 'M_n_SCm': M_n_SCm})

        results = []
        results.append(f"S_SCm-String = ∫d²⁶x √(-g) [R - ¼F² + ½ηρ_Av_UA²cos(πt_n) + L_phonon]")
        results.append(f"n_dim = {n_dim}, R_compact = {R_compact:.3e} m")
        results.append(f"V_compact = (2πR)^(d-4) = {V_compact:.6e}")
        results.append(f"T₀ = 1/(2πα') = {T_0:.6e} N")
        results.append(f"T_SCm = T₀·S_26^(3)·Φ_gauss = {T_SCm:.6e} N")
        results.append(f"S_26^(3)([SSq]={ssq}) = {S_26_3:.6e}")
        results.append(f"Φ_gauss(ω={omega:.3e}, Γ={Gamma:.3e}) = {Phi_gauss:.6f}")
        results.append(f"S_EH (Ricci) = {R_EH:.6e}")
        results.append(f"S_YM = {S_YM:.6e}")
        results.append(f"S_buoy = {S_buoy:.6e}")
        results.append(f"S_phonon = {S_phonon:.6e}")
        results.append(f"S_total (action density) = {S_total:.6e}")
        results.append(f"δS_brane (p={p_brane}) = {delta_S_brane:.6e}")
        results.append(f"E_net = ρ_phonon·S_26^(3) = {E_net:.6e} J/m³")
        for ms in mass_spectrum:
            results.append(f"  String level n={ms['n']}: M_n={ms['M_n']:.6e} kg, M_n_SCm={ms['M_n_SCm']:.6e} kg")

        return {
            'T_0': T_0,
            'T_SCm': T_SCm,
            'S_26_3': S_26_3,
            'Phi_gauss': Phi_gauss,
            'R_EH': R_EH,
            'S_YM': S_YM,
            'S_buoy': S_buoy,
            'S_phonon': S_phonon,
            'S_total': S_total,
            'delta_S_brane': delta_S_brane,
            'V_compact': V_compact,
            'mass_spectrum': mass_spectrum,
            'primary_equations': results,
        }


class UQFF26DGeometricFoldingOperatorCalculator:
    """UQFF 26D geometric folding operator and folded spacetime metric.

    Derives the 26D folding operator that maps fields from the full UQFF
    26-dimensional space to the effective 4D spacetime, paralleling Wolfram
    hypergraph rewriting:

        F_26(x) = x · (26!)^{-1/13} · S_26^(3)([SSq]) · Φ_1.25THz(ω, Γ)

    The folded metric is:
        ds²_folded = g_μν · F_26(r)

    where the folding operator acts as a conformal factor on the 4D slice
    of the full 26D metric.  Each of the 26 quantum layers is mapped to
    a node in the Wolfram-parallel hypergraph, with phonon resonances
    as edges.
    """

    HBAR = 1.054571817e-34
    C = 2.99792458e8
    G = 6.67430e-11
    SSQ = 0.57
    PHI_0 = 1.25e12
    KAPPA = 0.0005

    def compute(self, dataset: dict = None) -> dict:
        import math
        ds = dataset or {}
        x_field = ds.get('x_field', 1.0)              # Field amplitude/value
        r = ds.get('r', 1.0e10)                        # Radial coordinate (m)
        M = ds.get('M', 1.989e30)                      # Central mass (kg)
        omega = ds.get('omega', self.PHI_0)
        Gamma = ds.get('Gamma', 0.1e12)
        ssq = ds.get('SSq', self.SSQ)
        n_layers = ds.get('n_layers', 26)

        S_26_3 = (1.0 - ssq) ** 3

        # Gaussian phonon spectral function
        omega_SCm = 2.0 * math.pi * self.PHI_0
        omega_in = 2.0 * math.pi * omega
        Gamma_rad = 2.0 * math.pi * Gamma
        Phi_gauss = math.exp(-((omega_in - omega_SCm) ** 2) / (2.0 * Gamma_rad ** 2))

        # 26! factorial
        fact_26 = math.factorial(n_layers)

        # Folding normalization: (26!)^{-1/13}
        fold_norm = fact_26 ** (-1.0 / 13.0)

        # Folding operator: F_26(x)
        F_26 = x_field * fold_norm * S_26_3 * Phi_gauss

        # Folded metric components (Schwarzschild base)
        r_s = 2.0 * self.G * M / (self.C ** 2)        # Schwarzschild radius
        g_tt = -(1.0 - r_s / r) if r > r_s else -(1e-10)
        g_rr = 1.0 / (1.0 - r_s / r) if r > r_s else 1e10
        g_theta = r ** 2
        g_phi = r ** 2  # (sin²θ = 1 at equator)

        # Folded metric: ds²_folded = g_μν · F_26(r)
        F_26_r = r * fold_norm * S_26_3 * Phi_gauss    # F_26 evaluated at r
        g_tt_folded = g_tt * F_26_r
        g_rr_folded = g_rr * F_26_r
        g_theta_folded = g_theta * F_26_r
        g_phi_folded = g_phi * F_26_r

        # Effective Schwarzschild radius under folding
        r_s_folded = r_s * fold_norm * S_26_3 * Phi_gauss

        # Wolfram-parallel hypergraph properties
        # 26 nodes (quantum layers), edges = phonon resonance connections
        # Each layer i connects to layers i±1 and i±13 (modular topology)
        n_nodes = n_layers
        n_edges_nearest = n_layers                     # Ring topology (i → i+1)
        n_edges_skip = n_layers                        # Skip connections (i → i+13)
        n_edges_total = n_edges_nearest + n_edges_skip
        graph_density = n_edges_total / (n_nodes * (n_nodes - 1) / 2)

        # Layer-by-layer folding factors
        layer_factors = []
        for i in range(1, n_layers + 1):
            Q_i = 1.0 / (1.0 + self.KAPPA * i)
            F_i = Q_i * fold_norm * S_26_3 * Phi_gauss
            layer_factors.append({'layer': i, 'Q_i': Q_i, 'F_i': F_i})

        # Aggregate folding (product of all layer factors)
        F_product = 1.0
        for lf in layer_factors:
            F_product *= (1.0 + lf['F_i'])

        results = []
        results.append(f"F_26(x) = x · (26!)^{{-1/13}} · S_26^(3) · Φ_gauss")
        results.append(f"ds²_folded = g_μν · F_26(r)")
        results.append(f"26! = {fact_26:.6e}")
        results.append(f"(26!)^(-1/13) = {fold_norm:.6e}")
        results.append(f"S_26^(3)([SSq]={ssq}) = {S_26_3:.6e}")
        results.append(f"Φ_gauss(ω={omega:.3e}, Γ={Gamma:.3e}) = {Phi_gauss:.6f}")
        results.append(f"F_26(x={x_field}) = {F_26:.6e}")
        results.append(f"F_26(r={r:.3e}) = {F_26_r:.6e}")
        results.append(f"r_s = {r_s:.6e} m")
        results.append(f"r_s_folded = {r_s_folded:.6e} m")
        results.append(f"g_tt = {g_tt:.6e}, g_tt_folded = {g_tt_folded:.6e}")
        results.append(f"g_rr = {g_rr:.6e}, g_rr_folded = {g_rr_folded:.6e}")
        results.append(f"Wolfram hypergraph: {n_nodes} nodes, {n_edges_total} edges, density = {graph_density:.4f}")
        results.append(f"Layer product Π(1+F_i) = {F_product:.6e}")

        return {
            'F_26': F_26,
            'F_26_r': F_26_r,
            'fold_norm': fold_norm,
            'S_26_3': S_26_3,
            'Phi_gauss': Phi_gauss,
            'r_s': r_s,
            'r_s_folded': r_s_folded,
            'g_tt_folded': g_tt_folded,
            'g_rr_folded': g_rr_folded,
            'graph_density': graph_density,
            'F_product': F_product,
            'layer_factors': layer_factors,
            'primary_equations': results,
        }


class VDSDVPBHUnifiedNumberSystemCalculator:
    """Unified VDS/DVP/BH number system derivation engine.

    Computes all three UQFF number systems simultaneously with Ramanujan
    acceleration, cross-validating internal consistency:

    VDS (Vacuum Density Series):
        VDS = Σ_{n=1}^N [SSq]^n / n^26 · R_n^{(26,3)} ≈ Li_26([SSq])

    DVP (Dipole Vortex Primes):
        a(p) = [SSq]^{π(p)} / p^26,  p prime, p > 26

    BH (Buoyancy Harmonics):
        Ug2 = Σ_{m=1}^M H_m · (1 - exp(-[SSq]·m)) · cos(ω_Ug2 · t_n)

    The Ramanujan acceleration factor for VDS:
        R_n^{(26,3)} = (2π)^{n/6}/n! · [1 + Σ_{m=1}^3 1/n^{26m}
                       · Σ_{j=1}^{26} (-1)^{j+1} C(26,j) (26-j)!/n^j ]

    Cross-validation: VDS provides vacuum density, DVP encodes vortex
    quantization, BH drives Ug2 charge-reactivity gravity.
    """

    SSQ = 0.57
    PHI_0 = 1.25e12
    KAPPA = 0.0005

    def compute(self, dataset: dict = None) -> dict:
        import math
        ds = dataset or {}
        ssq = ds.get('SSq', self.SSQ)
        N_vds = ds.get('N_vds', 50)                    # VDS truncation order
        N_bh = ds.get('N_bh', 20)                       # BH harmonic order
        omega_Ug2 = ds.get('omega_Ug2', 1.0e6)        # Ug2 frequency (Hz)
        t_n = ds.get('t_n', 1.0)                        # Negative-time cycle

        # ─── VDS: Vacuum Density Series with Ramanujan acceleration ───
        def ramanujan_R(n, s=26, order=3):
            """Ramanujan acceleration factor R_n^{(s, order)}."""
            if n == 0:
                return 1.0
            base = ((2.0 * math.pi) ** (n / 6.0)) / math.factorial(min(n, 170))
            correction = 0.0
            for m in range(1, order + 1):
                inner = 0.0
                for j in range(1, min(s + 1, 20)):
                    sign = (-1.0) ** (j + 1)
                    binom = math.comb(s, j)
                    fact_term = math.factorial(s - j) if (s - j) <= 170 else 1e300
                    inner += sign * binom * fact_term / (n ** j)
                correction += inner / (n ** (s * m))
            return base * (1.0 + correction)

        vds_terms = []
        vds_sum = 0.0
        vds_plain = 0.0  # Without Ramanujan acceleration (plain Li_26)
        for n in range(1, N_vds + 1):
            plain_term = (ssq ** n) / (n ** 26)
            R_n = ramanujan_R(n)
            accel_term = plain_term * R_n
            vds_plain += plain_term
            vds_sum += accel_term
            if n <= 5 or n == N_vds:
                vds_terms.append({'n': n, 'plain': plain_term, 'R_n': R_n, 'accel': accel_term})

        # ─── DVP: Dipole Vortex Primes ───
        def is_prime(num):
            if num < 2:
                return False
            if num < 4:
                return True
            if num % 2 == 0 or num % 3 == 0:
                return False
            i = 5
            while i * i <= num:
                if num % i == 0 or num % (i + 2) == 0:
                    return False
                i += 6
            return True

        def prime_counting(p):
            """π(p) — count primes ≤ p."""
            count = 0
            for i in range(2, p + 1):
                if is_prime(i):
                    count += 1
            return count

        primes_above_26 = [p for p in range(27, 200) if is_prime(p)]
        dvp_values = []
        for p in primes_above_26[:15]:
            pi_p = prime_counting(p)
            a_p = (ssq ** pi_p) / (p ** 26)
            dvp_values.append({'p': p, 'pi_p': pi_p, 'a_p': a_p})

        # ─── BH: Buoyancy Harmonics ───
        def harmonic_H(m):
            """m-th partial harmonic sum H_m = Σ_{k=1}^m 1/k."""
            return sum(1.0 / k for k in range(1, m + 1))

        bh_terms = []
        Ug2_sum = 0.0
        for m in range(1, N_bh + 1):
            H_m = harmonic_H(m)
            saturation = 1.0 - math.exp(-ssq * m)
            oscillation = math.cos(omega_Ug2 * t_n * m)    # m-th harmonic
            term = H_m * saturation * oscillation
            Ug2_sum += term
            if m <= 5 or m == N_bh:
                bh_terms.append({'m': m, 'H_m': H_m, 'saturation': saturation, 'term': term})

        # ─── Cross-validation ───
        # VDS provides the vacuum density normalization
        # DVP encodes vortex quantization (proplyd orbital r_q)
        # BH drives Ug2 charge-reactivity in E_net
        # Consistency: VDS ≈ Li_26([SSq]) and DVP samples this at primes
        # DVP sum should be bounded by VDS
        dvp_total = sum(d['a_p'] for d in dvp_values)
        cross_ratio = dvp_total / vds_sum if vds_sum != 0 else 0.0
        # BH saturation at large m → 1 (consistent with [SSq] < 1)
        bh_saturation_limit = 1.0 - math.exp(-ssq * N_bh)

        results = []
        results.append(f"=== VDS: Vacuum Density Series (N={N_vds}) ===")
        results.append(f"VDS = Σ [SSq]^n/n^26 · R_n^(26,3) = Li_26([SSq]) accelerated")
        results.append(f"VDS_plain (no accel) = {vds_plain:.6e}")
        results.append(f"VDS_accelerated = {vds_sum:.6e}")
        results.append(f"Acceleration ratio = {vds_sum/vds_plain:.6f}" if vds_plain != 0 else "Acceleration ratio = N/A")
        for vt in vds_terms:
            results.append(f"  n={vt['n']}: plain={vt['plain']:.6e}, R_n={vt['R_n']:.6e}, accel={vt['accel']:.6e}")

        results.append(f"=== DVP: Dipole Vortex Primes (p > 26) ===")
        results.append(f"a(p) = [SSq]^π(p) / p^26")
        for dv in dvp_values[:8]:
            results.append(f"  p={dv['p']}: π(p)={dv['pi_p']}, a(p)={dv['a_p']:.6e}")
        results.append(f"DVP total (15 primes) = {dvp_total:.6e}")

        results.append(f"=== BH: Buoyancy Harmonics (M={N_bh}) ===")
        results.append(f"Ug2 = Σ H_m·(1-exp(-[SSq]·m))·cos(ω_Ug2·t_n)")
        results.append(f"ω_Ug2 = {omega_Ug2:.3e} Hz, t_n = {t_n}")
        for bt in bh_terms:
            results.append(f"  m={bt['m']}: H_m={bt['H_m']:.4f}, sat={bt['saturation']:.6f}, term={bt['term']:.6e}")
        results.append(f"Ug2_sum = {Ug2_sum:.6e}")
        results.append(f"BH saturation limit (m={N_bh}) = {bh_saturation_limit:.6f}")

        results.append(f"=== Cross-Validation ===")
        results.append(f"DVP/VDS ratio = {cross_ratio:.6e}")
        results.append(f"[SSq] = {ssq} (consistent across all three)")

        return {
            'VDS_plain': vds_plain,
            'VDS_accelerated': vds_sum,
            'DVP_values': dvp_values,
            'DVP_total': dvp_total,
            'Ug2_sum': Ug2_sum,
            'BH_terms': bh_terms,
            'cross_ratio': cross_ratio,
            'bh_saturation_limit': bh_saturation_limit,
            'primary_equations': results,
        }


# Registry for Session 225 Phase 5 -- SCm String Theory / 26D Folding / VDS-DVP-BH
SESSION_225_PHASE5_STRING_FOLDING_NUMBER_SYSTEMS_CALCULATORS = {
    'SCmStringTheory26DActionCalculator':            SCmStringTheory26DActionCalculator(),
    'UQFF26DGeometricFoldingOperatorCalculator':     UQFF26DGeometricFoldingOperatorCalculator(),
    'VDSDVPBHUnifiedNumberSystemCalculator':         VDSDVPBHUnifiedNumberSystemCalculator(),
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
last3 = final.strip().split("\n")[-3:]
for ln in last3:
    print(f"  | {ln}")
print("\nDone.")
