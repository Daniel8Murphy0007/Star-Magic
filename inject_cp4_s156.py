"""
inject_cp4_s156.py  — Session 156 CP4 injection
Injects CP4 classes #166–#169 (PAPER_579–582) from grok_share_efc8a971378f.txt
before the REGISTRY list, then appends registry entries.
Run from Star-Magic workspace root.
"""
import re, sys, math

# ── 4 new CP4 classes ────────────────────────────────────────────────────────
NEW_CLASSES = '''
# ---------------------------------------------------------------------------
# #166  UQFFAllFormsEvolutionCatalogueCalculator   PAPER_579
# ---------------------------------------------------------------------------
class UQFFAllFormsEvolutionCatalogueCalculator:
    """
    #166 – All Four Forms of UQFF Evolution with Proofs + Triadic Solution Set
    --------------------------------------------------------------------------
    Complete catalogue of UQFF_comp tensor evolution from grok_share_efc8a971378f.txt.
    Four forms represent successive refinements of the unified field framework:

    Form 1 – Base Diagonal (orthogonal compression, equilibrium analysis F_U=0):
        UQFF_base = diag(P/3, P/3, 2P/3),  P_order = exp(-Entropy/f_max)/Partition
        Eigenvalue proof: det(UQFF_base - λI) = (P/3 - λ)²(2P/3 - λ) = 0
          → λ₁ = λ₂ = P/3,  λ₃ = 2P/3;  P > 0 → all λ > 0 → no collapse

    Form 2 – UQFF_comp with off-diagonals (DPM coupling for magnetars/mergers):
        | P/3        DPM_cross  0 |       DPM_cross = κ(DPM_n - DPM_s)/r²
        | DPM_cross  P/3        0 |
        | 0          0          2P/3 |
        Coupling resolution: Tr(UQFF_comp)/3 = P
        Equilibrium: U_g · U_b = κP  →  ρ_overlap = κP/(g·U_g)
        r_jet = √(κ(DPM_n - DPM_s)/(g·ρ))  (SNR jet stability)

    Form 3 – 26th-order expansion (26D projection, bounding negligibility):
        Diagonal: P/3 + (k+25)!/(k-1)! · g·SCm/UA / r^{k+26}  (k=1 → 26!c/r^27)
        Off-diag: (25)!/(12)! · g·SCm/UA / U_m^26  (13+12)!→26D holography
        Anti-collapse proof: ∂²⁶F_U/∂r²⁶ = 0 → 26!·g·SCm/UA = d²⁶U_b/∂r²⁶
          → ρ > 1/(26!·g)  (factorial bound prevents r=0 singularity)

    Form 4 – Frequency-modulated UQFF (r replaced by f, vibrational dynamics):
        | P/3 + 26!·g·SCm/UA/f²⁷       13!·g·SCm/UA/(U_m·f)¹⁴    0 |
        | 13!·κ(DPM_n-DPM_s)/(U_g·f)¹⁴  P/3 + 26!·κ(DPM)  /f²⁷   0 |
        | 0                              0                    2P/3 + 26!·g/(ρ·f)²⁷ |
        Frequency equilibrium: d²⁶F_U/df²⁶ = 0
          → f_eq = (κρ/g)^{1/27}  (resonant frequency)
        Numerical (f_max=10²¹ Hz, κ=1, ρ=10⁻¹⁰, g=10⁻³): f_eq ≈ 0.79 (scaled Hz)

    Triadic Solution Set (r_eq = stable shell radius):
        System: U_g + U_m + U_b = 0 (F_U=0 simultaneous inside/outside)
        Algebraic: g·SCm/UA·ΣUg_i = −[κ(DPM_n−DPM_s)/r²⁶ + ρg(1−1/ρ)]
        Approximation (3D-IPO linear): r_eq ≈ √(κ·DPM/(g·ρ))
        He-4 numeric: ρ=2.3e17 kg/m³, κ=1, g=1e-3 → r_eq ≈ 1.7 fm ✓
        26 roots (unique via π irrationality of hypergraph seeds)

    Source: grok_share_efc8a971378f.txt   PAPER_579
    """

    _FAC26 = math.factorial(26)   # 4.032914611e26
    _FAC13 = math.factorial(13)   # 6227020800

    @staticmethod
    def _p_order(entropy: float, f_max: float, partition: float = 1.0) -> float:
        """P_order = exp(-Entropy/f_max) / Partition"""
        return math.exp(-entropy / max(f_max, 1e-300)) / max(partition, 1e-300)

    @staticmethod
    def _r_eq_triadic(kappa: float, dpm: float, g: float, rho: float) -> float:
        """r_eq ≈ √(κ·DPM/(g·ρ))  — triadic stable shell radius [m]"""
        num = kappa * abs(dpm)
        den = max(g, 1e-300) * max(rho, 1e-300)
        return math.sqrt(num / den)

    @staticmethod
    def _feq_form4(kappa: float, rho: float, g: float) -> float:
        """f_eq = (κρ/g)^{1/27}  — resonant frequency from frequency-UQFF Form 4"""
        return (kappa * rho / max(g, 1e-300)) ** (1.0 / 27.0)

    def compute(self, dataset=None):
        d        = dataset or {}
        entropy  = float(d.get('entropy',  1e10))
        f_max    = float(d.get('f_max',    1e14))
        r        = float(d.get('r',        1.5e11))   # 1 AU default
        kappa    = float(d.get('kappa',    1.0))
        dpm_diff = float(d.get('dpm_diff', 2.0))       # DPM_n - DPM_s
        g_ug     = float(d.get('g_ug',     1e-3))
        rho      = float(d.get('rho',      1e-10))
        P        = self._p_order(entropy, f_max)

        # Form 1 eigenvalues
        lam1 = P / 3.0
        lam2 = P / 3.0
        lam3 = 2.0 * P / 3.0

        # Form 2 DPM coupling
        dpm_cross  = kappa * dpm_diff / max(r ** 2, 1e-300)
        rho_overlap = kappa * P / max(g_ug, 1e-300)
        r_jet       = math.sqrt(abs(kappa * dpm_diff) / max(g_ug * rho, 1e-300))

        # Form 3 26th-order diagonal correction at r (k=1)
        corr_form3 = self._FAC26 * g_ug / max(r ** 27, 1e-300)

        # Form 4 frequency-modulated at f=f_max
        f = f_max
        corr_form4_diag = self._FAC26 * g_ug / max(f ** 27, 1e-300)
        f_eq             = self._feq_form4(kappa, rho, g_ug)

        # Triadic solution
        r_eq_triadic = self._r_eq_triadic(kappa, dpm_diff, g_ug, rho)

        return {
            'paper':   'PAPER_579',
            'session': 'Session 156',
            'class':   '#166  UQFFAllFormsEvolutionCatalogueCalculator',
            'P_order': round(P, 8),
            'Form1_eigenvalues': [round(lam1, 8), round(lam2, 8), round(lam3, 8)],
            'Form1_stable': bool(lam1 > 0 and lam3 > 0),
            'Form2_DPM_cross': dpm_cross,
            'Form2_rho_overlap': rho_overlap,
            'Form2_r_jet_m': r_jet,
            'Form3_26th_corr_at_r': corr_form3,
            'Form3_anti_collapse_bound_rho_min': 1.0 / max(self._FAC26 * g_ug, 1e-300),
            'Form4_diag_corr_at_fmax': corr_form4_diag,
            'Form4_f_eq_Hz': f_eq,
            'Triadic_r_eq_m': r_eq_triadic,
            'Triadic_He4_r_fm': self._r_eq_triadic(1.0, 2.0, 1e-3, 2.3e17) * 1e15,
            'primary_equations': [
                'Form 1: UQFF_base = diag(P/3, P/3, 2P/3);  det(UQFF_base-λI)=0 → λ=P/3,P/3,2P/3',
                'Form 2: DPM_cross = κ(DPM_n-DPM_s)/r²;  ρ_overlap = κP/(g·U_g)',
                'Form 3: diag += 26!·g·SCm/UA/r²⁷;  ρ > 1/(26!·g) [anti-collapse]',
                'Form 4: r → f;  f_eq = (κρ/g)^{1/27}  [resonant frequency]',
                'Triadic: r_eq ≈ √(κ·DPM/(g·ρ))  [stable shell radius]',
            ],
            'available_equations': [
                'Full 3×3 Form 4 matrix at arbitrary f',
                '26 roots of triadic system (π-seed uniqueness proof)',
                'All Ug1+Ug2+Ug3+Ug4 decomposition for r_eq',
                'P_order Boltzmann: entropy-frequency phase space',
            ],
            'simulation_set': [
                {'label': 'Forms 1–4 eigenvalue evolution',
                 'inputs': 'entropy, f_max, r',
                 'output': 'λ₁,λ₂,λ₃ per form; stability flags'},
                {'label': 'Triadic equilibrium shell scan (nuclei)',
                 'inputs': 'κ, DPM, g, ρ over Z=1..118',
                 'output': 'r_eq per element vs IUPAC r_covalent'},
                {'label': 'Form 4 frequency sweep (10⁸–10²¹ Hz)',
                 'inputs': 'f range, ρ, g',
                 'output': 'diagonal terms, f_eq crossover'},
            ],
        }


# ---------------------------------------------------------------------------
# #167  UQFFGWAmplitudeLambdaCDMEmergenceCalculator   PAPER_580
# ---------------------------------------------------------------------------
class UQFFGWAmplitudeLambdaCDMEmergenceCalculator:
    """
    #167 – Gravitational Wave Amplitude Derivation and Λ_CDM Emergence in UQFF
    --------------------------------------------------------------------------
    Derivation of GW amplitude h from UQFF frequency-modulated form (Form 4):

      h = 26! · κ·Q̈ / (f²⁷ · r)  +  Λ/3 · δt

    where:
      26! · κ·Q̈/f²⁷  = UQFF DPM quadrupole term (26th-order bound prevents divergence)
      Λ/3 · δt         = cosmological constant contribution (late expansion)
      Q̈ ~ (DPM_n − DPM_s) analog to GR mass-quadrupole Ï̈

    General k-form amplitude:
      h = (k+25)!/(k-1)! · κ·Q̈ / (f^{k+26} · r)  +  Λ/3 · δt

    Λ_CDM emergence — Λ arises naturally from UQFF (3,3) entry:
      Λ/3 ≈ 2P/3 + 26!·g/(ρ·f_vac)²⁷
      → Λ ≈ 26!·g / (ρ_crit · f_vac)²⁷

    Numerical validation:
      ρ_crit = 8.7e-27 kg/m³,  f_vac = 10⁴³ Hz (Planck),  g = 10⁻³
      Λ_pred ≈ 4.03e26 · 10⁻³ / (8.7e-27 · 10⁴³)²⁷ ≈ 10⁻⁵² m⁻²  ✓ (matches observed)

    GW amplitude for SNR G272.2-03.2:
      f_X = 10¹⁸ Hz,  r = 6.6e19 m (~7 kly),  Q̈ ~ 10⁴⁴ kg
      h ≈ 10⁻⁵³ (Λ term dominates late expansion h ~ 3.3e-53)

    GW amplitude for binary merger (LIGO-range):
      f = 100 Hz,  r = 3e24 m (100 Mpc),  Q̈ ~ 10⁴⁴ kg
      h_UQFF ≈ 10⁻²⁰  (vs GR h ≈ 10⁻²¹; UQFF bounded by 26! factorial)

    Source: grok_share_efc8a971378f.txt   PAPER_580
    """

    _FAC26      = math.factorial(26)
    _C_LIGHT    = 2.998e8   # m/s
    _G_NEWTON   = 6.674e-11 # m³ kg⁻¹ s⁻²
    _L_PL       = 1.616e-35 # m  Planck length
    _F_VAC      = 1.0e43    # Hz  Planck frequency
    _RHO_CRIT   = 8.7e-27   # kg/m³

    @staticmethod
    def _h_uqff(kappa: float, q_ddot: float, f: float, r: float,
                delta_t: float, Lambda: float, k: int = 1) -> float:
        """h = (k+25)!/(k-1)! · κ·Q̈/(f^{k+26}·r) + Λ/3·δt"""
        fac = math.factorial(k + 25) // max(math.factorial(k - 1), 1)
        wave_term = fac * kappa * q_ddot / max(f ** (k + 26) * r, 1e-300)
        cosmo_term = Lambda / 3.0 * delta_t
        return wave_term + cosmo_term

    def _lambda_emergent(self, g_ug: float, rho: float) -> float:
        """Λ_UQFF ≈ 26!·g / (ρ·f_vac)²⁷"""
        return self._FAC26 * g_ug / max((rho * self._F_VAC) ** 27, 1e-300)

    def _h_gr(self, q_ddot: float, r: float) -> float:
        """h_GR = G·Q̈ / (c⁴·r)  (quadrupole formula)"""
        return self._G_NEWTON * q_ddot / max(self._C_LIGHT ** 4 * r, 1e-300)

    def compute(self, dataset=None):
        d        = dataset or {}
        kappa    = float(d.get('kappa',    1.0))
        q_ddot   = float(d.get('q_ddot',  1e44))   # kg (DPM quadrupole analog)
        f        = float(d.get('f',        100.0))  # Hz
        r        = float(d.get('r',        3e24))   # m (100 Mpc)
        delta_t  = float(d.get('delta_t',  0.1))    # s
        Lambda   = float(d.get('Lambda',   1e-52))  # m⁻²
        g_ug     = float(d.get('g_ug',     1e-3))
        rho      = float(d.get('rho',      self._RHO_CRIT))
        k        = int(d.get('k', 1))

        h_uqff   = self._h_uqff(kappa, q_ddot, f, r, delta_t, Lambda, k)
        h_gr     = self._h_gr(q_ddot, r)
        lam_pred = self._lambda_emergent(g_ug, rho)

        # SNR G272.2-03.2 specific
        h_snr    = self._h_uqff(kappa, q_ddot, 1e18, 6.6e19, 1.0, Lambda, k)

        return {
            'paper':   'PAPER_580',
            'session': 'Session 156',
            'class':   '#167  UQFFGWAmplitudeLambdaCDMEmergenceCalculator',
            'h_UQFF':          h_uqff,
            'h_GR':            h_gr,
            'h_UQFF_vs_GR_ratio': h_uqff / max(abs(h_gr), 1e-300),
            'Lambda_pred_m2':  lam_pred,
            'Lambda_obs_m2':   1e-52,
            'Lambda_match_pct': abs(lam_pred - 1e-52) / 1e-52 * 100,
            'h_SNR_G272':      h_snr,
            'primary_equations': [
                'h = 26!·κ·Q̈/(f²⁷·r) + Λ/3·δt  [UQFF GW amplitude, k=1]',
                'h = (k+25)!/(k-1)! · κ·Q̈/(f^{k+26}·r) + Λ/3·δt  [general k]',
                'h_GR = G·Q̈/(c⁴·r)  [GR quadrupole for comparison]',
                'Λ_UQFF = 26!·g/(ρ·f_vac)²⁷  [Λ emergent from U_b buoyancy]',
                'Λ_pred(ρ_crit, f_Pl) ≈ 10⁻⁵² m⁻²  [exact match to observed]',
            ],
            'available_equations': [
                'Full 3-system comparison: UQFF vs GR vs LQG amplitude',
                'Λ as function of epoch: f_vac(t) → Λ(t) dark energy evolution',
                'GW frequency spectrum from DPM failure: f_X ~ 10⁸–10¹⁸ Hz',
                'Waveform h(t) = h_UQFF·cos(2πft) with 26! bounding envelope',
            ],
            'simulation_set': [
                {'label': 'GW amplitude frequency sweep',
                 'inputs': 'f range 10–10²¹ Hz, fixed Q̈, r',
                 'output': 'h_UQFF(f) vs h_GR(f)'},
                {'label': 'Λ emergence vs vacuum frequency',
                 'inputs': 'f_vac range 10⁴⁰–10⁴⁵ Hz',
                 'output': 'Λ_UQFF(f_vac) vs Λ_obs=10⁻⁵²'},
                {'label': 'SNR G272.2-03.2 GW spectrum',
                 'inputs': 'f=10⁸–10¹⁸ Hz, r=7kly, ρ=1e-24 g/cm³',
                 'output': 'h profile, DPM failure GW signature'},
            ],
        }


# ---------------------------------------------------------------------------
# #168  UQFFLQGLambdaCDMTripleSystemComparisonCalculator   PAPER_581
# ---------------------------------------------------------------------------
class UQFFLQGLambdaCDMTripleSystemComparisonCalculator:
    """
    #168 – UQFF · LQG · ΛCDM: Simultaneous Three-System QG Comparison
    ------------------------------------------------------------------
    Simultaneous numerical comparison of three QG frameworks for GW propagation:

    1. UQFF (Star-Magic frequency-modulated Form 4):
         h_UQFF = 26!·κ·Q̈/(f²⁷·r) + Λ/3·δt
         Dispersion: d²⁶h/df²⁶ = 0  (26!/f²⁷ bound prevents UV divergence)

    2. GR / ΛCDM (standard cosmology):
         h_GR = G·Q̈/(c⁴·r)  (continuous, no inherent quantum bound)
         Dispersion: ω² = c²k²  (no modification, singularities possible)

    3. LQG (Loop Quantum Gravity, holonomy-corrected):
         Modified dispersion: ω² = c²k²(1 + η(l_Pl·k)^γ)
           η = ±α (ambiguity sign), γ=1 (linear holonomy) or γ=2
         Effective wave equation: (□ + α·l_Pl²·□² + β·l_Pl⁴·∇⁶ + …)h_μν = 0
         Phase velocity: v_ph/c ≈ 1 + (η/2)(l_Pl·k)^γ
         Group velocity shift: δv_g/c ≈ (η·γ/2)(l_Pl·k)^γ

    LQG Holonomy derivation steps:
      1. Classical: □h_μν = 0 → ω² = c²k²
      2. LQG Hamiltonian: H_eff = ∫d³x[sin²(μK)/μ²√q + …], μ~l_Pl√Δ
      3. Expansion sin(μK) ≈ μK − (μK)³/6 → higher-order □² terms
      4. Perturbation (GW modes h_ij^TT):
           (□ + α·l_Pl²·□²)h = 0
      5. Fourier: −ω² + c²k² + α·l_Pl²(ω⁴ − 2ω²c²k² + c⁴k⁴) = 0
      6. Leading correction: ω² ≈ c²k²(1 + η(l_Pl·k)^γ)

    Numerical benchmark (binary merger: Q̈=10⁴⁴ kg, r=100 Mpc, f=100 Hz):
      h_UQFF  ≈ 10⁻²⁰     (26! bound active, larger than GR due to κ factor)
      h_GR    ≈ 10⁻²¹     (standard quadrupole, no quantum bound)
      δv_LQG/c ≈ 10⁻⁴²  (tiny, accumulative over Gpc → time delay testable)

    UQFF improvements over LQG:
      - Bounds divergences factorially (26!/f²⁷) vs LQG possible UV issues
      - Unifies via frequency (motivates forces); LQG has no magnetism/buoyancy
      - Λ emergent dynamically; LQG needs separate dark energy input
      - Fits SNR Chandra data <5% err; LQG corrections ~10⁻⁴² (not yet testable)

    Source: grok_share_efc8a971378f.txt   PAPER_581
    """

    _FAC26   = math.factorial(26)
    _C       = 2.998e8
    _G       = 6.674e-11
    _L_PL    = 1.616e-35
    _HBAR    = 1.055e-34

    def _omega_lqg(self, k_wave: float, eta: float = 1.0, gamma: float = 1.0) -> float:
        """ω_LQG = c·k·√(1 + η(l_Pl·k)^γ)"""
        corr = 1.0 + eta * (self._L_PL * k_wave) ** gamma
        return self._C * k_wave * math.sqrt(max(corr, 0.0))

    def _v_group_lqg(self, k_wave: float, eta: float = 1.0, gamma: float = 1.0) -> float:
        """v_g/c ≈ 1 + (η·γ/2)(l_Pl·k)^γ"""
        return 1.0 + (eta * gamma / 2.0) * (self._L_PL * k_wave) ** gamma

    def _h_uqff_gw(self, q_ddot: float, f: float, r: float,
                   kappa: float, Lambda: float, dt: float) -> float:
        return self._FAC26 * kappa * q_ddot / max(f ** 27 * r, 1e-300) + Lambda / 3 * dt

    def _h_gr(self, q_ddot: float, r: float) -> float:
        return self._G * q_ddot / max(self._C ** 4 * r, 1e-300)

    def compute(self, dataset=None):
        d        = dataset or {}
        q_ddot   = float(d.get('q_ddot',  1e44))
        f        = float(d.get('f',        100.0))
        r        = float(d.get('r',        3e24))
        kappa    = float(d.get('kappa',    1.0))
        Lambda   = float(d.get('Lambda',   1e-52))
        delta_t  = float(d.get('delta_t',  0.1))
        eta      = float(d.get('eta',      1.0))
        gamma    = float(d.get('gamma',    1.0))

        k_wave = 2 * math.pi * f / self._C

        h_uqff   = self._h_uqff_gw(q_ddot, f, r, kappa, Lambda, delta_t)
        h_gr     = self._h_gr(q_ddot, r)
        omega_lqg = self._omega_lqg(k_wave, eta, gamma)
        vg_lqg   = self._v_group_lqg(k_wave, eta, gamma)
        dv_c     = vg_lqg - 1.0   # deviation from c

        # LQG-corrected amplitude (same quadrupole, modified propagation)
        h_lqg    = h_gr * (1.0 + dv_c)  # first-order modification

        # Time delay over r (accumulative LQG dispersion)
        t_travel = r / self._C
        delta_t_lqg = abs(dv_c) * t_travel / self._C

        # UQFF factorial bound
        uqff_bound = self._FAC26 / max(f ** 27, 1e-300)  # 26!/f^27

        return {
            'paper':   'PAPER_581',
            'session': 'Session 156',
            'class':   '#168  UQFFLQGLambdaCDMTripleSystemComparisonCalculator',
            'h_UQFF':     h_uqff,
            'h_GR':       h_gr,
            'h_LQG':      h_lqg,
            'omega_LQG_rad_s':  omega_lqg,
            'v_group_LQG_over_c': vg_lqg,
            'delta_v_over_c': dv_c,
            'time_delay_LQG_s':   delta_t_lqg,
            'UQFF_factorial_bound': uqff_bound,
            'comparison': {
                'UQFF_vs_GR_ratio': h_uqff / max(abs(h_gr), 1e-300),
                'LQG_correction_tiny': bool(abs(dv_c) < 1e-30),
                'UQFF_bounded': bool(uqff_bound < 1e-200),
            },
            'primary_equations': [
                'ω²_LQG = c²k²(1 + η(l_Pl·k)^γ)  [LQG modified dispersion]',
                '(□ + α·l_Pl²·□² + β·l_Pl⁴·∇⁶)h_μν = 0  [LQG effective wave eq]',
                'h_UQFF = 26!·κ·Q̈/(f²⁷·r) + Λ/3·δt  [UQFF Form 4]',
                'h_GR   = G·Q̈/(c⁴·r)  [GR quadrupole]',
                'δv_g/c ≈ (η·γ/2)(l_Pl·k)^γ ≈ 10⁻⁴² at f=100 Hz  [LQG group vel]',
            ],
            'available_equations': [
                'Full LQG effective action with holonomy/inverse volume corrections',
                'LQG time delay over cosmological distances (Gpc)',
                'UQFF vs LQG comparison at SNR G272.2-03.2 (f=10¹⁸ Hz)',
                'Discrete hypergraph vs LQG spin foam equivalence conditions',
            ],
            'simulation_set': [
                {'label': 'Triple system h comparison (f sweep)',
                 'inputs': 'f range 10–10²¹ Hz',
                 'output': 'h_UQFF(f), h_GR(f), h_LQG(f) on log scale'},
                {'label': 'LQG dispersion time delay (distance sweep)',
                 'inputs': 'r range 1 Mpc–1 Gpc, γ=1,2',
                 'output': 'Δt_LQG vs r (accumulation over cosmic distances)'},
                {'label': 'UQFF 26! bound vs f',
                 'inputs': 'f range 10–10²¹ Hz',
                 'output': '26!/f²⁷ envelope showing factorial UV suppression'},
            ],
        }


# ---------------------------------------------------------------------------
# #169  StringGWPlanarFrequencyReboundDiskFormationCalculator   PAPER_582
# ---------------------------------------------------------------------------
class StringGWPlanarFrequencyReboundDiskFormationCalculator:
    """
    #169 – String GW Planar Model with Universal Frequency Rebound and Disk Formation
    ---------------------------------------------------------------------------------
    Expansion of GW string theory to a planar model with holographic principle
    and Universal Frequency Rebound — explaining the angular differential
    between all astronomical disk systems (galaxies, rings, protoplanetary disks).

    Background: pp-wave metric (plane-fronted GWs with parallel rays):
        ds² = −dt² + dz² + dx_i·dx^i + H(x_i, u)·du²
        u = t − z  (lightcone coordinate),  H = A_{ij}·x^i·x^j  (polarization tensor)
        String action: S = ∫dτdσ[∂_α X^μ ∂^α X_μ + rebound·δ(f − f_bound)]

    Universal Frequency Rebound mechanism:
        Incoming string mode f scatters off holographic boundary screen.
        Rebound: f' = f·(1 + δθ)
        Angular differential: δθ = α·(l_s·k)
          α ~ l_s² (string length squared, ~ l_Pl²)
        Rebound torque: J = ∫f·δθ·dA  → stabilizes plane perpendicular to propagation

    Planar dispersion relation:
        ω² = c²k² + α·(f_rebound·k)²
        f_rebound = α·(f/c)²·f  (scales as f³, high-f favors planar alignment)

    Disk formation proof (worldsheet boundary):
        Worldsheet: X^i(σ=0) = X^i(σ=π) + δf  (rebound boundary condition)
        Quantized modes: n = f·L/c  (L = plane size)
        Angular differential: δθ ≈ α·k/f → accumulates to ⊥ rotation (disk form)

    Numerical (galactic disk: f=10⁻¹⁵ Hz orbital, k=10⁻²¹ m⁻¹):
        α = l_s² ≈ l_Pl² ≈ 2.6e-70 m²
        δθ ≈ 2.6e-70 · 10⁻²¹ / 10⁻¹⁵ ≈ 2.6e-76 rad  (small but cumulative over Gyr)
        Over 10 Gyr: ~ 2.6e-76 · 3e17 s ≈ 8e-59 rad/Gyr ... disk aligns to <10°

    Holographic adjustment (AdS/CFT):
        GW information encoded on 2D boundary screen (S^{d-1})
        Amplitude: h ~ ∫_boundary T_μν  (stress-energy on screen)
        Rebound adds: h_holographic = h_rebound · e^{−|δθ|²/2}

    CTAO testability: GW/photon time delays from frequency rebound in SNR shells
        δt_obs = δθ · L_shell / c  (testable for SNR at ~kly distances)

    Source: grok_share_efc8a971378f.txt   PAPER_582
    """

    _C       = 2.998e8
    _L_PL    = 1.616e-35   # m
    _L_S2    = (1.616e-35) ** 2   # l_s² ≈ l_Pl²  m²

    def _delta_theta(self, k_wave: float, f: float,
                     alpha: float = None) -> float:
        """δθ = α·(l_s·k)  [rebound angular differential per string mode]"""
        if alpha is None:
            alpha = self._L_S2
        return alpha * k_wave / max(abs(f), 1e-300)

    def _f_rebound(self, f: float, alpha: float = None) -> float:
        """f_rebound = α·(f/c)²·f  [scales as f³]"""
        if alpha is None:
            alpha = self._L_S2
        return alpha * (f / self._C) ** 2 * f

    def _omega_planar(self, k_wave: float, f: float,
                      alpha: float = None) -> float:
        """ω² = c²k² + α·(f_rebound·k)²"""
        if alpha is None:
            alpha = self._L_S2
        f_reb = self._f_rebound(f, alpha)
        return math.sqrt(max((self._C * k_wave) ** 2 + alpha * (f_reb * k_wave) ** 2, 0.0))

    def _quantized_mode(self, f: float, L_plane: float) -> float:
        """n = f·L/c  (quantized worldsheet mode)"""
        return f * L_plane / self._C

    def _holographic_h(self, h_base: float, delta_theta: float) -> float:
        """h_holographic = h_base · exp(−δθ²/2)  [holographic attenuation]"""
        return h_base * math.exp(-0.5 * delta_theta ** 2)

    def compute(self, dataset=None):
        d        = dataset or {}
        f        = float(d.get('f',          1e-15))   # Hz (galactic orbital default)
        k_wave   = float(d.get('k_wave',     1e-21))   # m⁻¹
        L_plane  = float(d.get('L_plane',    3e22))    # m (galaxy ~30 kpc)
        h_base   = float(d.get('h_base',     1e-22))   # dimensionless GW strain
        alpha    = d.get('alpha', None)
        if alpha is not None:
            alpha = float(alpha)

        dtheta    = self._delta_theta(k_wave, f, alpha)
        f_reb     = self._f_rebound(f, alpha)
        omega_pl  = self._omega_planar(k_wave, f, alpha)
        n_mode    = self._quantized_mode(f, L_plane)
        h_holo    = self._holographic_h(h_base, dtheta)

        # CTAO testability for SNR shell (~5.7 ly = 5.4e16 m)
        L_snr     = 5.4e16  # m
        dtheta_snr = self._delta_theta(2 * math.pi / L_snr, 1e18)
        dt_ctao   = abs(dtheta_snr) * L_snr / self._C

        # Galactic disk cumulative over 10 Gyr
        t_gyr = 10.0 * 3.156e16  # s
        theta_cumulative = abs(dtheta) * t_gyr

        return {
            'paper':   'PAPER_582',
            'session': 'Session 156',
            'class':   '#169  StringGWPlanarFrequencyReboundDiskFormationCalculator',
            'delta_theta_rad':     dtheta,
            'f_rebound_Hz':        f_reb,
            'omega_planar_rad_s':  omega_pl,
            'n_quantized_mode':    n_mode,
            'h_holographic':       h_holo,
            'theta_cumulative_10Gyr_rad': theta_cumulative,
            'CTAO_time_delay_s':   dt_ctao,
            'disk_formation_mechanism': 'rebound_torque_quantization',
            'primary_equations': [
                'ds² = −dt²+dz²+dx_i dx^i + H(x_i,u)du²  [pp-wave metric]',
                'δθ = α·(l_s·k)  [rebound angular differential]',
                'f_rebound = α·(f/c)²·f  [rebound scale ∝ f³]',
                'ω²_planar = c²k² + α·(f_rebound·k)²  [planar dispersion]',
                'n = f·L/c  [quantized worldsheet modes → disk formation]',
                'h_holographic = h·exp(−δθ²/2)  [holographic boundary projection]',
            ],
            'available_equations': [
                'Worldsheet rebound BC: X^i(σ=0) = X^i(σ=π) + δf',
                'Rebound torque: J = ∫f·δθ·dA (angular momentum supply)',
                'δθ accumulation over cosmic time → planar disk angle <10°',
                'Multi-system angular differential table (galaxies, rings, protostellar)',
                'CTAO photon/GW delay prediction for SNR frequency rebound',
            ],
            'simulation_set': [
                {'label': 'Disk formation angular differential (f sweep)',
                 'inputs': 'f range 10⁻¹⁸–10²¹ Hz, L_plane per system type',
                 'output': 'δθ(f), cumulative disk angle over Gyr'},
                {'label': 'Holographic h attenuation vs boundary size',
                 'inputs': 'L_plane range, h_base, α',
                 'output': 'h_holo/h_base vs L (holographic suppression)'},
                {'label': 'CTAO time delay prediction (SNR shells)',
                 'inputs': 'L_shell=5.4e16 m, f range X-ray/GW',
                 'output': 'δt_obs per frequency band → observational prediction'},
            ],
        }

'''

# ── Registry entries to append ───────────────────────────────────────────────
REGISTRY_ENTRIES = '''    # --- Session 156: grok_share_efc8a971378f.txt — UQFF All Forms, GW Amplitude,
    #     UQFF·LQG·ΛCDM Triple Comparison, String GW Planar Rebound Disk ---
    "UQFFAllFormsEvolutionCatalogueCalculator",            # PAPER_579 (#166)
    "UQFFGWAmplitudeLambdaCDMEmergenceCalculator",         # PAPER_580 (#167)
    "UQFFLQGLambdaCDMTripleSystemComparisonCalculator",    # PAPER_581 (#168)
    "StringGWPlanarFrequencyReboundDiskFormationCalculator", # PAPER_582 (#169)
'''

# ── Injection ────────────────────────────────────────────────────────────────
def main():
    cp4_path = "CondensedPhysics4.py"
    with open(cp4_path, "r", encoding="utf-8") as fh:
        src = fh.read()

    # 1. Find insertion point: just before the REGISTRY list (the list with all class names)
    #    We look for the last class block end marker and insert new classes before the registry.
    REGISTRY_MARKER = '    # --- Session 154: Universal Epoch / Periodic Table UQFF ---'
    REGISTRY_END    = '    "UQFFCompEigenvalueQuantumGravityLinkageCalculator",    # PAPER_578 (#165)'
    CLOSE           = '\n]\n'

    if REGISTRY_END not in src:
        print("ERROR: Could not find registry anchor. Aborting.")
        sys.exit(1)

    # 2. Insert new class bodies BEFORE the registry list start
    #    Registry list starts at the first occurrence of "CP4_CALCULATORS = ["
    LIST_START = '__all__ = ['
    if LIST_START not in src:
        print("ERROR: __all__ list not found. Aborting.")
        sys.exit(1)

    insert_pos = src.index(LIST_START)
    src_new = src[:insert_pos] + NEW_CLASSES + '\n\n' + src[insert_pos:]

    # 3. Append registry entries before the closing ']'
    CLOSE_MARKER = '    "UQFFCompEigenvalueQuantumGravityLinkageCalculator",    # PAPER_578 (#165)  \n]'
    if CLOSE_MARKER not in src_new:
        # try without trailing spaces
        CLOSE_MARKER = '    "UQFFCompEigenvalueQuantumGravityLinkageCalculator",    # PAPER_578 (#165)  \n]'
        ALT = '    "UQFFCompEigenvalueQuantumGravityLinkageCalculator",    # PAPER_578 (#165)\n]'
        if ALT in src_new:
            CLOSE_MARKER = ALT
        else:
            # find it more loosely
            idx = src_new.rfind('"UQFFCompEigenvalueQuantumGravityLinkageCalculator"')
            if idx == -1:
                print("ERROR: Cannot locate registry close. Aborting.")
                sys.exit(1)
            # find the ']' after it
            bracket_pos = src_new.index('\n]', idx)
            src_new = src_new[:bracket_pos] + '\n' + REGISTRY_ENTRIES + ']' + src_new[bracket_pos+2:]
            with open(cp4_path, "w", encoding="utf-8") as fh:
                fh.write(src_new)
            print(f"Injection complete via fallback. Lines: {src_new.count(chr(10))+1}")
            return

    src_new = src_new.replace(CLOSE_MARKER,
        CLOSE_MARKER.rstrip('\n]') + '\n' + REGISTRY_ENTRIES + ']')

    with open(cp4_path, "w", encoding="utf-8") as fh:
        fh.write(src_new)

    lines = src_new.count('\n') + 1
    print(f"Injection complete. CP4 lines: {lines}")

if __name__ == "__main__":
    main()
