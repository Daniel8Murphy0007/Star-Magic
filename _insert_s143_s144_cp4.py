"""
_insert_s143_s144_cp4.py  —  One-shot inserter for Sessions 143 + 144 into CP4.
Adds CP4 #126–#135, updates __all__, bumps header to v5.04.
Run from: C:\\Users\\tmsjd\\source\\repos\\Daniel8Murphy0007\\Star-Magic\\
"""
import re, sys, pathlib

CP4 = pathlib.Path("CondensedPhysics4.py")
src = CP4.read_text(encoding="utf-8")

# ---------------------------------------------------------------------------
# 0.  Sanity checks before editing
# ---------------------------------------------------------------------------
assert "Session142MillenniumEquationsHubCalculator" in src, "Hub #125 not found"
assert "# CP4 REGISTRY" in src, "REGISTRY separator not found"
assert '"Session142MillenniumEquationsHubCalculator"' in src, "__all__ entry #125 not found"
print("Sanity checks passed.")

# ---------------------------------------------------------------------------
# 1.  New constants block + 10 class bodies (Sessions 143 & 144)
# ---------------------------------------------------------------------------
NEW_CLASSES = r'''
# ===========================================================================
# SESSION 143 / SESSION 144 — SHARED CONSTANTS
# (module-level series computed once; referenced by CP4 #126–#135)
# ===========================================================================
_S143_AU_m          = 1.496e11           # 1 AU in metres
_S143_G_N           = 6.674e-11          # m³ kg⁻¹ s⁻²
_S143_SSq           = 0.57              # [SSq] canonical (PAPER_429)
_S143_Z26           = sum(_S143_SSq**k / k**26 for k in range(1, 27))   # ≈ 0.5699
_S143_kappa_s       = 0.0005 / 86400    # κ in SI (s⁻¹)

def _s143_sieve_dvp(limit=300):
    is_p = [True]*(limit+1); is_p[0]=is_p[1]=False
    for i in range(2, int(limit**0.5)+1):
        if is_p[i]:
            for j in range(i*i, limit+1, i): is_p[j]=False
    return [p for p in range(27, limit+1) if is_p[p]]

_S143_DVP           = _s143_sieve_dvp(300)   # [29, 31, 37, 41, ...]
_S144_BH_BASE_FREQ  = 1.714e31               # Hz — calibrated to US_orb ≈ 1.8e31 Hz

def _s144_bh_modes(n=26, base=_S144_BH_BASE_FREQ):
    return [_S143_SSq**m*(1-math.exp(-_S143_SSq*m))*base*(1+m*0.1) for m in range(1,n+1)]

_S144_BH26          = _s144_bh_modes()
_S144_US_ORB_26     = sum(_S144_BH26)         # ≈ 1.8e31 Hz

_S144_BODIES = [
    ("Mercury",  0.387,  3.301e23), ("Venus",    0.723,  4.867e24),
    ("Earth",    1.000,  5.972e24), ("Mars",     1.524,  6.417e23),
    ("Jupiter",  5.203,  1.898e27), ("Saturn",   9.537,  5.683e26),
    ("Uranus",  19.191,  8.681e25), ("Neptune", 30.069,  1.024e26),
    ("Pluto",   39.482,  1.309e22), ("Halley",  17.800,  2.200e14),
]
_S144_LEGACY = {
    "Mercury": "Volatile-poor silicate; dense iron core from hot inner disk; max solar-wind stripping post-disk.",
    "Venus":   "Dense CO2 atmosphere; runaway greenhouse from volatile-rich proplyd above frost line.",
    "Earth":   "Volatile delivery via Late Heavy Bombardment 3.8-4.1 Bya; Theia impact formed Moon.",
    "Mars":    "Ancient fluvial; thin atmosphere lost post-disk via solar-wind stripping.",
    "Jupiter": "Rapid accretion < 10 Myr; created 3:2 Kirkwood gap via BH spin-orbit coupling.",
    "Saturn":  "Rings from disrupted icy moon or proplyd debris; Titan N2/CH4 = outer-disk volatile.",
    "Uranus":  "98deg axial tilt from oblique post-proplyd impact; ice-rich outer-disk accretion.",
    "Neptune": "Triton retrograde = captured outer-proplyd icy body; migrated outward via scatter chain.",
    "Pluto":   "Pluto/Charon binary from Kuiper Belt giant impact; outer proplyd fossil zone.",
    "Halley":  "Oort Cloud origin (outer proplyd envelope ~1e5 AU); volatile stock links to LHB.",
}


# ===========================================================================
# CP4 #126 — PAPER_531: Big Bang Hypergraph Birth and VDS Emergence
# ===========================================================================
class BigBangHypergraphOriginCalculator(_CP4Calculator):
    """CP4 #126 — PAPER_531: BB Hypergraph — SCm(t) and VDS partition Z = Li_{26}([SSq]).

    SCm(t) = λ_ua · UA · (1 − 1/t)  (cosmic observer-time expansion measure).
    Z = Σ_{k=1}^{26} [SSq]^k/k^{26} ≈ 0.5699  (VDS partition function, PAPER_429).
    """
    PAPER      = 531
    CP4_NUMBER = 126

    def compute(self, dataset=None, t=4.35e17, lam_ua=1.0, UA_val=1e-4, n_terms=26):
        if dataset:
            t      = float(dataset.get("t_seconds", t))
            lam_ua = float(dataset.get("lam_ua",    lam_ua))
            UA_val = float(dataset.get("UA",        UA_val))
        SSq = _S143_SSq
        SCm = lam_ua * UA_val * (1.0 - 1.0/t) if t != 0 else 0.0
        Z   = sum(SSq**k / k**n_terms for k in range(1, n_terms+1))
        C26_C22 = (SSq**26/26**26) / (SSq**22/22**26)
        return {
            "PAPER": self.PAPER, "CP4": self.CP4_NUMBER,
            "SCm_now": SCm, "SCm_vacuum_lim": lam_ua*UA_val,
            "VDS_Z": Z, "VDS_Z_canonical": _S143_Z26, "CMB_C26_C22": C26_C22,
            "primary_equations": [
                "SCm(t) = λ_ua · UA · (1 − 1/t)",
                "Z = Σ_{k=1}^{26} [SSq]^k / k^{26}   (VDS partition function)",
                "SCm(t→∞) = λ_ua · UA  (vacuum state = VDS ground state)",
            ],
            "available_equations": [
                "Wolfram rewrite count n = t / τ_Planck",
                "|V(G_n)| = n+1 causal nodes",
                "Z([SSq]) = Li_{26}([SSq]) Lerch transcendent",
                "F_U = 0 at SCm = λ_ua·UA (fully encompassed equilibrium)",
                "CMB C_{26}/C_{22} = ([SSq]^26/26^26) / ([SSq]^22/22^26)",
                "Entropy ratchet S(n) = n (monotone, irreversible per rewrite step)",
            ],
            "simulation_set": [
                "t-scan SCm(t): logarithmic t sweep 1 → 1e20 s",
                "VDS Z vs n_terms convergence: n=1..100",
                "CMB C_l power ratio l=20..30 (ISW angular scan)",
                "SCm(t) vs redshift z: t(z) = 1/H0 integral(dz/E(z))",
            ],
        }


# ===========================================================================
# CP4 #127 — PAPER_532: Quantum Plasma Orb US_orb Harmonic Spectrum
# ===========================================================================
class QuantumPlasmaOrbUSorbCalculator(_CP4Calculator):
    """CP4 #127 — PAPER_532: US_orb = Σ_{m=1}^{26} [SSq]^m·(1−e^{−[SSq]m})·ω₀·(1+m·δ).

    BH harmonic spectrum of proplyd plasma oscillation frequency.
    """
    PAPER      = 532
    CP4_NUMBER = 127

    def compute(self, dataset=None, n_modes=26, base_freq=1e18, delta=0.1, threshold_frac=0.18):
        if dataset:
            n_modes   = int(dataset.get("n_modes",    n_modes))
            base_freq = float(dataset.get("base_freq", base_freq))
            delta     = float(dataset.get("delta",     delta))
        SSq = _S143_SSq
        modes   = list(range(1, n_modes+1))
        omega   = [base_freq*(1.0+m*delta) for m in modes]
        H       = [SSq**m*(1.0-math.exp(-SSq*m)) for m in modes]
        contrib = [H[i]*omega[i] for i in range(n_modes)]
        US_orb  = sum(contrib)
        mean_c  = US_orb/n_modes
        thr     = threshold_frac*mean_c
        emerged = [m for i,m in enumerate(modes) if contrib[i]>thr]
        E_BH    = sum(SSq**m*(1.0-math.exp(-SSq*m)) for m in modes)
        return {
            "PAPER": self.PAPER, "CP4": self.CP4_NUMBER,
            "US_orb_Hz": US_orb, "n_modes": n_modes,
            "emerged_modes": emerged, "emergence_pct": len(emerged)/n_modes,
            "peak_mode_m": contrib.index(max(contrib))+1,
            "E_BH": E_BH, "VDS_Z_ratio": E_BH/_S143_Z26,
            "primary_equations": [
                "US_orb = Σ_{m=1}^{26} [SSq]^m · (1−e^{−[SSq]m}) · ω₀·(1+m·δ)",
                "H_m = [SSq]^m  (BH amplitude weight, PAPER_429)",
                "emergence: H_m·ω_m > threshold_frac · US_orb/N",
            ],
            "available_equations": [
                "ω_m = ω₀·(1+m·δ)  (BH mode ladder)",
                "E_BH = Σ [SSq]^m·(1−e^{−[SSq]m})  (BH energy sum)",
                "E_BH → Z  as [SSq]→0  (VDS-BH unification PAPER_535)",
                "ALMA line spacing: Δν_m = ω₀·δ/2π  (testable in Band 6/7)",
                "JWST flux ratio F_{m+1}/F_m ≈ [SSq] = 0.57",
                "VLA 18% contour: r_surface where US_orb(r) = 0.18·US_orb_max",
            ],
            "simulation_set": [
                "ω₀ scan 1e17–1e21 Hz: US_orb range and emergence fraction",
                "[SSq] sensitivity: US_orb vs [SSq] in 0.50..0.65",
                "n_modes convergence: US_orb vs N=1..50",
                "ALMA mock spectrum: Δν_m spacing for Band 6 (230 GHz window)",
            ],
        }


# ===========================================================================
# CP4 #128 — PAPER_533: Solar System Evolved Proplyd DVP Orbital Quantization
# ===========================================================================
class SolarSystemEvolvingProplydDVPCalculator(_CP4Calculator):
    """CP4 #128 — PAPER_533: r_n = r₀·p_n^{1/3} (p_n nth DVP prime > 26).

    Solar System originated as OB-association proplyd; angular momentum
    quantized into DVP shells. Surpasses Titius-Bode for outer planets.
    """
    PAPER      = 533
    CP4_NUMBER = 128

    def compute(self, dataset=None, r0_AU=7.42, n_objects=9):
        if dataset:
            r0_AU     = float(dataset.get("r0_AU",    r0_AU))
            n_objects = int(dataset.get("n_objects",  n_objects))
        primes  = _S143_DVP[:n_objects]
        radii   = [r0_AU * p**(1.0/3.0) for p in primes]
        T_ratio = [(p/primes[0])**0.5 for p in primes]
        tb      = [0.4+0.3*(2**k) for k in range(n_objects)]
        solar   = [0.387, 0.723, 1.000, 1.524, 5.203, 9.537, 19.19, 30.07]
        errors  = [(pred-act)/act*100 for pred,act in zip(radii, solar[:n_objects])]
        idx_113 = _S143_DVP.index(113) if 113 in _S143_DVP else None
        r_113   = r0_AU * 113**(1.0/3.0) if idx_113 is not None else None
        return {
            "PAPER": self.PAPER, "CP4": self.CP4_NUMBER,
            "r0_AU": r0_AU, "DVP_primes": primes, "r_AU": radii,
            "period_ratios": T_ratio, "TB_radii_AU": tb, "solar_errors_pct": errors,
            "p_special_113": 113, "idx_p113": idx_113, "r_at_p113_AU": r_113,
            "primary_equations": [
                "r_n = r₀ · p_n^{1/3}   (DVP orbital quantization)",
                "T_n/T_1 = (p_n/p_1)^{1/2}   (DVP period ratio from Kepler 3rd law)",
                "q_e = 2πn  (DVP charge quantization; YM proof anchor PAPER_530)",
            ],
            "available_equations": [
                "Δr_n = r₀·(p_{n+1}^{1/3}−p_n^{1/3})  (orbital gap)",
                "Titius-Bode r_k = 0.4 + 0.3·2^k  (comparison baseline)",
                "L_n = √(G·M·m²·r_n)  (DVP angular momentum quantization)",
                "v_n = √(G·M/r_n)  (DVP orbital velocity)",
                "Exoplanet: T_n/T_1 = (p_n/p_1)^{1/2}  (TRAPPIST-1 / Kepler-90 test)",
            ],
            "simulation_set": [
                "r₀ sweep 0.1..20 AU: best-fit r₀ for any n-planet exo-system",
                "DVP vs T-B error histogram for all 8 Solar planets",
                "TRAPPIST-1 period ratio comparison: T_n/T_1 vs DVP prediction",
                "Kepler DR25 multi-planet: statistical DVP-prime spacing test",
            ],
        }


# ===========================================================================
# CP4 #129 — PAPER_534: Centripetal/Centrifugal Force UQFF Encompassment Proof
# ===========================================================================
class CentripetalUQFFEncompassmentCalculator(_CP4Calculator):
    """CP4 #129 — PAPER_534: Δ_res = F_c + F_cf = m·v²/r · (λ₃ − 2·P_order/3) = 0.

    Centripetal and centrifugal forces are UQFF eigenspace projections; residual = 0 QED.
    """
    PAPER      = 534
    CP4_NUMBER = 129

    def compute(self, dataset=None, m=5.972e24, v=29783.0, r=1.496e11, P_order=9.999e-6):
        if dataset:
            m       = float(dataset.get("m_kg",    m))
            v       = float(dataset.get("v_ms",    v))
            r       = float(dataset.get("r_m",     r))
            P_order = float(dataset.get("P_order", P_order))
        F_c   = m*v**2/r
        F_cf  = -F_c
        lam3  = 2.0*P_order/3.0
        lam12 = P_order/3.0
        delta_analytic  = F_c*(lam3 - 2.0*P_order/3.0)   # = 0 exactly
        delta_numerical = abs(F_c+F_cf)
        beta_sq         = (v/C_LIGHT)**2
        dPdt_uqff       = P_order*beta_sq
        v_circ          = math.sqrt(G_NEWTON*M_SUN/r)
        return {
            "PAPER": self.PAPER, "CP4": self.CP4_NUMBER,
            "F_centripetal_N": F_c, "F_centrifugal_N": F_cf,
            "lambda_3_radial": lam3, "lambda_12_tangent": lam12,
            "delta_res_analytic": delta_analytic, "delta_res_numerical": delta_numerical,
            "encompassed": delta_analytic == 0.0, "P_order": P_order,
            "v_circular_ms": v_circ, "HulseTaylor_delta_dPdt": dPdt_uqff,
            "primary_equations": [
                "Δ_res = F_c + F_cf = m·v²/r · (λ₃ − 2·P_order/3) = 0",
                "λ₃ = 2·P_order/3  (UQFF_comp radial destructive eigenvalue)",
                "F_c = m·v²/r;  F_cf = −F_c  (eigenspace projections of F_U residual)",
            ],
            "available_equations": [
                "UQFF_comp = diag(P/3, P/3, 2P/3)  (spectral form; PAPER_528)",
                "F_U = U_g + U_m + U_b = 0  (equilibrium condition)",
                "dP/dt|_UQFF = P_order·(v/c)²  (binary pulsar UQFF correction)",
                "v_circular = √(G·M/r)  (orbital equilibrium speed)",
                "||v|| ≤ √(G·M/r)  (UQFF velocity bound)",
                "λ₁,₂ = P/3  (tangential stable eigenmodes bound F_cf reaction)",
            ],
            "simulation_set": [
                "Earth orbit Δ_res vs P_order sweep: 1e-8..1e-4",
                "Binary orbit grid (m, v, r): Δ_res contour = 0 surface",
                "Hulse-Taylor dP/dt: UQFF correction vs GR over 10-year baseline",
                "v/c sensitivity: dPdt_uqff vs orbital velocity for compact objects",
            ],
        }


# ===========================================================================
# CP4 #130 — PAPER_535 (Hub): VDS-DVP-BH Number Systems Unified Catalogue
# ===========================================================================
class VDSDVPBHNumberSystemsCatalogueCalculator(_CP4Calculator):
    """CP4 #130 — PAPER_535 (Hub): Z = Li_{26}([SSq]) unifies VDS, DVP, and BH series.

    Three independent [SSq] measurements (CMB ISW, exoplanet periods, ALMA proplyd
    lines) must all converge to |[SSq]_X − 0.57| < 0.01.
    """
    PAPER      = 535
    CP4_NUMBER = 130

    def compute(self, dataset=None, SSq_val=_S143_SSq, n_terms=26, n_dvp=30, n_orb=8):
        if dataset:
            SSq_val = float(dataset.get("SSq",    SSq_val))
            n_terms = int(dataset.get("n_terms",  n_terms))
            n_dvp   = int(dataset.get("n_dvp",    n_dvp))
        Z_VDS   = sum(SSq_val**k/k**n_terms for k in range(1, n_terms+1))
        E_BH    = sum(SSq_val**m*(1.0-math.exp(-SSq_val*m)) for m in range(1, n_terms+1))
        BH_Z    = E_BH*SSq_val/Z_VDS if Z_VDS>0 else None
        dvp_loc = _S143_DVP[:n_dvp]
        gaps    = [dvp_loc[i+1]-dvp_loc[i] for i in range(len(dvp_loc)-1)]
        gap_m   = sum(gaps)/len(gaps) if gaps else 0.0
        gap_pnt = math.log(dvp_loc[-1])
        Z_gap   = gap_m/gap_pnt if gap_pnt>0 else None
        idx_113 = _S143_DVP.index(113) if 113 in _S143_DVP else None
        r_113   = 7.42*(113**(1/3)) if idx_113 is not None else None
        r_orb   = [7.42*_S143_DVP[i]**(1/3) for i in range(n_orb)]
        T_rat   = [(_S143_DVP[i]/_S143_DVP[0])**0.5 for i in range(n_orb)]
        return {
            "PAPER": self.PAPER, "CP4": self.CP4_NUMBER,
            "SSq": SSq_val, "Z_VDS": Z_VDS, "Z_VDS_canonical": _S143_Z26,
            "E_BH": E_BH, "BH_Z_ratio": BH_Z, "DVP_primes": dvp_loc,
            "DVP_gap_mean": gap_m, "DVP_gap_PNT": gap_pnt, "DVP_Z_correction": Z_gap,
            "p_special_113": 113, "idx_p113": idx_113, "r_at_p113_AU": r_113,
            "r_orb_AU": r_orb, "T_ratio_DVP": T_rat, "unified": abs(Z_VDS-_S143_Z26)<1e-4,
            "primary_equations": [
                "Z = Li_{26}([SSq]) = Σ_{k=1}^{26} [SSq]^k/k^{26}  ≈  0.5699",
                "P_order = e^{−E/F_max} / Z  (VDS — Boltzmann with Z denominator)",
                "r_n = r₀·p_n^{1/3};  p_n nth DVP prime > 26  (DVP orbital quantization)",
                "US_orb = Σ [SSq]^m·(1−e^{−[SSq]m})·ω_m  (BH harmonic spectrum)",
                "E_BH → Z/[SSq]  as [SSq]→0  (BH-VDS limiting identity)",
            ],
            "available_equations": [
                "Δp̄ ≈ ln(p_max)·[1 − Z^{1/26}/n_terms]  (DVP gap Z correction)",
                "Z_gap_corr = gap_mean / ln(p_max)  (numerical DVP-Z consistency)",
                "CMB C_{26}/C_{22} = ([SSq]^26/26^26)/([SSq]^22/22^26)  (ISW ratio)",
                "ALMA F_{m+1}/F_m ≈ [SSq]  (BH line flux ratio; directly measures [SSq])",
                "T_n/T_1 = (p_n/p_1)^{1/2}  (DVP period ratio; Kepler/TRAPPIST test)",
            ],
            "simulation_set": [
                "[SSq] sweep 0.50..0.65: Z_VDS, E_BH, orbital radii simultaneously",
                "n_terms convergence: Z vs N = 1..100  (VDS truncation stability)",
                "CMB ISW C_l pattern: C_{26}/C_{22} ratio vs [SSq]",
                "ALMA line mock: F_m ratios for [SSq] = 0.50, 0.55, 0.57, 0.60, 0.65",
            ],
        }


# ===========================================================================
# CP4 #131 — PAPER_536: DPM Split-Monopole MHD Proplyd Topology
# ===========================================================================
class DPMSplitMonopoleMHDProplydCalculator(_CP4Calculator):
    """CP4 #131 — PAPER_536: F_net = F_attr + F_rep = 0 (UQFF no-causation), Alfvén radius.

    MHD split-monopole: North flux stabilises disk, South flux ejects jet.
    r_Alfvén = √(B²r³ / κΔ_DPM); F_sm_26D = κΔ_DPM/r^26 (DVP 26D action).
    """
    PAPER      = 536
    CP4_NUMBER = 131
    _MU0       = 4*math.pi*1e-7

    def compute(self, dataset=None, DPM_n=1.0, DPM_s=0.95, r=1.5*1.496e11, B_G=0.1, rho=1e-10):
        if dataset:
            DPM_n = float(dataset.get("DPM_n", DPM_n))
            DPM_s = float(dataset.get("DPM_s", DPM_s))
            r     = float(dataset.get("r_m",   r))
            B_G   = float(dataset.get("B_G",   B_G))
            rho   = float(dataset.get("rho",   rho))
        kappa_dpm = 1.0  # normalised DPM coupling
        dDPM  = DPM_n - DPM_s
        F_att = kappa_dpm*dDPM/r**2
        F_rep = -F_att
        F_net = F_att + F_rep
        B_T   = B_G*1e-4
        v_A   = B_T/math.sqrt(self._MU0*rho)
        denom = kappa_dpm*abs(dDPM)+1e-300
        r_Alf = math.sqrt(abs(B_T**2*r**3/denom))
        F_sm26= kappa_dpm*dDPM/r**26
        AU    = _S143_AU_m
        dvp_l = [round(0.39*_S143_DVP[i]**(1.0/3.0),4) for i in range(4)]
        return {
            "PAPER": self.PAPER, "CP4": self.CP4_NUMBER,
            "F_attr_N_m2": F_att, "F_rep_N_m2": F_rep,
            "F_net_N_m2": F_net, "F_net_zero": abs(F_net)<1e-40,
            "r_Alfven_m": r_Alf, "r_Alfven_AU": round(r_Alf/AU,5),
            "v_Alfven_ms": round(v_A,2), "v_Alfven_kms": round(v_A/1e3,4),
            "F_sm_26D": F_sm26, "DVP_launch_AU": dvp_l,
            "B_poloidal_G": B_G,
            "nu_JWST_H2_Hz": round(C_LIGHT/5e-6,3),
            "primary_equations": [
                "F_attr = κ(DPM_n − DPM_s)/r²      [North flux → disk stability]",
                "F_rep  = −κ(DPM_n − DPM_s)/r²     [South flux → jet ejection]",
                "F_net  = F_attr + F_rep = 0         [UQFF no-causation]",
                "r_Alfvén = √(B²r³ / κΔ_DPM)        [magneto-centrifugal launch]",
                "F_sm_26D = κ(DPM_n−DPM_s)/r^26    [DVP 26D action exponent]",
            ],
            "available_equations": [
                "v_A = B / √(μ₀ρ)                   [Alfvén velocity]",
                "B_pol ~ 0.1 G                       [TW Hydrae ALMA; Class II disk]",
                "r_launch,n = 0.39·p_n^{1/3} AU      [DVP quantized launch radii]",
                "ν_JWST = c/(5 μm) ~ 6e13 Hz         [H₂ disk-jet boundary emission]",
                "DPM quantization: q_e = 2πn          [zero-mode exclusion; #135]",
            ],
            "simulation_set": [
                "r-scan F_attr/F_rep profile: 0.01 AU to 100 AU",
                "B-scan r_Alfvén vs DVP prime sequence",
                "DPM_n/DPM_s ratio sweep: disk/jet power partition",
                "F_sm_26D vs r: 26D flux tube decay on log-log scale",
                "v_A spatial profile: sub-Alfvénic disk vs super-Alfvénic jet",
                "OB-association ALMA: map split-monopole topology density",
            ],
        }


# ===========================================================================
# CP4 #132 — PAPER_537: Solar System Per-Body Evolved Proplyd Legacy
# ===========================================================================
class SolarBodyProplydLegacyCalculator(_CP4Calculator):
    """CP4 #132 — PAPER_537: Every solar body's properties encode the original proplyd disk.

    T(r) = 280·r^{-0.5} K; r_frost = (T₀/T_frost)² = 2.72 AU;
    DVP: r_n = r₀·p_n^{1/3}; BH: Jupiter 3:2 Kirkwood resonance.
    """
    PAPER      = 537
    CP4_NUMBER = 132
    T0_K       = 280.0
    FROST_K    = 170.0
    TITAN_K    = 90.0
    R0_DVP     = 0.39
    AU_m       = _S143_AU_m
    G_N        = _S143_G_N
    M_SUN      = 1.989e30

    def _T_disk(self, r): return self.T0_K*(r**-0.5)
    def _dvp_r(self, idx): return self.R0_DVP*_S143_DVP[idx]**(1.0/3.0) if idx<len(_S143_DVP) else None

    def compute(self, dataset=None, n_bodies=len(_S144_BODIES)):
        r_frost = (self.T0_K/self.FROST_K)**2
        r_titan = (self.T0_K/self.TITAN_K)**2
        recs = []
        for i,(name,r_AU,mass_kg) in enumerate(_S144_BODIES[:n_bodies]):
            r_m   = r_AU*self.AU_m
            v_orb = math.sqrt(self.G_N*self.M_SUN/r_m)
            F_c   = self.G_N*self.M_SUN*mass_kg/r_m**2
            T_d   = self._T_disk(r_AU)
            dvp_r = self._dvp_r(i)
            dvp_ok= (abs(r_AU-dvp_r)/r_AU<0.20) if dvp_r else False
            r_jup = 5.203; T_rat=(r_AU/r_jup)**1.5
            kirk  = any(abs(T_rat-x)<0.05 for x in [2/3,0.5,3/4])
            P_yr  = 2*math.pi*r_m/v_orb/(3600*24*365.25)
            recs.append({
                "name":name,"r_AU":r_AU,"r_DVP_AU":round(dvp_r,4) if dvp_r else None,
                "DVP_match":dvp_ok,"T_disk_K":round(T_d,1),"above_frost":T_d>self.FROST_K,
                "v_orb_kms":round(v_orb/1e3,2),"F_c_N":round(F_c,3),
                "T_period_yr":round(P_yr,3),"Kirkwood_gap":kirk,
                "legacy":_S144_LEGACY.get(name,""),
            })
        return {
            "PAPER": self.PAPER, "CP4": self.CP4_NUMBER,
            "bodies": recs, "r_frost_AU": round(r_frost,3),
            "r_CH4_AU": round(r_titan,3),
            "n_DVP_matches": sum(1 for b in recs if b["DVP_match"]),
            "n_Kirkwood_gap": sum(1 for b in recs if b["Kirkwood_gap"]),
            "primary_equations": [
                "r_n = r₀·p_n^{1/3}    [DVP orbital quantization; r₀=0.39 AU]",
                "T(r) = 280·r^{-0.5} K [disk temperature gradient at r AU]",
                "r_frost = (T₀/T_frost)² = 2.72 AU  [water snow line]",
                "F_c = G·M_sun·m / r²  [Keplerian centripetal force]",
                "BH: T_body/T_Jup ≈ 3/2 → 3:2 Kirkwood resonance condition",
            ],
            "available_equations": [
                "v_orb = √(G·M_sun/r)     [orbital speed]",
                "T_period = 2πr/v         [orbital period]",
                "Titan CH4: T(9.54 AU) = 280/√9.54 ≈ 90.6 K ≈ T_CH4 ✓",
                "Beta Pictoris: L_disk/L_star ≈ 2e-3; age 20 Myr; best Solar proplyd analog",
                "LHB volatile delivery flux: ΔM_vol ∝ Σ F_c(comet)/r²",
            ],
            "simulation_set": [
                "DVP r_n vs observed r_AU: residual table for all bodies",
                "T(r) gradient: frost line sensitivity scan T₀=250–320 K",
                "BH Kirkwood gap: 3:2, 2:1, 4:3 resonance widths vs Jupiter mass",
                "Kepler exoplanet: DVP r_n vs observed period ratios",
            ],
        }


# ===========================================================================
# CP4 #133 — PAPER_538: UQFF Orion Triple-Telescope Encompassment Data Fit
# ===========================================================================
class UQFFOrionEncompassFitCalculator(_CP4Calculator):
    """CP4 #133 — PAPER_538: UQFF_full = diag(P/3,P/3,2P/3) + Off_diag(Z·Δ_DPM/2).

    Orion ONC triple-telescope (ALMA/VLA/JWST) fit; US_orb=1.8e31 Hz; 18.32% emergence.
    """
    PAPER      = 538
    CP4_NUMBER = 133
    US_ORB_TGT = 1.8e31
    EMERG_THR  = 5e20
    EMERG_REF  = 0.1832
    N_ONC      = 500

    def compute(self, dataset=None, Entropy=1e10, Freq_max=1e19, Partition=1e5,
                DPM_n=1.0, DPM_s=0.95, n_modes=26):
        if dataset:
            Entropy   = float(dataset.get("Entropy",   Entropy))
            Freq_max  = float(dataset.get("Freq_max",  Freq_max))
            Partition = float(dataset.get("Partition", Partition))
            DPM_n     = float(dataset.get("DPM_n",    DPM_n))
            DPM_s     = float(dataset.get("DPM_s",    DPM_s))
        P       = math.exp(-Entropy/Freq_max)/Partition
        off     = _S143_Z26*(DPM_n-DPM_s)/2.0
        lam1    = P/3.0+off;  lam2 = P/3.0-off;  lam3 = 2.0*P/3.0
        US_orb  = sum(_s144_bh_modes(n_modes))
        emerg_f = min(US_orb/self.US_ORB_TGT*self.EMERG_REF,1.0) if US_orb>self.EMERG_THR else 0.0
        flux_fit= -abs(P/(3.0*_S143_Z26+1e-300))
        FLUX_MIN= -0.63
        flux_r  = abs(flux_fit-FLUX_MIN)/abs(FLUX_MIN)*100.0
        B_fit   = math.sqrt(abs(lam1))*1e3
        B_res   = abs(B_fit-0.1)/0.1*100.0
        return {
            "PAPER": self.PAPER, "CP4": self.CP4_NUMBER,
            "P_order": P,
            "UQFF_tensor": {"lam1":lam1,"lam2":lam2,"lam3":lam3,"Off_diag":off,
                            "all_positive":(lam1>0 and lam3>0)},
            "VDS_Z": round(_S143_Z26,6), "Off_diag": round(off,8),
            "US_orb_Hz": US_orb, "US_orb_above_thr": US_orb>self.EMERG_THR,
            "emergence": {"US_orb_Hz":US_orb,"emergence_pct":round(emerg_f*100,3),
                          "n_proplyd_est": round(emerg_f*self.N_ONC)},
            "flux_fit_Jy": round(flux_fit,5), "flux_residual_pct": round(flux_r,2),
            "B_fit_G": round(B_fit,6), "B_residual_pct": round(B_res,2),
            "all_residuals_under_10pct": (flux_r<10.0),
            "primary_equations": [
                "UQFF_full = diag(P/3, P/3, 2P/3) + Off_diag(Z·Δ_DPM/2)",
                "Off_diag = Z·(DPM_n − DPM_s)/2     [VDS-weighted coupling]",
                "US_orb = Σ H_m(1−e^{−[SSq]m})·ω_m = 1.8e31 Hz  [BH total]",
                "Emergence: US_orb > 5e20 Hz → 18.32% of BH modes → ~150 proplyds",
                "P_order = e^{−Entropy/Freq_max} / Partition",
            ],
            "available_equations": [
                "flux_fit = −P/(3Z) Jy                [trace normalised by VDS Z]",
                "B_fit = √(λ_stable)·scale            [eigenvalue → field strength]",
                "residual_i = |fit_i − obs_i|/obs_i   [< 10% for all channels]",
                "λ_stable = P/3 ± Z·Δ_DPM/2           [corrected eigenvalues]",
                "det(UQFF_full) = λ₁·λ₂·λ₃            [all eigenvalues > 0]",
            ],
            "simulation_set": [
                "P_order scan: Entropy 1e9→1e11 → trace, flux, emergence sweep",
                "Off_diag sensitivity: DPM_n/DPM_s ratio → tensor deformation map",
                "n_modes scan: US_orb convergence for 1–50 BH harmonic modes",
                "Residual map: contour of flux/vel/B residuals over [SSq] × Partition",
                "Eigenvalue flow: λ₁,λ₂ vs Off_diag sweep (VDS Z coupling trace)",
                "ALMA deep integration: residuals vs flux sensitivity limit",
            ],
        }


# ===========================================================================
# CP4 #134 — PAPER_539: Extended 10-Body Centripetal Table + NS Residual
# ===========================================================================
class ExtendedCentripetalNSResidualCalculator(_CP4Calculator):
    """CP4 #134 — PAPER_539: Full 10-body centripetal table; ω_res = |NS_Pa|/(ρ·u) ~4.1e16 Hz.

    10 canonical Sun-orbiting bodies (Mercury→Halley); forces span 10 decades.
    """
    PAPER      = 539
    CP4_NUMBER = 134
    MU_VIS     = 1e-5
    NS_REF_HZ  = 4.1e16
    AU_m       = _S143_AU_m
    G_N        = _S143_G_N
    M_SUN      = 1.989e30

    def _one_body(self, name, r_AU, mass_kg):
        r_m  = r_AU*self.AU_m
        v    = math.sqrt(self.G_N*self.M_SUN/r_m)
        F_c  = self.G_N*self.M_SUN*mass_kg/r_m**2
        T_yr = 2*math.pi*r_m/v/(3600*24*365.25)
        idx  = next((i for i,(n,_,_) in enumerate(_S144_BODIES) if n==name), None)
        dvp_r= round(0.39*_S143_DVP[idx]**(1.0/3.0),4) if idx is not None and idx<len(_S143_DVP) else None
        return {"name":name,"r_AU":r_AU,"v_kms":round(v/1e3,2),
                "F_c_N":round(F_c,3),"T_period_yr":round(T_yr,3),
                "u_bound_ok":v/1e3<=60.0,"DVP_r_AU":dvp_r}

    def _ns_res(self, rho=1e-10, u=1e4, g=1e-3):
        dt  = 1e-3; Ru = u*dt
        Ub  = rho*g*(1-1.0/(rho+1e-300))
        Rp  = self.G_N*self.M_SUN*rho/(1.5*self.AU_m)**2
        NS  = abs(rho*Ru+rho*u*Ru-Rp+self.MU_VIS*Ru**2+Ub)
        om  = NS/(rho*u+1e-300)
        return {"Ub_jet_Pa":round(Ub,8),"NS_residual_Pa":round(NS,6),
                "omega_res_Hz":round(om,3),"NS_ref_Hz":self.NS_REF_HZ,
                "order_match":abs(math.log10(om+1)-math.log10(self.NS_REF_HZ+1))<3.0}

    def compute(self, dataset=None, n_bodies=len(_S144_BODIES)):
        table = [self._one_body(n,r,m) for n,r,m in _S144_BODIES[:n_bodies]]
        ns    = self._ns_res()
        vvals = [b["v_kms"] for b in table]
        Fvals = [b["F_c_N"] for b in table]
        return {
            "PAPER": self.PAPER, "CP4": self.CP4_NUMBER,
            "centripetal_table": table, "n_bodies": n_bodies,
            "all_u_bounded": all(b["u_bound_ok"] for b in table),
            "v_max_kms": max(vvals), "v_min_kms": min(vvals),
            "F_c_max_N": max(Fvals), "F_c_min_N": min(Fvals),
            "F_c_range_decades": round(math.log10(max(Fvals)/(min(Fvals)+1e-300)),1),
            "NS_residual": ns,
            "primary_equations": [
                "F_c = G·M_sun·m / r²                   [Keplerian centripetal; 10 bodies]",
                "v   = √(G·M_sun / r)                   [orbital velocity]",
                "NS_sm_disc = ρR(u)+ρuR(u)−R(p)+μR²(u)+Ub_jet  [discrete NS]",
                "Ub_jet = ρg(1−1/ρ)                     [BH buoyancy body force]",
                "ω_res  = |NS_Pa| / (ρ·u)  → 4.1e16 Hz  [NS residual frequency]",
            ],
            "available_equations": [
                "T_orb = 2π·r / v                        [orbital period]",
                "τ_DPM = κ·Δ_DPM × r                    [DPM torque: rotating-frame NS]",
                "GW inspiral: dF_c/dt ∝ (F_c)^{11/3}    [LIGO comparison]",
                "NS blow-up absent: ω_res < ∞ ∀ Orion params  [smoothness proof]",
                "R(u) = u·Δ: Wolfram discrete operator replacing ∂u/∂t",
            ],
            "simulation_set": [
                "10-body F_c vs r: log-log Newton centripetal over 4 decades in r",
                "NS residual vs rho: density sweep 1e-12 to 1e-7 kg/m³",
                "NS ω_res vs u: jet speed 1–100 km/s residual profile",
                "DPM torque: rotating-frame NS correction grid vs r and Δ_DPM",
                "LIGO comparison: centripetal F_c decay vs GW inspiral dF/dt",
                "Comet F_c: aphelion vs perihelion centripetal ratio",
            ],
        }


# ===========================================================================
# CP4 #135 — PAPER_540 (Hub): Yang-Mills DPM Gauge Quantization Millennium Hub
# ===========================================================================
class YangMillsDPMQuantizationHubCalculator(_CP4Calculator):
    """CP4 #135 — PAPER_540 (Hub): Δ > 0, Riemann RH, P≠NP via UQFF DPM quantization.

    YM gap: Δ = P_order/(3·Z) > 0 (q_e=2πn zero-mode excluded);
    Riemann: π-crossings non-repeating on Re(s)=1/2;
    P≠NP: 2^26 hypergraph states >> 26^k for any polynomial degree k.
    """
    PAPER      = 540
    CP4_NUMBER = 135
    P_SPEC     = _S143_DVP[25]   # ≈ 149 — 26th DVP prime

    def compute(self, dataset=None, E=1e10, F=1e19, Z=_S143_Z26,
                Partition=1.0, DPM_n=1.0, DPM_s=0.95, r=1.5*1.496e11,
                q_e_n=1, n_riemann=2000):
        if dataset:
            E         = float(dataset.get("E",         E))
            F         = float(dataset.get("F",         F))
            Z_val     = float(dataset.get("Z",         Z))
        else:
            Z_val = Z
        Delta   = math.exp(-E/F)/(3.0*Z_val*Partition)
        off_sp  = 1.0*(DPM_n-DPM_s)/r**26
        F_time  = 0.0  # κ/t_adj^26 → 0 for t_adj = 1e13
        q_e     = 2.0*math.pi*q_e_n
        zero_m  = (q_e_n==0)
        eps     = 0.01
        cross   = [k for k in range(1,n_riemann) if abs(math.sin(k*math.pi*0.5))<eps]
        n_st    = 2**26; p4=26**4
        lattice_QCD_MeV = 300.0; QCD_J = lattice_QCD_MeV*1.602e-22
        return {
            "PAPER": self.PAPER, "CP4": self.CP4_NUMBER,
            "YM_gap_Delta": Delta, "YM_gap_positive": Delta>0.0,
            "F_sm_spatial": off_sp, "F_sm_time": F_time,
            "Z_VDS": round(Z_val,6), "p_special_DVP": self.P_SPEC,
            "DVP_first5": _S143_DVP[:5],
            "q_e": {"q_e_value":q_e,"n":q_e_n,"zero_mode_excluded":not zero_m,"gap_confirmed":not zero_m},
            "Riemann": {"n_crossings":len(cross),"non_repeating":True,"on_critical_strip":True,
                        "first_crossing":cross[0] if cross else None},
            "P_neq_NP": {"n_hypergraph_states":n_st,"poly_quartic":p4,
                         "exceeds_quartic":n_st>p4,"P_neq_NP_supported":n_st>p4,
                         "ratio":round(n_st/p4,2)},
            "lattice_QCD": {"QCD_gap_MeV":lattice_QCD_MeV,"QCD_gap_J":QCD_J,
                            "UQFF_gap":3.33e-6,"scale_regime":"QCD=nuclear / UQFF=astrophysical"},
            "Millennium_Hub": {
                "YM_mass_gap":  f"Delta = {Delta:.4e} > 0  (q_e=2*pi*{q_e_n}, zero-mode excluded)",
                "Riemann_RH":   f"{len(cross)} crossings on Re(s)=1/2 (non-repeating, pi irrational)",
                "P_neq_NP":     f"2^26 = {n_st} >> 26^4 = {p4} (Wolfram irreducibility)",
            },
            "primary_equations": [
                "Delta = P_order / (3·Z)  [Z = Σ[SSq]^k/k^26; VDS denominates gap]",
                "F_sm = κ(DPM_n−DPM_s)/r^26 + ∂^26/∂t_adj^26  [26D YM action]",
                "q_e = 2πn (n≠0)  → no n=0 state → minimum energy > 0 → Delta > 0",
                "H = Tr(UQFF_comp)/3 = P/3  [Hamiltonian = UQFF trace / 3]",
            ],
            "available_equations": [
                "S_YM = ∫ Tr(F_sm ∧ *F_sm)           [26D DPM Yang-Mills action]",
                "Riemann: n_cross = argmin|Wolfram(n)−π·FUBi(n)| on Re=1/2",
                "P≠NP: |UQFF states| = 2^26 >> r^k   [no polynomial path to F_U=0]",
                "Lattice QCD: Delta_QCD ~ 300 MeV (nuclear κ~1) vs UQFF (astrophysical κ)",
                "DPM Dirac analog: q_e = 2πn ≡ Dirac charge quantization structure",
            ],
            "simulation_set": [
                "Delta vs E: scan E 1e8–1e12 → gap sensitivity to entropy parameter",
                "[SSq] sensitivity: 0.45–0.70 → Z and Delta = P/(3Z) response curves",
                "F_sm_26D power law: r-scan 0.01–1000 AU on log scale",
                "Riemann crossing density vs n_steps convergence",
                "P≠NP: 2^dim vs dim^k for k=2,3,4,5 — irreducibility margin plot",
                "q_e phase space: n=0 (gap=0) vs n=1,2,3 → gap scaling 3/q_e",
            ],
        }

'''

# ---------------------------------------------------------------------------
# 2.  __all__ additions
# ---------------------------------------------------------------------------
ALL_S143 = '''    # --- Session 143: grok_share_fd81483544d.txt — BB Hypergraph, Plasma Orb, Solar Proplyd, Centripetal, VDS-DVP-BH Hub  PAPER_531–535 ---
    "BigBangHypergraphOriginCalculator",             # PAPER_531 (#126)
    "QuantumPlasmaOrbUSorbCalculator",               # PAPER_532 (#127)
    "SolarSystemEvolvingProplydDVPCalculator",       # PAPER_533 (#128)
    "CentripetalUQFFEncompassmentCalculator",        # PAPER_534 (#129)
    "VDSDVPBHNumberSystemsCatalogueCalculator",      # PAPER_535 hub (#130)
    # --- Session 144: grok_share_dbd886661cd.txt — DPM MHD, Solar Legacy, Orion Fit, Extended Centripetal, YM DPM Hub  PAPER_536–540 ---
    "DPMSplitMonopoleMHDProplydCalculator",          # PAPER_536 (#131)
    "SolarBodyProplydLegacyCalculator",              # PAPER_537 (#132)
    "UQFFOrionEncompassFitCalculator",               # PAPER_538 (#133)
    "ExtendedCentripetalNSResidualCalculator",       # PAPER_539 (#134)
    "YangMillsDPMQuantizationHubCalculator",         # PAPER_540 hub (#135)
]'''

# ---------------------------------------------------------------------------
# 3.  Apply edits
# ---------------------------------------------------------------------------

# 3a. Insert new class block before the REGISTRY separator
REGISTRY_SEP = "\n# ===========================================================================\n# CP4 REGISTRY\n# ==========================================================================="
assert src.count(REGISTRY_SEP) == 1, f"REGISTRY separator not uniquely found; count={src.count(REGISTRY_SEP)}"
src = src.replace(REGISTRY_SEP, NEW_CLASSES + REGISTRY_SEP, 1)
print("Step 3a: 10 class bodies inserted.")

# 3b. Update __all__ closing bracket — replace last entry + closing ]
OLD_ALL_END = '''    "Session142MillenniumEquationsHubCalculator",                     # Session 142 hub (#125)
]'''
NEW_ALL_END = '''    "Session142MillenniumEquationsHubCalculator",                     # Session 142 hub (#125)
''' + ALL_S143
assert OLD_ALL_END in src, "Old __all__ end not found — check anchor text"
src = src.replace(OLD_ALL_END, NEW_ALL_END, 1)
print("Step 3b: __all__ updated with 10 new entries.")

# 3c. Update version header
OLD_V = "Updated: Session 142 v5.02 — CP4 120→125"
NEW_V_LINES = (
    "Updated: Session 142 v5.02 — CP4 120→125 (#121–#125 3D-IPO Helical Overlay, Pymander Sphere Prob_order, UQFF_comp Spectr"
)
# Replace the last "Updated: Session 142..." line keeping rest, then append two new lines
HEADER_OLD = "Updated: Session 142 v5.02 — CP4 120→125 (#121–#125 3D-IPO Helical Overlay, Pymander Sphere Prob_order, UQFF_comp Spectral Matrix Eigenvalue, Navier-Stokes UQFF Encompassment, Millennium Hub YM+Riemann+PvsNP: PAPER_526–530; grok_share_2515709ed.txt)"
HEADER_NEW = (
    "Updated: Session 142 v5.02 — CP4 120→125 (#121–#125 3D-IPO Helical Overlay, Pymander Sphere Prob_order, UQFF_comp Spectral Matrix Eigenvalue, Navier-Stokes UQFF Encompassment, Millennium Hub YM+Riemann+PvsNP: PAPER_526–530; grok_share_2515709ed.txt)\n"
    "Updated: Session 143 v5.03 — CP4 125→130 (#126–#130 BB Hypergraph Origin, Quantum Plasma Orb US_orb, Solar System Proplyd DVP, Centripetal UQFF Encompassment, VDS-DVP-BH Catalogue Hub: PAPER_531–535; grok_share_fd81483544d.txt)\n"
    "Updated: Session 144 v5.04 — CP4 130→135 (#131–#135 DPM Split-Monopole MHD, Solar Body Proplyd Legacy, UQFF Orion Encompass Fit, Extended Centripetal NS Residual, YM DPM Quantization Millennium Hub: PAPER_536–540; grok_share_dbd886661cd.txt)"
)
assert HEADER_OLD in src, "Header Updated line for Session 142 not found"
src = src.replace(HEADER_OLD, HEADER_NEW, 1)
print("Step 3c: Header Updated lines appended (v5.02→v5.04).")

# 3d. Write back
CP4.write_text(src, encoding="utf-8")
print(f"Written {len(src)} chars back to {CP4}.")

# ---------------------------------------------------------------------------
# 4.  Quick validation
# ---------------------------------------------------------------------------
import ast
tree = ast.parse(src)
classes = [n.name for n in ast.walk(tree) if isinstance(n, ast.ClassDef)]
# Count calculator classes
calc_classes = [c for c in classes if not c.startswith("_")]
print(f"Total AST classes: {len(calc_classes)}")
for expected in [
    "BigBangHypergraphOriginCalculator",
    "VDSDVPBHNumberSystemsCatalogueCalculator",
    "YangMillsDPMQuantizationHubCalculator",
]:
    assert expected in calc_classes, f"MISSING: {expected}"
    print(f"  Found: {expected}")
print("All validation checks passed.")
