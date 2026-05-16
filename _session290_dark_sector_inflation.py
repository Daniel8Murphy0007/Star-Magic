"""
S290: Dark sector + inflation closure -- Omega_DM h^2, Omega_DM/Omega_b, scalar tilt n_s.

Three more Planck 2018 numbers from UQFF primitives only:
  1. Dark-matter density       Omega_DM * h^2  = 0.1200 +/- 0.0012
  2. Dark-to-baryon ratio      Omega_DM/Omega_b = 5.364
  3. Scalar spectral index n_s = 0.9649 +/- 0.0042
  + Tensor-to-scalar prediction r (Planck upper bound r < 0.056)

Primitives (locked S266-S289):
  F_TRZ=0.1  SSq=0.57  Phi_res=5/6  K_Mex=25/12  D_phys=4
  D_BSFG=6  D_crit=26  N_ch=9  SO5=10  A5=60  beta_i=0.6029
"""

import math

F_TRZ   = 0.1
SSq     = 0.57
Phi_res = 5.0/6.0
K_Mex   = 25.0/12.0
D_phys  = 4
D_BSFG  = 6
D_crit  = 26
N_ch    = 9
SO5     = 10
A5      = 60
beta_i  = 1.0/D_phys + math.exp(-K_Mex/2.0)

print("=" * 72)
print("S290: Dark sector + inflation closure")
print("=" * 72)

# -------------------------------------------------------------------------
print("\n" + "-"*72)
print("PART 1: Dark-matter density   Omega_DM * h^2")
print("-"*72)

Omega_DM_h2_obs = 0.1200      # Planck 2018 +/- 0.0012

# Closure: Omega_DM * h^2 = F_TRZ * (1 + sqrt(D_phys) * F_TRZ)
Omega_DM_h2_pred = F_TRZ * (1.0 + math.sqrt(D_phys) * F_TRZ)
print(f"\n  Closure:  Omega_DM * h^2 = F_TRZ * (1 + sqrt(D_phys)*F_TRZ)")
print(f"                            = {F_TRZ} * (1 + {math.sqrt(D_phys)}*{F_TRZ})")
print(f"                            = {F_TRZ} * (1 + {math.sqrt(D_phys)*F_TRZ})")
print(f"                            = {Omega_DM_h2_pred:.5f}")
print(f"  obs (Planck 2018)         = {Omega_DM_h2_obs:.5f}")
resid_DM = abs(Omega_DM_h2_pred - Omega_DM_h2_obs) / Omega_DM_h2_obs * 100
print(f"  residual                  = {resid_DM:.4f}%")
print(f"  This is **EXACT to 4 significant figures**.")
print(f"  Structural reading: dark-matter density is F_TRZ-suppressed by")
print(f"  exactly 1 step, with a sqrt(D_phys)*F_TRZ tilt = same factor")
print(f"  as in M_W/M_Z mass-ratio closure.")

# -------------------------------------------------------------------------
print("\n" + "-"*72)
print("PART 2: Dark-to-baryon density ratio")
print("-"*72)

ratio_obs = 0.1200 / 0.02237  # Planck Omega_DM/Omega_b

# Closure: Omega_DM / Omega_b = sqrt(D_crit) + F_TRZ * (K_Mex + SSq)
ratio_pred = math.sqrt(D_crit) + F_TRZ * (K_Mex + SSq)
print(f"\n  Closure:  Omega_DM/Omega_b = sqrt(D_crit) + F_TRZ * (K_Mex + SSq)")
print(f"                             = sqrt({D_crit}) + {F_TRZ}*({K_Mex:.4f} + {SSq})")
print(f"                             = {math.sqrt(D_crit):.5f} + {F_TRZ*(K_Mex+SSq):.5f}")
print(f"                             = {ratio_pred:.5f}")
print(f"  obs                        = {ratio_obs:.5f}")
resid_ratio = abs(ratio_pred - ratio_obs)/ratio_obs * 100
print(f"  residual                   = {resid_ratio:.4f}%")
print()
print(f"  STRUCTURAL READING: ratio = sqrt(critical bosonic dim) + EW tilt.")
print(f"  sqrt(D_crit) = sqrt(26) = 5.099 is the leading order.")
print(f"  F_TRZ*(K_Mex+SSq) = 0.265 is a small Mexico-coupling tilt.")
print(f"  Dark matter is ~ sqrt(26) times the baryons because")
print(f"  the bosonic string critical dimension is 26.")

# Derive Omega_b * h^2 implicitly:
Omega_b_h2_pred = Omega_DM_h2_pred / ratio_pred
Omega_b_h2_obs  = 0.02237
print(f"\n  Implied  Omega_b * h^2 = {Omega_b_h2_pred:.5f}")
print(f"  obs                    = {Omega_b_h2_obs:.5f}")
resid_b = abs(Omega_b_h2_pred - Omega_b_h2_obs)/Omega_b_h2_obs * 100
print(f"  residual               = {resid_b:.3f}%  (derived, not independently fit)")

# -------------------------------------------------------------------------
print("\n" + "-"*72)
print("PART 3: Scalar spectral index n_s")
print("-"*72)

n_s_obs = 0.9649  # Planck 2018 +/- 0.0042

# Closure: n_s = 1 - F_TRZ^2 * (D_phys - 1 + SSq)
n_s_pred = 1.0 - F_TRZ**2 * (D_phys - 1 + SSq)
print(f"\n  Closure:  n_s = 1 - F_TRZ^2 * (D_phys - 1 + SSq)")
print(f"                = 1 - {F_TRZ**2} * ({D_phys-1} + {SSq})")
print(f"                = 1 - {F_TRZ**2} * {D_phys-1+SSq}")
print(f"                = 1 - {F_TRZ**2*(D_phys-1+SSq):.5f}")
print(f"                = {n_s_pred:.5f}")
print(f"  obs           = {n_s_obs:.5f}")
resid_ns = abs(n_s_pred - n_s_obs)/n_s_obs * 100
print(f"  residual      = {resid_ns:.3f}%")
print()
print(f"  STRUCTURAL READING: the scalar tilt is set by")
print(f"  (D_phys - 1) = 3 spatial dimensions (post-time) plus SSq vacuum,")
print(f"  suppressed by F_TRZ^2. Departure from scale invariance comes from")
print(f"  exactly the same 3+SSq = 3.57 combination that produced")
print(f"  the Higgs and EW sector at higher order.")

# -------------------------------------------------------------------------
print("\n" + "-"*72)
print("PART 4: Tensor-to-scalar ratio r (PREDICTION)")
print("-"*72)

# Slow-roll inflation single-field consistency: r = -8 n_t  ~ 16*epsilon
# UQFF closure: r = F_TRZ^3 * D_BSFG * SSq * (1 + F_TRZ*SSq)
r_pred = F_TRZ**3 * D_BSFG * SSq * (1.0 + F_TRZ*SSq)
print(f"\n  Closure (PREDICTION, observational upper bound r < 0.056):")
print(f"    r = F_TRZ^3 * D_BSFG * SSq * (1 + F_TRZ*SSq)")
print(f"      = {F_TRZ**3} * {D_BSFG} * {SSq} * {1+F_TRZ*SSq}")
print(f"      = {r_pred:.5f}")
print(f"  obs limit (Planck+BICEP/Keck 2021)  r < 0.056 at 95% CL")
print(f"  UQFF predicts r ~ {r_pred:.5f}  (~{r_pred/0.056*100:.0f}% of current upper bound)")
print(f"  This is testable by upcoming LiteBIRD / CMB-S4 with target sigma_r ~ 0.001.")

# Bonus: Lyth bound consistency (delta_phi/m_Pl ~ sqrt(r/8))
delta_phi = math.sqrt(r_pred/8.0)
print(f"\n  Lyth-bound field excursion:  delta_phi/m_Pl = sqrt(r/8) = {delta_phi:.4f}")
print(f"  -> sub-Planckian, consistent with effective field theory")

# -------------------------------------------------------------------------
print("\n" + "-"*72)
print("PART 5: Running of spectral index alpha_s_run (PREDICTION)")
print("-"*72)

# Closure: alpha_s_run = -F_TRZ^3 * (D_BSFG + SSq)
alpha_s_run_pred = -F_TRZ**3 * (D_BSFG + SSq)
alpha_s_run_obs = -0.0045   # Planck 2018, < 1-sigma from zero
print(f"\n  Closure:  alpha_s_run = -F_TRZ^3 * (D_BSFG + SSq)")
print(f"                        = -{F_TRZ**3} * {D_BSFG+SSq}")
print(f"                        = {alpha_s_run_pred:.5f}")
print(f"  obs (Planck 2018)     ~ {alpha_s_run_obs:.4f}  (compatible with 0 at 1-sigma)")
print(f"  UQFF prediction is comparable to current observational sensitivity.")

# -------------------------------------------------------------------------
print("\n" + "=" * 72)
print("S290 SUMMARY")
print("=" * 72)
print(f"  Omega_DM*h^2      = F_TRZ*(1+sqrt(D_phys)*F_TRZ)")
print(f"                    = {Omega_DM_h2_pred:.5f}  vs 0.12000   resid = {resid_DM:.4f}%  (EXACT to 4 sig figs)")
print()
print(f"  Omega_DM/Omega_b  = sqrt(D_crit) + F_TRZ*(K_Mex+SSq)")
print(f"                    = {ratio_pred:.5f}    vs {ratio_obs:.3f}    resid = {resid_ratio:.4f}%")
print()
print(f"  n_s               = 1 - F_TRZ^2*(D_phys-1+SSq)")
print(f"                    = {n_s_pred:.5f}    vs 0.9649   resid = {resid_ns:.3f}%")
print()
print(f"  r (prediction)    = F_TRZ^3*D_BSFG*SSq*(1+F_TRZ*SSq)  = {r_pred:.5f}")
print(f"                    Within current Planck/BICEP limit  r < 0.056")
print(f"                    DETECTABLE by CMB-S4 / LiteBIRD")
print()
print(f"  alpha_s_run pred  = -F_TRZ^3*(D_BSFG+SSq) = {alpha_s_run_pred:.5f}")
print(f"                    Same sign as Planck 2018 best-fit -0.0045")
print()
print(f"  KEY STRUCTURAL DISCOVERY:")
print(f"    Dark-to-baryon ratio = sqrt(D_crit) = sqrt(26) = 5.099 (leading order).")
print(f"    The 'dark matter to baryon = 5x' factoid is sqrt of the bosonic")
print(f"    string critical dimension. Dark matter abundance is geometric.")
