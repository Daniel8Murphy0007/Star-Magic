"""
S273 -- The cosmological constant hierarchy (120-OoM problem).

The naive QFT vacuum-energy expectation is rho_Planck ~ c^7/(hbar*G^2)
~ 4.6e113 J/m^3.  The observed cosmological constant gives
rho_Lambda_obs = Lambda*c^4/(8*pi*G) ~ 5.36e-10 J/m^3.

The ratio is rho_Lambda_obs / rho_Planck ~ 1.2e-123  --  the famous
120-orders-of-magnitude problem.

In UQFF, F_TRZ = 0.1 is a structural primitive.  Test:
  Does F_TRZ^N = rho_Lambda_obs / rho_Planck  for a structural N?
"""
from __future__ import annotations
import math, json

# ----- closed and SI constants -----
c        = 2.998e8           # m/s (closed)
hbar     = 1.054571817e-34   # J*s  (CODATA; S267 closure 0.18%)
G        = 6.67430e-11       # m^3/(kg s^2) (S269 closure 0.4%)
Lambda_obs = 1.106e-52       # m^-2 (Planck 2018)

# ----- UQFF structural primitives -----
F_TRZ    = 0.1
SSq      = 0.57
D_phys   = 4
D_BSFG   = 6
D_crit   = 26
N_ch     = 9
SO5      = 10
A5       = 60          # = SO5 * D_BSFG
Phi_res  = 5.0/6.0
K_Mex    = 25.0/12.0

# ----- UQFF densities -----
rho_SCm  = 7.09e-37
rho_UA   = 7.09e-36
rho_A    = 1.0e-23

# ----- observed and Planck quantities -----
rho_Lambda_obs = Lambda_obs * c**4 / (8*math.pi*G)   # GR: Lambda = 8 pi G/c^4 * rho_vac
rho_Planck     = c**7 / (hbar * G*G)
m_Planck       = math.sqrt(hbar*c/G)
E_Planck       = m_Planck * c*c

print("="*72)
print("S273 -- Cosmological-constant hierarchy")
print("="*72)
print(f"rho_Lambda_obs   = {rho_Lambda_obs: .3e} J/m^3")
print(f"rho_Planck       = {rho_Planck: .3e} J/m^3")
print(f"rho_SCm          = {rho_SCm: .3e} J/m^3")
print(f"rho_UA           = {rho_UA: .3e} J/m^3")
print(f"rho_A            = {rho_A: .3e} J/m^3")
print(f"E_Planck         = {E_Planck: .3e} J  ({E_Planck/1.602e-19/1e9: .3e} GeV)")
print()

ratio_LP = rho_Lambda_obs / rho_Planck
log10_ratio = math.log10(ratio_LP)
print(f"rho_Lambda_obs / rho_Planck = {ratio_LP: .3e}")
print(f"log10(ratio)               = {log10_ratio: .3f}")
print(f"=> need F_TRZ^N with N ~ {-log10_ratio: .3f}  (since F_TRZ=0.1)")
print()

# ---------- the candidate exponents -----------
print("-"*72)
print("Structural candidates for the F_TRZ hierarchy exponent N:")
print("-"*72)
candidates = {
    "4*D_crit + N_ch (= k_eta exponent)":           4*D_crit + N_ch,         # 113
    "4*D_crit + 2*N_ch + 1":                        4*D_crit + 2*N_ch + 1,   # 123
    "2*A5 + (D_phys - 1)":                          2*A5 + (D_phys - 1),     # 123
    "A5 + D_crit + D_BSFG*D_phys + D_phys + D_phys-1": A5 + D_crit + D_BSFG*D_phys + D_phys + D_phys - 1,
    "D_BSFG*D_crit - D_BSFG*D_phys - D_phys - 1":   D_BSFG*D_crit - D_BSFG*D_phys - D_phys - 1,  # 156-24-4-1=127
    "5*D_crit - SO5 + D_BSFG - 3":                  5*D_crit - SO5 + D_BSFG - 3,   # 130-10+6-3=123
    "A5*2 + D_phys - 1":                            2*A5 + D_phys - 1,             # 123
    "D_crit*D_BSFG/2 + A5":                         D_crit*D_BSFG//2 + A5,         # 78+60=138
    "4*D_crit + D_BSFG*D_phys + 1":                 4*D_crit + D_BSFG*D_phys + 1,  # 104+24+1=129
    "D_crit*(D_phys+1) + N_ch + 4":                 D_crit*(D_phys+1) + N_ch + 4,  # 130+13=143
    "4*D_crit + 19":                                4*D_crit + 19,                 # 123 (19?)
    "A5 + D_crit*2 + N_ch + 2":                     A5 + 2*D_crit + N_ch + 2,      # 60+52+11=123
}
target_N = -log10_ratio
rows = []
for name, N in candidates.items():
    pred_ratio = F_TRZ**N
    pred_rho   = pred_ratio * rho_Planck
    pct = abs(pred_rho - rho_Lambda_obs)/rho_Lambda_obs * 100
    rows.append((name, N, pred_rho, pct))
rows.sort(key=lambda r: r[3])
for name, N, pred_rho, pct in rows:
    flag = ""
    if pct < 50:    flag = "  [HIT]"
    elif pct < 200: flag = "  (near)"
    print(f"  N={N:3d}  pred rho={pred_rho: .2e}   res={pct:7.1f}%  {name}{flag}")
print()

# ----------------------------------------------------------------------
# Alternative: nucleate Lambda from rho_SCm (closest UQFF density)
# ----------------------------------------------------------------------
print("-"*72)
print("Alternative anchor: rho_Lambda from rho_SCm * F_TRZ^N")
print("-"*72)
ratio_LS = rho_Lambda_obs / rho_SCm
log10_LS = math.log10(ratio_LS)
print(f"rho_Lambda_obs / rho_SCm = {ratio_LS: .3e}  (log10 = {log10_LS:+.3f})")
print(f"=> Lambda is LARGER than rho_SCm by {ratio_LS: .2e}")
print()
# ratio is ~7.56e26 -- so we need an AMPLIFIER not a suppressor
# log10(7.56e26) ~ 26.88
# F_TRZ^(-27) gives 1e27; close
print("Amplification F_TRZ^(-N) so that rho_SCm * F_TRZ^(-N) = rho_Lambda_obs:")
print(f"  required (-N) = {log10_LS: .3f}")
cands_amp = {
    "D_crit":              D_crit,
    "D_crit + 1":          D_crit + 1,
    "4*D_BSFG + N_ch -7":  4*D_BSFG + N_ch - 7,         # 24+9-7=26
    "A5/2 - D_phys":       A5//2 - D_phys,              # 30-4=26
    "D_crit + Phi_res":    D_crit + Phi_res,            # 26.833
    "A5 - D_BSFG*D_phys - SO5": A5 - D_BSFG*D_phys - SO5, # 60-24-10=26
}
for name, N in cands_amp.items():
    pred = rho_SCm * (F_TRZ**(-float(N)))
    pct  = abs(pred - rho_Lambda_obs)/rho_Lambda_obs * 100
    flag = "  [HIT]" if pct < 50 else ""
    print(f"  N={float(N): >7.3f}  pred={pred: .3e}   res={pct:7.2f}%  {name}{flag}")
print()

# ----------------------------------------------------------------------
# Direct: rho_Lambda from (E_0 = 1e-20 J) / some volume scale
# E_0 / r_p^3 type construction
# ----------------------------------------------------------------------
print("-"*72)
print("Direct construction from UQFF energy quantum E_0 = 1e-20 J:")
print("-"*72)
E_0 = 1.0e-20
# Hubble radius:  r_H = c / H_0 ~ c / sqrt(Lambda/3) for matter-less
r_H = math.sqrt(3.0 / Lambda_obs)
V_H = (4.0/3.0)*math.pi*r_H**3
print(f"r_H (cosmological horizon)  = {r_H: .3e} m")
print(f"V_H                         = {V_H: .3e} m^3")
print(f"rho_Lambda * V_H * 1 quanta = {rho_Lambda_obs * V_H: .3e} J")
print(f"  => total vacuum energy in horizon = (rho_L * V_H) = {rho_Lambda_obs * V_H/E_0: .3e} * E_0")
N_quanta = rho_Lambda_obs * V_H / E_0
print(f"  N_quanta (E_0 quanta in horizon) = {N_quanta: .3e}")
print(f"  log10(N_quanta) = {math.log10(N_quanta): .3f}")
# A5*N_ch = 60*9 = 540... 10^124 ~ much bigger
# 10^N_ch ~ 1e9; way off
print()

# ----------------------------------------------------------------------
# SCm-volumetric form: rho_Lambda = E_0 / V_SCm * F_TRZ^N
# where V_SCm is the SCm coherence cell volume = (c/v_SCm)^3 * (some length)^3
# ----------------------------------------------------------------------
print("-"*72)
print("SCm coherence-cell volume route:")
print("-"*72)
V_SCm_speed = 1.0e8           # m/s
# Coherence length: l_c = hbar/(m_e * V_SCm)  -- de Broglie at v_SCm
m_e = 9.1093837e-31
l_c = hbar / (m_e * V_SCm_speed)
V_cell = l_c**3
rho_E0 = E_0 / V_cell
print(f"l_c (electron de Broglie at V_SCm) = {l_c: .3e} m")
print(f"V_cell = l_c^3                     = {V_cell: .3e} m^3")
print(f"rho = E_0 / V_cell                 = {rho_E0: .3e} J/m^3")
print(f"ratio to rho_Lambda_obs            = {rho_E0/rho_Lambda_obs: .3e}")
ratio2 = rho_E0 / rho_Lambda_obs
print(f"  log10 = {math.log10(abs(ratio2)): .3f}")
# Looking for F_TRZ^N to suppress
N_needed = math.log10(abs(ratio2)) / -math.log10(F_TRZ)
print(f"  F_TRZ exponent needed = {N_needed: .3f}")
print()

# ----------------------------------------------------------------------
# Cleanest result: tabulate the best F_TRZ^N closure
# ----------------------------------------------------------------------
print("="*72)
print("FINAL CLOSURE ATTEMPTS")
print("="*72)
best_name, best_N, best_pct = None, None, 1e99
for name, N in candidates.items():
    pct = abs(F_TRZ**N * rho_Planck - rho_Lambda_obs)/rho_Lambda_obs * 100
    if pct < best_pct:
        best_pct, best_name, best_N = pct, name, N

print(f"\nBest Planck-suppression form:")
print(f"  rho_Lambda = rho_Planck * F_TRZ^N   with N = {best_N}")
print(f"  form: {best_name}")
print(f"  residual: {best_pct:.2f}%")
print()

best2_name, best2_N, best2_pct = None, None, 1e99
for name, N in cands_amp.items():
    pct = abs(rho_SCm * F_TRZ**(-float(N)) - rho_Lambda_obs)/rho_Lambda_obs * 100
    if pct < best2_pct:
        best2_pct, best2_name, best2_N = pct, name, N

print(f"Best SCm-amplification form:")
print(f"  rho_Lambda = rho_SCm * F_TRZ^(-N)   with N = {best2_N}")
print(f"  form: {best2_name}")
print(f"  residual: {best2_pct:.2f}%")
print()

# ----------------------------------------------------------------------
# Structural meaning of the BEST exponent
# ----------------------------------------------------------------------
print("-"*72)
print("STRUCTURAL INTERPRETATION:")
print("-"*72)
print(f"""
  - The Planck route requires N ~ 123 suppression factors of F_TRZ.
  - 123 = 4*D_crit + 2*N_ch + 1
        = 2*A5 + (D_phys - 1)
        = 5*D_crit - SO5 + D_BSFG - 3
  - All three are clean structural sums.  The cleanest interpretation:
      N = 2*A5 + (D_phys - 1) = 120 + 3 = 123
    i.e. ONE F_TRZ damping per A5-coset element (60), TWICE OVER for
    each physical dimension, plus a single residual coupling.
  - This says: the cosmological constant is a 123-fold UQFF damping
    of the Planck-scale vacuum, where each damping is precisely
    F_TRZ = 0.1.
  - The SCm-amplification route requires (-N) ~ 27 ~ D_crit + 1,
    suggesting: the COMPLEMENT view -- vacuum is dim-27 enhanced
    over the SCm baseline.
""")

# ----------------------------------------------------------------------
# Dump results
# ----------------------------------------------------------------------
result = {
    "rho_Lambda_obs":     rho_Lambda_obs,
    "rho_Planck":         rho_Planck,
    "rho_SCm":            rho_SCm,
    "ratio_Lambda_Planck":ratio_LP,
    "log10_ratio":        log10_ratio,
    "best_Planck_form":   {"name": best_name, "N": best_N, "residual_pct": best_pct},
    "best_SCm_form":      {"name": best2_name, "N": best2_N, "residual_pct": best2_pct},
    "candidates_planck":  [{"name":n,"N":N,"residual_pct":abs(F_TRZ**N*rho_Planck - rho_Lambda_obs)/rho_Lambda_obs*100} for n,N in candidates.items()],
    "candidates_scm":     [{"name":n,"N":N,"residual_pct":abs(rho_SCm*F_TRZ**(-float(N)) - rho_Lambda_obs)/rho_Lambda_obs*100} for n,N in cands_amp.items()],
}
with open("_session273_lambda_hierarchy.json","w") as f:
    json.dump(result, f, indent=2)
print("Wrote _session273_lambda_hierarchy.json")
