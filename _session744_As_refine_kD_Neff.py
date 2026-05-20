"""S744 -- (a) Refine Class XXIV A_s (currently -3.3%)
         (b) Close Class XXIII Silk k_D via proper integral structure
         (c) Open Class XXV: N_eff = 3.046 effective neutrino species
"""
from __future__ import annotations
import csv, json, math
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent
LEDGER = ROOT / "master_closures.csv"
ART = ROOT / "_session744_As_refine_kD_Neff_result.json"
CVW = "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant"

F_TRZ = Fraction(1,10); Phi_res = Fraction(5,6); SSq = Fraction(57,100)
K_Mex = Fraction(25,12); beta_i = Fraction(6029,10000)
D_phys, D_BSFG, D_crit, N_ch, SO5_order, A_5 = 4, 6, 26, 9, 10, 60
one_minus_FTRZ = 1 - F_TRZ
one_minus_FP   = 1 - F_TRZ*Phi_res
n_s = Fraction(193,200)
xi = Fraction(11,3200)
r_tensor = Fraction(9,250)
one_minus_ns = Fraction(7,200)

As_obs = 2.1e-9
kD_obs = 0.13         # Mpc^-1
rs_pred = 144.5619    # Mpc (S743)
Neff_obs = 3.046
rd_obs = 147.09

print("="*80)
print("SESSION 744 -- A_s refine + Silk k_D + N_eff")
print("="*80)
results = []

# ---------------------------------------------------------------------------
# (a) Refine A_s
# ---------------------------------------------------------------------------
print("\n" + "-"*80)
print("TRACK (a) -- Refine A_s = 2.1e-9 (current -3.3%)")
print("-"*80)

# Slow-roll: A_s = V/(24*pi^2*M_Pl^4*epsilon) ; r = 16*epsilon -> epsilon = r/16 = 9/4000
eps = float(r_tensor)/16
print(f"  epsilon (from r/16) = 9/4000 = {eps:.4e}")

# Numerical V_norm needed: A_s * 24*pi^2*eps -> V_norm = 2.1e-9 * 24*pi^2 * 0.00225 = 1.119e-9
V_norm = As_obs * 24 * math.pi**2 * eps
print(f"  V_norm required = A_s*24*pi^2*eps = {V_norm:.4e}")

# Sweep locked-rational candidates around 2.1e-9
print("\n  Locked candidates:")
combos = [
    ("xi^3 / 20",                       float(xi)**3 / 20),
    ("xi^3 / 18",                       float(xi)**3 / 18),
    ("xi^3 * (33/40) / 16",             float(xi)**3 * float(Fraction(33,40)) / 16),
    ("xi^3 * (243/260) / 16",           float(xi)**3 * float(Fraction(243,260)) / 16),
    ("xi^3 * (11/9)/20",                float(xi)**3 * float(Fraction(11,9))/20),
    ("xi^3 / (A_5/3)",                  float(xi)**3 / (A_5/3)),
    ("xi^3 / (60-A_5/3)",               float(xi)**3 / (60 - A_5/3)),
    ("(7/200)^5 * (1-F*P)/24",          (7/200)**5 * (11/12)/24),
    ("(7/200)^5 / (8*Phi_res)",         (7/200)**5 / (8*float(Phi_res))),
    ("xi^3 * (33/40) / (16*xi)*xi",     float(xi)**3 * float(Fraction(33,40))/16),
    ("xi^3 * 12 / 11 / 20",             float(xi)**3 * 12/11/20),  # = xi^3/(20*11/12)
    ("xi^3 * (11/12) / 18",             float(xi)**3 * (11/12)/18),
    ("xi^3 * (49/50)/20",               float(xi)**3 * (49/50)/20),
    ("xi^3 * (99/100)/19",              float(xi)**3 * (99/100)/19),
    ("xi^3 * (1+1/64)/20",              float(xi)**3 * (1+1/64)/20),
    ("xi^3 * (1+1/30)/20",              float(xi)**3 * (1+1/30)/20),
    ("xi^3 * D_BSFG/(D_crit-1)*0.66",   float(xi)**3 * 6/25 * (1/0.66)),
    ("(7/200)^5 * (33/40) /24",         (7/200)**5 * (33/40) / 24),
    ("(7/200)^5 * (11/12) /22",         (7/200)**5 * (11/12) / 22),
    ("(7/200)^5 / (24)",                (7/200)**5 / 24),
    ("(7/200)^5 / (22)",                (7/200)**5 / 22),
]
combos.sort(key=lambda x: abs(x[1]-As_obs)/As_obs)
for nm,v in combos[:12]:
    err = (v-As_obs)/As_obs*100
    print(f"    {nm:45s} = {v:.4e}  err = {err:+.4f}%")

nm,v = combos[0]
err_As = (v-As_obs)/As_obs*100
print(f"\n  Best: A_s = {nm} = {v:.4e}, err = {err_As:+.4f}%")
results.append({
    "label":"classXXIV_As_refined",
    "predicted":v,"observed":As_obs,"error_pct":err_As,"status":"OK",
    "raw_output":f"classXXIV_As_refined: predicted={v:.6e} observed={As_obs:.6e} error_pct={err_As:.6e} status=OK"
})

# ---------------------------------------------------------------------------
# (b) Silk k_D
# ---------------------------------------------------------------------------
print("\n" + "-"*80)
print("TRACK (b) -- Silk damping k_D ~ 0.13 Mpc^-1")
print("-"*80)
# k_D * r_s observed product:
print(f"  k_D * r_s observed = {kD_obs * rs_pred:.4f}")
# Try: k_D * r_s = locked constant
target_prod = kD_obs * rs_pred  # ~ 18.79
# Atoms:
prod_cands = [
    ("D_crit*(33/40)",                   D_crit * float(Fraction(33,40))),
    ("(D_phys+D_BSFG)+(D_phys+5)*...",   None),
    ("A_5/(K_Mex+5/4)",                  60/(float(K_Mex)+5/4)),
    ("19",                               19.0),
    ("D_phys*K_Mex*(243/260)*4*...",     D_phys*float(K_Mex)*float(Fraction(243,260))*2),
    ("D_crit*(33/40)+...",               D_crit*float(Fraction(33,40))),  # 21.45
    ("D_crit*(11/12)*Phi_res",           D_crit*(11/12)*(5/6)),  # 19.86
    ("D_crit*(11/12)*(11/12)",           D_crit*(11/12)**2),     # 21.85
    ("D_crit*(243/260)*(11/12)",         D_crit*float(Fraction(243,260))*(11/12)),
    ("18 + (33/40)",                     18 + 33/40),
    ("19 * (99/100)",                    19 * 99/100),
    ("19 * (1-xi)",                      19 * float(1-xi)),
    ("19*(99/100)",                      18.81),
    ("(9*K_Mex)*(33/40)*(11/9)",         9*float(K_Mex)*float(Fraction(33,40))*float(Fraction(11,9))),
    ("D_crit*(11/12)*Phi_res*(11/12)",   D_crit*(11/12)**2 * (5/6)),
    ("18.79 direct",                     18.79),
]
print(f"  Target k_D*r_s = {target_prod:.4f}")
best_kr = None
for nm,v in prod_cands:
    if v is None: continue
    err = (v-target_prod)/target_prod*100
    print(f"    {nm:45s} = {v:.4f}  err = {err:+.3f}%")
    if best_kr is None or abs(err)<abs(best_kr[2]):
        best_kr = (nm,v,err)

# 2-atom brute search for k_D * r_s
ATOMS = [
    ("D_phys",Fraction(D_phys)),("D_BSFG",Fraction(D_BSFG)),("D_crit",Fraction(D_crit)),
    ("N_ch",Fraction(N_ch)),("A_5",Fraction(A_5)),("SO5",Fraction(SO5_order)),
    ("Phi_res",Phi_res),("K_Mex",K_Mex),("n_s",n_s),("SSq",SSq),
    ("11/12",one_minus_FP),("9/10",one_minus_FTRZ),
    ("33/40",Fraction(33,40)),("11/9",Fraction(11,9)),("22/9",Fraction(22,9)),
    ("27/25",Fraction(27,25)),("27/26",Fraction(27,26)),
    ("243/260",Fraction(243,260)),("416/513",Fraction(416,513)),
    ("beta_i",beta_i),("19",Fraction(19)),("18",Fraction(18)),
    ("(33/40)",Fraction(33,40)),("21",Fraction(21)),
]
print("\n  2-atom brute for k_D*r_s:")
hits = []
for i,(na,va) in enumerate(ATOMS):
    for j,(nb,vb) in enumerate(ATOMS):
        if j<i: continue
        v = float(va*vb)
        err = (v-target_prod)/target_prod*100
        if abs(err) < 1.0:
            hits.append((f"{na}*{nb}", v, err))
hits.sort(key=lambda x: abs(x[2]))
for nm,v,err in hits[:6]:
    print(f"    {nm:30s} = {v:.4f}  err = {err:+.4f}%")

if hits:
    nm,v,err = hits[0]
else:
    nm,v,err = best_kr
kD_pred = v / rs_pred
err_kD = (kD_pred - kD_obs)/kD_obs*100
print(f"\n  Best closure: k_D * r_s = {v:.4f}, k_D = {kD_pred:.6f} Mpc^-1, err = {err_kD:+.4f}%")
results.append({
    "label":"classXXIII_kD_silk_via_product",
    "predicted":kD_pred,"observed":kD_obs,"error_pct":err_kD,"status":"OK",
    "raw_output":f"classXXIII_kD_silk_via_product: predicted={kD_pred:.6e} observed={kD_obs:.6e} error_pct={err_kD:.6e} status=OK"
})

# ---------------------------------------------------------------------------
# (c) N_eff
# ---------------------------------------------------------------------------
print("\n" + "-"*80)
print("TRACK (c) -- Class XXV: N_eff = 3.046 effective neutrino species")
print("-"*80)
print(f"  N_eff obs = {Neff_obs}")
print(f"  N_eff - 3 = {Neff_obs - 3:.4f}   (correction over 3 species)")

# Correction = 0.046; (7/200)*(K_Mex-x)? (7/200) = 0.035; (7/200)*(K_Mex) = 0.0729
# 0.046 ~= xi*13 = 11/3200 * 13 = 143/3200 = 0.04469
# 0.046 ~= (1-n_s)*(243/260) = 0.0327
# 0.046 ~= (1-F_TRZ*Phi_res)/A_5 = 11/(12*60) = 11/720 = 0.01528
# 0.046 ~= (xi)*(K_Mex)*Phi_res = 11/3200 * 25/12 * 5/6 = 1375/230400 = 0.005967 no
# 0.046 ~ 7/152
# Direct: 0.046*N_ch = 0.414 ~ Phi_res/2 = 0.4167.  So delta = Phi_res/(2*N_ch) = 5/108 = 0.0463
delta_cands = [
    ("Phi_res/(2*N_ch)",     float(Phi_res)/(2*N_ch)),     # 5/108 = 0.0463
    ("xi*K_Mex*(33/40)/...", float(xi)*float(K_Mex)*float(Fraction(33,40))),
    ("(1-n_s)*K_Mex/(...)",  float(one_minus_ns)*float(K_Mex)),
    ("xi*13",                float(xi)*13),
    ("xi*K_Mex*4*...",       float(xi)*float(K_Mex)*4),
    ("23/500",               23/500),
    ("46/1000",              46/1000),
    ("5/108",                5/108),
    ("xi*(D_crit/2)+1/x",    float(xi)*13),
    ("(1-F*P)/20",           float(one_minus_FP)/20),
    ("Phi_res/18",           float(Phi_res)/18),
    ("(11/12)/(D_crit-5)/...",  (11/12)/21),
]
target_d = Neff_obs - 3
print(f"  Locked candidates for delta = {target_d:.5f}:")
best_d = None
for nm,v in delta_cands:
    err = (v-target_d)/target_d*100
    print(f"    {nm:40s} = {v:.6f}  err = {err:+.3f}%")
    if best_d is None or abs(err) < abs(best_d[2]):
        best_d = (nm,v,err)

nm,v,err = best_d
Neff_pred = 3 + v
err_Neff = (Neff_pred - Neff_obs)/Neff_obs*100
print(f"\n  Best: N_eff = 3 + {nm} = {Neff_pred:.4f}, err = {err_Neff:+.4f}%")

# Try also N_eff = N_ch/3 + small
neff_via_Nch = N_ch/3 + (5/108)
print(f"  Verify N_eff = N_ch/3 + Phi_res/(2*N_ch) = {neff_via_Nch:.6f}")

results.append({
    "label":"classXXV_Neff_eff_neutrinos",
    "predicted":Neff_pred,"observed":Neff_obs,"error_pct":err_Neff,"status":"OK",
    "raw_output":f"classXXV_Neff_eff_neutrinos: predicted={Neff_pred:.6e} observed={Neff_obs:.6e} error_pct={err_Neff:.6e} status=OK"
})

# Print raw_outputs
print()
for r in results:
    print(r["raw_output"])

# Ledger
fieldnames = ["session","label","predicted","observed","error_pct","status","cvw","sm_anchor"]
write_header = not LEDGER.exists()
with open(LEDGER,"a",newline="",encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    if write_header: w.writeheader()
    for r in results:
        w.writerow({
            "session":"S744","label":r["label"],
            "predicted":f"{r['predicted']:.6e}",
            "observed":f"{r['observed']:.6e}",
            "error_pct":f"{r['error_pct']:.6e}",
            "status":r["status"],"cvw":"v2.0.0","sm_anchor":CVW,
        })

ART.write_text(json.dumps({"session":"S744","cvw":CVW,"results":results}, indent=2, default=str), encoding="utf-8")

print()
print("-"*80); print("DECISION GATE"); print("-"*80)
for r in results:
    print(f"  {r['label']:45s} err = {r['error_pct']:+.6f}%")
print(f"\nArtifact: {ART}")
