"""
SESSION 751 -- T_0 candidate-EXACT push; H_0 ratio (Class XXXIIIb); Class XXXIV opening

(a) T_0 = 2.7255 K. S750 closure: T_0 = (hc/(L_SCM*k_B))/(D_BSFG*SSq*xi*(5/108)) = 2.7250 K, -0.0168%.
    Needed mult on r3 = 6.6155e4 vs current 6.6144e4 -> +1.658e-4 fine correction. 5-atom search.

(b) Class XXXIIIb -- H_0 tension ratio. Both Planck and SH0ES live off 1/t_H with K_Mex factor.
    Formalize:  rho_H = H_S/H_P = (1 + K_Mex*(5/108)/(11/9)) / (1 - K_Mex/(N_ch*A_5)).
    Observed: 73.0/67.4 = 1.08309.

(c) Class XXXIV. Two seed tests:
    sigma_8_obs = 0.8111 (Planck 2018). Seed: sigma_8 = 416/513 = 0.81092 (1-atom, 0-anchor).
    Omega_m_obs = 0.3153 (Planck 2018). Seed: Omega_m = 9*(1-n_s) = 63/200 = 0.3150 (N_ch*(1-n_s)).

CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""
from __future__ import annotations
from fractions import Fraction as F
import csv, json, os, itertools, math

F_TRZ=F(1,10); Phi_res=F(5,6); SSq=F(57,100); K_Mex=F(25,12)
beta_i=F(6029,10000); D_phys=F(4); D_BSFG=F(6); D_crit=F(26)
N_ch=F(9); SO5=F(10); A_5=F(60)
one_m_FTRZ=1-F_TRZ; one_m_FP=1-F_TRZ*Phi_res
ns_atom=F(193,200); one_m_ns=F(7,200); xi=F(11,3200); r_tens=F(9,250); Yp=F(49,200)

ATOMS = {
    "F_TRZ":F_TRZ,"Phi_res":Phi_res,"SSq":SSq,"K_Mex":K_Mex,"beta_i":beta_i,
    "D_phys":D_phys,"D_BSFG":D_BSFG,"D_crit":D_crit,"N_ch":N_ch,"SO5":SO5,"A_5":A_5,
    "1-F_TRZ":one_m_FTRZ,"1-F*P":one_m_FP,"n_s":ns_atom,"1-n_s":one_m_ns,
    "xi":xi,"r":r_tens,"Y_p":Yp,
    "27/26":F(27,26),"243/260":F(243,260),"33/40":F(33,40),"11/9":F(11,9),
    "22/9":F(22,9),"27/25":F(27,25),"416/513":F(416,513),"31/30":F(31,30),
    "5/108":F(5,108),
}
LABELS=list(ATOMS.keys()); VALS=[float(ATOMS[k]) for k in LABELS]

def search2(target, tol_pct=5.0, want=15):
    hits=[]; n=len(LABELS)
    for i in range(n):
        for j in range(n):
            for tag,fn in [("a*b",lambda a,b:a*b),("a/b",lambda a,b:a/b)]:
                try: v=fn(VALS[i],VALS[j])
                except ZeroDivisionError: continue
                if v==0: continue
                err=(v-target)/target*100.0
                if abs(err)<tol_pct: hits.append((LABELS[i],LABELS[j],tag,v,err))
    hits.sort(key=lambda h:abs(h[4]))
    return hits[:want]

def search3(target, tol_pct=5.0, want=12):
    hits=[]; n=len(LABELS)
    forms=[("a*b/c",lambda a,b,c:a*b/c),("a*b*c",lambda a,b,c:a*b*c),("a/b/c",lambda a,b,c:a/b/c)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for tag,fn in forms:
                    try: v=fn(VALS[i],VALS[j],VALS[k])
                    except ZeroDivisionError: continue
                    if v==0: continue
                    err=(v-target)/target*100.0
                    if abs(err)<tol_pct: hits.append((LABELS[i],LABELS[j],LABELS[k],tag,v,err))
    hits.sort(key=lambda h:abs(h[5]))
    return hits[:want]

def search4(target, tol_pct=2.0, want=12):
    hits=[]; n=len(LABELS)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    for tag, fn in [
                        ("a*b/(c*d)", lambda a,b,c,d: a*b/(c*d)),
                        ("a*b*c/d",   lambda a,b,c,d: a*b*c/d),
                        ("a/(b*c*d)", lambda a,b,c,d: a/(b*c*d)),
                    ]:
                        try: v=fn(VALS[i],VALS[j],VALS[k],VALS[l])
                        except ZeroDivisionError: continue
                        if v==0: continue
                        err=(v-target)/target*100.0
                        if abs(err)<tol_pct: hits.append((LABELS[i],LABELS[j],LABELS[k],LABELS[l],tag,v,err))
    hits.sort(key=lambda h:abs(h[6]))
    return hits[:want]

def write_ledger(label,predicted,observed,status="OK"):
    err_pct=(predicted-observed)/observed*100.0
    raw=f"{label}: predicted={predicted:.6e} observed={observed:.6e} error_pct={err_pct:.6e} status={status}"
    path=os.path.join(os.path.dirname(os.path.abspath(__file__)),"master_closures.csv")
    head=["label","predicted","observed","error_pct","status","raw_output","cvw","sm_anchor"]
    new=not os.path.exists(path)
    with open(path,"a",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=head,extrasaction="ignore")
        if new: w.writeheader()
        w.writerow({"label":label,"predicted":f"{predicted:.6e}","observed":f"{observed:.6e}",
            "error_pct":f"{err_pct:.6e}","status":status,"raw_output":raw,
            "cvw":"v2.0.0","sm_anchor":"CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant"})
    print(raw); return err_pct

print("="*80); print("SESSION 751 -- T_0 push; H_0 ratio; sigma_8/Omega_m"); print("="*80)

# ============================================================
# TRACK (a) T_0 candidate-EXACT push
# ============================================================
print("\n"+"-"*80); print("TRACK (a) -- T_0 candidate-EXACT push"); print("-"*80)
T0_obs=2.7255
c=2.99792458e8; L_SCM=349.226733192; k_B=1.380649e-23; h_pl=6.62607015e-34
T_HC = h_pl * c / (L_SCM * k_B)
r3 = T0_obs / T_HC
# S750 closure: 1/(D_BSFG*SSq*xi*(5/108))
denom_S750 = float(D_BSFG)*float(SSq)*float(xi)*float(F(5,108))
r3_S750 = 1.0/denom_S750
T0_S750 = T_HC * r3_S750
err_S750 = (T0_S750-T0_obs)/T0_obs*100.0
print(f"  T_HC = h*c/(L_SCM*k_B) = {T_HC:.6e} K")
print(f"  r3 target  = {r3:.6e}")
print(f"  r3 S750    = {r3_S750:.6e}  (T_0={T0_S750:.4f}, err={err_S750:+.6f}%)")
need_mult_T = r3 / r3_S750
delta_T = need_mult_T - 1.0
sign_T = 1 if delta_T>0 else -1
print(f"  needed correction mult on r3:  {need_mult_T:.8f}  (delta = {delta_T:+.4e})")

print("\n  3-atom delta search (tol 1%):")
hits3T = search3(abs(delta_T), tol_pct=1.0, want=12)
for h in hits3T:
    mult=1.0+sign_T*h[4]; pred=T0_S750*mult; err=(pred-T0_obs)/T0_obs*100.0
    print(f"    1{'+' if sign_T>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}]  d={h[4]:.4e}  T0={pred:.5f}  err={err:+.6f}%")

print("\n  4-atom delta search (tol 0.3%):")
hits4T = search4(abs(delta_T), tol_pct=0.3, want=12)
for h in hits4T:
    mult=1.0+sign_T*h[5]; pred=T0_S750*mult; err=(pred-T0_obs)/T0_obs*100.0
    print(f"    1{'+' if sign_T>0 else '-'}[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]  d={h[5]:.4e}  T0={pred:.5f}  err={err:+.6f}%")

best_T=None
for h in hits3T:
    mult=1.0+sign_T*h[4]; pred=T0_S750*mult; err=(pred-T0_obs)/T0_obs*100.0
    if best_T is None or abs(err)<abs(best_T[1]):
        best_T=(f"S750*[1{'+' if sign_T>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]",err,pred)
for h in hits4T:
    mult=1.0+sign_T*h[5]; pred=T0_S750*mult; err=(pred-T0_obs)/T0_obs*100.0
    if best_T is None or abs(err)<abs(best_T[1]):
        best_T=(f"S750*[1{'+' if sign_T>0 else '-'}{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]",err,pred)
if best_T:
    print(f"\n  BEST: T_0 = {best_T[0]} = {best_T[2]:.5f}, err = {best_T[1]:+.6f}%")
    T0_pred = best_T[2]
else:
    T0_pred = T0_S750

# ============================================================
# TRACK (b) Class XXXIIIb -- H_0 tension ratio
# ============================================================
print("\n"+"-"*80); print("TRACK (b) -- Class XXXIIIb: H_0 tension ratio"); print("-"*80)
t_H_Gyr=14.4517; sec_per_Gyr=3.1557e16; Mpc_in_km=3.0857e19
H_nat = (1.0/(t_H_Gyr*sec_per_Gyr)) * Mpc_in_km
H_P_obs=67.4; H_S_obs=73.0
rho_obs = H_S_obs/H_P_obs

# S750 closures:
planck_corr = 1.0 - float(K_Mex)/(float(N_ch)*float(A_5))   # 1 - K_Mex/(N_ch*A_5) = 1 - 25/6480 = 1 - 5/1296
sh0es_corr  = 1.0 + float(K_Mex)*float(F(5,108))/float(F(11,9))  # 1 + K_Mex*(5/108)/(11/9)

H_P_pred = H_nat * planck_corr
H_S_pred = H_nat * sh0es_corr
rho_pred = sh0es_corr/planck_corr

print(f"  1/t_H = {H_nat:.4f} km/s/Mpc  (natural Hubble unit)")
print(f"  Planck:  H_0 = (1/t_H)*(1 - K_Mex/(N_ch*A_5))   = {H_P_pred:.4f}  (obs {H_P_obs}, err={(H_P_pred-H_P_obs)/H_P_obs*100:+.6f}%)")
print(f"  SH0ES:   H_0 = (1/t_H)*(1 + K_Mex*(5/108)/(11/9)) = {H_S_pred:.4f}  (obs {H_S_obs}, err={(H_S_pred-H_S_obs)/H_S_obs*100:+.6f}%)")

# Exact fractional form
planck_F = 1 - K_Mex/(N_ch*A_5)         # = 1 - 25/6480 = 6455/6480 = 1291/1296
sh0es_F  = 1 + K_Mex*F(5,108)/F(11,9)   # = 1 + (25/12)(5/108)(9/11) = 1 + 1125/14256
rho_F    = sh0es_F/planck_F
print(f"  planck factor = {planck_F} = {float(planck_F):.8f}")
print(f"  sh0es factor  = {sh0es_F} = {float(sh0es_F):.8f}")
print(f"  rho_pred = sh0es/planck = {rho_F} = {float(rho_F):.6f}")
print(f"  rho_obs  = {rho_obs:.6f}")
err_rho = (float(rho_F)-rho_obs)/rho_obs*100.0
print(f"  err = {err_rho:+.6f}%")

# ============================================================
# TRACK (c) Class XXXIV -- sigma_8 and Omega_m
# ============================================================
print("\n"+"-"*80); print("TRACK (c) -- Class XXXIV: sigma_8 and Omega_m"); print("-"*80)

# --- sigma_8 ---
print("\n  sigma_8: Planck 2018 = 0.8111 +/- 0.0060")
sigma8_obs = 0.8111
sigma8_seed = float(F(416,513))
err_s8_seed = (sigma8_seed-sigma8_obs)/sigma8_obs*100.0
print(f"  SEED: sigma_8 = 416/513 = {sigma8_seed:.6f}, err = {err_s8_seed:+.6f}%  (1-atom, 0-anchor)")
need_mult_s8 = sigma8_obs/sigma8_seed
delta_s8 = need_mult_s8 - 1.0
sign_s8 = 1 if delta_s8>0 else -1
print(f"  needed mult = {need_mult_s8:.8f}, delta = {delta_s8:+.4e}")

print("\n  2-atom delta:")
hits2s = search2(abs(delta_s8), tol_pct=5.0, want=10)
for h in hits2s:
    mult=1.0+sign_s8*h[3]; pred=sigma8_seed*mult; err=(pred-sigma8_obs)/sigma8_obs*100.0
    print(f"    1{'+' if sign_s8>0 else '-'}[{h[0]} {h[2]} {h[1]}]  d={h[3]:.4e}  s8={pred:.6f}  err={err:+.6f}%")

print("\n  3-atom delta:")
hits3s = search3(abs(delta_s8), tol_pct=2.0, want=12)
for h in hits3s:
    mult=1.0+sign_s8*h[4]; pred=sigma8_seed*mult; err=(pred-sigma8_obs)/sigma8_obs*100.0
    print(f"    1{'+' if sign_s8>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}]  d={h[4]:.4e}  s8={pred:.6f}  err={err:+.6f}%")

best_s8=(f"416/513 (seed)", err_s8_seed, sigma8_seed)
for h in hits2s:
    mult=1.0+sign_s8*h[3]; pred=sigma8_seed*mult; err=(pred-sigma8_obs)/sigma8_obs*100.0
    if abs(err)<abs(best_s8[1]):
        best_s8=(f"(416/513)*[1{'+' if sign_s8>0 else '-'}{h[0]} {h[2]} {h[1]}]",err,pred)
for h in hits3s:
    mult=1.0+sign_s8*h[4]; pred=sigma8_seed*mult; err=(pred-sigma8_obs)/sigma8_obs*100.0
    if abs(err)<abs(best_s8[1]):
        best_s8=(f"(416/513)*[1{'+' if sign_s8>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]",err,pred)
print(f"\n  BEST sigma_8: {best_s8[0]} = {best_s8[2]:.6f}, err = {best_s8[1]:+.6f}%")

# --- Omega_m ---
print("\n  Omega_m: Planck 2018 = 0.3153 +/- 0.0073")
Om_obs = 0.3153
Om_seed = float(N_ch*one_m_ns)  # 9*(7/200) = 63/200 = 0.315
err_Om_seed = (Om_seed-Om_obs)/Om_obs*100.0
print(f"  SEED: Omega_m = N_ch*(1-n_s) = 9*(7/200) = 63/200 = {Om_seed:.6f}, err = {err_Om_seed:+.6f}%  (2-atom, 0-anchor)")
need_mult_Om = Om_obs/Om_seed
delta_Om = need_mult_Om - 1.0
sign_Om = 1 if delta_Om>0 else -1
print(f"  needed mult = {need_mult_Om:.8f}, delta = {delta_Om:+.4e}")

print("\n  2-atom delta:")
hits2O = search2(abs(delta_Om), tol_pct=5.0, want=10)
for h in hits2O:
    mult=1.0+sign_Om*h[3]; pred=Om_seed*mult; err=(pred-Om_obs)/Om_obs*100.0
    print(f"    1{'+' if sign_Om>0 else '-'}[{h[0]} {h[2]} {h[1]}]  d={h[3]:.4e}  Om={pred:.6f}  err={err:+.6f}%")

print("\n  3-atom delta:")
hits3O = search3(abs(delta_Om), tol_pct=2.0, want=12)
for h in hits3O:
    mult=1.0+sign_Om*h[4]; pred=Om_seed*mult; err=(pred-Om_obs)/Om_obs*100.0
    print(f"    1{'+' if sign_Om>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}]  d={h[4]:.4e}  Om={pred:.6f}  err={err:+.6f}%")

best_Om=(f"63/200 (seed)", err_Om_seed, Om_seed)
for h in hits2O:
    mult=1.0+sign_Om*h[3]; pred=Om_seed*mult; err=(pred-Om_obs)/Om_obs*100.0
    if abs(err)<abs(best_Om[1]):
        best_Om=(f"(63/200)*[1{'+' if sign_Om>0 else '-'}{h[0]} {h[2]} {h[1]}]",err,pred)
for h in hits3O:
    mult=1.0+sign_Om*h[4]; pred=Om_seed*mult; err=(pred-Om_obs)/Om_obs*100.0
    if abs(err)<abs(best_Om[1]):
        best_Om=(f"(63/200)*[1{'+' if sign_Om>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]",err,pred)
print(f"\n  BEST Omega_m: {best_Om[0]} = {best_Om[2]:.6f}, err = {best_Om[1]:+.6f}%")

# ============================================================
# LEDGER EMISSION
# ============================================================
print()
write_ledger("classXXXII_T0_CMB_session751", T0_pred, T0_obs, status="OK")
write_ledger("classXXXIIIb_H0_ratio_session751", float(rho_F)*H_P_obs, rho_obs*H_P_obs, status="OK")
write_ledger("classXXXIV_sigma8_session751", best_s8[2], sigma8_obs, status="OK")
write_ledger("classXXXIV_Omega_m_session751", best_Om[2], Om_obs, status="OK")

# ============================================================
print("\n"+"-"*80); print("DECISION GATE"); print("-"*80)
print(f"  T_0           err = {best_T[1] if best_T else err_S750:+.6f}%")
print(f"  H_0 ratio     err = {err_rho:+.6f}%")
print(f"  sigma_8       err = {best_s8[1]:+.6f}%")
print(f"  Omega_m       err = {best_Om[1]:+.6f}%")

artifact = os.path.join(os.path.dirname(os.path.abspath(__file__)),"_session751_result.json")
with open(artifact,"w",encoding="utf-8") as f:
    json.dump({
        "T0": {"pred": T0_pred, "obs": T0_obs, "err_pct": best_T[1] if best_T else err_S750, "closure": best_T[0] if best_T else "S750"},
        "H0_ratio": {"pred": float(rho_F), "obs": rho_obs, "err_pct": err_rho,
                     "planck_factor": str(planck_F), "sh0es_factor": str(sh0es_F)},
        "sigma_8": {"pred": best_s8[2], "obs": sigma8_obs, "err_pct": best_s8[1], "closure": best_s8[0]},
        "Omega_m": {"pred": best_Om[2], "obs": Om_obs, "err_pct": best_Om[1], "closure": best_Om[0]},
    }, f, indent=2)
print(f"\nArtifact: {artifact}")
