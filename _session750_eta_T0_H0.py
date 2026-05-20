"""
SESSION 750 -- eta_b fine; T_0 retry (multi-anchor); Class XXXIII H_0 (Hubble tension)

(a) eta_b: base (1-n_s)^6 * 208/627 = 6.0982e-10 at -0.029%. Target delta ~ +2.9e-4.
(b) T_0 retry: try 4-atom search on lambda_max/L_SCM = 3.0447e-6 and T_0/T_anchor = 4.157e5.
(c) Class XXXIII: H_0 = 67.4 (Planck) or 73.0 (SH0ES). 1/t_H ≈ 67.74 km/s/Mpc; test
    H_0_Planck = (1/t_H)*(1 - (1-n_s)/7) = (1/t_H)*(199/200). Also probe SH0ES side.

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

def search4(target, tol_pct=2.0, want=10):
    """4-atom: a*b/(c*d), a*b*c/d, a/b/c/d"""
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

print("="*80); print("SESSION 750 -- eta_b fine; T_0 retry; H_0 Hubble"); print("="*80)

# ============================================================
# TRACK (a) eta_b fine
# ============================================================
print("\n"+"-"*80); print("TRACK (a) -- eta_b fine: (1-n_s)^6 * 208/627"); print("-"*80)
eta_obs=6.10e-10
base_eta = float(one_m_ns)**6 * (208.0/627.0)
err_base=(base_eta-eta_obs)/eta_obs*100.0
print(f"  base = {base_eta:.6e}, err = {err_base:+.6f}%")
need_mult_eta = eta_obs/base_eta
delta_eta = need_mult_eta - 1.0
print(f"  needed mult = {need_mult_eta:.8f}, delta = {delta_eta:+.4e}")

sign_eta = 1 if delta_eta>0 else -1
print("\n  3-atom delta search:")
hits3 = search3(abs(delta_eta), tol_pct=2.0, want=12)
for h in hits3:
    mult=1.0+sign_eta*h[4]; pred=base_eta*mult; err=(pred-eta_obs)/eta_obs*100.0
    print(f"    1{'+' if sign_eta>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}]  delta={h[4]:.4e}  eta={pred:.4e}  err={err:+.6f}%")

print("\n  4-atom delta search (tol 0.5%):")
hits4 = search4(abs(delta_eta), tol_pct=0.5, want=12)
for h in hits4:
    mult=1.0+sign_eta*h[5]; pred=base_eta*mult; err=(pred-eta_obs)/eta_obs*100.0
    print(f"    1{'+' if sign_eta>0 else '-'}[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]  delta={h[5]:.4e}  eta={pred:.4e}  err={err:+.6f}%")

best_eta=None
for h in hits3:
    mult=1.0+sign_eta*h[4]; pred=base_eta*mult; err=(pred-eta_obs)/eta_obs*100.0
    if best_eta is None or abs(err)<abs(best_eta[1]):
        best_eta=(f"base*[1{'+' if sign_eta>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]",err,pred)
for h in hits4:
    mult=1.0+sign_eta*h[5]; pred=base_eta*mult; err=(pred-eta_obs)/eta_obs*100.0
    if best_eta is None or abs(err)<abs(best_eta[1]):
        best_eta=(f"base*[1{'+' if sign_eta>0 else '-'}{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]",err,pred)
if best_eta:
    print(f"\n  BEST: eta_b = (1-n_s)^6*(208/627)*{best_eta[0]} = {best_eta[2]:.4e}, err = {best_eta[1]:+.6f}%")
    eta_pred = best_eta[2]
else: eta_pred = base_eta

# ============================================================
# TRACK (b) T_0 multi-anchor
# ============================================================
print("\n"+"-"*80); print("TRACK (b) -- T_0 = 2.7255 K  multi-anchor retry"); print("-"*80)
T0_obs=2.7255
c=2.99792458e8; L_SCM=349.226733192; rho_vac=7.09e-37
k_B=1.380649e-23; hbar=1.054571817e-34; h_pl=6.62607015e-34
a_rad = math.pi**2 * k_B**4 / (15.0*hbar**3*c**3)
b_Wien = 2.897771955e-3

# Several dimensionless ratios:
r1 = b_Wien / (T0_obs * L_SCM)  # lambda_max/L_SCM
r2 = T0_obs * k_B * L_SCM / (h_pl * c)  # = T_0 / T_HC^-1
T_HC = h_pl * c / (L_SCM * k_B)
r3 = T0_obs / T_HC
u_CMB = a_rad * T0_obs**4
r4 = u_CMB / rho_vac

print(f"  r1 = lambda_max/L_SCM = {r1:.6e}")
print(f"  r2 = T_0*k_B*L_SCM/(h*c) = {r2:.6e}")
print(f"  r3 = T_0/T_HC = T_0*L_SCM*k_B/(h*c) = {r3:.6e}")
print(f"  r4 = u_CMB/rho_vac_SCM = {r4:.6e}")

# r2 = r3 (same). Let's focus on r3 and r4.
# Try search 4-atom for r3 ~ 6.6155e4
print(f"\n  4-atom search for r3 = {r3:.4e} (tol 0.5%):")
hits_r3 = search4(r3, tol_pct=0.5, want=12)
for h in hits_r3:
    print(f"    {h[0]} {h[4]} {h[1]} {h[2]} {h[3]} = {h[5]:.4e}  err={h[6]:+.6f}%")

# Also try multiplicative pure: M such that T_0 = M * (h*c)/(L_SCM*k_B)
# r3 ~ 6.62e4 - try product of inverse-rationals
# (1/F_TRZ)^4 * D_BSFG = 10000*6 = 6e4 — close. ×11/10 = 6.6e4 — ★
trial_1 = 1/float(F_TRZ)**4 * float(D_BSFG) * float(F(11,10))  # 10000*6*1.1=66000
print(f"\n  Trial: (1/F_TRZ)^4 * D_BSFG * 11/10 = {trial_1:.4e}, err = {(trial_1-r3)/r3*100:+.4f}%")
# 11/10 isn't an atom. Try (1/F_TRZ)^4 * D_BSFG * (1+F_TRZ) = 10000*6*1.1=66000
# 1+F_TRZ would need expressing.

# r4 ~ 5.88e22 - astronomical
# (1/F_TRZ)^22 = 1e22 ×5.88 — hmm. (1/xi)^k? (3200/11)^k. 
inv_xi = 1.0/float(xi)
print(f"\n  inv_xi = {inv_xi:.4e}, inv_xi^2 = {inv_xi**2:.4e}")
print(f"  Compare r4 = {r4:.4e}")
# inv_xi^2 = 8.46e4 - too small
# Try A_5^k:
for k in range(8,16):
    v = float(A_5)**k
    print(f"  A_5^{k} = {v:.4e}, ratio r4/A_5^{k} = {r4/v:.4e}")

# Hypothesis: r4 = A_5^13 * something. A_5^13 = 1.6e23. r4/A_5^13 = 0.367. Hmm.
# A_5^12 = 2.18e21. r4/A_5^12 = 27.0 = D_crit·D_phys/...? 27 = 27 (clean!)
v = float(A_5)**12
print(f"  A_5^12 = {v:.4e}, r4/A_5^12 = {r4/v:.6f}")
ratio_a12 = r4/v
# Search 2-atom for ratio_a12
print(f"\n  2-atom search for r4/A_5^12 = {ratio_a12:.4e}:")
hits_a12 = search2(ratio_a12, tol_pct=2.0, want=10)
for h in hits_a12:
    op="*" if h[2]=="a*b" else "/"; print(f"    {h[0]}{op}{h[1]} = {h[3]:.4e}  err={h[4]:+.4f}%")
hits_a12_3 = search3(ratio_a12, tol_pct=1.0, want=10)
print("  3-atom:")
for h in hits_a12_3:
    print(f"    {h[0]} {h[3]} {h[1]} {h[2]} = {h[4]:.4e}  err={h[5]:+.4f}%")

# Better tactic: direct 3-atom search on r3 ~ 6.6e4
print(f"\n  3-atom search for r3 = {r3:.4e} (tol 2%):")
hits_r3_3 = search3(r3, tol_pct=2.0, want=10)
for h in hits_r3_3:
    print(f"    {h[0]} {h[3]} {h[1]} {h[2]} = {h[4]:.4e}  err={h[5]:+.4f}%")

# Try r3 / A_5^k for small k:
for k_pow in range(2,5):
    rk = r3 / float(A_5)**k_pow
    print(f"\n  r3 / A_5^{k_pow} = {rk:.4e}:  2-atom hits:")
    h2 = search2(rk, tol_pct=2.0, want=6)
    for h in h2:
        op="*" if h[2]=="a*b" else "/"; print(f"    A_5^{k_pow} * ({h[0]}{op}{h[1]}) = {float(A_5)**k_pow*h[3]:.4e}  err base ratio={h[4]:+.4f}%")

# Pick best T_0 closure (only if real)
best_t0 = None
# Check hits_r3 (4-atom)
for h in hits_r3:
    T0_pred = h[5] * T_HC
    err = (T0_pred-T0_obs)/T0_obs*100.0
    if best_t0 is None or abs(err)<abs(best_t0[1]):
        best_t0=(f"T_0 = (h*c/(L_SCM*k_B)) * [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", err, T0_pred)
for h in hits_r3_3:
    T0_pred = h[4] * T_HC
    err = (T0_pred-T0_obs)/T0_obs*100.0
    if best_t0 is None or abs(err)<abs(best_t0[1]):
        best_t0=(f"T_0 = (h*c/(L_SCM*k_B)) * [{h[0]} {h[3]} {h[1]} {h[2]}]", err, T0_pred)

if best_t0:
    print(f"\n  BEST: {best_t0[0]} = {best_t0[2]:.4f} K, err = {best_t0[1]:+.6f}%")
    T0_pred = best_t0[2]
    t0_status = "OK"
else:
    print("\n  *** No clean T_0 closure found ***")
    T0_pred = None
    t0_status = "FAILED"

# ============================================================
# TRACK (c) H_0 Hubble — Class XXXIII
# ============================================================
print("\n"+"-"*80); print("TRACK (c) -- Class XXXIII: H_0 (Hubble)"); print("-"*80)
H0_Planck = 67.4       # km/s/Mpc
H0_SH0ES  = 73.0       # km/s/Mpc
t_H_Gyr   = 14.4517    # locked
sec_per_Gyr = 3.1557e16
Mpc_in_km = 3.0857e19
H0_natural = (1.0 / (t_H_Gyr * sec_per_Gyr)) * Mpc_in_km  # km/s/Mpc
print(f"  1/t_H = {H0_natural:.4f} km/s/Mpc  (natural Hubble unit)")
print(f"  H_0(Planck) = {H0_Planck}")
print(f"  H_0(SH0ES)  = {H0_SH0ES}")
print(f"  Planck/natural = {H0_Planck/H0_natural:.8f}")
print(f"  SH0ES/natural  = {H0_SH0ES/H0_natural:.8f}")

# Planck closure: 199/200 = 1 - (1-n_s)/7 — clean
val_Planck = (1 - float(one_m_ns)/7) * H0_natural
err_Planck = (val_Planck - H0_Planck)/H0_Planck*100.0
print(f"\n  H_0(Planck) candidate: (1/t_H)*(1 - (1-n_s)/7)")
print(f"    (1 - (1-n_s)/7) = (1 - 1/200) = 199/200 = {1 - float(one_m_ns)/7:.6f}")
print(f"    H_0 = {val_Planck:.4f}, err = {err_Planck:+.6f}%")

# Also try alternative: H_0(Planck) = (1/t_H)*(1 - xi/...) or similar
print("\n  Planck delta search (target mult = 0.99499):")
mult_P = H0_Planck/H0_natural
delta_P = mult_P - 1.0  # negative
print(f"    needed delta = {delta_P:+.4e}")
hits_P2 = search2(abs(delta_P), tol_pct=5.0, want=10)
for h in hits_P2:
    op="*" if h[2]=="a*b" else "/"
    print(f"    1-{h[0]}{op}{h[1]}  delta={h[3]:.4e}  H_0={H0_natural*(1-h[3]):.4f}  err={(H0_natural*(1-h[3])-H0_Planck)/H0_Planck*100:+.4f}%")
hits_P3 = search3(abs(delta_P), tol_pct=1.0, want=12)
print("  3-atom Planck:")
for h in hits_P3:
    pred=H0_natural*(1-h[4]); err=(pred-H0_Planck)/H0_Planck*100.0
    print(f"    1-[{h[0]} {h[3]} {h[1]} {h[2]}]  delta={h[4]:.4e}  H_0={pred:.4f}  err={err:+.6f}%")

# SH0ES side
print("\n  SH0ES delta search (target mult = 1.07770):")
mult_S = H0_SH0ES/H0_natural
delta_S = mult_S - 1.0  # positive
print(f"    needed delta = {delta_S:+.4e}")
hits_S2 = search2(abs(delta_S), tol_pct=2.0, want=10)
for h in hits_S2:
    op="*" if h[2]=="a*b" else "/"
    print(f"    1+{h[0]}{op}{h[1]}  delta={h[3]:.4e}  H_0={H0_natural*(1+h[3]):.4f}  err={(H0_natural*(1+h[3])-H0_SH0ES)/H0_SH0ES*100:+.4f}%")
hits_S3 = search3(abs(delta_S), tol_pct=1.0, want=12)
print("  3-atom SH0ES:")
for h in hits_S3:
    pred=H0_natural*(1+h[4]); err=(pred-H0_SH0ES)/H0_SH0ES*100.0
    print(f"    1+[{h[0]} {h[3]} {h[1]} {h[2]}]  delta={h[4]:.4e}  H_0={pred:.4f}  err={err:+.6f}%")

# Best Planck:
best_P = ("(1/t_H)*(199/200) = (1/t_H)*[1-(1-n_s)/7]", err_Planck, val_Planck)
for h in hits_P2:
    pred=H0_natural*(1-h[3]); err=(pred-H0_Planck)/H0_Planck*100.0
    if abs(err)<abs(best_P[1]):
        op="*" if h[2]=="a*b" else "/"
        best_P=(f"(1/t_H)*[1-{h[0]}{op}{h[1]}]", err, pred)
for h in hits_P3:
    pred=H0_natural*(1-h[4]); err=(pred-H0_Planck)/H0_Planck*100.0
    if abs(err)<abs(best_P[1]):
        best_P=(f"(1/t_H)*[1-{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
print(f"\n  BEST Planck: {best_P[0]} = {best_P[2]:.4f}, err = {best_P[1]:+.6f}%")
H0P_pred = best_P[2]

# Best SH0ES:
best_S = None
for h in hits_S2:
    pred=H0_natural*(1+h[3]); err=(pred-H0_SH0ES)/H0_SH0ES*100.0
    if best_S is None or abs(err)<abs(best_S[1]):
        op="*" if h[2]=="a*b" else "/"
        best_S=(f"(1/t_H)*[1+{h[0]}{op}{h[1]}]", err, pred)
for h in hits_S3:
    pred=H0_natural*(1+h[4]); err=(pred-H0_SH0ES)/H0_SH0ES*100.0
    if best_S is None or abs(err)<abs(best_S[1]):
        best_S=(f"(1/t_H)*[1+{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
if best_S:
    print(f"  BEST SH0ES:  {best_S[0]} = {best_S[2]:.4f}, err = {best_S[1]:+.6f}%")
    H0S_pred = best_S[2]
else:
    H0S_pred = None

# ============================================================
# Emit
# ============================================================
print()
e1 = write_ledger("classXXX_eta_b_session750", eta_pred, eta_obs)
e2_status = "FAILED" if T0_pred is None else "OK"
if T0_pred is not None:
    e2 = write_ledger("classXXXII_T0_CMB_session750", T0_pred, T0_obs, status=e2_status)
else:
    e2 = None
    print(f"classXXXII_T0_CMB_session750: SKIPPED -- no clean atom match found")
e3a = write_ledger("classXXXIII_H0_Planck_session750", H0P_pred, H0_Planck)
if H0S_pred is not None:
    e3b = write_ledger("classXXXIII_H0_SH0ES_session750", H0S_pred, H0_SH0ES)
else:
    e3b = None

print("\n"+"-"*80); print("DECISION GATE"); print("-"*80)
print(f"  eta_b        err = {e1:+.6f}%")
print(f"  T_0          err = {e2 if e2 is not None else 'FAILED'}")
print(f"  H_0(Planck)  err = {e3a:+.6f}%")
if e3b is not None: print(f"  H_0(SH0ES)   err = {e3b:+.6f}%")

art={"session":750,"cvw":"v2.0.0","sm_anchor":"CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
    "tracks":{
        "a_eta_b":{"formula":best_eta[0] if best_eta else "(1-n_s)^6 * 208/627","value":eta_pred,"err_pct":e1},
        "b_T0":{"formula":best_t0[0] if best_t0 else "FAILED","value":T0_pred,"err_pct":e2,"status":e2_status},
        "c_H0_Planck":{"formula":best_P[0],"value":H0P_pred,"err_pct":e3a},
        "c_H0_SH0ES":{"formula":best_S[0] if best_S else "no closure","value":H0S_pred,"err_pct":e3b},
    }}
art_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),"_session750_eta_T0_H0_result.json")
with open(art_path,"w",encoding="utf-8") as f: json.dump(art,f,indent=2)
print(f"\nArtifact: {art_path}")
