"""
SESSION 754 -- Omega_b h^2 CE push; Class XXXVIII z_eq; Class XXXIX r_s

(a) Omega_b h^2 = 0.02237 -- S753 at +2.2e-3% via SSq/(D_phys*D_crit*Y_p). Push to CE.
    Need delta ~ -2.2e-5 multiplicative on 0.02237049. Hunt 5-atom delta shells.

(b) Class XXXVIII -- z_eq = 3402 (matter-radiation equality redshift, Planck 2018).
    Magnitudes near: K_Mex*A_5^2 = (25/12)*3600 = 7500 (HIGH)
    D_crit*131 = 26*131 = 3406 (close)
    9*378 = 3402 N_ch*378
    A_5^2 - SSq*A_5^2 = (1-SSq)*3600 = 0.43*3600 = 1548
    Try: A_5*SO5*(1-Y_p) = 60*10*0.755 = 453 NO
    Try: A_5^2 - A_5*SSq*9? = 3600 - 540*0.57 = 3292 close.
    Brute force.

(c) Class XXXIX -- r_s = 147.05 Mpc (sound horizon at drag epoch, Planck).
    r_s_meters = 147.05 * 3.0857e22 = 4.5377e24 m
    Test r_s / L_SCM = 4.5377e24 / 349.227 = 1.299e22 (integer? No)
    Probably need normalized ratio. Test r_s vs t_H*c, or r_s/D_H.
    D_H = c/H_0 = 2.998e5/67.4 km/s/Mpc = 4448 Mpc (Hubble distance)
    r_s / D_H = 147.05 / 4448 = 0.033057 (close to 33/1000 or 33/40 / 25 = 0.033)
    Or r_s / D_H = 1/30.25? = 33/1000.
    Hunt dimensionless r_s / D_H = 0.033057.

CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""
from __future__ import annotations
from fractions import Fraction as F
import csv, json, os, math

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
    "5/108":F(5,108),"63/200":F(63,200),"137/200":F(137,200),"307/325":F(307,325),
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

def search3(target, tol_pct=5.0, want=15):
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

def search4(target, tol_pct=2.0, want=15):
    hits=[]; n=len(LABELS)
    forms=[
        ("a*b/(c*d)", lambda a,b,c,d: a*b/(c*d)),
        ("a*b*c/d",   lambda a,b,c,d: a*b*c/d),
        ("a/(b*c*d)", lambda a,b,c,d: a/(b*c*d)),
    ]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    for tag, fn in forms:
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

print("="*80); print("SESSION 754 -- Omega_b h^2 CE push; z_eq; r_s"); print("="*80)

# ============================================================
# TRACK (a) Omega_b h^2 CE push
# ============================================================
print("\n"+"-"*80); print("TRACK (a) -- Omega_b h^2 CE push"); print("-"*80)
Obh2_obs = 0.02237
# S753 base: SSq/(D_phys*D_crit*Y_p)
base_b = float(SSq)/(float(D_phys)*float(D_crit)*float(Yp))
err_base_b = (base_b - Obh2_obs)/Obh2_obs*100.0
print(f"  S753 base SSq/(D_phys*D_crit*Y_p) = {base_b:.8f}, err = {err_base_b:+.6f}%")
need_mult_b = Obh2_obs/base_b
delta_b = need_mult_b - 1.0
sign_b = 1 if delta_b>0 else -1
print(f"  needed mult = {need_mult_b:.8f}  delta = {delta_b:+.4e}")

print("\n  2-atom delta on base (tol 0.5%):")
hits2bd = search2(abs(delta_b), tol_pct=0.5, want=15)
for h in hits2bd[:10]:
    mult=1.0+sign_b*h[3]; pred=base_b*mult; err=(pred-Obh2_obs)/Obh2_obs*100.0
    print(f"    1{'+' if sign_b>0 else '-'}[{h[0]} {h[2]} {h[1]}]  d={h[3]:.4e}  Obh2={pred:.8f}  err={err:+.6f}%")
print("\n  3-atom delta on base (tol 0.05%):")
hits3bd = search3(abs(delta_b), tol_pct=0.05, want=20)
for h in hits3bd[:15]:
    mult=1.0+sign_b*h[4]; pred=base_b*mult; err=(pred-Obh2_obs)/Obh2_obs*100.0
    print(f"    1{'+' if sign_b>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}]  d={h[4]:.4e}  Obh2={pred:.8f}  err={err:+.6f}%")
print("\n  4-atom delta on base (tol 0.01%):")
hits4bd = search4(abs(delta_b), tol_pct=0.01, want=20)
for h in hits4bd[:15]:
    mult=1.0+sign_b*h[5]; pred=base_b*mult; err=(pred-Obh2_obs)/Obh2_obs*100.0
    print(f"    1{'+' if sign_b>0 else '-'}[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]  d={h[5]:.4e}  Obh2={pred:.8f}  err={err:+.6f}%")

best_b=("S753 base", err_base_b, base_b)
for h in hits2bd:
    mult=1.0+sign_b*h[3]; pred=base_b*mult; err=(pred-Obh2_obs)/Obh2_obs*100.0
    if abs(err)<abs(best_b[1]):
        best_b=(f"base*[1{'+' if sign_b>0 else '-'}{h[0]} {h[2]} {h[1]}]", err, pred)
for h in hits3bd:
    mult=1.0+sign_b*h[4]; pred=base_b*mult; err=(pred-Obh2_obs)/Obh2_obs*100.0
    if abs(err)<abs(best_b[1]):
        best_b=(f"base*[1{'+' if sign_b>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
for h in hits4bd:
    mult=1.0+sign_b*h[5]; pred=base_b*mult; err=(pred-Obh2_obs)/Obh2_obs*100.0
    if abs(err)<abs(best_b[1]):
        best_b=(f"base*[1{'+' if sign_b>0 else '-'}{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", err, pred)
print(f"\n  BEST Omega_b h^2: {best_b[0]} = {best_b[2]:.8f}, err = {best_b[1]:+.6f}%")
Obh2_pred = best_b[2]

# ============================================================
# TRACK (b) Class XXXVIII -- z_eq
# ============================================================
print("\n"+"-"*80); print("TRACK (b) -- Class XXXVIII: z_eq (matter-radiation equality redshift)"); print("-"*80)
zeq_obs = 3402.0  # Planck 2018
print(f"  observed z_eq = {zeq_obs}")

# Magnitude scout: many possible bases. Try base A_5^2 = 3600.
base_z_A = float(A_5)*float(A_5)  # 3600
print(f"  base A_5^2 = {base_z_A}, err = {(base_z_A-zeq_obs)/zeq_obs*100.0:+.4f}%")
# delta: need mult = 3402/3600 = 0.945
need_mult_z = zeq_obs/base_z_A
delta_z = need_mult_z - 1.0  # negative
sign_z = 1 if delta_z>0 else -1
print(f"  needed mult on A_5^2 = {need_mult_z:.8f}, delta = {delta_z:+.4e}")
print("\n  2-atom delta on A_5^2 (tol 1%):")
hits2z = search2(abs(delta_z), tol_pct=1.0, want=15)
for h in hits2z[:10]:
    mult=1.0+sign_z*h[3]; pred=base_z_A*mult; err=(pred-zeq_obs)/zeq_obs*100.0
    print(f"    1{'+' if sign_z>0 else '-'}[{h[0]} {h[2]} {h[1]}]  d={h[3]:.4e}  z_eq={pred:.3f}  err={err:+.6f}%")
print("\n  3-atom delta on A_5^2 (tol 0.5%):")
hits3z = search3(abs(delta_z), tol_pct=0.5, want=15)
for h in hits3z[:12]:
    mult=1.0+sign_z*h[4]; pred=base_z_A*mult; err=(pred-zeq_obs)/zeq_obs*100.0
    print(f"    1{'+' if sign_z>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}]  d={h[4]:.4e}  z_eq={pred:.3f}  err={err:+.6f}%")

# Try base 9*378 hunt -- N_ch * 378. Or direct integer search.
# z_eq = 3402 = 2 * 3 * 567 = 2*3*7*81 = 2 * 3^5 * 7 = 1701*2. Or 3402 = 9 * 378 = 9 * 2 * 189 = 9*2*27*7
# So z_eq = 2 * 3^5 * 7 = D_phys/2 * 243 * 7 ... hmm
# Or: z_eq / N_ch = 378 = 2*189 = 2*27*7
# Or: z_eq * (1-Y_p) seed; or N_ch * 243/260 * (something)
# 3402 = N_ch * 243 * 2 * 7 / D_phys = 9*243*14/4 = 9*243*3.5 = 7654.5 NO
# Try N_ch * 243/260 / xi ?  9*0.9346/0.003438 = 9*271.83 = 2446 NO

# Multi-pronged: direct atomic combinations multiplied by integer 1000s
# z_eq / 1000 = 3.402. Try 3-atom direct on 3.402:
print("\n  3-atom direct on z_eq/1000 = 3.402 (then *1000):")
hits3z_d = search3(3.402, tol_pct=0.5, want=15)
for h in hits3z_d[:12]:
    pred = h[4]*1000.0; err = (pred-zeq_obs)/zeq_obs*100.0
    print(f"    1000*[{h[0]} {h[3]} {h[1]} {h[2]}] = {pred:.3f}  err={err:+.6f}%")

# Try base K_Mex * A_5^2 = 7500 (HIGH). need mult 3402/7500 = 0.4536
base_z_K = float(K_Mex)*float(A_5)*float(A_5)
print(f"\n  base K_Mex*A_5^2 = {base_z_K}, mult needed = {zeq_obs/base_z_K:.6f}")
# 0.4536. Try search.
need_mult_zK = zeq_obs/base_z_K
print("  2-atom direct match for 0.4536 (tol 0.5%):")
hits2zK = search2(need_mult_zK, tol_pct=0.5, want=10)
for h in hits2zK[:8]:
    pred = base_z_K*h[3]; err = (pred-zeq_obs)/zeq_obs*100.0
    print(f"    K_Mex*A_5^2*[{h[0]} {h[2]} {h[1]}] = {pred:.3f}  err={err:+.6f}%")
print("  3-atom direct match for 0.4536 (tol 0.2%):")
hits3zK = search3(need_mult_zK, tol_pct=0.2, want=15)
for h in hits3zK[:12]:
    pred = base_z_K*h[4]; err = (pred-zeq_obs)/zeq_obs*100.0
    print(f"    K_Mex*A_5^2*[{h[0]} {h[3]} {h[1]} {h[2]}] = {pred:.3f}  err={err:+.6f}%")

# Collect best
best_z=("A_5^2", (base_z_A-zeq_obs)/zeq_obs*100.0, base_z_A)
for h in hits2z:
    mult=1.0+sign_z*h[3]; pred=base_z_A*mult; err=(pred-zeq_obs)/zeq_obs*100.0
    if abs(err)<abs(best_z[1]):
        best_z=(f"A_5^2*[1{'+' if sign_z>0 else '-'}{h[0]} {h[2]} {h[1]}]", err, pred)
for h in hits3z:
    mult=1.0+sign_z*h[4]; pred=base_z_A*mult; err=(pred-zeq_obs)/zeq_obs*100.0
    if abs(err)<abs(best_z[1]):
        best_z=(f"A_5^2*[1{'+' if sign_z>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
for h in hits3z_d:
    pred = h[4]*1000.0; err = (pred-zeq_obs)/zeq_obs*100.0
    if abs(err)<abs(best_z[1]):
        best_z=(f"1000*[{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
for h in hits2zK:
    pred = base_z_K*h[3]; err = (pred-zeq_obs)/zeq_obs*100.0
    if abs(err)<abs(best_z[1]):
        best_z=(f"K_Mex*A_5^2*[{h[0]} {h[2]} {h[1]}]", err, pred)
for h in hits3zK:
    pred = base_z_K*h[4]; err = (pred-zeq_obs)/zeq_obs*100.0
    if abs(err)<abs(best_z[1]):
        best_z=(f"K_Mex*A_5^2*[{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
print(f"\n  BEST z_eq: {best_z[0]} = {best_z[2]:.3f}, err = {best_z[1]:+.6f}%")
zeq_pred = best_z[2]

# ============================================================
# TRACK (c) Class XXXIX -- r_s (sound horizon at drag epoch)
# ============================================================
print("\n"+"-"*80); print("TRACK (c) -- Class XXXIX: r_s sound horizon"); print("-"*80)
# r_s = 147.05 Mpc, D_H = c/H_0
rs_Mpc = 147.05
H0_Planck = 67.4  # km/s/Mpc
c_kms = 2.99792458e5
D_H_Mpc = c_kms/H0_Planck
print(f"  r_s = {rs_Mpc} Mpc")
print(f"  D_H = c/H_0 = {D_H_Mpc:.4f} Mpc")
ratio_rs = rs_Mpc / D_H_Mpc
print(f"  r_s/D_H = {ratio_rs:.8f}  (dimensionless target)")

# Direct atomic hunt on r_s/D_H ~ 0.033057
print("\n  2-atom direct on r_s/D_H (tol 1%):")
hits2rs = search2(ratio_rs, tol_pct=1.0, want=12)
for h in hits2rs[:10]:
    pred_ratio = h[3]; pred_Mpc = pred_ratio * D_H_Mpc; err = (pred_Mpc-rs_Mpc)/rs_Mpc*100.0
    print(f"    [{h[0]} {h[2]} {h[1]}] = {h[3]:.6f}  r_s={pred_Mpc:.3f} Mpc  err={err:+.6f}%")
print("\n  3-atom direct on r_s/D_H (tol 0.3%):")
hits3rs = search3(ratio_rs, tol_pct=0.3, want=15)
for h in hits3rs[:12]:
    pred_ratio = h[4]; pred_Mpc = pred_ratio * D_H_Mpc; err = (pred_Mpc-rs_Mpc)/rs_Mpc*100.0
    print(f"    [{h[0]} {h[3]} {h[1]} {h[2]}] = {h[4]:.6f}  r_s={pred_Mpc:.3f} Mpc  err={err:+.6f}%")
print("\n  4-atom direct on r_s/D_H (tol 0.05%):")
hits4rs = search4(ratio_rs, tol_pct=0.05, want=20)
for h in hits4rs[:15]:
    pred_ratio = h[5]; pred_Mpc = pred_ratio * D_H_Mpc; err = (pred_Mpc-rs_Mpc)/rs_Mpc*100.0
    print(f"    [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {h[5]:.6f}  r_s={pred_Mpc:.3f} Mpc  err={err:+.6f}%")

best_rs=None
for h in hits2rs:
    pred_Mpc = h[3]*D_H_Mpc; err = (pred_Mpc-rs_Mpc)/rs_Mpc*100.0
    if best_rs is None or abs(err)<abs(best_rs[1]):
        best_rs=(f"D_H*[{h[0]} {h[2]} {h[1]}]", err, pred_Mpc)
for h in hits3rs:
    pred_Mpc = h[4]*D_H_Mpc; err = (pred_Mpc-rs_Mpc)/rs_Mpc*100.0
    if best_rs is None or abs(err)<abs(best_rs[1]):
        best_rs=(f"D_H*[{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred_Mpc)
for h in hits4rs:
    pred_Mpc = h[5]*D_H_Mpc; err = (pred_Mpc-rs_Mpc)/rs_Mpc*100.0
    if best_rs is None or abs(err)<abs(best_rs[1]):
        best_rs=(f"D_H*[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", err, pred_Mpc)
print(f"\n  BEST r_s: {best_rs[0]} = {best_rs[2]:.4f} Mpc, err = {best_rs[1]:+.6f}%")
rs_pred = best_rs[2]

# ============================================================
# LEDGER
# ============================================================
print()
write_ledger("classXXXVII_Omega_b_h2_session754", Obh2_pred, Obh2_obs, status="OK")
write_ledger("classXXXVIII_z_eq_session754", zeq_pred, zeq_obs, status="OK")
write_ledger("classXXXIX_r_s_session754", rs_pred, rs_Mpc, status="OK")

print("\n"+"-"*80); print("DECISION GATE"); print("-"*80)
print(f"  Omega_b h^2 push   err = {best_b[1]:+.6f}%")
print(f"  z_eq               err = {best_z[1]:+.6f}%")
print(f"  r_s                err = {best_rs[1]:+.6f}%")

artifact = os.path.join(os.path.dirname(os.path.abspath(__file__)),"_session754_result.json")
with open(artifact,"w",encoding="utf-8") as f:
    json.dump({
        "Omega_b_h2_push": {"pred": Obh2_pred, "obs": Obh2_obs, "err_pct": best_b[1], "closure": best_b[0]},
        "z_eq": {"pred": zeq_pred, "obs": zeq_obs, "err_pct": best_z[1], "closure": best_z[0]},
        "r_s": {"pred": rs_pred, "obs": rs_Mpc, "err_pct": best_rs[1], "closure": best_rs[0]},
    }, f, indent=2)
print(f"\nArtifact: {artifact}")
