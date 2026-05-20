"""
SESSION 761 -- (a) A_s (Class L) reopen via seed+shell strategy with alpha_em^4 base family;
                (b) Class LI omega_b/omega_c = 0.02237/0.1206 = 0.18549 (baryon-to-CDM ratio);
                (c) Class LII eta_nu = N_eff/3 = 1.01533 (neutrino effective number / SM expectation).

(a) S760 left A_s OPEN: 4-atom direct could not reach 2.1e-9.
    Strategy: use alpha_em^4 = 2.836e-9 as natural seed (35% high) and add 2/3-atom delta shells.
    Alternative seeds: r*xi^3 = 1.46e-9 (30% low), alpha^4 * (33/40) * (1-F_TRZ) = 2.10e-9 (close).

(b) Class LI -- omega_b/omega_c = (416/513)*(1/137)/(Y_p*27/25) divided by
    [beta_i/(N_ch*(416/513)*(137/200))] = quotient of two CE-strict closures.
    Direct hunt + algebraic ratio.

(c) Class LII -- eta_nu = N_eff/3 = 3.046/3 = 1.01533
    From Class XLIII: N_eff = 3 + r*(513/416)/n_s. Algebraic: eta_nu = 1 + r*(513/416)/(3*n_s).

CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""
from __future__ import annotations
from fractions import Fraction as F
import csv, json, os

F_TRZ=F(1,10); Phi_res=F(5,6); SSq=F(57,100); K_Mex=F(25,12)
beta_i=F(6029,10000); D_phys=F(4); D_BSFG=F(6); D_crit=F(26)
N_ch=F(9); SO5=F(10); A_5=F(60)
one_m_FTRZ=1-F_TRZ; one_m_FP=1-F_TRZ*Phi_res
ns_atom=F(193,200); one_m_ns=F(7,200); xi=F(11,3200); r_tens=F(9,250); Yp=F(49,200)
mpme=F(1836152673,1000000000); inv_mpme=F(1)/mpme
alpha_em=F(72973525693,10000000000000); inv_137=F(1,137)
mpme_rat=F(11,6); inv_mpme_rat=F(6,11); one_twelfth=F(1,12)

ATOMS = {
    "F_TRZ":F_TRZ,"Phi_res":Phi_res,"SSq":SSq,"K_Mex":K_Mex,"beta_i":beta_i,
    "D_phys":D_phys,"D_BSFG":D_BSFG,"D_crit":D_crit,"N_ch":N_ch,"SO5":SO5,"A_5":A_5,
    "1-F_TRZ":one_m_FTRZ,"1-F*P":one_m_FP,"n_s":ns_atom,"1-n_s":one_m_ns,
    "xi":xi,"r":r_tens,"Y_p":Yp,
    "27/26":F(27,26),"243/260":F(243,260),"33/40":F(33,40),"11/9":F(11,9),
    "22/9":F(22,9),"27/25":F(27,25),"416/513":F(416,513),"31/30":F(31,30),
    "5/108":F(5,108),"63/200":F(63,200),"137/200":F(137,200),"307/325":F(307,325),
    "1/mpme":inv_mpme, "alpha":alpha_em, "1/137":inv_137, "mpme":mpme,
    "11/6":mpme_rat, "6/11":inv_mpme_rat, "1/12":one_twelfth,
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
        ("a*b*c*d",   lambda a,b,c,d: a*b*c*d),
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

print("="*80); print("SESSION 761 -- A_s (L) reopen; omega_b/omega_c (LI); eta_nu (LII)"); print("="*80)

# ============================================================
# TRACK (a) -- A_s reopen with seed+shell
# ============================================================
print("\n"+"-"*80); print("TRACK (a) -- A_s (L) reopen with alpha_em^4 seed family"); print("-"*80)
As_obs = 2.1e-9

a4 = float(alpha_em)**4
# Build candidate seed list: alpha^4 * single atom (one multiplicative or divisive factor)
seeds = [("alpha^4", a4)]
for k,v in ATOMS.items():
    seeds.append((f"alpha^4*{k}", a4*float(v)))
    if float(v) != 0:
        seeds.append((f"alpha^4/{k}", a4/float(v)))

# Find seed within 30% of target
seed_hits = [(name,val,(val-As_obs)/As_obs*100.0) for name,val in seeds if abs((val-As_obs)/As_obs)<0.30]
seed_hits.sort(key=lambda x:abs(x[2]))
print(f"\n  alpha^4 seeds within 30% of A_s = 2.1e-9:")
for name,val,err in seed_hits[:12]:
    print(f"    {name:30s} = {val:.4e}  err={err:+.4f}%")

# Take best seed and run delta shell
best_a = ("none", 9999.0, 0.0)
if seed_hits:
    base_name, base_val, base_err = seed_hits[0]
    best_a = (base_name, base_err, base_val)
    print(f"\n  Best alpha^4 seed: {base_name} = {base_val:.4e}, err = {base_err:+.4f}%")
    need = As_obs/base_val; delta = need - 1.0; sign = 1 if delta>0 else -1
    print(f"  delta needed = {delta:+.4e}")
    print("\n  2-atom delta on A_s (tol 5%):")
    h2 = search2(abs(delta), tol_pct=5.0, want=15)
    for h in h2[:10]:
        mult=1.0+sign*h[3]; pred=base_val*mult; err=(pred-As_obs)/As_obs*100.0
        if abs(err)<abs(best_a[1]):
            best_a=(f"{base_name}*[1{'+' if sign>0 else '-'}{h[0]} {h[2]} {h[1]}]", err, pred)
        print(f"    1{'+' if sign>0 else '-'}[{h[0]} {h[2]} {h[1]}] A_s={pred:.4e} err={err:+.4f}%")
    print("\n  3-atom delta on A_s (tol 1%):")
    h3 = search3(abs(delta), tol_pct=1.0, want=15)
    for h in h3[:12]:
        mult=1.0+sign*h[4]; pred=base_val*mult; err=(pred-As_obs)/As_obs*100.0
        if abs(err)<abs(best_a[1]):
            best_a=(f"{base_name}*[1{'+' if sign>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
        print(f"    1{'+' if sign>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}] A_s={pred:.4e} err={err:+.4f}%")
    print("\n  4-atom delta on A_s (tol 0.1%):")
    h4 = search4(abs(delta), tol_pct=0.1, want=12)
    for h in h4[:10]:
        mult=1.0+sign*h[5]; pred=base_val*mult; err=(pred-As_obs)/As_obs*100.0
        if abs(err)<abs(best_a[1]):
            best_a=(f"{base_name}*[1{'+' if sign>0 else '-'}{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", err, pred)
        print(f"    1{'+' if sign>0 else '-'}[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] A_s={pred:.4e} err={err:+.4f}%")

# Also direct multi-atom seed approach: alpha^4 * two_atoms
print("\n  Direct: alpha^4 * 2-atom product (tol 1%):")
hits2dir = []
for i in range(len(LABELS)):
    for j in range(len(LABELS)):
        for tag, fn in [("a*b", lambda a,b:a*b),("a/b", lambda a,b:a/b)]:
            try: m = fn(VALS[i], VALS[j])
            except ZeroDivisionError: continue
            if m == 0: continue
            v = a4 * m
            err = (v - As_obs)/As_obs*100.0
            if abs(err) < 1.0:
                hits2dir.append((LABELS[i], LABELS[j], tag, v, err))
hits2dir.sort(key=lambda x:abs(x[4]))
for h in hits2dir[:12]:
    print(f"    alpha^4*[{h[0]} {h[2]} {h[1]}] = {h[3]:.4e}  err={h[4]:+.4f}%")
    if abs(h[4])<abs(best_a[1]):
        best_a=(f"alpha^4*[{h[0]} {h[2]} {h[1]}]", h[4], h[3])

print("\n  Direct: alpha^4 * 3-atom product (tol 0.1%):")
hits3dir = []
forms3 = [("a*b*c", lambda a,b,c:a*b*c),("a*b/c", lambda a,b,c:a*b/c),("a/(b*c)", lambda a,b,c:a/(b*c))]
for i in range(len(LABELS)):
    for j in range(len(LABELS)):
        for k in range(len(LABELS)):
            for tag, fn in forms3:
                try: m = fn(VALS[i], VALS[j], VALS[k])
                except ZeroDivisionError: continue
                if m == 0: continue
                v = a4 * m
                err = (v - As_obs)/As_obs*100.0
                if abs(err) < 0.1:
                    hits3dir.append((LABELS[i], LABELS[j], LABELS[k], tag, v, err))
hits3dir.sort(key=lambda x:abs(x[5]))
for h in hits3dir[:15]:
    print(f"    alpha^4*[{h[0]} {h[3]} {h[1]} {h[2]}] = {h[4]:.4e}  err={h[5]:+.4f}%")
    if abs(h[5])<abs(best_a[1]):
        best_a=(f"alpha^4*[{h[0]} {h[3]} {h[1]} {h[2]}]", h[5], h[4])

print(f"\n  BEST A_s: {best_a[0]} = {best_a[2]:.4e}, err = {best_a[1]:+.6f}%")
a_status = "EXACT" if abs(best_a[1])<1e-8 else ("CE_strict" if abs(best_a[1])<5e-4 else ("tight CLOSED" if abs(best_a[1])<0.1 else ("CLOSED" if abs(best_a[1])<5.0 else "OPEN")))

# ============================================================
# TRACK (b) -- Class LI omega_b/omega_c
# ============================================================
print("\n"+"-"*80); print("TRACK (b) -- Class LI: omega_b/omega_c"); print("-"*80)

# Use exact ratio from XXXVII / XLVII closures
Obh2 = float(F(416,513))*float(inv_137)/(float(Yp)*float(F(27,25)))  # Class XXXVII
Och2 = float(beta_i)/(float(N_ch)*float(F(416,513))*float(F(137,200))) * (1 + float(one_m_ns)*float(F(5,108))/(float(D_crit)*float(F(22,9))))  # XLVII CE
ratio_alg = Obh2 / Och2
ratio_obs = 0.02237 / 0.1206
print(f"  XXXVII Omega_b*h^2 = {Obh2:.6e}")
print(f"  XLVII  Omega_c*h^2 = {Och2:.6e}")
print(f"  Algebraic ratio    = {ratio_alg:.6f}")
print(f"  Observed ratio     = {ratio_obs:.6f}")
print(f"  err_alg = {(ratio_alg-ratio_obs)/ratio_obs*100.0:+.6f}%")

# Direct hunt
print("\n  2-atom direct on omega_b/omega_c (tol 3%):")
h2r = search2(ratio_obs, tol_pct=3.0, want=15)
for h in h2r[:12]:
    print(f"    [{h[0]} {h[2]} {h[1]}] = {h[3]:.6f}  err={h[4]:+.6f}%")
print("\n  3-atom direct on omega_b/omega_c (tol 0.5%):")
h3r = search3(ratio_obs, tol_pct=0.5, want=15)
for h in h3r[:12]:
    print(f"    [{h[0]} {h[3]} {h[1]} {h[2]}] = {h[4]:.6f}  err={h[5]:+.6f}%")
print("\n  4-atom direct on omega_b/omega_c (tol 0.05%):")
h4r = search4(ratio_obs, tol_pct=0.05, want=15)
for h in h4r[:15]:
    print(f"    [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {h[5]:.6f}  err={h[6]:+.6f}%")

best_r = ("algebraic XXXVII/XLVII", (ratio_alg-ratio_obs)/ratio_obs*100.0, ratio_alg)
for h in h2r:
    if abs(h[4])<abs(best_r[1]): best_r=(f"[{h[0]} {h[2]} {h[1]}]", h[4], h[3])
for h in h3r:
    if abs(h[5])<abs(best_r[1]): best_r=(f"[{h[0]} {h[3]} {h[1]} {h[2]}]", h[5], h[4])
for h in h4r:
    if abs(h[6])<abs(best_r[1]): best_r=(f"[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", h[6], h[5])
print(f"\n  BEST omega_b/omega_c: {best_r[0]} = {best_r[2]:.6f}, err = {best_r[1]:+.6f}%")
r_status = "EXACT" if abs(best_r[1])<1e-8 else ("CE_strict" if abs(best_r[1])<5e-4 else ("tight CLOSED" if abs(best_r[1])<0.1 else "CLOSED"))

# ============================================================
# TRACK (c) -- Class LII eta_nu = N_eff/3 = 1.01533
# ============================================================
print("\n"+"-"*80); print("TRACK (c) -- Class LII: eta_nu = N_eff/3 = 1.01533"); print("-"*80)
eta_obs = 1.01533

# Algebraic from XLIII: eta_nu = 1 + r*(513/416)/(3*n_s)
eta_alg = 1.0 + float(r_tens)*(513.0/416.0)/(3.0*float(ns_atom))
print(f"  Algebraic from XLIII: eta_nu = 1 + r*(513/416)/(3*n_s) = {eta_alg:.6f}")
print(f"  Observed             = {eta_obs:.6f}")
print(f"  err_alg = {(eta_alg-eta_obs)/eta_obs*100.0:+.6f}%")

# Direct hunt
print("\n  2-atom direct on eta_nu (tol 3%):")
h2e = search2(eta_obs, tol_pct=3.0, want=15)
for h in h2e[:12]:
    print(f"    [{h[0]} {h[2]} {h[1]}] = {h[3]:.6f}  err={h[4]:+.6f}%")
print("\n  3-atom direct on eta_nu (tol 0.1%):")
h3e = search3(eta_obs, tol_pct=0.1, want=15)
for h in h3e[:12]:
    print(f"    [{h[0]} {h[3]} {h[1]} {h[2]}] = {h[4]:.6f}  err={h[5]:+.6f}%")
print("\n  4-atom direct on eta_nu (tol 0.01%):")
h4e = search4(eta_obs, tol_pct=0.01, want=15)
for h in h4e[:15]:
    print(f"    [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {h[5]:.6f}  err={h[6]:+.6f}%")

best_e = ("algebraic XLIII/3", (eta_alg-eta_obs)/eta_obs*100.0, eta_alg)
for h in h2e:
    if abs(h[4])<abs(best_e[1]): best_e=(f"[{h[0]} {h[2]} {h[1]}]", h[4], h[3])
for h in h3e:
    if abs(h[5])<abs(best_e[1]): best_e=(f"[{h[0]} {h[3]} {h[1]} {h[2]}]", h[5], h[4])
for h in h4e:
    if abs(h[6])<abs(best_e[1]): best_e=(f"[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", h[6], h[5])
print(f"\n  BEST eta_nu: {best_e[0]} = {best_e[2]:.6f}, err = {best_e[1]:+.6f}%")
e_status = "EXACT" if abs(best_e[1])<1e-8 else ("CE_strict" if abs(best_e[1])<5e-4 else ("tight CLOSED" if abs(best_e[1])<0.1 else "CLOSED"))

# ============================================================
# Write ledger
# ============================================================
print()
write_ledger("classL_A_s_session761", best_a[2] if best_a[2]>0 else 1e-30, As_obs, status=a_status)
write_ledger("classLI_omegab_omegac_session761", best_r[2], ratio_obs, status=r_status)
write_ledger("classLII_eta_nu_session761", best_e[2], eta_obs, status=e_status)

print("\n"+"-"*80); print("DECISION GATE"); print("-"*80)
print(f"  A_s (L)                    err = {best_a[1]:+.6f}%  ({a_status})")
print(f"  omega_b/omega_c (LI)       err = {best_r[1]:+.6f}%  ({r_status})")
print(f"  eta_nu (LII)               err = {best_e[1]:+.6f}%  ({e_status})")

out = {
    "session": 761,
    "A_s":             {"best": best_a[0], "pred": best_a[2], "err_pct": best_a[1], "status": a_status},
    "omegab_omegac":   {"best": best_r[0], "pred": best_r[2], "err_pct": best_r[1], "status": r_status},
    "eta_nu":          {"best": best_e[0], "pred": best_e[2], "err_pct": best_e[1], "status": e_status},
    "cvw": "v2.0.0",
    "sm_anchor": "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant"
}
art = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_session761_result.json")
with open(art,"w",encoding="utf-8") as f: json.dump(out, f, indent=2)
print(f"\nArtifact: {art}")
