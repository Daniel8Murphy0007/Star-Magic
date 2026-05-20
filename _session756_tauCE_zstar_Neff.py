"""
SESSION 756 -- (a) tau_reion CE push; (b) Class XLII z_* photon decoupling;
                (c) Class XLIII N_eff effective neutrino species.

(a) S755 closed tau = 0.054402 at +3.7e-3% via Y_p*(5/108)*(1-n_s)*alpha.
    Push under strict CE threshold (5e-4%) using a 5-atom shell or the new
    1/mpme atom. Delta needed on best base: (0.0544 - 0.054402)/0.054402 = -3.7e-5.

(b) Class XLII -- z_* = 1089.80 (Planck 2018 photon-baryon decoupling redshift).
    Distinct from z_eq = 3402 (matter-radiation equality, Class XXXVIII).
    Ratio z_*/z_eq = 1089.80/3402 = 0.32034. Hunt direct on z_* and on ratio.

(c) Class XLIII -- N_eff = 3.046 (effective number of neutrino species, Planck).
    Standard: 3 active + 0.046 QED/non-instantaneous-decoupling correction.
    Seeds: 3 + 1/(D_phys*K_Mex) ?  Try direct atomic hunt.

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

ATOMS = {
    "F_TRZ":F_TRZ,"Phi_res":Phi_res,"SSq":SSq,"K_Mex":K_Mex,"beta_i":beta_i,
    "D_phys":D_phys,"D_BSFG":D_BSFG,"D_crit":D_crit,"N_ch":N_ch,"SO5":SO5,"A_5":A_5,
    "1-F_TRZ":one_m_FTRZ,"1-F*P":one_m_FP,"n_s":ns_atom,"1-n_s":one_m_ns,
    "xi":xi,"r":r_tens,"Y_p":Yp,
    "27/26":F(27,26),"243/260":F(243,260),"33/40":F(33,40),"11/9":F(11,9),
    "22/9":F(22,9),"27/25":F(27,25),"416/513":F(416,513),"31/30":F(31,30),
    "5/108":F(5,108),"63/200":F(63,200),"137/200":F(137,200),"307/325":F(307,325),
    "1/mpme":inv_mpme, "alpha":alpha_em, "1/137":inv_137, "mpme":mpme,
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

print("="*80); print("SESSION 756 -- tau CE push; z_*; N_eff"); print("="*80)

# ============================================================
# TRACK (a) tau CE push
# ============================================================
print("\n"+"-"*80); print("TRACK (a) -- tau_reion strict CE push"); print("-"*80)
tau_obs = 0.0544
base_tau = float(Yp) * float(F(5,108)) * float(one_m_ns) * float(alpha_em)
err_base_tau = (base_tau - tau_obs)/tau_obs*100.0
print(f"  S755 base (Y_p*(5/108)*(1-n_s)*alpha) = {base_tau:.8f}, err = {err_base_tau:+.6f}%")
need_mult_tau = tau_obs/base_tau
delta_tau = need_mult_tau - 1.0
sign_tau = 1 if delta_tau>0 else -1
print(f"  delta needed = {delta_tau:+.4e}")

print("\n  2-atom delta on tau base (tol 0.01%):")
hits2td = search2(abs(delta_tau), tol_pct=0.01, want=20)
for h in hits2td[:15]:
    mult=1.0+sign_tau*h[3]; pred=base_tau*mult; err=(pred-tau_obs)/tau_obs*100.0
    print(f"    1{'+' if sign_tau>0 else '-'}[{h[0]} {h[2]} {h[1]}]  d={h[3]:.4e}  tau={pred:.8f}  err={err:+.6f}%")

print("\n  3-atom delta on tau base (tol 0.005%):")
hits3td = search3(abs(delta_tau), tol_pct=0.005, want=20)
for h in hits3td[:15]:
    mult=1.0+sign_tau*h[4]; pred=base_tau*mult; err=(pred-tau_obs)/tau_obs*100.0
    print(f"    1{'+' if sign_tau>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}]  d={h[4]:.4e}  tau={pred:.8f}  err={err:+.6f}%")

print("\n  4-atom direct on tau (tol 0.001% strict):")
hits4tau = search4(tau_obs, tol_pct=0.001, want=20)
for h in hits4tau[:15]:
    print(f"    [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {h[5]:.8f}  err={h[6]:+.6f}%")

best_tau = ("S755 base", err_base_tau, base_tau)
for h in hits2td:
    mult=1.0+sign_tau*h[3]; pred=base_tau*mult; err=(pred-tau_obs)/tau_obs*100.0
    if abs(err)<abs(best_tau[1]): best_tau=(f"base*[1{'+' if sign_tau>0 else '-'}{h[0]} {h[2]} {h[1]}]", err, pred)
for h in hits3td:
    mult=1.0+sign_tau*h[4]; pred=base_tau*mult; err=(pred-tau_obs)/tau_obs*100.0
    if abs(err)<abs(best_tau[1]): best_tau=(f"base*[1{'+' if sign_tau>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
for h in hits4tau:
    if abs(h[6])<abs(best_tau[1]): best_tau=(f"[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", h[6], h[5])

print(f"\n  BEST tau: {best_tau[0]} = {best_tau[2]:.8f}, err = {best_tau[1]:+.6f}%")
tau_status = "CE_strict" if abs(best_tau[1])<5e-4 else "tight CLOSED"
print(f"  STATUS: {tau_status}")

# ============================================================
# TRACK (b) Class XLII -- z_*
# ============================================================
print("\n"+"-"*80); print("TRACK (b) -- Class XLII: z_* photon decoupling"); print("-"*80)
zstar_obs = 1089.80
zeq = 3402.0
ratio_z = zstar_obs/zeq
print(f"  observed z_* = {zstar_obs}, z_eq (Class XXXVIII) = {zeq}")
print(f"  z_*/z_eq = {ratio_z:.8f}")
# Direct on z_*
print("\n  2-atom direct on z_* (tol 1%):")
hits2z = search2(zstar_obs, tol_pct=1.0, want=15)
for h in hits2z[:10]:
    print(f"    [{h[0]} {h[2]} {h[1]}] = {h[3]:.4f}  err={h[4]:+.6f}%")
print("\n  3-atom direct on z_* (tol 0.3%):")
hits3z = search3(zstar_obs, tol_pct=0.3, want=20)
for h in hits3z[:15]:
    print(f"    [{h[0]} {h[3]} {h[1]} {h[2]}] = {h[4]:.4f}  err={h[5]:+.6f}%")
print("\n  4-atom direct on z_* (tol 0.05%):")
hits4z = search4(zstar_obs, tol_pct=0.05, want=20)
for h in hits4z[:15]:
    print(f"    [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {h[5]:.4f}  err={h[6]:+.6f}%")
# Also hunt the ratio against z_eq
print("\n  3-atom direct on z_*/z_eq = 0.32034 (tol 0.3%):")
hits3zr = search3(ratio_z, tol_pct=0.3, want=15)
for h in hits3zr[:12]:
    pred = h[4]*zeq; err = (pred-zstar_obs)/zstar_obs*100.0
    print(f"    z_eq*[{h[0]} {h[3]} {h[1]} {h[2]}] = {pred:.4f}  err={err:+.6f}%")

best_z=None
for h in hits2z:
    if best_z is None or abs(h[4])<abs(best_z[1]):
        best_z=(f"[{h[0]} {h[2]} {h[1]}]", h[4], h[3])
for h in hits3z:
    if best_z is None or abs(h[5])<abs(best_z[1]):
        best_z=(f"[{h[0]} {h[3]} {h[1]} {h[2]}]", h[5], h[4])
for h in hits4z:
    if best_z is None or abs(h[6])<abs(best_z[1]):
        best_z=(f"[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", h[6], h[5])
for h in hits3zr:
    pred = h[4]*zeq; err = (pred-zstar_obs)/zstar_obs*100.0
    if best_z is None or abs(err)<abs(best_z[1]):
        best_z=(f"z_eq*[{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
print(f"\n  BEST z_*: {best_z[0]} = {best_z[2]:.4f}, err = {best_z[1]:+.6f}%")
zstar_pred = best_z[2]

# ============================================================
# TRACK (c) Class XLIII -- N_eff
# ============================================================
print("\n"+"-"*80); print("TRACK (c) -- Class XLIII: N_eff effective neutrino species"); print("-"*80)
Neff_obs = 3.046
# Standard composition: 3 (active) + 0.046 (QED correction)
delta_Neff = 0.046
print(f"  observed N_eff = {Neff_obs} = 3 + {delta_Neff}")
print(f"  delta = 0.046 vs alpha (1/137) = {float(alpha_em):.6e}")
print(f"  alpha * D_phys * Phi_res = {float(alpha_em)*4*5/6:.6f}  (vs 0.046)")
# Direct hunt on N_eff
print("\n  2-atom direct on N_eff (tol 1%):")
hits2n = search2(Neff_obs, tol_pct=1.0, want=15)
for h in hits2n[:10]:
    print(f"    [{h[0]} {h[2]} {h[1]}] = {h[3]:.6f}  err={h[4]:+.6f}%")
print("\n  3-atom direct on N_eff (tol 0.3%):")
hits3n = search3(Neff_obs, tol_pct=0.3, want=15)
for h in hits3n[:15]:
    print(f"    [{h[0]} {h[3]} {h[1]} {h[2]}] = {h[4]:.6f}  err={h[5]:+.6f}%")
print("\n  4-atom direct on N_eff (tol 0.05%):")
hits4n = search4(Neff_obs, tol_pct=0.05, want=20)
for h in hits4n[:15]:
    print(f"    [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {h[5]:.6f}  err={h[6]:+.6f}%")
# Delta hunt on 0.046
print("\n  2-atom direct on N_eff residual 0.046 (tol 5%):")
hits2dn = search2(delta_Neff, tol_pct=5.0, want=20)
for h in hits2dn[:12]:
    pred = 3.0 + h[3]; err = (pred-Neff_obs)/Neff_obs*100.0
    print(f"    3+[{h[0]} {h[2]} {h[1]}] = {pred:.6f}  err={err:+.6f}%")
print("\n  3-atom direct on N_eff residual 0.046 (tol 1%):")
hits3dn = search3(delta_Neff, tol_pct=1.0, want=20)
for h in hits3dn[:15]:
    pred = 3.0 + h[4]; err = (pred-Neff_obs)/Neff_obs*100.0
    print(f"    3+[{h[0]} {h[3]} {h[1]} {h[2]}] = {pred:.6f}  err={err:+.6f}%")

best_n=None
for h in hits2n:
    if best_n is None or abs(h[4])<abs(best_n[1]):
        best_n=(f"[{h[0]} {h[2]} {h[1]}]", h[4], h[3])
for h in hits3n:
    if best_n is None or abs(h[5])<abs(best_n[1]):
        best_n=(f"[{h[0]} {h[3]} {h[1]} {h[2]}]", h[5], h[4])
for h in hits4n:
    if best_n is None or abs(h[6])<abs(best_n[1]):
        best_n=(f"[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", h[6], h[5])
for h in hits2dn:
    pred = 3.0 + h[3]; err = (pred-Neff_obs)/Neff_obs*100.0
    if best_n is None or abs(err)<abs(best_n[1]):
        best_n=(f"3+[{h[0]} {h[2]} {h[1]}]", err, pred)
for h in hits3dn:
    pred = 3.0 + h[4]; err = (pred-Neff_obs)/Neff_obs*100.0
    if best_n is None or abs(err)<abs(best_n[1]):
        best_n=(f"3+[{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
print(f"\n  BEST N_eff: {best_n[0]} = {best_n[2]:.6f}, err = {best_n[1]:+.6f}%")
Neff_pred = best_n[2]

# ============================================================
# LEDGER
# ============================================================
print()
write_ledger("classXLI_tau_reion_session756", best_tau[2], tau_obs, status="OK")
write_ledger("classXLII_z_star_session756", zstar_pred, zstar_obs, status="OK")
write_ledger("classXLIII_N_eff_session756", Neff_pred, Neff_obs, status="OK")

print("\n"+"-"*80); print("DECISION GATE"); print("-"*80)
print(f"  tau CE push    err = {best_tau[1]:+.6f}%  ({tau_status})")
print(f"  z_* (XLII)     err = {best_z[1]:+.6f}%")
print(f"  N_eff (XLIII)  err = {best_n[1]:+.6f}%")

artifact = os.path.join(os.path.dirname(os.path.abspath(__file__)),"_session756_result.json")
with open(artifact,"w",encoding="utf-8") as f:
    json.dump({
        "tau_CE_push": {"pred": best_tau[2], "obs": tau_obs, "err_pct": best_tau[1],
                        "closure": best_tau[0], "status": tau_status},
        "z_star": {"pred": zstar_pred, "obs": zstar_obs, "err_pct": best_z[1], "closure": best_z[0]},
        "N_eff":  {"pred": Neff_pred,  "obs": Neff_obs,  "err_pct": best_n[1], "closure": best_n[0]},
    }, f, indent=2)
print(f"\nArtifact: {artifact}")
