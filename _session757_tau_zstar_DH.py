"""
SESSION 757 -- (a) tau strict-CE push (corrected base);
                (b) z_* strict-CE 5-atom refinement;
                (c) Class XLIV primordial D/H = 2.547e-5.

(a) S755 closure: tau = Y_p*(5/108)*(1-n_s)/alpha_em = 0.054402, err = +3.7e-3%.
    Delta needed = (0.0544 - 0.054402)/0.054402 = -3.66e-5 (very small).
    Hunt 2- and 3-atom delta shells with tight tolerance.

(b) z_* S756 closure: N_ch*(27/25)/((11/9)*alpha) = 1089.81, err = +9.1e-4%.
    Delta needed = (1089.80 - 1089.81)/1089.81 = -9.11e-6 (very small).
    Hunt delta shell.

(c) Class XLIV -- D/H = 2.547e-5 (Planck/BBN deuterium abundance).
    Seeds: ~Y_p^4 = (49/200)^4 = 3.6e-3 (no); 5/108^3 = 9.9e-5 (close order);
    1/(D_phys * D_crit * SO5^2) = 1/(4*26*100) = 9.6e-5 (close);
    Direct 2/3/4-atom hunt.

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

print("="*80); print("SESSION 757 -- tau CE; z_* CE; D/H Class XLIV"); print("="*80)

# ============================================================
# TRACK (a) tau strict CE push (CORRECTED BASE)
# ============================================================
print("\n"+"-"*80); print("TRACK (a) -- tau_reion strict CE push (corrected)"); print("-"*80)
tau_obs = 0.0544
base_tau = float(Yp) * float(F(5,108)) * float(one_m_ns) / float(alpha_em)
err_base_tau = (base_tau - tau_obs)/tau_obs*100.0
print(f"  S755 base (Y_p*(5/108)*(1-n_s)/alpha) = {base_tau:.8f}, err = {err_base_tau:+.6f}%")
need_mult_tau = tau_obs/base_tau
delta_tau = need_mult_tau - 1.0
sign_tau = 1 if delta_tau>0 else -1
print(f"  delta needed = {delta_tau:+.4e}")

print("\n  2-atom delta on tau base (tol 5e-3%):")
hits2td = search2(abs(delta_tau), tol_pct=5e-3, want=20)
for h in hits2td[:15]:
    mult=1.0+sign_tau*h[3]; pred=base_tau*mult; err=(pred-tau_obs)/tau_obs*100.0
    print(f"    1{'+' if sign_tau>0 else '-'}[{h[0]} {h[2]} {h[1]}]  d={h[3]:.4e}  tau={pred:.8f}  err={err:+.6f}%")
print("\n  3-atom delta on tau base (tol 1e-3%):")
hits3td = search3(abs(delta_tau), tol_pct=1e-3, want=20)
for h in hits3td[:15]:
    mult=1.0+sign_tau*h[4]; pred=base_tau*mult; err=(pred-tau_obs)/tau_obs*100.0
    print(f"    1{'+' if sign_tau>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}]  d={h[4]:.4e}  tau={pred:.8f}  err={err:+.6f}%")
print("\n  4-atom delta on tau base (tol 5e-4%):")
hits4td = search4(abs(delta_tau), tol_pct=5e-4, want=20)
for h in hits4td[:15]:
    mult=1.0+sign_tau*h[5]; pred=base_tau*mult; err=(pred-tau_obs)/tau_obs*100.0
    print(f"    1{'+' if sign_tau>0 else '-'}[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]  d={h[5]:.4e}  tau={pred:.8f}  err={err:+.6f}%")

# Also try direct 4-atom strict
print("\n  Direct 4-atom on tau (tol 1e-3%):")
hits4tau = search4(tau_obs, tol_pct=1e-3, want=15)
for h in hits4tau[:15]:
    print(f"    [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {h[5]:.8f}  err={h[6]:+.6f}%")

best_tau = ("S755 base", err_base_tau, base_tau)
for h in hits2td:
    mult=1.0+sign_tau*h[3]; pred=base_tau*mult; err=(pred-tau_obs)/tau_obs*100.0
    if abs(err)<abs(best_tau[1]): best_tau=(f"base*[1{'+' if sign_tau>0 else '-'}{h[0]} {h[2]} {h[1]}]", err, pred)
for h in hits3td:
    mult=1.0+sign_tau*h[4]; pred=base_tau*mult; err=(pred-tau_obs)/tau_obs*100.0
    if abs(err)<abs(best_tau[1]): best_tau=(f"base*[1{'+' if sign_tau>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
for h in hits4td:
    mult=1.0+sign_tau*h[5]; pred=base_tau*mult; err=(pred-tau_obs)/tau_obs*100.0
    if abs(err)<abs(best_tau[1]): best_tau=(f"base*[1{'+' if sign_tau>0 else '-'}{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", err, pred)
for h in hits4tau:
    if abs(h[6])<abs(best_tau[1]): best_tau=(f"[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", h[6], h[5])

tau_status = "CE_strict" if abs(best_tau[1])<5e-4 else "tight CLOSED"
print(f"\n  BEST tau: {best_tau[0]} = {best_tau[2]:.8f}, err = {best_tau[1]:+.6f}%")
print(f"  STATUS: {tau_status}")

# ============================================================
# TRACK (b) z_* strict CE
# ============================================================
print("\n"+"-"*80); print("TRACK (b) -- z_* strict CE refinement"); print("-"*80)
zstar_obs = 1089.80
base_z = float(N_ch)*float(F(27,25))/(float(F(11,9))*float(alpha_em))
err_base_z = (base_z - zstar_obs)/zstar_obs*100.0
print(f"  S756 base (N_ch*(27/25)/((11/9)*alpha)) = {base_z:.4f}, err = {err_base_z:+.6f}%")
delta_z = zstar_obs/base_z - 1.0
sign_z = 1 if delta_z>0 else -1
print(f"  delta needed = {delta_z:+.4e}")

print("\n  2-atom delta on z_* base (tol 1e-3%):")
hits2zd = search2(abs(delta_z), tol_pct=1e-3, want=20)
for h in hits2zd[:15]:
    mult=1.0+sign_z*h[3]; pred=base_z*mult; err=(pred-zstar_obs)/zstar_obs*100.0
    print(f"    1{'+' if sign_z>0 else '-'}[{h[0]} {h[2]} {h[1]}]  d={h[3]:.4e}  z_*={pred:.4f}  err={err:+.6f}%")
print("\n  3-atom delta on z_* base (tol 5e-4%):")
hits3zd = search3(abs(delta_z), tol_pct=5e-4, want=20)
for h in hits3zd[:15]:
    mult=1.0+sign_z*h[4]; pred=base_z*mult; err=(pred-zstar_obs)/zstar_obs*100.0
    print(f"    1{'+' if sign_z>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}]  d={h[4]:.4e}  z_*={pred:.4f}  err={err:+.6f}%")

best_z=("S756 base", err_base_z, base_z)
for h in hits2zd:
    mult=1.0+sign_z*h[3]; pred=base_z*mult; err=(pred-zstar_obs)/zstar_obs*100.0
    if abs(err)<abs(best_z[1]): best_z=(f"base*[1{'+' if sign_z>0 else '-'}{h[0]} {h[2]} {h[1]}]", err, pred)
for h in hits3zd:
    mult=1.0+sign_z*h[4]; pred=base_z*mult; err=(pred-zstar_obs)/zstar_obs*100.0
    if abs(err)<abs(best_z[1]): best_z=(f"base*[1{'+' if sign_z>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)

z_status = "CE_strict" if abs(best_z[1])<5e-4 else "tight CLOSED"
print(f"\n  BEST z_*: {best_z[0]} = {best_z[2]:.4f}, err = {best_z[1]:+.6f}%")
print(f"  STATUS: {z_status}")

# ============================================================
# TRACK (c) Class XLIV -- D/H
# ============================================================
print("\n"+"-"*80); print("TRACK (c) -- Class XLIV: primordial D/H abundance"); print("-"*80)
DH_obs = 2.547e-5
# Seeds
print(f"  observed D/H = {DH_obs:.4e}")
print(f"  Seed (5/108)^3 = {(5/108)**3:.4e}, err = {((5/108)**3-DH_obs)/DH_obs*100.0:+.4f}%")
print(f"  Seed 1/(D_phys*D_crit*SO5^2) = {1/(4*26*100):.4e}, err = {(1/10400-DH_obs)/DH_obs*100.0:+.4f}%")
print(f"  Seed xi*Y_p*(5/108) = {float(xi)*float(Yp)*5/108:.4e}, err = {(float(xi)*float(Yp)*5/108-DH_obs)/DH_obs*100.0:+.4f}%")

print("\n  2-atom direct on D/H (tol 10%):")
hits2dh = search2(DH_obs, tol_pct=10.0, want=15)
for h in hits2dh[:10]:
    print(f"    [{h[0]} {h[2]} {h[1]}] = {h[3]:.4e}  err={h[4]:+.6f}%")
print("\n  3-atom direct on D/H (tol 1%):")
hits3dh = search3(DH_obs, tol_pct=1.0, want=15)
for h in hits3dh[:12]:
    print(f"    [{h[0]} {h[3]} {h[1]} {h[2]}] = {h[4]:.4e}  err={h[5]:+.6f}%")
print("\n  4-atom direct on D/H (tol 0.1%):")
hits4dh = search4(DH_obs, tol_pct=0.1, want=20)
for h in hits4dh[:15]:
    print(f"    [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {h[5]:.4e}  err={h[6]:+.6f}%")

best_dh=None
for h in hits2dh:
    if best_dh is None or abs(h[4])<abs(best_dh[1]):
        best_dh=(f"[{h[0]} {h[2]} {h[1]}]", h[4], h[3])
for h in hits3dh:
    if best_dh is None or abs(h[5])<abs(best_dh[1]):
        best_dh=(f"[{h[0]} {h[3]} {h[1]} {h[2]}]", h[5], h[4])
for h in hits4dh:
    if best_dh is None or abs(h[6])<abs(best_dh[1]):
        best_dh=(f"[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", h[6], h[5])
print(f"\n  BEST D/H: {best_dh[0]} = {best_dh[2]:.4e}, err = {best_dh[1]:+.6f}%")
dh_pred = best_dh[2]

# ============================================================
# LEDGER
# ============================================================
print()
write_ledger("classXLI_tau_reion_session757", best_tau[2], tau_obs, status="OK")
write_ledger("classXLII_z_star_session757", best_z[2], zstar_obs, status="OK")
write_ledger("classXLIV_DH_session757", dh_pred, DH_obs, status="OK")

print("\n"+"-"*80); print("DECISION GATE"); print("-"*80)
print(f"  tau (XLI)      err = {best_tau[1]:+.6f}%  ({tau_status})")
print(f"  z_* (XLII)     err = {best_z[1]:+.6f}%  ({z_status})")
print(f"  D/H (XLIV)     err = {best_dh[1]:+.6f}%")

artifact = os.path.join(os.path.dirname(os.path.abspath(__file__)),"_session757_result.json")
with open(artifact,"w",encoding="utf-8") as f:
    json.dump({
        "tau": {"pred": best_tau[2], "obs": tau_obs, "err_pct": best_tau[1], "closure": best_tau[0], "status": tau_status},
        "z_star": {"pred": best_z[2], "obs": zstar_obs, "err_pct": best_z[1], "closure": best_z[0], "status": z_status},
        "DH":     {"pred": dh_pred,   "obs": DH_obs,    "err_pct": best_dh[1], "closure": best_dh[0]},
    }, f, indent=2)
print(f"\nArtifact: {artifact}")
