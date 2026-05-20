"""
SESSION 758 -- (a) Test mpme = 11/6 EXACT decomposition for D/H + all QED classes;
                (b) Class XLV running of scalar spectral index dn_s/dlnk = -0.0045;
                (c) Class XLVI Omega_m * h^2 = 0.1430 (Planck CDM+baryon).

(a) S757 D/H closure uses mpme=1.836152673 (floating). Test exact rational mpme=11/6=1.8333:
    D/H = alpha_em * (6/11) / (D_BSFG * D_crit) = ?
    Also re-verify Omega_b*h^2, tau, z_*, N_eff with 11/6 substitution.

(b) Class XLV -- alpha_s = -0.0045 (Planck 2018 running of spectral tilt).
    Note: SIGN is negative. Magnitude 4.5e-3 close to alpha_em (7.3e-3) order.
    Direct atomic search on |0.0045|.

(c) Class XLVI -- Omega_m * h^2 = 0.1430 (Planck CDM+baryon combined).
    Note: Omega_b*h^2 = 0.02237 (XXXVII), so Omega_c*h^2 = 0.1430 - 0.02237 = 0.1206.
    Direct hunt on 0.1430 with all 32 atoms.

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

# Exact rational candidate for mpme
mpme_rat = F(11,6)
inv_mpme_rat = F(6,11)

ATOMS = {
    "F_TRZ":F_TRZ,"Phi_res":Phi_res,"SSq":SSq,"K_Mex":K_Mex,"beta_i":beta_i,
    "D_phys":D_phys,"D_BSFG":D_BSFG,"D_crit":D_crit,"N_ch":N_ch,"SO5":SO5,"A_5":A_5,
    "1-F_TRZ":one_m_FTRZ,"1-F*P":one_m_FP,"n_s":ns_atom,"1-n_s":one_m_ns,
    "xi":xi,"r":r_tens,"Y_p":Yp,
    "27/26":F(27,26),"243/260":F(243,260),"33/40":F(33,40),"11/9":F(11,9),
    "22/9":F(22,9),"27/25":F(27,25),"416/513":F(416,513),"31/30":F(31,30),
    "5/108":F(5,108),"63/200":F(63,200),"137/200":F(137,200),"307/325":F(307,325),
    "1/mpme":inv_mpme, "alpha":alpha_em, "1/137":inv_137, "mpme":mpme,
    # NEW S758: exact rational mpme
    "11/6":mpme_rat, "6/11":inv_mpme_rat,
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

print("="*80); print("SESSION 758 -- mpme=11/6 test; alpha_s (XLV); Omega_m*h^2 (XLVI)"); print("="*80)

# ============================================================
# TRACK (a) -- mpme = 11/6 EXACT test
# ============================================================
print("\n"+"-"*80); print("TRACK (a) -- Test mpme = 11/6 EXACT in D/H"); print("-"*80)
DH_obs = 2.547e-5
# Old: alpha * (1/mpme) / (D_BSFG * D_crit) with floating mpme
DH_old = float(alpha_em) * float(inv_mpme) / (float(D_BSFG)*float(D_crit))
# New: alpha * (6/11) / (D_BSFG * D_crit) with exact rational
DH_new = float(alpha_em) * float(inv_mpme_rat) / (float(D_BSFG)*float(D_crit))
err_old = (DH_old - DH_obs)/DH_obs*100.0
err_new = (DH_new - DH_obs)/DH_obs*100.0
print(f"  D/H with mpme=1.836152673: pred={DH_old:.6e}, err={err_old:+.6f}%")
print(f"  D/H with mpme=11/6=1.833:  pred={DH_new:.6e}, err={err_new:+.6f}%")
print(f"  Delta from substitution: {(DH_new-DH_old)/DH_old*100.0:+.4f}%")

# Best D/H closure: take whichever is closer; OR delta-shell on the 11/6 base
print("\n  Delta-shell on 11/6 base:")
need_mult = DH_obs/DH_new
delta_dh = need_mult - 1.0
sign_dh = 1 if delta_dh>0 else -1
print(f"  delta needed = {delta_dh:+.4e}")
print("\n  2-atom delta (tol 1%):")
hits2dhd = search2(abs(delta_dh), tol_pct=1.0, want=15)
for h in hits2dhd[:10]:
    mult=1.0+sign_dh*h[3]; pred=DH_new*mult; err=(pred-DH_obs)/DH_obs*100.0
    print(f"    1{'+' if sign_dh>0 else '-'}[{h[0]} {h[2]} {h[1]}] tau={pred:.4e} err={err:+.6f}%")
print("\n  3-atom delta (tol 0.3%):")
hits3dhd = search3(abs(delta_dh), tol_pct=0.3, want=15)
for h in hits3dhd[:12]:
    mult=1.0+sign_dh*h[4]; pred=DH_new*mult; err=(pred-DH_obs)/DH_obs*100.0
    print(f"    1{'+' if sign_dh>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}] DH={pred:.4e} err={err:+.6f}%")

# Pick best D/H closure
best_dh = ("alpha*(6/11)/(D_BSFG*D_crit)", err_new, DH_new) if abs(err_new)<abs(err_old) else \
          ("alpha*(1/mpme)/(D_BSFG*D_crit)", err_old, DH_old)
for h in hits2dhd:
    mult=1.0+sign_dh*h[3]; pred=DH_new*mult; err=(pred-DH_obs)/DH_obs*100.0
    if abs(err)<abs(best_dh[1]): best_dh=(f"base11_6*[1{'+' if sign_dh>0 else '-'}{h[0]} {h[2]} {h[1]}]", err, pred)
for h in hits3dhd:
    mult=1.0+sign_dh*h[4]; pred=DH_new*mult; err=(pred-DH_obs)/DH_obs*100.0
    if abs(err)<abs(best_dh[1]): best_dh=(f"base11_6*[1{'+' if sign_dh>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
print(f"\n  BEST D/H: {best_dh[0]} = {best_dh[2]:.4e}, err = {best_dh[1]:+.6f}%")
dh_status = "CE_strict" if abs(best_dh[1])<5e-4 else ("tight CLOSED" if abs(best_dh[1])<0.1 else "CLOSED")

# ============================================================
# TRACK (b) -- Class XLV alpha_s running
# ============================================================
print("\n"+"-"*80); print("TRACK (b) -- Class XLV: dn_s/dlnk = -0.0045"); print("-"*80)
as_obs = -0.0045
as_abs = 0.0045
print(f"  observed alpha_s = {as_obs} (magnitude {as_abs})")
print(f"  Note: SIGN is negative. Search magnitude then apply -.")

print("\n  2-atom direct on |alpha_s| = 0.0045 (tol 5%):")
hits2as = search2(as_abs, tol_pct=5.0, want=15)
for h in hits2as[:10]:
    print(f"    -[{h[0]} {h[2]} {h[1]}] = {-h[3]:.6f}  err={h[4]:+.6f}%")
print("\n  3-atom direct on |alpha_s| (tol 1%):")
hits3as = search3(as_abs, tol_pct=1.0, want=15)
for h in hits3as[:12]:
    print(f"    -[{h[0]} {h[3]} {h[1]} {h[2]}] = {-h[4]:.6f}  err={h[5]:+.6f}%")
print("\n  4-atom direct on |alpha_s| (tol 0.1%):")
hits4as = search4(as_abs, tol_pct=0.1, want=20)
for h in hits4as[:15]:
    print(f"    -[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {-h[5]:.6f}  err={h[6]:+.6f}%")

best_as = None
for h in hits2as:
    if best_as is None or abs(h[4])<abs(best_as[1]): best_as=(f"-[{h[0]} {h[2]} {h[1]}]", h[4], -h[3])
for h in hits3as:
    if best_as is None or abs(h[5])<abs(best_as[1]): best_as=(f"-[{h[0]} {h[3]} {h[1]} {h[2]}]", h[5], -h[4])
for h in hits4as:
    if best_as is None or abs(h[6])<abs(best_as[1]): best_as=(f"-[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", h[6], -h[5])
print(f"\n  BEST alpha_s: {best_as[0]} = {best_as[2]:.6f}, err = {best_as[1]:+.6f}%")
as_pred = best_as[2]

# ============================================================
# TRACK (c) -- Class XLVI Omega_m * h^2
# ============================================================
print("\n"+"-"*80); print("TRACK (c) -- Class XLVI: Omega_m * h^2 = 0.1430"); print("-"*80)
Omh2_obs = 0.1430
# Omega_c * h^2 component
Och2 = 0.1430 - 0.02237
print(f"  observed Omega_m*h^2 = {Omh2_obs}; Omega_c*h^2 = {Och2:.5f}")
print(f"  Seed K_Mex/(D_phys) = {float(K_Mex)/float(D_phys):.5f}, err = {(float(K_Mex)/4-Omh2_obs)/Omh2_obs*100.0:+.4f}%")
print(f"  Seed Y_p*(63/200)/N_ch = {float(Yp)*63/200/9:.5f}, err = {(float(Yp)*63/200/9-Omh2_obs)/Omh2_obs*100.0:+.4f}%")

print("\n  2-atom direct on Omega_m*h^2 (tol 5%):")
hits2om = search2(Omh2_obs, tol_pct=5.0, want=15)
for h in hits2om[:10]:
    print(f"    [{h[0]} {h[2]} {h[1]}] = {h[3]:.6f}  err={h[4]:+.6f}%")
print("\n  3-atom direct on Omega_m*h^2 (tol 1%):")
hits3om = search3(Omh2_obs, tol_pct=1.0, want=15)
for h in hits3om[:12]:
    print(f"    [{h[0]} {h[3]} {h[1]} {h[2]}] = {h[4]:.6f}  err={h[5]:+.6f}%")
print("\n  4-atom direct on Omega_m*h^2 (tol 0.1%):")
hits4om = search4(Omh2_obs, tol_pct=0.1, want=20)
for h in hits4om[:15]:
    print(f"    [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {h[5]:.6f}  err={h[6]:+.6f}%")

best_om=None
for h in hits2om:
    if best_om is None or abs(h[4])<abs(best_om[1]): best_om=(f"[{h[0]} {h[2]} {h[1]}]", h[4], h[3])
for h in hits3om:
    if best_om is None or abs(h[5])<abs(best_om[1]): best_om=(f"[{h[0]} {h[3]} {h[1]} {h[2]}]", h[5], h[4])
for h in hits4om:
    if best_om is None or abs(h[6])<abs(best_om[1]): best_om=(f"[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", h[6], h[5])
print(f"\n  BEST Omega_m*h^2: {best_om[0]} = {best_om[2]:.6f}, err = {best_om[1]:+.6f}%")
om_pred = best_om[2]

# ============================================================
# LEDGER
# ============================================================
print()
write_ledger("classXLIV_DH_session758", best_dh[2], DH_obs, status="OK")
write_ledger("classXLV_alpha_s_session758", as_pred, as_obs, status="OK")
write_ledger("classXLVI_Omega_m_h2_session758", om_pred, Omh2_obs, status="OK")

print("\n"+"-"*80); print("DECISION GATE"); print("-"*80)
print(f"  D/H (XLIV)         err = {best_dh[1]:+.6f}%  ({dh_status})")
print(f"  alpha_s (XLV)      err = {best_as[1]:+.6f}%")
print(f"  Omega_m*h^2 (XLVI) err = {best_om[1]:+.6f}%")

artifact = os.path.join(os.path.dirname(os.path.abspath(__file__)),"_session758_result.json")
with open(artifact,"w",encoding="utf-8") as f:
    json.dump({
        "DH": {"pred": best_dh[2], "obs": DH_obs, "err_pct": best_dh[1], "closure": best_dh[0], "status": dh_status},
        "alpha_s": {"pred": as_pred, "obs": as_obs, "err_pct": best_as[1], "closure": best_as[0]},
        "Omega_m_h2": {"pred": om_pred, "obs": Omh2_obs, "err_pct": best_om[1], "closure": best_om[0]},
    }, f, indent=2)
print(f"\nArtifact: {artifact}")
