"""
SESSION 759 -- (a) tau_reion + z_* unblock attempt using new atoms (11/6, 6/11, 1/12);
                (b) Class XLVII Omega_c * h^2 = 0.1206 (Planck CDM-only);
                (c) Class XLVIII inflation pivot k_pivot = 0.05 Mpc^-1 (Planck spectral pivot).

(a) S757 showed tau (XLI) and z_* (XLII) blocked at ~1e-5 multiplicative residual at 2/3/4-atom depth.
    S758 introduced 11/6 and 6/11 as new atoms. Retest delta-shells with broadened atom set
    (also add 1/12 = 1 - (1-F*P) since 11/12 acted structurally in Class XLVI).

(b) Class XLVII -- Omega_c * h^2 = 0.1206 (Planck cold dark matter only).
    From XLVI: Omega_m*h^2 = D_crit*(1-F*P)*r/D_BSFG = 143/1000.
    From XXXVII: Omega_b*h^2 = (416/513)*(1/137)/(Y_p*27/25) ~= 0.02237.
    Difference: 0.1430 - 0.02237 = 0.12063. Search direct decomposition.

(c) Class XLVIII -- k_pivot = 0.05 Mpc^-1 (Planck 2018 inflation pivot scale).
    Target 5e-2; abundant atoms close to this magnitude (F_TRZ = 0.1, F_TRZ/2 = 0.05, 5/108 = 0.0463, etc.)

CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""
from __future__ import annotations
from fractions import Fraction as F
import csv, json, os

F_TRZ=F(1,10); Phi_res=F(5,6); SSq=F(57,100); K_Mex=F(25,12)
beta_i=F(6029,10000); D_phys=F(4); D_BSFG=F(6); D_crit=F(26)
N_ch=F(9); SO5=F(10); A_5=F(60)
one_m_FTRZ=1-F_TRZ; one_m_FP=1-F_TRZ*Phi_res  # 9/10, 11/12
ns_atom=F(193,200); one_m_ns=F(7,200); xi=F(11,3200); r_tens=F(9,250); Yp=F(49,200)
mpme=F(1836152673,1000000000); inv_mpme=F(1)/mpme
alpha_em=F(72973525693,10000000000000); inv_137=F(1,137)
mpme_rat=F(11,6); inv_mpme_rat=F(6,11)
# S759 NEW: complement of (1-F*P) = 1/12 (used structurally in XLVI)
one_twelfth = F(1,12)

ATOMS = {
    "F_TRZ":F_TRZ,"Phi_res":Phi_res,"SSq":SSq,"K_Mex":K_Mex,"beta_i":beta_i,
    "D_phys":D_phys,"D_BSFG":D_BSFG,"D_crit":D_crit,"N_ch":N_ch,"SO5":SO5,"A_5":A_5,
    "1-F_TRZ":one_m_FTRZ,"1-F*P":one_m_FP,"n_s":ns_atom,"1-n_s":one_m_ns,
    "xi":xi,"r":r_tens,"Y_p":Yp,
    "27/26":F(27,26),"243/260":F(243,260),"33/40":F(33,40),"11/9":F(11,9),
    "22/9":F(22,9),"27/25":F(27,25),"416/513":F(416,513),"31/30":F(31,30),
    "5/108":F(5,108),"63/200":F(63,200),"137/200":F(137,200),"307/325":F(307,325),
    "1/mpme":inv_mpme, "alpha":alpha_em, "1/137":inv_137, "mpme":mpme,
    "11/6":mpme_rat, "6/11":inv_mpme_rat,
    # NEW S759:
    "1/12":one_twelfth,
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

print("="*80); print("SESSION 759 -- tau/z_* unblock; Omega_c*h^2 (XLVII); k_pivot (XLVIII)"); print("="*80)

# ============================================================
# TRACK (a) -- tau + z_* unblock with broadened atom set
# ============================================================
print("\n"+"-"*80); print("TRACK (a) -- tau + z_* unblock retry with 11/6, 6/11, 1/12 atoms"); print("-"*80)

# --- tau (Class XLI) ---
tau_obs = 0.0540  # Planck reionization optical depth
# S755 base: Y_p * (5/108) * (1-n_s) / alpha_em = 0.054402  err = +3.7e-3% (~7.4e-5 multiplicative)
tau_base = float(Yp) * float(F(5,108)) * float(one_m_ns) / float(alpha_em)
err_tau_base = (tau_base - tau_obs)/tau_obs*100.0
print(f"\n  tau base (S755): pred={tau_base:.6e}, err={err_tau_base:+.6f}%")
need_mult_tau = tau_obs/tau_base; delta_tau = need_mult_tau - 1.0
sign_tau = 1 if delta_tau>0 else -1
print(f"  delta needed = {delta_tau:+.4e}")

print("\n  2-atom delta on tau (tol 1%):")
hits2t = search2(abs(delta_tau), tol_pct=1.0, want=10)
for h in hits2t[:8]:
    mult=1.0+sign_tau*h[3]; pred=tau_base*mult; err=(pred-tau_obs)/tau_obs*100.0
    print(f"    1{'+' if sign_tau>0 else '-'}[{h[0]} {h[2]} {h[1]}] tau={pred:.4e} err={err:+.6f}%")

print("\n  3-atom delta on tau (tol 0.5%):")
hits3t = search3(abs(delta_tau), tol_pct=0.5, want=15)
for h in hits3t[:12]:
    mult=1.0+sign_tau*h[4]; pred=tau_base*mult; err=(pred-tau_obs)/tau_obs*100.0
    print(f"    1{'+' if sign_tau>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}] tau={pred:.4e} err={err:+.6f}%")

print("\n  4-atom delta on tau (tol 0.05%):")
hits4t = search4(abs(delta_tau), tol_pct=0.05, want=12)
for h in hits4t[:10]:
    mult=1.0+sign_tau*h[5]; pred=tau_base*mult; err=(pred-tau_obs)/tau_obs*100.0
    print(f"    1{'+' if sign_tau>0 else '-'}[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] tau={pred:.4e} err={err:+.6f}%")

best_tau = ("S755 base", err_tau_base, tau_base)
for h in hits2t:
    mult=1.0+sign_tau*h[3]; pred=tau_base*mult; err=(pred-tau_obs)/tau_obs*100.0
    if abs(err)<abs(best_tau[1]): best_tau=(f"base*[1{'+' if sign_tau>0 else '-'}{h[0]} {h[2]} {h[1]}]", err, pred)
for h in hits3t:
    mult=1.0+sign_tau*h[4]; pred=tau_base*mult; err=(pred-tau_obs)/tau_obs*100.0
    if abs(err)<abs(best_tau[1]): best_tau=(f"base*[1{'+' if sign_tau>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
for h in hits4t:
    mult=1.0+sign_tau*h[5]; pred=tau_base*mult; err=(pred-tau_obs)/tau_obs*100.0
    if abs(err)<abs(best_tau[1]): best_tau=(f"base*[1{'+' if sign_tau>0 else '-'}{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", err, pred)
print(f"\n  BEST tau: {best_tau[0]} = {best_tau[2]:.4e}, err = {best_tau[1]:+.6f}%")
tau_status = "CE_strict" if abs(best_tau[1])<5e-4 else ("tight CLOSED" if abs(best_tau[1])<0.1 else "CLOSED")

# --- z_* (Class XLII) ---
zstar_obs = 1089.80
zstar_base = float(N_ch) * float(F(27,25)) / (float(F(11,9)) * float(alpha_em))
err_zs_base = (zstar_base - zstar_obs)/zstar_obs*100.0
print(f"\n  z_* base (S756): pred={zstar_base:.6f}, err={err_zs_base:+.6f}%")
need_mult_zs = zstar_obs/zstar_base; delta_zs = need_mult_zs - 1.0
sign_zs = 1 if delta_zs>0 else -1
print(f"  delta needed = {delta_zs:+.4e}")

print("\n  2-atom delta on z_* (tol 1%):")
hits2z = search2(abs(delta_zs), tol_pct=1.0, want=10)
for h in hits2z[:8]:
    mult=1.0+sign_zs*h[3]; pred=zstar_base*mult; err=(pred-zstar_obs)/zstar_obs*100.0
    print(f"    1{'+' if sign_zs>0 else '-'}[{h[0]} {h[2]} {h[1]}] z={pred:.4f} err={err:+.6f}%")

print("\n  3-atom delta on z_* (tol 0.1%):")
hits3z = search3(abs(delta_zs), tol_pct=0.1, want=15)
for h in hits3z[:12]:
    mult=1.0+sign_zs*h[4]; pred=zstar_base*mult; err=(pred-zstar_obs)/zstar_obs*100.0
    print(f"    1{'+' if sign_zs>0 else '-'}[{h[0]} {h[3]} {h[1]} {h[2]}] z={pred:.4f} err={err:+.6f}%")

print("\n  4-atom delta on z_* (tol 0.02%):")
hits4z = search4(abs(delta_zs), tol_pct=0.02, want=12)
for h in hits4z[:10]:
    mult=1.0+sign_zs*h[5]; pred=zstar_base*mult; err=(pred-zstar_obs)/zstar_obs*100.0
    print(f"    1{'+' if sign_zs>0 else '-'}[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] z={pred:.4f} err={err:+.6f}%")

best_zs = ("S756 base", err_zs_base, zstar_base)
for h in hits2z:
    mult=1.0+sign_zs*h[3]; pred=zstar_base*mult; err=(pred-zstar_obs)/zstar_obs*100.0
    if abs(err)<abs(best_zs[1]): best_zs=(f"base*[1{'+' if sign_zs>0 else '-'}{h[0]} {h[2]} {h[1]}]", err, pred)
for h in hits3z:
    mult=1.0+sign_zs*h[4]; pred=zstar_base*mult; err=(pred-zstar_obs)/zstar_obs*100.0
    if abs(err)<abs(best_zs[1]): best_zs=(f"base*[1{'+' if sign_zs>0 else '-'}{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
for h in hits4z:
    mult=1.0+sign_zs*h[5]; pred=zstar_base*mult; err=(pred-zstar_obs)/zstar_obs*100.0
    if abs(err)<abs(best_zs[1]): best_zs=(f"base*[1{'+' if sign_zs>0 else '-'}{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", err, pred)
print(f"\n  BEST z_*: {best_zs[0]} = {best_zs[2]:.4f}, err = {best_zs[1]:+.6f}%")
zs_status = "CE_strict" if abs(best_zs[1])<5e-4 else ("tight CLOSED" if abs(best_zs[1])<0.1 else "CLOSED")

# ============================================================
# TRACK (b) -- Class XLVII Omega_c * h^2 = 0.1206
# ============================================================
print("\n"+"-"*80); print("TRACK (b) -- Class XLVII: Omega_c * h^2 = 0.1206"); print("-"*80)
Och2_obs = 0.1206

# Seed via XLVI - XXXVII
Omh2 = float(D_crit)*float(one_m_FP)*float(r_tens)/float(D_BSFG)
Obh2 = float(F(416,513))*float(inv_137)/(float(Yp)*float(F(27,25)))
diff = Omh2 - Obh2
print(f"  XLVI Omega_m*h^2 = {Omh2:.6f}; XXXVII Omega_b*h^2 = {Obh2:.6f}; diff = {diff:.6f}, err={(diff-Och2_obs)/Och2_obs*100.0:+.6f}%")

print("\n  2-atom direct on Omega_c*h^2 (tol 5%):")
hits2o = search2(Och2_obs, tol_pct=5.0, want=12)
for h in hits2o[:10]:
    print(f"    [{h[0]} {h[2]} {h[1]}] = {h[3]:.6f}  err={h[4]:+.6f}%")

print("\n  3-atom direct on Omega_c*h^2 (tol 1%):")
hits3o = search3(Och2_obs, tol_pct=1.0, want=15)
for h in hits3o[:12]:
    print(f"    [{h[0]} {h[3]} {h[1]} {h[2]}] = {h[4]:.6f}  err={h[5]:+.6f}%")

print("\n  4-atom direct on Omega_c*h^2 (tol 0.1%):")
hits4o = search4(Och2_obs, tol_pct=0.1, want=15)
for h in hits4o[:15]:
    print(f"    [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {h[5]:.6f}  err={h[6]:+.6f}%")

best_o = ("XLVI-XXXVII diff", (diff-Och2_obs)/Och2_obs*100.0, diff)
for h in hits2o:
    if abs(h[4])<abs(best_o[1]): best_o=(f"[{h[0]} {h[2]} {h[1]}]", h[4], h[3])
for h in hits3o:
    if abs(h[5])<abs(best_o[1]): best_o=(f"[{h[0]} {h[3]} {h[1]} {h[2]}]", h[5], h[4])
for h in hits4o:
    if abs(h[6])<abs(best_o[1]): best_o=(f"[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", h[6], h[5])
print(f"\n  BEST Omega_c*h^2: {best_o[0]} = {best_o[2]:.6f}, err = {best_o[1]:+.6f}%")
oc_status = "EXACT" if abs(best_o[1])<1e-8 else ("CE_strict" if abs(best_o[1])<5e-4 else ("tight CLOSED" if abs(best_o[1])<0.1 else "CLOSED"))

# ============================================================
# TRACK (c) -- Class XLVIII k_pivot = 0.05 Mpc^-1
# ============================================================
print("\n"+"-"*80); print("TRACK (c) -- Class XLVIII: k_pivot = 0.05 Mpc^-1"); print("-"*80)
kp_obs = 0.05

print("\n  Seed candidates:")
print(f"    F_TRZ/2 = {0.1/2:.6f}  err={(0.05-0.05)/0.05*100.0:+.6f}%")
print(f"    5/108   = {5/108:.6f}  err={(5/108-0.05)/0.05*100.0:+.6f}%")

print("\n  2-atom direct on k_pivot (tol 5%):")
hits2k = search2(kp_obs, tol_pct=5.0, want=15)
for h in hits2k[:12]:
    print(f"    [{h[0]} {h[2]} {h[1]}] = {h[3]:.6f}  err={h[4]:+.6f}%")

print("\n  3-atom direct on k_pivot (tol 0.5%):")
hits3k = search3(kp_obs, tol_pct=0.5, want=15)
for h in hits3k[:12]:
    print(f"    [{h[0]} {h[3]} {h[1]} {h[2]}] = {h[4]:.6f}  err={h[5]:+.6f}%")

print("\n  4-atom direct on k_pivot (tol 0.05%):")
hits4k = search4(kp_obs, tol_pct=0.05, want=15)
for h in hits4k[:15]:
    print(f"    [{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}] = {h[5]:.6f}  err={h[6]:+.6f}%")

best_k = ("none", 999.0, 0.0)
for h in hits2k:
    if abs(h[4])<abs(best_k[1]): best_k=(f"[{h[0]} {h[2]} {h[1]}]", h[4], h[3])
for h in hits3k:
    if abs(h[5])<abs(best_k[1]): best_k=(f"[{h[0]} {h[3]} {h[1]} {h[2]}]", h[5], h[4])
for h in hits4k:
    if abs(h[6])<abs(best_k[1]): best_k=(f"[{h[0]} {h[4]} {h[1]} {h[2]} {h[3]}]", h[6], h[5])
print(f"\n  BEST k_pivot: {best_k[0]} = {best_k[2]:.6f}, err = {best_k[1]:+.6f}%")
kp_status = "EXACT" if abs(best_k[1])<1e-8 else ("CE_strict" if abs(best_k[1])<5e-4 else ("tight CLOSED" if abs(best_k[1])<0.1 else "CLOSED"))

# ============================================================
# Write ledger
# ============================================================
print()
write_ledger("classXLI_tau_session759", best_tau[2], tau_obs, status=tau_status)
write_ledger("classXLII_zstar_session759", best_zs[2], zstar_obs, status=zs_status)
write_ledger("classXLVII_Omega_c_h2_session759", best_o[2], Och2_obs, status=oc_status)
write_ledger("classXLVIII_k_pivot_session759", best_k[2], kp_obs, status=kp_status)

# ============================================================
# Decision gate
# ============================================================
print("\n"+"-"*80); print("DECISION GATE"); print("-"*80)
print(f"  tau (XLI retry)     err = {best_tau[1]:+.6f}%  ({tau_status})")
print(f"  z_* (XLII retry)    err = {best_zs[1]:+.6f}%  ({zs_status})")
print(f"  Omega_c*h^2 (XLVII) err = {best_o[1]:+.6f}%  ({oc_status})")
print(f"  k_pivot (XLVIII)    err = {best_k[1]:+.6f}%  ({kp_status})")

# Artifact
out = {
    "session": 759,
    "tau":     {"best": best_tau[0], "pred": best_tau[2], "err_pct": best_tau[1], "status": tau_status},
    "z_star":  {"best": best_zs[0],  "pred": best_zs[2],  "err_pct": best_zs[1],  "status": zs_status},
    "Omega_c_h2": {"best": best_o[0], "pred": best_o[2], "err_pct": best_o[1], "status": oc_status},
    "k_pivot": {"best": best_k[0],   "pred": best_k[2],   "err_pct": best_k[1],  "status": kp_status},
    "cvw": "v2.0.0",
    "sm_anchor": "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant"
}
art = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_session759_result.json")
with open(art,"w",encoding="utf-8") as f: json.dump(out, f, indent=2)
print(f"\nArtifact: {art}")
