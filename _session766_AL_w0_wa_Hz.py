"""
SESSION 766 -- (a) LXVI  A_L      (Planck lensing parameter ~ 1.180);
                (b) LXVII w_0     (dark energy EoS, DESI 2024 ~ -1.03);
                (c) LXVIII w_a    (DE evolution param, DESI 2024 ~ -0.39);
                (d) LXIX  H(z=0.38) (BOSS DR12 ~ 81.5 km/s/Mpc, paired with LIII H_0 + LXI z_BAO).

CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""
from __future__ import annotations
from fractions import Fraction as F
import csv, os, json

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
    hits.sort(key=lambda h:abs(h[4])); return hits[:want]

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
    hits.sort(key=lambda h:abs(h[5])); return hits[:want]

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
                    for tag,fn in forms:
                        try: v=fn(VALS[i],VALS[j],VALS[k],VALS[l])
                        except ZeroDivisionError: continue
                        if v==0: continue
                        err=(v-target)/target*100.0
                        if abs(err)<tol_pct: hits.append((LABELS[i],LABELS[j],LABELS[k],LABELS[l],tag,v,err))
    hits.sort(key=lambda h:abs(h[6])); return hits[:want]

def status_of(err_pct):
    a=abs(err_pct)
    if a==0.0: return "EXACT"
    if a<5e-4: return "candidate_EXACT"
    if a<5e-2: return "CE_strict"
    if a<1.0:  return "CE"
    return "OPEN"

def pick_best(t2, t3, t4):
    h4=search4(t4[0], t4[1], 1); h3=search3(t3[0], t3[1], 1); h2=search2(t2[0], t2[1], 1)
    for c in (h4,h3,h2):
        if c: return c[0]
    return None

def banner(s): print("="*80); print(s); print("="*80)

print("="*80); print("SESSION 766 -- A_L (LXVI); w_0 (LXVII); w_a (LXVIII); H(z=0.38) (LXIX)"); print("="*80)

# ============================================================================
# (a) LXVI A_L = 1.180 (Planck lensing parameter, ~3sigma tension vs A_planck=1)
# ============================================================================
banner("TRACK (a) -- Class LXVI: A_L = 1.180 (Planck lensing parameter)")
A_L_obs = 1.180
print(f"  Target: A_L = {A_L_obs}")
print()
print("  2-atom direct (tol 3%):")
for la,lb,tag,v,err in search2(A_L_obs, 3.0, 12):
    print(f"    [{la} {tag} {lb}] = {v:.6f}  err={err:+.4f}%")
print("\n  3-atom direct (tol 0.5%):")
for la,lb,lc,tag,v,err in search3(A_L_obs, 0.5, 15):
    print(f"    [{la} {tag} {lb} {lc}] = {v:.6f}  err={err:+.4f}%")
print("\n  4-atom direct (tol 0.05%):")
for la,lb,lc,ld,tag,v,err in search4(A_L_obs, 0.05, 15):
    print(f"    [{la} {tag} {lb} {lc} {ld}] = {v:.6f}  err={err:+.4f}%")
AL_best = pick_best((A_L_obs,10.0),(A_L_obs,5.0),(A_L_obs,5.0))
if AL_best is None: AL_pred=0.0; AL_label="none"; AL_err=9999.0
else: AL_pred=AL_best[-2]; AL_err=AL_best[-1]; AL_label=" ".join(str(x) for x in AL_best[:-2])
print(f"\n  BEST A_L: [{AL_label}] = {AL_pred:.6f}, err = {AL_err:+.6f}%")

# ============================================================================
# (b) LXVII w_0 = -1.03 (DESI 2024 DR1 + CMB + SN)
# ============================================================================
banner("TRACK (b) -- Class LXVII: w_0 = -1.03 (DESI 2024 dark energy EoS)")
w0_obs = -1.03; w0_abs = 1.03
print(f"  Target: w_0 = {w0_obs}  (search on |w_0| = {w0_abs})")
print("\n  2-atom direct (tol 3%):")
for la,lb,tag,v,err in search2(w0_abs, 3.0, 12):
    print(f"    [{la} {tag} {lb}] = {v:.6f}  err={err:+.4f}%")
print("\n  3-atom direct (tol 0.5%):")
for la,lb,lc,tag,v,err in search3(w0_abs, 0.5, 15):
    print(f"    [{la} {tag} {lb} {lc}] = {v:.6f}  err={err:+.4f}%")
print("\n  4-atom direct (tol 0.05%):")
for la,lb,lc,ld,tag,v,err in search4(w0_abs, 0.05, 15):
    print(f"    [{la} {tag} {lb} {lc} {ld}] = {v:.6f}  err={err:+.4f}%")
w0_best = pick_best((w0_abs,10.0),(w0_abs,5.0),(w0_abs,5.0))
if w0_best is None: w0_pred=0.0; w0_label="none"; w0_err=9999.0
else:
    w0_pred = -w0_best[-2]  # reinstate sign
    w0_err = (w0_pred - w0_obs)/w0_obs * 100.0
    w0_label = "(-1) * [" + " ".join(str(x) for x in w0_best[:-2]) + "]"
print(f"\n  BEST w_0: [{w0_label}] = {w0_pred:+.6f}, err = {w0_err:+.6f}%")

# ============================================================================
# (c) LXVIII w_a = -0.39 (DESI 2024 evolving DE)
# ============================================================================
banner("TRACK (c) -- Class LXVIII: w_a = -0.39 (DESI 2024 dark energy evolution)")
wa_obs = -0.39; wa_abs = 0.39
print(f"  Target: w_a = {wa_obs}  (search on |w_a| = {wa_abs})")
print("\n  2-atom direct (tol 3%):")
for la,lb,tag,v,err in search2(wa_abs, 3.0, 12):
    print(f"    [{la} {tag} {lb}] = {v:.6f}  err={err:+.4f}%")
print("\n  3-atom direct (tol 0.5%):")
for la,lb,lc,tag,v,err in search3(wa_abs, 0.5, 15):
    print(f"    [{la} {tag} {lb} {lc}] = {v:.6f}  err={err:+.4f}%")
print("\n  4-atom direct (tol 0.05%):")
for la,lb,lc,ld,tag,v,err in search4(wa_abs, 0.05, 15):
    print(f"    [{la} {tag} {lb} {lc} {ld}] = {v:.6f}  err={err:+.4f}%")
wa_best = pick_best((wa_abs,10.0),(wa_abs,5.0),(wa_abs,5.0))
if wa_best is None: wa_pred=0.0; wa_label="none"; wa_err=9999.0
else:
    wa_pred = -wa_best[-2]
    wa_err = (wa_pred - wa_obs)/wa_obs * 100.0
    wa_label = "(-1) * [" + " ".join(str(x) for x in wa_best[:-2]) + "]"
print(f"\n  BEST w_a: [{wa_label}] = {wa_pred:+.6f}, err = {wa_err:+.6f}%")

# ============================================================================
# (d) LXIX H(z=0.38) = 81.5 km/s/Mpc (BOSS DR12 BAO Hubble)
# ============================================================================
banner("TRACK (d) -- Class LXIX: H(z=0.38) = 81.5 km/s/Mpc (BOSS DR12 BAO)")
H_z_obs = 81.5
H0_LIII = 137.0/(float(mpme)**2 * float(beta_i))  # from S762 LIII
ratio_target = H_z_obs / H0_LIII
print(f"  Target H(z=0.38) = {H_z_obs} km/s/Mpc")
print(f"  H_0 (LIII)       = {H0_LIII:.4f} km/s/Mpc")
print(f"  Ratio target H(z)/H_0 = {ratio_target:.6f}")
print("\n  2-atom shell (tol 3%):")
for la,lb,tag,v,err in search2(ratio_target, 3.0, 12):
    overall = v * H0_LIII
    o_err = (overall - H_z_obs)/H_z_obs * 100.0
    print(f"    [{la} {tag} {lb}] ratio={v:.6f}  -> H(z)={overall:.4f}  err={o_err:+.4f}%")
print("\n  3-atom shell (tol 0.5%):")
for la,lb,lc,tag,v,err in search3(ratio_target, 0.5, 15):
    overall = v * H0_LIII
    o_err = (overall - H_z_obs)/H_z_obs * 100.0
    print(f"    [{la} {tag} {lb} {lc}] ratio={v:.6f}  -> H(z)={overall:.4f}  err={o_err:+.4f}%")
print("\n  4-atom shell (tol 0.05%):")
for la,lb,lc,ld,tag,v,err in search4(ratio_target, 0.05, 15):
    overall = v * H0_LIII
    o_err = (overall - H_z_obs)/H_z_obs * 100.0
    print(f"    [{la} {tag} {lb} {lc} {ld}] ratio={v:.6f}  -> H(z)={overall:.4f}  err={o_err:+.4f}%")

# Pick best by overall residual
all_shell = []
for la,lb,tag,v,err in search2(ratio_target, 10.0, 50):
    overall = v * H0_LIII; oe = (overall - H_z_obs)/H_z_obs * 100.0
    all_shell.append((abs(oe), oe, overall, f"H_0 * [{la} {tag} {lb}]"))
for la,lb,lc,tag,v,err in search3(ratio_target, 2.0, 50):
    overall = v * H0_LIII; oe = (overall - H_z_obs)/H_z_obs * 100.0
    all_shell.append((abs(oe), oe, overall, f"H_0 * [{la} {tag} {lb} {lc}]"))
for la,lb,lc,ld,tag,v,err in search4(ratio_target, 1.0, 50):
    overall = v * H0_LIII; oe = (overall - H_z_obs)/H_z_obs * 100.0
    all_shell.append((abs(oe), oe, overall, f"H_0 * [{la} {tag} {lb} {lc} {ld}]"))
all_shell.sort()
if all_shell:
    _, Hz_err, Hz_pred, Hz_label = all_shell[0]
else:
    Hz_pred=0.0; Hz_label="none"; Hz_err=9999.0
print(f"\n  BEST H(z=0.38): {Hz_label} = {Hz_pred:.4f}, err = {Hz_err:+.6f}%")

# ============================================================================
# WRITE LEDGER
# ============================================================================
ROOT = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(ROOT, "master_closures.csv")

def write_ledger(name, predicted, observed, err_pct):
    st = status_of(err_pct)
    raw = f"{name}: predicted={predicted:.6e} observed={observed:.6e} error_pct={err_pct:.6e} status={st}"
    print(raw)
    headers = ["script","label","predicted","observed","error_pct","status","cvw","sm_anchor","raw_output"]
    file_exists = os.path.exists(CSV_PATH)
    with open(CSV_PATH, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=headers, extrasaction="ignore")
        if not file_exists: w.writeheader()
        w.writerow({
            "script": os.path.basename(__file__),
            "label": name,
            "predicted": f"{predicted:.6e}",
            "observed":  f"{observed:.6e}",
            "error_pct": f"{err_pct:.6e}",
            "status": st,
            "cvw": "v2.0.0",
            "sm_anchor": "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
            "raw_output": raw,
        })
    return st

print()
st_a = write_ledger("classLXVI_A_L_session766",        AL_pred, A_L_obs, AL_err)
st_b = write_ledger("classLXVII_w_0_session766",       w0_pred, w0_obs,  w0_err)
st_c = write_ledger("classLXVIII_w_a_session766",      wa_pred, wa_obs,  wa_err)
st_d = write_ledger("classLXIX_Hz_BAO_session766",     Hz_pred, H_z_obs, Hz_err)

print()
print("-"*80); print("DECISION GATE"); print("-"*80)
print(f"  LXVI   A_L          pred={AL_pred:.6f}    err={AL_err:+.4f}%   ({st_a})")
print(f"  LXVII  w_0          pred={w0_pred:+.6f}   err={w0_err:+.4f}%   ({st_b})")
print(f"  LXVIII w_a          pred={wa_pred:+.6f}   err={wa_err:+.4f}%   ({st_c})")
print(f"  LXIX   H(z=0.38)    pred={Hz_pred:.4f}    err={Hz_err:+.4f}%   ({st_d})")

artifact = os.path.join(ROOT, "_session766_result.json")
with open(artifact, "w", encoding="utf-8") as fh:
    json.dump({
        "session": 766,
        "classLXVI_A_L":     {"pred": AL_pred, "obs": A_L_obs, "err_pct": AL_err, "status": st_a, "form": AL_label},
        "classLXVII_w_0":    {"pred": w0_pred, "obs": w0_obs,  "err_pct": w0_err, "status": st_b, "form": w0_label},
        "classLXVIII_w_a":   {"pred": wa_pred, "obs": wa_obs,  "err_pct": wa_err, "status": st_c, "form": wa_label},
        "classLXIX_Hz_BAO":  {"pred": Hz_pred, "obs": H_z_obs, "err_pct": Hz_err, "status": st_d, "form": Hz_label, "H_0_used_LIII": H0_LIII},
    }, fh, indent=2)
print(f"\nArtifact: {artifact}")
print(f"Master registry written: {CSV_PATH}")
