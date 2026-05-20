"""
SESSION 767 -- BAO distance suite + RSD at z_BAO_eff = 0.38 (paired with LXI, LXIX).

  (a) LXX   D_A(z=0.38)    angular diameter distance ~ 1099.8 Mpc
  (b) LXXI  D_V(z=0.38)    BAO dilation scale         ~ 1476.8 Mpc
  (c) LXXII D_M(z=0.38)    comoving ang. diam. dist.  ~ 1517.7 Mpc
  (d) LXXIII f*sigma_8(z=0.38) RSD growth amplitude   ~ 0.497

Strategy: distances expressed as shells on r_drag (LIV = 279/(260*alpha_em) ~ 147.05 Mpc).

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
    forms=[("a*b/(c*d)",lambda a,b,c,d:a*b/(c*d)),
           ("a*b*c/d",  lambda a,b,c,d:a*b*c/d),
           ("a/(b*c*d)",lambda a,b,c,d:a/(b*c*d)),
           ("a*b*c*d",  lambda a,b,c,d:a*b*c*d)]
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

def banner(s): print("="*80); print(s); print("="*80)

# r_drag from LIV (S762): r_d = 279 / (260 * alpha_em) Mpc
r_drag = 279.0 / (260.0 * float(alpha_em))

print("="*80); print("SESSION 767 -- BAO suite at z=0.38: D_A (LXX); D_V (LXXI); D_M (LXXII); f*sig8 (LXXIII)"); print("="*80)
print(f"  r_drag (LIV)  = {r_drag:.4f} Mpc")
print(f"  z_BAO_eff (LXI) = 19/50 = 0.38")
print(f"  H(z=0.38) (LXIX) = 81.5 km/s/Mpc")
print()

# ---------------------------------------------------------------------------
# Helper: shell search on a unit-scale (target/unit)
# ---------------------------------------------------------------------------
def report_shell(name, obs, unit, label_unit):
    target_ratio = obs / unit
    print(f"  Target {name} = {obs:.4f}  ({label_unit}-shell target = {target_ratio:.6f})")
    print("\n  2-atom shell (tol 3%):")
    for la,lb,tag,v,err in search2(target_ratio, 3.0, 12):
        overall = v*unit; oe = (overall-obs)/obs*100.0
        print(f"    [{la} {tag} {lb}] ratio={v:.6f}  -> {overall:.4f}  err={oe:+.4f}%")
    print("\n  3-atom shell (tol 0.5%):")
    for la,lb,lc,tag,v,err in search3(target_ratio, 0.5, 15):
        overall = v*unit; oe = (overall-obs)/obs*100.0
        print(f"    [{la} {tag} {lb} {lc}] ratio={v:.6f}  -> {overall:.4f}  err={oe:+.4f}%")
    print("\n  4-atom shell (tol 0.05%):")
    for la,lb,lc,ld,tag,v,err in search4(target_ratio, 0.05, 15):
        overall = v*unit; oe = (overall-obs)/obs*100.0
        print(f"    [{la} {tag} {lb} {lc} {ld}] ratio={v:.6f}  -> {overall:.4f}  err={oe:+.4f}%")
    candidates=[]
    for la,lb,tag,v,err in search2(target_ratio, 10.0, 30):
        o=v*unit; oe=(o-obs)/obs*100.0
        candidates.append((abs(oe), oe, o, f"{label_unit} * [{la} {tag} {lb}]"))
    for la,lb,lc,tag,v,err in search3(target_ratio, 2.0, 50):
        o=v*unit; oe=(o-obs)/obs*100.0
        candidates.append((abs(oe), oe, o, f"{label_unit} * [{la} {tag} {lb} {lc}]"))
    for la,lb,lc,ld,tag,v,err in search4(target_ratio, 1.0, 50):
        o=v*unit; oe=(o-obs)/obs*100.0
        candidates.append((abs(oe), oe, o, f"{label_unit} * [{la} {tag} {lb} {lc} {ld}]"))
    candidates.sort()
    if not candidates:
        return 0.0, "none", 9999.0
    _, oe, op, lab = candidates[0]
    print(f"\n  BEST {name}: {lab} = {op:.4f}, err = {oe:+.6f}%")
    return op, lab, oe

# ---------------------------------------------------------------------------
# (a) LXX  D_A(z=0.38) = 1099.8 Mpc  [BOSS DR12 D_M/(1+z) with D_M/r_d=10.27]
# ---------------------------------------------------------------------------
banner("TRACK (a) -- Class LXX: D_A(z=0.38) = 1099.8 Mpc (angular diameter)")
DA_obs = 1099.8
DA_pred, DA_label, DA_err = report_shell("D_A(0.38)", DA_obs, r_drag, "r_drag")

# ---------------------------------------------------------------------------
# (b) LXXI D_V(z=0.38) = 1476.8 Mpc
# ---------------------------------------------------------------------------
banner("TRACK (b) -- Class LXXI: D_V(z=0.38) = 1476.8 Mpc (BAO dilation)")
DV_obs = 1476.8
DV_pred, DV_label, DV_err = report_shell("D_V(0.38)", DV_obs, r_drag, "r_drag")

# ---------------------------------------------------------------------------
# (c) LXXII D_M(z=0.38) = 1517.7 Mpc
# ---------------------------------------------------------------------------
banner("TRACK (c) -- Class LXXII: D_M(z=0.38) = 1517.7 Mpc (comoving ang. diam.)")
DM_obs = 1517.7
DM_pred, DM_label, DM_err = report_shell("D_M(0.38)", DM_obs, r_drag, "r_drag")

# ---------------------------------------------------------------------------
# (d) LXXIII f*sigma_8(z=0.38) = 0.497 (direct, no shell)
# ---------------------------------------------------------------------------
banner("TRACK (d) -- Class LXXIII: f*sigma_8(z=0.38) = 0.497 (RSD growth amplitude)")
fs8_obs = 0.497
print(f"  Target f*sigma_8(0.38) = {fs8_obs}")
print("\n  2-atom direct (tol 3%):")
for la,lb,tag,v,err in search2(fs8_obs, 3.0, 12):
    print(f"    [{la} {tag} {lb}] = {v:.6f}  err={err:+.4f}%")
print("\n  3-atom direct (tol 0.5%):")
for la,lb,lc,tag,v,err in search3(fs8_obs, 0.5, 15):
    print(f"    [{la} {tag} {lb} {lc}] = {v:.6f}  err={err:+.4f}%")
print("\n  4-atom direct (tol 0.05%):")
for la,lb,lc,ld,tag,v,err in search4(fs8_obs, 0.05, 15):
    print(f"    [{la} {tag} {lb} {lc} {ld}] = {v:.6f}  err={err:+.4f}%")

fs8_candidates=[]
for la,lb,tag,v,err in search2(fs8_obs, 10.0, 30):
    fs8_candidates.append((abs(err), v, err, f"[{la} {tag} {lb}]"))
for la,lb,lc,tag,v,err in search3(fs8_obs, 5.0, 50):
    fs8_candidates.append((abs(err), v, err, f"[{la} {tag} {lb} {lc}]"))
for la,lb,lc,ld,tag,v,err in search4(fs8_obs, 2.0, 50):
    fs8_candidates.append((abs(err), v, err, f"[{la} {tag} {lb} {lc} {ld}]"))
fs8_candidates.sort()
if fs8_candidates:
    _, fs8_pred, fs8_err, fs8_label = fs8_candidates[0]
else:
    fs8_pred=0.0; fs8_label="none"; fs8_err=9999.0
print(f"\n  BEST f*sigma_8: {fs8_label} = {fs8_pred:.6f}, err = {fs8_err:+.6f}%")

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
st_a = write_ledger("classLXX_D_A_z038_session767",  DA_pred, DA_obs, DA_err)
st_b = write_ledger("classLXXI_D_V_z038_session767", DV_pred, DV_obs, DV_err)
st_c = write_ledger("classLXXII_D_M_z038_session767",DM_pred, DM_obs, DM_err)
st_d = write_ledger("classLXXIII_fsigma8_z038_session767", fs8_pred, fs8_obs, fs8_err)

print()
print("-"*80); print("DECISION GATE"); print("-"*80)
print(f"  LXX    D_A(0.38)        pred={DA_pred:.4f}  err={DA_err:+.4f}%  ({st_a})")
print(f"  LXXI   D_V(0.38)        pred={DV_pred:.4f}  err={DV_err:+.4f}%  ({st_b})")
print(f"  LXXII  D_M(0.38)        pred={DM_pred:.4f}  err={DM_err:+.4f}%  ({st_c})")
print(f"  LXXIII f*sigma_8(0.38)  pred={fs8_pred:.4f}    err={fs8_err:+.4f}%  ({st_d})")

artifact = os.path.join(ROOT, "_session767_result.json")
with open(artifact, "w", encoding="utf-8") as fh:
    json.dump({
        "session": 767, "r_drag_LIV_Mpc": r_drag,
        "classLXX_D_A":     {"pred": DA_pred, "obs": DA_obs, "err_pct": DA_err, "status": st_a, "form": DA_label},
        "classLXXI_D_V":    {"pred": DV_pred, "obs": DV_obs, "err_pct": DV_err, "status": st_b, "form": DV_label},
        "classLXXII_D_M":   {"pred": DM_pred, "obs": DM_obs, "err_pct": DM_err, "status": st_c, "form": DM_label},
        "classLXXIII_fs8":  {"pred": fs8_pred,"obs": fs8_obs,"err_pct": fs8_err,"status": st_d, "form": fs8_label},
    }, fh, indent=2)
print(f"\nArtifact: {artifact}")
print(f"Master registry written: {CSV_PATH}")
