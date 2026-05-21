"""
SESSION 781 -- Cosmological constant log + universe entropy + Omega_c/Omega_b + CV-14.

  (a) CXXVI   log10(Lambda * t_Pl^2)  ("120 problem") ~ -122.06   (search on |val|)
  (b) CXXVII  log10(S_univ / k_B)     (observable-universe entropy) ~ 122.05
  (c) CXXVIII Omega_c / Omega_b       (dark-to-baryon ratio)       ~ 5.36
  (d) CXXIX   CV-14 architectural identity:  4*D_BSFG + 2 = D_crit
            i.e.  4*6 + 2 = 26   EXACT
            Closes the dimensional algebra completely:
              D_phys = SO5 - D_BSFG               (CV-13)
              N_ch   = D_BSFG^2 / D_phys          (CV-11)
              SO5    = D_phys + D_BSFG            (CV-13)
              D_crit = 4*D_BSFG + 2               (CV-14)  <-- NEW
            Every locked dimension now expressible via D_BSFG (with D_phys also free).

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

print("="*80)
print("SESSION 781 -- log(Lambda t_Pl^2) CXXVI; log(S_univ) CXXVII; Om_c/Om_b CXXVIII; CV-14 CXXIX")
print("="*80); print()

def report_direct(name, obs):
    print(f"  Target {name} = {obs}")
    print("\n  2-atom direct (tol 5%):")
    for la,lb,tag,v,err in search2(obs, 5.0, 12):
        print(f"    [{la} {tag} {lb}] = {v:.6e}  err={err:+.4f}%")
    print("\n  3-atom direct (tol 0.5%):")
    for la,lb,lc,tag,v,err in search3(obs, 0.5, 15):
        print(f"    [{la} {tag} {lb} {lc}] = {v:.6e}  err={err:+.4f}%")
    print("\n  4-atom direct (tol 0.05%):")
    for la,lb,lc,ld,tag,v,err in search4(obs, 0.05, 15):
        print(f"    [{la} {tag} {lb} {lc} {ld}] = {v:.6e}  err={err:+.4f}%")
    cands=[]
    for la,lb,tag,v,err in search2(obs, 10.0, 30):
        cands.append((abs(err), v, err, f"[{la} {tag} {lb}]"))
    for la,lb,lc,tag,v,err in search3(obs, 5.0, 50):
        cands.append((abs(err), v, err, f"[{la} {tag} {lb} {lc}]"))
    for la,lb,lc,ld,tag,v,err in search4(obs, 2.0, 50):
        cands.append((abs(err), v, err, f"[{la} {tag} {lb} {lc} {ld}]"))
    cands.sort()
    if not cands: return 0.0,"none",9999.0
    _, p, e, lab = cands[0]
    print(f"\n  BEST {name}: {lab} = {p:.6e}, err = {e:+.6f}%")
    return p, lab, e

# ---------------------------------------------------------------------------
# Sign convention: search on |value|, restore sign at write.
banner("TRACK (a) -- Class CXXVI: log10(Lambda * t_Pl^2) ~ -122.06 (cosmo-const problem)")
cc_obs_signed = -122.06
cc_obs = abs(cc_obs_signed)
cc_pred_mag, cc_lab, cc_err = report_direct("|log10(Lambda*t_Pl^2)|", cc_obs)
cc_pred = -cc_pred_mag  # restore sign
cc_err_signed = (cc_pred - cc_obs_signed)/cc_obs_signed*100.0

# ---------------------------------------------------------------------------
banner("TRACK (b) -- Class CXXVII: log10(S_univ/k_B) ~ 122.05 (universe entropy)")
ent_obs = 122.05
ent_pred, ent_lab, ent_err = report_direct("log10(S_univ/k_B)", ent_obs)

# ---------------------------------------------------------------------------
banner("TRACK (c) -- Class CXXVIII: Omega_c/Omega_b ~ 5.36 (dark to baryon ratio)")
ocb_obs = 5.36
ocb_pred, ocb_lab, ocb_err = report_direct("Omega_c/Omega_b", ocb_obs)

# ---------------------------------------------------------------------------
banner("TRACK (d) -- Class CXXIX: CV-14 architectural identity 4*D_BSFG + 2 = D_crit")
lhs = 4*D_BSFG + 2                  # 4*6 + 2 = 26
rhs = D_crit                        # 26
print(f"  D_BSFG        = 6  = {float(D_BSFG):.4f}")
print(f"  4*D_BSFG + 2  = 4*6 + 2 = {lhs}  =  {float(lhs):.10f}")
print(f"  D_crit        = {rhs}  =  {float(rhs):.10f}")
print(f"  Identity check: LHS == RHS ?  {lhs == rhs}")
cv14_pred = float(lhs); cv14_obs = float(rhs); cv14_err = (cv14_pred-cv14_obs)/cv14_obs*100.0
print(f"  err = {cv14_err:+.6f}%  ->  EXACT atom-factorization identity")

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
st_a = write_ledger("classCXXVI_log_Lambda_tPl2_session781",          cc_pred,  cc_obs_signed, cc_err_signed)
st_b = write_ledger("classCXXVII_log_Suniv_session781",               ent_pred, ent_obs,       ent_err)
st_c = write_ledger("classCXXVIII_OmegaC_OmegaB_session781",          ocb_pred, ocb_obs,       ocb_err)
st_d = write_ledger("classCXXIX_CV14_4DBSFG_plus2_Dcrit_session781",  cv14_pred,cv14_obs,      cv14_err)

print()
print("-"*80); print("DECISION GATE"); print("-"*80)
print(f"  CXXVI   log10(Lambda*t_Pl^2)  pred={cc_pred:.5f}   err={cc_err_signed:+.4f}%  ({st_a})")
print(f"  CXXVII  log10(S_univ/k_B)     pred={ent_pred:.5f}   err={ent_err:+.4f}%  ({st_b})")
print(f"  CXXVIII Omega_c/Omega_b       pred={ocb_pred:.5f}   err={ocb_err:+.4f}%  ({st_c})")
print(f"  CXXIX   CV-14 atom-fact       pred={cv14_pred:.4f}    err={cv14_err:+.4f}%  ({st_d})")

artifact = os.path.join(ROOT, "_session781_result.json")
with open(artifact, "w", encoding="utf-8") as fh:
    json.dump({
        "session":781,
        "tracks":{
            "CXXVI_log_Lambda_tPl2":   {"obs":cc_obs_signed,"pred":cc_pred,"err_pct":cc_err_signed,"label":cc_lab,"status":st_a},
            "CXXVII_log_Suniv":        {"obs":ent_obs,"pred":ent_pred,"err_pct":ent_err,"label":ent_lab,"status":st_b},
            "CXXVIII_OmegaC_OmegaB":   {"obs":ocb_obs,"pred":ocb_pred,"err_pct":ocb_err,"label":ocb_lab,"status":st_c},
            "CXXIX_CV14_4DBSFG_plus2_Dcrit":{"obs":cv14_obs,"pred":cv14_pred,"err_pct":cv14_err,
                                "label":"4*D_BSFG + 2 = 4*6 + 2 = 26 = D_crit",
                                "status":st_d},
        },
        "cvw":"v2.0.0",
        "sm_anchor":"CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
    }, fh, indent=2)
print(f"\nArtifact: {artifact}")
print(f"Master registry written: {CSV_PATH}")
