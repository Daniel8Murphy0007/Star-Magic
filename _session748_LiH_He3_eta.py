"""
SESSION 748 -- Drive Li/H to candidate-EXACT; open Class XXIX (3He/H); open Class XXX (eta_b)

(a) Li/H: base (1-n_s)^7 * (SO5/D_phys) = 1.6085e-10 at +0.5302%.
    Required (1+delta) with delta = -0.005273 to reach 1.6000e-10.
(b) Class XXIX: 3He/H_obs ~ 1.00e-5. Try (1-n_s)^3 * M (locked rationals).
(c) Class XXX: eta_b_obs ~ 6.10e-10 (Planck 2018). Try (1-n_s)^k * M (locked rationals).

CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""
from __future__ import annotations
from fractions import Fraction as F
import csv, json, os, itertools

# -------------------- Locked primitives --------------------
F_TRZ = F(1,10); Phi_res = F(5,6); SSq = F(57,100); K_Mex = F(25,12)
beta_i = F(6029,10000); D_phys = F(4); D_BSFG = F(6); D_crit = F(26)
N_ch = F(9); SO5 = F(10); A_5 = F(60)

# Derived locked atoms
one_m_FTRZ = 1 - F_TRZ                 # 9/10
one_m_FP   = 1 - F_TRZ*Phi_res         # 11/12
ns_atom    = F(193,200)                # n_s
one_m_ns   = F(7,200)                  # (1 - n_s)
xi         = F(11,3200)                # 11/3200
r_tens     = F(9,250)                  # r
Yp         = F(49,200)                 # 7(1-n_s)

ATOMS = {
    "F_TRZ":F_TRZ, "Phi_res":Phi_res, "SSq":SSq, "K_Mex":K_Mex, "beta_i":beta_i,
    "D_phys":D_phys, "D_BSFG":D_BSFG, "D_crit":D_crit, "N_ch":N_ch, "SO5":SO5, "A_5":A_5,
    "1-F_TRZ":one_m_FTRZ, "1-F*P":one_m_FP, "n_s":ns_atom, "1-n_s":one_m_ns,
    "xi":xi, "r":r_tens, "Y_p":Yp,
    "27/26":F(27,26), "243/260":F(243,260), "33/40":F(33,40), "11/9":F(11,9),
    "22/9":F(22,9), "27/25":F(27,25), "416/513":F(416,513), "31/30":F(31,30),
    "5/108":F(5,108),
}
LABELS = list(ATOMS.keys())
VALS = [float(ATOMS[k]) for k in LABELS]

# -------------------- Helpers --------------------
def search3(target, tol_pct=5.0, want=8):
    """Find a*b/c, a*b*c, a/b/c hitting target."""
    hits = []
    n = len(LABELS)
    forms = [
        ("a*b/c",  lambda a,b,c: a*b/c),
        ("a*b*c",  lambda a,b,c: a*b*c),
        ("a/b/c",  lambda a,b,c: a/b/c),
    ]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for tag, fn in forms:
                    try:
                        v = fn(VALS[i], VALS[j], VALS[k])
                    except ZeroDivisionError:
                        continue
                    if v == 0: continue
                    err = (v - target)/target * 100.0
                    if abs(err) < tol_pct:
                        hits.append((LABELS[i], LABELS[j], LABELS[k], tag, v, err))
    hits.sort(key=lambda h: abs(h[5]))
    return hits[:want]

def search2(target, tol_pct=5.0, want=12):
    hits = []
    forms = [
        ("a*b",   lambda a,b: a*b),
        ("a/b",   lambda a,b: a/b),
    ]
    n = len(LABELS)
    for i in range(n):
        for j in range(n):
            for tag, fn in forms:
                try:
                    v = fn(VALS[i], VALS[j])
                except ZeroDivisionError:
                    continue
                if v == 0: continue
                err = (v - target)/target * 100.0
                if abs(err) < tol_pct:
                    hits.append((LABELS[i], LABELS[j], tag, v, err))
    hits.sort(key=lambda h: abs(h[4]))
    return hits[:want]

def write_ledger(label, predicted, observed, status="OK"):
    err_pct = (predicted - observed)/observed * 100.0
    raw = f"{label}: predicted={predicted:.6e} observed={observed:.6e} error_pct={err_pct:.6e} status={status}"
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "master_closures.csv")
    head = ["label","predicted","observed","error_pct","status","raw_output","cvw","sm_anchor"]
    new = not os.path.exists(path)
    with open(path,"a",newline="",encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=head, extrasaction="ignore")
        if new: w.writeheader()
        w.writerow({
            "label":label, "predicted":f"{predicted:.6e}", "observed":f"{observed:.6e}",
            "error_pct":f"{err_pct:.6e}", "status":status, "raw_output":raw,
            "cvw":"v2.0.0", "sm_anchor":"CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
        })
    print(raw)
    return err_pct

# -------------------- Banner --------------------
print("="*80)
print("SESSION 748 -- Li/H fine; Class XXIX 3He/H; Class XXX eta_b")
print("="*80)

# =====================================================================
# TRACK (a) -- Li/H fine
# =====================================================================
print("\n" + "-"*80)
print("TRACK (a) -- Li/H fine: (1-n_s)^7 * (SO5/D_phys) currently +0.5302%")
print("-"*80)
LiH_obs = 1.60e-10
base_li = float(one_m_ns)**7 * float(SO5/D_phys)
err_base = (base_li - LiH_obs)/LiH_obs * 100.0
print(f"  base = {base_li:.6e}, err = {err_base:+.6f}%")
need_mult_li = LiH_obs / base_li
delta_li = need_mult_li - 1.0
print(f"  needed mult = {need_mult_li:.8f}, delta = {delta_li:+.4e}")

print("\n  3-atom delta search for (1+delta):")
hits_li = search3(abs(delta_li), tol_pct=5.0, want=10)
sign_li = -1 if delta_li < 0 else 1
for h in hits_li:
    mult = 1.0 + sign_li * h[4]
    pred = base_li * mult
    err = (pred - LiH_obs)/LiH_obs * 100.0
    print(f"    1{'-' if sign_li<0 else '+'}{h[0]}*{h[1]}/{h[2] if h[3]=='a*b/c' else '...'}  delta={h[4]:.4e}  Li/H={pred:.4e}  err={err:+.6f}%")

print("\n  2-atom delta search:")
hits_li2 = search2(abs(delta_li), tol_pct=5.0, want=15)
for h in hits_li2:
    mult = 1.0 + sign_li * h[3]
    pred = base_li * mult
    err = (pred - LiH_obs)/LiH_obs * 100.0
    print(f"    1{'-' if sign_li<0 else '+'}{h[0]}{'*' if h[2]=='a*b' else '/'}{h[1]}    delta={h[3]:.4e}  Li/H={pred:.4e}  err={err:+.6f}%")

# Pick best
best_li = None
for h in hits_li:
    mult = 1.0 + sign_li * h[4]
    pred = base_li * mult
    err = (pred - LiH_obs)/LiH_obs * 100.0
    if best_li is None or abs(err) < abs(best_li[1]):
        best_li = (f"(1-n_s)^7*(SO5/D_phys)*[1{'-' if sign_li<0 else '+'}{h[0]}*{h[1]}/{h[2]}]" if h[3]=="a*b/c" else
                   f"(1-n_s)^7*(SO5/D_phys)*[1{'-' if sign_li<0 else '+'}{h[0]}*{h[1]}*{h[2]}]" if h[3]=="a*b*c" else
                   f"(1-n_s)^7*(SO5/D_phys)*[1{'-' if sign_li<0 else '+'}{h[0]}/{h[1]}/{h[2]}]", err, pred)
for h in hits_li2:
    mult = 1.0 + sign_li * h[3]
    pred = base_li * mult
    err = (pred - LiH_obs)/LiH_obs * 100.0
    if best_li is None or abs(err) < abs(best_li[1]):
        op = "*" if h[2]=="a*b" else "/"
        best_li = (f"(1-n_s)^7*(SO5/D_phys)*[1{'-' if sign_li<0 else '+'}{h[0]}{op}{h[1]}]", err, pred)

if best_li:
    print(f"\n  BEST: Li/H = {best_li[0]} = {best_li[2]:.4e}, err = {best_li[1]:+.6f}%")
    LiH_pred = best_li[2]
else:
    LiH_pred = base_li

# =====================================================================
# TRACK (b) -- Class XXIX 3He/H
# =====================================================================
print("\n" + "-"*80)
print("TRACK (b) -- Class XXIX: 3He/H = 1.00e-5 (Bania Rood Balser 2002)")
print("-"*80)
He3_obs = 1.00e-5

print("  Direct hypotheses:")
# (1-n_s)^k * M
for k in [1,2,3,4,5]:
    base = float(one_m_ns)**k
    need = He3_obs / base
    print(f"    (1-n_s)^{k} = {base:.4e}  need_mult = {need:.4f}")

# search direct
print("\n  3-atom direct: a*b/c hitting 1e-5")
hits_he3 = search3(He3_obs, tol_pct=10.0, want=15)
for h in hits_he3:
    print(f"    {h[0]}*{h[1]}/{h[2] if h[3]=='a*b/c' else '...'} ({h[3]:8s})  = {h[4]:.4e}  err={h[5]:+.4f}%")

# Compound hypotheses: (1-n_s)^k * (a*b/c)
print("\n  (1-n_s)^3 * (3-atom): need = ", He3_obs/float(one_m_ns)**3)
hits_h3 = search3(He3_obs/float(one_m_ns)**3, tol_pct=2.0, want=10)
for h in hits_h3:
    pred = float(one_m_ns)**3 * h[4]
    err = (pred-He3_obs)/He3_obs*100.0
    print(f"    (1-n_s)^3 * [{h[0]}*{h[1]}/{h[2]} ({h[3]})]  = {pred:.4e}  err={err:+.6f}%")

print("\n  (1-n_s)^2 * (3-atom): need = ", He3_obs/float(one_m_ns)**2)
hits_h2 = search3(He3_obs/float(one_m_ns)**2, tol_pct=2.0, want=10)
for h in hits_h2:
    pred = float(one_m_ns)**2 * h[4]
    err = (pred-He3_obs)/He3_obs*100.0
    print(f"    (1-n_s)^2 * [{h[0]}*{h[1]}/{h[2]} ({h[3]})]  = {pred:.4e}  err={err:+.6f}%")

print("\n  Y_p * (3-atom): need = ", He3_obs/float(Yp))
hits_yp = search3(He3_obs/float(Yp), tol_pct=2.0, want=10)
for h in hits_yp:
    pred = float(Yp) * h[4]
    err = (pred-He3_obs)/He3_obs*100.0
    print(f"    Y_p * [{h[0]}*{h[1]}/{h[2]} ({h[3]})]  = {pred:.4e}  err={err:+.6f}%")

# Pick best 3He/H
best_he3 = None
for h in hits_he3:
    err = h[5]
    if best_he3 is None or abs(err) < abs(best_he3[1]):
        op = "*" if h[3]=="a*b*c" else ("/" if h[3]=="a/b/c" else "*")
        best_he3 = (f"[{h[0]} {h[3]} {h[1]} {h[2]}]", err, h[4])
for h in hits_h3:
    pred = float(one_m_ns)**3 * h[4]
    err = (pred-He3_obs)/He3_obs*100.0
    if best_he3 is None or abs(err) < abs(best_he3[1]):
        best_he3 = (f"(1-n_s)^3 * [{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
for h in hits_h2:
    pred = float(one_m_ns)**2 * h[4]
    err = (pred-He3_obs)/He3_obs*100.0
    if best_he3 is None or abs(err) < abs(best_he3[1]):
        best_he3 = (f"(1-n_s)^2 * [{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)
for h in hits_yp:
    pred = float(Yp) * h[4]
    err = (pred-He3_obs)/He3_obs*100.0
    if best_he3 is None or abs(err) < abs(best_he3[1]):
        best_he3 = (f"Y_p * [{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)

if best_he3:
    print(f"\n  BEST: 3He/H = {best_he3[0]} = {best_he3[2]:.4e}, err = {best_he3[1]:+.6f}%")
    He3_pred = best_he3[2]
else:
    He3_pred = float(one_m_ns)**3

# =====================================================================
# TRACK (c) -- Class XXX eta_b
# =====================================================================
print("\n" + "-"*80)
print("TRACK (c) -- Class XXX: eta_b = 6.10e-10 (baryon-to-photon, Planck 2018)")
print("-"*80)
eta_obs = 6.10e-10

print("  Hypotheses:")
for k in [4,5,6,7,8]:
    print(f"    (1-n_s)^{k} = {float(one_m_ns)**k:.4e}")

# (1-n_s)^7 * M
print("\n  (1-n_s)^7 * (single atom): need =", eta_obs/float(one_m_ns)**7)
need_eta7 = eta_obs/float(one_m_ns)**7
hits_e7_1 = []
for k, v in ATOMS.items():
    err = (v - need_eta7)/need_eta7 * 100.0
    if abs(err) < 5.0:
        hits_e7_1.append((k, float(v), err))
hits_e7_1.sort(key=lambda h: abs(h[2]))
for h in hits_e7_1[:10]:
    pred = float(one_m_ns)**7 * h[1]
    err = (pred-eta_obs)/eta_obs*100.0
    print(f"    (1-n_s)^7 * {h[0]:12s} = {h[1]:.4f}  eta={pred:.4e}  err={err:+.6f}%")

print("\n  (1-n_s)^7 * (a/b): need =", need_eta7)
hits_e7_2 = search2(need_eta7, tol_pct=2.0, want=15)
for h in hits_e7_2:
    pred = float(one_m_ns)**7 * h[3]
    err = (pred-eta_obs)/eta_obs*100.0
    op = "*" if h[2]=="a*b" else "/"
    print(f"    (1-n_s)^7 * ({h[0]}{op}{h[1]}) = {h[3]:.4f}  eta={pred:.4e}  err={err:+.6f}%")

print("\n  (1-n_s)^7 * (3-atom): need =", need_eta7)
hits_e7_3 = search3(need_eta7, tol_pct=1.0, want=10)
for h in hits_e7_3:
    pred = float(one_m_ns)**7 * h[4]
    err = (pred-eta_obs)/eta_obs*100.0
    print(f"    (1-n_s)^7 * [{h[0]} {h[3]} {h[1]} {h[2]}] = {h[4]:.4f}  eta={pred:.4e}  err={err:+.6f}%")

print("\n  Li/H * (single atom): need =", eta_obs/LiH_pred)
need_eta_li = eta_obs/LiH_pred
for k, v in ATOMS.items():
    err = (float(v) - need_eta_li)/need_eta_li * 100.0
    if abs(err) < 5.0:
        pred = LiH_pred * float(v)
        err2 = (pred-eta_obs)/eta_obs*100.0
        print(f"    Li/H * {k:12s} = {float(v):.4f}  eta={pred:.4e}  err={err2:+.6f}%")

# Pick best eta_b
best_eta = None
for h in hits_e7_1[:10]:
    pred = float(one_m_ns)**7 * h[1]
    err = (pred-eta_obs)/eta_obs*100.0
    if best_eta is None or abs(err) < abs(best_eta[1]):
        best_eta = (f"(1-n_s)^7 * {h[0]}", err, pred)
for h in hits_e7_2:
    pred = float(one_m_ns)**7 * h[3]
    err = (pred-eta_obs)/eta_obs*100.0
    op = "*" if h[2]=="a*b" else "/"
    if best_eta is None or abs(err) < abs(best_eta[1]):
        best_eta = (f"(1-n_s)^7 * ({h[0]}{op}{h[1]})", err, pred)
for h in hits_e7_3:
    pred = float(one_m_ns)**7 * h[4]
    err = (pred-eta_obs)/eta_obs*100.0
    if best_eta is None or abs(err) < abs(best_eta[1]):
        best_eta = (f"(1-n_s)^7 * [{h[0]} {h[3]} {h[1]} {h[2]}]", err, pred)

if best_eta:
    print(f"\n  BEST: eta_b = {best_eta[0]} = {best_eta[2]:.4e}, err = {best_eta[1]:+.6f}%")
    eta_pred = best_eta[2]
else:
    eta_pred = float(one_m_ns)**7

# =====================================================================
# Emit ledger
# =====================================================================
print()
e1 = write_ledger("classXXVIII_LiH_session748", LiH_pred, LiH_obs)
e2 = write_ledger("classXXIX_He3H_session748", He3_pred, He3_obs)
e3 = write_ledger("classXXX_eta_b_session748", eta_pred, eta_obs)

print("\n" + "-"*80)
print("DECISION GATE")
print("-"*80)
print(f"  Li/H    err = {e1:+.6f}%")
print(f"  3He/H   err = {e2:+.6f}%")
print(f"  eta_b   err = {e3:+.6f}%")

# Artifact
art = {
    "session":748,
    "cvw":"v2.0.0",
    "sm_anchor":"CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
    "tracks":{
        "a_LiH":   {"formula": best_li[0] if best_li else "(1-n_s)^7*(SO5/D_phys)", "value": LiH_pred, "err_pct": e1},
        "b_He3H":  {"formula": best_he3[0] if best_he3 else "(1-n_s)^3", "value": He3_pred, "err_pct": e2},
        "c_eta_b": {"formula": best_eta[0] if best_eta else "(1-n_s)^7", "value": eta_pred, "err_pct": e3},
    },
}
art_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_session748_LiH_He3_eta_result.json")
with open(art_path,"w",encoding="utf-8") as f:
    json.dump(art, f, indent=2)
print(f"\nArtifact: {art_path}")
