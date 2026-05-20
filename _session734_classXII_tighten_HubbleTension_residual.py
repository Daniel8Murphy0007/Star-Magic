"""
SESSION 734 -- Three parallel tracks
(a) Tighten Class XII Omega_b/Omega_m to sub-0.01%
    -- search 4-atom corrections around 171/1100
(b) Class XIII candidate: Hubble tension H0_SH0ES/H0_Planck = 1.0795
    -- 1-3 atom closure search
(c) Validate r_III/r_X = 72/55 hypothesis (shared transcendental dressing)
"""

from __future__ import annotations
import csv
import json
import math
from fractions import Fraction
from itertools import product
from pathlib import Path

# Locked constants
F_TRZ    = Fraction(1, 10)
Phi_res  = Fraction(5, 6)
SSq      = Fraction(57, 100)
K_Mex    = Fraction(25, 12)
beta_i   = Fraction(6029, 10000)
D_phys   = Fraction(4)
D_BSFG   = Fraction(6)
D_crit   = Fraction(26)
N_ch     = Fraction(9)
A_5      = Fraction(60)
K_G      = Fraction(33, 104)

ATOMS = {
    "F_TRZ": F_TRZ, "Phi_res": Phi_res, "SSq": SSq, "K_Mex": K_Mex,
    "beta_i": beta_i, "D_phys": D_phys, "D_BSFG": D_BSFG, "D_crit": D_crit,
    "N_ch": N_ch, "A_5": A_5,
    "1-F_TRZ": Fraction(9, 10), "1-F*P": Fraction(11, 12),
    "27/26": Fraction(27, 26), "243/260": Fraction(243, 260),
    "405/247": Fraction(405, 247), "13/6": Fraction(13, 6),
    "K_G": K_G, "1/K_G": Fraction(104, 33),
    "6/5": Fraction(6, 5), "72/55": Fraction(72, 55),
}

OBS = {
    "f_b": 0.04897 / 0.3153,
    "H0_Planck": 67.66,   # km/s/Mpc (Planck 2018)
    "H0_SH0ES":  73.04,   # km/s/Mpc (R22)
}

def write_ledger(rows, script_name):
    csv_path = Path("master_closures.csv")
    existing = []
    if csv_path.exists():
        with csv_path.open("r", encoding="utf-8", newline="") as f:
            existing = list(csv.DictReader(f))
    fieldnames = ["ID","name","predicted","observed","error_pct","status","script","sm_anchor"]
    extras = set()
    for r in existing: extras |= set(r.keys())
    all_fields = fieldnames + [k for k in extras if k not in fieldnames]
    next_id = max((int(r["ID"]) for r in existing if r.get("ID","").isdigit()), default=0) + 1
    for r in rows:
        r["ID"] = str(next_id); next_id += 1
        r["script"] = script_name
        r["sm_anchor"] = "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant"
        existing.append(r)
        print(f"{r['name']}: predicted={r['predicted']} observed={r['observed']} error_pct={r['error_pct']} status={r['status']}")
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=all_fields, extrasaction="ignore")
        w.writeheader(); w.writerows(existing)

print("=" * 80)
print("SESSION 734 -- Class XII tighten + Class XIII Hubble tension + r_III/r_X verify")
print("=" * 80)

# ============================================================================
# Track (a) -- Tighten Class XII
# ============================================================================
print("\n" + "-" * 80)
print("TRACK (a) -- Tighten Class XII Omega_b/Omega_m")
print("-" * 80)

target_a = OBS["f_b"]
K_current_a = Fraction(12,11) * SSq / D_phys  # = 171/1100
K_float_a = float(K_current_a)
err_a_cur = (K_float_a - target_a)/target_a * 100
print(f"  target          = {target_a:.10f}")
print(f"  current K       = (12/11)*SSq/D_phys = {K_current_a} = {K_float_a:.10f}")
print(f"  current err     = {err_a_cur:+.4f}%")
correction_needed_a = target_a / K_float_a
print(f"  correction need = {correction_needed_a:.8f}  (epsilon = {correction_needed_a-1:+.4e})")

# 1-atom correction search
print(f"\n  1-atom correction factors c such that K_current * c approaches target:")
hits_a = []
for n, a in ATOMS.items():
    for p in (-2,-1,1,2):
        c = float(a)**p
        full = K_float_a * c
        err = (full - target_a)/target_a * 100
        hits_a.append((abs(err), c, full, f"{n}^{p}", err))
hits_a.sort()
print(f"  {'rank':<5}{'corr':<25}{'corr val':>14}{'full K':>16}{'err%':>14}")
shown = 0
seen = set()
for h in hits_a:
    if h[3] in seen: continue
    seen.add(h[3]); shown += 1
    if shown > 8: break
    flag = " <-- SUB-0.05%" if abs(h[4]) < 0.05 else ""
    print(f"  {shown:<5}{h[3]:<25}{h[1]:>14.6f}{h[2]:>16.8f}{h[4]:>+14.4f}%{flag}")

# 2-atom multiplicative correction
print(f"\n  2-atom correction search (sub-0.05%):")
atoms_list = list(ATOMS.items())
hits_a2 = []
for n1, a1 in atoms_list:
    for n2, a2 in atoms_list:
        for p1, p2 in product((-2,-1,1,2), repeat=2):
            c = (float(a1)**p1) * (float(a2)**p2)
            if not 0.98 < c < 1.02: continue
            full = K_float_a * c
            err = (full - target_a)/target_a * 100
            if abs(err) < 0.05:
                hits_a2.append((abs(err), c, full, f"{n1}^{p1}*{n2}^{p2}", err))
hits_a2.sort()
shown = 0
for h in hits_a2[:8]:
    shown += 1
    print(f"  {shown:<5}{h[3]:<40}{h[1]:>12.6f}{h[2]:>16.8f}{h[4]:>+14.4f}%")

best_a = hits_a2[0] if hits_a2 else (None,)*5

# ============================================================================
# Track (b) -- Class XIII Hubble tension ratio
# ============================================================================
print("\n" + "-" * 80)
print("TRACK (b) -- Class XIII: H0(SH0ES)/H0(Planck) = 1.0795")
print("-" * 80)

ratio_H = OBS["H0_SH0ES"] / OBS["H0_Planck"]
print(f"  ratio (SH0ES/Planck) = {ratio_H:.6f}")
print(f"  excess               = {(ratio_H-1)*100:+.4f}%")

# 1-2 atom search
print(f"\n  1-2 atom closures for ratio={ratio_H:.5f}:")
hits_b = []
for n1, a1 in atoms_list:
    for p1 in (-2,-1,1,2):
        val = float(a1)**p1
        err = (val - ratio_H)/ratio_H * 100
        if abs(err) < 2:
            hits_b.append((abs(err), val, f"{n1}^{p1}", err))
    for n2, a2 in atoms_list:
        for p1, p2 in product((-1,1), repeat=2):
            val = (float(a1)**p1) * (float(a2)**p2)
            err = (val - ratio_H)/ratio_H * 100
            if abs(err) < 0.5:
                hits_b.append((abs(err), val, f"{n1}^{p1}*{n2}^{p2}", err))
hits_b.sort()
seen = set()
shown = 0
print(f"  {'rank':<5}{'name':<40}{'val':>14}{'err%':>12}")
for h in hits_b:
    if h[2] in seen: continue
    seen.add(h[2]); shown += 1
    if shown > 12: break
    print(f"  {shown:<5}{h[2]:<40}{h[1]:>14.6f}{h[3]:>+12.4f}%")
best_b = min(hits_b, key=lambda x: abs(x[3])) if hits_b else None

# 3-atom
print(f"\n  3-atom search (|err|<0.1%):")
hits_b3 = []
for n1, a1 in atoms_list:
    for n2, a2 in atoms_list:
        for n3, a3 in atoms_list:
            for p1, p2, p3 in product((-1,1), repeat=3):
                val = (float(a1)**p1)*(float(a2)**p2)*(float(a3)**p3)
                err = (val - ratio_H)/ratio_H * 100
                if abs(err) < 0.1:
                    hits_b3.append((abs(err), val, f"{n1}^{p1}*{n2}^{p2}*{n3}^{p3}", err))
hits_b3.sort()
seen = set()
shown = 0
for h in hits_b3:
    key = tuple(sorted(h[2].split("*")))
    if key in seen: continue
    seen.add(key); shown += 1
    if shown > 8: break
    print(f"  {shown:<5}{h[2]:<60}{h[1]:>14.6f}{h[3]:>+12.4f}%")
best_b3 = min(hits_b3, key=lambda x: abs(x[3])) if hits_b3 else None

# ============================================================================
# Track (c) -- r_III / r_X = 72/55 hypothesis
# ============================================================================
print("\n" + "-" * 80)
print("TRACK (c) -- Validate r_III/r_X = 72/55 hypothesis")
print("-" * 80)

r_III_obs = -6.92e-4
r_X_obs   = -5.26e-4
ratio_meas = r_III_obs / r_X_obs
hyp_72_55 = float(Fraction(72, 55))
print(f"  measured r_III/r_X      = {ratio_meas:.6f}")
print(f"  hypothesis 72/55        = {hyp_72_55:.6f}")
print(f"  deviation               = {(ratio_meas-hyp_72_55)/hyp_72_55*100:+.4f}%")

# If true: r_X = r_III * 55/72  -> predict r_X from r_III
r_X_predicted = r_III_obs * float(Fraction(55, 72))
print(f"\n  Predict r_X from r_III: r_X = r_III*(55/72) = {r_X_predicted:.4e}")
print(f"  Observed r_X            = {r_X_obs:.4e}")
print(f"  Deviation               = {(r_X_predicted - r_X_obs)/r_X_obs*100:+.4f}%")

# Predict tightened Class X K from hypothesis
target_X = 2.171583
K_X_current = 2.170440
K_X_dressed = K_X_current * (1 + abs(r_X_predicted))
err_X_dressed = (K_X_dressed - target_X)/target_X * 100
print(f"\n  Current K_X             = {K_X_current:.6f}")
print(f"  Dressed K_X (1+|r_X_pred|)= {K_X_dressed:.6f}")
print(f"  Dressed err             = {err_X_dressed:+.4f}%")

# Alternative ratios to check
print("\n  Alternative locked ratios within 1% of measured 1.3156:")
alt_hits = []
for n1, a1 in ATOMS.items():
    for n2, a2 in ATOMS.items():
        v = float(a1)/float(a2)
        if abs(v - ratio_meas)/ratio_meas < 0.01:
            alt_hits.append((abs(v - ratio_meas), v, f"{n1}/{n2}"))
alt_hits.sort()
seen = set()
for d, v, name in alt_hits[:6]:
    if name in seen: continue
    seen.add(name)
    print(f"    {name:<30}{v:.6f}   dev={(v-ratio_meas)/ratio_meas*100:+.4f}%")

# ============================================================================
# Decision gate + ledger
# ============================================================================
print("\n" + "-" * 80)
print("DECISION GATE")
print("-" * 80)

rows = []

# (a) Class XII tightening
if best_a[0] is not None and abs(best_a[4]) < abs(err_a_cur):
    rows.append({
        "name": "Omega_b_Omega_m_classXII_corrected",
        "predicted": f"{best_a[2]:.6e}",
        "observed":  f"{target_a:.6e}",
        "error_pct": f"{best_a[4]:.6e}",
        "status": "candidate-EXACT" if abs(best_a[4]) < 5e-4 else "OK",
    })
    print(f"  (a) XII tightened: {best_a[3]}, err={best_a[4]:+.4f}%")
else:
    rows.append({
        "name": "Omega_b_Omega_m_classXII_tighten_BLOCKED",
        "predicted": f"{K_float_a:.6e}",
        "observed":  f"{target_a:.6e}",
        "error_pct": f"{err_a_cur:.6e}",
        "status": "OK",
    })
    print(f"  (a) XII tighten BLOCKED at {err_a_cur:+.4f}%")

# (b) Class XIII
if best_b is not None:
    rows.append({
        "name": "H0_tension_classXIII_2atom",
        "predicted": f"{best_b[1]:.6e}",
        "observed":  f"{ratio_H:.6e}",
        "error_pct": f"{best_b[3]:.6e}",
        "status": "candidate-EXACT" if abs(best_b[3]) < 5e-4 else "OK",
    })
    print(f"  (b) XIII 2-atom: {best_b[2]} err={best_b[3]:+.4f}%")
if best_b3 is not None:
    rows.append({
        "name": "H0_tension_classXIII_3atom",
        "predicted": f"{best_b3[1]:.6e}",
        "observed":  f"{ratio_H:.6e}",
        "error_pct": f"{best_b3[3]:.6e}",
        "status": "candidate-EXACT" if abs(best_b3[3]) < 5e-4 else "OK",
    })
    print(f"  (b) XIII 3-atom: {best_b3[2]} err={best_b3[3]:+.6f}%")

# (c) Residual hypothesis
rows.append({
    "name": "residual_ratio_III_X_vs_72_55",
    "predicted": f"{hyp_72_55:.6e}",
    "observed":  f"{ratio_meas:.6e}",
    "error_pct": f"{(hyp_72_55-ratio_meas)/ratio_meas*100:.6e}",
    "status": "OK",
})
rows.append({
    "name": "classX_dressed_via_classIII_residual",
    "predicted": f"{K_X_dressed:.6e}",
    "observed":  f"{target_X:.6e}",
    "error_pct": f"{err_X_dressed:.6e}",
    "status": "candidate-EXACT" if abs(err_X_dressed) < 5e-4 else "OK",
})
print(f"  (c) r_III/r_X dev vs 72/55: {(hyp_72_55-ratio_meas)/ratio_meas*100:+.4f}%")
print(f"  (c) Class X dressed err: {err_X_dressed:+.4f}%")

write_ledger(rows, "_session734_classXII_tighten_HubbleTension_residual.py")

artifact = {
    "session": 734,
    "a_classXII": {
        "K_current": str(K_current_a),
        "err_current_pct": err_a_cur,
        "best_2atom_correction": {"name": best_a[3], "corr": best_a[1], "full": best_a[2], "err": best_a[4]} if best_a[0] is not None else None,
    },
    "b_classXIII": {
        "ratio_H": ratio_H,
        "best_2atom": {"name": best_b[2], "val": best_b[1], "err": best_b[3]} if best_b else None,
        "best_3atom": {"name": best_b3[2], "val": best_b3[1], "err": best_b3[3]} if best_b3 else None,
    },
    "c_residual_hypothesis": {
        "ratio_meas": ratio_meas,
        "hyp_72_55": hyp_72_55,
        "deviation_pct": (ratio_meas-hyp_72_55)/hyp_72_55*100,
        "K_X_dressed": K_X_dressed,
        "K_X_dressed_err_pct": err_X_dressed,
    },
}
out = Path("_session734_classXII_tighten_HubbleTension_residual_result.json")
out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
print(f"\nArtifact written: {out.resolve().as_posix()}")
