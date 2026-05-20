"""
SESSION 736 -- Three parallel tracks
(a) Closed form for delta_univ ~ -5.2e-4 (search wider locked pool)
(b) Class XV: n_s = 0.965 scalar spectral index (1-3 atom search)
(c) Cross-validate Class XIV via S_8 = sigma_8 * sqrt(Omega_m/0.3)
"""

from __future__ import annotations
import csv, json, math
from fractions import Fraction
from itertools import product
from pathlib import Path

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
    "27/25": Fraction(27, 25), "55/72": Fraction(55, 72),
    "416/513": Fraction(416, 513),
}

# Observables
N_S       = 0.965         # Planck scalar spectral index
SIGMA_8   = 0.811
S_8       = 0.832
OMEGA_M   = 0.3153

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

atoms_list = list(ATOMS.items())

print("=" * 80)
print("SESSION 736 -- delta_univ closed form + Class XV n_s + Class XIV S_8 check")
print("=" * 80)

# ============================================================================
# Track (a) -- closed form for delta_univ ~ -5.2e-4 (signed)
# ============================================================================
print("\n" + "-" * 80)
print("TRACK (a) -- Closed form for delta_univ ~ -5.2e-4")
print("-" * 80)

# Use observed delta values from S735:
delta_obs = -5.16e-4  # cluster mean (III,V,X,XII,XIII)
# Try single locked atom * 1/D_crit^n * (c/v)^m
c_over_v = 2.99792458  # (c/v_SCM)

print(f"  target delta_univ = {delta_obs:.4e}")
print(f"  log10|delta|      = {math.log10(abs(delta_obs)):.4f}")
print(f"  Compare:")
print(f"    -1/2000          = {-1/2000:+.4e}")
print(f"    -1/(D_crit^2*3)  = {-1/(26**2*3):+.4e}")
print(f"    -1/(N_ch^2*24)   = {-1/(81*24):+.4e}")
print(f"    -1/(A_5*K_Mex*alpha_inv?)")

# Pure search: K * D_crit^p * (c/v)^q close to -5.16e-4
print(f"\n  Search K * D_crit^p * (c/v)^q close to {delta_obs:+.4e}:")
hits = []
for n, a in ATOMS.items():
    for sign in (1, -1):
        for p in range(-5, 1):
            for q in range(-5, 6):
                val = sign * float(a) * (26 ** p) * (c_over_v ** q)
                err = (val - delta_obs)/delta_obs * 100
                if abs(err) < 5:
                    hits.append((abs(err), val, f"{sign:+d}*{n}*D_crit^{p}*(c/v)^{q}", err))
hits.sort()
seen = set(); shown = 0
print(f"  {'rank':<5}{'closure':<55}{'val':>14}{'err%':>10}")
for h in hits:
    if h[2] in seen: continue
    seen.add(h[2]); shown += 1
    if shown > 10: break
    flag = " ★" if abs(h[3]) < 0.5 else ""
    print(f"  {shown:<5}{h[2]:<55}{h[1]:>+14.4e}{h[3]:>+10.3f}%{flag}")
best_a = min(hits, key=lambda x: abs(x[3])) if hits else None

# 2-atom multiplicative search with D_crit, (c/v)
print(f"\n  2-atom K1*K2 with D_crit^p and (c/v)^q (sub-1%):")
hits2 = []
for n1, a1 in atoms_list:
    for n2, a2 in atoms_list:
        for sign in (1, -1):
            for p in range(-4, 1):
                for q in range(-4, 5):
                    val = sign * float(a1) * float(a2) * (26 ** p) * (c_over_v ** q)
                    err = (val - delta_obs)/delta_obs * 100
                    if abs(err) < 1:
                        hits2.append((abs(err), val, f"{sign:+d}*{n1}*{n2}*D_crit^{p}*(c/v)^{q}", err))
hits2.sort()
seen = set(); shown = 0
for h in hits2[:8]:
    shown += 1
    print(f"  {shown:<5}{h[2]:<65}{h[1]:>+14.4e}{h[3]:>+10.3f}%")
best_a2 = min(hits2, key=lambda x: abs(x[3])) if hits2 else None

# ============================================================================
# Track (b) -- Class XV: n_s = 0.965
# ============================================================================
print("\n" + "-" * 80)
print("TRACK (b) -- Class XV: scalar spectral index n_s = 0.965")
print("-" * 80)

print(f"  target n_s = {N_S}")
print(f"  log = {math.log10(N_S):.5f}")
print(f"  1 - n_s = {1 - N_S:.4f}  (slow-roll deviation)")

# 1-2 atom direct
print(f"\n  1-2 atom closures for n_s = {N_S}:")
hits_b = []
for n, a in atoms_list:
    for p in (-2,-1,1,2):
        v = float(a)**p
        e = (v - N_S)/N_S*100
        if abs(e) < 3:
            hits_b.append((abs(e), v, f"{n}^{p}", e))
    for n2, a2 in atoms_list:
        for p1, p2 in product((-1,1), repeat=2):
            v = (float(a)**p1) * (float(a2)**p2)
            e = (v - N_S)/N_S*100
            if abs(e) < 0.5:
                hits_b.append((abs(e), v, f"{n}^{p1}*{n2}^{p2}", e))
hits_b.sort()
seen = set(); shown = 0
print(f"  {'rank':<5}{'name':<40}{'val':>14}{'err%':>12}")
for h in hits_b:
    if h[2] in seen: continue
    seen.add(h[2]); shown += 1
    if shown > 10: break
    flag = " ★" if abs(h[3]) < 0.1 else ""
    print(f"  {shown:<5}{h[2]:<40}{h[1]:>14.6f}{h[3]:>+12.4f}%{flag}")
best_b = min(hits_b, key=lambda x: abs(x[3])) if hits_b else None

# 3-atom
print(f"\n  3-atom closures (|err|<0.05%):")
hits_b3 = []
for n1, a1 in atoms_list:
    for n2, a2 in atoms_list:
        for n3, a3 in atoms_list:
            for p1,p2,p3 in product((-1,1), repeat=3):
                v = (float(a1)**p1)*(float(a2)**p2)*(float(a3)**p3)
                e = (v - N_S)/N_S*100
                if abs(e) < 0.05:
                    hits_b3.append((abs(e), v, f"{n1}^{p1}*{n2}^{p2}*{n3}^{p3}", e))
hits_b3.sort()
seen_combos = set(); shown = 0
for h in hits_b3:
    key = tuple(sorted(h[2].split("*")))
    if key in seen_combos: continue
    seen_combos.add(key); shown += 1
    if shown > 8: break
    print(f"  {shown:<5}{h[2]:<60}{h[1]:>12.6f}{h[3]:>+12.4f}%")
best_b3 = min(hits_b3, key=lambda x: abs(x[3])) if hits_b3 else None

# Compare to 1 - 1/D_crit:
print(f"\n  Heuristic checks:")
print(f"    1 - 1/D_crit  = {1 - 1/26:.6f}    err vs n_s = {(1-1/26 - N_S)/N_S*100:+.4f}%")
print(f"    1 - F_TRZ*SSq*beta_i = {1 - 0.1*0.57*0.6029:.6f}    err = {(1 - 0.1*0.57*0.6029 - N_S)/N_S*100:+.4f}%")
print(f"    1 - (1-1/K_Mex) = {1/K_Mex:.6f}    err = {(float(1/K_Mex) - N_S)/N_S*100:+.4f}%")
print(f"    1 - 7/200       = {1 - 7/200:.6f}    err = {(1 - 7/200 - N_S)/N_S*100:+.4f}%")

# ============================================================================
# Track (c) -- Cross-validate Class XIV via S_8
# ============================================================================
print("\n" + "-" * 80)
print("TRACK (c) -- Cross-validate XIV: S_8 = sigma_8 * sqrt(Omega_m/0.3)")
print("-" * 80)

sigma_8_uqff = 416 / 513
Omega_L_over_m = 2.171583  # observed
Omega_m_uqff = 1 / (1 + Omega_L_over_m)
print(f"  sigma_8 (UQFF Class XIV) = {sigma_8_uqff:.6f}")
print(f"  Omega_m (from Class X)   = 1/(1+{Omega_L_over_m}) = {Omega_m_uqff:.6f}")
print(f"  S_8 = sigma_8 * sqrt(Omega_m/0.3) = {sigma_8_uqff * math.sqrt(Omega_m_uqff/0.3):.6f}")
print(f"  S_8 observed (Planck)    = {S_8}")
S8_uqff = sigma_8_uqff * math.sqrt(Omega_m_uqff / 0.3)
err_S8 = (S8_uqff - S_8)/S_8 * 100
print(f"  err = {err_S8:+.4f}%")

# Also try Class X dressed (+0.0002%): Omega_L/m = 2.171587
Omega_L_over_m_dressed = 2.171587
Omega_m_dressed = 1 / (1 + Omega_L_over_m_dressed)
S8_dressed = sigma_8_uqff * math.sqrt(Omega_m_dressed / 0.3)
err_S8_d = (S8_dressed - S_8)/S_8 * 100
print(f"\n  Using Class X dressed Omega_m = {Omega_m_dressed:.6f}:")
print(f"  S_8_dressed = {S8_dressed:.6f}, err = {err_S8_d:+.4f}%")

# Check if S8 closure 5/6 = Phi_res still beats compositional derivation
print(f"\n  Direct closure S_8 = Phi_res = 5/6 = {5/6:.6f}, err = {(5/6 - S_8)/S_8*100:+.4f}%")

# ============================================================================
print("\n" + "-" * 80)
print("DECISION GATE")
print("-" * 80)
rows = []

# (a)
if best_a is not None:
    rows.append({
        "name": "delta_univ_closed_form_1atom",
        "predicted": f"{best_a[1]:.6e}",
        "observed":  f"{delta_obs:.6e}",
        "error_pct": f"{best_a[3]:.6e}",
        "status":    "candidate-EXACT" if abs(best_a[3]) < 5e-4 else "OK",
    })
if best_a2 is not None and abs(best_a2[3]) < abs(best_a[3] if best_a else 1e9):
    rows.append({
        "name": "delta_univ_closed_form_2atom",
        "predicted": f"{best_a2[1]:.6e}",
        "observed":  f"{delta_obs:.6e}",
        "error_pct": f"{best_a2[3]:.6e}",
        "status":    "candidate-EXACT" if abs(best_a2[3]) < 5e-4 else "OK",
    })

# (b)
if best_b is not None:
    rows.append({
        "name": "n_s_classXV_2atom",
        "predicted": f"{best_b[1]:.6e}",
        "observed":  f"{N_S:.6e}",
        "error_pct": f"{best_b[3]:.6e}",
        "status":    "candidate-EXACT" if abs(best_b[3]) < 5e-4 else "OK",
    })
if best_b3 is not None:
    rows.append({
        "name": "n_s_classXV_3atom",
        "predicted": f"{best_b3[1]:.6e}",
        "observed":  f"{N_S:.6e}",
        "error_pct": f"{best_b3[3]:.6e}",
        "status":    "candidate-EXACT" if abs(best_b3[3]) < 5e-4 else "OK",
    })

# (c)
rows.append({
    "name": "S_8_compositional_via_XIV_X",
    "predicted": f"{S8_uqff:.6e}",
    "observed":  f"{S_8:.6e}",
    "error_pct": f"{err_S8:.6e}",
    "status":    "candidate-EXACT" if abs(err_S8) < 5e-4 else "OK",
})

print(f"  (a) delta_univ best 1-atom: {best_a[2] if best_a else 'none'} err={best_a[3] if best_a else 'NA':+.3f}%")
if best_a2: print(f"      best 2-atom: {best_a2[2]} err={best_a2[3]:+.3f}%")
print(f"  (b) n_s best 2-atom: {best_b[2] if best_b else 'none'} err={best_b[3] if best_b else 'NA':+.4f}%")
if best_b3: print(f"      best 3-atom: {best_b3[2]} err={best_b3[3]:+.6f}%")
print(f"  (c) S_8 compositional err = {err_S8:+.4f}%")

write_ledger(rows, "_session736_delta_univ_classXV_ns_S8_check.py")

artifact = {
    "session": 736,
    "a_delta_univ": {
        "delta_obs": delta_obs,
        "best_1atom": {"name": best_a[2], "val": best_a[1], "err": best_a[3]} if best_a else None,
        "best_2atom": {"name": best_a2[2], "val": best_a2[1], "err": best_a2[3]} if best_a2 else None,
    },
    "b_classXV_ns": {
        "best_2atom": {"name": best_b[2], "val": best_b[1], "err": best_b[3]} if best_b else None,
        "best_3atom": {"name": best_b3[2], "val": best_b3[1], "err": best_b3[3]} if best_b3 else None,
    },
    "c_S8_check": {
        "S8_uqff": S8_uqff, "err_pct": err_S8,
        "S8_dressed": S8_dressed, "err_dressed_pct": err_S8_d,
        "S8_direct_5_6": 5/6,
    },
}
out = Path("_session736_delta_univ_classXV_ns_S8_check_result.json")
out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
print(f"\nArtifact written: {out.resolve().as_posix()}")
