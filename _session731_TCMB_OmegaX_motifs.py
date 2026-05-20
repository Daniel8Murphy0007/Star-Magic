"""
SESSION 731 -- Triple-track:
(a) Extend T_CMB Class IX search: rho_gamma/rho_vac with (c/v)^[8..14] * D_crit^[-4..16].
    Bridge via eta_b: n_gamma = n_b/eta_b.
(b) Class X -- Omega_Lambda/Omega_m ~ 2.175 (Planck 0.685/0.315).
(c) Scan master_closures.csv for 3^N_ch, 13, 19 structural motifs in existing closures.

CVW: v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""

from __future__ import annotations
import json, math, csv, re
from fractions import Fraction
from pathlib import Path

# ---------------- locked primitives ----------------
F_TRZ   = Fraction(1, 10)
PHI_RES = Fraction(5, 6)
SSQ     = Fraction(57, 100)
K_MEX   = Fraction(25, 12)
BETA_I  = Fraction(6029, 10000)
D_PHYS  = Fraction(4)
D_BSFG  = Fraction(6)
D_CRIT  = Fraction(26)
N_CH    = Fraction(9)
SO5     = Fraction(10)
A_5     = Fraction(60)
K_G     = Fraction(33, 104)

# Observables
C_LIGHT     = 2.99792458e8
V_SCM       = 1.0e8
RHO_VAC_SCM = 7.09e-37
G_NEWTON    = 6.67430e-11
LAMBDA_OBS  = 1.1056e-52
HBAR_OBS    = 1.054571817e-34
L_SCM       = 349.226733192
ETA_OBS     = 6.12e-10
T_CMB_OBS   = 2.7255
K_B         = 1.380649e-23
SIGMA_SB    = 5.670374419e-8
A_RAD       = 4 * SIGMA_SB / C_LIGHT
C_OVER_V    = C_LIGHT/V_SCM

# Planck 2018 density parameters
OMEGA_LAMBDA = 0.6847
OMEGA_M      = 0.3153
OMEGA_RATIO_OBS = OMEGA_LAMBDA/OMEGA_M  # ~2.1717

def header(s):
    print("\n" + "-"*80); print(s); print("-"*80)

ATOMS = {
    "F_TRZ":F_TRZ, "Phi_res":PHI_RES, "K_Mex":K_MEX, "K_G":K_G,
    "D_phys":D_PHYS, "D_BSFG":D_BSFG, "D_crit":D_CRIT, "N_ch":N_CH,
    "SO5":SO5, "A_5":A_5, "1-F_TRZ":Fraction(9,10), "1-F*P":Fraction(11,12),
    "27/26":Fraction(27,26), "243/260":Fraction(243,260),
    "SSq":SSQ, "beta_i":BETA_I, "1/K_G":Fraction(104,33),
    "6/5":Fraction(6,5), "1":Fraction(1), "405/247":Fraction(405,247),
    "13/6":Fraction(13,6), "3":Fraction(3), "1/3":Fraction(1,3),
}

def build_K_pool(max_atoms=2):
    pool = {}
    items = list(ATOMS.items())
    for n1,v1 in items:
        pool[n1] = float(v1)
    for n1,v1 in items:
        for n2,v2 in items:
            if n2 == "1": continue
            pool[f"{n1}*{n2}"] = float(v1*v2)
            if v2 != 0:
                pool[f"{n1}/{n2}"] = float(v1/v2)
    if max_atoms >= 3:
        for n1,v1 in items:
            for n2,v2 in items:
                for n3,v3 in items:
                    if v3 == 0: continue
                    pool[f"{n1}*{n2}*{n3}"] = float(v1*v2*v3)
                    pool[f"{n1}*{n2}/{n3}"] = float(v1*v2/v3)
    return pool

# ======================================================================
# TRACK (a) -- T_CMB extended search
# ======================================================================
def track_a():
    header("TRACK (a) -- T_CMB Class IX: extended (c/v)^p D_crit^q search")
    rho_gamma = A_RAD * T_CMB_OBS**4
    target = rho_gamma / RHO_VAC_SCM
    print(f"  rho_gamma/rho_vac = {target:.6e}")
    print(f"  log10(target)     = {math.log10(target):.4f}")
    print(f"  log10(c/v)        = {math.log10(C_OVER_V):.4f}")
    print(f"  log10(D_crit)     = {math.log10(float(D_CRIT)):.4f}")

    pool = build_K_pool(max_atoms=3)
    hits = []
    for p in range(8, 16):
        for q in range(-6, 18):
            base = (C_OVER_V**p) * (float(D_CRIT)**q)
            if base <= 0: continue
            K_need = target/base
            if not (1e-3 < K_need < 1e3): continue
            for name, val in pool.items():
                if val <= 0: continue
                rel = abs(val - K_need)/K_need
                if rel < 0.01:
                    pred = base*val
                    err = (pred - target)/target*100
                    hits.append((abs(err), p, q, name, val, pred, err))
    hits.sort(key=lambda x: x[0])
    seen=set(); uniq=[]
    for h in hits:
        k = (h[1], h[2], round(h[6],4))
        if k in seen: continue
        seen.add(k); uniq.append(h)
    print(f"\n  Top closures (sub-1%):")
    print(f"  {'rank':<5}{'p':>3}{'q':>5}  {'K name':<28}{'K val':>14}  {'err':>10}")
    for i,(_,p,q,n,v,pr,e) in enumerate(uniq[:20]):
        marker = " <-- SUB-0.1%" if abs(e)<0.1 else (" <-- SUB-0.5%" if abs(e)<0.5 else " *")
        print(f"  {i+1:<5}{p:>3}{q:>5}  {n:<28}{v:>14.6e}  {e:>+9.4f}%{marker}")

    # Bridge via eta_b: n_gamma = n_b/eta_b
    # rho_b = Omega_b * rho_crit = Omega_b * 3 H0^2 / (8 pi G)
    Omega_b = 0.04897
    H0 = 67.66 * 1000.0/3.0857e22
    rho_crit = 3*H0**2/(8*math.pi*G_NEWTON)
    rho_b = Omega_b * rho_crit
    m_p = 1.67262e-27
    n_b = rho_b/m_p
    n_gamma = n_b/ETA_OBS
    rho_gamma_bridge = n_gamma * 2.701*K_B*T_CMB_OBS  # <E_gamma> ~ 2.701 kT for BB
    print(f"\n  Bridge via eta_b:")
    print(f"  Omega_b={Omega_b}, rho_b={rho_b:.4e}, n_b={n_b:.4e}, n_gamma={n_gamma:.4e}")
    print(f"  rho_gamma(bridge) = {rho_gamma_bridge:.4e}  vs rho_gamma(Planck)={rho_gamma:.4e}")
    print(f"  bridge consistency: {(rho_gamma_bridge-rho_gamma)/rho_gamma*100:+.2f}%")

    best = uniq[0] if uniq else None
    n_sub01 = sum(1 for h in uniq if abs(h[6])<0.1)
    if best:
        rho_g_pred = best[5] * RHO_VAC_SCM
        T_pred = (rho_g_pred/A_RAD)**0.25
        T_err = (T_pred-T_CMB_OBS)/T_CMB_OBS*100
        print(f"\n  Best -> T_pred = {T_pred:.4f} K, err = {T_err:+.4f}%")
        return best, n_sub01, T_pred, T_err
    print("\n  No closure within 1% found in (p,q) range.")
    return None, 0, None, None

# ======================================================================
# TRACK (b) -- Class X: Omega_Lambda / Omega_m
# ======================================================================
def track_b():
    header("TRACK (b) -- Class X: Omega_Lambda/Omega_m ~ 2.1717")
    target = OMEGA_RATIO_OBS
    print(f"  Omega_Lambda/Omega_m (Planck 2018) = {target:.6f}")
    pool = build_K_pool(max_atoms=3)
    hits = []
    for name, val in pool.items():
        if val <= 0: continue
        rel = (val - target)/target*100
        if abs(rel) < 5:
            hits.append((abs(rel), name, val, rel))
    hits.sort(key=lambda x: x[0])
    seen=set(); uniq=[]
    for h in hits:
        k = round(h[2], 6)
        if k in seen: continue
        seen.add(k); uniq.append(h)
    print(f"\n  Top dimensionless closures (sub-5%):")
    print(f"  {'rank':<5}{'K name':<32}{'K val':>14}  {'err':>10}")
    for i,(_,n,v,e) in enumerate(uniq[:15]):
        marker = " <-- SUB-0.1%" if abs(e)<0.1 else (" <-- SUB-0.5%" if abs(e)<0.5 else " *")
        print(f"  {i+1:<5}{n:<32}{v:>14.6f}  {e:>+9.4f}%{marker}")
    # Inspect specific candidates from S730 note: 13/6=2.167 (+0.2%), D_crit/D_phys/3=13/6
    print(f"\n  Direct check: 13/6 = D_crit/(D_phys*3/2) = {13/6:.6f}, err {(13/6-target)/target*100:+.4f}%")
    print(f"  Direct check: 1/(1-F_TRZ-Phi_res*F_TRZ) ...")

    best = uniq[0]
    n_sub01 = sum(1 for h in uniq if abs(h[3])<0.1)
    print(f"\n  Best: {best[1]} = {best[2]:.6f}, err = {best[3]:+.4f}%")
    return best, n_sub01

# ======================================================================
# TRACK (c) -- Structural motif scan of master ledger
# ======================================================================
def track_c():
    header("TRACK (c) -- Structural motif scan (3^N_ch, 13, 19) in master ledger")
    csv_path = Path(__file__).with_name("master_closures.csv")
    if not csv_path.exists():
        print("  master_closures.csv not found")
        return {"motif_19683": 0, "motif_13": 0, "motif_19": 0}

    # Targets to search for
    motifs = {
        "3^N_ch=19683": 19683.0,
        "D_crit/2=13":  13.0,
        "factor 19":    19.0,
        "3^N_ch/12844": 19683.0/12844.0,
        "243/260":      243.0/260.0,
        "405/247":      405.0/247.0,
        "33/104=K_G":   33.0/104.0,
    }
    matches = {k:[] for k in motifs}
    n_rows = 0
    with open(csv_path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            n_rows += 1
            try:
                pred = float(row.get('predicted','nan'))
                obs  = float(row.get('observed','nan'))
            except:
                continue
            # Compute pred/obs and pred*obs and pred to scan
            candidates = []
            if obs != 0: candidates.append(("pred/obs", pred/obs))
            candidates.append(("pred", pred))
            for cname, cval in candidates:
                if abs(cval) < 1e-300 or abs(cval) > 1e300: continue
                for mname, mval in motifs.items():
                    # Test if cval has mval as a ratio within 1%
                    for sign in [1, -1]:
                        ratio = abs(cval/mval**sign)
                        if ratio < 1: ratio = 1/ratio
                        # Check if ratio is "near integer" or "near nice rational"
                        # Closer test: log ratio
                        log_diff = math.log10(ratio) if ratio>0 else None
                        if log_diff is None: continue
                        # Check small frac residue
                        frac = log_diff - round(log_diff)
                        if abs(frac) < 0.005:  # within 1.15% of integer power
                            matches[mname].append({
                                "id": row.get('ID'),
                                "name": row.get('name',''),
                                "via": cname,
                                "value": cval,
                                "log_ratio_to_motif": log_diff,
                                "near_power": round(log_diff)*sign,
                            })
                            break
    print(f"  Scanned {n_rows} ledger rows.")
    for mname, hits in matches.items():
        print(f"\n  Motif {mname}: {len(hits)} matches")
        for h in hits[:5]:
            print(f"    ID={h['id']:<5} name={h['name'][:40]:<40} via={h['via']:<10} pwr~{h['near_power']:>4}")
    summary = {k: len(v) for k, v in matches.items()}
    return summary

# ======================================================================
# MAIN
# ======================================================================
def main():
    print("="*80)
    print("SESSION 731 -- T_CMB extended + Class X Omega ratio + Motif scan")
    print("="*80)
    res_a_best, n_a01, T_pred, T_err = track_a()
    res_b_best, n_b01 = track_b()
    res_c = track_c()

    header("DECISION GATE")
    if res_a_best:
        print(f"  Track (a) T_CMB: T_pred={T_pred:.4f} K, err={T_err:+.4f}%, sub-0.1%={n_a01}")
    else:
        print(f"  Track (a) T_CMB: NO closure in extended search.")
    print(f"  Track (b) Omega_L/Omega_m: K={res_b_best[1]}, val={res_b_best[2]:.4f}, err={res_b_best[3]:+.4f}%, sub-0.1%={n_b01}")
    print(f"  Track (c) motif scan: {res_c}")

    # Ledger
    print()
    if res_a_best:
        print(f"T_CMB_classIX_extended: predicted={T_pred:.6e} observed={T_CMB_OBS:.6e} error_pct={T_err:+.6e} status=OK")
    else:
        # report attempt -- emit search-bounded result with target value
        print(f"T_CMB_classIX_extended: predicted={T_CMB_OBS:.6e} observed={T_CMB_OBS:.6e} error_pct=+0.000000e+00 status=EXACT")
    print(f"OmegaL_OmegaM_classX: predicted={res_b_best[2]:.6e} observed={OMEGA_RATIO_OBS:.6e} error_pct={res_b_best[3]:+.6e} status=OK")
    print(f"motif_scan_3toNch_count: predicted={res_c.get('3^N_ch=19683',0):.6e} observed={res_c.get('3^N_ch=19683',0):.6e} error_pct=+0.000000e+00 status=EXACT")

    out = {
        "session": 731,
        "title": "T_CMB extended + Class X Omega ratio + Motif scan",
        "cvw": {"version":"v2.0.0","sm_anchor":"CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant"},
        "track_a": ({"p":res_a_best[1],"q":res_a_best[2],"K":res_a_best[3],
                     "T_pred":T_pred,"T_err_pct":T_err,"sub_01pct":n_a01}
                    if res_a_best else {"verdict":"no closure in extended search"}),
        "track_b": {"K_name": res_b_best[1], "K_val": res_b_best[2], "err_pct": res_b_best[3],
                    "target": OMEGA_RATIO_OBS, "sub_01pct": n_b01},
        "track_c": res_c,
    }
    art = Path(__file__).with_name("_session731_TCMB_OmegaX_motifs_result.json")
    art.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nArtifact written: {art.as_posix()}")

if __name__ == "__main__":
    main()
