"""
S277 -- structural hunt for beta_i = 0.603

The Lambda hierarchy (S274) requires the fractional correction in the exponent
to be beta_i * F_TRZ where beta_i = 0.6028 (solved value).
S270 had calibrated beta_i ~ 0.603 from buoyancy-force vacuum coupling.

Goal: find the deepest structural form of 0.6028 using only UQFF primitives.

Target: beta_i = 0.6028 (from S274 closure of Lambda)
"""
import math, itertools, json

# ----- target -----
beta_i_target = 0.6028     # solved from S274 (Lambda closure exponent)
TOL = 0.005                # accept residual < 0.5%

# ----- primitives -----
P = {
    "F_TRZ":   0.1,
    "SSq":     0.57,
    "Phi_res": 5/6,
    "K_Mex":   25/12,
    "D_phys":  4,
    "D_BSFG":  6,
    "D_crit":  26,
    "N_ch":    9,
    "SO5":     10,
    "A5":      60,
    "pi":      math.pi,
    "e":       math.e,
    "sqrt2":   math.sqrt(2),
    "sqrt3":   math.sqrt(3),
    "sqrt5":   math.sqrt(5),
    "1":       1.0,
    "2":       2.0,
    "3":       3.0,
    "ln2":     math.log(2),
    "ln3":     math.log(3),
}

# Generate all unary forms (value, name)
def unary_forms():
    out = []
    for n,v in P.items():
        if v <= 0: continue
        out.append((v, n))
        out.append((1/v, f"1/{n}"))
        if v < 10:
            out.append((math.sqrt(v), f"sqrt({n})"))
            out.append((v*v, f"{n}^2"))
            out.append((math.exp(-v), f"exp(-{n})"))
            out.append((math.exp(-v/2), f"exp(-{n}/2)"))
            out.append((math.log(v) if v>1 else None, f"log({n})"))
            out.append((math.log(v+1), f"log({n}+1)"))
    return [(v,n) for v,n in out if v is not None]

unary = unary_forms()
print(f"Generated {len(unary)} unary forms.")

# Find single-term matches
print("\n=== UNARY MATCHES (single primitive) ===")
single = []
for v,n in unary:
    if v <= 0 or v > 5: continue
    res = abs(v - beta_i_target)/beta_i_target
    if res < 0.02:
        single.append((res, v, n))
single.sort()
for res,v,n in single[:15]:
    print(f"  beta_i ~= {n:30s} = {v:.6f}   residual = {res*100:.3f}%")

# Binary combinations: a op b
print("\n=== BINARY MATCHES (a OP b) ===")
binary = []
ops = [
    ("+",  lambda a,b: a+b),
    ("-",  lambda a,b: a-b),
    ("*",  lambda a,b: a*b),
    ("/",  lambda a,b: a/b if abs(b)>1e-12 else None),
]
for (va,na),(vb,nb) in itertools.product(unary, repeat=2):
    if va<=0 or vb<=0: continue
    for sym,op in ops:
        try: r = op(va,vb)
        except: continue
        if r is None or r<=0 or r>5: continue
        res = abs(r - beta_i_target)/beta_i_target
        if res < TOL:
            binary.append((res, r, f"{na} {sym} {nb}"))
binary.sort()
seen = set()
print("Top 25 unique binary forms:")
for res,r,form in binary:
    if form in seen: continue
    seen.add(form)
    if len(seen) > 25: break
    print(f"  {form:50s} = {r:.6f}   residual = {res*100:.3f}%")

# Highlight especially clean small-integer forms
print("\n=== CLEAN BINARY (small integers / clean primitives) ===")
clean_names = {"SSq","Phi_res","F_TRZ","K_Mex","D_phys","D_BSFG","D_crit","N_ch","SO5","A5","e","pi","2","3","1"}
clean_binary = []
for res,r,form in binary:
    parts = form.replace("(",  " ").replace(")"," ").replace("/"," ").replace("*"," ").replace("+"," ").replace("-"," ").split()
    tokens = [p for p in parts if not p.replace(".","").isdigit()]
    if all(any(c in t for c in clean_names) or t.isdigit() for t in tokens):
        clean_binary.append((res,r,form))
clean_binary.sort()
seen=set()
for res,r,form in clean_binary:
    if form in seen: continue
    seen.add(form)
    if len(seen) > 12: break
    print(f"  {form:50s} = {r:.6f}   residual = {res*100:.3f}%")

# Targeted hypothesis tests:
print("\n=== TARGETED HYPOTHESES ===")
hyps = [
    ("exp(-1/2)",                  math.exp(-0.5)),
    ("SSq + Phi_res/D_crit",       0.57 + (5/6)/26),
    ("SSq + F_TRZ*(D_BSFG-2*pi*F_TRZ)", 0.57 + 0.1*(6 - 2*math.pi*0.1)),
    ("1 - SSq*Phi_res - F_TRZ",    1 - 0.57*(5/6) - 0.1),
    ("Phi_res - F_TRZ*(K_Mex+...)",0.833 - 0.1*(25/12 + 0.196)),
    ("SSq + (Phi_res-SSq)*F_TRZ",  0.57 + (5/6 - 0.57)*0.1),
    ("SSq + F_TRZ*Phi_res/(D_phys-1)*K_Mex", 0.57 + 0.1*(5/6)/3*(25/12)),
    ("(D_BSFG-SSq*D_phys)/(D_BSFG+SSq*D_BSFG)", (6-0.57*4)/(6+0.57*6)),
    ("(D_crit-N_ch-N_ch)/(D_crit+SO5-D_crit/D_BSFG*pi)", None),
    ("3/5 + F_TRZ*F_TRZ/3",        0.6 + 0.01/3),
    ("3/5 + Phi_res*F_TRZ/sqrt(5)",0.6 + (5/6)*0.1/math.sqrt(5)),
    ("3/5 + F_TRZ*Phi_res/...",    0.6 + 0.1*(5/6)*0.0335),
    ("(1+SSq)/(SSq+pi)",           (1+0.57)/(0.57+math.pi)),
    ("Phi_res*Phi_res*K_Mex/(D_phys/2)", (5/6)**2*(25/12)/2),
    ("K_Mex - 3/2 - F_TRZ/sqrt(D_BSFG)",   25/12 - 1.5 - 0.1/math.sqrt(6)),
    ("(D_BSFG+2*F_TRZ)/(K_Mex*5)", (6+0.2)/(25/12*5)),
    ("SSq + Phi_res*K_Mex*F_TRZ^2",0.57 + (5/6)*(25/12)*0.01),
    ("1 - K_Mex*Phi_res*F_TRZ*2",  1 - (25/12)*(5/6)*0.1*2),
    ("Phi_res*Phi_res*K_Mex*...",  None),
    ("1 - pi/N_ch + F_TRZ/...",    None),
    ("(D_crit-D_BSFG*sqrt(pi))/D_crit*...", None),
    ("Phi_res - sqrt(F_TRZ)/...",  None),
    ("3*Phi_res - K_Mex*F_TRZ - Phi_res", None),
    # SI-Lambda inspired:
    ("Phi_res*Phi_res - pi*F_TRZ/sqrt(D_crit)", (5/6)**2 - math.pi*0.1/math.sqrt(26)),
    ("Phi_res - F_TRZ*pi*Phi_res*Phi_res", (5/6) - 0.1*math.pi*(5/6)**2),
    ("3/5 + (pi-3)*F_TRZ/sqrt(2)", 0.6 + (math.pi-3)*0.1/math.sqrt(2)),
    ("3/5 + (pi-3)/(K_Mex*D_BSFG)",0.6 + (math.pi-3)/(25/12*6)),
]
for name, val in hyps:
    if val is None: continue
    res = abs(val - beta_i_target)/beta_i_target
    flag = "  <<<<" if res < 0.003 else ""
    print(f"  {name:55s} = {val:.6f}  residual = {res*100:6.3f}%{flag}")

# ------- Top result -------
all_results = [(r,v,n) for r,v,n in binary]
all_results.sort()
if all_results:
    best = all_results[0]
    print("\n=== BEST OVERALL ===")
    print(f"  beta_i ~= {best[2]} = {best[1]:.6f}   residual = {best[0]*100:.4f}%")

# Save
out = {
    "target_beta_i": beta_i_target,
    "tolerance":     TOL,
    "n_binary_matches": len(binary),
    "best_unary":    [{"form":n, "value":v, "residual_pct":res*100} for res,v,n in single[:10]],
    "best_binary":   [{"form":f, "value":r, "residual_pct":res*100} for res,r,f in binary[:30]],
    "clean_forms":   [{"form":f, "value":r, "residual_pct":res*100} for res,r,f in clean_binary[:15]],
}
with open("_session277_beta_i_hunt.json","w") as f:
    json.dump(out, f, indent=2)
print(f"\nWrote _session277_beta_i_hunt.json  ({len(binary)} binary candidates)")
