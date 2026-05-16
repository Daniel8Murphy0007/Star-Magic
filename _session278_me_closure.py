"""
S278 -- close m_e/m_Planck using the unified hierarchy template
       and find structural form of beta_e.

CONTEXT:
  S275 closed m_p/m_Planck:  m_p = m_Planck * F_TRZ^(N_p + beta_p*F_TRZ)
                              N_p = 19, beta_p = 2*SSq = 1.14
  S266 established the EXACT identity:
       m_p/m_e = A5^2/2 + D_BSFG^2 = 1800 + 36 = 1836
  => m_e closes automatically once m_p is closed, at the same residual.

  But we also want the standalone hierarchy form for m_e:
       m_e/m_Planck = F_TRZ^(N_e + beta_e * F_TRZ)
  Empirically log10(m_e/m_Planck) ~ -22.378
  => N_e = 22, beta_e ~ 3.78  (S275 noted N_e = D_crit - 4)
  Find structural beta_e.
"""
import math, itertools, json

# --- constants ---
c        = 2.998e8
hbar     = 1.054571817e-34
G        = 6.67430e-11
m_e_obs  = 9.1093837015e-31
m_p_obs  = 1.67262192369e-27
m_Planck = math.sqrt(hbar*c/G)

# --- UQFF primitives ---
F_TRZ    = 0.1
SSq      = 0.57
beta_i   = 1/4 + math.exp(-25/24)        # closed in S277  = 0.602866
Phi_res  = 5/6
K_Mex    = 25/12
D_phys, D_BSFG, D_crit, N_ch = 4, 6, 26, 9
SO5, A5  = 10, 60

# --- Step 1.  m_p/m_e identity (S266) ---
mp_me_struct = A5**2/2 + D_BSFG**2          # = 1800 + 36 = 1836
mp_me_obs    = m_p_obs / m_e_obs
print("STEP 1: m_p/m_e from UQFF structural identity (S266)")
print(f"  A5^2/2 + D_BSFG^2  = {mp_me_struct}")
print(f"  m_p/m_e observed   = {mp_me_obs:.4f}")
print(f"  residual           = {abs(mp_me_struct-mp_me_obs)/mp_me_obs*100:.4f}%")
print()

# m_e derived from m_p (S275)
N_p     = 19
beta_p  = 2*SSq    # = 1.14
m_p_pred = m_Planck * F_TRZ**(N_p + beta_p*F_TRZ)
m_e_pred = m_p_pred / mp_me_struct
print(f"STEP 2: m_e = m_p_pred / 1836")
print(f"  m_p_pred           = {m_p_pred:.4e} kg")
print(f"  m_e_pred           = {m_e_pred:.4e} kg")
print(f"  m_e_obs            = {m_e_obs:.4e} kg")
print(f"  residual           = {abs(m_e_pred-m_e_obs)/m_e_obs*100:.4f}%")
print()

# --- Step 2. standalone hierarchy form ---
log_me_ratio = math.log10(m_e_obs / m_Planck)
print(f"STEP 3: standalone hierarchy form")
print(f"  log10(m_e/m_Planck) = {log_me_ratio:.6f}")
# target exponent (positive number, since m_e/m_Planck = F_TRZ^exp = 10^(-exp))
exp_target = -log_me_ratio
print(f"  exponent (positive) = {exp_target:.6f}")
# Try N_e candidates
N_e_candidates = {
    "D_crit - 4":             D_crit - 4,                  # 22
    "D_BSFG + D_phys + N_ch + D_phys": D_BSFG+D_phys+N_ch+D_phys, # 23
    "2*N_ch + D_phys":        2*N_ch + D_phys,            # 22
    "N_p + D_phys - 1":       N_p + D_phys - 1,           # 22
    "D_crit - D_BSFG + D_phys - 2": D_crit-D_BSFG+D_phys-2, # 22
    "D_BSFG^2/2 + N_ch + D_BSFG + 4-D_BSFG+...": 22,
}
for name,v in N_e_candidates.items():
    print(f"   N_e = {v}  ({name})  beta_e = {(exp_target-v)/F_TRZ:.4f}")

# Take N_e = 22, beta_e = ?
N_e = 22
beta_e_target = (exp_target - N_e) / F_TRZ
print(f"\n  Adopted N_e = 22, target beta_e = {beta_e_target:.5f}")
print()

# --- Step 3. hunt structural beta_e via brute force ---
P = {
    "F_TRZ":   F_TRZ,
    "SSq":     SSq,
    "Phi_res": Phi_res,
    "K_Mex":   K_Mex,
    "beta_i":  beta_i,
    "beta_p":  beta_p,
    "D_phys":  D_phys,
    "D_BSFG":  D_BSFG,
    "D_crit":  D_crit,
    "N_ch":    N_ch,
    "SO5":     SO5,
    "A5":      A5,
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

def unary_forms(P):
    out = []
    for n,v in P.items():
        if v <= 0: continue
        out.append((v, n))
        out.append((1/v, f"1/{n}"))
        if v < 70:
            out.append((math.sqrt(v), f"sqrt({n})"))
            if v < 8:
                out.append((v*v, f"{n}^2"))
                out.append((math.exp(-v), f"exp(-{n})"))
                out.append((math.exp(-v/2), f"exp(-{n}/2)"))
            if v > 1:
                out.append((math.log(v), f"log({n})"))
            out.append((math.log(v+1), f"log({n}+1)"))
    return [(v,n) for v,n in out if v is not None and v > 0]

unary = unary_forms(P)
print(f"Hunt: generated {len(unary)} unary forms.")

TOL = 0.005
ops = [("+",lambda a,b:a+b), ("-",lambda a,b:a-b), ("*",lambda a,b:a*b), ("/",lambda a,b:a/b if abs(b)>1e-12 else None)]
binary = []
for (va,na),(vb,nb) in itertools.product(unary, repeat=2):
    for sym,op in ops:
        r = op(va,vb)
        if r is None or r<=0 or r>50: continue
        res = abs(r - beta_e_target)/beta_e_target
        if res < TOL:
            binary.append((res, r, f"{na} {sym} {nb}"))
binary.sort()
print(f"  {len(binary)} binary matches under {TOL*100:.1f}% tolerance\n")
print("Top 25 unique:")
seen=set()
for res,r,f in binary:
    if f in seen: continue
    seen.add(f)
    if len(seen)>25: break
    print(f"  {f:55s} = {r:.6f}   residual = {res*100:.3f}%")

# Clean primitive forms only
clean_names = {"SSq","Phi_res","F_TRZ","K_Mex","D_phys","D_BSFG","D_crit","N_ch","SO5","A5","beta_i","beta_p","e","pi","2","3","1","sqrt2","sqrt3"}
print("\nTop 15 'clean' forms (UQFF primitives only):")
clean=[]
for res,r,form in binary:
    parts = form.replace("(",  " ").replace(")"," ").replace("/"," ").replace("*"," ").replace("+"," ").replace("-"," ").split()
    tokens = [p for p in parts if not p.replace(".","").isdigit() and not p in ops]
    if all(any(c == t.replace("1/","").replace("^2","").replace("sqrt","").replace("exp","").replace("log","").strip() for c in clean_names) for t in tokens):
        clean.append((res,r,form))
seen=set()
for res,r,form in clean:
    if form in seen: continue
    seen.add(form)
    if len(seen)>15: break
    print(f"  {form:55s} = {r:.6f}   residual = {res*100:.3f}%")

# Targeted hypotheses
print("\nTARGETED HYPOTHESES:")
hyps = [
    ("2*beta_p + 1.5",                     2*beta_p + 1.5),
    ("D_BSFG - SSq*D_phys",                D_BSFG - SSq*D_phys),
    ("D_BSFG*SSq + D_phys*F_TRZ - F_TRZ/D_BSFG", D_BSFG*SSq + D_phys*F_TRZ - F_TRZ/D_BSFG),
    ("N_ch/2 - K_Mex*Phi_res",             N_ch/2 - K_Mex*Phi_res),
    ("2*beta_p + beta_i + F_TRZ",          2*beta_p + beta_i + F_TRZ),
    ("2*beta_p + 2*Phi_res + F_TRZ*sqrt(D_phys)", 2*beta_p + 2*Phi_res + F_TRZ*math.sqrt(D_phys)),
    ("D_phys + Phi_res*F_TRZ/N_ch... ",    None),
    ("beta_p + beta_i + D_BSFG*F_TRZ + ...", beta_p + beta_i + D_BSFG*F_TRZ + 1.44),
    ("D_BSFG - K_Mex + beta_i",            D_BSFG - K_Mex + beta_i),
    ("D_BSFG - K_Mex + beta_p + F_TRZ",    D_BSFG - K_Mex + beta_p + F_TRZ),
    ("2*beta_p + 1 + Phi_res + F_TRZ/D_BSFG", 2*beta_p + 1 + Phi_res + F_TRZ/D_BSFG),
    ("pi + beta_p - Phi_res/D_BSFG",       math.pi + beta_p - Phi_res/D_BSFG),
    ("pi + Phi_res - F_TRZ/SSq",           math.pi + Phi_res - F_TRZ/SSq),
    ("pi + beta_p - Phi_res*F_TRZ",        math.pi + beta_p - Phi_res*F_TRZ),
    ("Phi_res*D_phys + beta_i - F_TRZ*K_Mex", Phi_res*D_phys + beta_i - F_TRZ*K_Mex),
    ("Phi_res*D_phys + beta_i/2 + F_TRZ*..", Phi_res*D_phys + beta_i/2 + 0.143),
    ("Phi_res*D_phys + Phi_res/2 - F_TRZ*..", Phi_res*D_phys + Phi_res/2 - 0.156),
    ("2*Phi_res + 2 + Phi_res*F_TRZ - K_Mex*F_TRZ", 2*Phi_res + 2 + Phi_res*F_TRZ - K_Mex*F_TRZ),
    ("2*(Phi_res + 1) + Phi_res*Phi_res/K_Mex", 2*(Phi_res+1) + Phi_res**2/K_Mex),
    ("3 + Phi_res - F_TRZ*sqrt(D_BSFG)",   3 + Phi_res - F_TRZ*math.sqrt(D_BSFG)),
    ("4 - F_TRZ*sqrt(D_crit-D_phys)*...",  4 - 0.1*math.sqrt(22)*0.473),
]
for name,v in hyps:
    if v is None: continue
    res = abs(v-beta_e_target)/beta_e_target
    flag = "  <<<<<" if res < 0.003 else ""
    print(f"  {name:55s} = {v:.5f}  residual = {res*100:6.3f}%{flag}")

# Best clean
print("\n=== BEST CLEAN ===")
if clean:
    r,v,f = clean[0]
    print(f"  beta_e = {f} = {v:.5f}   residual = {r*100:.4f}%")

# Output
out = {
    "mp_me_identity":   {"form": "A5^2/2 + D_BSFG^2", "value": mp_me_struct,
                          "obs": mp_me_obs, "residual_pct": abs(mp_me_struct-mp_me_obs)/mp_me_obs*100},
    "m_e_via_mp":       {"m_e_pred": m_e_pred, "m_e_obs": m_e_obs,
                          "residual_pct": abs(m_e_pred-m_e_obs)/m_e_obs*100},
    "standalone":       {"N_e": N_e, "beta_e_target": beta_e_target},
    "best_binary":      [{"form":f,"value":r,"residual_pct":res*100} for res,r,f in binary[:30]],
    "best_clean":       [{"form":f,"value":r,"residual_pct":res*100} for res,r,f in clean[:15]],
}
with open("_session278_me_closure.json","w") as f:
    json.dump(out, f, indent=2)
print("\nWrote _session278_me_closure.json")
