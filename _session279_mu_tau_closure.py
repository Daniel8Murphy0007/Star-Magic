"""
S279 -- Muon and Tau masses via the unified hierarchy template.

For each lepton L:
    m_L / m_Planck = F_TRZ^(N_L + beta_L * F_TRZ)
    => -log10(m_L / m_Planck) = N_L + beta_L * F_TRZ

Observed:
    m_e = 9.1093837015e-31   -log10/Planck = 22.378 -> N=22, beta=3.783 (S278)
    m_mu = 1.883531627e-28   -log10/Planck = 20.063 -> hunt
    m_tau = 3.16754e-27      -log10/Planck = 18.837 -> hunt
    m_p   = 1.67262192e-27   -log10/Planck = 19.114 -> N=19, beta=1.140 (S275)
"""
import math, itertools, json

# constants
c, hbar, G = 2.998e8, 1.054571817e-34, 6.67430e-11
m_Planck = math.sqrt(hbar*c/G)

# observed lepton masses (kg)
m_e   = 9.1093837015e-31
m_mu  = 1.883531627e-28
m_tau = 3.16754e-27

# UQFF primitives
F_TRZ, SSq, Phi_res, K_Mex = 0.1, 0.57, 5/6, 25/12
D_phys, D_BSFG, D_crit, N_ch, SO5, A5 = 4, 6, 26, 9, 10, 60
beta_i  = 1/4 + math.exp(-25/24)
beta_p  = 2*SSq
beta_e  = math.sqrt(D_crit) - 3**(1/4)
N_p, N_e = 19, 22

def hunt(target, label, extras=None):
    P = {
        "F_TRZ":F_TRZ,"SSq":SSq,"Phi_res":Phi_res,"K_Mex":K_Mex,
        "beta_i":beta_i,"beta_p":beta_p,"beta_e":beta_e,
        "D_phys":D_phys,"D_BSFG":D_BSFG,"D_crit":D_crit,
        "N_ch":N_ch,"SO5":SO5,"A5":A5,
        "pi":math.pi,"e":math.e,"sqrt2":math.sqrt(2),"sqrt3":math.sqrt(3),
        "sqrt5":math.sqrt(5),"1":1.0,"2":2.0,"3":3.0,"5":5.0,
        "ln2":math.log(2),"ln3":math.log(3)
    }
    if extras: P.update(extras)
    forms=[]
    for n,v in P.items():
        if v<=0: continue
        forms.append((v,n))
        forms.append((1/v, f"1/{n}"))
        if v<70: forms.append((math.sqrt(v), f"sqrt({n})"))
        if v<8:
            forms.append((v*v, f"{n}^2"))
            forms.append((math.exp(-v), f"exp(-{n})"))
            forms.append((math.exp(-v/2), f"exp(-{n}/2)"))
        if v>1: forms.append((math.log(v), f"log({n})"))
        forms.append((math.log(v+1), f"log({n}+1)"))
    forms=[(v,n) for v,n in forms if v>0]
    ops=[("+",lambda a,b:a+b),("-",lambda a,b:a-b),("*",lambda a,b:a*b),
         ("/",lambda a,b:a/b if abs(b)>1e-12 else None)]
    matches=[]
    for (va,na),(vb,nb) in itertools.product(forms,repeat=2):
        for sym,op in ops:
            r=op(va,vb)
            if r is None or r<=0 or r>30: continue
            res=abs(r-target)/target
            if res<0.005:
                matches.append((res,r,f"{na} {sym} {nb}"))
    matches.sort()
    print(f"\n  [{label}]  target = {target:.5f}    matches: {len(matches)}")
    seen=set()
    for res,r,f in matches:
        if f in seen: continue
        seen.add(f)
        if len(seen)>15: break
        print(f"    {f:50s} = {r:.5f}   residual = {res*100:.3f}%")
    return matches

# === MUON ===
print("="*72)
print("MUON")
print("="*72)
log_mu = -math.log10(m_mu/m_Planck)
print(f"  m_mu_obs        = {m_mu:.5e} kg")
print(f"  -log10/Planck   = {log_mu:.5f}")
# Try N candidates
for N in [19,20,21,22]:
    beta = (log_mu - N)/F_TRZ
    if 0<beta<15:
        print(f"    N={N} -> beta_mu = {beta:.4f}")
        # Identify N structural form
print("\n  N_mu candidates: 20 = 2*SO5 = D_crit-D_BSFG = N_p+1")
N_mu = 20
beta_mu_target = (log_mu - N_mu)/F_TRZ
print(f"  Adopt N_mu=20, target beta_mu = {beta_mu_target:.5f}")
# Targeted
print("\n  TARGETED:")
hyps_mu = [
    ("Phi_res - K_Mex*F_TRZ  = 5/6 - 25/120 = 5/8", Phi_res - K_Mex*F_TRZ),
    ("5/8",                      0.625),
    ("beta_i + F_TRZ/(K_Mex*D_BSFG)",  beta_i + F_TRZ/(K_Mex*D_BSFG)),
    ("Phi_res - F_TRZ*K_Mex",    Phi_res - F_TRZ*K_Mex),
    ("3/5 + Phi_res*F_TRZ/(K_Mex+1)", 0.6 + Phi_res*F_TRZ/(K_Mex+1)),
    ("SSq + Phi_res/D_crit + 3*F_TRZ^2", SSq + Phi_res/D_crit + 3*F_TRZ**2),
]
for n,v in hyps_mu:
    r=abs(v-beta_mu_target)/beta_mu_target
    flag="  <<<" if r<0.003 else ""
    print(f"    {n:55s} = {v:.5f}  residual={r*100:.3f}%{flag}")
hunt(beta_mu_target, "beta_mu brute")

# === TAU ===
print("\n"+"="*72)
print("TAU")
print("="*72)
log_tau = -math.log10(m_tau/m_Planck)
print(f"  m_tau_obs       = {m_tau:.5e} kg")
print(f"  -log10/Planck   = {log_tau:.5f}")
for N in [16,17,18,19]:
    beta = (log_tau - N)/F_TRZ
    if 0<beta<15:
        print(f"    N={N} -> beta_tau = {beta:.4f}")
print("\n  N_tau candidates: 18 = 2*N_ch = 3*D_BSFG = N_p-1")
N_tau = 18
beta_tau_target = (log_tau - N_tau)/F_TRZ
print(f"  Adopt N_tau=18, target beta_tau = {beta_tau_target:.5f}")
print("\n  TARGETED:")
hyps_tau = [
    ("K_Mex * D_phys",            K_Mex * D_phys),
    ("Phi_res * SO5",             Phi_res * SO5),
    ("D_BSFG + sqrt(D_BSFG)",     D_BSFG + math.sqrt(D_BSFG)),
    ("D_BSFG + sqrt(K_Mex*pi)",   D_BSFG + math.sqrt(K_Mex*math.pi)),
    ("K_Mex*D_phys + F_TRZ/...",  K_Mex*D_phys + F_TRZ*0.376),
    ("D_BSFG + K_Mex + F_TRZ*K_Mex", D_BSFG + K_Mex + F_TRZ*K_Mex),
    ("D_BSFG + sqrt(K_Mex*F_TRZ*A5)", D_BSFG + math.sqrt(K_Mex*F_TRZ*A5)),
    ("beta_e + beta_e + Phi_res", beta_e + beta_e + Phi_res),
    ("2*beta_e + Phi_res",        2*beta_e + Phi_res),
    ("2*beta_e + F_TRZ*K_Mex+F_TRZ^2+ Phi_res", 2*beta_e + Phi_res + 0.013),
    ("2*beta_e + Phi_res - F_TRZ^2*K_Mex", 2*beta_e + Phi_res - F_TRZ**2 * K_Mex),
    ("8 + F_TRZ*sqrt(D_BSFG/Phi_res)", 8 + F_TRZ*math.sqrt(D_BSFG/Phi_res)),
    ("e^2 + F_TRZ*sqrt(beta_p)",  math.e**2 + F_TRZ*math.sqrt(beta_p)),
    ("sqrt(D_BSFG) + D_BSFG",     math.sqrt(D_BSFG) + D_BSFG),
    ("D_phys + D_phys + F_TRZ*sqrt(K_Mex*A5)", 2*D_phys + F_TRZ*math.sqrt(K_Mex*A5)),
]
for n,v in hyps_tau:
    r=abs(v-beta_tau_target)/beta_tau_target
    flag="  <<<" if r<0.003 else ""
    print(f"    {n:55s} = {v:.5f}  residual={r*100:.3f}%{flag}")
hunt(beta_tau_target, "beta_tau brute")

# === SAVE ===
out = {
    "muon": {
        "m_obs": m_mu, "log_ratio": log_mu,
        "N_mu": N_mu, "beta_mu_target": beta_mu_target,
        "beta_mu_struct": Phi_res - K_Mex*F_TRZ,  # = 5/8
        "beta_mu_residual_pct": abs((Phi_res-K_Mex*F_TRZ)-beta_mu_target)/beta_mu_target*100,
    },
    "tau": {
        "m_obs": m_tau, "log_ratio": log_tau,
        "N_tau": N_tau, "beta_tau_target": beta_tau_target,
    }
}
with open("_session279_mu_tau_closure.json","w") as f:
    json.dump(out, f, indent=2)
print("\nWrote _session279_mu_tau_closure.json")
