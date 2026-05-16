"""
S283: Deep fine-tune pass on the four worst residuals.

Targets (current best residual, mass space):
  m_tau:  0.186%   (was: 2*beta_e + Phi_res - K_Mex*F_TRZ^2)
  m_t  :  0.130%   (was: v_EW/sqrt(2))            -- different route
  m_u  :  0.107%   (was: K_Mex^2*sqrt(3))
  m_c  :  0.106%   (was: Phi_res - 1)

Method: expand token table, allow TRIPLE combinations (a op1 b op2 c),
restrict to UQFF-primitive-only forms, choose tightest sub-0.05% matches.
"""
import math, itertools, json

F_TRZ      = 0.1
SSq        = 0.57
Phi_res    = 5/6
K_Mex      = 25/12
D_phys, D_BSFG, D_crit = 4, 6, 26
N_ch, SO5, A5 = 9, 10, 60
beta_i  = 1/D_phys + math.exp(-K_Mex/2)
beta_p  = 2*SSq
beta_e  = math.sqrt(D_crit) - 3**0.25
beta_mu = beta_i**2 + 1/beta_e
beta_u  = K_Mex**2 * math.sqrt(3)
beta_b  = K_Mex * math.sqrt(5)

m_Planck_kg = 2.176434e-8
eV_to_kg    = 1.78266192e-36

# Targets
TARGETS = {
    'tau' : dict(m_obs=1.77686e9*eV_to_kg, N_list=[18, 19],
                 prev_resid=0.186),
    't'   : dict(m_obs=172.69e9*eV_to_kg, N_list=[16, 17, 18],
                 prev_resid=0.130),
    'u'   : dict(m_obs=2.16e6*eV_to_kg,    N_list=[20, 21, 22],
                 prev_resid=0.107),
    'c'   : dict(m_obs=1.27e9*eV_to_kg,    N_list=[18, 19, 20],
                 prev_resid=0.106),
}

# Restrict to ONLY pure UQFF primitives + low integers (no transcendentals
# unless they appeared elsewhere structurally)
TOKENS = {
    # primitives
    'F_TRZ': F_TRZ, 'SSq': SSq, 'Phi_res': Phi_res, 'K_Mex': K_Mex,
    'D_phys': D_phys, 'D_BSFG': D_BSFG, 'D_crit': D_crit,
    'N_ch': N_ch, 'SO5': SO5, 'A5': A5,
    # reciprocals & roots
    '1/Phi_res': 1/Phi_res, '1/K_Mex': 1/K_Mex, '1/D_BSFG': 1/D_BSFG,
    '1/D_phys': 1/D_phys, '1/SSq': 1/SSq,
    'sqrtD_phys': D_phys**0.5, 'sqrtD_BSFG': D_BSFG**0.5,
    'sqrtD_crit': D_crit**0.5, 'sqrtN_ch': N_ch**0.5, 'sqrtSO5': SO5**0.5,
    'sqrtA5': A5**0.5, 'sqrtK_Mex': K_Mex**0.5,
    # squares
    'D_phys^2': 16, 'D_BSFG^2': 36, 'K_Mex^2': K_Mex**2, 'Phi_res^2': Phi_res**2,
    'SSq^2': SSq**2, 'F_TRZ^2': F_TRZ**2,
    # transcendentals (used in prior closures: pi^2 in nu_3, e in beta_i)
    'pi': math.pi, 'pi^2': math.pi**2, '2pi': 2*math.pi,
    'e': math.e, 'e^2': math.e**2,
    # primitive roots (used in S278, S279)
    'sqrt2': 2**0.5, 'sqrt3': 3**0.5, 'sqrt5': 5**0.5,
    '3^(1/4)': 3**0.25, '5^(1/4)': 5**0.25,
    # derived betas (cross-sector links)
    'beta_i': beta_i, 'beta_p': beta_p, 'beta_e': beta_e,
    'beta_mu': beta_mu, 'beta_u': beta_u, 'beta_b': beta_b,
    '1/beta_e': 1/beta_e, '1/beta_i': 1/beta_i,
    'beta_e^2': beta_e**2, 'beta_i^2': beta_i**2,
    # small ints
    '0': 0.0, '1': 1.0, '2': 2.0, '3': 3.0, '4': 4.0,
    '5': 5.0, '6': 6.0, '7': 7.0, '8': 8.0,
}

def hunt_triple(name, m_obs, N_list, prev_resid):
    total = -math.log10(m_obs / m_Planck_kg)
    print(f"\n--- {name}  log10(m_Planck/m)={total:.4f}  prev_best={prev_resid}% ---")
    items = list(TOKENS.items())
    for N in N_list:
        target_beta = (total - N) / F_TRZ
        scored = []
        # singles
        for n,v in items:
            scored.append((n, v))
        # binaries
        for (n1,v1),(n2,v2) in itertools.product(items, repeat=2):
            for op,fn in [('+',lambda a,b:a+b),('-',lambda a,b:a-b),
                          ('*',lambda a,b:a*b)]:
                try:
                    scored.append((f"{n1}{op}{n2}", fn(v1,v2)))
                except: pass
        # triples a + b*c  (additive correction to product)
        prim_small = [(k,v) for k,v in items if abs(v) < 50]
        for (n1,v1) in prim_small:
            for (n2,v2),(n3,v3) in itertools.product(prim_small, repeat=2):
                try:
                    val_p = v2*v3
                    if abs(val_p) > 1000: continue
                    scored.append((f"{n1}+{n2}*{n3}", v1 + val_p))
                    scored.append((f"{n1}-{n2}*{n3}", v1 - val_p))
                except: pass

        ranked = []
        for nm, val in scored:
            try:
                m_pred = m_Planck_kg * 10**(-(N + val*F_TRZ))
                mr = abs(m_pred - m_obs)/m_obs*100
                ranked.append((mr, nm, val, m_pred))
            except: pass
        ranked.sort()
        print(f"  N={N}  target_beta={target_beta:.5f}")
        seen_vals = set()
        shown = 0
        for mr, nm, val, m_pred in ranked:
            key = round(val, 5)
            if key in seen_vals: continue
            seen_vals.add(key)
            if mr > prev_resid * 1.5: break
            print(f"     beta={val:>10.5f}  '{nm:<32}'  mass_resid={mr:>7.4f}%")
            shown += 1
            if shown >= 8: break

for name, t in TARGETS.items():
    hunt_triple(name, t['m_obs'], t['N_list'], t['prev_resid'])
