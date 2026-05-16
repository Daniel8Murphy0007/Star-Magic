"""
S280: Quark mass closure via unified hierarchy template.

Target: m_q = m_Planck * F_TRZ^(N_q + beta_q * F_TRZ)
     => N_q + beta_q*F_TRZ = log10(m_Planck/m_q)

Observed quark masses (PDG 2024):
  m_u ~ 2.16 MeV (MSbar @ 2 GeV)
  m_d ~ 4.67 MeV
  m_s ~ 93.4 MeV
  m_c ~ 1.27 GeV (MSbar)
  m_b ~ 4.18 GeV (MSbar)
  m_t ~ 172.69 GeV  [already closed S272 via v_EW]

Strategy: same brute-force-then-clean approach used for leptons.
Hypothesis: integer ladder descends in steps; gen 1 -> N=21, gen 2 -> N=20/19,
gen 3 -> N=18 (matches lepton tower pattern).
"""

import math, json, itertools, sys

# --- locked UQFF primitives ---
F_TRZ      = 0.1
SSq        = 0.57
Phi_res    = 5/6
K_Mex      = 25/12
D_phys     = 4
D_BSFG     = 6
D_crit     = 26
N_ch       = 9
SO5        = 10
A5         = 60
# derived betas locked from prior sessions
beta_i     = 1/D_phys + math.exp(-K_Mex/2)            # S277
beta_p     = 2*SSq                                     # S275
beta_e     = math.sqrt(D_crit) - 3**0.25               # S278
beta_mu    = beta_i**2 + 1/beta_e                      # S279
beta_tau   = 2*beta_e + Phi_res - K_Mex*F_TRZ**2       # S279

# --- constants ---
m_Planck_kg = 2.176434e-8
eV_to_kg    = 1.78266192e-36

QUARKS = {
    'u': 2.16e6  * eV_to_kg,
    'd': 4.67e6  * eV_to_kg,
    's': 93.4e6  * eV_to_kg,
    'c': 1.27e9  * eV_to_kg,
    'b': 4.18e9  * eV_to_kg,
}

# --- primitive token table ---
TOKENS = {
    'F_TRZ': F_TRZ, 'SSq': SSq, 'Phi_res': Phi_res, 'K_Mex': K_Mex,
    'D_phys': D_phys, 'D_BSFG': D_BSFG, 'D_crit': D_crit,
    'N_ch': N_ch, 'SO5': SO5, 'A5': A5,
    '1/Phi_res': 1/Phi_res, '1/K_Mex': 1/K_Mex,
    'sqrtD_phys': D_phys**0.5, 'sqrtD_BSFG': D_BSFG**0.5, 'sqrtD_crit': D_crit**0.5,
    'sqrtN_ch': N_ch**0.5, 'sqrtSO5': SO5**0.5,
    'D_phys^2': 16, 'D_BSFG^2': 36, 'K_Mex^2': K_Mex**2,
    'pi': math.pi, 'pi^2': math.pi**2, 'e': math.e, 'e^2': math.e**2,
    '2pi': 2*math.pi, 'sqrt2': 2**0.5, 'sqrt3': 3**0.5, 'sqrt5': 5**0.5,
    '3^(1/4)': 3**0.25, '5^(1/4)': 5**0.25, '2^(1/4)': 2**0.25,
    'beta_i': beta_i, 'beta_p': beta_p, 'beta_e': beta_e,
    'beta_mu': beta_mu, 'beta_tau': beta_tau,
    '1/beta_i': 1/beta_i, '1/beta_e': 1/beta_e, '1/beta_p': 1/beta_p,
    'beta_e^2': beta_e**2, 'beta_i^2': beta_i**2,
    'K_Mex*D_BSFG': K_Mex*D_BSFG, 'SSq*D_BSFG': SSq*D_BSFG,
    '2K_Mex': 2*K_Mex, '2SSq': 2*SSq, '2Phi_res': 2*Phi_res,
    'F_TRZ^2': F_TRZ**2, '10F_TRZ': 1.0,
    'D_crit-D_BSFG': 20, 'D_crit-D_phys': 22, 'D_crit-N_ch': 17,
    'log_pi': math.log(math.pi), 'ln10': math.log(10),
    '0': 0.0, '1': 1.0, '2': 2.0, '3': 3.0, '4': 4.0, '5': 5.0, '6': 6.0,
}

# integer ladder candidates per quark
N_CANDIDATES = {
    'u': [20, 21, 22],
    'd': [20, 21, 22],
    's': [19, 20, 21],
    'c': [18, 19, 20],
    'b': [17, 18, 19],
}


def hunt(name, m_obs):
    target_total = -math.log10(m_obs / m_Planck_kg)  # = N + beta*F_TRZ
    results = []
    for N in N_CANDIDATES[name]:
        target_beta = (target_total - N) / F_TRZ
        # search single tokens and binary combos
        candidates = []
        for tname, tval in TOKENS.items():
            candidates.append((tname, tval))
        items = list(TOKENS.items())
        for (n1, v1), (n2, v2) in itertools.product(items, repeat=2):
            for op, fn in [('+', lambda a,b: a+b),
                           ('-', lambda a,b: a-b),
                           ('*', lambda a,b: a*b)]:
                try:
                    val = fn(v1, v2)
                    candidates.append((f"{n1}{op}{n2}", val))
                except Exception:
                    pass
        # rank by residual
        scored = []
        for cname, cval in candidates:
            if cval is None: continue
            try:
                resid = abs(cval - target_beta) / abs(target_beta) if target_beta != 0 else abs(cval)
                scored.append((resid, cname, cval))
            except Exception:
                pass
        scored.sort()
        for resid, cname, cval in scored[:6]:
            mass_resid = abs(10**(-(N + cval*F_TRZ)) - m_obs/m_Planck_kg) / (m_obs/m_Planck_kg)
            results.append({
                'N': N, 'target_beta': target_beta,
                'form': cname, 'value': cval, 'beta_resid_pct': resid*100,
                'mass_resid_pct': mass_resid*100,
            })
    return target_total, results


report = {}
print(f"{'q':>3}  {'log10(M/m)':>10}  best closures")
print("="*90)
for q, m in QUARKS.items():
    total, res = hunt(q, m)
    print(f"\n--- {q}-quark  log10(m_Planck/m_{q}) = {total:.4f} ---")
    # show top 8 across all N
    res.sort(key=lambda r: r['mass_resid_pct'])
    for r in res[:8]:
        print(f"  N={r['N']:>2}  beta={r['value']:>9.4f}  form='{r['form']:<28}'  "
              f"mass_resid={r['mass_resid_pct']:>7.3f}%")
    report[q] = {'log10_ratio': total, 'top': res[:8]}

with open('_session280_quark_closure.json', 'w') as f:
    json.dump(report, f, indent=2, default=str)
print("\nWrote _session280_quark_closure.json")
