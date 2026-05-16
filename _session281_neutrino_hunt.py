"""
S281: Neutrino mass scale closure.

Oscillation observables (PDG 2024, normal ordering, m_1 ~ 0 limit):
  dm21^2 = 7.42e-5 eV^2  ->  m_2 ~ sqrt(dm21^2) = 8.61 meV
  dm31^2 = 2.51e-3 eV^2  ->  m_3 ~ sqrt(dm31^2) = 50.10 meV
  Sigma_m_nu cosmological best fit (Planck+DESI) ~ 0.06 eV  ->  consistent

Apply unified hierarchy:
  m_nu / m_Planck = F_TRZ^(N + beta*F_TRZ)
  log10(m_Planck/m_nu) = N + beta*F_TRZ
"""
import math, json, itertools

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
beta_tau = 2*beta_e + Phi_res - K_Mex*F_TRZ**2

m_Planck_kg = 2.176434e-8
eV_to_kg    = 1.78266192e-36

TARGETS = {
    'nu_2': dict(m_eV=8.61e-3,
                 description='m_2 = sqrt(dm21^2) = sqrt(7.42e-5) eV'),
    'nu_3': dict(m_eV=5.010e-2,
                 description='m_3 = sqrt(dm31^2) = sqrt(2.51e-3) eV'),
    'sum':  dict(m_eV=0.060,
                 description='Sigma m_nu cosmological best fit'),
}

TOKENS = {
    'F_TRZ': F_TRZ, 'SSq': SSq, 'Phi_res': Phi_res, 'K_Mex': K_Mex,
    'D_phys': D_phys, 'D_BSFG': D_BSFG, 'D_crit': D_crit,
    'N_ch': N_ch, 'SO5': SO5, 'A5': A5,
    '1/Phi_res': 1/Phi_res, '1/K_Mex': 1/K_Mex,
    'sqrtD_phys': D_phys**0.5, 'sqrtD_BSFG': D_BSFG**0.5, 'sqrtD_crit': D_crit**0.5,
    'sqrtN_ch': N_ch**0.5, 'sqrtSO5': SO5**0.5, 'sqrtA5': A5**0.5,
    'D_phys^2': 16, 'D_BSFG^2': 36, 'K_Mex^2': K_Mex**2,
    'pi': math.pi, 'pi^2': math.pi**2, 'e': math.e, 'e^2': math.e**2,
    '2pi': 2*math.pi, 'sqrt2': 2**0.5, 'sqrt3': 3**0.5, 'sqrt5': 5**0.5,
    '3^(1/4)': 3**0.25, '5^(1/4)': 5**0.25,
    'beta_i': beta_i, 'beta_p': beta_p, 'beta_e': beta_e,
    'beta_mu': beta_mu, 'beta_tau': beta_tau,
    '1/beta_e': 1/beta_e,  '1/beta_tau': 1/beta_tau,
    'K_Mex^2*sqrt3': K_Mex**2 * math.sqrt(3),
    '0': 0.0, '1': 1.0, '2': 2.0, '3': 3.0, '4': 4.0, '5': 5.0, '6': 6.0,
    '7': 7.0, '8': 8.0, '9': 9.0, '10': 10.0,
}

def hunt(name, m_eV, N_list):
    m_kg = m_eV * eV_to_kg
    total = -math.log10(m_kg / m_Planck_kg)
    print(f"\n--- {name}  m={m_eV*1e3:.3f} meV  log10(m_Planck/m)={total:.4f} ---")
    items = list(TOKENS.items())
    for N in N_list:
        target_beta = (total - N) / F_TRZ
        scored = []
        for (n1,v1) in items:
            scored.append((abs(v1-target_beta)/abs(target_beta) if target_beta else abs(v1),
                           n1, v1))
        for (n1,v1),(n2,v2) in itertools.product(items, repeat=2):
            for op,fn in [('+',lambda a,b:a+b),('-',lambda a,b:a-b),('*',lambda a,b:a*b)]:
                try:
                    val = fn(v1,v2)
                    r = abs(val-target_beta)/abs(target_beta) if target_beta else abs(val)
                    scored.append((r, f"{n1}{op}{n2}", val))
                except Exception: pass
        scored.sort()
        print(f"  N={N}  target_beta={target_beta:.4f}")
        for r,nm,v in scored[:5]:
            m_pred = m_Planck_kg * 10**(-(N + v*F_TRZ))
            mr = abs(m_pred - m_kg)/m_kg*100
            print(f"     beta={v:>9.4f}  '{nm:<28}'  mass_resid={mr:>7.3f}%")

hunt('nu_2', TARGETS['nu_2']['m_eV'], [29, 30, 31])
hunt('nu_3', TARGETS['nu_3']['m_eV'], [28, 29, 30])
hunt('sum',  TARGETS['sum']['m_eV'],  [28, 29, 30])
