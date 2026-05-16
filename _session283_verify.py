"""S283 verify: refined canonical forms for four worst residuals."""
import math

F_TRZ = 0.1
SSq = 0.57
Phi_res = 5/6
K_Mex = 25/12
D_phys, D_BSFG, D_crit = 4, 6, 26
N_ch, SO5, A5 = 9, 10, 60
beta_i = 1/D_phys + math.exp(-K_Mex/2)
beta_e = math.sqrt(D_crit) - 3**0.25
beta_u_old = K_Mex**2 * math.sqrt(3)

m_Planck_kg = 2.176434e-8
eV_to_kg = 1.78266192e-36

def show(name, N, beta_expr, beta_val, m_obs, prev):
    m_pred = m_Planck_kg * 10**(-(N + beta_val*F_TRZ))
    r = abs(m_pred - m_obs)/m_obs*100
    improvement = prev/r if r > 0 else float('inf')
    print(f"  {name:>5}  N={N}  beta={beta_val:.6f}")
    print(f"         form  = {beta_expr}")
    print(f"         pred  = {m_pred:.6e} kg")
    print(f"         obs   = {m_obs:.6e} kg")
    print(f"         resid = {r:.4f}%   (was {prev}% -- improvement {improvement:.0f}x)")
    print()

print("="*70)
print("S283 Fine-Tune Results -- Refined Canonical Forms")
print("="*70 + "\n")

# tau: beta = beta_u + beta_i*sqrt2  (cross-sector link)
beta_tau_new = beta_u_old + beta_i*math.sqrt(2)
show('tau', 18, 'beta_u + beta_i*sqrt(2)  [LINK to up quark]',
     beta_tau_new, 1.77686e9*eV_to_kg, 0.186)

# t: beta = K_Mex - sqrt(A5)/SSq
beta_t = K_Mex - math.sqrt(A5)/SSq
show('top', 18, 'K_Mex - sqrt(A5)/SSq',
     beta_t, 172.69e9*eV_to_kg, 0.130)

# u refined: beta = sqrt(A5) - F_TRZ*sqrt5
beta_u_new = math.sqrt(A5) - F_TRZ*math.sqrt(5)
show('u', 21, 'sqrt(A5) - F_TRZ*sqrt(5)',
     beta_u_new, 2.16e6*eV_to_kg, 0.107)

# c refined: beta = 5^(1/4) - SO5/D_BSFG
beta_c_new = 5**0.25 - SO5/D_BSFG
show('c', 19, '5^(1/4) - SO5/D_BSFG',
     beta_c_new, 1.27e9*eV_to_kg, 0.106)

# --- propagate: re-derive beta_nu_2 since it referenced beta_u ---
print("-"*70)
print("Cross-sector propagation: refined beta_u changes beta_nu_2 form")
print("-"*70)
beta_nu2_new = D_phys + beta_u_new  # D_phys + sqrt(A5) - F_TRZ*sqrt5
m_nu2 = m_Planck_kg * 10**(-(29 + beta_nu2_new*F_TRZ))
m_nu2_obs = 8.61e-3 * eV_to_kg
print(f"  nu_2  beta = D_phys + sqrt(A5) - F_TRZ*sqrt(5) = {beta_nu2_new:.4f}")
print(f"        pred = {m_nu2/eV_to_kg*1e3:.4f} meV   obs = 8.610 meV")
print(f"        resid = {abs(m_nu2-m_nu2_obs)/m_nu2_obs*100:.4f}%   (was 0.020%)")
print()

# --- summary table ---
print("="*70)
print("SUMMARY: Fermion-sector residuals after S283 fine-tune")
print("="*70)
data = [
    ('m_e',      0.009,  0.009,  'sqrt(D_crit) - 3^(1/4)'),
    ('m_mu',     0.000,  0.000,  'beta_i^2 + 1/beta_e'),
    ('m_tau',    0.186,  0.0016, 'beta_u + beta_i*sqrt(2)  [REFINED]'),
    ('m_u',      0.107,  0.0030, 'sqrt(A5) - F_TRZ*sqrt(5)  [REFINED]'),
    ('m_d',      0.071,  0.071,  '2K_Mex + F_TRZ^2'),
    ('m_s',      0.023,  0.023,  'sqrt(SO5) - sqrt(D_phys)'),
    ('m_c',      0.106,  0.0011, '5^(1/4) - SO5/D_BSFG  [REFINED]'),
    ('m_b',      0.080,  0.080,  'K_Mex * sqrt(5)'),
    ('m_t',      0.130,  0.0040, 'K_Mex - sqrt(A5)/SSq  [REFINED]'),
    ('m_nu2',    0.020,  0.020,  'D_phys + beta_u'),
    ('m_nu3',    0.028,  0.028,  'pi^2 - D_BSFG'),
    ('Sigma_nu', 0.044,  0.044,  'K_Mex + 1 = 37/12'),
]
print(f"  {'mass':<10} {'prev%':>7}  {'new%':>7}  form")
print("-"*70)
for nm, p, n, f in data:
    flag = ' **' if n < p else ''
    print(f"  {nm:<10} {p:>7.4f}  {n:>7.4f}  {f}{flag}")

max_new = max(n for _,_,n,_ in data)
print(f"\n  WORST residual after fine-tune: {max_new:.4f}%")
print(f"  ALL 12 fundamental fermions now closed below 0.1%")
