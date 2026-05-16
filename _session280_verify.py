"""S280 verification: canonical clean quark closures."""
import math

F_TRZ      = 0.1
SSq        = 0.57
Phi_res    = 5/6
K_Mex      = 25/12
D_phys, D_BSFG, D_crit = 4, 6, 26
N_ch, SO5 = 9, 10
m_Planck_kg = 2.176434e-8
eV_to_kg    = 1.78266192e-36

# canonical forms (chosen for UQFF-primitive cleanliness)
QUARKS = {
    'u': dict(m_obs=2.16e6*eV_to_kg, N=21,
              beta_form='K_Mex^2 * sqrt(3)',
              beta_val = K_Mex**2 * math.sqrt(3)),
    'd': dict(m_obs=4.67e6*eV_to_kg, N=21,
              beta_form='2*K_Mex + F_TRZ^2',
              beta_val = 2*K_Mex + F_TRZ**2),
    's': dict(m_obs=93.4e6*eV_to_kg, N=20,
              beta_form='sqrt(SO5) - sqrt(D_phys)',
              beta_val = math.sqrt(SO5) - math.sqrt(D_phys)),
    'c': dict(m_obs=1.27e9*eV_to_kg, N=19,
              beta_form='Phi_res - 1  =  -1/D_BSFG',
              beta_val = Phi_res - 1),
    'b': dict(m_obs=4.18e9*eV_to_kg, N=18,
              beta_form='K_Mex * sqrt(5)',
              beta_val = K_Mex * math.sqrt(5)),
}

print(f"{'q':>3} {'N':>3} {'beta':>10}  {'form':<28}  {'m_pred[kg]':>14}  {'m_obs[kg]':>14}  resid")
print("="*110)
for q, d in QUARKS.items():
    m_pred = m_Planck_kg * 10**(-(d['N'] + d['beta_val']*F_TRZ))
    resid  = abs(m_pred - d['m_obs'])/d['m_obs']*100
    print(f"{q:>3} {d['N']:>3} {d['beta_val']:>10.4f}  {d['beta_form']:<28}  "
          f"{m_pred:>14.4e}  {d['m_obs']:>14.4e}  {resid:>6.3f}%")

print("\n--- mass ratios (predicted vs observed) ---")
m_pred = {q: m_Planck_kg * 10**(-(d['N'] + d['beta_val']*F_TRZ)) for q, d in QUARKS.items()}
m_obs  = {q: d['m_obs'] for q, d in QUARKS.items()}
pairs = [('d','u'),('s','d'),('c','s'),('b','c')]
for a,b in pairs:
    pr = m_pred[a]/m_pred[b]; ob = m_obs[a]/m_obs[b]
    print(f"  m_{a}/m_{b}:  pred={pr:>9.4f}   obs={ob:>9.4f}   resid={abs(pr-ob)/ob*100:.3f}%")

# integer ladder summary
print("\n--- integer ladder (mass-ordered) ---")
print("  u : N=21 = D_crit - SO5/2     (= D_phys+D_BSFG+N_ch+2)")
print("  d : N=21    (same as u; isospin doublet)")
print("  s : N=20 = 2*SO5               (same as muon N_mu)")
print("  c : N=19 = N_p                 (same as proton N_p)")
print("  b : N=18 = 2*N_ch              (same as tau N_tau)")
print("  t : closed via v_EW/sqrt(2)   (S272)")
