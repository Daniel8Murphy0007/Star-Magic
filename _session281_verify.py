"""S281 verify: canonical clean neutrino mass closures."""
import math

F_TRZ      = 0.1
K_Mex      = 25/12
D_phys, D_BSFG, D_crit = 4, 6, 26
m_Planck_kg = 2.176434e-8
eV_to_kg    = 1.78266192e-36

# canonical structural forms
def predict(N, beta):
    return m_Planck_kg * 10**(-(N + beta*F_TRZ))

# beta_u from S280 (up quark)
beta_u = K_Mex**2 * math.sqrt(3)
print(f"beta_u (S280, up quark) = K_Mex^2 * sqrt(3) = {beta_u:.6f}")

# nu_2  -  link to up quark!
beta_nu2 = D_phys + beta_u  # = D_phys + K_Mex^2*sqrt(3)
m_nu2 = predict(29, beta_nu2)
print(f"\nm_nu_2:  N=29, beta = D_phys + beta_u = {beta_nu2:.4f}")
print(f"   pred = {m_nu2/eV_to_kg*1e3:.4f} meV   obs = 8.610 meV  resid = {abs(m_nu2/eV_to_kg*1e3-8.610)/8.610*100:.3f}%")

# nu_3  -  pure UQFF primitive
beta_nu3 = math.pi**2 - D_BSFG
m_nu3 = predict(29, beta_nu3)
print(f"\nm_nu_3:  N=29, beta = pi^2 - D_BSFG = {beta_nu3:.4f}")
print(f"   pred = {m_nu3/eV_to_kg*1e3:.4f} meV   obs = 50.10 meV  resid = {abs(m_nu3/eV_to_kg*1e3-50.10)/50.10*100:.3f}%")

# Sigma m_nu cosmological
beta_sum = K_Mex + 1  # = 25/12 + 1 = 37/12
m_sum = predict(29, beta_sum)
print(f"\nSigma m_nu:  N=29, beta = K_Mex + 1 = 37/12 = {beta_sum:.4f}")
print(f"   pred = {m_sum/eV_to_kg*1e3:.4f} meV   obs ~ 60 meV  resid = {abs(m_sum/eV_to_kg*1e3-60)/60*100:.3f}%")

print("\n--- mass-squared differences (oscillation observables) ---")
dm21_pred = (m_nu2/eV_to_kg)**2  # m_1 ~ 0
dm31_pred = (m_nu3/eV_to_kg)**2
print(f"dm21^2:  pred = {dm21_pred*1e5:.3f}e-5 eV^2   obs = 7.420e-5 eV^2   resid = {abs(dm21_pred*1e5-7.420)/7.420*100:.3f}%")
print(f"dm31^2:  pred = {dm31_pred*1e3:.3f}e-3 eV^2   obs = 2.510e-3 eV^2   resid = {abs(dm31_pred*1e3-2.510)/2.510*100:.3f}%")

# mass ratio - this is parameter-free
ratio_pred = m_nu3/m_nu2
ratio_obs  = math.sqrt(2.51e-3/7.42e-5)
print(f"\nm_3/m_2:  pred = {ratio_pred:.4f}   obs = {ratio_obs:.4f}   resid = {abs(ratio_pred-ratio_obs)/ratio_obs*100:.3f}%")
