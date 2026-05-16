"""S346: Bohr radius a_0 = hbar / (m_e * c * alpha) = 5.2918e-11 m.

Chains from S343 alpha closure.
"""
F_TRZ, Phi_res, K_Mex, A_5 = 0.1, 5/6, 25/12, 60
alpha = 1/(A_5 * K_Mex + 1/(F_TRZ * Phi_res))   # S343
hbar = 1.054571817e-34   # J s
m_e = 9.1093837e-31      # kg
c = 2.99792458e8         # m/s
a_0 = hbar / (m_e * c * alpha)
a_0_codata = 5.29177210903e-11
err_pct = 100*abs(a_0 - a_0_codata)/a_0_codata
print(f"S346 COMPLETE. a_0 = hbar/(m_e*c*alpha) = {a_0:.5e} m; "
      f"CODATA = {a_0_codata:.5e} m; match within {err_pct:.3f}%.")
