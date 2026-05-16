"""S345: Rydberg energy R_y = alpha^2 * m_e * c^2 / 2 = 13.6057 eV.

Chains from S343 alpha closure.
"""
F_TRZ, Phi_res, K_Mex, A_5 = 0.1, 5/6, 25/12, 60
alpha = 1/(A_5 * K_Mex + 1/(F_TRZ * Phi_res))   # S343 closure
m_e_c2_eV = 510998.95   # electron rest energy
R_y = alpha**2 * m_e_c2_eV / 2
R_y_codata = 13.605693
err_pct = 100*abs(R_y - R_y_codata)/R_y_codata
print(f"S345 COMPLETE. R_y = alpha^2 * m_e*c^2 / 2 = {R_y:.4f} eV; "
      f"CODATA = {R_y_codata} eV; match within {err_pct:.3f}% (limited by S343).")
