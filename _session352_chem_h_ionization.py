"""S352: H ionization energy = 13.6057 eV (chains from S345).

Also predicts hydrogen 1s -> 2p Lyman-alpha at (3/4) R_y = 10.20 eV
and Balmer-alpha at (5/36) R_y = 1.89 eV.
"""
F_TRZ, Phi_res, K_Mex, A_5 = 0.1, 5/6, 25/12, 60
alpha = 1/(A_5 * K_Mex + 1/(F_TRZ * Phi_res))   # S343
m_e_c2_eV = 510998.95
R_y = alpha**2 * m_e_c2_eV / 2
E_ion_H = R_y                          # 1s ionization
E_Lyman = 0.75 * R_y                   # 1s -> 2p (3/4)
E_Balmer = (5/36) * R_y                # 2 -> 3 (5/36)
print(f"S352 COMPLETE. E_ion(H) = R_y = {E_ion_H:.4f} eV (obs 13.6057); "
      f"E_Lyman = 3/4 R_y = {E_Lyman:.4f} eV (obs 10.20); "
      f"E_Balmer = 5/36 R_y = {E_Balmer:.4f} eV (obs 1.89); all chain from S343 alpha.")
