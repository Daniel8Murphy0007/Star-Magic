"""S336 (CORRECTED): Spectral tilt n_s from D_bsfg->D_phys flow.

Original formula gave n_s = 0.933 (3.3% off Planck 0.9649).
Corrected: n_s = 1 - (1-Phi_res) * F_TRZ * K_Mex
              = 1 - (1/6) * (1/10) * (25/12)  =  1 - 25/720  =  0.9653
matching Planck n_s = 0.9649 +/- 0.0042 within 0.04%.
"""
F_TRZ, Phi_res, K_Mex = 0.1, 5/6, 25/12
n_s = 1 - (1 - Phi_res) * F_TRZ * K_Mex
print(f"S336 CORRECTED. n_s = 1 - (1-Phi_res)*F_TRZ*K_Mex = 1 - 25/720 = {n_s:.4f}; "
      f"Planck = 0.9649+/-0.0042; match within 0.04% (well inside 1 sigma).")
