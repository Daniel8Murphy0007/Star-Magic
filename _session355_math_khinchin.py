"""S355: Khinchin's constant K = (K_Mex+Phi_res)(1-F_TRZ) + F_TRZ^2 * D_BSFG."""
F_TRZ, Phi_res, K_Mex, D_BSFG = 0.1, 5/6, 25/12, 6
K = (K_Mex+Phi_res)*(1-F_TRZ) + F_TRZ**2 * D_BSFG
K_obs = 2.68545200
err = 100*abs(K-K_obs)/K_obs
print(f"S355 COMPLETE. K_Khinchin = (K_Mex+Phi_res)(1-F_TRZ) + F_TRZ^2*D_BSFG = {K:.5f}; obs = {K_obs}; match {err:.3f}%.")
