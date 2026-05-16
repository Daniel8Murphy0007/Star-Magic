"""S365: Matter clustering amplitude sigma_8 = (1 + Phi_res - F_TRZ*K_Mex)/2."""
F_TRZ, Phi_res, K_Mex = 0.1, 5/6, 25/12
sigma_8 = (1 + Phi_res - F_TRZ*K_Mex)/2
obs = 0.811
err = 100*abs(sigma_8-obs)/obs
print(f"S365 COMPLETE. sigma_8 = (1 + Phi_res - F_TRZ*K_Mex)/2 = {sigma_8:.4f}; obs = {obs}; match {err:.3f}%.")
