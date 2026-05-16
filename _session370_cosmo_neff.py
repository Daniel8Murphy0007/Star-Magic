"""S370: Effective neutrino species N_eff = D_phys - Phi_res - F_TRZ*K_Mex*SSq."""
F_TRZ, Phi_res, K_Mex, SSq, D_phys = 0.1, 5/6, 25/12, 0.57, 4
N_eff = D_phys - Phi_res - F_TRZ * K_Mex * SSq
obs = 3.046
err = 100*abs(N_eff-obs)/obs
print(f"S370 COMPLETE. N_eff = D_phys - Phi_res - F_TRZ*K_Mex*SSq = {N_eff:.4f}; obs (SM) = {obs}; match {err:.3f}%.")
