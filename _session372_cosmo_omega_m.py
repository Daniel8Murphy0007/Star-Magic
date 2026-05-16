"""S372: Matter density Omega_m = SSq - K_Mex*F_TRZ - Phi_res*F_TRZ*SSq."""
F_TRZ, Phi_res, K_Mex, SSq = 0.1, 5/6, 25/12, 0.57
Om = SSq - K_Mex*F_TRZ - Phi_res*F_TRZ*SSq
obs = 0.315
err = 100*abs(Om-obs)/obs
print(f"S372 COMPLETE. Omega_m = SSq - K_Mex*F_TRZ - Phi_res*F_TRZ*SSq = {Om:.4f}; obs (Planck) = {obs}; match {err:.3f}%.")
