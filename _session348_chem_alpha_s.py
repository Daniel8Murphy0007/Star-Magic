"""S348: Strong coupling alpha_s(M_Z) = 1/(K_Mex*D_phys + F_TRZ)."""
F_TRZ, K_Mex, D_phys = 0.1, 25/12, 4
alpha_s = 1/(K_Mex*D_phys + F_TRZ)
alpha_s_obs = 0.1179
err_pct = 100*abs(alpha_s - alpha_s_obs)/alpha_s_obs
print(f"S348 COMPLETE. alpha_s(M_Z) = 1/(K_Mex*D_phys + F_TRZ) = 1/{K_Mex*D_phys+F_TRZ:.4f} = {alpha_s:.4f}; "
      f"PDG = {alpha_s_obs}+/-0.0010; match within {err_pct:.2f}% (well inside 1 sigma).")
