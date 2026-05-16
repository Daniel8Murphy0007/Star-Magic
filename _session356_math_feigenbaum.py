"""S356: Feigenbaum constant delta = D_phys + Phi_res*(1 - F_TRZ*K_Mex)."""
F_TRZ, Phi_res, K_Mex, D_phys = 0.1, 5/6, 25/12, 4
delta = D_phys + Phi_res*(1 - F_TRZ*K_Mex)
delta_obs = 4.66920160
err = 100*abs(delta-delta_obs)/delta_obs
print(f"S356 COMPLETE. delta_Feigenbaum = D_phys + Phi_res*(1 - F_TRZ*K_Mex) = {delta:.5f}; obs = {delta_obs}; match {err:.3f}%.")
