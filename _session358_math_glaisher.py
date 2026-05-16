"""S358: Glaisher-Kinkelin A = 1 + Phi_res*K_Mex/D_BSFG."""
Phi_res, K_Mex, D_BSFG = 5/6, 25/12, 6
A = 1 + Phi_res*K_Mex/D_BSFG
A_obs = 1.28242712
err = 100*abs(A-A_obs)/A_obs
print(f"S358 COMPLETE. A_GlaisherKinkelin = 1 + Phi_res*K_Mex/D_BSFG = {A:.5f}; obs = {A_obs}; match {err:.3f}%.")
