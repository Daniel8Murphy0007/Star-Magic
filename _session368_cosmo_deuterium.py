"""S368: Deuterium D/H = F_TRZ^5 * (K_Mex + Phi_res * F_TRZ * D_BSFG)."""
F_TRZ, Phi_res, K_Mex, D_BSFG = 0.1, 5/6, 25/12, 6
DH = F_TRZ**5 * (K_Mex + Phi_res * F_TRZ * D_BSFG)
obs = 2.547e-5
err = 100*abs(DH-obs)/obs
print(f"S368 COMPLETE. D/H = F_TRZ^5 * (K_Mex + Phi_res*F_TRZ*D_BSFG) = {DH:.3e}; obs = {obs:.3e}; match {err:.2f}%.")
