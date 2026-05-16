"""S371: CMB temperature T_CMB = Phi_res * (D_phys - Phi_res + F_TRZ) Kelvin."""
F_TRZ, Phi_res, D_phys = 0.1, 5/6, 4
T = Phi_res * (D_phys - Phi_res + F_TRZ)
obs = 2.7255
err = 100*abs(T-obs)/obs
print(f"S371 COMPLETE. T_CMB = Phi_res*(D_phys-Phi_res+F_TRZ) = {T:.4f} K; obs = {obs} K; match {err:.3f}%.")
