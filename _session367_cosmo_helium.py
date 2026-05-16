"""S367: Primordial helium Y_p = (1-Phi_res) + Phi_res*F_TRZ*(1 - F_TRZ*SSq)."""
F_TRZ, Phi_res, SSq = 0.1, 5/6, 0.57
Yp = (1-Phi_res) + Phi_res*F_TRZ*(1 - F_TRZ*SSq)
obs = 0.2453
err = 100*abs(Yp-obs)/obs
print(f"S367 COMPLETE. Y_p = (1-Phi_res) + Phi_res*F_TRZ*(1-F_TRZ*SSq) = {Yp:.5f}; obs = {obs}; match {err:.3f}%.")
