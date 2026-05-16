"""S363: Baryon-to-photon ratio eta_B = D_BSFG * F_TRZ^10 = 6e-10."""
F_TRZ, D_BSFG = 0.1, 6
eta = D_BSFG * F_TRZ**10
eta_obs = 6.14e-10
err = 100*abs(eta-eta_obs)/eta_obs
print(f"S363 COMPLETE. eta_B = D_BSFG * F_TRZ^10 = {eta:.3e}; obs = {eta_obs:.3e}; match {err:.2f}%.")
