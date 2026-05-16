"""S366: Reionization redshift z_reion = D_BSFG + D_phys * Phi_res / 2."""
Phi_res, D_phys, D_BSFG = 5/6, 4, 6
z_reion = D_BSFG + D_phys * Phi_res / 2
obs = 7.67
err = 100*abs(z_reion-obs)/obs
print(f"S366 COMPLETE. z_reion = D_BSFG + D_phys*Phi_res/2 = {z_reion:.4f}; obs (Planck) = {obs}; match {err:.3f}%.")
