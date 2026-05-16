"""S375: Atmospheric neutrino mixing sin^2(theta_23)."""
F_TRZ=0.1; D_phys=4; SSq=0.57
val = SSq * (1 - F_TRZ**2 * D_phys)
obs = 0.55
print(f"S375 COMPLETE. sin^2(theta_23) = SSq*(1-F_TRZ^2*D_phys) = {val:.4f}; obs (NuFit 5.2) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
