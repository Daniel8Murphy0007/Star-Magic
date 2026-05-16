"""S392: Solar Schwarzschild radius / Solar radius ratio R_s/R_sun."""
F_TRZ=0.1; D_phys=4; SSq=0.57
val = F_TRZ**6 * D_phys * (1 + F_TRZ*SSq)
obs = 4.24e-6
print(f"S392 COMPLETE. R_s/R_sun = F_TRZ^6*D_phys*(1+F_TRZ*SSq) = {val:.3e}; obs (2.95 km / 6.96e5 km) = {obs:.2e}; match {abs(val-obs)/obs*100:.3f}%.")
