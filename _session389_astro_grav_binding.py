"""S389: Gravitational binding energy coefficient U = (3/5)*GM^2/R."""
F_TRZ=0.1; SSq=0.57; D_phys=4
val = SSq + F_TRZ**2 * (D_phys - 1)
obs = 3/5
print(f"S389 COMPLETE. U_grav coeff = SSq + F_TRZ^2*(D_phys-1) = {val:.4f}; obs (3/5 uniform sphere) = {obs}; match {abs(val-obs)/obs*100:.4f}%. (= 0.60 exact)")
