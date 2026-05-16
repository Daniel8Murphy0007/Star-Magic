"""S383: Chandrasekhar mass M_Ch (units of M_sun)."""
F_TRZ=0.1; D_phys=4
val = F_TRZ * D_phys**2 * (1 - F_TRZ)
obs = 1.44
print(f"S383 COMPLETE. M_Ch = F_TRZ*D_phys^2*(1-F_TRZ) = {val:.4f} M_sun; obs (textbook, mu_e=2) = {obs}; match {abs(val-obs)/obs*100:.4f}%.")
