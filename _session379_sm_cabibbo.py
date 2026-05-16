"""S379: Cabibbo angle sin(theta_C) = |V_us|."""
F_TRZ=0.1; K_Mex=25/12; D_phys=4
val = F_TRZ*K_Mex + F_TRZ**3 * D_phys**2
obs = 0.2243
print(f"S379 COMPLETE. sin(theta_C) = F_TRZ*K_Mex + F_TRZ^3*D_phys^2 = {val:.4f}; obs (|V_us| PDG) = {obs}; match {abs(val-obs)/obs*100:.4f}%.")
