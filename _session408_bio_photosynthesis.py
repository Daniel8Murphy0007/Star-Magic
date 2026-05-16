"""S408 Photosynthesis quantum requirement 8 photons/O2 -> yield 1/8 = 0.125"""
F_TRZ=0.1; K_Mex=25/12; D_phys=4
val = F_TRZ + F_TRZ*F_TRZ*K_Mex + F_TRZ*F_TRZ*F_TRZ*D_phys
obs = 0.125
print(f"S408 COMPLETE. PS quantum yield = F_TRZ+F_TRZ^2*K_Mex+F_TRZ^3*D_phys = {val:.4f}; obs (1/8) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
