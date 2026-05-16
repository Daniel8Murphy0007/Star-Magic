"""S390: Salpeter initial mass function slope alpha."""
F_TRZ=0.1; K_Mex=25/12; Phi_res=5/6; D_BSFG=6; D_phys=4
val = K_Mex + Phi_res - F_TRZ*D_BSFG + F_TRZ**2 * (D_phys - Phi_res)
obs = 2.35
print(f"S390 COMPLETE. alpha_Salpeter = K_Mex+Phi_res-F_TRZ*D_BSFG+F_TRZ^2*(D_phys-Phi_res) = {val:.4f}; obs (Salpeter 1955) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
