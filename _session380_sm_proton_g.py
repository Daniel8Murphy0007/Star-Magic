"""S380: Proton g-factor (magnetic moment in nuclear magnetons * 2)."""
F_TRZ=0.1; D_BSFG=6; Phi_res=5/6; D_phys=4
val = D_BSFG - Phi_res + F_TRZ*D_phys
obs = 5.5857
print(f"S380 COMPLETE. g_p = D_BSFG - Phi_res + F_TRZ*D_phys = {val:.4f}; obs (CODATA) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
