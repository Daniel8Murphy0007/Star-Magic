"""S374: Jarlskog CP invariant J."""
F_TRZ=0.1; D_BSFG=6; K_Mex=25/12; SSq=0.57
val = F_TRZ**5 * D_BSFG * SSq * (1 - F_TRZ*K_Mex*SSq)
obs = 3.00e-5
print(f"S374 COMPLETE. J_CP = F_TRZ^5*D_BSFG*SSq*(1-F_TRZ*K_Mex*SSq) = {val:.3e}; obs = {obs:.2e}; match {abs(val-obs)/obs*100:.3f}%.")
