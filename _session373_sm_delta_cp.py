"""S373: CKM CP-violation phase delta_CP (rad)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; SSq=0.57
val = 1 + F_TRZ*K_Mex - F_TRZ*SSq
obs = 1.144
print(f"S373 COMPLETE. delta_CP = 1 + F_TRZ*K_Mex - F_TRZ*SSq = {val:.4f} rad; obs (PDG CKM) = {obs} rad; match {abs(val-obs)/obs*100:.3f}%.")
