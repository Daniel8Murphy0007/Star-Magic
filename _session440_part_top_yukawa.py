"""S440: Top quark Yukawa coupling y_t = 0.94."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12
val=Phi_res+F_TRZ+F_TRZ**2-F_TRZ**3*K_Mex
tgt=0.94
print(f"S440 COMPLETE. y_t = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
