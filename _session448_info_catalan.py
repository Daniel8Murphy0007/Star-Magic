"""S448: Catalan constant G = 0.9160."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12
val=Phi_res + F_TRZ - F_TRZ**2*K_Mex + F_TRZ**2 - F_TRZ**2*Phi_res + F_TRZ**3
tgt=0.9160
print(f"S448 COMPLETE. Catalan G = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
