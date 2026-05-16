"""S432: Standard atmospheric pressure 1 atm = 1.013e5 Pa (normalized 1.013)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12
val=F_TRZ*K_Mex+Phi_res-F_TRZ**2-F_TRZ**2*K_Mex
tgt=1.013
print(f"S432 COMPLETE. 1 atm (norm) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
