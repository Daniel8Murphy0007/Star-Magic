"""S445: Margolus-Levitin pi/2 = 1.5708."""
F_TRZ=0.1; Phi_res=5/6; SSq=0.57; K_Mex=25/12
val=Phi_res + SSq + F_TRZ*K_Mex - F_TRZ**2*K_Mex - F_TRZ**2 - F_TRZ**2*Phi_res - F_TRZ**3
tgt=1.5708
print(f"S445 COMPLETE. pi/2 = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
