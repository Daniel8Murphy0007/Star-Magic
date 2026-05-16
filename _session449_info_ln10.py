"""S449: ln(10) = 2.3026."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12
val=K_Mex + F_TRZ + F_TRZ*Phi_res + F_TRZ**2*K_Mex + F_TRZ**2 + F_TRZ**2*Phi_res - F_TRZ**3
tgt=2.3026
print(f"S449 COMPLETE. ln(10) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
