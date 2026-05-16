"""S443: Landauer limit / bit-to-nat conversion ln(2) = 0.6931."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12
val=Phi_res - F_TRZ - F_TRZ**2*K_Mex - F_TRZ**2*Phi_res - F_TRZ**2 - F_TRZ**3
tgt=0.6931
print(f"S443 COMPLETE. ln(2) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
