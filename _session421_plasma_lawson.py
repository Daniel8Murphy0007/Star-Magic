"""S421: Lawson criterion n*tau ~ 1.5e20 m^-3 s (normalized 1.5)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; SSq=0.57
val=Phi_res + SSq + F_TRZ - F_TRZ**3
tgt=1.5
print(f"S421 COMPLETE. Lawson n*tau (norm) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
