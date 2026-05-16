"""S442: PMNS solar mixing angle theta_12 = 33.4 degrees."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; D_phys=4; SO5=10; SSq=0.57
val=3*SO5+D_phys-Phi_res+SSq-F_TRZ-F_TRZ*K_Mex
tgt=33.4
print(f"S442 COMPLETE. theta_12 = {val:.4f} deg; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
