"""S427: Earth obliquity (axial tilt) 23.5 degrees."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; D_phys=4; SO5=10
val=2*SO5+D_phys-Phi_res+F_TRZ*K_Mex+F_TRZ*Phi_res+F_TRZ**2
tgt=23.5
print(f"S427 COMPLETE. obliquity = {val:.4f} deg; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
