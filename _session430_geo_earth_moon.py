"""S430: Earth-Moon semimajor axis / Earth radius = 60.34."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; A_5=60
val=A_5+F_TRZ*K_Mex+F_TRZ*Phi_res+F_TRZ**2
tgt=60.34
print(f"S430 COMPLETE. a_Moon/R_E = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
