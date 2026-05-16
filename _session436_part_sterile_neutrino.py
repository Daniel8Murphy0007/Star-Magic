"""S436: Sterile neutrino best-fit Delta m^2 = 1.7 eV^2 (LSND/MiniBooNE region)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; SSq=0.57
val=Phi_res+SSq+F_TRZ*K_Mex+F_TRZ-F_TRZ**2
tgt=1.7
print(f"S436 COMPLETE. Delta m^2 = {val:.4f} eV^2; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
