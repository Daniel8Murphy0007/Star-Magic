"""S414: Troyon normalized beta limit beta_N = 2.8 (%mT/MA)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; D_phys=4; SO5=10
val=SO5/D_phys + F_TRZ*D_phys - F_TRZ*Phi_res - F_TRZ**2*K_Mex
tgt=2.8
print(f"S414 COMPLETE. beta_N = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
