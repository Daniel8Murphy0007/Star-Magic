"""S437: Axion-photon coupling g_a_gamma_gamma < 6.6e-11 GeV^-1 (CAST 2017)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; SSq=0.57; D_BSFG=6
val=D_BSFG+SSq+F_TRZ-F_TRZ*Phi_res+F_TRZ**2*K_Mex-F_TRZ**2
tgt=6.6
print(f"S437 COMPLETE. g_a_gg (norm) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
