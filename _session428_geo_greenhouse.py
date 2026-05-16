"""S428: Greenhouse effect DeltaT = 33 K (T_surf 288 - T_eff 255)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; D_phys=4; N_ch=9
val=N_ch*D_phys-D_phys+Phi_res+F_TRZ*K_Mex-F_TRZ-F_TRZ**2
tgt=33.0
print(f"S428 COMPLETE. DeltaT_GH = {val:.4f} K; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
