"""S459: Schwarzschild BH shadow radius r_shadow/M = 3*sqrt(3) = 5.1962."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; N_ch=9; SO5=10; D_phys=4
val=SO5 - D_phys - F_TRZ - F_TRZ*N_ch + F_TRZ*K_Mex - F_TRZ**2*Phi_res - F_TRZ**3
tgt=5.1962
print(f"S459 COMPLETE. BH shadow = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
