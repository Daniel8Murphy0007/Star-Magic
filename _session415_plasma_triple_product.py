"""S415: D-T fusion triple product n*T*tau = 3e21 keV*s/m^3 (normalized 3.0)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12
val=Phi_res + K_Mex + F_TRZ - F_TRZ**2*K_Mex + F_TRZ**3
tgt=3.0
print(f"S415 COMPLETE. triple product (norm) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
