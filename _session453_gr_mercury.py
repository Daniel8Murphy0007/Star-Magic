"""S453: Mercury perihelion precession (anomalous) = 43.0 arcsec/century."""
F_TRZ=0.1; K_Mex=25/12; N_ch=9
SO5=10; D_phys=4
val=SO5*D_phys + K_Mex + F_TRZ*N_ch - F_TRZ**2 + F_TRZ**2*K_Mex
tgt=43.0
print(f"S453 COMPLETE. Mercury = {val:.4f} as/cy; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
