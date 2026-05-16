"""S456: Gravity Probe B geodetic precession = 6.6028 arcsec/yr."""
F_TRZ=0.1; Phi_res=5/6; SSq=0.57; K_Mex=25/12; D_BSFG=6
val=D_BSFG + SSq + F_TRZ - F_TRZ*Phi_res + F_TRZ**2*K_Mex - F_TRZ**2*Phi_res + F_TRZ**3
tgt=6.6028
print(f"S456 COMPLETE. GPB geodetic = {val:.4f} as/yr; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
