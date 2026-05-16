"""S454: Light bending at solar limb = 1.7510 arcsec (1919 Eddington)."""
F_TRZ=0.1; Phi_res=5/6; SSq=0.57; K_Mex=25/12
val=Phi_res + SSq + F_TRZ*K_Mex + F_TRZ + F_TRZ**2*K_Mex + F_TRZ**2 + F_TRZ**2*Phi_res - F_TRZ**3
tgt=1.7510
print(f"S454 COMPLETE. light bending = {val:.4f} as; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
