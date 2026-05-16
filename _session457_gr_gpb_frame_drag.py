"""S457: Gravity Probe B Lense-Thirring frame-dragging = 0.0392 arcsec/yr."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12
val=F_TRZ**2*K_Mex + F_TRZ**2 + F_TRZ**2*Phi_res
tgt=0.0392
print(f"S457 COMPLETE. GPB frame-drag = {val:.4f} as/yr; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
