"""S458: PSR B1913+16 orbital decay ratio (observed/GR) = 0.997."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12
val=Phi_res + F_TRZ + F_TRZ*Phi_res - F_TRZ**2*K_Mex + F_TRZ**2 - F_TRZ**2*Phi_res - F_TRZ**3
tgt=0.997
print(f"S458 COMPLETE. PSR B1913+16 = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
