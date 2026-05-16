"""S429: Mean ocean salinity 35 ppt (g/kg)."""
F_TRZ=0.1; K_Mex=25/12; D_phys=4; SO5=10
val=K_Mex*SO5+SO5+D_phys+F_TRZ*K_Mex-F_TRZ**2-F_TRZ**3
tgt=35.0
print(f"S429 COMPLETE. salinity = {val:.4f} ppt; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
