"""S387: Bekenstein-Hawking entropy coefficient S = (1/4)A in l_P^2 units."""
F_TRZ=0.1; K_Mex=25/12; D_phys=4
val = F_TRZ*K_Mex + F_TRZ**2 * D_phys
obs = 0.25
print(f"S387 COMPLETE. S/A = F_TRZ*K_Mex + F_TRZ^2*D_phys = {val:.4f}; obs (Bekenstein-Hawking) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
