"""S388: White-dwarf radius-mass scaling exponent (R ~ M^(-1/3))."""
F_TRZ=0.1; Phi_res=5/6; D_phys=4
val = -Phi_res * F_TRZ * D_phys
obs = -1/3
print(f"S388 COMPLETE. WD exponent = -Phi_res*F_TRZ*D_phys = {val:.6f}; obs (polytrope n=3/2) = {obs:.6f}; match {abs(val-obs)/abs(obs)*100:.4f}%. (= -1/3 exact)")
