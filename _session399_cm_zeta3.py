"""S399 Apery constant zeta(3) = 1.2020569 (appears in lattice sums, Sommerfeld expansions)"""
F_TRZ=0.1; Phi_res=5/6; D_phys=4
val = Phi_res + F_TRZ*D_phys - F_TRZ*F_TRZ*Phi_res*D_phys
obs = 1.2020569
print(f"S399 COMPLETE. zeta(3) closure = Phi_res+F_TRZ*D_phys-F_TRZ^2*Phi_res*D_phys = {val:.4f}; obs (Apery) = {obs:.4f}; match {abs(val-obs)/obs*100:.3f}%.")
