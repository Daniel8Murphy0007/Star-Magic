"""S400 BCS isotope effect coefficient alpha = 1/2 (exact)"""
F_TRZ=0.1; Phi_res=5/6; D_phys=4
val = Phi_res - F_TRZ*Phi_res*D_phys     # 5/6 - 4*5/60 = 5/6 - 1/3 = 1/2
obs = 0.5
print(f"S400 COMPLETE. alpha_iso = Phi_res-F_TRZ*Phi_res*D_phys = {val:.6f}; obs (BCS) = {obs}; match {abs(val-obs)/obs*100:.4f}%. (= 1/2 exact)")
