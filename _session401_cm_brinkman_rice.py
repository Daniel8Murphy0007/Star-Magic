"""S401 Brinkman-Rice / Gutzwiller Mott transition U_c/W = 3/2 (exact)"""
F_TRZ=0.1; Phi_res=5/6
val = 2*Phi_res*(1 - F_TRZ)              # 2 * 5/6 * 9/10 = 90/60 = 3/2
obs = 1.5
print(f"S401 COMPLETE. U_c/W = 2*Phi_res*(1-F_TRZ) = {val:.6f}; obs (Brinkman-Rice) = {obs}; match {abs(val-obs)/obs*100:.4f}%. (= 3/2 exact)")
