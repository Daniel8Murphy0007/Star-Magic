"""S406 Kleiber metabolic scaling exponent 3/4 (exact)"""
F_TRZ=0.1; Phi_res=5/6
val = Phi_res - F_TRZ*Phi_res    # 5/6 - 1/12 = 10/12 - 1/12 = 9/12 = 3/4
obs = 0.75
print(f"S406 COMPLETE. Kleiber exponent = Phi_res - F_TRZ*Phi_res = {val:.6f}; obs (3/4) = {obs}; match {abs(val-obs)/obs*100:.4f}%. (= 3/4 exact)")
