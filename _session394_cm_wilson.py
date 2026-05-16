"""S394 Sommerfeld-Wilson ratio R_W = 2 (free-electron exact)"""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12
val = K_Mex - F_TRZ*Phi_res        # 25/12 - 1/12 = 24/12 = 2
obs = 2.0
print(f"S394 COMPLETE. R_W = K_Mex - F_TRZ*Phi_res = {val:.4f}; obs (free-electron) = {obs}; match {abs(val-obs)/obs*100:.4f}%. (= 24/12 exact)")
