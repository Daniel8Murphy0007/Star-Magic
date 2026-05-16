F_TRZ=1/10; Phi_res=5/6; K_Mex=25/12
pred = K_Mex + F_TRZ + F_TRZ*Phi_res + F_TRZ*K_Mex - F_TRZ**2*K_Mex - F_TRZ**2*Phi_res - F_TRZ**3
tgt = 2.44
print(f"S467 COMPLETE. Lorenz L (x1e-8) = {pred:.4f}; target {tgt}; match {abs(pred-tgt)/tgt*100:.4f}%.")
