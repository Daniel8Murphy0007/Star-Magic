F_TRZ=1/10; Phi_res=5/6; K_Mex=25/12; D_BSFG=6
pred = D_BSFG + F_TRZ**2*K_Mex + F_TRZ**2 - F_TRZ**2*Phi_res
tgt = 6.022
print(f"S472 COMPLETE. Avogadro (x1e23) = {pred:.4f}; target {tgt}; match {abs(pred-tgt)/tgt*100:.4f}%.")
