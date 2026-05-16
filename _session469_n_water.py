F_TRZ=1/10; Phi_res=5/6; SSq=0.57; K_Mex=25/12
pred = Phi_res + SSq - F_TRZ + F_TRZ**2*K_Mex + F_TRZ**2*Phi_res + F_TRZ**3
tgt = 1.333
print(f"S469 COMPLETE. n_water = {pred:.4f}; target {tgt}; match {abs(pred-tgt)/tgt*100:.4f}%.")
