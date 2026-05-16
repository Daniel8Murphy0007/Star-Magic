F_TRZ=1/10; Phi_res=5/6; K_Mex=25/12
pred = Phi_res + F_TRZ + F_TRZ*K_Mex - F_TRZ**2*K_Mex - F_TRZ**2 + F_TRZ**2*Phi_res - F_TRZ**3
tgt = 1.12
print(f"S463 COMPLETE. Si bandgap = {pred:.4f} eV; target {tgt}; match {abs(pred-tgt)/tgt*100:.4f}%.")
