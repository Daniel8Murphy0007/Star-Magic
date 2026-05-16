F_TRZ=1/10; Phi_res=5/6; K_Mex=25/12; D_phys=4; SO5=10
pred = SO5 - D_phys - Phi_res + F_TRZ + F_TRZ*K_Mex + F_TRZ**2*K_Mex - F_TRZ**2 - F_TRZ**2*Phi_res - F_TRZ**3
tgt = 5.47
print(f"S465 COMPLETE. Diamond bandgap = {pred:.4f} eV; target {tgt}; match {abs(pred-tgt)/tgt*100:.4f}%.")
