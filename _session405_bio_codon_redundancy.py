"""S405 codon redundancy 64 codons / 20 amino acids = 3.2"""
F_TRZ=0.1; Phi_res=5/6; SSq=0.57; K_Mex=25/12; D_phys=4
val = K_Mex + SSq + F_TRZ*D_phys + SSq*F_TRZ + F_TRZ*Phi_res
obs = 64.0/20.0
print(f"S405 COMPLETE. codons/aa = K_Mex+SSq+F_TRZ*D_phys+SSq*F_TRZ+F_TRZ*Phi_res = {val:.4f}; obs (64/20) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
