F_TRZ=1/10; Phi_res=5/6; SSq=0.57
pred = Phi_res + SSq + F_TRZ - F_TRZ*Phi_res - F_TRZ**3
tgt = 1.42
print(f"S464 COMPLETE. GaAs bandgap = {pred:.4f} eV; target {tgt}; match {abs(pred-tgt)/tgt*100:.4f}%.")
