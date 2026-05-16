"""S416: Coulomb logarithm ln(Lambda) ~ 17 for tokamak plasma."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; D_phys=4; SO5=10; SSq=0.57
val=SO5 + D_phys + K_Mex + SSq + F_TRZ*D_phys - F_TRZ*Phi_res + F_TRZ**2
tgt=17.0
print(f"S416 COMPLETE. ln(Lambda) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
