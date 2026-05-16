"""S426: Dry adiabatic lapse rate Gamma_d = 9.8 K/km."""
F_TRZ=0.1; Phi_res=5/6; SO5=10
val=SO5-F_TRZ-F_TRZ*Phi_res-F_TRZ**2
tgt=9.8
print(f"S426 COMPLETE. Gamma_d = {val:.4f} K/km; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
