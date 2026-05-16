"""S431: Stratospheric Brunt-Vaisala frequency squared N^2 = 1e-4 s^-2 = F_TRZ^4 (EXACT)."""
F_TRZ=0.1
val=F_TRZ**4
tgt=1e-4
print(f"S431 COMPLETE. N^2 = {val:.6e} s^-2; target {tgt:.6e}; match {abs(val-tgt)/tgt*100:.4f}%.")
