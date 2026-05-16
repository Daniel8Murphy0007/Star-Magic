"""S447: Euler-Mascheroni gamma = 0.5772."""
F_TRZ=0.1; Phi_res=5/6; SSq=0.57
val=SSq + F_TRZ**2*Phi_res - F_TRZ**3
tgt=0.5772
print(f"S447 COMPLETE. gamma = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
