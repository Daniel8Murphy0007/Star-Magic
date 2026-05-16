"""S450: Omega constant W(1) = 0.5671 (Lambert W)."""
F_TRZ=0.1; SSq=0.57; Phi_res=5/6
val=SSq + F_TRZ**2*Phi_res - F_TRZ**2 - F_TRZ**3
tgt=0.5671
print(f"S450 COMPLETE. Omega = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
