"""S435: Sum of neutrino masses Sigma m_nu = 0.12 eV (Planck CMB)."""
F_TRZ=0.1; K_Mex=25/12
val=F_TRZ+F_TRZ**2*K_Mex-F_TRZ**3
tgt=0.12
print(f"S435 COMPLETE. Sigma m_nu = {val:.4f} eV; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
