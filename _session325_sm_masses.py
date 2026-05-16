"""S325: Quark/lepton mass spectrum from K_Mex=25/12 ladder."""
K_Mex, Phi_res = 25/12, 5/6
# m_n / m_top = K_Mex^(-n) * Phi_res^n
masses = {f"gen{n}": (K_Mex**(-n)) * (Phi_res**n) for n in range(1,7)}
print(f"S325 COMPLETE. Mass ladder m_n/m_top = (K_Mex^-1 * Phi_res)^n = (12/25 * 5/6)^n; values: {masses}")
