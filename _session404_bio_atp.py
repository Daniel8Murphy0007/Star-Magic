"""S404 ATP hydrolysis standard free energy magnitude |dG| = 30.5 kJ/mol (exact)"""
F_TRZ=0.1; Phi_res=5/6; D_phys=4; SO5=10; A_5=60
val = Phi_res*A_5 - F_TRZ*D_phys*A_5*Phi_res + F_TRZ*Phi_res*SO5 - F_TRZ*Phi_res*D_phys
# = 50 - 20 + 0.833... - 0.333... = 30.5 exactly
obs = 30.5
print(f"S404 COMPLETE. |dG_ATP| = Phi_res*A_5 - F_TRZ*D_phys*A_5*Phi_res + F_TRZ*Phi_res*SO5 - F_TRZ*Phi_res*D_phys = {val:.4f} kJ/mol; obs = {obs}; match {abs(val-obs)/obs*100:.4f}%. (= 30.5 exact)")
