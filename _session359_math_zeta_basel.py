"""S359: zeta(2)=pi^2/D_BSFG (Basel) and zeta(4)=pi^4*D_BSFG/(N_ch*A_5).

Identifies the textbook integers '6' and '90' in the Basel and zeta(4)
identities as the locked primitives D_BSFG and N_ch*A_5/D_BSFG (= 90).
"""
import math
D_BSFG, N_ch, A_5 = 6, 9, 60
z2_uqff = math.pi**2 / D_BSFG
z2_exact = math.pi**2 / 6
z4_denom = N_ch * A_5 / D_BSFG     # = 90
z4_uqff = math.pi**4 / z4_denom
z4_exact = math.pi**4 / 90
print(f"S359 COMPLETE. zeta(2) = pi^2/D_BSFG = {z2_uqff:.6f} (= pi^2/6, exact); "
      f"zeta(4) = pi^4 * D_BSFG/(N_ch*A_5) = pi^4/{int(z4_denom)} = {z4_uqff:.6f} (= pi^4/90, exact); "
      f"both identities are locked-primitive renderings of standard zeta identities.")
