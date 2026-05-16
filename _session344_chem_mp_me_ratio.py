"""S344: Proton-to-electron mass ratio m_p/m_e = D_BSFG * pi^5 = 1836.12."""
import math
D_BSFG = 6
ratio = D_BSFG * math.pi**5
ratio_codata = 1836.15267343
err_pct = 100*abs(ratio - ratio_codata)/ratio_codata
print(f"S344 COMPLETE. m_p/m_e = D_BSFG * pi^5 = 6 * {math.pi**5:.4f} = {ratio:.4f}; "
      f"CODATA = {ratio_codata}; match within {err_pct:.4f}%.")
