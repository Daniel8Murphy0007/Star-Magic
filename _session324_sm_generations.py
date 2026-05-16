"""S324 (CORRECTED): Three generations from Phi_res^k generation ladder."""
import math
F_TRZ, Phi_res = 0.1, 5/6
# Generations = largest n such that Phi_res^n > 1/2 (mixing-stability threshold)
# Phi_res^3 = 0.579 > 0.5; Phi_res^4 = 0.482 < 0.5 => exactly 3
n_gen = int(math.log(0.5)/math.log(Phi_res))   # floor, no spurious +1
print(f"S324 CORRECTED. Generation count = floor(log(1/2)/log(Phi_res)) = {n_gen}; "
      f"Phi_res^3 = {Phi_res**3:.4f} > 0.5 stable, Phi_res^4 = {Phi_res**4:.4f} < 0.5 decouples; "
      f"exactly 3 mass-stable families.")
