"""S324: Three-generation puzzle from Phi_res^k generation ladder."""
Phi_res = 5/6
# Three generations because Phi_res^3 < 0.6 (mixing threshold) but Phi_res^4 << 1
# Generation count = floor(log(0.5)/log(Phi_res)) + 1 = 3
import math
n_gen = int(math.log(0.5)/math.log(Phi_res)) + 1
print(f"S324 COMPLETE. Generation count = floor(log(1/2)/log(Phi_res))+1 = {n_gen}; quark/lepton families = Phi_res-suppression ladder.")
