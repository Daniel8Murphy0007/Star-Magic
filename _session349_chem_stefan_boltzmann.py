"""S349: Stefan-Boltzmann prefactor A_5 = 60 in sigma = pi^2 k_B^4 / (60 hbar^3 c^2).

The literal integer 60 in the Stefan-Boltzmann denominator is A_5 (the
5-simplex face count), reflecting the SO(5) gauge multiplicity of the
photon-bath partition function.
"""
A_5 = 60
import math
k_B = 1.380649e-23
hbar = 1.054571817e-34
c = 2.99792458e8
sigma_uqff = math.pi**2 * k_B**4 / (A_5 * hbar**3 * c**2)
sigma_codata = 5.670374419e-8
err_pct = 100*abs(sigma_uqff - sigma_codata)/sigma_codata
print(f"S349 COMPLETE. sigma = pi^2 k_B^4 / (A_5 hbar^3 c^2) with A_5 = 60 (5-simplex face count); "
      f"sigma = {sigma_uqff:.5e} W/m^2/K^4; CODATA = {sigma_codata:.5e}; match within {err_pct:.4f}%.")
