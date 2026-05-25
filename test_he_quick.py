#!/usr/bin/env python3
from superposition_pair_solver import SuperpositionPairStateCalculator, PhysicalConstants
from scipy import constants

const = PhysicalConstants()
solver = SuperpositionPairStateCalculator(const)

# Test Helium
result = solver.solve_pair_system(Z=2, n=1, M_nucleus=4*constants.u, E_neutrino=0.0)
print('Helium (Z=2, n=1):')
print(f'  Single electron energy: {result["single_electron_energy_eV"]:.3f} eV')
print(f'  Pair total energy: {result["pair_total_energy_eV"]:.3f} eV')
print(f'  Expected: -79.0 eV')
print(f'  Error: {abs(result["pair_total_energy_eV"] - (-79.0)):.3f} eV')
print(f'  Pass: {abs(result["pair_total_energy_eV"] - (-79.0)) < 1.0}')
