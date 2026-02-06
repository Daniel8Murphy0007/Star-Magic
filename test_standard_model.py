#!/usr/bin/env python3
"""Test Standard Model UQFF Integration in CondensedPhysics.py."""

from CondensedPhysics import STANDARD_MODEL_UQFF, CONSTANTS

print('=' * 70)
print('STANDARD MODEL UQFF INTEGRATION VALIDATION')
print('=' * 70)

# Test 1: Quark calculations
print('\n1. QUARK CALCULATIONS')
print('-' * 50)
for flavor in ['u', 'd', 'c', 's', 't', 'b']:
    U_q, _ = STANDARD_MODEL_UQFF.compute_quark_energy_density(flavor, t=0, t_n=0, n=1)
    E_q, _ = STANDARD_MODEL_UQFF.compute_quark_mass_energy(flavor, n=1)
    print(f'  {flavor} quark: U = {U_q:.2e} J/m³, E = {E_q:.2e} J')

# Test 2: Lepton calculations
print('\n2. LEPTON CALCULATIONS')
print('-' * 50)
for lepton in ['e', 'mu', 'tau']:
    U_l, _ = STANDARD_MODEL_UQFF.compute_lepton_energy_density(lepton, t=0, t_n=0, n=2)
    print(f'  {lepton}: U = {U_l:.2e} J/m³')

# Test 3: Gluon calculations
print('\n3. GLUON CALCULATIONS')
print('-' * 50)
U_g, _ = STANDARD_MODEL_UQFF.compute_gluon_energy_density(t=0, n=1)
print(f'  Gluon: U = {U_g:.2e} J/m³')
E_color, _ = STANDARD_MODEL_UQFF.compute_color_charge_field(r=1e-15)
print(f'  Color field at 1 fm: {E_color:.2e}')

# Test 4: Higgs calculations
print('\n4. HIGGS CALCULATIONS')
print('-' * 50)
U_H, _ = STANDARD_MODEL_UQFF.compute_higgs_energy_density(t=0, n=1)
print(f'  Higgs: U = {U_H:.2e} J/m³')
V_H, _ = STANDARD_MODEL_UQFF.compute_higgs_potential(phi_H=246)
print(f'  Higgs potential at VEV: V = {V_H:.2e} GeV⁴')
m_t, _ = STANDARD_MODEL_UQFF.compute_mass_from_yukawa('t')
print(f'  Top mass from Yukawa: m_t = {m_t:.1f} GeV/c²')

# Test 5: Hadron calculations
print('\n5. HADRON CALCULATIONS')
print('-' * 50)
U_p, _ = STANDARD_MODEL_UQFF.compute_proton_energy_density(t=0, t_n=0, n=1)
U_n, _ = STANDARD_MODEL_UQFF.compute_neutron_energy_density(t=0, t_n=0, n=1)
print(f'  Proton: U = {U_p:.2e} J/m³')
print(f'  Neutron: U = {U_n:.2e} J/m³')
U_K, _ = STANDARD_MODEL_UQFF.compute_kaon_energy_density('K_plus', t=0, t_n=0, n=3)
print(f'  K⁺ kaon: U = {U_K:.2e} J/m³')

# Test 6: Quantum state effects
print('\n6. QUANTUM STATE EFFECTS (26 STATES)')
print('-' * 50)
E_n1, _ = STANDARD_MODEL_UQFF.compute_mass_energy_modified(CONSTANTS['m_e'], n=1)
E_n13, _ = STANDARD_MODEL_UQFF.compute_mass_energy_modified(CONSTANTS['m_e'], n=13)
E_n26, _ = STANDARD_MODEL_UQFF.compute_mass_energy_modified(CONSTANTS['m_e'], n=26)
print(f'  Electron at n=1:  E = {E_n1:.4e} J (factor 1.038)')
print(f'  Electron at n=13: E = {E_n13:.4e} J (factor 1.649)')
print(f'  Electron at n=26: E = {E_n26:.4e} J (factor 2.718 = e)')

# Test 7: Non-local jump probability
print('\n7. NON-LOCAL INTERACTIONS')
print('-' * 50)
P_jump, _ = STANDARD_MODEL_UQFF.compute_nonlocal_jump_probability(t=0, n=1)
print(f'  Jump probability (t=0, n=1): P = {P_jump:.4e}')
t_minus, _ = STANDARD_MODEL_UQFF.compute_negative_time_operator(t_n=1.0)
print(f'  Negative time operator (t_n=1): t⁻ = {t_minus:.4e}')

# Test 8: Neutrino oscillation
print('\n8. NEUTRINO OSCILLATION')
print('-' * 50)
P_osc, _ = STANDARD_MODEL_UQFF.compute_neutrino_oscillation_probability(
    L=1000e3,  # 1000 km baseline
    E_nu=1e-12,  # 1 GeV in J (approximate)
    flavor_from='e', flavor_to='mu'
)
print(f'  P(νe → νμ) at 1000 km, 1 GeV: P = {P_osc:.4f}')

# Test 9: Decay rates
print('\n9. DECAY RATES')
print('-' * 50)
for particle in ['neutron', 'K_plus', 'muon']:
    Gamma, _ = STANDARD_MODEL_UQFF.compute_decay_rate(particle, t=0)
    tau = STANDARD_MODEL_UQFF.decay_lifetimes[particle]
    print(f'  {particle}: τ = {tau:.2e} s, Γ(0) = {Gamma:.2e} s⁻¹')

# Validation tests
print('\n' + '=' * 70)
print('VALIDATION TESTS')
print('=' * 70)
results = STANDARD_MODEL_UQFF.validate_standard_model_physics()
for test, status, value in results['results']:
    print(f'  {test}: {status} ({value})')
print(f'\nTotal: {results["tests_passed"]}/{results["total_tests"]} tests PASS')
print(f'Pass rate: {results["pass_rate"]:.0f}%')

# Summary
print('\n' + '=' * 70)
print('STANDARD MODEL PARAMETERS LOADED')
print('=' * 70)
print(f'  Quark masses: {len(STANDARD_MODEL_UQFF.quark_masses)} flavors')
print(f'  Lepton masses: {len(STANDARD_MODEL_UQFF.lepton_masses)} types')
print(f'  Boson masses: {len(STANDARD_MODEL_UQFF.boson_masses)} bosons')
print(f'  Hadron masses: {len(STANDARD_MODEL_UQFF.hadron_masses)} hadrons')
print(f'  Decay lifetimes: {len(STANDARD_MODEL_UQFF.decay_lifetimes)} particles')
print(f'  Quantum states: 1-{STANDARD_MODEL_UQFF.n_quantum_states}')
