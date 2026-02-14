#!/usr/bin/env python3
"""
QCalc_test.py - Unit Tests for Wolfram Extracted Physics Functions
===================================================================

Pytest unit tests for Wolfram physics functions extracted from UQFF framework.

Test Coverage:
- SOURCE14: 12 magnetar physics terms (SGR 0501+4516)
- SOURCE15: 15 SMBH physics terms (Sagittarius A*)
- SOURCE16: 3 star formation terms (Tapestry Nebula)
- SOURCE17: 2 cluster terms (Westerlund 2)
- SOURCE26: 3 cosmological terms (HUDF)
- Integration: QCalc.py UnifiedFieldSolver Wolfram integration

Architecture Validation:
- NO hardcoded system data (all via InputParameters)
- Generic function names (not system-specific)
- EquationResult metadata completeness
- Realistic physics values and units
- m/s² vacuum proof unit consistency

Author: Daniel T. Murphy
Date: February 13, 2026
Updated: February 13, 2026 (Added SOURCE16,17,26 tests)
"""

import pytest
import numpy as np
from IPData import create_manual_input, InputParameters
from QCalc_Wolfram_Extensions import *
from QCalc import UnifiedFieldSolver, ComputeParams, CONSTANTS


class TestSource14MagnetarPhysics:
    """Test 12 magnetar physics terms from source14_wolfram.cpp."""
    
    @pytest.fixture
    def magnetar_params(self):
        """SGR 0501+4516 test parameters."""
        return create_manual_input(
            "SGR 0501+4516 Test",
            M=1.4 * CONSTANTS['M_sun'],  # 1.4 solar masses
            r=20e3,                       # 20 km radius
            B=1e10,                       # 10^10 Tesla
            tau_B=4000 * 3.156e7,         # 4000 years magnetic decay
            tau_Omega=10000 * 3.156e7,    # 10,000 years spin-down
            P=5.0,                        # 5 second rotation period
            rho=1e17,                     # 10^17 kg/m³ density
            v_surf=1e6,                   # 1,000 km/s surface velocity
            delta_x=1e-3,                 # 1 mm position uncertainty
            delta_p=1e-20,                # 10^-20 kg·m/s momentum uncertainty
            psi_integral=1.0,             # Normalized wavefunction
            M_halo=1e29                   # Dark matter halo mass
        )
    
    def test_base_gravity_hubble_magnetic(self, magnetar_params):
        """Test base gravity with Hubble expansion and magnetic suppression."""
        t = 1e8  # 100 million seconds (~3 years)
        result = calculate_base_gravity_hubble_magnetic(magnetar_params, t)
        
        assert result.result > 0, "Gravity should be positive"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert 'Hubble' in result.name, "Name should reference Hubble"
        assert 1e11 < result.result < 1e12, f"Magnetar gravity should be ~10^11 m/s², got {result.result:.3e}"
        assert 'G' in result.parameters_used, "Should track G constant"
        assert 'M' in result.parameters_used, "Should track mass"
    
    def test_uqff_unification_time_reversal(self, magnetar_params):
        """Test UQFF unification with time-reversal factor."""
        Ug1, Ug2, Ug3, Ug4 = 1e9, 1e8, 1e7, 1e6
        result = calculate_uqff_unification_time_reversal(magnetar_params, Ug1, Ug2, Ug3, Ug4)
        
        assert result.result > 0, "Unified field should be positive"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        # F_U = (Ug1+Ug2+Ug3+Ug4) × (1 + f_TRZ) where f_TRZ = 0.1
        expected = (Ug1 + Ug2 + Ug3 + Ug4) * 1.1
        assert result.result == pytest.approx(expected, rel=1e-6), "UQFF unification formula incorrect"
    
    def test_cosmological_constant_acceleration(self, magnetar_params):
        """Test cosmological constant acceleration term."""
        result = calculate_cosmological_constant_acceleration(magnetar_params)
        
        assert result.result > 0, "Lambda acceleration should be positive"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert result.result < 1e-30, f"Lambda term should be tiny (~10^-36), got {result.result:.3e}"
        assert 'Lambda' in result.latex or 'Λ' in result.latex, "Should include Lambda in equation"
    
    def test_em_acceleration_vacuum_corrected(self, magnetar_params):
        """Test EM acceleration with vacuum correction."""
        t = 1e8
        result = calculate_em_acceleration_vacuum_corrected(magnetar_params, t)
        
        assert result.result > 0, "EM acceleration should be positive"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert result.result > 1e10, f"Magnetar EM should be strong (>10^10 m/s²), got {result.result:.3e}"
        assert 'v_surf' in result.parameters_used or 'Bt' in result.parameters_used, "Should use velocity or magnetic field"
    
    def test_gravitational_wave_spin_down(self, magnetar_params):
        """Test GW emission from spin-down."""
        t = 1e8
        result = calculate_gravitational_wave_spin_down(magnetar_params, t)
        
        assert result.result >= 0, "GW strain should be non-negative"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert result.result < 1e-5, f"GW term should be small (<10^-5), got {result.result:.3e}"
        assert 'dOmega_dt' in result.parameters_used, "Should track spin-down rate"
    
    def test_quantum_uncertainty_heisenberg(self, magnetar_params):
        """Test quantum uncertainty contribution."""
        result = calculate_quantum_uncertainty_heisenberg(magnetar_params)
        
        assert result.result > 0, "Quantum term should be positive"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert result.result < 1e-35, f"Quantum term should be tiny (<10^-35), got {result.result:.3e}"
        assert 'hbar' in result.parameters_used, "Should use Planck constant"
        assert 'delta_x' in result.parameters_used, "Should use position uncertainty"
    
    def test_fluid_density_coupling(self, magnetar_params):
        """Test fluid density coupling term."""
        result = calculate_fluid_density_coupling(magnetar_params)
        
        assert result.result > 0, "Fluid coupling should be positive"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert 'rho' in result.parameters_used or 'rho_fluid' in result.parameters_used, "Should use density"
        assert result.result > 1e10, f"Fluid coupling should be significant, got {result.result:.3e}"
    
    def test_oscillatory_wave_superposition(self, magnetar_params):
        """Test wave superposition (standing + traveling)."""
        t, x = 1e8, 1e4
        result = calculate_oscillatory_wave_superposition(magnetar_params, t, x)
        
        # Can be positive or negative (oscillatory)
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert abs(result.result) < 1e12, f"Wave amplitude should be reasonable (<10^12), got {result.result:.3e}"
        assert 'A_osc' in result.parameters_used, "Should track amplitude"
    
    def test_dark_matter_perturbation(self, magnetar_params):
        """Test dark matter perturbation."""
        result = calculate_dark_matter_perturbation(magnetar_params)
        
        assert result.result > 0, "DM perturbation should be positive"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert result.result < 1e10, f"DM perturbation should be small (<10^10), got {result.result:.3e}"
        assert 'M_DM' in result.parameters_used or 'M_halo' in result.parameters_used, "Should use DM mass"
    
    def test_magnetic_field_decay(self, magnetar_params):
        """Test exponential magnetic field decay."""
        # Test at t=0 (initial field)
        result_t0 = calculate_magnetic_field_decay(magnetar_params, t=0)
        B0 = 1e10  # Initial field from params
        assert result_t0.result == pytest.approx(B0, rel=1e-6), "At t=0, B(t) should equal B0"
        
        # Test at t=tau_B (1/e decay)
        tau_B = 4000 * 3.156e7
        result_tau = calculate_magnetic_field_decay(magnetar_params, t=tau_B)
        expected_tau = B0 / np.e
        assert result_tau.result == pytest.approx(expected_tau, rel=0.01), "At t=τ_B, B(t) should equal B0/e"
        assert result_tau.unit == 'T', "Unit should be Tesla"
    
    def test_spin_evolution_angular_velocity(self, magnetar_params):
        """Test spin angular velocity evolution."""
        # Test at t=0 (initial spin)
        result_t0 = calculate_spin_evolution_angular_velocity(magnetar_params, t=0)
        P = 5.0  # 5 second period
        Omega_0 = 2 * np.pi / P
        assert result_t0.result == pytest.approx(Omega_0, rel=1e-6), "At t=0, Ω(t) should equal Ω0"
        
        # Test at t=tau_Omega (1/e decay)
        tau_Omega = 10000 * 3.156e7
        result_tau = calculate_spin_evolution_angular_velocity(magnetar_params, t=tau_Omega)
        expected_tau = Omega_0 / np.e
        assert result_tau.result == pytest.approx(expected_tau, rel=0.01), "At t=τ_Ω, Ω(t) should equal Ω0/e"
        assert result_tau.unit == 'rad/s', "Unit should be rad/s"
    
    def test_time_reversal_factor(self, magnetar_params):
        """Test time-reversal factor constant."""
        result = calculate_time_reversal_factor(magnetar_params)
        
        assert result.result == pytest.approx(0.1, rel=1e-6), "f_TRZ should be 0.1"
        assert result.unit == '(dimensionless)', "Should be dimensionless"


class TestSource15SMBHPhysics:
    """Test 15 SMBH physics terms from source15_wolfram.cpp."""
    
    @pytest.fixture
    def smbh_params(self):
        """Sagittarius A* test parameters."""
        return create_manual_input(
            "Sagittarius A* Test",
            M=4.3e6 * CONSTANTS['M_sun'],  # 4.3 million solar masses
            r=1.27e10,                     # Schwarzschild radius (~12.7 million km)
            B=1e4,                         # 10^4 Gauss (1 Tesla)
            tau_B=1e6 * 3.156e7,           # 1 Myr magnetic decay
            tau_Omega=1e9 * 3.156e7,       # 1 Gyr spin-down
            tau_acc=9e9 * 3.156e7,         # 9 Gyr accretion timescale
            M_dot=0.01,                    # 1% dimensionless accretion rate
            rho=1e-10,                     # Low-density accretion disk
            v_surf=1e5,                    # 100 km/s accretion velocity
            delta_x=1e6,                   # 1,000 km uncertainty
            delta_p=1e-15,                 # Momentum uncertainty
            psi_integral=1.0,              # Normalized wavefunction
            M_halo=4.3e4 * CONSTANTS['M_sun'],  # 1% DM halo
            precession_angle=30.0 * np.pi / 180  # 30° in radians
        )
    
    def test_smbh_time_dependent_mass(self, smbh_params):
        """Test SMBH time-dependent mass M(t) with accretion."""
        # Test at t=0 (initial mass with accretion boost)
        result_t0 = calculate_smbh_time_dependent_mass(smbh_params, t=0)
        M0 = 4.3e6 * CONSTANTS['M_sun']
        M_dot_0 = 0.01
        expected_t0 = M0 * (1.0 + M_dot_0)
        assert result_t0.result == pytest.approx(expected_t0, rel=1e-6), "At t=0, M(t) = M0(1+Ṁ0)"
        
        # Test at t=tau_acc (1/e decay of accretion)
        tau_acc = 9e9 * 3.156e7
        result_tau = calculate_smbh_time_dependent_mass(smbh_params, t=tau_acc)
        expected_tau = M0 * (1.0 + M_dot_0 / np.e)
        assert result_tau.result == pytest.approx(expected_tau, rel=0.01), "At t=τ_acc, accretion decays by e"
        assert result_tau.unit == 'kg', "Unit should be kg"
    
    def test_smbh_base_gravity_mass_evolution(self, smbh_params):
        """Test SMBH base gravity with M(t) evolution."""
        t = 1e12  # 1 trillion seconds
        result = calculate_smbh_base_gravity_mass_evolution(smbh_params, t)
        
        assert result.result > 0, "Gravity should be positive"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert 1e6 < result.result < 1e7, f"SMBH gravity should be ~10^6 m/s², got {result.result:.3e}"
        assert 'Mt' in result.parameters_used, "Should track time-dependent mass"
    
    def test_smbh_uqff_unification(self, smbh_params):
        """Test SMBH UQFF unification (same formula as magnetar)."""
        Ug1, Ug2, Ug3, Ug4 = 1e5, 1e4, 1e3, 1e2
        result = calculate_smbh_uqff_unification(smbh_params, Ug1, Ug2, Ug3, Ug4)
        
        expected = (Ug1 + Ug2 + Ug3 + Ug4) * 1.1
        assert result.result == pytest.approx(expected, rel=1e-6), "SMBH UQFF formula incorrect"
        assert 'SMBH' in result.name, "Name should indicate SMBH"
    
    def test_smbh_cosmological_constant(self, smbh_params):
        """Test SMBH cosmological constant (same as magnetar)."""
        result = calculate_smbh_cosmological_constant(smbh_params)
        
        assert result.result > 0, "Lambda acceleration should be positive"
        assert result.result < 1e-30, "Lambda term should be tiny"
    
    def test_smbh_em_acceleration(self, smbh_params):
        """Test SMBH EM acceleration (accretion disk)."""
        t = 1e12
        result = calculate_smbh_em_acceleration(smbh_params, t)
        
        assert result.result > 0, "EM acceleration should be positive"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert result.result > 1e10, f"SMBH EM should be strong, got {result.result:.3e}"
    
    def test_smbh_gravitational_wave(self, smbh_params):
        """Test SMBH GW emission with M(t) dependence."""
        t = 1e12
        result = calculate_smbh_gravitational_wave(smbh_params, t)
        
        assert result.result >= 0, "GW strain should be non-negative"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert result.result < 1e-10, f"SMBH GW should be tiny, got {result.result:.3e}"
    
    def test_smbh_quantum_uncertainty(self, smbh_params):
        """Test SMBH quantum uncertainty (SMBH scale)."""
        result = calculate_smbh_quantum_uncertainty(smbh_params)
        
        assert result.result > 0, "Quantum term should be positive"
        assert result.result < 1e-40, f"SMBH quantum term should be tiny, got {result.result:.3e}"
    
    def test_smbh_fluid_density(self, smbh_params):
        """Test SMBH fluid density with M(t) dependence."""
        t = 1e12
        result = calculate_smbh_fluid_density(smbh_params, t)
        
        assert result.result > 0, "Fluid coupling should be positive"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert 'Mt' in result.parameters_used, "Should use time-dependent mass"
    
    def test_smbh_oscillatory_wave_orbital(self, smbh_params):
        """Test SMBH orbital oscillations (light-crossing time)."""
        t, x = 1e12, 1e9
        result = calculate_smbh_oscillatory_wave_orbital(smbh_params, t, x)
        
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert abs(result.result) < 1e8, f"SMBH wave amplitude reasonable (<10^8), got {result.result:.3e}"
        assert 'light_crossing_time' in result.parameters_used, "Should use light-crossing time"
    
    def test_smbh_dark_matter_precession(self, smbh_params):
        """Test SMBH DM with 30° precession angle."""
        result = calculate_smbh_dark_matter_precession(smbh_params)
        
        assert result.result > 0, "DM perturbation should be positive"
        assert result.unit == 'm/s²', "Unit should be m/s²"
        assert 'precession_factor' in result.parameters_used, "Should use precession sin(30°)"
        assert 0.4 < result.parameters_used['precession_factor'] < 0.6, "sin(30°) ≈ 0.5"
    
    def test_smbh_magnetic_decay_gauss_conversion(self, smbh_params):
        """Test SMBH magnetic decay with Gauss→Tesla conversion."""
        # Test at t=0
        result_t0 = calculate_smbh_magnetic_decay_gauss_conversion(smbh_params, t=0)
        B0_gauss = 1e4  # 10^4 Gauss
        B0_tesla = B0_gauss * 1e-4  # 1 Tesla
        assert result_t0.result == pytest.approx(B0_tesla, rel=1e-6), "At t=0, B(t) = B0 in Tesla"
        assert result_t0.unit == 'T', "Unit should be Tesla"
    
    def test_smbh_spin_evolution_relativistic(self, smbh_params):
        """Test SMBH spin with relativistic initial velocity."""
        result_t0 = calculate_smbh_spin_evolution_relativistic(smbh_params, t=0)
        
        # Ω0 = 0.3c/r
        c = CONSTANTS['c']
        r = 1.27e10
        spin_factor = 0.3
        Omega_0_expected = spin_factor * c / r
        assert result_t0.result == pytest.approx(Omega_0_expected, rel=1e-6), "Ω0 = 0.3c/r"
        assert result_t0.unit == 'rad/s', "Unit should be rad/s"
    
    def test_smbh_precession_factor(self, smbh_params):
        """Test SMBH precession factor sin(30°) = 0.5."""
        result = calculate_smbh_precession_factor(smbh_params)
        
        assert result.result == pytest.approx(0.5, rel=1e-6), "sin(30°) = 0.5"
        assert result.unit == '(dimensionless)', "Should be dimensionless"
    
    def test_smbh_accretion_rate(self, smbh_params):
        """Test SMBH accretion rate exponential decay."""
        # Test at t=0
        result_t0 = calculate_smbh_accretion_rate(smbh_params, t=0)
        M_dot_0 = 0.01
        assert result_t0.result == pytest.approx(M_dot_0, rel=1e-6), "At t=0, Ṁ(t) = Ṁ0"
        
        # Test at t=tau_acc
        tau_acc = 9e9 * 3.156e7
        result_tau = calculate_smbh_accretion_rate(smbh_params, t=tau_acc)
        expected_tau = M_dot_0 / np.e
        assert result_tau.result == pytest.approx(expected_tau, rel=0.01), "At t=τ_acc, Ṁ decays by e"
    
    def test_smbh_schwarzschild_radius(self, smbh_params):
        """Test SMBH Schwarzschild radius calculation."""
        result = calculate_smbh_schwarzschild_radius(smbh_params)
        
        G = CONSTANTS['G']
        M = 4.3e6 * CONSTANTS['M_sun']
        c = CONSTANTS['c']
        r_s_expected = (2 * G * M) / (c ** 2)
        
        assert result.result == pytest.approx(r_s_expected, rel=1e-6), "r_s = 2GM/c²"
        assert result.unit == 'm', "Unit should be meters"
        assert 1e10 < result.result < 2e10, f"Sgr A* r_s ≈ 1.27×10^10 m, got {result.result:.3e}"


class TestSource16StarFormation:
    """Test 3 star formation functions from SOURCE16."""
    
    @pytest.fixture
    def star_formation_params(self):
        """Tapestry Nebula star formation test parameters."""
        return create_manual_input(
            "Tapestry Star Formation Test",
            M=5e3 * CONSTANTS['M_sun'],  # 5,000 solar masses
            r=5e16,                      # ~1.5 pc radius
            SFR=0.05,                    # 0.05 M_sun/year star formation rate
            rho=1e-18,                   # 10^-18 kg/m³ cloud density
            v_surf=1e4                   # 10 km/s stellar wind velocity
        )
    
    def test_star_formation_mass_growth(self, star_formation_params):
        """Test star formation mass growth rate."""
        t = 1e6 * 3.156e7  # 1 Myr
        result = calculate_star_formation_mass_growth(star_formation_params, t)
        
        assert result.result != 0, "Mass growth should be non-zero"
        assert result.unit == 'm/s²', "Unit should be m/s² (vacuum proof)"
        assert 'SFR' in result.parameters_used or 'M' in result.parameters_used, "Should track SFR or mass"
    
    def test_stellar_wind_ram_pressure(self, star_formation_params):
        """Test stellar wind ram pressure contribution."""
        result = calculate_stellar_wind_ram_pressure(star_formation_params)
        
        assert result.result != 0, "Wind pressure should be non-zero"
        assert result.unit == 'm/s²', "Unit should be m/s² (vacuum proof)"
    
    def test_tapestry_radiation_pressure(self, star_formation_params):
        """Test Tapestry radiation pressure term."""
        result = calculate_tapestry_radiation_pressure(star_formation_params)
        
        assert result.result != 0, "Radiation pressure should be non-zero"
        assert result.unit == 'm/s²', "Unit should be m/s² (vacuum proof)"


class TestSource17Clusters:
    """Test 2 cluster functions from SOURCE17."""
    
    @pytest.fixture
    def cluster_params(self):
        """Westerlund 2 cluster test parameters."""
        return create_manual_input(
            "Westerlund 2 Test",
            M=1e4 * CONSTANTS['M_sun'],  # 10,000 solar masses
            r=2e16,                      # ~0.6 pc radius
            age=2e6 * 3.156e7            # 2 Myr age
        )
    
    def test_cluster_mass_evolution(self, cluster_params):
        """Test cluster mass evolution over time - returns mass M(t), not acceleration."""
        t = 1e6 * 3.156e7  # 1 Myr
        result = calculate_cluster_mass_evolution(cluster_params, t)
        
        assert result.result != 0, "Cluster evolution should produce non-zero mass"
        assert result.unit == 'kg', "Unit should be kg (mass evolution)"
        assert result.result > 0, "Mass should be positive"
    
    def test_westerlund2_composite_muge(self, cluster_params):
        """Test Westerlund2 composite MUGE calculation."""
        t = 1e6 * 3.156e7  # 1 Myr
        result = calculate_westerlund2_composite_muge(cluster_params, t)
        
        assert result.result != 0, "Composite MUGE should be non-zero"
        assert result.unit == 'm/s²', "Unit should be m/s² (vacuum proof)"


class TestSource26Cosmological:
    """Test 3 HUDF cosmological functions from SOURCE26."""
    
    @pytest.fixture
    def hudf_params(self):
        """Hubble Ultra Deep Field test parameters."""
        return create_manual_input(
            "HUDF Galaxy Test",
            M=1e11 * CONSTANTS['M_sun'],  # 100 billion solar masses
            r=1e22,                        # ~3 Mpc
            z=6.5,                         # High redshift
            SFR=100.0                      # 100 M_sun/year
        )
    
    def test_hudf_star_formation_mass(self, hudf_params):
        """Test HUDF star formation mass contribution - returns M(t), not acceleration."""
        t = 1e9 * 3.156e7  # 1 Gyr
        result = calculate_hudf_star_formation_mass(hudf_params, t)
        
        assert result.result != 0, "HUDF SFR contribution should be non-zero"
        assert result.unit == 'kg', "Unit should be kg (mass from star formation)"
        assert result.result > 0, "Mass should be positive"
    
    def test_hudf_intergalaxy_interaction(self, hudf_params):
        """Test HUDF intergalaxy interaction term - returns coupling strength I(t)."""
        result = calculate_hudf_intergalaxy_interaction(hudf_params)
        
        assert result.result != 0, "Intergalaxy interaction should be non-zero"
        assert result.unit == 'dimensionless', "Unit should be dimensionless (coupling strength)"
        assert 0 <= result.result <= 1.0, "Coupling strength should be normalized"
    
    def test_hudf_complete_muge(self, hudf_params):
        """Test HUDF complete MUGE calculation."""
        t = 1e9 * 3.156e7  # 1 Gyr
        result = calculate_hudf_complete_muge(hudf_params, t)
        
        assert result.result != 0, "Complete MUGE should be non-zero"
        assert result.unit == 'm/s²', "Unit should be m/s² (vacuum proof)"


class TestQCalcIntegration:
    """Test integration of Wolfram functions into QCalc.py UnifiedFieldSolver."""
    
    def test_wolfram_integration_magnetar(self):
        """Test that QCalc.solve() calls Wolfram functions for magnetar."""
        solver = UnifiedFieldSolver()
        
        # Magnetar parameters
        params = ComputeParams(
            query_name="Magnetar Integration Test",
            M=1.4 * CONSTANTS['M_sun'],
            r=20e3,
            B=1e10,
            t=1e8
        )
        
        result = solver.solve(params)
        
        # Check that some Wolfram results are present
        eq_names = [eq['name'] for eq in result['long_form_equations']]
        
        # Should have some magnetar terms
        assert any('Base Gravity' in name and 'Hubble' in name for name in eq_names), \
            "Should include Wolfram base gravity term"
    
    def test_wolfram_integration_smbh(self):
        """Test that QCalc.solve() calls Wolfram functions for SMBH."""
        solver = UnifiedFieldSolver()
        
        # SMBH parameters
        params = ComputeParams(
            query_name="SMBH Integration Test",
            M=4.3e6 * CONSTANTS['M_sun'],
            r=1.27e10,
            B=1e4,
            t=1e12
        )
        
        result = solver.solve(params)
        
        # Check that some SMBH Wolfram results are present
        eq_names = [eq['name'] for eq in result['long_form_equations']]
        
        # Should have some SMBH terms
        assert any('SMBH' in name for name in eq_names), \
            "Should include SMBH Wolfram terms"
    
    def test_wolfram_graceful_failure(self):
        """Test that QCalc handles missing Wolfram dependencies gracefully."""
        solver = UnifiedFieldSolver()
        
        # Minimal parameters (may not have all Wolfram requirements)
        params = ComputeParams(
            query_name="Minimal Test",
            M=1e30,
            r=1e6
        )
        
        # Should not raise exception
        result = solver.solve(params)
        
        # Should have some results even if Wolfram terms fail
        assert len(result['long_form_equations']) > 0, "Should compute something"
    
    def test_equation_result_structure(self):
        """Test that all Wolfram functions return properly structured EquationResults."""
        magnetar_params = create_manual_input(
            "Structure Test",
            M=1.4 * CONSTANTS['M_sun'],
            r=20e3,
            B=1e10
        )
        
        result = calculate_base_gravity_hubble_magnetic(magnetar_params, t=1e8)
        
        # Validate EquationResult structure
        assert hasattr(result, 'name'), "Should have name"
        assert hasattr(result, 'latex'), "Should have latex equation"
        assert hasattr(result, 'substituted'), "Should have substituted equation"
        assert hasattr(result, 'result'), "Should have numerical result"
        assert hasattr(result, 'unit'), "Should have unit"
        assert hasattr(result, 'parameters_used'), "Should have parameters_used dict"
        assert hasattr(result, 'notes'), "Should have notes"
        
        assert isinstance(result.name, str), "Name should be string"
        assert isinstance(result.result, (int, float, np.number)), "Result should be numeric"
        assert isinstance(result.parameters_used, dict), "Parameters should be dict"


# ═══════════════════════════════════════════════════════════════════════════════
# RUN TESTS
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # Run pytest with verbose output
    import sys
    sys.exit(pytest.main([__file__, '-v', '--tb=short']))
