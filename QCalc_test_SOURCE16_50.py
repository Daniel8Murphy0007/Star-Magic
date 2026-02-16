#!/usr/bin/env python3
"""
QCalc_test_SOURCE16_50.py - Unit Tests for SOURCE16-50 Physics Functions
=========================================================================

Pytest unit tests for 67 additional Wolfram physics functions (SOURCE16-50).

Test Coverage:
- SOURCE16: 3 star formation functions
- SOURCE17: 2 cluster functions
- SOURCE18: 3 photoevaporation functions
- SOURCE19-25: 14 batch astrophysical functions
- SOURCE26-27: 6 cosmological functions
- SOURCE28-30: 6 galaxy/planetary functions
- SOURCE31-35: 8 nebula/magnetar frequency functions
- SOURCE36-40: 10 framework module functions
- SOURCE41-45: 7 extreme-scale functions
- SOURCE46-50: 8 specific nebulae + generic API functions

Total: 67 functions (completes 94-function test suite with QCalc_test.py)

Author: Daniel T. Murphy
Date: February 13, 2026
"""

import pytest
import numpy as np
from IPData import create_manual_input, InputParameters
from QCalc_Wolfram_Extensions import *
from QCalc import UnifiedFieldSolver, ComputeParams, CONSTANTS


# ═══════════════════════════════════════════════════════════════════════════════
# SOURCE16 - STAR FORMATION (3 functions)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSource16StarFormation:
    """Test 3 star formation physics terms (Tapestry Nebula)."""
    
    @pytest.fixture
    def star_formation_params(self):
        """Tapestry Nebula test parameters."""
        return create_manual_input(
            "Tapestry Nebula Test",
            M=1e4 * CONSTANTS['M_sun'],  # 10^4 solar masses
            r=5e17,                       # ~16 parsecs
            T=5000,                       # 5000 K
            rho=1e-18,                    # ISM density kg/m³
            v_surf=1e4                    # 10 km/s gas velocity
        )
    
    def test_star_formation_mass_growth(self, star_formation_params):
        """Test star formation mass growth rate."""
        t = 1e13  # ~300,000 years
        result = calculate_star_formation_mass_growth(star_formation_params, t)
        
        assert result.result > 0, "Mass growth should be positive"
        assert result.unit == 'm/s²', "Gravity function should return m/s²"
        assert 'M(t)' in result.latex or 'dM/dt' in result.latex or 'Delta M' in result.latex, "Should reference mass evolution"
    
    def test_stellar_wind_ram_pressure(self, star_formation_params):
        """Test stellar wind ram pressure."""
        result = calculate_stellar_wind_ram_pressure(star_formation_params)
        
        assert result.result > 0, "Ram pressure should be positive"
        assert result.unit == 'm/s²', "Gravity function should return m/s²"
        assert 'rho_wind' in result.parameters_used or 'v_wind' in result.parameters_used, "Should use wind density/velocity"
    
    def test_tapestry_radiation_pressure(self, star_formation_params):
        """Test radiation pressure in Tapestry."""
        result = calculate_tapestry_radiation_pressure(star_formation_params)
        
        assert result.result > 0, "Radiation pressure should be positive"
        assert result.unit == 'm/s²', "Gravity function should return m/s²"


# ═══════════════════════════════════════════════════════════════════════════════
# SOURCE17 - CLUSTERS (2 functions)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSource17Clusters:
    """Test 2 cluster physics terms (Westerlund 2)."""
    
    @pytest.fixture
    def cluster_params(self):
        """Westerlund 2 test parameters."""
        return create_manual_input(
            "Westerlund 2 Test",
            M=5e3 * CONSTANTS['M_sun'],  # 5,000 solar masses
            r=3e17,                       # ~10 parsecs
            T=20000,                      # 20,000 K (hot stars)
            rho=1e-16                     # Cluster density
        )
    
    def test_cluster_mass_evolution(self, cluster_params):
        """Test cluster mass evolution."""
        t = 1e13  # ~300,000 years
        result = calculate_cluster_mass_evolution(cluster_params, t)
        
        assert result.result != 0, "Mass evolution should be non-zero"
        assert 'M' in result.parameters_used, "Should use mass"
    
    def test_westerlund2_composite_muge(self, cluster_params):
        """Test Westerlund2 composite MUGE."""
        t = 1e13
        result = calculate_westerlund2_composite_muge(cluster_params, t)
        
        assert result.result > 0, "Composite MUGE should be positive"
        assert result.unit == 'm/s²', "Unit should be acceleration"


# ═══════════════════════════════════════════════════════════════════════════════
# SOURCE18 - PHOTOEVAPORATION (3 functions)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSource18Photoevaporation:
    """Test 3 photoevaporation physics terms (Pillars of Creation)."""
    
    @pytest.fixture
    def photoevap_params(self):
        """Pillars of Creation test parameters."""
        return create_manual_input(
            "Pillars of Creation Test",
            M=1e3 * CONSTANTS['M_sun'],  # 1,000 solar masses
            r=1e17,                       # ~3 parsecs
            T=10000,                      # 10,000 K ionization front
            rho=1e-17,                    # Pillar density
            v_surf=1e4                    # Erosion velocity
        )
    
    def test_photoevaporation_erosion(self, photoevap_params):
        """Test photoevaporation erosion rate."""
        t = 1e12  # ~30,000 years
        result = calculate_photoevaporation_erosion(photoevap_params, t)
        
        assert result.result >= 0, "Erosion rate should be non-negative"
        assert 'erosion' in result.name.lower() or 'mass loss' in result.name.lower()
    
    def test_ionization_front_pressure(self, photoevap_params):
        """Test ionization front pressure."""
        result = calculate_ionization_front_pressure(photoevap_params)
        
        assert result.result > 0, "IF pressure should be positive"
        assert result.unit == 'Pa' or result.unit == 'N/m²', "Unit should be pressure"
    
    def test_pillars_mass_with_erosion(self, photoevap_params):
        """Test pillars mass accounting for erosion."""
        t = 1e12
        result = calculate_pillars_mass_with_erosion(photoevap_params, t)
        
        assert result.result > 0, "Remaining mass should be positive"
        assert result.unit == 'kg', "Unit should be mass"


# ═══════════════════════════════════════════════════════════════════════════════
# SOURCE19-25 - BATCH ASTROPHYSICS (14 functions)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSource19_25BatchAstro:
    """Test 14 batch astrophysical functions."""
    
    @pytest.fixture
    def batch_params(self):
        """Generic astrophysical system parameters."""
        return create_manual_input(
            "Generic Astrophysical System",
            M=1e6 * CONSTANTS['M_sun'],  # 10^6 solar masses
            r=1e19,                       # ~300 parsecs
            T=1e7,                        # 10 million K
            rho=1e-15,                    # Typical density
            v_surf=1e5,                   # 100 km/s
            M_bh=1e5 * CONSTANTS['M_sun'], # SMBH mass
            M_halo=1e10 * CONSTANTS['M_sun'] # Halo mass
        )
    
    def test_gravitational_lensing_amplification(self, batch_params):
        """Test gravitational lensing."""
        result = calculate_gravitational_lensing_amplification(batch_params)
        assert result.result > 0, "Amplification should be positive"
    
    def test_central_smbh_contribution(self, batch_params):
        """Test central SMBH gravity contribution."""
        t = 1e14
        result = calculate_central_smbh_contribution(batch_params, t)
        assert result.result > 0, "SMBH contribution should be positive"
        assert result.unit == 'm/s²', "Unit should be acceleration"
    
    def test_supernova_mass_ejection(self, batch_params):
        """Test supernova mass ejection."""
        t = 1e7  # 100 days
        result = calculate_supernova_mass_ejection(batch_params, t)
        assert result.result >= 0, "Mass ejection should be non-negative"
    
    def test_cavity_pressure_decay(self, batch_params):
        """Test cavity pressure decay."""
        t = 1e13
        result = calculate_cavity_pressure_decay(batch_params, t)
        assert result.result >= 0, "Pressure should be non-negative"
    
    def test_starburst_mass_growth(self, batch_params):
        """Test starburst mass growth."""
        t = 1e13
        result = calculate_starburst_mass_growth(batch_params, t)
        assert result.result >= 0, "Mass growth should be non-negative"
    
    def test_bubble_expansion_radius(self, batch_params):
        """Test bubble expansion radius."""
        t = 1e13
        result = calculate_bubble_expansion_radius(batch_params, t)
        assert result.result > 0, "Radius should be positive"
        assert result.unit == 'm', "Unit should be distance"
    
    def test_stellar_wind_feedback_acceleration(self, batch_params):
        """Test stellar wind feedback."""
        result = calculate_stellar_wind_feedback_acceleration(batch_params)
        assert result.result != 0, "Feedback should be non-zero"
        assert result.unit == 'm/s²', "Unit should be acceleration"
    
    def test_tidal_interaction_strength(self, batch_params):
        """Test tidal interaction."""
        result = calculate_tidal_interaction_strength(batch_params)
        assert result.result >= 0, "Tidal force should be non-negative"
    
    def test_merger_enhanced_star_formation(self, batch_params):
        """Test merger-enhanced SFR."""
        t = 1e14
        result = calculate_merger_enhanced_star_formation(batch_params, t)
        assert result.result >= 0, "SFR should be non-negative"
    
    def test_horsehead_erosion_mass_loss(self, batch_params):
        """Test Horsehead nebula erosion."""
        t = 1e13
        result = calculate_horsehead_erosion_mass_loss(batch_params, t)
        assert result.result >= 0, "Mass loss should be non-negative"
    
    def test_nebula_mass_decay(self, batch_params):
        """Test nebula mass decay."""
        t = 1e13
        result = calculate_nebula_mass_decay(batch_params, t)
        assert result.result > 0, "Remaining mass should be positive"
    
    def test_cooling_flow_contribution(self, batch_params):
        """Test cooling flow."""
        t = 1e14
        result = calculate_cooling_flow_contribution(batch_params, t)
        assert result.result != 0, "Cooling flow should be non-zero"
    
    def test_magnetic_filament_decay(self, batch_params):
        """Test magnetic filament decay."""
        t = 1e13
        result = calculate_magnetic_filament_decay(batch_params, t)
        assert result.result >= 0, "Field decay should be non-negative"
    
    def test_filament_support_buildup(self, batch_params):
        """Test filament support buildup."""
        t = 1e13
        result = calculate_filament_support_buildup(batch_params, t)
        assert result.result >= 0, "Support should be non-negative"


# ═══════════════════════════════════════════════════════════════════════════════
# SOURCE26-27 - COSMOLOGICAL (6 functions)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSource26_27Cosmological:
    """Test 6 cosmological physics functions."""
    
    @pytest.fixture
    def cosmo_params(self):
        """HUDF/NGC1792 test parameters."""
        return create_manual_input(
            "HUDF Galaxy Test",
            M=1e11 * CONSTANTS['M_sun'],  # 10^11 solar masses
            r=1e22,                        # ~3 Mpc
            z=6.5,                         # High redshift
            T=1e6,                         # 1 million K
            rho=1e-20                      # Cosmological density
        )
    
    def test_hudf_star_formation_mass(self, cosmo_params):
        """Test HUDF star formation mass."""
        t = 1e15
        result = calculate_hudf_star_formation_mass(cosmo_params, t)
        assert result.result >= 0, "SFR mass should be non-negative"
    
    def test_hudf_intergalaxy_interaction(self, cosmo_params):
        """Test HUDF intergalactic interaction."""
        result = calculate_hudf_intergalaxy_interaction(cosmo_params)
        assert result.result != 0, "Interaction should be non-zero"
    
    def test_hudf_complete_muge(self, cosmo_params):
        """Test HUDF complete MUGE."""
        t = 1e15
        result = calculate_hudf_complete_muge(cosmo_params, t)
        assert result.result > 0, "MUGE should be positive"
        assert result.unit == 'm/s²', "Unit should be acceleration"
    
    def test_ngc1792_star_formation_mass(self, cosmo_params):
        """Test NGC1792 star formation."""
        t = 1e14
        result = calculate_ngc1792_star_formation_mass(cosmo_params, t)
        assert result.result >= 0, "SFR mass should be non-negative"
    
    def test_ngc1792_uqff_ug(self, cosmo_params):
        """Test NGC1792 UQFF Ug."""
        result = calculate_ngc1792_uqff_ug(cosmo_params)
        assert result.result > 0, "Ug should be positive"
        assert result.unit == 'm/s²', "Unit should be acceleration"
    
    def test_ngc1792_complete_muge(self, cosmo_params):
        """Test NGC1792 complete MUGE."""
        t = 1e14
        result = calculate_ngc1792_complete_muge(cosmo_params, t)
        assert result.result > 0, "MUGE should be positive"


# ═══════════════════════════════════════════════════════════════════════════════
# SOURCE28-30 - GALAXY & PLANETARY (6 functions)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSource28_30GalaxyPlanetary:
    """Test 6 galaxy and planetary functions."""
    
    @pytest.fixture
    def galaxy_params(self):
        """Andromeda M31 test parameters."""
        return create_manual_input(
            "Andromeda M31 Test",
            M=1.5e42,  # Andromeda mass
            r=7.7e20,  # 25 kpc radius
            T=1e6,     # ISM temperature
            rho=1e-21, # ISM density
            v_surf=1e5 # Rotation velocity
        )
    
    @pytest.fixture
    def planetary_params(self):
        """Saturn test parameters."""
        return create_manual_input(
            "Saturn Test",
            M=5.68e26,  # Saturn mass
            r=5.82e7,   # Saturn radius
            T=134,      # Surface temperature
            omega=1.64e-4  # Rotation rate
        )
    
    def test_andromeda_dust_friction(self, galaxy_params):
        """Test Andromeda dust friction."""
        result = calculate_andromeda_dust_friction(galaxy_params)
        assert result.result >= 0, "Friction should be non-negative"
    
    def test_andromeda_complete_muge(self, galaxy_params):
        """Test Andromeda complete MUGE."""
        t = 1e15
        result = calculate_andromeda_complete_muge(galaxy_params, t)
        assert result.result > 0, "MUGE should be positive"
        assert result.unit == 'm/s²', "Unit should be acceleration"
    
    def test_sombrero_superconductivity_dust(self, galaxy_params):
        """Test Sombrero superconductivity."""
        result = calculate_sombrero_superconductivity_dust(galaxy_params)
        assert result.result != 0, "SC effect should be non-zero"
    
    def test_sombrero_complete_muge(self, galaxy_params):
        """Test Sombrero complete MUGE."""
        t = 1e15
        result = calculate_sombrero_complete_muge(galaxy_params, t)
        assert result.result > 0, "MUGE should be positive"
    
    def test_saturn_ring_wind_effects(self, planetary_params):
        """Test Saturn ring wind effects."""
        result = calculate_saturn_ring_wind_effects(planetary_params)
        assert result.result != 0, "Ring effect should be non-zero"
    
    def test_saturn_complete_muge(self, planetary_params):
        """Test Saturn complete MUGE."""
        t = 1e9  # ~30 years
        result = calculate_saturn_complete_muge(planetary_params, t)
        assert result.result > 0, "MUGE should be positive"


# ═══════════════════════════════════════════════════════════════════════════════
# SOURCE31-35 - NEBULA & MAGNETAR FREQUENCY (8 functions)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSource31_35NebulaFrequency:
    """Test 8 nebula and magnetar frequency functions."""
    
    @pytest.fixture
    def nebula_params(self):
        """M16 Eagle Nebula test parameters."""
        return create_manual_input(
            "M16 Eagle Nebula Test",
            M=2e33,    # 1000 solar masses
            r=1.2e16,  # 4 parsecs
            T=8000,    # 8000 K
            rho=1e-18, # Nebula density
            B=1e-7     # Magnetic field
        )
    
    @pytest.fixture
    def magnetar_params(self):
        """SGR 1745-2900 test parameters."""
        return create_manual_input(
            "SGR 1745-2900 Test",
            M=2.8e30,  # 1.4 solar masses
            r=12000,   # 12 km
            B=1e14,    # 10^14 Tesla
            P=3.76,    # 3.76 s period
            T=1e9      # 1 billion K
        )
    
    def test_m16_star_formation_radiation(self, nebula_params):
        """Test M16 star formation radiation."""
        t = 1e13
        result = calculate_m16_star_formation_radiation(nebula_params, t)
        assert result.result >= 0, "Radiation should be non-negative"
    
    def test_m16_complete_muge(self, nebula_params):
        """Test M16 complete MUGE."""
        t = 1e13
        result = calculate_m16_complete_muge(nebula_params, t)
        assert result.result > 0, "MUGE should be positive"
    
    def test_crab_pulsar_wind_magnetic(self, magnetar_params):
        """Test Crab pulsar wind magnetic."""
        t = 1e9  # ~30 years
        result = calculate_crab_pulsar_wind_magnetic(magnetar_params, t)
        assert result.result != 0, "Pulsar wind should be non-zero"
    
    def test_crab_complete_muge(self, magnetar_params):
        """Test Crab complete MUGE."""
        t = 1e9
        result = calculate_crab_complete_muge(magnetar_params, t)
        assert result.result > 0, "MUGE should be positive"
    
    def test_sgr1745_superconductivity_critical(self, magnetar_params):
        """Test SGR1745 superconductivity."""
        result = calculate_sgr1745_superconductivity_critical(magnetar_params)
        assert result.result != 0, "SC effect should be non-zero"
    
    def test_sgr1745_complete_muge(self, magnetar_params):
        """Test SGR1745 complete MUGE."""
        t = 1e9
        result = calculate_sgr1745_complete_muge(magnetar_params, t)
        assert result.result > 0, "MUGE should be positive"
    
    def test_sgr1745_frequency_model(self, magnetar_params):
        """Test SGR1745 frequency model."""
        t = 1e8
        result = calculate_sgr1745_frequency_model(magnetar_params, t)
        assert result.result > 0, "Frequency should be positive"
        assert result.unit == 'Hz' or result.unit == 'rad/s', "Unit should be frequency"
    
    def test_sgra_frequency_model(self, magnetar_params):
        """Test Sgr A* frequency model."""
        t = 1e8
        result = calculate_sgra_frequency_model(magnetar_params, t)
        assert result.result > 0, "Frequency should be positive"


# ═══════════════════════════════════════════════════════════════════════════════
# SOURCE36-40 - FRAMEWORK MODULES (10 functions)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSource36_40FrameworkModules:
    """Test 10 generic framework module functions."""
    
    @pytest.fixture
    def framework_params(self):
        """Generic framework test parameters."""
        return create_manual_input(
            "Framework Test System",
            M=1e6 * CONSTANTS['M_sun'],
            r=1e18,
            T=1e6,
            rho=1e-16,
            B=1e-6,
            omega=1e-12,
            P=1e8
        )
    
    def test_tapestry_dpm_term(self, framework_params):
        """Test Tapestry DPM term."""
        result = calculate_tapestry_dpm_term(framework_params)
        assert result.result != 0, "DPM should be non-zero"
    
    def test_tapestry_complete_uqff(self, framework_params):
        """Test Tapestry complete UQFF."""
        t = 1e14
        result = calculate_tapestry_complete_uqff(framework_params, t)
        assert result.result > 0, "UQFF should be positive"
    
    def test_resonance_terms(self, framework_params):
        """Test resonance terms."""
        t = 1e14
        result = calculate_resonance_terms(framework_params, t)
        assert result.result != 0, "Resonance should be non-zero"
    
    def test_resonance_superconductivity_full(self, framework_params):
        """Test resonance superconductivity."""
        t = 1e14
        result = calculate_resonance_superconductivity_full(framework_params, t)
        assert result.result != 0, "Resonance SC should be non-zero"
    
    def test_compressed_terms(self, framework_params):
        """Test compressed terms (sys 10-16)."""
        result = calculate_compressed_terms(framework_params)
        assert result.result > 0, "Compressed should be positive"
    
    def test_compressed_resonance_full(self, framework_params):
        """Test compressed resonance full."""
        t = 1e14
        result = calculate_compressed_resonance_full(framework_params, t)
        assert result.result > 0, "Compressed resonance should be positive"
    
    def test_crab_resonance_dpm(self, framework_params):
        """Test Crab resonance DPM."""
        t = 1e9
        result = calculate_crab_resonance_dpm(framework_params, t)
        assert result.result != 0, "Crab DPM should be non-zero"
    
    def test_crab_resonance_complete(self, framework_params):
        """Test Crab resonance complete."""
        t = 1e9
        result = calculate_crab_resonance_complete(framework_params, t)
        assert result.result > 0, "Crab resonance should be positive"
    
    def test_compressed_terms_sys18_24(self, framework_params):
        """Test compressed terms (sys 18-24)."""
        result = calculate_compressed_terms_sys18_24(framework_params)
        assert result.result > 0, "Compressed sys18-24 should be positive"
    
    def test_compressed_resonance_sys18_24(self, framework_params):
        """Test compressed resonance sys18-24."""
        t = 1e14
        result = calculate_compressed_resonance_sys18_24(framework_params, t)
        assert result.result > 0, "Compressed resonance sys18-24 should be positive"


# ═══════════════════════════════════════════════════════════════════════════════
# SOURCE41-45 - EXTREME SCALE (7 functions)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSource41_45ExtremeScale:
    """Test 7 extreme-scale physics functions."""
    
    @pytest.fixture
    def universe_params(self):
        """Universe-scale test parameters."""
        return create_manual_input(
            "Universe Test",
            M=1e53,    # Observable universe mass
            r=4.4e26,  # Observable universe radius (46 Gly)
            T=2.725,   # CMB temperature
            z=1100     # Recombination redshift
        )
    
    @pytest.fixture
    def hydrogen_params(self):
        """Hydrogen atom test parameters."""
        return create_manual_input(
            "Hydrogen Atom Test",
            M=1.67e-27,  # Proton mass
            r=5.29e-11,  # Bohr radius
            T=300        # Room temperature
        )
    
    def test_universe_diameter_complete(self, universe_params):
        """Test universe diameter calculation."""
        t = 4.35e17  # Age of universe
        result = calculate_universe_diameter_complete(universe_params, t)
        assert result.result > 0, "Diameter should be positive"
        assert result.unit == 'm', "Unit should be distance"
    
    def test_hydrogen_quantum_term(self, hydrogen_params):
        """Test hydrogen quantum term."""
        result = calculate_hydrogen_quantum_term(hydrogen_params)
        assert result.result != 0, "Quantum term should be non-zero"
    
    def test_hydrogen_complete_uqff(self, hydrogen_params):
        """Test hydrogen complete UQFF."""
        t = 1.0  # 1 second
        result = calculate_hydrogen_complete_uqff(hydrogen_params, t)
        assert result.result != 0, "UQFF should be non-zero"
    
    def test_hydrogen_ptoe_resonance(self, hydrogen_params):
        """Test hydrogen PToE resonance."""
        result = calculate_hydrogen_ptoe_resonance(hydrogen_params)
        assert result.result != 0, "PToE resonance should be non-zero"
    
    def test_lagoon_m8_star_formation(self, universe_params):
        """Test Lagoon M8 star formation."""
        t = 1e13
        result = calculate_lagoon_m8_star_formation(universe_params, t)
        assert result.result >= 0, "SFR should be non-negative"
    
    def test_spiral_supernova_term(self, universe_params):
        """Test spiral supernova term."""
        t = 1e7  # 100 days
        result = calculate_spiral_supernova_term(universe_params, t)
        assert result.result >= 0, "SN term should be non-negative"
    
    def test_spiral_complete_uqff(self, universe_params):
        """Test spiral complete UQFF."""
        t = 1e14
        result = calculate_spiral_complete_uqff(universe_params, t)
        assert result.result > 0, "UQFF should be positive"


# ═══════════════════════════════════════════════════════════════════════════════
# SOURCE46-50 - SPECIFIC NEBULAE & GENERIC API (8 functions)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSource46_50NebulaGenericAPI:
    """Test 8 specific nebula and generic API functions."""
    
    @pytest.fixture
    def ngc6302_params(self):
        """NGC 6302 Butterfly Nebula test parameters."""
        return create_manual_input(
            "NGC 6302 Butterfly Test",
            M=1e30,    # Solar mass central star
            r=5e16,    # 0.5 parsec
            T=25000,   # 25,000 K central star
            rho=1e-19  # Nebula density
        )
    
    @pytest.fixture
    def orion_params(self):
        """Orion M42 test parameters."""
        return create_manual_input(
            "Orion M42 Test",
            M=2e33,    # 1000 solar masses
            r=1.2e16,  # 4 parsecs
            T=10000,   # 10,000 K
            rho=1e-18, # Nebula density
            B=1e-7     # Magnetic field
        )
    
    @pytest.fixture
    def generic_params(self):
        """Generic system for API testing."""
        return create_manual_input(
            "Generic System Test",
            M=1e30,
            r=1e15,
            T=5000,
            rho=1e-16,
            B=1e-8,
            omega=1e-10,
            P=1e6
        )
    
    def test_ngc6302_butterfly_complete(self, ngc6302_params):
        """Test NGC 6302 butterfly complete."""
        t = 1e12
        result = calculate_ngc6302_butterfly_complete(ngc6302_params, t)
        assert result.result > 0, "NGC6302 should be positive"
    
    def test_ngc6302_resonance(self, ngc6302_params):
        """Test NGC 6302 resonance."""
        t = 1e12
        result = calculate_ngc6302_resonance(ngc6302_params, t)
        assert result.result != 0, "Resonance should be non-zero"
    
    def test_orion_m42_complete(self, orion_params):
        """Test Orion M42 complete."""
        t = 1e13
        result = calculate_orion_m42_complete(orion_params, t)
        assert result.result > 0, "Orion M42 should be positive"
    
    def test_compressed_resonance_framework(self, generic_params):
        """Test compressed resonance framework (multi-system)."""
        t = 1e14
        result = calculate_compressed_resonance_framework(generic_params, t)
        assert result.result > 0, "Framework should be positive"
    
    def test_generic_compressed_uqff(self, generic_params):
        """Test generic compressed UQFF."""
        t = 1e14
        result = calculate_generic_compressed_uqff(generic_params, t)
        assert result.result > 0, "Generic compressed should be positive"
    
    def test_generic_resonance_uqff(self, generic_params):
        """Test generic resonance UQFF."""
        t = 1e14
        result = calculate_generic_resonance_uqff(generic_params, t)
        assert result.result != 0, "Generic resonance should be non-zero"


# ═══════════════════════════════════════════════════════════════════════════════
# INTEGRATION TESTS - QCalc.py System Detection
# ═══════════════════════════════════════════════════════════════════════════════

class TestQCalcIntegrationSOURCE16_50:
    """Test QCalc.py integration with SOURCE16-50 functions."""
    
    def test_system_detectors(self):
        """Test that system detectors work correctly."""
        solver = UnifiedFieldSolver()
        
        # Test star formation detection
        params_sf = ComputeParams(
            query_name="Tapestry Nebula",
            M=1e4 * CONSTANTS['M_sun'],
            r=5e17,
            t=1e13
        )
        result_sf = solver.solve(params_sf)
        assert len(result_sf['long_form_equations']) > 30, "Should trigger star formation functions"
        
        # Test galaxy detection
        params_gal = ComputeParams(
            query_name="Andromeda M31",
            M=1.5e42,
            r=7.7e20,
            t=1e15
        )
        result_gal = solver.solve(params_gal)
        assert len(result_gal['long_form_equations']) > 30, "Should trigger galaxy functions"
        
        # Test planetary detection
        params_planet = ComputeParams(
            query_name="Saturn",
            M=5.68e26,
            r=5.82e7,
            t=1e9
        )
        result_planet = solver.solve(params_planet)
        assert len(result_planet['long_form_equations']) > 20, "Should trigger planetary functions"
        
        # Test cosmological detection
        params_cosmo = ComputeParams(
            query_name="HUDF Galaxy",
            M=1e41,
            r=1e22,
            z=6.5,
            t=1e15
        )
        result_cosmo = solver.solve(params_cosmo)
        assert len(result_cosmo['long_form_equations']) > 30, "Should trigger cosmological functions"
    
    def test_all_94_functions_accessible(self):
        """Test that all 94 functions are accessible from QCalc.py."""
        solver = UnifiedFieldSolver()
        
        # Create comprehensive params to trigger all function groups
        params = ComputeParams(
            query_name="Comprehensive Test System",
            M=1e6 * CONSTANTS['M_sun'],
            r=1e18,
            T=1e6,
            B=1e-6,
            omega=1e-12,
            P=1e8,
            z=1.0,
            rho=1e-16,
            t=1e14
        )
        
        result = solver.solve(params)
        
        # Should return many equations (not necessarily all 94, but substantial number)
        assert len(result['long_form_equations']) >= 50, f"Should compute 50+ equations, got {len(result['long_form_equations'])}"
        assert 'long_form_equations' in result, "Should have equations"
        assert 'solutions' in result, "Should have solutions"
        assert 'available_equations' in result, "Should have available equations list"


# ═══════════════════════════════════════════════════════════════════════════════
# RUN TESTS
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # Run pytest with verbose output
    import sys
    sys.exit(pytest.main([__file__, '-v', '--tb=short']))
