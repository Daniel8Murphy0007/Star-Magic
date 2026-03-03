#!/usr/bin/env python3
"""
test_priority2_integration.py - Priority 2 Integration Test Suite
=================================================================

Tests UQFFSystemsDatabase integration with CondensedPhysics2.py and
validates database functionality for UQFF calculations.

Priority 2 Validation:
1. UQFFSystemsDatabase initialization (36+ systems)
2. System parameter access methods
3. Category filtering and queries
4. CondensedPhysics2 import integration
5. Combined UQFF+Database calculations

Author: Daniel T. Murphy
Created: March 3, 2026
"""

import sys
import traceback
from typing import Dict, List


def test_database_initialization():
    """Test 1: UQFFSystemsDatabase initialization."""
    print("\n" + "="*70)
    print("TEST 1: UQFFSystemsDatabase Initialization")
    print("="*70)
    
    try:
        from UQFFSystemsDatabase import UQFFSystemsDatabase, UQFFPhysicalConstants
        
        db = UQFFSystemsDatabase()
        systems_count = len(db.list_systems())
        categories_count = len(db.list_categories())
        
        print(f"✓ Database initialized successfully")
        print(f"✓ Total systems: {systems_count}")
        print(f"✓ Total categories: {categories_count}")
        
        assert systems_count >= 36, f"Expected >= 36 systems, got {systems_count}"
        print(f"✓ System count validation passed (>= 36)")
        
        return True
    except Exception as e:
        print(f"✗ FAILED: {e}")
        traceback.print_exc()
        return False


def test_system_access_methods():
    """Test 2: System parameter access methods."""
    print("\n" + "="*70)
    print("TEST 2: System Access Methods")
    print("="*70)
    
    try:
        from UQFFSystemsDatabase import UQFFSystemsDatabase
        
        db = UQFFSystemsDatabase()
        
        # Test get_system
        m87 = db.get_system('M87')
        assert m87 is not None, "M87 system not found"
        print(f"✓ get_system('M87'): {m87.name}")
        
        sgr = db.get_system('SGR1745')
        assert sgr is not None, "SGR1745 system not found"
        print(f"✓ get_system('SGR1745'): {sgr.name}")
        
        at2024 = db.get_system('AT2024tvd')
        assert at2024 is not None, "AT2024tvd system not found"
        print(f"✓ get_system('AT2024tvd'): {at2024.name}")
        
        # Verify key parameters
        assert m87.mass == 6.5e9, f"M87 mass incorrect: {m87.mass}"
        print(f"✓ M87 SMBH mass: {m87.mass:.2e} M☉")
        
        assert sgr.B == 2e14, f"SGR1745 B-field incorrect: {sgr.B}"
        print(f"✓ SGR1745 B-field: {sgr.B:.2e} Gauss")
        
        assert at2024.L_X == 5e43, f"AT2024tvd L_X incorrect: {at2024.L_X}"
        print(f"✓ AT2024tvd X-ray luminosity: {at2024.L_X:.2e} erg/s")
        
        return True
    except Exception as e:
        print(f"✗ FAILED: {e}")
        traceback.print_exc()
        return False


def test_category_filtering():
    """Test 3: Category filtering and queries."""
    print("\n" + "="*70)
    print("TEST 3: Category Filtering")
    print("="*70)
    
    try:
        from UQFFSystemsDatabase import UQFFSystemsDatabase
        
        db = UQFFSystemsDatabase()
        
        # Test category queries
        agn_systems = db.get_systems_by_category('AGN')
        print(f"✓ AGN systems: {len(agn_systems)} found")
        agn_names = [s.name for s in agn_systems]
        print(f"  Systems: {', '.join(agn_names)}")
        
        tde_systems = db.get_systems_by_category('TDE')
        print(f"✓ TDE systems: {len(tde_systems)} found")
        tde_names = [s.name for s in tde_systems]
        print(f"  Systems: {', '.join(tde_names)}")
        
        wr_systems = db.get_systems_by_category('Wolf-Rayet')
        print(f"✓ Wolf-Rayet systems: {len(wr_systems)} found")
        wr_names = [s.name for s in wr_systems]
        print(f"  Systems: {', '.join(wr_names)}")
        
        # Test list methods
        all_systems = db.list_systems()
        all_categories = db.list_categories()
        print(f"✓ list_systems(): {len(all_systems)} systems")
        print(f"✓ list_categories(): {len(all_categories)} categories")
        
        return True
    except Exception as e:
        print(f"✗ FAILED: {e}")
        traceback.print_exc()
        return False


def test_condensedphysics2_import():
    """Test 4: CondensedPhysics2 integration."""
    print("\n" + "="*70)
    print("TEST 4: CondensedPhysics2 Import Integration")
    print("="*70)
    
    try:
        from CondensedPhysics2 import (
            UQFFPhysicalConstants,
            AstrophysicalSystem,
            UQFFSystemsDatabase,
        )
        
        print("✓ UQFFPhysicalConstants imported from CondensedPhysics2")
        print("✓ AstrophysicalSystem imported from CondensedPhysics2")
        print("✓ UQFFSystemsDatabase imported from CondensedPhysics2")
        
        # Verify constants
        print(f"  c = {UQFFPhysicalConstants.c:.3e} m/s")
        print(f"  G = {UQFFPhysicalConstants.G:.3e} m³/kg·s²")
        print(f"  M_sun = {UQFFPhysicalConstants.M_sun:.3e} kg")
        
        # Initialize database via CP2 import
        db = UQFFSystemsDatabase()
        print(f"✓ Database initialized via CondensedPhysics2: {len(db.list_systems())} systems")
        
        return True
    except Exception as e:
        print(f"✗ FAILED: {e}")
        traceback.print_exc()
        return False


def test_combined_uqff_database():
    """Test 5: Combined UQFF calculations with database parameters."""
    print("\n" + "="*70)
    print("TEST 5: Combined UQFF + Database Calculations")
    print("="*70)
    
    try:
        from UQFFSystemsDatabase import UQFFSystemsDatabase
        from GrokThreadUQFFExtensions import (
            SystemParams, 
            UniversalMagnetismCalculator,
            ResonanceGravityCalculator
        )
        from BuoyancyProofVariants import FUBiiVirialXray
        
        db = UQFFSystemsDatabase()
        
        # Test 1: M87 with UniversalMagnetismCalculator
        print("\n--- M87 Universal Magnetism ---")
        m87 = db.get_system('M87')
        um_calc = UniversalMagnetismCalculator()
        
        # Convert M87 to SystemParams
        m87_params = SystemParams(
            M=m87.mass,  # 6.5e9 M☉
            r_obs=m87.distance * 3.086e16,  # pc to meters
            B=m87.B if m87.B else 1e-3,  # Gauss
            T=1e6,  # Assume 1 MK accretion disk
            v=m87.v_rot * 1e3 if m87.v_rot else 550e3,  # km/s to m/s
            z=m87.redshift,
            SFR=0.0,  # AGN, not star-forming
            tau_c=1e6  # years
        )
        
        um_result = um_calc.compute_Um_with_metadata(m87_params, r=1e15, t=1e6, t_n=0.5)
        print(f"✓ M87 Um = {um_result['Um']:.3e} (normalized)")
        print(f"  LENR enhancement factor: {um_result['LENR_enhancement']:.3e}")
        print(f"  Quasi-particle modulation: {um_result['quasi_particle_mod']:.6f}")
        
        # Test 2: Perseus Cluster with FUBiiVirialXray
        print("\n--- Perseus Cluster Buoyancy (Virial X-ray) ---")
        perseus = db.get_system('PerseusCluster')
        buoy_calc = FUBiiVirialXray()
        
        # Use cluster velocity dispersion and X-ray data
        sigma_X = perseus.additional_params.get('velocity_dispersion', 1200) * 1e3  # km/s to m/s
        T_X = perseus.T  # Already in Kelvin
        n_e = perseus.n_e * 1e6  # cm⁻³ to m⁻³
        r = 1e22  # 1 Mpc in meters
        
        F_UBii_virx = buoy_calc.compute(sigma_X, T_X, n_e, r)
        print(f"✓ Perseus F_UBii_virx = {F_UBii_virx:.3e}")
        print(f"  σ_X = {sigma_X:.3e} m/s")
        print(f"  T_X = {T_X:.3e} K")
        print(f"  n_e = {n_e:.3e} m⁻³")
        
        # Test 3: SGR1745 with ResonanceGravityCalculator
        print("\n--- SGR1745 Resonance Gravity ---")
        sgr = db.get_system('SGR1745')
        g_res_calc = ResonanceGravityCalculator()
        
        sgr_params = SystemParams(
            M=sgr.mass,  # 1.4 M☉
            r_obs=sgr.distance * 3.086e16,  # pc to meters
            B=sgr.B,  # 2e14 Gauss
            T=sgr.T,  # 5e6 K
            v=1e5,  # Assume 100 km/s surface velocity
            z=0.0,  # Galactic
            SFR=0.0,
            tau_c=1e3  # Magnetar cooling timescale (years)
        )
        
        g_res_result = g_res_calc.compute_g_res_with_metadata(sgr_params, r=1e4, t=1e9)
        print(f"✓ SGR1745 g_res = {g_res_result['g_res']:.3e} m/s²")
        print(f"  Quantum freq contribution: {g_res_result['a_quantum_freq']:.3e}")
        print(f"  Aether freq contribution: {g_res_result['a_Aether_freq']:.3e}")
        print(f"  Total resonance modes: {len(g_res_result['resonance_modes'])}")
        
        # Test 4: AT2024tvd TDE parameters
        print("\n--- AT2024tvd TDE Event ---")
        tde = db.get_system('AT2024tvd')
        print(f"✓ AT2024tvd retrieved:")
        print(f"  L_bol (peak) = {tde.L_bol:.3e} erg/s")
        print(f"  L_X = {tde.L_X:.3e} erg/s")
        print(f"  T_BB = {tde.T:.3e} K")
        print(f"  Peak time = {tde.additional_params['peak_time']:.1f} days")
        print(f"  Recurrence = {tde.additional_params['recurrence_time']:.1f} days")
        
        return True
    except Exception as e:
        print(f"✗ FAILED: {e}")
        traceback.print_exc()
        return False


def test_system_summary_output():
    """Test 6: System summary formatting."""
    print("\n" + "="*70)
    print("TEST 6: System Summary Output")
    print("="*70)
    
    try:
        from UQFFSystemsDatabase import UQFFSystemsDatabase
        
        db = UQFFSystemsDatabase()
        
        # Test summary for WR 124
        summary = db.get_system_summary('WR124')
        print(summary)
        
        assert 'WR 124' in summary, "WR 124 not in summary"
        assert 'Wolf-Rayet' in summary, "Category not in summary"
        print("✓ Summary format validated")
        
        return True
    except Exception as e:
        print(f"✗ FAILED: {e}")
        traceback.print_exc()
        return False


def main():
    """Run all Priority 2 integration tests."""
    print("\n" + "="*70)
    print("PRIORITY 2 INTEGRATION TEST SUITE")
    print("UQFFSystemsDatabase + CondensedPhysics2 Validation")
    print("="*70)
    
    tests = [
        ("Database Initialization", test_database_initialization),
        ("System Access Methods", test_system_access_methods),
        ("Category Filtering", test_category_filtering),
        ("CondensedPhysics2 Import", test_condensedphysics2_import),
        ("Combined UQFF+Database", test_combined_uqff_database),
        ("System Summary Output", test_system_summary_output),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"\n✗ CRITICAL ERROR in {test_name}: {e}")
            traceback.print_exc()
            results.append((test_name, False))
    
    # Final summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    print(f"\n{passed}/{total} tests passed ({100*passed/total:.1f}%)")
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED! Priority 2 integration complete.")
        return 0
    else:
        print(f"\n⚠️  {total - passed} test(s) failed. Review output above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
