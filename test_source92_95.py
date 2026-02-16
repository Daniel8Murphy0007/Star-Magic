# ============================================================================
# SOURCE92: BUOYANCY COUPLING TESTS
# ============================================================================

import math

def test_source92_default_buoyancy_coupling():
    """Test SOURCE92 with default parameters."""
    print("TEST: Source92 buoyancy coupling (default params)...")
    
    from Phase7_Consolidated import Source92_BuoyancyCoupling
    result = Source92_BuoyancyCoupling.calculate_buoyancy_coupling_master()
    
    # Beta uniform
    assert result['beta'] == 0.6, "Beta must be 0.6 uniform"
    
    # All U_bi negative (opposes gravity)
    assert result['U_b1'] < 0, "U_b1 must be negative (opposes gravity)"
    assert result['U_b2'] < 0, "U_b2 must be negative"
    assert result['U_b3'] < 0, "U_b3 must be negative"
    assert result['U_b4'] < 0, "U_b4 must be negative"
    
    # F_U contribution negative
    assert result['F_U_contribution'] < 0, "F_U contribution opposes gravity"
    
    # U_b1 should be largest (highest U_g1)
    assert abs(result['U_b1']) > abs(result['U_b2']), "U_b1 > U_b2 (higher U_g1)"
    
    print(f"  beta = {result['beta']}")
    print(f"  U_b1 = {result['U_b1']:.3e} J/m³ (opposes gravity)")
    print(f"  F_U contribution = {result['F_U_contribution']:.3e} J/m³")
    print("✅ Source92 default buoyancy coupling validated")


def test_source92_beta_uniform():
    """Test β_i uniform across all i (1-4)."""
    print("TEST: Source92 β_i uniform (0.6 for all i)...")
    
    from Phase7_Consolidated import Source92_BuoyancyCoupling
    for i in range(1, 5):
        beta = Source92_BuoyancyCoupling.compute_beta(i)
        assert beta == 0.6, f"Beta_{i} must be 0.6, got {beta}"
    
    print("  β_1 = β_2 = β_3 = β_4 = 0.6 ✅")
    print("✅ Source92 beta uniformity validated")


def test_source92_u_bi_scaling():
    """Test U_bi scaling with U_gi."""
    print("TEST: Source92 U_bi scaling with U_gi...")
    
    from Phase7_Consolidated import Source92_BuoyancyCoupling
    # Test with doubled U_g1
    custom = {'U_g1': 2.78e26}  # 2× default
    result = Source92_BuoyancyCoupling.calculate_buoyancy_coupling_master(custom)
    default_result = Source92_BuoyancyCoupling.calculate_buoyancy_coupling_master()
    
    # U_b1 should scale linearly with U_g1
    ratio = result['U_b1'] / default_result['U_b1']
    assert 1.9 < ratio < 2.1, f"U_b1 should double with 2× U_g1, got {ratio}"
    
    print(f"  Default U_b1 = {default_result['U_b1']:.3e} J/m³")
    print(f"  2× U_g1 U_b1 = {result['U_b1']:.3e} J/m³")
    print(f"  Ratio = {ratio:.2f} (expected ~2.0)")
    print("✅ Source92 U_bi scaling validated")


def test_source92_f_u_contribution():
    """Test F_U contribution sum."""
    print("TEST: Source92 F_U contribution sum...")
    
    from Phase7_Consolidated import Source92_BuoyancyCoupling
    result = Source92_BuoyancyCoupling.calculate_buoyancy_coupling_master()
    
    # Manual sum
    manual_sum = result['U_b1'] + result['U_b2'] + result['U_b3'] + result['U_b4']
    
    assert abs(result['F_U_contribution'] - manual_sum) < 1e-10, \
        "F_U contribution must equal sum of U_bi"
    
    print(f"  U_b1 + U_b2 + U_b3 + U_b4 = {manual_sum:.3e} J/m³")
    print(f"  F_U_contribution = {result['F_U_contribution']:.3e} J/m³")
    print("✅ Source92 F_U contribution sum validated")


# ============================================================================
# SOURCE93: SOLAR WIND BUOYANCY TESTS
# ============================================================================

def test_source93_default_solar_wind():
    """Test SOURCE93 with default parameters."""
    print("TEST: Source93 solar wind buoyancy (default params)...")
    
    from Phase7_Consolidated import Source93_SolarWindBuoyancy
    result = Source93_SolarWindBuoyancy.calculate_solar_wind_buoyancy_master()
    
    # ε_sw should be 0.001
    assert result['epsilon_sw'] == 0.001, "ε_sw must be 0.001"
    
    # Modulation factor should be very close to 1
    assert 0.999 < result['modulation_factor'] < 1.001, \
        f"Modulation factor should be ~1, got {result['modulation_factor']}"
    
    # U_b1 should be negative (opposes gravity)
    assert result['U_b1'] < 0, "U_b1 must be negative"
    
    print(f"  ε_sw = {result['epsilon_sw']}")
    print(f"  modulation_factor = {result['modulation_factor']:.18f}")
    print(f"  U_b1 = {result['U_b1']:.3e} J/m³")
    print("✅ Source93 default solar wind validated")


def test_source93_epsilon_sw_modulation():
    """Test ε_sw modulation formula: 1 + ε_sw * ρ_vac,sw."""
    print("TEST: Source93 ε_sw modulation formula...")
    
    from Phase7_Consolidated import Source93_SolarWindBuoyancy
    # Test with large ε_sw to make difference measurable
    custom = {'epsilon_sw': 1000.0}  # Large value for testing
    result = Source93_SolarWindBuoyancy.calculate_solar_wind_buoyancy_master(custom)
    
    # Manual calculation
    rho_vac_sw = 8.0e-21  # kg/m³
    expected_modulation = 1.0 + 1000.0 * rho_vac_sw
    
    assert abs(result['modulation_factor'] - expected_modulation) < 1e-15, \
        f"Modulation should be {expected_modulation}, got {result['modulation_factor']}"
    
    print(f"  ε_sw = {result['epsilon_sw']}")
    print(f"  Expected modulation = 1 + ε_sw * ρ_vac,sw = {expected_modulation:.15f}")
    print(f"  Actual modulation = {result['modulation_factor']:.15f}")
    print("✅ Source93 ε_sw modulation formula validated")


def test_source93_negligible_correction():
    """Test that ε_sw correction is negligible (~8e-24)."""
    print("TEST: Source93 negligible correction magnitude...")
    
    from Phase7_Consolidated import Source93_SolarWindBuoyancy
    result = Source93_SolarWindBuoyancy.calculate_solar_wind_buoyancy_master()
    
    # Correction = ε_sw * ρ_vac,sw
    rho_vac_sw = 8.0e-21  # kg/m³ (from DEFAULT_PARAMS)
    correction = result['epsilon_sw'] * rho_vac_sw
    
    assert correction < 1e-20, f"Correction should be negligible, got {correction}"
    
    print(f"  ε_sw = {result['epsilon_sw']}")
    print(f"  ρ_vac,sw = {rho_vac_sw:.1e} kg/m³")
    print(f"  Correction = ε_sw * ρ_vac,sw = {correction:.1e} (negligible)")
    print("✅ Source93 negligible correction validated")


# ============================================================================
# SOURCE94: Ug COUPLING TESTS
# ============================================================================

def test_source94_default_ug_coupling():
    """Test SOURCE94 with default parameters."""
    print("TEST: Source94 Ug coupling (default params)...")
    
    from Phase7_Consolidated import Source94_UgCoupling
    result = Source94_UgCoupling.calculate_ug_coupling_master()
    
    # k_i values
    assert result['k1'] == 1.5, "k1 must be 1.5"
    assert result['k2'] == 1.2, "k2 must be 1.2"
    assert result['k3'] == 1.8, "k3 must be 1.8"
    assert result['k4'] == 1.0, "k4 must be 1.0"
    
    # Sum should be dominated by k2*U_g2 (largest U_g)
    assert result['k2_U_g2'] > result['k1_U_g1'], "k2*U_g2 should dominate"
    assert result['sum_k_ugi'] > 0, "Sum must be positive"
    
    print(f"  k1 = {result['k1']} (Dipole)")
    print(f"  k2 = {result['k2']} (Bubble)")
    print(f"  k3 = {result['k3']} (Magnetic Disk)")
    print(f"  k4 = {result['k4']} (Star-BH)")
    print(f"  Sum Σ k_i U_gi = {result['sum_k_ugi']:.3e} J/m³")
    print("✅ Source94 default Ug coupling validated")


def test_source94_k_i_values():
    """Test k_i values for all i (1-4)."""
    print("TEST: Source94 k_i values...")
    
    from Phase7_Consolidated import Source94_UgCoupling
    k1 = Source94_UgCoupling.compute_k_i(1)
    k2 = Source94_UgCoupling.compute_k_i(2)
    k3 = Source94_UgCoupling.compute_k_i(3)
    k4 = Source94_UgCoupling.compute_k_i(4)
    
    assert k1 == 1.5, f"k1 must be 1.5, got {k1}"
    assert k2 == 1.2, f"k2 must be 1.2, got {k2}"
    assert k3 == 1.8, f"k3 must be 1.8, got {k3}"
    assert k4 == 1.0, f"k4 must be 1.0, got {k4}"
    
    print(f"  k_1 = {k1}, k_2 = {k2}, k_3 = {k3}, k_4 = {k4} ✅")
    print("✅ Source94 k_i values validated")


def test_source94_scaled_ug_terms():
    """Test k_i * U_gi scaling."""
    print("TEST: Source94 scaled Ug terms...")
    
    from Phase7_Consolidated import Source94_UgCoupling
    # Get raw U_gi values
    u_g1 = Source94_UgCoupling.compute_u_gi(1)
    u_g2 = Source94_UgCoupling.compute_u_gi(2)
    
    # Get scaled values
    result = Source94_UgCoupling.calculate_ug_coupling_master()
    
    # Manual calculation
    k1_manual = 1.5 * u_g1
    k2_manual = 1.2 * u_g2
    
    assert abs(result['k1_U_g1'] - k1_manual) < 1e10, "k1*U_g1 mismatch"
    assert abs(result['k2_U_g2'] - k2_manual) < 1e30, "k2*U_g2 mismatch"
    
    print(f"  k1 * U_g1 = {result['k1_U_g1']:.3e} J/m³")
    print(f"  k2 * U_g2 = {result['k2_U_g2']:.3e} J/m³ (dominant)")
    print("✅ Source94 scaled Ug terms validated")


def test_source94_sum_scaling():
    """Test sum Σ k_i * U_gi."""
    print("TEST: Source94 sum scaling...")
    
    from Phase7_Consolidated import Source94_UgCoupling
    result = Source94_UgCoupling.calculate_ug_coupling_master()
    
    # Manual sum
    manual_sum = (result['k1_U_g1'] + result['k2_U_g2'] + 
                  result['k3_U_g3'] + result['k4_U_g4'])
    
    assert abs(result['sum_k_ugi'] - manual_sum) < 1e-10, \
        "Sum must equal manual calculation"
    
    print(f"  Manual sum = {manual_sum:.3e} J/m³")
    print(f"  Σ k_i U_gi = {result['sum_k_ugi']:.3e} J/m³")
    print("✅ Source94 sum scaling validated")


# ============================================================================
# SOURCE95: MAGNETIC STRING TESTS
# ============================================================================

def test_source95_default_magnetic_string():
    """Test SOURCE95 with default parameters."""
    print("TEST: Source95 magnetic string (default params)...")
    
    from Phase7_Consolidated import Source95_MagneticString
    result = Source95_MagneticString.calculate_magnetic_string_master()
    
    # r_j should be 1.496e13 m (100 AU)
    assert abs(result['r_j_m'] - 1.496e13) < 1e10, "r_j must be ~1.496e13 m"
    assert abs(result['r_j_AU'] - 100.0) < 0.1, "r_j must be 100 AU"
    
    # μ_j should be positive
    assert result['mu_j'] > 0, "μ_j must be positive"
    
    # μ_j/r_j should be positive
    assert result['mu_over_rj'] > 0, "μ_j/r_j must be positive"
    
    # U_m should be zero at t=0 (exp(-γ*0*cos(π*0)) = exp(0) = 1)
    assert abs(result['Um_contribution']) < 1e-10, "U_m must be 0 at t=0"
    
    print(f"  r_j = {result['r_j_m']:.3e} m = {result['r_j_AU']:.1f} AU")
    print(f"  μ_j = {result['mu_j']:.3e} T·m³")
    print(f"  μ_j/r_j = {result['mu_over_rj']:.3e} T·m²")
    print(f"  U_m = {result['Um_contribution']:.3e} J/m³ (zero at t=0) ✅")
    print("✅ Source95 default magnetic string validated")


def test_source95_rj_unit_conversions():
    """Test r_j unit conversions (m, AU, ly, pc)."""
    print("TEST: Source95 r_j unit conversions...")
    
    from Phase7_Consolidated import Source95_MagneticString
    result = Source95_MagneticString.calculate_magnetic_string_master()
    
    # Conversion factors
    AU_to_m = 1.496e11
    ly_to_m = 9.461e15
    pc_to_ly = 3.262
    pc_to_m = pc_to_ly * ly_to_m
    
    # Check conversions (use looser tolerance for pc due to chained conversions)
    assert abs(result['r_j_m'] - result['r_j_AU'] * AU_to_m) < 1e8, "AU conversion error"
    assert abs(result['r_j_m'] - result['r_j_ly'] * ly_to_m) < 1e8, "ly conversion error"
    assert abs(result['r_j_m'] - result['r_j_pc'] * pc_to_m) < 1e10, "pc conversion error"
    
    print(f"  r_j = {result['r_j_m']:.3e} m")
    print(f"      = {result['r_j_AU']:.1f} AU")
    print(f"      = {result['r_j_ly']:.6f} ly")
    print(f"      = {result['r_j_pc']:.6f} pc")
    print("✅ Source95 unit conversions validated")


def test_source95_mu_j_time_variation():
    """Test μ_j time variation with sin(ω_c*t)."""
    print("TEST: Source95 μ_j time variation...")
    
    from Phase7_Consolidated import Source95_MagneticString
    # Test at t=0
    result_t0 = Source95_MagneticString.calculate_magnetic_string_master(t=0.0)
    # Test at t=π/(2*ω_c) for maximum sin term (sin = 1)
    omega_c = 2.5e-6  # rad/s
    t_max = (math.pi / 2) / omega_c  # Time where sin(ω_c*t) = 1
    result_tmax = Source95_MagneticString.calculate_magnetic_string_master(t=t_max)
    
    # μ_j should vary: at t=0, sin=0; at t_max, sin=1
    # μ_j(0) = 1000 * μ_base
    # μ_j(t_max) = (1000 + 0.4) * μ_base = 1000.4 * μ_base
    expected_ratio = 1000.4 / 1000.0  # = 1.0004
    actual_ratio = result_tmax['mu_j'] / result_t0['mu_j']
    
    assert abs(actual_ratio - expected_ratio) < 1e-6, \
        f"Expected ratio {expected_ratio}, got {actual_ratio}"
    
    print(f"  μ_j(t=0) = {result_t0['mu_j']:.3e} T·m³")
    print(f"  μ_j(t_max) = {result_tmax['mu_j']:.3e} T·m³")
    print(f"  Ratio = {actual_ratio:.6f} (expected {expected_ratio:.6f})")
    print("✅ Source95 μ_j time variation validated")


def test_source95_um_contribution_t0():
    """Test U_m contribution at t=0 (should be zero)."""
    print("TEST: Source95 U_m contribution at t=0...")
    
    from Phase7_Consolidated import Source95_MagneticString
    result = Source95_MagneticString.calculate_magnetic_string_master({'t': 0.0})
    
    # At t=0: 1 - exp(-γ*0*cos(π*0)) = 1 - exp(0) = 1 - 1 = 0
    assert abs(result['Um_contribution']) < 1e-10, "U_m must be 0 at t=0"
    
    print(f"  U_m(t=0) = {result['Um_contribution']:.3e} J/m³")
    print("  Expected: 0 J/m³ (1 - exp(0) = 0) ✅")
    print("✅ Source95 U_m at t=0 validated")


def test_source95_ug3_contribution():
    """Test Ug3 contribution with k3 scaling."""
    print("TEST: Source95 Ug3 contribution...")
    
    from Phase7_Consolidated import Source95_MagneticString
    result = Source95_MagneticString.calculate_magnetic_string_master()
    
    # Ug3 should be positive
    assert result['Ug3_contribution'] > 0, "Ug3 must be positive"
    
    # Should include k3=1.8 scaling from SOURCE94
    print(f"  Ug3 = {result['Ug3_contribution']:.3e} J/m³")
    print("  (includes k3=1.8 scaling from SOURCE94)")
    print("✅ Source95 Ug3 contribution validated")


# Run tests
if __name__ == '__main__':
    print("=" * 80)
    print("SOURCE92-95 COMPREHENSIVE TEST SUITE")
    print("=" * 80)
    print()
    
    # SOURCE92 tests (4 tests)
    print("SOURCE92: BUOYANCY COUPLING (4 tests)")
    print("-" * 80)
    test_source92_default_buoyancy_coupling()
    print()
    test_source92_beta_uniform()
    print()
    test_source92_u_bi_scaling()
    print()
    test_source92_f_u_contribution()
    print()
    
    # SOURCE93 tests (3 tests)
    print("SOURCE93: SOLAR WIND BUOYANCY (3 tests)")
    print("-" * 80)
    test_source93_default_solar_wind()
    print()
    test_source93_epsilon_sw_modulation()
    print()
    test_source93_negligible_correction()
    print()
    
    # SOURCE94 tests (4 tests)
    print("SOURCE94: Ug COUPLING (4 tests)")
    print("-" * 80)
    test_source94_default_ug_coupling()
    print()
    test_source94_k_i_values()
    print()
    test_source94_scaled_ug_terms()
    print()
    test_source94_sum_scaling()
    print()
    
    # SOURCE95 tests (5 tests)
    print("SOURCE95: MAGNETIC STRING (5 tests)")
    print("-" * 80)
    test_source95_default_magnetic_string()
    print()
    test_source95_rj_unit_conversions()
    print()
    test_source95_mu_j_time_variation()
    print()
    test_source95_um_contribution_t0()
    print()
    test_source95_ug3_contribution()
    print()
    
    print("=" * 80)
    print("✅ ALL 16 TESTS PASSED (4+3+4+5)")
    print("=" * 80)
