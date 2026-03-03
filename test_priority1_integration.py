"""
Test Priority 1 Integration - New Grok Thread Physics
"""

print("=" * 70)
print("TESTING PRIORITY 1 INTEGRATION (Grok Thread 98b2e77d)")
print("=" * 70)

# Test 1: BuoyancyProofVariants import
print("\n1. Testing BuoyancyProofVariants...")
try:
    from BuoyancyProofVariants import BuoyancyProofVariantsCalculator
    calc = BuoyancyProofVariantsCalculator()
    variants = calc.list_variants()
    print(f"   ✓ SUCCESS: {len(variants)} buoyancy variants loaded")
    print(f"   Variants: {', '.join(variants)}")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

# Test 2: GrokThreadUQFFExtensions new classes
print("\n2. Testing GrokThreadUQFFExtensions new classes...")
try:
    from GrokThreadUQFFExtensions import (
        UniversalMagnetismCalculator,
        AetherMetricTensor,
        UnifiedFieldCalculator
    )
    print("   ✓ UniversalMagnetismCalculator imported")
    print("   ✓ AetherMetricTensor imported")
    print("   ✓ UnifiedFieldCalculator imported")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

# Test 3: CondensedPhysics2 integration
print("\n3. Testing CondensedPhysics2 integration...")
try:
    import CondensedPhysics2 as cp2
    print("   ✓ CondensedPhysics2 imports successfully")
    
    # Check if new classes are accessible
    has_buoyancy = hasattr(cp2, 'BuoyancyProofVariantsCalculator')
    has_um = hasattr(cp2, 'UniversalMagnetismCalculator')
    has_aether = hasattr(cp2, 'AetherMetricTensor')
    has_unified = hasattr(cp2, 'UnifiedFieldCalculator')
    
    print(f"   ✓ BuoyancyProofVariantsCalculator accessible: {has_buoyancy}")
    print(f"   ✓ UniversalMagnetismCalculator accessible: {has_um}")
    print(f"   ✓ AetherMetricTensor accessible: {has_aether}")
    print(f"   ✓ UnifiedFieldCalculator accessible: {has_unified}")
    
except Exception as e:
    print(f"   ✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

# Test 4: Quick functional test
print("\n4. Running functional tests...")
try:
    from BuoyancyProofVariants import FUBiiVirialXray
    virx_calc = FUBiiVirialXray()
    F_virx = virx_calc.compute(sigma_X=1000e3, r_h=2e22, Q_wave=1.0)
    print(f"   ✓ F_UBii_virx = {F_virx:.3e} N (virial X-ray)")
    
    from GrokThreadUQFFExtensions import UniversalMagnetismCalculator
    um_calc = UniversalMagnetismCalculator()
    dipoles = [{'mu': 1e30, 'r': 1e3, 'phi_hat': 1.0}]
    Um = um_calc.compute_Um(dipoles, gamma_t=1.0, P_SCm=1e-10, E_react=5.0, f_Heav=1.0)
    print(f"   ✓ Um (LENR) = {Um:.3e} (magnetism with Heaviside)")
    
    from GrokThreadUQFFExtensions import AetherMetricTensor
    aether = AetherMetricTensor()
    metric = aether.compute_full_metric_tensor(M=1.989e30, J=1e42, r=1e9)
    print(f"   ✓ A^00 (temporal) = {metric['A00_temporal']:.6f}")
    print(f"   ✓ A^11 (spatial) = {metric['A11_spatial_x']:.6f}")
    
except Exception as e:
    print(f"   ✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("INTEGRATION TEST COMPLETE")
print("=" * 70)
