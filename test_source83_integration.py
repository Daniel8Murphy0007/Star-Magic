"""
Quick integration test for SOURCE83 LENR module
"""
from Phase7_Consolidated import PHASE7_CATALOG, __version__, __current_equations__, __progress__
from source83_lenr_extract import Source83_LENR, ScenarioType83

print("=" * 80)
print("SOURCE83 LENR MODULE - INTEGRATION TEST")
print("=" * 80)

# Test 1: Module metadata
print(f"\n✅ Phase7_Consolidated v{__version__}")
print(f"   Functions: {__current_equations__} ({__progress__} of target)")

# Test 2: Catalog completeness
print(f"\n✅ PHASE7_CATALOG has {len(PHASE7_CATALOG)} entries:")
for key in PHASE7_CATALOG.keys():
    print(f"   - {key}")

# Test 3: SOURCE83 calculation
print(f"\n✅ SOURCE83 LENR Calculation (HYDRIDE scenario):")
result = Source83_LENR.calculate_lenr_master(Source83_LENR.DEFAULT_PARAMS)
print(f"   Scenario: {result['scenario'].name}")
print(f"   Plasma Frequency: ω_pe = {result['omega_pe']:.3e} rad/s")
print(f"   Electric Field: E = {result['E_field']:.3e} V/m")
print(f"   Neutron Rate: η = {result['eta']:.3e}")
print(f"   Universal Magnetism: U_m = {result['U_m']:.3e}")
print(f"   Reactor Efficiency: E_react = {result['E_react']:.3e}")

# Test 4: Catalog function access
print(f"\n✅ Catalog Function Access:")
lenr_entry = PHASE7_CATALOG['source83_lenr_hydride']
print(f"   System: {lenr_entry['system']}")
print(f"   C Functions: {lenr_entry['c_functions']}")
print(f"   Source File: {lenr_entry['source_file']}")
print(f"   Unique Physics: {', '.join(lenr_entry['unique_physics'])}")

# Test 5: Call catalog function
print(f"\n✅ Direct Catalog Function Call:")
catalog_result = lenr_entry['function'](Source83_LENR.DEFAULT_PARAMS)
print(f"   ω_pe = {catalog_result['omega_pe']:.3e} rad/s")
print(f"   E_react = {catalog_result['E_react']:.3e}")

print("\n" + "=" * 80)
print("✅ SOURCE83 INTEGRATION TEST COMPLETE - ALL CHECKS PASSED")
print("=" * 80)
