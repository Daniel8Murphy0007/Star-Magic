"""Test script for NGC3603StarClusterModel - Original 10-Term Document Equation."""

from CondensedPhysics import NGC3603StarClusterModel

# Run tests
result = NGC3603StarClusterModel.run_tests()
print(f"NGC3603StarClusterModel: {result['passed']}/{result['total']} tests passed")
print()
for t in result['tests']:
    status = 'PASS' if t['passed'] else 'FAIL'
    error = f" - {t.get('error', '')}" if 'error' in t and not t['passed'] else ''
    print(f"  {status}: {t['name']}{error}")

print()

# Test the master gravity equation with original 10 terms
model = NGC3603StarClusterModel()
t = 0.5e6 * 3.156e7  # 0.5 Myr

print("=== ORIGINAL 10-TERM EQUATION (May 2025 Document) ===")
print()

# Individual non-extension terms
g1, _ = model.compute_Ug1(t)
print(f"Ug1 (base gravity):        {g1:.4e} m/s2")

g4, _ = model.compute_Ug4(t)
print(f"Ug4 (superconductive):     {g4:.4e} m/s2")

g_uqff, _ = model.compute_Ug_total_with_TRZ(t, include_dynamics=False)
print(f"UQFF total with TRZ:       {g_uqff:.4e} m/s2 (Ug2=Ug3=0 per document)")

g_em, _ = model.compute_electromagnetic_with_UA()
print(f"EM with UA:                {g_em:.4e} m/s2")

g_dm, _ = model.compute_dark_matter_perturbation()
print(f"Dark Matter:               {g_dm:.4e} m/s2")

g_wind, _ = model.compute_stellar_wind_feedback()
print(f"Stellar Wind:              {g_wind:.4e} m/s2")

g_buoy, _ = model.compute_fluid_buoyancy()
print(f"Fluid Buoyancy:            {g_buoy:.4e} m/s2")

print()
g_total, _ = model.compute_master_gravity(t, include_all_terms=True, include_extensions=False)
print(f"TOTAL g_NGC3603 (10 terms): {g_total:.4e} m/s2")
print(f"DOCUMENT REFERENCE:         ~1.053e-3 m/s2")

print()
print("=== WITH RESEARCH EXTENSIONS (Terms 11-14) ===")
g_ext, _ = model.compute_master_gravity(t, include_all_terms=True, include_extensions=True)
print(f"TOTAL with extensions:      {g_ext:.4e} m/s2")
print("(Extensions are NOT part of original document equation)")
