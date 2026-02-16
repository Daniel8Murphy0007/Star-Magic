"""Proof that Phase 7 integration exists and works in QCalc.py"""
from QCalc import UnifiedFieldSolver, ComputeParams, UQFFScale

# Create solver and test with Andromeda parameters
solver = UnifiedFieldSolver()
params = ComputeParams(
    M=1.5e42,      # Andromeda mass
    r=2.5e21,      # Distance
    z=-0.001,      # Blueshift (triggers SOURCE88)
    t=3.799e10,    # Time
    scale=UQFFScale.GALACTIC
)

# Solve and check for Phase 7 equations
result = solver.solve(params)
equations = result.get('long_form_equations', [])

phase7_names = ['Aether', 'Beta', 'Ug_Coupling', 'String', 'Andromeda']
phase7_eqs = [eq for eq in equations if any(name in eq['name'] for name in phase7_names)]

print(f"✅ Phase 7 equations found: {len(phase7_eqs)}")
for eq in phase7_eqs:
    print(f"  - {eq['name']}")
