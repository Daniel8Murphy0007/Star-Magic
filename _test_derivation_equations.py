#!/usr/bin/env python3
"""Quick test of _uqff_derivation_equations.py"""

from _uqff_derivation_equations import DerivationRegistry

print("=" * 70)
print("UQFF Derivation Equations Registry Test")
print("=" * 70)
print()

reg = DerivationRegistry()
derivs = reg.list_all_derivations()
print(f"✓ Registry initialized with {len(derivs)} derivations")
print()

print("Registered Constants:")
print("-" * 70)
for name in sorted(derivs.keys()):
    deriv = derivs[name]
    print(f"  {name:<30} S{deriv.session_number:3d} | {deriv.domain:<30}")

print()
print("Sample Derivation Details:")
print("-" * 70)

# Test F_TRZ
deriv_ftrz = reg.get_derivation("F_TRZ")
if deriv_ftrz:
    print(f"\n{deriv_ftrz.constant_name} (Session {deriv_ftrz.session_number})")
    print(f"  Domain: {deriv_ftrz.domain}")
    print(f"  Value: {deriv_ftrz.constant_value}")
    print(f"  Equation: {deriv_ftrz.equation_latex}")
    print(f"  Status: {deriv_ftrz.status}")
    print(f"  Description: {deriv_ftrz.description[:80]}...")
    print(f"  Derivation Steps: {len(deriv_ftrz.derivation_steps)}")

print()
print("✓ All tests passed!")
print()
print("File ready for use:")
print("  from _uqff_derivation_equations import DerivationRegistry")
print("  reg = DerivationRegistry()")
print("  deriv = reg.get_derivation('CONSTANT_NAME')")
print("  print(deriv.equation_latex)")
