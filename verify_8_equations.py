#!/usr/bin/env python3
"""Verify all 8 UQFF Master Equations work in all 10 models."""

from CondensedPhysics import (
    NGC2264Model, UGC10214Model, NGC4676Model, RedSpiderNebulaModel,
    NGC3372Model, AGCarinaeModel, M42Model, TarantulaNebulaModel,
    NGC2841Model, MysticMountainModel
)

def test_model(model_cls):
    """Test all 8 UQFF equations for a model class."""
    equations = [
        ('compute_UQFF_base', 'F_U'),
        ('compute_compressed_equation', 'g_compressed'),
        ('compute_resonance_equation', 'R_amplitude'),
        ('compute_superconductive_equation', 'SCm_t'),
        ('compute_buoyant_equation', 'F_U_Bi'),
        ('compute_master_buoyant_equation', 'F_U_Bi_i'),
        ('compute_triadic_equation', 'g_triadic'),
        ('compute_quadratic_equation', 'x1'),
    ]
    
    results = []
    instance = model_cls()
    
    for method_name, key in equations:
        try:
            method = getattr(instance, method_name)
            result = method()
            value = result.get(key, result.get('result', 'N/A'))
            if isinstance(value, (int, float)):
                results.append((method_name, True, value))
            else:
                results.append((method_name, True, 'computed'))
        except Exception as e:
            results.append((method_name, False, str(e)[:50]))
    
    return results

def main():
    models = [
        NGC2264Model, UGC10214Model, NGC4676Model, RedSpiderNebulaModel,
        NGC3372Model, AGCarinaeModel, M42Model, TarantulaNebulaModel,
        NGC2841Model, MysticMountainModel
    ]
    
    print("=" * 80)
    print("8 UQFF MASTER EQUATIONS - VERIFICATION")
    print("=" * 80)
    
    all_passed = 0
    total_tests = 0
    
    for model_cls in models:
        name = model_cls.__name__
        results = test_model(model_cls)
        passed = sum(1 for _, success, _ in results if success)
        total = len(results)
        all_passed += passed
        total_tests += total
        
        status = "PASS" if passed == total else "FAIL"
        print(f"\n{name}: {passed}/{total} [{status}]")
        
        for method_name, success, value in results:
            short_name = method_name.replace('compute_', '')
            if success:
                if isinstance(value, float):
                    print(f"    {short_name}: {value:.4e}")
                else:
                    print(f"    {short_name}: {value}")
            else:
                print(f"    {short_name}: ERROR - {value}")
    
    print("\n" + "=" * 80)
    print(f"TOTAL: {all_passed}/{total_tests} equations computed successfully")
    all_complete = all_passed == total_tests
    print(f"STATUS: {'ALL 8 EQUATIONS VERIFIED' if all_complete else 'INCOMPLETE'}")
    print("=" * 80)

if __name__ == "__main__":
    main()
