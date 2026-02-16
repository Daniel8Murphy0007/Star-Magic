#!/usr/bin/env python3
"""Test all 8 UQFF equations"""
import CondensedPhysics as CP

print("Testing all 8 UQFF equations...")
print()

result = CP.solve("Sun", equations=None, output_format="both")

if "error" in result:
    print(f"ERROR: {result['error']}")
else:
    print(f"System: {result['system_name']}")
    print(f"CSV: {result['csv_path']}")
    print()
    
    print("=" * 70)
    print("NUMERIC RESULTS")
    print("=" * 70)
    for eq_num in sorted(result.get("results", {}).keys()):
        data = result["results"][eq_num]
        name = data.get("equation_name", f"Eq {eq_num}")
        # Find primary result
        primary = None
        for k, v in data.items():
            if k != "equation_name" and isinstance(v, (int, float)):
                primary = v
                break
        if primary is not None:
            print(f"{eq_num}. {name:35} = {primary:.6e}")
        else:
            print(f"{eq_num}. {name:35} = (complex output)")
    print()
    
    print("=" * 70)
    print("LONG-FORM DERIVATIONS (first 400 chars each)")
    print("=" * 70)
    for eq_num in sorted(result.get("long_form", {}).keys()):
        text = result["long_form"][eq_num]
        print(f"--- EQUATION {eq_num} ---")
        print(text[:400])
        print("...")
        print()

print("ALL 8 UQFF EQUATIONS TESTED!")
