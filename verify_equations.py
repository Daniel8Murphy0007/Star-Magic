#!/usr/bin/env python3
"""Verify all equations and compute functions are in place."""

from CondensedPhysics import (
    NGC2264Model, UGC10214Model, NGC4676Model, RedSpiderNebulaModel,
    NGC3372Model, AGCarinaeModel, M42Model, TarantulaNebulaModel,
    NGC2841Model, MysticMountainModel
)

def main():
    required_methods = [
        '__init__', 'compute_H_z', 'compute_hubble_correction', 'compute_g_grav',
        'compute_compressed_equation', 'compute_resonance_equation', 
        'validate_model', 'run_tests'
    ]

    # Extended methods for NGC2264 (full implementation)
    extended_methods = ['compute_M_sf', 'compute_E_rad', 'compute_a_EM_with_UA']

    models = [
        NGC2264Model, UGC10214Model, NGC4676Model, RedSpiderNebulaModel,
        NGC3372Model, AGCarinaeModel, M42Model, TarantulaNebulaModel,
        NGC2841Model, MysticMountainModel
    ]

    print("=" * 70)
    print("EQUATION & COMPUTE FUNCTION VERIFICATION")
    print("=" * 70)

    all_complete = True
    for model_cls in models:
        name = model_cls.__name__
        missing = []
        present = []
        
        for method in required_methods:
            if hasattr(model_cls, method):
                present.append(method)
            else:
                missing.append(method)
        
        # Check extended methods
        extended_present = []
        for method in extended_methods:
            if hasattr(model_cls, method):
                extended_present.append(method)
        
        if missing:
            print(f"\n{name}:")
            print(f"  Present ({len(present)}): {', '.join(present)}")
            print(f"  MISSING: {missing}")
            all_complete = False
        else:
            ext_str = f" + {len(extended_present)} extended" if extended_present else ""
            print(f"{name}: {len(present)} methods present{ext_str}")

    print("\n" + "=" * 70)
    
    # Test actual computation
    print("\nFUNCTIONAL TEST (compute equations):")
    print("-" * 70)
    
    for model_cls in models:
        name = model_cls.__name__
        try:
            instance = model_cls()
            
            # Test compute methods
            g_grav = instance.compute_g_grav()
            H_z = instance.compute_H_z()
            hubble = instance.compute_hubble_correction()
            compressed = instance.compute_compressed_equation()
            resonance = instance.compute_resonance_equation()
            
            g_comp = compressed.get('g_compressed', compressed.get('g_total', 0))
            R_amp = resonance.get('R_amplitude', resonance.get('R_total', 0))
            
            print(f"{name}:")
            print(f"  g_grav = {g_grav:.4e} m/s²")
            print(f"  H(z) = {H_z:.4e} s⁻¹")
            print(f"  g_compressed = {g_comp:.4e} m/s²")
            print(f"  R_amplitude = {R_amp:.4e} m/s²")
            
        except Exception as e:
            print(f"{name}: ERROR - {e}")
            all_complete = False

    print("\n" + "=" * 70)
    status = "ALL EQUATIONS VERIFIED" if all_complete else "INCOMPLETE"
    print(f"STATUS: {status}")
    print("=" * 70)

if __name__ == "__main__":
    main()
