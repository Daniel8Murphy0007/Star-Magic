#!/usr/bin/env python3
"""
extraction_demo.py - Demonstration of UQFF Extraction Layer
============================================================

Shows how to use the complete pipeline from user queries to UQFF calculations.

Examples:
    1. Single object query
    2. Batch processing
    3. Quick formatted output
    4. Custom parameter requirements
"""

from ExtractionLayer import compute_for_object, compute_batch, quick_query
from QCalc import UQFFScale

def demo_single_object():
    """Demonstrate single object processing."""
    print("\n" + "="*80)
    print("DEMO 1: Single Object Processing")
    print("="*80 + "\n")
    
    # Query a well-known object
    result = compute_for_object(
       "Betelgeuse",
        required_params=['M', 'r', 'T'],
        scale=UQFFScale.STELLAR,
        export_csv=True,
        verbose=True
    )
    
    # Access results
    print("\n" + "-"*80)
    print("Result Summary:")
    print("-"*80)
    print(f"Object: {result['object_name']}")
    print(f"Sources: {', '.join(result['sources'])}")
    print(f"Fetch Time: {result['fetch_time_seconds']:.2f}s")
    print(f"Equations Computed: {len(result['equations'])}")
    print(f"CSV Path: {result['csv_path']}")
    
    # Show key solutions
    print("\nKey Solutions:")
    for name in ['Ug', 'UQFF_Triadic', 'UQFF_Compressed']:
        if name in result['solutions']:
            print(f"  {name}: {result['solutions'][name]:.4e}")


def demo_batch_processing():
    """Demonstrate batch processing of multiple objects."""
    print("\n" + "="*80)
    print("DEMO 2: Batch Processing")
    print("="*80 + "\n")
    
    objects = [
        "Sirius",
        "Vega",
        "Arcturus"
    ]
    
    results = compute_batch(
        objects,
        required_params=['M', 'r', 'T'],
        scale=UQFFScale.STELLAR,
        export_csv=True,
        delay_seconds=1.0
    )
    
    # Summary
    print("\n" + "-"*80)
    print("Batch Results Summary:")
    print("-"*80)
    
    for idx, result in enumerate(results, 1):
        if 'error' in result:
            print(f"{idx}. {result['object_name']}: ✗ ERROR - {result['error']}")
        else:
            num_eq = len(result['equations'])
            print(f"{idx}. {result['object_name']}: ✓ {num_eq} equations computed")


def demo_quick_query():
    """Demonstrate quick query with formatted output."""
    print("\n" + "="*80)
    print("DEMO 3: Quick Formatted Query")
    print("="*80 + "\n")
    
    quick_query("Rigel")


def demo_custom_requirements():
    """Demonstrate processing with custom parameter requirements."""
    print("\n" + "="*80)
    print("DEMO 4: Custom Parameter Requirements")
    print("="*80 + "\n")
    
    result = compute_for_object(
        "Crab Nebula",
        required_params=['M', 'r', 'T', 'B', 'omega'],  # Include magnetic field and rotation
        scale=UQFFScale.GALACTIC,
        export_csv=True,
        verbose=True
    )
    
    print("\n" + "-"*80)
    print("Parameters Retrieved:")
    print("-"*80)
    
    params = result['input_parameters']
    for param, value in params.items():
        if param not in ['query_id', 'query_name', 'timestamp', 'sources'] and value is not None:
            if isinstance(value, float) and abs(value) > 1e6:
                print(f"  {param}: {value:.4e}")
            else:
                print(f"  {param}: {value}")


if __name__ == "__main__":
    import sys
    
    print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║                UQFF Extraction Layer - Demonstration Suite                   ║
╚══════════════════════════════════════════════════════════════════════════════╝

This demonstration shows the complete UQFF data extraction and computation pipeline:

  User Query → APIFetch (SIMBAD/NED/Grok) → IPData → QCalc → OPData → Results

Running demos:
  1. Single object processing
  2. Batch processing (multiple objects)
  3. Quick formatted query
  4. Custom parameter requirements

Press Ctrl+C to skip any demo.
""")
    
    demos = [
        ("Single Object", demo_single_object),
        ("Batch Processing", demo_batch_processing),
        ("Quick Query", demo_quick_query),
        ("Custom Requirements", demo_custom_requirements)
    ]
    
    for idx, (name, demo_func) in enumerate(demos, 1):
        try:
            print(f"\n\n{'#'*80}")
            print(f"Running Demo {idx}/{len(demos)}: {name}")
            print(f"{'#'*80}")
            
            demo_func()
            
            if idx < len(demos):
                input("\nPress Enter to continue to next demo (or Ctrl+C to exit)...")
                
        except KeyboardInterrupt:
            print(f"\n\nDemo {idx} interrupted. Exiting...")
            break
        except Exception as e:
            print(f"\n✗ Demo {idx} error: {e}")
            print("This may be normal if API keys are not configured.")
            print("Note: Full functionality requires NASA_API_KEY and XAI_API_KEY environment variables.")
            
            if idx < len(demos):
                try:
                    input("\nPress Enter to continue to next demo (or Ctrl+C to exit)...")
                except KeyboardInterrupt:
                    print("\nExiting...")
                    break
    
    print("\n" + "="*80)
    print("Demonstration Complete!")
    print("="*80)
    print("""
To use in your own code:

    from ExtractionLayer import compute_for_object, compute_batch, quick_query
    
    # Single object
    result = compute_for_object("M87")
    
    # Batch
    results = compute_batch(["NGC 1365", "NGC 4258", "M51"])
    
    # Quick query
    quick_query("Andromeda")

""")
