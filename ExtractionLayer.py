#!/usr/bin/env python3
"""
ExtractionLayer.py - UQFF Complete Data Extraction & Computation Pipeline
===========================================================================

Orchestrates the complete data flow from user queries to UQFF calculations:

    User Query > APIFetch > IPData > QCalc > OPData > Results

FEATURES:
    - Single query or batch processing
    - Multi-source API fallback (SIMBAD > NED > Grok)
    - Automatic CSV export (bodies_YYYYMMDD_HHMMSS.csv)
    - Progress tracking with ETA
    - Error handling with graceful degradation
    - Complete audit trail (sources, timestamps, parameters)

PHASE 6: QCalc automatically includes Phase 6 galaxy physics (M51, NGC1316,
SMBH binaries) when appropriate parameters are detected.

PHASE 7: QCalc automatically includes Phase 7 cosmological systems (SOURCE81-95:
Andromeda, SMBH M-σ, NGC346, LENR, Aether, DPM Birth, coupling constants) with
14 complete modules and 110 functions.

USAGE:
    # Single object
    >>> from ExtractionLayer import compute_for_object
    >>> result = compute_for_object("Sagittarius A*")
    
    # Batch processing
    >>> from ExtractionLayer import compute_batch
    >>> results = compute_batch(["Betelgeuse", "M87", "NGC 1365"])
    
    # With parameter requirements
    >>> result = compute_for_object("Vela Pulsar", required_params=['M', 'r', 'B'])

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import os
import time
import csv
from typing import Dict, List, Any, Tuple, Optional
from datetime import datetime
from dataclasses import asdict

# Import UQFF modules
from APIFetch import UnifiedFetcher, manual_input
from IPData import InputParameters, InputDataStore
from QCalc import UnifiedFieldSolver, ComputeParams, UQFFScale
from OPData import OutputDataStore

# ═══════════════════════════════════════════════════════════════════════════════
# EXTRACTION LAYER ORCHESTRATOR
# ═══════════════════════════════════════════════════════════════════════════════

class ExtractionPipeline:
    """
    Complete extraction and computation pipeline.
    
    Architecture:
        1. User provides object name(s)
        2. APIFetch queries SIMBAD, NED, NASA, Grok
        3. Parameters stored in IPData
        4. QCalc computes all applicable UQFF equations
        5. Results stored in OPData
        6. CSV export for source2.cpp compatibility
    """
    
    def __init__(self, verbose: bool = True):
        """
        Initialize extraction pipeline.
        
        Args:
            verbose: Print progress messages
        """
        self.fetcher = UnifiedFetcher()
        self.solver = UnifiedFieldSolver()
        self.input_store = InputDataStore()
        self.output_store = OutputDataStore()
        self.verbose = verbose
    
    def compute_for_object(
        self, 
        object_name: str,
        required_params: List[str] = None,
        scale: UQFFScale = UQFFScale.GALACTIC,
        export_csv: bool = True
    ) -> Dict[str, Any]:
        """
        Complete pipeline for a single astronomical object.
        
        Args:
            object_name: Name to query (e.g., "Sagittarius A*", "Betelgeuse")
            required_params: List of required parameter names ['M', 'r', 'T', etc.]
            scale: UQFF scale (ATOMIC, STELLAR, GALACTIC, COSMOLOGICAL)
            export_csv: Save results to CSV file
            
        Returns:
            Complete result dictionary with equations, solutions, metadata
        """
        if self.verbose:
            print(f"\n{'='*80}")
            print(f"UQFF Extraction Pipeline: {object_name}")
            print(f"{'='*80}\n")
        
        # ───────────────────────────────────────────────────────────────────────
        # STEP 1: Fetch parameters from APIs
        # ───────────────────────────────────────────────────────────────────────
        if self.verbose:
            print("[1/4] Fetching parameters from APIs...")
        
        start_time = time.time()
        
        if required_params is None:
            required_params = ['M', 'r', 'T']  # Default requirements
        
        # Fetch from all sources
        raw_data = self.fetcher.fetch(object_name, required_params)
        
        fetch_time = time.time() - start_time
        
        if self.verbose:
            sources = raw_data.get('sources', ['No data found'])
            print(f"   ✓ Data retrieved from: {', '.join(sources)}")
            print(f"   ✓ Fetch time: {fetch_time:.2f}s")
            print(f"   ✓ Parameters found: {len([k for k in raw_data.keys() if k not in ['name', 'sources']])}")
        
        # ───────────────────────────────────────────────────────────────────────
        # STEP 2: Store in IPData
        # ───────────────────────────────────────────────────────────────────────
        if self.verbose:
            print("\n[2/4] Storing parameters in IPData...")
        
        # Create InputParameters object
        input_params = InputParameters(
            query_name=object_name,
            sources=raw_data.get('sources', []),
            # Core parameters
            M=raw_data.get('mass'),
            r=raw_data.get('distance'),
            R=raw_data.get('radius'),
            T=raw_data.get('temperature'),
            L=raw_data.get('luminosity'),
            # Additional parameters
            B=raw_data.get('magnetic_field'),
            z=raw_data.get('redshift'),
            omega=raw_data.get('rotation_rate'),
            P=raw_data.get('period'),
            v_rad=raw_data.get('radial_velocity'),
            v_disp=raw_data.get('velocity_dispersion'),
            parallax=raw_data.get('parallax'),
            spectral_type=raw_data.get('spectral_type'),
            mag_V=raw_data.get('magnitude_V'),
        )
        
        # Store in IPData
        query_id = self.input_store.store(input_params)
        
        if self.verbose:
            print(f"   ✓ Stored with Query ID: {query_id}")
            print(f"   ✓ Available parameters: {len(input_params.get_available_params())}")
        
        # ───────────────────────────────────────────────────────────────────────
        # STEP 3: Compute UQFF equations with QCalc
        # ───────────────────────────────────────────────────────────────────────
        if self.verbose:
            print("\n[3/4] Computing UQFF equations with QCalc...")
        
        # Convert InputParameters to ComputeParams
        compute_params = ComputeParams(
            query_name=object_name,
            scale=scale,
            M=input_params.M,
            r=input_params.r,
            R=input_params.R,
            T=input_params.T,
            L=input_params.L,
            B=input_params.B,
            z=input_params.z,
            omega=input_params.omega,
            P=input_params.P,
            t=4.5e9 * 365.25 * 86400  # Default: 4.5 billion years
        )
        
        # Solve all applicable equations
        calc_result = self.solver.solve(compute_params)
        
        if self.verbose:
            num_equations = len(calc_result['long_form_equations'])
            num_available = len(calc_result['available_equations'])
            print(f"   ✓ Computed: {num_equations} equations")
            print(f"   ✓ Available methods: {num_available}")
            
            # Show key results
            if 'Ug' in calc_result['solutions']:
                print(f"   ✓ Unified Gravity (Ug): {calc_result['solutions']['Ug']:.4e} m/s²")
            if 'UQFF_Triadic' in calc_result['solutions']:
                print(f"   ✓ 26-Layer Triadic: {calc_result['solutions']['UQFF_Triadic']:.4e} m/s²")
        
        # ───────────────────────────────────────────────────────────────────────
        # STEP 4: Store results in OPData
        # ───────────────────────────────────────────────────────────────────────
        if self.verbose:
            print("\n[4/4] Storing results in OPData...")
        
        # Store computation results
        result_query_id = self.output_store.store(calc_result)
        
        if self.verbose:
            print(f"   ✓ Results stored with ID: {result_query_id}")
        
        # ───────────────────────────────────────────────────────────────────────
        # STEP 5: Export to CSV if requested
        # ───────────────────────────────────────────────────────────────────────
        csv_path = None
        if export_csv:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            csv_filename = f"bodies_{timestamp}.csv"
            csv_path = self._export_to_csv(raw_data, calc_result, csv_filename)
            
            if self.verbose:
                print(f"\n[5/5] CSV export complete")
                print(f"   ✓ File: {csv_path}")
        
        # Build complete result
        complete_result = {
            'query_id': query_id,
            'result_query_id': result_query_id,
            'object_name': object_name,
            'sources': raw_data.get('sources', []),
            'fetch_time_seconds': fetch_time,
            'input_parameters': input_params.to_dict(),
            'equations': calc_result['long_form_equations'],
            'solutions': calc_result['solutions'],
            'available_equations': calc_result['available_equations'],
            'csv_path': csv_path,
            'timestamp': datetime.now().isoformat()
        }
        
        if self.verbose:
            print(f"\n{'='*80}")
            print(f"✓ Pipeline complete for {object_name}")
            print(f"{'='*80}\n")
        
        return complete_result
    
    def compute_batch(
        self,
        object_names: List[str],
        required_params: List[str] = None,
        scale: UQFFScale = UQFFScale.GALACTIC,
        export_csv: bool = True,
        delay_seconds: float = 1.0
    ) -> List[Dict[str, Any]]:
        """
        Process multiple objects in batch.
        
        Args:
            object_names: List of object names to process
            required_params: Required parameter names
            scale: UQFF scale
            export_csv: Export each result to CSV
            delay_seconds: Delay between API calls (rate limiting)
            
        Returns:
            List of complete result dictionaries
        """
        results = []
        total = len(object_names)
        
        print(f"\n{'='*80}")
        print(f"UQFF Batch Processing: {total} objects")
        print(f"{'='*80}\n")
        
        for idx, obj_name in enumerate(object_names, 1):
            print(f"\nProcessing [{idx}/{total}]: {obj_name}")
            print("-" * 80)
            
            try:
                result = self.compute_for_object(
                    obj_name,
                    required_params=required_params,
                    scale=scale,
                    export_csv=export_csv
                )
                results.append(result)
                
            except Exception as e:
                print(f"✗ Error processing {obj_name}: {e}")
                results.append({
                    'object_name': obj_name,
                    'error': str(e),
                    'timestamp': datetime.now().isoformat()
                })
            
            # Rate limiting
            if idx < total:
                time.sleep(delay_seconds)
        
        # Summary
        successful = len([r for r in results if 'error' not in r])
        failed = total - successful
        
        print(f"\n{'='*80}")
        print(f"Batch Processing Complete")
        print(f"{'='*80}")
        print(f"Successful: {successful}/{total}")
        print(f"Failed: {failed}/{total}")
        print(f"{'='*80}\n")
        
        return results
    
    def _export_to_csv(
        self,
        raw_data: Dict[str, Any],
        calc_result: Dict[str, Any],
        filename: str
    ) -> str:
        """
        Export combined API data and calculation results to CSV.
        
        Format compatible with source2.cpp bodies_YYYYMMDD_HHMMSS.csv.
        
        Args:
            raw_data: Raw API fetch results
            calc_result: QCalc computation results
            filename: Output filename
            
        Returns:
            Full path to exported CSV file
        """
        # Combine all data
        combined = {
            'name': raw_data.get('name', 'unknown'),
            'timestamp': datetime.now().isoformat(),
            'sources': ';'.join(raw_data.get('sources', [])),
            **{k: v for k, v in raw_data.items() if k not in ['name', 'sources']},
            **{f'computed_{k}': v for k, v in calc_result['solutions'].items() if isinstance(v, (int, float))}
        }
        
        # Write to CSV
        with open(filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=combined.keys())
            writer.writeheader()
            writer.writerow(combined)
        
        return os.path.abspath(filename)


# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL PIPELINE INSTANCE
# ═══════════════════════════════════════════════════════════════════════════════

PIPELINE = ExtractionPipeline()


# ═══════════════════════════════════════════════════════════════════════════════
# CONVENIENCE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def compute_for_object(
    object_name: str,
    required_params: List[str] = None,
    scale: UQFFScale = UQFFScale.GALACTIC,
    export_csv: bool = True,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    One-line computation for an astronomical object.
    
    Orchestrates complete pipeline: API fetch → IPData → QCalc → OPData → CSV
    
    Args:
        object_name: Object to query (e.g., "Sagittarius A*")
        required_params: Required parameter names (default: ['M', 'r', 'T'])
        scale: UQFF scale (default: GALACTIC)
        export_csv: Export to CSV (default: True)
        verbose: Print progress messages (default: True)
    
    Returns:
        Complete result dictionary
    
    Example:
        >>> result = compute_for_object("Betelgeuse")
        >>> print(f"Computed {len(result['equations'])} equations")
        >>> print(f"Ug = {result['solutions']['Ug']:.4e} m/s²")
    """
    pipeline = ExtractionPipeline(verbose=verbose)
    return pipeline.compute_for_object(object_name, required_params, scale, export_csv)


def compute_batch(
    object_names: List[str],
    required_params: List[str] = None,
    scale: UQFFScale = UQFFScale.GALACTIC,
    export_csv: bool = True,
    delay_seconds: float = 1.0
) -> List[Dict[str, Any]]:
    """
    Batch processing for multiple astronomical objects.
    
    Args:
        object_names: List of object names
        required_params: Required parameter names
        scale: UQFF scale
        export_csv: Export each to CSV
        delay_seconds: Delay between API calls
    
    Returns:
        List of result dictionaries
    
    Example:
        >>> results = compute_batch(["Betelgeuse", "Rigel", "Sirius"])
        >>> successful = [r for r in results if 'error' not in r]
        >>> print(f"Successfully processed {len(successful)} objects")
    """
    pipeline = ExtractionPipeline(verbose=True)
    return pipeline.compute_batch(object_names, required_params, scale, export_csv, delay_seconds)


def quick_query(object_name: str) -> None:
    """
    Quick query with formatted output to console.
    
    Perfect for interactive exploration.
    
    Args:
        object_name: Object to query
    
    Example:
        >>> quick_query("M87")
    """
    result = compute_for_object(object_name, verbose=True)
    
    print(f"\n{'='*80}")
    print(f"Quick Query Results: {object_name}")
    print(f"{'='*80}\n")
    
    print(f"Sources: {', '.join(result['sources'])}")
    print(f"Query ID: {result['query_id']}")
    print(f"\nKey Solutions:")
    
    key_solutions = ['Ug', 'UQFF_Compressed', 'UQFF_Resonant', 'UQFF_Triadic', 
                     'UQFF_Superconductive', 'F_U_Bi', 'F_U_Bi_i']
    
    for sol_name in key_solutions:
        if sol_name in result['solutions']:
            val = result['solutions'][sol_name]
            if isinstance(val, (int, float)):
                print(f"  {sol_name}: {val:.4e}")
    
    print(f"\nTotal equations computed: {len(result['equations'])}")
    print(f"CSV exported to: {result['csv_path']}")
    print(f"{'='*80}\n")


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE TEST
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import sys
    
    print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║                 UQFF Extraction Layer - Interactive Demo                     ║
╚══════════════════════════════════════════════════════════════════════════════╝

Complete pipeline from user query to UQFF calculations:
  User Query → APIFetch → IPData → QCalc → OPData → Results

Available functions:
  1. compute_for_object(name)  - Single object processing
  2. compute_batch([names])    - Batch processing
  3. quick_query(name)         - Quick formatted output

""")
    
    # Demo with a well-known object
    demo_object = "Betelgeuse"
    
    print(f"Running demo query for: {demo_object}\n")
    print("-" * 80)
    
    try:
        result = compute_for_object(demo_object, verbose=True)
        
        print("\n" + "="*80)
        print("Demo complete! Try these examples:")
        print("="*80)
        print('1. from ExtractionLayer import compute_for_object')
        print('   result = compute_for_object("Sagittarius A*")')
        print()
        print('2. from ExtractionLayer import compute_batch')
        print('   results = compute_batch(["M87", "NGC 1365", "Vela Pulsar"])')
        print()
        print('3. from ExtractionLayer import quick_query')
        print('   quick_query("Crab Nebula")')
        print("="*80 + "\n")
        
    except KeyboardInterrupt:
        print("\n\nDemo interrupted by user.")
        sys.exit(0)
    except Exception as e:
        print(f"\n✗ Demo error: {e}")
        print("\nThis is normal if API keys are not configured.")
        print("Set NASA_API_KEY and XAI_API_KEY environment variables to enable full functionality.")
