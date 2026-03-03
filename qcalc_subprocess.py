#!/usr/bin/env python3
"""
qcalc_subprocess.py
===================
Fast subprocess wrapper for QCalc.py - Phase 0 IPC Integration

QCalc.py advantages:
- Fast import: ~1.1s (vs CondensedPhysics.py 30s+)
- Has UnifiedFieldSolver interface
- 9,149 lines (vs CondensedPhysics.py 168,494 lines)
- Production-ready 8 UQFF master equations

Usage:
    echo '{"object_name": "SGR 1745+29", "M": 2.0, "r": 1e6}' | python qcalc_subprocess.py

Input JSON format:
    {
        "object_name": "SGR 1745+29",
        "M": 2.0,          // Solar masses (optional)
        "r": 1000000.0,    // Meters (optional)
        "z": 0.001,        // Redshift (optional)
        "B": 1e13,         // Tesla (optional)
        "T": 300,          // Kelvin (optional)
        "SFR": 10          // M☉/yr (optional)
    }

Output JSON format:
    {
        "success": true,
        "query_id": "UQFF_2026-03-02",
        "timestamp": "2026-03-02T12:34:56",
        "input_params": {...},
        "equations": ["g = G·M/r² = 1.334e-08 m/s²", ...],
        "solutions": {"g_newton": 1.334e-08, ...},
        "available_equations": ["g_newton", "g_triadic", ...],
        "compute_time_ms": 123.45
    }

Author: Phase 0 IPC Unification
Date: March 2, 2026
"""

import sys
import json
import time

def process_calculation(input_data):
    """
    Process a single QCalc calculation request
    
    Args:
        input_data: Dictionary with query parameters
    
    Returns:
        Dictionary with computation results
    """
    try:
        # Import QCalc (fast: ~1.1s)
        from QCalc import UnifiedFieldSolver, ComputeParams
    except ImportError as e:
        return {
            "success": False,
            "error": f"QCalc module not available: {str(e)}",
            "solutions": {},
            "equations": [],
            "compute_time_ms": 0
        }
    
    calc_start = time.time()
    
    try:
        # Convert dict to ComputeParams object
        params = ComputeParams()
        params.query_name = input_data.get('object_name', 'UNKNOWN')
        params.M = input_data.get('M')
        params.r = input_data.get('r')
        params.z = input_data.get('z')
        params.B = input_data.get('B')
        params.T = input_data.get('T')
        params.SFR = input_data.get('SFR')
        params.L = input_data.get('L')
        params.v = input_data.get('v')
        params.P = input_data.get('P')
        params.R = input_data.get('R')
        params.d = input_data.get('d')
        params.omega = input_data.get('omega')
        
        # Instantiate solver
        solver = UnifiedFieldSolver()
        
        # Call solve() with ComputeParams object
        result = solver.solve(params)
        
        # Ensure required fields exist
        if "success" not in result:
            result["success"] = True
        if "compute_time_ms" not in result:
            result["compute_time_ms"] = (time.time() - calc_start) * 1000
        
        return result
        
    except Exception as e:
        import traceback
        return {
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc(),
            "solutions": {},
            "equations": [],
            "compute_time_ms": (time.time() - calc_start) * 1000
        }

def main():
    """Main subprocess entry point - reads JSON from stdin, writes JSON to stdout"""
    overall_start = time.time()
    
    try:
        # Log to stderr (won't interfere with stdout JSON)
        print("[QCalc Subprocess] Starting...", file=sys.stderr, flush=True)
        
        # Read input from stdin
        input_json = sys.stdin.read()
        
        if not input_json.strip():
            error_result = {
                "success": False,
                "error": "No input data received on stdin",
                "solutions": {},
                "equations": [],
                "compute_time_ms": 0
            }
            print(json.dumps(error_result), flush=True)
            return 1
        
        #Parse JSON
        try:
            input_data = json.loads(input_json)
        except json.JSONDecodeError as e:
            error_result = {
                "success": False,
                "error": f"Invalid JSON input: {str(e)}",
                "solutions": {},
                "equations": [],
                "compute_time_ms": 0
            }
            print(json.dumps(error_result), flush=True)
            return 1
        
        # Log request
        object_name = input_data.get('object_name', 'UNKNOWN')
        print(f"[QCalc Subprocess] Processing: {object_name}", file=sys.stderr, flush=True)
        
        # Process calculation
        result = process_calculation(input_data)
        
        # Add overall timing
        overall_ms = (time.time() - overall_start) * 1000
        result["total_time_ms"] = overall_ms
        
        # Log completion
        if result.get("success", False):
            print(f"[QCalc Subprocess] ✓ Completed in {overall_ms:.1f}ms", 
                  file=sys.stderr, flush=True)
        else:
            print(f"[QCalc Subprocess] ✗ Failed: {result.get('error', 'Unknown error')}", 
                  file=sys.stderr, flush=True)
        
        # Write result to stdout (with custom JSON encoder for enums)
        def json_encoder(obj):
            """Custom JSON encoder to handle enums and other non-serializable types"""
            if hasattr(obj, '__name__'):  # Handle enums
                return str(obj)
            elif hasattr(obj, 'value'):  # Handle enum values
                return obj.value
            elif hasattr(obj, '__dict__'):  # Handle dataclasses/objects
                return str(obj)
            else:
                return str(obj)
        
        output_json = json.dumps(result, default=json_encoder)
        print(output_json, flush=True)
        
        # Return exit code based on success
        return 0 if result.get("success", False) else 1
        
    except Exception as e:
        import traceback
        error_result = {
            "success": False,
            "error": f"Subprocess error: {str(e)}",
            "traceback": traceback.format_exc(),
            "solutions": {},
            "equations": [],
            "compute_time_ms": 0,
            "total_time_ms": (time.time() - overall_start) * 1000
        }
        print(json.dumps(error_result), flush=True)
        print(f"[QCalc Subprocess] ✗ Exception: {e}", file=sys.stderr, flush=True)
        return 1

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
