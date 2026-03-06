#!/usr/bin/env python3
"""
condensed_physics_subprocess_FAST.py
=====================================
Optimized subprocess wrapper for CondensedPhysics with lazy imports.

This version delays the slow CondensedPhysics import until after confirming
stdin is available, avoiding the 30s+ import overhead for every subprocess call.

Usage:
    echo '{"object_name": "SGR 1745+29", "M": 2.0}' | python condensed_physics_subprocess_FAST.py
"""

import sys
import json
import time

def process_calculation(input_data):
    """
    Process a single calculation request
    
    Args:
        input_data: Dictionary with keys: object_name, M, r, z, B, T, SFR (all optional except object_name)
    
    Returns:
        Dictionary with computation results
    """
    # LAZY IMPORT: Only import CondensedPhysics when actually needed
    # This allows the subprocess to start quickly and fail fast if stdin is unavailable
    try:
        from CondensedPhysics import UnifiedFieldSolver
    except ImportError as e:
        return {
            "success": False,
            "error": f"CondensedPhysics module not available: {str(e)}",
            "solutions": [],
            "equations": [],
            "compute_time_ms": 0
        }

    # Import CondensedPhysics2 calculators via Aggregator (CP2 wiring, Problem 4)
    try:
        from CondensedPhysicsAggregator import CondensedPhysicsAggregator
        CP2_AVAILABLE = True
    except ImportError:
        CP2_AVAILABLE = False
    
    start_time = time.time()
    
    try:
        # Instantiate solver
        solver = UnifiedFieldSolver()
        
        # Call solve() with input parameters
        result = solver.solve(input_data)
        
        # Add success flag and timing
        result["success"] = True
        result["compute_time_ms"] = (time.time() - start_time) * 1000
        
        return result
        
    except Exception as e:
        import traceback
        return {
            "success": False,
            "error": str(e),
            "traceback": traceback.format_exc(),
            "solutions": [],
            "equations": [],
            "compute_time_ms": (time.time() - start_time) * 1000
        }

def main():
    """Main subprocess entry point - reads JSON from stdin, writes JSON to stdout"""
    try:
        # Read input from stdin
        input_json = sys.stdin.read()
        
        if not input_json.strip():
            error_result = {
                "success": False,
                "error": "No input data received on stdin",
                "solutions": [],
                "equations": [],
                "compute_time_ms": 0
            }
            print(json.dumps(error_result), flush=True)
            return 1
        
        # Parse JSON
        try:
            input_data = json.loads(input_json)
        except json.JSONDecodeError as e:
            error_result = {
                "success": False,
                "error": f"Invalid JSON input: {str(e)}",
                "solutions": [],
                "equations": [],
                "compute_time_ms": 0
            }
            print(json.dumps(error_result), flush=True)
            return 1
        
        # Log to stderr (won't interfere with stdout JSON)
        print(f"[Subprocess] Processing request for: {input_data.get('object_name', 'UNKNOWN')}", 
              file=sys.stderr, flush=True)
        
        # Process calculation (this is where CondensedPhysics gets imported)
        result = process_calculation(input_data)
        
        # Write result to stdout
        output_json = json.dumps(result)
        print(output_json, flush=True)
        
        # Return exit code based on success
        return 0 if result.get("success", False) else 1
        
    except Exception as e:
        import traceback
        error_result = {
            "success": False,
            "error": f"Subprocess error: {str(e)}",
            "traceback": traceback.format_exc(),
            "solutions": [],
            "equations": [],
            "compute_time_ms": 0
        }
        print(json.dumps(error_result), flush=True)
        return 1

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
