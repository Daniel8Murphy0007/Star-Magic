#!/usr/bin/env python3
"""
condensed_physics_subprocess.py - Subprocess Wrapper for CondensedPhysics.py
=============================================================================

Called by source2(HEAD PROGRAM).cpp via QProcess to perform physics calculations.

Input:  JSON via stdin with InputData parameters
Output: JSON via stdout with computation results

Usage:
    python condensed_physics_subprocess.py

Data Flow:
    source2(HEAD).cpp → QProcess → this script → CondensedPhysics.solve() → stdout → source2(HEAD).cpp

Author: Daniel T. Murphy
Date: March 2, 2026
Phase: 0 - Unification (IPC Wiring)
"""

import sys
import json
import traceback
from typing import Dict, Any

# Import CondensedPhysics calculator
try:
    from CondensedPhysics import UnifiedFieldSolver
    CONDENSED_PHYSICS_AVAILABLE = True
except ImportError:
    CONDENSED_PHYSICS_AVAILABLE = False
    print(json.dumps({"error": "CondensedPhysics.py not available"}), file=sys.stderr, flush=True)

# Import CondensedPhysics2 calculators via Aggregator (CP2 wiring, Problem 4)
# CP2_AVAILABLE flags availability of all 529 CP2 calculators for routing.
# Callers of process_calculation() may check CP2_AVAILABLE to enable CP2 paths.
try:
    from CondensedPhysicsAggregator import CondensedPhysicsAggregator
    CP2_AVAILABLE = True
except ImportError:
    CP2_AVAILABLE = False

def process_calculation(input_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process physics calculation request using CondensedPhysics.
    
    Args:
        input_data: Dictionary with keys:
            - object_name: str (e.g., "SGR 1745+29")
            - M: float (mass in solar masses, optional)
            - r: float (radius in meters, optional)
            - z: float (redshift, optional)
            - B: float (magnetic field in Tesla, optional)
            - T: float (temperature in K, optional)
            - SFR: float (star formation rate, optional)
            ... (additional InputData fields)
    
    Returns:
        Dictionary with:
            - success: bool
            - long_form_equations: List[str] (if success)
            - solutions: Dict (if success)
            - available_equations: List[str] (if success)
            - query_id: str (if success)
            - error: str (if failure)
            - compute_time_ms: float
    """
    import time
    start_time = time.time()
    
    try:
        if not CONDENSED_PHYSICS_AVAILABLE:
            return {
                "success": False,
                "error": "CondensedPhysics.py module not available",
                "compute_time_ms": (time.time() - start_time) * 1000
            }
        
        # Create solver instance
        solver = UnifiedFieldSolver()
        
        # Call solve() with input parameters
        result = solver.solve(input_data)
        
        # Add success flag and timing
        result["success"] = True
        result["compute_time_ms"] = (time.time() - start_time) * 1000
        
        return result
        
    except Exception as e:
        # Capture any errors
        error_traceback = traceback.format_exc()
        return {
            "success": False,
            "error": str(e),
            "traceback": error_traceback,
            "compute_time_ms": (time.time() - start_time) * 1000
        }

def main():
    """
    Main subprocess loop: reads JSON from stdin, processes, writes JSON to stdout.
    """
    try:
        # Read input JSON from stdin
        input_line = sys.stdin.readline()
        
        if not input_line:
            print(json.dumps({
                "success": False,
                "error": "No input received on stdin"
            }), flush=True)
            sys.exit(1)
        
        # Parse JSON input
        try:
            input_data = json.loads(input_line)
        except json.JSONDecodeError as e:
            print(json.dumps({
                "success": False,
                "error": f"Invalid JSON input: {str(e)}"
            }), flush=True)
            sys.exit(1)
        
        # Process calculation
        result = process_calculation(input_data)
        
        # Output result as JSON to stdout
        print(json.dumps(result), flush=True)
        
        # Exit successfully
        sys.exit(0 if result["success"] else 1)
        
    except Exception as e:
        # Catch-all error handler
        error_result = {
            "success": False,
            "error": f"Subprocess error: {str(e)}",
            "traceback": traceback.format_exc()
        }
        print(json.dumps(error_result), flush=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
