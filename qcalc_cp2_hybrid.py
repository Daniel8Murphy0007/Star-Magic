#!/usr/bin/env python3
"""
qcalc_cp2_hybrid.py
====================
Phase 2+3: Hybrid QCalc + CondensedPhysics2 with Production Polish

Intelligent routing:
- Standard UQFF queries → QCalc.UnifiedFieldSolver (fast, 1.1s import)
- Experimental queries → CondensedPhysics2 calculators (advanced physics)

Phase 3 enhancements:
- Error handling with exponential backoff retry
- Circuit breaker pattern for cascading failure prevention
- LRU caching with TTL (1-hour default)
- Progress tracking with status updates

CP2 triggers:
- Query contains: "Orb10", "Orb11", "Orb12", "Orb13", "Orb14", "Orb15", "Orb16"
- Query contains: "Red Mercury", "Silver Mercury", "Plasmoid", "UFEQFET"
- Query contains: "Exp2", "Plasma Flow", "Cosmic Wind"

Performance:
- QCalc path: ~1.1s import + 920ms compute (with cache: <10ms)
- CP2 path: ~2.5s import + variable compute (advanced calculus)

Author: Phase 2+3 (Production)
Date: March 3, 2026
"""

import sys
import json
import time
import re

# Phase 3: Import error handling, caching, and progress modules
try:
    from qcalc_error_handler import (
        retry_with_backoff, classify_error, format_error_response,
        CIRCUIT_BREAKER, ErrorCategory
    )
    from qcalc_cache import GLOBAL_CACHE
    from qcalc_progress import ProgressTracker, CalculationStage, DummyProgressTracker
    PHASE3_ENABLED = True
except ImportError as e:
    print(f"[Warning] Phase 3 modules not available: {str(e)}", file=sys.stderr)
    print("[Warning] Running without error handling, caching, and progress tracking", file=sys.stderr)
    PHASE3_ENABLED = False
    GLOBAL_CACHE = None
    CIRCUIT_BREAKER = None
    
    # Fallback implementations when Phase 3 modules not available
    class DummyProgressTracker:
        """No-op progress tracker for Phase 2 compatibility."""
        def start(self): pass
        def update(self, stage, percent): pass
        def complete(self, result): pass
    
    class CalculationStage:
        """Dummy enum for Phase 2 compatibility."""
        INIT = "INIT"
        IMPORT_MODULES = "IMPORT_MODULES"
        PARSE_INPUT = "PARSE_INPUT"
        VALIDATE_PARAMS = "VALIDATE_PARAMS"
        COMPUTE_FU = "COMPUTE_FU"
        COMPUTE_UG1 = "COMPUTE_UG1"
        COMPUTE_UG2 = "COMPUTE_UG2"
        COMPUTE_UG3 = "COMPUTE_UG3"
        COMPUTE_UG4 = "COMPUTE_UG4"
        COMPUTE_UM = "COMPUTE_UM"
        COMPUTE_UBI = "COMPUTE_UBI"
        FORMAT_RESULTS = "FORMAT_RESULTS"
        COMPLETE = "COMPLETE"

# ═══════════════════════════════════════════════════════════════════════════════
# CP2 QUERY DETECTION
# ═══════════════════════════════════════════════════════════════════════════════

CP2_KEYWORDS = [
    # Orb Analysis series
    r'\bOrb\s*1[0-6]\b', r'\bOrb\s*Analysis\b',
    
    # Experimental physics
    r'\bRed\s*Mercury\b', r'\bSilver\s*Mercury\b',
    r'\bPlasmoid\b', r'\bUFEQFET\b',
    r'\bExp\s*2\b', r'\bExperiment\s*2\b',
    
    # Advanced modules
    r'\bPlasma\s*Flow\b', r'\bCosmic\s*Wind\b',
    r'\bNavier[-\s]?Stokes\b',
    r'\bMagnetic\s*Bubble\b',
    r'\bQuantum\s*Overlay\b',
    r'\bSuperconductor\b', r'\bPropulsion\b',
    
    # Frame sequences
    r'\b(42|44|47|500)\s*Frame\b',
    
    # Energy calculators
    r'\bTotal\s*Energy\s*Budget\b',
    r'\bEnergy\s*Efficiency\b',
    r'\bEnergy\s*Progression\b',
]

def is_cp2_query(query_name: str, params: dict) -> bool:
    """
    Detect if query should route to CondensedPhysics2.
    
    Args:
        query_name: Object name or query identifier
        params: Full parameter dictionary
        
    Returns:
        True if CP2 calculators should handle this query
    """
    # Check query name against CP2 keywords
    query_text = query_name.lower()
    
    for pattern in CP2_KEYWORDS:
        if re.search(pattern, query_text, re.IGNORECASE):
            return True
    
    # Check for CP2-specific parameters
    cp2_params = [
        'orb_number', 'frame_sequence', 'plasmoid_type',
        'mercury_type', 'plasma_viscosity', 'cosmic_wind_velocity'
    ]
    
    if any(key in params for key in cp2_params):
        return True
    
    return False


# ═══════════════════════════════════════════════════════════════════════════════
# CP2 CALCULATOR ROUTER (with Phase 3 progress tracking)
# ═══════════════════════════════════════════════════════════════════════════════

def route_to_cp2(input_data: dict, progress=None) -> dict:
    """
    Route query to appropriate CondensedPhysics2 calculator.
    
    Args:
        input_data: Dictionary with query parameters
        progress: ProgressTracker instance (optional)
        
    Returns:
        Dictionary with computation results
    """
    if progress is None:
        progress = DummyProgressTracker()
    
    calc_start = time.time()
    
    try:
        # Phase 3: Import stage (CP2 is heavy: ~2.5s)
        progress.update(CalculationStage.IMPORT_MODULES, 0)
        import CondensedPhysics2 as CP2
        progress.update(CalculationStage.IMPORT_MODULES, 100)
        
        query_name = input_data.get('object_name', 'UNKNOWN')
        
        # Phase 3: Parse input
        progress.update(CalculationStage.PARSE_INPUT, 50)
        
        # Route to specific calculator based on query
        calculator = None
        result = {}
        
        progress.update(CalculationStage.PARSE_INPUT, 100)
        
        # Phase 3: Computation stage (simplified tracking for CP2)
        progress.update(CalculationStage.COMPUTE_FU, 0)
        
        # Orb Analysis routing
        if re.search(r'\bOrb\s*10\b', query_name, re.IGNORECASE):
            calculator = CP2.MagneticBubbleConfinementOrb10Calculator()
            result = calculator.compute(input_data)
            
        elif re.search(r'\bOrb\s*11\b', query_name, re.IGNORECASE):
            if 'intelligent' in query_name.lower() or 'plasmoid' in query_name.lower():
                calculator = CP2.IntelligentPlasmoidBehaviorOrb11Calculator()
            else:
                calculator = CP2.TotalEnergyBudgetOrb11Calculator()
            result = calculator.compute(input_data)
            
        elif re.search(r'\bOrb\s*12\b', query_name, re.IGNORECASE):
            if '42' in query_name or 'frame' in query_name.lower():
                calculator = CP2.FortyTwoFrameSequenceCalculator()
            elif 'convection' in query_name.lower():
                calculator = CP2.CyclicalConvectionOrb12Calculator()
            else:
                calculator = CP2.Orb12RefinedFUCalculator()
            result = calculator.compute(input_data)
            
        elif re.search(r'\bOrb\s*13\b', query_name, re.IGNORECASE):
            if '44' in query_name:
                calculator = CP2.FortyFourFrameSequenceCalculator()
            elif 'diagonal' in query_name.lower():
                calculator = CP2.DiagonalShiftOrb13Calculator()
            else:
                calculator = CP2.Orb13RefinedFUCalculator()
            result = calculator.compute(input_data)
            
        elif re.search(r'\bOrb\s*14\b', query_name, re.IGNORECASE):
            if '47' in query_name:
                calculator = CP2.FortySevenFrameSequenceCalculator()
            elif 'convection' in query_name.lower():
                calculator = CP2.CyclicalConvectionOrb14Calculator()
            else:
                calculator = CP2.Orb14RefinedFUCalculator()
            result = calculator.compute(input_data)
            
        elif re.search(r'\bOrb\s*15\b', query_name, re.IGNORECASE):
            if '500' in query_name:
                calculator = CP2.FiveHundredFrameDatasetCalculator()
            elif 'spin' in query_name.lower():
                calculator = CP2.PlasmoidSpinRateCalculator()
            else:
                calculator = CP2.FinalizedFURefinementCalculator()
            result = calculator.compute(input_data)
            
        # Experimental physics routing
        elif re.search(r'\bRed\s*Mercury\b', query_name, re.IGNORECASE):
            calculator = CP2.RedMercurySuperconductorCalculator()
            result = calculator.compute(input_data)
            
        elif re.search(r'\bSilver\s*Mercury\b', query_name, re.IGNORECASE):
            calculator = CP2.SilverMercuryPropulsionCalculator()
            result = calculator.compute(input_data)
            
        elif re.search(r'\bNavier[-\s]?Stokes\b', query_name, re.IGNORECASE):
            calculator = CP2.NavierStokesPlasmaFlowCalculator()
            result = calculator.compute(input_data)
            
        elif re.search(r'\bCosmic\s*Wind\b', query_name, re.IGNORECASE):
            if 'disk' in query_name.lower() or 'stability' in query_name.lower():
                calculator = CP2.CosmicWindDiskStabilityCalculator()
            else:
                calculator = CP2.CosmicWindInteractionCalculator()
            result = calculator.compute(input_data)
            
        elif re.search(r'\bPlasmoid\b', query_name, re.IGNORECASE):
            if 'frame' in query_name.lower():
                calculator = CP2.PlasmoidFrameAnalysisCalculator()
            elif 'intelligence' in query_name.lower():
                calculator = CP2.PlasmoidIntelligenceMetricsCalculator()
            elif 'species' in query_name.lower():
                calculator = CP2.PlasmoidSpeciesClassifierCalculator()
            else:
                calculator = CP2.IntelligentQuantumPlasmoidCalculator()
            result = calculator.compute(input_data)
            
        elif re.search(r'\bUFEQFET\b', query_name, re.IGNORECASE):
            calculator = CP2.UFEQFETenComponentCalculator()
            result = calculator.compute(input_data)
            
        else:
            # Fallback: Return info about CP2 capabilities
            result = {
                'success': False,
                'error': f'CP2 query detected but no matching calculator found',
                'query': query_name,
                'available_cp2_modules': [
                    'Orb Analysis 10-16',
                    'Red Mercury Superconductor',
                    'Silver Mercury Propulsion',
                    'Navier-Stokes Plasma Flow',
                    'Cosmic Wind Interactions',
                    'Plasmoid Intelligence Metrics',
                    'UFEQFET Ten-Component'
                ],
                'suggestion': 'Use QCalc.UnifiedFieldSolver for standard UQFF queries'
            }
        
        # Phase 3: Mark computation complete
        progress.update(CalculationStage.COMPUTE_FU, 100)
        
        # Phase 3: Format results
        progress.update(CalculationStage.FORMAT_RESULTS, 50)
        
        # Add timing
        result['compute_time_ms'] = (time.time() - calc_start) * 1000
        result['calculator_used'] = calculator.__class__.__name__ if calculator else 'None'
        result['source'] = 'CondensedPhysics2'
        
        progress.update(CalculationStage.FORMAT_RESULTS, 100)
        
        return result
        
    except ImportError as e:
        return {
            "success": False,
            "error": f"CondensedPhysics2 module not available: {str(e)}",
            "solutions": {},
            "equations": [],
            "compute_time_ms": (time.time() - calc_start) * 1000,
            "source": "CP2 (import failed)"
        }
    except Exception as e:
        return {
            "success": False,
            "error": f"CP2 calculation failed: {str(e)}",
            "error_type": type(e).__name__,
            "query": input_data.get('object_name', 'UNKNOWN'),
            "compute_time_ms": (time.time() - calc_start) * 1000,
            "source": "CP2 (runtime error)"
        }


# ═══════════════════════════════════════════════════════════════════════════════
# QCALC STANDARD PATH (with Phase 3 progress tracking)
# ═══════════════════════════════════════════════════════════════════════════════

def route_to_qcalc(input_data: dict, progress=None) -> dict:
    """
    Route query to QCalc.UnifiedFieldSolver (standard UQFF).
    
    Args:
        input_data: Dictionary with query parameters
        progress: ProgressTracker instance (optional)
        
    Returns:
        Dictionary with computation results
    """
    if progress is None:
        progress = DummyProgressTracker()
    
    calc_start = time.time()
    
    try:
        # Phase 3: Import stage
        progress.update(CalculationStage.IMPORT_MODULES, 0)
        from QCalc import UnifiedFieldSolver, ComputeParams
        progress.update(CalculationStage.IMPORT_MODULES, 100)
        
    except ImportError as e:
        return {
            "success": False,
            "error": f"QCalc module not available: {str(e)}",
            "solutions": {},
            "equations": [],
            "compute_time_ms": 0,
            "source": "QCalc (import failed)"
        }
    
    try:
        # Phase 3: Parse input stage
        progress.update(CalculationStage.PARSE_INPUT, 0)
        
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
        
        progress.update(CalculationStage.PARSE_INPUT, 100)
        
        # Phase 3: Validate parameters stage
        progress.update(CalculationStage.VALIDATE_PARAMS, 50)
        
        # Instantiate solver
        solver = UnifiedFieldSolver()
        
        progress.update(CalculationStage.VALIDATE_PARAMS, 100)
        
        # Phase 3: Computation stages (simplified tracking)
        progress.update(CalculationStage.COMPUTE_FU, 0)
        
        # Call solve() with ComputeParams object
        result = solver.solve(params)
        
        # Mark computation complete
        progress.update(CalculationStage.COMPUTE_FU, 100)
        progress.update(CalculationStage.COMPUTE_UG1, 100)
        progress.update(CalculationStage.COMPUTE_UG2, 100)
        progress.update(CalculationStage.COMPUTE_UG3, 100)
        progress.update(CalculationStage.COMPUTE_UG4, 100)
        progress.update(CalculationStage.COMPUTE_UM, 100)
        progress.update(CalculationStage.COMPUTE_UBI, 100)
        
        # Phase 3: Format results stage
        progress.update(CalculationStage.FORMAT_RESULTS, 50)
        
        # Ensure required fields exist
        if "success" not in result:
            result["success"] = True
        if "compute_time_ms" not in result:
            result["compute_time_ms"] = (time.time() - calc_start) * 1000
        
        result['source'] = 'QCalc'
        
        progress.update(CalculationStage.FORMAT_RESULTS, 100)
        
        return result
        
    except Exception as e:
        return {
            "success": False,
            "error": f"QCalc calculation failed: {str(e)}",
            "error_type": type(e).__name__,
            "query": input_data.get('object_name', 'UNKNOWN'),
            "compute_time_ms": (time.time() - calc_start) * 1000,
            "source": "QCalc (runtime error)"
        }


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN HYBRID ROUTER (with Phase 3 enhancements)
# ═══════════════════════════════════════════════════════════════════════════════

def process_calculation(input_data: dict) -> dict:
    """
    Intelligent router: QCalc vs CP2 with caching, error handling, and progress.
    
    Phase 3 enhancements:
    1. Check cache before calculation
    2. Wrap in circuit breaker
    3. Retry with exponential backoff on transient failures
    4. Track progress
    
    Args:
        input_data: Dictionary with query parameters
        
    Returns:
        Dictionary with computation results
    """
    query_name = input_data.get('object_name', 'UNKNOWN')
    
    # Phase 3: Check cache first
    if PHASE3_ENABLED and GLOBAL_CACHE:
        cached_result = GLOBAL_CACHE.get(input_data)
        if cached_result:
            print(f"[Cache HIT] {query_name} (age: {cached_result.get('__cache_age_seconds__', 0):.1f}s)", 
                  file=sys.stderr)
            return cached_result
        else:
            print(f"[Cache MISS] {query_name}", file=sys.stderr)
    
    # Initialize progress tracker
    if PHASE3_ENABLED:
        progress = ProgressTracker(enable_stderr_output=True)
        progress.start()
    else:
        progress = DummyProgressTracker()
    
    # Route decision
    if is_cp2_query(query_name, input_data):
        print(f"[ROUTER] Detected CP2 query: {query_name}", file=sys.stderr)
        calc_func = lambda: route_to_cp2(input_data, progress)
    else:
        print(f"[ROUTER] Standard UQFF query: {query_name}", file=sys.stderr)
        calc_func = lambda: route_to_qcalc(input_data, progress)
    
    # Phase 3: Execute with error handling and circuit breaker
    if PHASE3_ENABLED and CIRCUIT_BREAKER:
        try:
            # Wrap in circuit breaker
            result = CIRCUIT_BREAKER.call(
                retry_with_backoff,
                calc_func,
                max_retries=3,
                initial_delay=1.0,
                max_delay=10.0
            )
            
            # Store in cache
            if result.get('success') and GLOBAL_CACHE:
                GLOBAL_CACHE.put(input_data, result)
                print(f"[Cache] Stored result for {query_name}", file=sys.stderr)
            
            progress.complete(result)
            return result
            
        except Exception as e:
            # Format error response
            category = classify_error(e)
            error_response = format_error_response(e, category, retries_attempted=3)
            error_response['source'] = 'hybrid_router'
            
            print(f"[ERROR] {type(e).__name__}: {str(e)}", file=sys.stderr)
            
            progress.complete(error_response)
            return error_response
    else:
        # No Phase 3, run directly
        try:
            result = calc_func()
            progress.complete(result)
            return result
        except Exception as e:
            error_response = {
                'success': False,
                'error': str(e),
                'error_type': type(e).__name__,
                'source': 'hybrid_router'
            }
            progress.complete(error_response)
            return error_response


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # Read JSON from stdin
    try:
        input_json = sys.stdin.read()
        input_data = json.loads(input_json)
    except json.JSONDecodeError as e:
        result = {
            "success": False,
            "error": f"Invalid JSON input: {str(e)}",
            "compute_time_ms": 0
        }
        print(json.dumps(result, indent=2))
        sys.exit(1)
    
    # Process calculation
    result = process_calculation(input_data)
    
    # Output JSON to stdout
    print(json.dumps(result, indent=2))
    sys.exit(0 if result.get('success', False) else 1)
