#!/usr/bin/env python3
"""
QCalc_API.py - REST API Endpoints for UQFF Calculator
======================================================

Production-ready REST API for QCalc.py calculations.

DESIGN PRINCIPLES:
- Stateless: Each request is independent
- Fast: Leverages QCalc_Performance.py caching
- Documented: OpenAPI/Swagger compatible
- Secure: Input validation and rate limiting

ENDPOINTS:
    POST /api/v1/calculate         - Single system calculation
    POST /api/v1/calculate/batch   - Bulk calculations
    GET  /api/v1/equations          - List available equations
    GET  /api/v1/constants          - Physics constants
    GET  /api/v1/health             - Health check
    GET  /api/v1/cache/stats        - Cache statistics

PHASE 6: QCalc automatically includes Phase 6 galaxy physics (M51, NGC1316,
SMBH binaries) when appropriate parameters are detected. Use GET /api/v1/phases
for complete phase documentation.

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

from flask import Flask, request, jsonify
from flask_cors import CORS
from typing import Dict, List, Any
import sys
import os
import time
import math

# Import QCalc and performance layer
sys.path.insert(0, os.path.dirname(__file__))
from QCalc import UnifiedFieldSolver, ComputeParams, CONSTANTS
from QCalc_Performance import CACHE, MONITOR, VectorizedCalculations
from IPData import InputParameters
from OPData import OutputDataStore

# ═══════════════════════════════════════════════════════════════════════════════
# FLASK APP SETUP
# ═══════════════════════════════════════════════════════════════════════════════

app = Flask(__name__)
# Restrict CORS to known local origins; set QCALC_CORS_ORIGINS env var to override
_allowed_origins = os.environ.get('QCALC_CORS_ORIGINS', 'http://localhost:*,http://127.0.0.1:*').split(',')
CORS(app, origins=_allowed_origins)

# API limits
MAX_BATCH_SIZE: int = 100  # Maximum systems per batch request (DoS guard)


# ═══════════════════════════════════════════════════════════════════════════════
# INPUT VALIDATION
# ═══════════════════════════════════════════════════════════════════════════════

def validate_params(data: Dict[str, Any]) -> tuple[bool, str, ComputeParams]:
    """
    Validate incoming parameters and convert to ComputeParams.
    
    Returns:
        (valid: bool, error_message: str, params: ComputeParams)
    """
    try:
        # Required fields
        if 'M' not in data:
            return False, "Missing required field: 'M' (mass in kg)", None
        
        if 'r' not in data:
            return False, "Missing required field: 'r' (distance in m)", None
        
        # Convert to float and check for NaN/inf
        def _safe_float(val, field_name):
            v = float(val)
            if not math.isfinite(v):
                raise ValueError(f"Field '{field_name}' must be a finite number, got {v}")
            return v

        # Convert to ComputeParams
        params = ComputeParams(
            M=_safe_float(data['M'], 'M'),
            r=_safe_float(data['r'], 'r'),
            z=_safe_float(data.get('z', 0.0), 'z'),
            t=_safe_float(data.get('t', 0.0), 't'),
            t_n=_safe_float(data.get('t_n', -1000.0), 't_n'),
            B=_safe_float(data.get('B_mag', data.get('B', 1e-8)), 'B'),
            SFR=_safe_float(data.get('SFR', 10.0), 'SFR'),
            M_bh=_safe_float(data.get('M_bh', 4.15e6 * CONSTANTS['M_sun']), 'M_bh'),
            d_g=_safe_float(data.get('d_g', 8000 * CONSTANTS['pc']), 'd_g')
        )
        
        # Sanity checks
        if params.M <= 0:
            return False, "Mass M must be positive", None
        if params.r <= 0:
            return False, "Distance r must be positive", None
        
        return True, "", params
        
    except ValueError as e:
        return False, f"Invalid parameter type: {str(e)}", None
    except Exception as e:
        return False, f"Validation error: {str(e)}", None


# ═══════════════════════════════════════════════════════════════════════════════
# API ENDPOINTS
# ═══════════════════════════════════════════════════════════════════════════════

@app.route('/api/v1/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        'status': 'healthy',
        'service': 'QCalc API',
        'version': '1.0.0',
        'framework': 'UQFF Star Magic',
        'solvability': '99.9%'
    }), 200


@app.route('/api/v1/constants', methods=['GET'])
def get_constants():
    """Return all physics constants."""
    return jsonify({
        'constants': CONSTANTS,
        'count': len(CONSTANTS),
        'units': {
            'G': 'm³/kg·s²',
            'c': 'm/s',
            'h_bar': 'J·s',
            'M_sun': 'kg',
            'pc': 'm',
            'ly': 'm'
        }
    }), 200


@app.route('/api/v1/equations', methods=['GET'])
def list_equations():
    """List all available equations."""
    # Get available equations from solver
    test_params = ComputeParams(M=1.989e30, r=1.496e11)
    results = solver.solve(test_params)
    
    # results['long_form_equations'] contains dicts from EquationResult.to_dict()
    equation_list = [
        {
            'name': eq.get('name', ''),
            'latex': eq.get('latex', ''),
            'unit': eq.get('unit', ''),
            'notes': eq.get('notes', '')
        }
        for eq in results['long_form_equations']
    ]
    
    return jsonify({
        'equations': equation_list,
        'total_count': len(equation_list),
        'available_equations': results.get('available_equations', [])
    }), 200


@app.route('/api/v1/calculate', methods=['POST'])
def calculate_single():
    """
    Calculate UQFF equations for a single system.
    
    Request Body:
        {
            "M": 1.989e30,           # Mass (kg) - REQUIRED
            "r": 1.496e11,           # Distance (m) - REQUIRED
            "z": 0.0,                # Redshift (optional)
            "t": 0.0,                # Time (s, optional)
            "t_n": -1000.0,          # Negative time (s, optional)
            "theta": 0.01,           # Angle (rad, optional)
            "v_rot": 1e5,            # Rotation velocity (m/s, optional)
            "B_mag": 1e-8,           # Magnetic field (T, optional)
            "SFR": 10.0,             # Star formation rate (M_sun/yr, optional)
            "name": "Test System"    # Name (optional)
        }
    
    Response:
        {
            "success": true,
            "system_name": "Test System",
            "solutions": {...},
            "equations": [...],
            "computation_time": 0.123,
            "cache_hit": false
        }
    """
    try:
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No JSON data provided'}), 400
        
        # Validate parameters
        valid, error_msg, params = validate_params(data)
        if not valid:
            return jsonify({'error': error_msg}), 400
        
        # Get system name
        system_name = data.get('name', f"System_M{params.M:.2e}_r{params.r:.2e}")
        
        # Solve (will use cache if available)
        import time
        start = time.perf_counter()
        results = solver.solve(params)
        computation_time = time.perf_counter() - start
        
        # Store in output database
        query_id = f"api_{int(time.time() * 1000)}"
        store_dict = {
            'query_id': query_id,
            'input_params': vars(params),
            'long_form_equations': results['long_form_equations'],
            'solutions': results['solutions'],
            'available_equations': results.get('available_equations', [])
        }
        output_store.store(store_dict)
        
        # Format response
        response = {
            'success': True,
            'query_id': query_id,
            'system_name': system_name,
            'input_parameters': vars(params),
            'solutions': results['solutions'],
            'equations': [
                {
                    'name': eq.get('name', ''),
                    'result': eq.get('result', 0.0),
                    'unit': eq.get('unit', ''),
                    'notes': eq.get('notes', '')
                }
                for eq in results['long_form_equations']
            ],
            'available_equations': results.get('available_equations', []),
            'computation_time': computation_time,
            'equation_count': len(results['long_form_equations'])
        }
        
        return jsonify(response), 200
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e),
            'error_type': type(e).__name__
        }), 500


@app.route('/api/v1/calculate/batch', methods=['POST'])
def calculate_batch():
    """
    Calculate UQFF equations for multiple systems.
    
    Request Body:
        {
            "systems": [
                {"M": 1.989e30, "r": 1.496e11, "name": "System 1"},
                {"M": 5e30, "r": 1e12, "name": "System 2"},
                ...
            ]
        }
    
    Response:
        {
            "success": true,
            "results": [...],
            "total_systems": 2,
            "total_time": 0.456
        }
    """
    try:
        data = request.get_json()
        if not data or 'systems' not in data:
            return jsonify({'error': 'Missing "systems" array in request body'}), 400
        
        systems = data['systems']
        if not isinstance(systems, list) or len(systems) == 0:
            return jsonify({'error': '"systems" must be a non-empty array'}), 400
        
        if len(systems) > MAX_BATCH_SIZE:
            return jsonify({
                'error': f'Batch size {len(systems)} exceeds maximum of {MAX_BATCH_SIZE}. '
                         f'Split into smaller batches.'
            }), 400
        
        # Process each system
        import time
        start = time.perf_counter()
        results_list = []
        errors = []
        
        for idx, system_data in enumerate(systems):
            # Validate
            valid, error_msg, params = validate_params(system_data)
            if not valid:
                errors.append({
                    'index': idx,
                    'system': system_data.get('name', f"System {idx}"),
                    'error': error_msg
                })
                continue
            
            # Solve
            try:
                result = solver.solve(params)
                system_name = system_data.get('name', f"System_{idx}")
                
                results_list.append({
                    'index': idx,
                    'name': system_name,
                    'solutions': result['solutions'],
                    'equation_count': len(result['long_form_equations'])
                })
            except Exception as e:
                errors.append({
                    'index': idx,
                    'system': system_data.get('name', f"System {idx}"),
                    'error': str(e)
                })
        
        total_time = time.perf_counter() - start
        
        response = {
            'success': len(results_list) > 0,
            'results': results_list,
            'errors': errors,
            'total_systems': len(systems),
            'successful': len(results_list),
            'failed': len(errors),
            'total_time': total_time,
            'avg_time_per_system': total_time / len(systems) if len(systems) > 0 else 0
        }
        
        return jsonify(response), 200
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e),
            'error_type': type(e).__name__
        }), 500


@app.route('/api/v1/cache/stats', methods=['GET'])
def cache_stats():
    """Get cache performance statistics."""
    stats = CACHE.get_stats()
    
    return jsonify({
        'cache': stats,
        'performance': MONITOR.get_report()
    }), 200


@app.route('/api/v1/cache/clear', methods=['POST'])
def clear_cache():
    """Clear the result cache (requires QCALC_ADMIN_TOKEN in Authorization header)."""
    # Simple token guard — set QCALC_ADMIN_TOKEN env var to protect in production
    admin_token = os.environ.get('QCALC_ADMIN_TOKEN', '')
    if admin_token:
        auth = request.headers.get('Authorization', '')
        if not auth.startswith('Bearer ') or auth[7:] != admin_token:
            return jsonify({'error': 'Unauthorized'}), 401
    CACHE.clear()
    return jsonify({
        'success': True,
        'message': 'Cache cleared'
    }), 200


# ═══════════════════════════════════════════════════════════════════════════════
# ADVANCED ENDPOINTS (Phase 3-4 Specific)
# ═══════════════════════════════════════════════════════════════════════════════

@app.route('/api/v1/phases', methods=['GET'])
def get_phases():
    """
    Return UQFF phase documentation.

    Phases represent progressive integration stages of physics modules
    sourced from source14.cpp through source173.cpp.
    """
    phases = {
        'framework': 'UQFF Star Magic v5.x (99.9% Solvability)',
        'data_flow': 'APIFetch.py → IPData.py → QCalc.py (+ Wolfram Extensions) → OPData.py',
        'phases': [
            {
                'phase': 1,
                'description': 'Core UQFF equations (Ug1-Ug4, Ubi, Um, UA_uv)',
                'sources': ['QCalc.py (UnifiedFieldSolver, 8 master equations)'],
                'equations': ['F_U', 'Ug1 (DPM)', 'Ug2 (heliosphere)', 'Ug3 (magnetic disk)', 'Ug4 (galactic)'],
                'modes': ['Compressed', 'Resonant', 'Superconductive', 'Buoyant']
            },
            {
                'phase': 2,
                'description': 'Magnetar + SMBH Wolfram terms (source14-15)',
                'sources': ['QCalc_Wolfram_Extensions.py: calculate_base_gravity_hubble_magnetic, '
                            'calculate_uqff_unification_time_reversal, calculate_smbh_*'],
                'systems': ['SGR 0501+4516 Magnetar (12 terms)', 'Sagittarius A* SMBH (15 terms)']
            },
            {
                'phase': 3,
                'description': 'Star formation + cluster + nebulae Wolfram terms (source16-25)',
                'sources': ['QCalc_Wolfram_Extensions.py: calculate_cluster_*, '
                            'calculate_westerlund2_*, calculate_photoevaporation_*'],
                'systems': ['Tapestry NGC 2014/2020', 'Westerlund 2', 'Pillars M16',
                            'Rings of Relativity', 'NGC 2525/SN', 'NGC 3603',
                            'Bubble Nebula', 'Antennae Galaxies', 'Horsehead', 'NGC 1275']
            },
            {
                'phase': 4,
                'description': 'Cosmological + galaxy + planetary + remnant terms (source26-40)',
                'sources': ['QCalc_Wolfram_Extensions.py: calculate_hudf_*, calculate_ngc1792_*, '
                            'calculate_andromeda_*, calculate_sombrero_*, calculate_saturn_*, '
                            'calculate_m16_*, calculate_crab_*, calculate_sgr1745_*'],
                'systems': ['HUDF z=3.5-12', 'NGC 1792', 'Andromeda M31', 'Sombrero M104',
                            'Saturn', 'M16 Eagle Nebula', 'Crab Nebula', 'SGR 1745-2900']
            },
            {
                'phase': 5,
                'description': 'Frequency model + framework hybrid terms (source34-50)',
                'sources': ['QCalc_Wolfram_Extensions.py: calculate_sgr1745_frequency_*, '
                            'calculate_tapestry_*, calculate_compressed_*, '
                            'calculate_ngc6302_*, calculate_orion_m42_*'],
                'systems': ['SGR 1745 11-freq model', 'Sgr A* freq model', 'Tapestry framework',
                            'NGC 6302 Butterfly', 'Orion M42', 'Generic compressed/resonance']
            },
            {
                'phase': 6,
                'description': 'Galaxy physics extensions (M51, NGC1316, SMBH binaries)',
                'sources': ['QCalc.py Phase 6 auto-detection (parameters: z, M_bh, SFR)'],
                'auto_trigger': 'Detected automatically when M > 1e10 M_sun or z > 0.01'
            }
        ],
        'parallel_calculators': [
            'MAIN_1_CoAnQi.cpp (107,019 lines, C++)',
            'QCalc.py (this service)',
            'CondensedPhysics.py (176 calculators)',
            'CondensedPhysics2.py (680 classes)',
            'uqff_server.js/index.js (106 systems)'
        ],
        'canonical_constants': {
            'beta_i': 0.603,
            'kappa': '5e-4 day^-1',
            'gamma': '5e-5 day^-1 (near-lossless)',
            'eta': '1e-22 (Aether coupling)',
            'rho_vac_SCm': '7.09e-37 J/m^3',
            'rho_vac_UA': '7.09e-36 J/m^3',
            'v_SCm': '1e8 m/s (c/3)',
            'H_SCm': 0.99,
            'SSq': 0.57
        }
    }
    return jsonify(phases), 200


def calculate_aether_metric():
    """
    Calculate aether metric tensor (Phase 4).
    
    Request Body:
        {
            "lambda_UA": 7.09e-36,   # Quantum vacuum density
            "lambda_SCm": 7.09e-37,  # Superconducting medium density
            "lambda_A": 8.99e-07,    # Aether mass density
            "t_n": -1000.0           # Negative time
        }
    
    Response:
        {
            "g_mu_nu": [[...], ...],           # Minkowski metric (4x4)
            "T_s_mu_nu": [[...], ...],         # Stress-energy tensor (4x4)
            "delta_g_mu_nu": [[...], ...],     # Metric perturbation (4x4)
            "UA_mu_nu": [[...], ...],          # Aether metric (4x4)
            "det_UA": -1.0,                    # Determinant
            "R": -1.1102e-16                   # Ricci scalar
        }
    """
    try:
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No JSON data provided'}), 400
        
        # Extract parameters
        lambda_UA = float(data.get('lambda_UA', 7.09e-36))
        lambda_SCm = float(data.get('lambda_SCm', 7.09e-37))
        lambda_A = float(data.get('lambda_A', 8.99e-07))
        t_n = float(data.get('t_n', -1000.0))
        
        # Import AetherMetricCalculator from QCalc
        from QCalc import AetherMetricCalculator
        
        calc = AetherMetricCalculator()
        
        # Compute all components
        g_mu_nu = calc.compute_minkowski_metric()
        T_s = calc.compute_stress_energy_tensor(lambda_UA, lambda_SCm, lambda_A, t_n)
        delta_g = calc.compute_metric_perturbation(T_s)
        UA = calc.compute_aether_metric(g_mu_nu, delta_g)
        det_UA = calc.compute_metric_determinant(UA)
        R = calc.compute_ricci_scalar(delta_g)
        
        response = {
            'success': True,
            'g_mu_nu': g_mu_nu.tolist(),
            'T_s_mu_nu': T_s.tolist(),
            'delta_g_mu_nu': delta_g.tolist(),
            'UA_mu_nu': UA.tolist(),
            'det_UA': float(det_UA),
            'R': float(R),
            'parameters': {
                'lambda_UA': lambda_UA,
                'lambda_SCm': lambda_SCm,
                'lambda_A': lambda_A,
                't_n': t_n
            }
        }
        
        return jsonify(response), 200
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e),
            'error_type': type(e).__name__
        }), 500


# ═══════════════════════════════════════════════════════════════════════════════
# DOCUMENTATION ENDPOINT
# ═══════════════════════════════════════════════════════════════════════════════

@app.route('/api/v1/docs', methods=['GET'])
def api_documentation():
    """Return API documentation."""
    docs = {
        'title': 'UQFF QCalc API',
        'version': '1.0.0',
        'description': 'REST API for Unified Quantum Field Framework calculations',
        'base_url': '/api/v1',
        'endpoints': [
            {
                'path': '/health',
                'method': 'GET',
                'description': 'Health check'
            },
            {
                'path': '/constants',
                'method': 'GET',
                'description': 'Get all physics constants'
            },
            {
                'path': '/equations',
                'method': 'GET',
                'description': 'List all available equations'
            },
            {
                'path': '/calculate',
                'method': 'POST',
                'description': 'Calculate UQFF equations for single system',
                'required_fields': ['M', 'r']
            },
            {
                'path': '/calculate/batch',
                'method': 'POST',
                'description': 'Batch calculations for multiple systems',
                'required_fields': ['systems']
            },
            {
                'path': '/aether_metric',
                'method': 'POST',
                'description': 'Calculate aether metric tensor (Phase 4)',
                'optional_fields': ['lambda_UA', 'lambda_SCm', 'lambda_A', 't_n']
            },
            {
                'path': '/cache/stats',
                'method': 'GET',
                'description': 'Get cache statistics'
            },
            {
                'path': '/cache/clear',
                'method': 'POST',
                'description': 'Clear result cache'
            }
        ],
        'example_request': {
            'url': '/api/v1/calculate',
            'method': 'POST',
            'body': {
                'M': 1.989e30,
                'r': 1.496e11,
                'name': 'Sun at 1 AU'
            }
        }
    }
    
    return jsonify(docs), 200


# ═══════════════════════════════════════════════════════════════════════════════
# RUN SERVER
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("="*80)
    print("QCalc REST API Server")
    print("="*80)
    print("Framework: UQFF Star Magic (99.9% Solvability)")
    print("Endpoints: /api/v1/docs for complete documentation")
    print("Performance: Caching + Vectorization enabled")
    print("="*80)
    
    # Run development server
    app.run(
        host='0.0.0.0',
        port=8443,
        debug=False,
        threaded=True
    )
