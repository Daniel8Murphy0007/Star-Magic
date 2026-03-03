#!/usr/bin/env python3
"""
QCalc_Start_API.py - Flask REST API Server Launcher
====================================================

Quick start script for QCalc_API.py that:
1. Verifies dependencies (Flask, QCalc, CP2)
2. Sets appropriate PYTHONPATH
3. Starts server on port 5000 (or custom port)

Phase 2: Now includes CondensedPhysics2 routing for experimental queries

Usage:
    python QCalc_Start_API.py            # Port 5000
    python QCalc_Start_API.py --port 8443  # Custom port

Author: Phase 2 Extension
Date: March 3, 2026
"""

import sys
import os
import argparse

def check_dependencies():
    """Verify all required modules are available"""
    missing = []
    
    try:
        import flask
        print("[✓] Flask installed:", flask.__version__)
    except ImportError:
        missing.append("flask")
        print("[✗] Flask not found - install with: pip install flask flask-cors")
    
    try:
        from QCalc import UnifiedFieldSolver, ComputeParams
        print("[✓] QCalc module available (9,149 lines)")
    except ImportError:
        missing.append("QCalc")
        print("[✗] QCalc.py not found in current directory")
    
    try:
        import CondensedPhysics2 as CP2
        print(f"[✓] CondensedPhysics2 module available ({CP2.CP2_CLASS_COUNT} classes)")
    except ImportError:
        print("[!] CondensedPhysics2.py not found (optional for experimental queries)")
    
    if missing:
        print("\n[ERROR] Missing required dependencies:")
        for mod in missing:
            print(f"  - {mod}")
        return False
    
    return True


def test_qcalc():
    """Quick sanity test of QCalc"""
    try:
        from QCalc import UnifiedFieldSolver, ComputeParams
        
        solver = UnifiedFieldSolver()
        params = ComputeParams(M=1.989e30, r=1.496e11)  # Sun at 1 AU
        
        result = solver.solve(params)
        
        if result.get('success'):
            print("[✓] QCalc test successful")
            print(f"    F_U = {result['solutions'].get('F_U', 'N/A')} N/kg")
            return True
        else:
            print("[✗] QCalc test failed:", result.get('error', 'Unknown error'))
            return False
            
    except Exception as e:
        print(f"[✗] QCalc test error: {str(e)}")
        return False


def start_server(port=5000, debug=True):
    """Start Flask API server"""
    print("="*80)
    print("QCalc REST API Server (Phase 2: QCalc + CP2 Hybrid)")
    print("="*80)
    print(f"Port: {port}")
    print(f"Debug mode: {debug}")
    print(f"API docs: http://localhost:{port}/api/v1/docs")
    print("="*80)
    
    # Import and run Flask app
    from QCalc_API import app
    
    app.run(
        host='0.0.0.0',
        port=port,
        debug=debug,
        threaded=True
    )


def main():
    parser = argparse.ArgumentParser(description='Start QCalc REST API Server')
    parser.add_argument('--port', type=int, default=5000, help='Port (default: 5000)')
    parser.add_argument('--no-debug', action='store_true', help='Disable debug mode')
    parser.add_argument('--skip-tests', action='store_true', help='Skip dependency checks')
    
    args = parser.parse_args()
    
    print("\n[Phase 2] QCalc REST API Launcher\n")
    
    if not args.skip_tests:
        print("1. Checking dependencies...")
        if not check_dependencies():
            print("\n[ERROR] Dependency check failed. Fix errors and try again.")
            sys.exit(1)
        
        print("\n2. Testing QCalc...")
        if not test_qcalc():
            print("\n[ERROR] QCalc test failed. Check installation.")
            sys.exit(1)
        
        print("\n3. Starting server...")
    
    try:
        start_server(port=args.port, debug=not args.no_debug)
    except KeyboardInterrupt:
        print("\n\n[Server] Shutting down gracefully...")
        sys.exit(0)
    except Exception as e:
        print(f"\n[ERROR] Server crashed: {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    main()
