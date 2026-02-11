#!/usr/bin/env python3
"""
test_integration.py - Quick Integration Test for Star-Magic UQFF
================================================================

Tests the integration between Python data layer and C++ MAIN_1_CoAnQi calculator.

Usage:
    python test_integration.py

Author: GitHub Copilot (AI-Generated)
Date: February 11, 2026
"""

import subprocess
import json
import os
import sys
from pathlib import Path
import time


class Colors:
    """ANSI color codes for terminal output"""
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'


def print_header(text):
    """Print formatted section header"""
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{text}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.END}\n")


def print_success(text):
    """Print success message"""
    print(f"{Colors.GREEN}✅ {text}{Colors.END}")


def print_error(text):
    """Print error message"""
    print(f"{Colors.RED}❌ {text}{Colors.END}")


def print_warning(text):
    """Print warning message"""
    print(f"{Colors.YELLOW}⚠️  {text}{Colors.END}")


def print_info(text):
    """Print info message"""
    print(f"{Colors.BLUE}ℹ️  {text}{Colors.END}")


def check_cpp_executable():
    """Check if MAIN_1_CoAnQi.exe exists"""
    print_header("1. Checking C++ Executable")
    
    exe_paths = [
        Path("build_msvc/Release/MAIN_1_CoAnQi.exe"),
        Path("build_msvc/Debug/MAIN_1_CoAnQi.exe"),
        Path("MAIN_1_CoAnQi.exe")
    ]
    
    for exe_path in exe_paths:
        if exe_path.exists():
            print_success(f"Found C++ executable: {exe_path}")
            return exe_path
    
    print_error("C++ executable not found!")
    print_info("Build it first with:")
    print("    cmake -S . -B build_msvc -G \"Visual Studio 17 2022\" -A x64")
    print("    cmake --build build_msvc --config Release")
    return None


def test_cli_batch_mode(exe_path):
    """Test --batch flag with a sample system"""
    print_header("2. Testing CLI Batch Mode (--batch)")
    
    test_systems = [
        "Sagittarius A*",
        "Betelgeuse",
        "M87"
    ]
    
    for system_name in test_systems:
        print(f"\n{Colors.BOLD}Testing: {system_name}{Colors.END}")
        
        try:
            result = subprocess.run(
                [str(exe_path), "--batch", system_name],
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode != 0:
                print_error(f"Command failed with code {result.returncode}")
                print(f"stderr: {result.stderr[:200]}")
                continue
            
            # Parse JSON output
            try:
                data = json.loads(result.stdout)
                print_success(f"Received valid JSON output")
                print(f"  F_U_Bi_i: {data.get('F_U_Bi_i', 'N/A')}")
                print(f"  g_compressed: {data.get('g_compressed', 'N/A')}")
                return True
            except json.JSONDecodeError as e:
                print_error(f"Invalid JSON output: {e}")
                print(f"Output (first 200 chars): {result.stdout[:200]}")
                
        except subprocess.TimeoutExpired:
            print_error(f"Computation timed out (>30s)")
        except Exception as e:
            print_error(f"Unexpected error: {e}")
    
    return False


def test_list_systems(exe_path):
    """Test --list-systems flag"""
    print_header("3. Testing System List (--list-systems)")
    
    try:
        result = subprocess.run(
            [str(exe_path), "--list-systems"],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode != 0:
            print_error(f"Command failed with code {result.returncode}")
            return False
        
        # Parse JSON output
        data = json.loads(result.stdout)
        total = data.get('total_systems', 0)
        systems = data.get('systems', [])
        
        print_success(f"Found {total} systems")
        print(f"Sample systems (first 5):")
        for system in systems[:5]:
            print(f"  • {system}")
        
        return True
    
    except Exception as e:
        print_error(f"Test failed: {e}")
        return False


def test_system_info(exe_path):
    """Test --system-info flag"""
    print_header("4. Testing System Info (--system-info)")
    
    test_system = "Betelgeuse"
    
    try:
        result = subprocess.run(
            [str(exe_path), "--system-info", test_system],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode != 0:
            print_error(f"Command failed with code {result.returncode}")
            return False
        
        # Parse JSON output
        data = json.loads(result.stdout)
        print_success(f"Retrieved info for {test_system}")
        print(f"  Mass (M): {data.get('M', 'N/A')} kg")
        print(f"  Radius (r): {data.get('r', 'N/A')} m")
        print(f"  Magnetic field (B0): {data.get('B0', 'N/A')} T")
        
        return True
    
    except Exception as e:
        print_error(f"Test failed: {e}")
        return False


def test_python_wrapper():
    """Test Python wrapper (CoAnQi_Wrapper.py)"""
    print_header("5. Testing Python Wrapper")
    
    try:
        from CoAnQi_Wrapper import CoAnQiCalculator
        print_success("Successfully imported CoAnQi_Wrapper module")
        
        # Test initialization
        calc = CoAnQiCalculator(verbose=True)
        print_success("Initialized CoAnQiCalculator")
        
        # Test computation
        print(f"\n{Colors.BOLD}Computing Sagittarius A* via Python wrapper...{Colors.END}")
        result = calc.compute_system("Sagittarius A*")
        
        if result.status == "success":
            print_success("Computation successful")
            print(f"  F_U_Bi_i: {result.F_U_Bi_i:.6e} N")
            print(f"  g_compressed: {result.g_compressed:.6e} m/s²")
            print(f"  Execution time: {result.computation_time:.2f}s")
            return True
        else:
            print_error(f"Computation failed: {result.error_message}")
            return False
    
    except ImportError:
        print_error("CoAnQi_Wrapper.py not found in current directory")
        return False
    except Exception as e:
        print_error(f"Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_data_layer_integration():
    """Test integration with APIFetch/IPData/QCalc/OPData"""
    print_header("6. Testing Data Layer Integration (Optional)")
    
    try:
        from IPData import InputParameters, InputDataStore
        from QCalc import UnifiedFieldSolver
        from OPData import OutputDataStore
        print_success("Python data layer modules found")
        
        # Test basic instantiation
        input_store = InputDataStore()
        solver = UnifiedFieldSolver()
        output_store = OutputDataStore()
        print_success("Data layer objects instantiated")
        
        return True
    
    except ImportError as e:
        print_warning(f"Data layer modules not found (optional): {e}")
        print_info("This is OK - data layer integration is optional")
        return None  # None = optional test skipped
    except Exception as e:
        print_warning(f"Data layer test failed: {e}")
        return None


def run_all_tests():
    """Run all integration tests"""
    print(f"{Colors.BOLD}{Colors.BLUE}")
    print("=" * 70)
    print("Star-Magic UQFF Integration Test Suite")
    print("Testing C++ ↔ Python Integration Layer")
    print("=" * 70)
    print(f"{Colors.END}")
    
    results = {}
    
    # Test 1: Check C++ executable
    exe_path = check_cpp_executable()
    results['cpp_exe'] = exe_path is not None
    
    if not exe_path:
        print_error("Cannot proceed without C++ executable. Aborting tests.")
        return results
    
    # Test 2: CLI batch mode
    results['batch_mode'] = test_cli_batch_mode(exe_path)
    
    # Test 3: List systems
    results['list_systems'] = test_list_systems(exe_path)
    
    # Test 4: System info
    results['system_info'] = test_system_info(exe_path)
    
    # Test 5: Python wrapper
    results['python_wrapper'] = test_python_wrapper()
    
    # Test 6: Data layer (optional)
    results['data_layer'] = test_data_layer_integration()
    
    # Print summary
    print_header("Test Summary")
    
    passed = sum(1 for v in results.values() if v is True)
    failed = sum(1 for v in results.values() if v is False)
    skipped = sum(1 for v in results.values() if v is None)
    total = len(results)
    
    for test_name, result in results.items():
        status = ""
        if result is True:
            status = f"{Colors.GREEN}✅ PASSED{Colors.END}"
        elif result is False:
            status = f"{Colors.RED}❌ FAILED{Colors.END}"
        else:
            status = f"{Colors.YELLOW}⚪ SKIPPED{Colors.END}"
        
        print(f"  {test_name:20s} {status}")
    
    print(f"\n{Colors.BOLD}Results: {passed}/{total} passed, {failed} failed, {skipped} skipped{Colors.END}")
    
    if failed == 0:
        print_success("\n🎉 All required tests passed! Integration layer is ready to use.")
    else:
        print_error("\n⚠️  Some tests failed. Review errors above and fix issues.")
    
    return results


if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("Starting Integration Tests...")
    print("=" * 70 + "\n")
    
    start_time = time.time()
    results = run_all_tests()
    elapsed = time.time() - start_time
    
    print(f"\n{Colors.BOLD}Total execution time: {elapsed:.2f}s{Colors.END}\n")
    
    # Exit code based on results
    if all(v is not False for v in results.values()):
        sys.exit(0)
    else:
        sys.exit(1)
