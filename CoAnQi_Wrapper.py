"""
CoAnQi_Wrapper.py - Python Interface to MAIN_1_CoAnQi.cpp C++ Calculator
========================================================================

Provides subprocess-based integration between Python data layer and C++ computational core.

Usage:
    from CoAnQi_Wrapper import CoAnQiCalculator
    
    calc = CoAnQiCalculator()
    result = calc.compute_system("Sagittarius A*")
    print(f"F_U_Bi_i = {result.F_U_Bi_i} N")

Author: GitHub Copilot (AI-Generated)
Date: February 11, 2026
"""

import subprocess
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Union
from dataclasses import dataclass, asdict
import time


@dataclass
class CoAnQiResult:
    """Results from C++ MAIN_1_CoAnQi calculator"""
    system_name: str
    F_U_Bi_i: float                 # Universal buoyancy force (N)
    g_compressed: float             # Compressed gravity (m/s²)
    F_jet_rel: Optional[float] = None      # Relativistic jet force (N)
    E_acc_rel: Optional[float] = None      # Relativistic accretion energy (J)
    F_drag_rel: Optional[float] = None     # Relativistic drag force (N)
    F_gw_rel: Optional[float] = None       # Gravitational wave force (N)
    computation_time: Optional[float] = None  # Execution time (seconds)
    status: str = "success"         # "success", "error", "timeout"
    error_message: Optional[str] = None
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization"""
        return asdict(self)
    
    def to_json(self) -> str:
        """Convert to JSON string"""
        return json.dumps(self.to_dict(), indent=2)


class CoAnQiCalculator:
    """
    Python interface to MAIN_1_CoAnQi.cpp C++ calculator
    
    Integration Methods:
    1. CLI batch mode (--batch flag) - RECOMMENDED
    2. Interactive menu automation (stdin piping)
    3. REST API (requires server mode in C++)
    
    Attributes:
        exe_path (Path): Path to MAIN_1_CoAnQi.exe
        timeout (float): Max execution time in seconds
        verbose (bool): Enable debug output
    """
    
    def __init__(
        self, 
        exe_path: str = "./build_msvc/Release/MAIN_1_CoAnQi.exe",
        timeout: float = 60.0,
        verbose: bool = False
    ):
        """
        Initialize CoAnQi calculator interface
        
        Args:
            exe_path: Path to MAIN_1_CoAnQi.exe (default: ./build_msvc/Release/MAIN_1_CoAnQi.exe)
            timeout: Maximum execution time in seconds (default: 60.0)
            verbose: Enable verbose output (default: False)
        
        Raises:
            FileNotFoundError: If C++ calculator executable not found
        """
        self.exe_path = Path(exe_path)
        self.timeout = timeout
        self.verbose = verbose
        
        if not self.exe_path.exists():
            raise FileNotFoundError(
                f"C++ calculator not found: {exe_path}\n"
                f"Build it first with: cmake --build build_msvc --config Release"
            )
        
        if self.verbose:
            print(f"[CoAnQi_Wrapper] Initialized with exe: {self.exe_path}")
    
    def compute_system(
        self, 
        system_name: str, 
        mode: str = "batch"
    ) -> CoAnQiResult:
        """
        Compute UQFF for a given astronomical system
        
        Args:
            system_name: Name of system (must exist in C++ systems database)
            mode: 'batch' (CLI) or 'interactive' (menu automation)
        
        Returns:
            CoAnQiResult with F_U_Bi_i, g_compressed, and auxiliary forces
        
        Raises:
            ValueError: If mode is invalid
            RuntimeError: If C++ calculator returns error
            subprocess.TimeoutExpired: If computation exceeds timeout
        """
        start_time = time.time()
        
        if mode == "batch":
            result = self._batch_compute(system_name)
        elif mode == "interactive":
            result = self._interactive_compute(system_name)
        else:
            raise ValueError(f"Unknown mode: {mode}. Use 'batch' or 'interactive'")
        
        result.computation_time = time.time() - start_time
        
        if self.verbose:
            print(f"[CoAnQi_Wrapper] Computed {system_name} in {result.computation_time:.2f}s")
        
        return result
    
    def _batch_compute(self, system_name: str) -> CoAnQiResult:
        """
        Execute via CLI --batch flag (requires modification to C++ main())
        
        Args:
            system_name: Astronomical system name
        
        Returns:
            CoAnQiResult with computation results
        
        Raises:
            RuntimeError: If C++ calculator returns error
        """
        try:
            result = subprocess.run(
                [str(self.exe_path), "--batch", system_name],
                capture_output=True,
                text=True,
                timeout=self.timeout
            )
            
            if result.returncode != 0:
                return CoAnQiResult(
                    system_name=system_name,
                    F_U_Bi_i=0.0,
                    g_compressed=0.0,
                    status="error",
                    error_message=f"C++ calculator error (code {result.returncode}): {result.stderr}"
                )
            
            # Parse JSON output
            try:
                data = json.loads(result.stdout)
                return CoAnQiResult(
                    system_name=data['system_name'],  # Fixed: C++ outputs 'system_name' not 'system'
                    F_U_Bi_i=data['F_U_Bi_i'],
                    g_compressed=data['g_compressed'],
                    F_jet_rel=data.get('F_jet_rel'),
                    E_acc_rel=data.get('E_acc_rel'),
                    F_drag_rel=data.get('F_drag_rel'),
                    F_gw_rel=data.get('F_gw_rel'),
                    status="success"
                )
            except json.JSONDecodeError as e:
                return CoAnQiResult(
                    system_name=system_name,
                    F_U_Bi_i=0.0,
                    g_compressed=0.0,
                    status="error",
                    error_message=f"JSON parse error: {e}\nOutput: {result.stdout[:200]}"
                )
        
        except subprocess.TimeoutExpired:
            return CoAnQiResult(
                system_name=system_name,
                F_U_Bi_i=0.0,
                g_compressed=0.0,
                status="timeout",
                error_message=f"Computation exceeded {self.timeout}s timeout"
            )
    
    def _interactive_compute(self, system_name: str) -> CoAnQiResult:
        """
        Execute via automated menu navigation (works with current C++ code)
        
        Args:
            system_name: Astronomical system name
        
        Returns:
            CoAnQiResult with computation results
        
        Note:
            This method automates the interactive menu:
            1. Select menu option 1 (Calculate system single)
            2. Navigate category selection (simplified)
            3. Select system by name
            4. Parse console output
            5. Exit menu (option 18)
        """
        # Menu navigation commands
        commands = [
            "1",           # Menu option 1: Calculate system (single)
            system_name,   # System selection (simplified - assumes direct name match)
            "18"           # Exit
        ]
        
        try:
            process = subprocess.Popen(
                [str(self.exe_path)],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            stdout, stderr = process.communicate(
                input="\n".join(commands), 
                timeout=self.timeout
            )
            
            # Parse console output
            return self._parse_console_output(stdout, system_name)
        
        except subprocess.TimeoutExpired:
            process.kill()
            return CoAnQiResult(
                system_name=system_name,
                F_U_Bi_i=0.0,
                g_compressed=0.0,
                status="timeout",
                error_message=f"Interactive computation exceeded {self.timeout}s timeout"
            )
    
    def _parse_console_output(self, output: str, system_name: str) -> CoAnQiResult:
        """
        Extract computation results from console output
        
        Args:
            output: Raw stdout from C++ calculator
            system_name: Astronomical system name
        
        Returns:
            CoAnQiResult with parsed values
        
        Example console output:
            === RESULTS: Sagittarius A* ===
            F_U_Bi_i:           1.234567e+45 N
            g_compressed:       9.876543e+12 m/s^2
            F_jet_rel:          3.456789e+43 N
        """
        result = CoAnQiResult(
            system_name=system_name, 
            F_U_Bi_i=0.0, 
            g_compressed=0.0
        )
        
        lines = output.split('\n')
        for line in lines:
            # Parse F_U_Bi_i
            if "F_U_Bi_i:" in line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        result.F_U_Bi_i = float(parts[1])
                    except ValueError:
                        pass
            
            # Parse g_compressed
            elif "g_compressed:" in line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        result.g_compressed = float(parts[1])
                    except ValueError:
                        pass
            
            # Parse auxiliary forces
            elif "F_jet_rel:" in line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        result.F_jet_rel = float(parts[1])
                    except ValueError:
                        pass
            
            elif "E_acc_rel:" in line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        result.E_acc_rel = float(parts[1])
                    except ValueError:
                        pass
            
            elif "F_drag_rel:" in line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        result.F_drag_rel = float(parts[1])
                    except ValueError:
                        pass
            
            elif "F_gw_rel:" in line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        result.F_gw_rel = float(parts[1])
                    except ValueError:
                        pass
        
        # Validate results
        if result.F_U_Bi_i == 0.0 and result.g_compressed == 0.0:
            result.status = "error"
            result.error_message = "Failed to parse results from console output"
        else:
            result.status = "success"
        
        return result
    
    def list_available_systems(self) -> List[str]:
        """
        Query C++ calculator for available systems
        
        Returns:
            List of system names
        
        Note:
            Requires --list-systems CLI flag (TO BE IMPLEMENTED in C++)
        """
        # TODO: Add --list-systems flag to MAIN_1_CoAnQi.cpp
        raise NotImplementedError(
            "list_available_systems() requires --list-systems CLI flag\n"
            "Add to MAIN_1_CoAnQi.cpp main() function"
        )
    
    def batch_compute_all(
        self, 
        system_names: List[str],
        parallel: bool = False
    ) -> List[CoAnQiResult]:
        """
        Compute UQFF for multiple systems
        
        Args:
            system_names: List of system names
            parallel: Use parallel processing (not yet implemented)
        
        Returns:
            List of CoAnQiResult for each system
        """
        results = []
        
        for system_name in system_names:
            if self.verbose:
                print(f"[CoAnQi_Wrapper] Computing {system_name}...")
            
            result = self.compute_system(system_name)
            results.append(result)
        
        return results
    
    def export_results_csv(
        self, 
        results: List[CoAnQiResult], 
        output_file: str = "uqff_results.csv"
    ):
        """
        Export computation results to CSV
        
        Args:
            results: List of CoAnQiResult objects
            output_file: Output CSV filename
        """
        import csv
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(
                f, 
                fieldnames=[
                    'system_name', 'F_U_Bi_i', 'g_compressed',
                    'F_jet_rel', 'E_acc_rel', 'F_drag_rel', 'F_gw_rel',
                    'computation_time', 'status', 'error_message'
                ]
            )
            writer.writeheader()
            
            for result in results:
                writer.writerow(result.to_dict())
        
        if self.verbose:
            print(f"[CoAnQi_Wrapper] Exported {len(results)} results to {output_file}")


# ============================================================================
# INTEGRATION WITH EXISTING PYTHON DATA LAYER
# ============================================================================

def integrate_with_qcalc():
    """
    Example integration with QCalc.py UnifiedFieldSolver
    
    Demonstrates dual-engine architecture:
    - Python UQFF (QCalc.py) for 8 Master Equations
    - C++ UQFF (MAIN_1_CoAnQi.cpp) for 492 PhysicsTerms
    - Cross-validation for scientific integrity
    """
    try:
        from IPData import InputParameters, InputDataStore
        from QCalc import UnifiedFieldSolver, ComputeParams
        from OPData import OutputDataStore, QueryResult
    except ImportError:
        print("WARNING: Python data layer not found (IPData, QCalc, OPData)")
        print("This is a demonstration function - full integration requires data layer")
        return None
    
    # Initialize both calculators
    python_solver = UnifiedFieldSolver()
    cxx_calculator = CoAnQiCalculator(verbose=True)
    
    # Example system
    system_name = "Sagittarius A*"
    params = ComputeParams(system_name=system_name)
    
    # Python UQFF
    print("\n=== Python UQFF Computation ===")
    python_result = python_solver.solve(params)
    print(f"Python F_U_Bi_i: {python_result.F_U_Bi_i}")
    
    # C++ UQFF
    print("\n=== C++ UQFF Computation ===")
    cxx_result = cxx_calculator.compute_system(system_name)
    print(f"C++ F_U_Bi_i: {cxx_result.F_U_Bi_i} N")
    print(f"C++ g_compressed: {cxx_result.g_compressed} m/s²")
    
    # Cross-validate
    print("\n=== Cross-Validation ===")
    if python_result.F_U_Bi_i != 0:
        relative_error = abs(
            python_result.F_U_Bi_i - cxx_result.F_U_Bi_i
        ) / abs(python_result.F_U_Bi_i)
        print(f"Relative Error: {relative_error * 100:.2f}%")
        
        if relative_error < 0.01:
            print("✅ Cross-validation PASSED (error < 1%)")
        else:
            print("⚠️ Cross-validation WARNING (error >= 1%)")
    
    return {
        'python': python_result,
        'cxx': cxx_result,
        'cross_validation': relative_error < 0.01 if python_result.F_U_Bi_i != 0 else False
    }


# ============================================================================
# COMMAND-LINE INTERFACE
# ============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Python wrapper for MAIN_1_CoAnQi.cpp C++ calculator"
    )
    parser.add_argument(
        "system", 
        type=str, 
        help="Astronomical system name (e.g., 'Sagittarius A*')"
    )
    parser.add_argument(
        "--exe", 
        type=str, 
        default="./build_msvc/Release/MAIN_1_CoAnQi.exe",
        help="Path to MAIN_1_CoAnQi.exe"
    )
    parser.add_argument(
        "--mode", 
        type=str, 
        default="batch",
        choices=["batch", "interactive"],
        help="Execution mode (batch or interactive)"
    )
    parser.add_argument(
        "--timeout", 
        type=float, 
        default=60.0,
        help="Timeout in seconds (default: 60.0)"
    )
    parser.add_argument(
        "--verbose", 
        action="store_true",
        help="Enable verbose output"
    )
    parser.add_argument(
        "--json", 
        action="store_true",
        help="Output as JSON"
    )
    
    args = parser.parse_args()
    
    # Initialize calculator
    try:
        calc = CoAnQiCalculator(
            exe_path=args.exe,
            timeout=args.timeout,
            verbose=args.verbose
        )
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        exit(1)
    
    # Compute system
    result = calc.compute_system(args.system, mode=args.mode)
    
    # Output results
    if args.json:
        print(result.to_json())
    else:
        print(f"\n=== UQFF Computation Results: {result.system_name} ===")
        print(f"Status: {result.status}")
        if result.status == "success":
            print(f"F_U_Bi_i:       {result.F_U_Bi_i:.6e} N")
            print(f"g_compressed:   {result.g_compressed:.6e} m/s²")
            if result.F_jet_rel is not None:
                print(f"F_jet_rel:      {result.F_jet_rel:.6e} N")
            if result.E_acc_rel is not None:
                print(f"E_acc_rel:      {result.E_acc_rel:.6e} J")
            if result.F_drag_rel is not None:
                print(f"F_drag_rel:     {result.F_drag_rel:.6e} N")
            if result.F_gw_rel is not None:
                print(f"F_gw_rel:       {result.F_gw_rel:.6e} N")
            print(f"Execution Time: {result.computation_time:.2f}s")
        else:
            print(f"Error: {result.error_message}")
    
    exit(0 if result.status == "success" else 1)
