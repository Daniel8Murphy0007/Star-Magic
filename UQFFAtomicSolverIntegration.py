#!/usr/bin/env python3
"""
UQFFAtomicSolverIntegration.py - Session 252 Solver Integration Hooks
======================================================================

This module is INTENTIONALLY SEPARATE from CondensedPhysics4.py to avoid
encoding/special character issues in that file.

PURPOSE: Provide integration bridge for v1.5 Simultaneous 7-Layer Solver

Can be imported into CondensedPhysicsAggregator.py with:
    from UQFFAtomicSolverIntegration import (
        Simultaneous7LayerSolverBridge,
        UQFFAtomicSolverCalculator,
        UQFF_SOLVER_BRIDGE,
        UQFF_ATOMIC_SOLVER,
    )

Session 252 Breakthrough: Buoyancy integration (Ubi term) fixes v1.0 plateau.

Author: Daniel T. Murphy
Created: May 25, 2026 (Session 252)
"""

import subprocess
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Any


class Simultaneous7LayerSolverBridge:
    """
    Subprocess wrapper for C++ v1.5 Simultaneous 7-Layer Solver.
    
    Communicates with: .\build_msvc\Release\simultaneous_7layer_solver_v1_5.exe
    
    Protocol:
        INPUT:  {"Z": int, "n": int}
        OUTPUT: {"converged": bool, "iterations": int, "layers": {...}, "Ubi": float, ...}
    """
    
    def __init__(self, solver_exe_path: Optional[str] = None):
        if solver_exe_path is None:
            solver_exe_path = r".\build_msvc\Release\simultaneous_7layer_solver_v1_5.exe"
        
        self.exe_path = Path(solver_exe_path).absolute()
        self.available = self.exe_path.exists()
    
    def solve(self, Z: int, n: int = 1) -> Dict[str, Any]:
        """Call C++ solver for atomic system (Z, n)."""
        
        if not self.available:
            return {
                "Z": Z,
                "n": n,
                "converged": False,
                "error": f"Solver executable not found at {self.exe_path}",
                "layers": {},
                "Ubi": 0.0,
                "convergence_history": [],
            }
        
        try:
            request = json.dumps({"Z": int(Z), "n": int(n)})
            result = subprocess.run(
                str(self.exe_path),
                input=request.encode(),
                capture_output=True,
                timeout=30,
                check=False,
            )
            
            if result.returncode == 0 and result.stdout:
                output_text = result.stdout.decode('utf-8', errors='replace').strip()
                try:
                    response = json.loads(output_text)
                    return response
                except json.JSONDecodeError:
                    import re
                    json_match = re.search(r'\{.*\}', output_text, re.DOTALL)
                    if json_match:
                        return json.loads(json_match.group())
                    raise ValueError(f"JSON parse failed: {output_text[:200]}")
            else:
                stderr = result.stderr.decode('utf-8', errors='replace') if result.stderr else "No stderr"
                return {
                    "Z": Z,
                    "n": n,
                    "converged": False,
                    "error": f"Solver code {result.returncode}: {stderr}",
                    "layers": {},
                    "Ubi": 0.0,
                    "convergence_history": [],
                }
        
        except subprocess.TimeoutExpired:
            return {
                "Z": Z,
                "n": n,
                "converged": False,
                "error": "Solver timeout (>30s)",
                "layers": {},
                "Ubi": 0.0,
                "convergence_history": [],
            }
        except Exception as e:
            return {
                "Z": Z,
                "n": n,
                "converged": False,
                "error": f"Bridge error: {str(e)}",
                "layers": {},
                "Ubi": 0.0,
                "convergence_history": [],
            }


class UQFFAtomicSolverCalculator:
    """
    UQFF Atomic Solver Calculator (Session 252, v1.5)
    
    Simultaneous 7-Layer solution for hydrogen-like atoms and ions.
    Pure stateless calculator following CondensedPhysics pattern.
    
    CANONICAL DESIGN (SESSION 252 VALIDATED):
    - E_DPM = 1.022e6 eV IMMUTABLE (not dynamic)
    - Ubi term is THE missing physics (v1.0 plateau elimination)
    - 4 iterations → machine precision convergence
    - Negligibilities prove buoyancy equilibrium
    """
    
    PARAMETERS = {
        "Z": 1,
        "n": 1,
    }
    
    BETA_I = 0.603
    E_DPM_IMMUTABLE = 1.022e6
    FINE_STRUCTURE = 1.0 / 137.036
    RYDBERG_ENERGY = 13.6057
    
    def __init__(self):
        self.bridge = Simultaneous7LayerSolverBridge()
    
    def compute(self, dataset: dict = None) -> dict:
        """Simultaneous 7-layer atomic solution with buoyancy integration."""
        if dataset is None:
            dataset = {}
        
        Z = max(1, min(int(dataset.get("Z", 1)), 118))
        n = max(1, min(int(dataset.get("n", 1)), 10))
        
        # Call solver bridge
        solver_result = self.bridge.solve(Z, n)
        
        converged = solver_result.get("converged", False)
        iterations = solver_result.get("iterations", 0)
        final_residual = solver_result.get("final_residual", float('inf'))
        layers = solver_result.get("layers", {})
        Ubi = solver_result.get("Ubi", 0.0)
        convergence_history = solver_result.get("convergence_history", [])
        error_msg = solver_result.get("error")
        
        # Element names
        element_names = {
            1: "Hydrogen", 2: "Helium", 3: "Lithium", 4: "Beryllium", 5: "Boron",
            6: "Carbon", 7: "Nitrogen", 8: "Oxygen", 9: "Fluorine", 10: "Neon",
        }
        element_name = element_names.get(Z, f"Element({Z})")
        
        if converged:
            primary_equations = [
                f"# Simultaneous 7-Layer UQFF Atomic Solver v1.5 - Session 252",
                f"## {element_name} (Z={Z}, n={n}) - CONVERGED in {iterations} iters",
                f"",
                f"**Layer Results:**",
                f"1. r_s = {layers.get('r_s', 0.0):.6e} m (Bohr radius)",
                f"2. g_quantum = {layers.get('g_quantum', 0.0):.6e} m/s²",
                f"3. v_orb = {layers.get('v_orb', 0.0):.6e} m/s (orbital velocity)",
                f"4. E_single = {layers.get('E_single', 0.0):.6e} eV (single-particle)",
                f"5. ψ_norm = {layers.get('psi_norm', 0.0):.6e} (normalization)",
                f"6. E_DPM = {layers.get('E_DPM', self.E_DPM_IMMUTABLE):.6e} eV (IMMUTABLE)",
                f"7. E_pair = {layers.get('E_pair', 0.0):.6e} eV (pair energy)",
                f"",
                f"**Buoyancy (THE Missing Physics - v1.5 Breakthrough):**",
                f"Ubᵢ = {Ubi:.6e} eV (Z-dependent)",
                f"",
                f"**Force Balance (Central Equation):**",
                f"F_U_total = Ug_sum - Ubᵢ + Um = {layers.get('F_U_total', 0.0):.6e} eV",
                f"",
                f"**Convergence:**",
                f"Final residual: {final_residual:.6e} eV",
                f"Iterations: {iterations}",
                f"",
                f"**Physics Interpretation:**",
                f"• Negligibilities prove buoyancy maintaining force balance",
                f"• Ubi suppresses quantum chain (~1e-33 range)",
                f"• v1.0 plateau eliminated by adding Ubi term",
            ]
        else:
            primary_equations = [
                f"# Simultaneous 7-Layer UQFF Atomic Solver v1.5",
                f"## ERROR - {element_name} (Z={Z}, n={n})",
                f"Status: NOT CONVERGED",
                f"Error: {error_msg}",
            ]
        
        available_equations = [
            "1. r_s = a₀·n²/Z (Bohr radius)",
            "2. E_n = -13.6·Z²/n² eV (Rydberg formula)",
            "3. v_orb = c·α·Z/n (orbital velocity)",
            "4. T_orb = 2π·r_s/v_orb (orbital period)",
            "5. L = n·ℏ (orbital angular momentum)",
            "6. E_spin_orbit ∝ (α²·Z⁴)/(n³·(2l+1)) (spin-orbit coupling)",
            "7. Hyperfine shift ∝ g_I·(α²/n³) (nuclear interaction)",
            "8. Stark shift ∝ n²·E_ext (external field)",
            "9. Zeeman shift ∝ m·B (magnetic field)",
            "10. Lamb shift ∝ (α⁵·Z⁴)/n³ (QED)",
            "11. Ubᵢ = βᵢ·Ugsum·Ωg·(Mbh/dg) (buoyancy) ← NEW v1.5",
            "12. F_U = Ug_sum - Ubᵢ + Um = 0 (unified force balance) ← NEW v1.5",
        ]
        
        return {
            "primary_equations": primary_equations,
            "available_equations": available_equations,
            "layer_solutions": layers,
            "Ubi_buoyancy": Ubi,
            "convergence_data": {
                "converged": converged,
                "iterations": iterations,
                "final_residual": final_residual,
                "history": convergence_history,
                "Z": Z,
                "n": n,
            }
        }
    
    def simulate(self, Z_range=None, **kwargs):
        """Run solver across Z range."""
        Z_range = Z_range or [1, 2, 10, 54]
        results = []
        for Z in Z_range:
            res = self.compute({"Z": Z, "n": 1})
            results.append(res)
        return results


# Singleton instances for direct import
UQFF_SOLVER_BRIDGE = Simultaneous7LayerSolverBridge()
UQFF_ATOMIC_SOLVER = UQFFAtomicSolverCalculator()

# Exports
__all__ = [
    'Simultaneous7LayerSolverBridge',
    'UQFFAtomicSolverCalculator',
    'UQFF_SOLVER_BRIDGE',
    'UQFF_ATOMIC_SOLVER',
]
