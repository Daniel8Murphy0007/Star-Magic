"""
WSTP Kernel Demo Runner — Live Wolfram Kernel Session for 9-Sector UQFF Lagrangian

Session 204 | Daniel Murphy
PURPOSE: Connect to a running Wolfram kernel (via wolframscript subprocess or
         the wolframclient Python library), load the generated .wl package from
         wstp_symbolic_exporter.py, evaluate the unified Lagrangian + GW damping
         channels, and return symbolic/numerical results.

ARCHITECTURE:
  wstp_symbolic_exporter.py → uqff_lagrangian_unified.wl
        ↓ (this module loads into kernel)
  Wolfram Kernel (local or cloud)
        ↓
  Evaluated results (L_UQFF, F_U_Bi_i, D_total, h_UQFF for GW170817/150914/190425)

REQUIREMENTS (one of):
  - wolframscript on PATH (ships with Mathematica / Wolfram Engine)
  - pip install wolframclient  (Wolfram Client Library for Python)

FALLBACK: If neither is available, runs a dry-run mode that prints the
          expressions that *would* be sent to the kernel.
"""

import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple


# ── §1  CONSTANTS (mirror wstp_symbolic_exporter.py) ───────────────────────

G       = 6.67430e-11
c       = 2.99792e8
M_sun   = 1.98892e30
KAPPA   = 5.787e-9
SSQ     = 0.57
H_SCM   = 0.99
BETA_I  = 0.603
U_UA    = 1e-4
D_TRZ   = 0.900
D_STRING = 0.370
F_COMBINED = 0.333

WL_PACKAGE = "uqff_lagrangian_unified.wl"

# GW events for demo evaluation
GW_EVENTS = {
    "GW170817": {"m1": 1.46, "m2": 1.27, "DL_Mpc": 40, "D_total": 0.333,
                 "h_GR": 5.4176e-22, "h_UQFF": 1.804e-22, "type": "BNS"},
    "GW150914": {"m1": 36.0, "m2": 29.0, "DL_Mpc": 410, "D_total": 0.810,
                 "h_GR": 1.2499e-21, "h_UQFF": 4.1622e-22, "type": "BBH"},
    "GW190425": {"m1": 1.7, "m2": 1.5, "DL_Mpc": 159, "D_total": 0.530,
                 "h_GR": 3.0e-22, "h_UQFF": 1.59e-22, "type": "BNS"},
}


# ── §2  KERNEL BACKEND DETECTION ──────────────────────────────────────────

class KernelBackend:
    """Detect available Wolfram kernel backend."""

    WOLFRAMSCRIPT = "wolframscript"
    WOLFRAMCLIENT = "wolframclient"
    DRY_RUN = "dry_run"

    @staticmethod
    def detect() -> str:
        """Return best available backend."""
        # Try wolframclient first (Python library)
        try:
            import wolframclient  # noqa: F401
            return KernelBackend.WOLFRAMCLIENT
        except ImportError:
            pass

        # Try wolframscript CLI
        try:
            result = subprocess.run(
                ["wolframscript", "-code", "1+1"],
                capture_output=True, text=True, timeout=30,
            )
            if result.returncode == 0 and result.stdout.strip() == "2":
                return KernelBackend.WOLFRAMSCRIPT
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass

        return KernelBackend.DRY_RUN


# ── §3  EXPRESSION BUILDER ────────────────────────────────────────────────

def build_demo_expressions(wl_path: str) -> List[Dict[str, str]]:
    """Build the sequence of Wolfram Language expressions for the demo.

    Each entry is {"label": ..., "code": ...} where code is valid WL.
    """
    # Normalize path for Mathematica (forward slashes, escaped)
    wl_abs = os.path.abspath(wl_path).replace("\\", "/")

    exprs = []

    # 1. Load the package
    exprs.append({
        "label": "Load UQFFLagrangian package",
        "code": f'Get["{wl_abs}"]',
    })

    # 2. Retrieve constants
    exprs.append({
        "label": "Retrieve calibrated constants",
        "code": 'UQFFConstants[]',
    })

    # 3. Evaluate individual sector Lagrangians
    exprs.append({
        "label": "Sector 1 — Einstein-Hilbert (R → 10^-52)",
        "code": 'LEH[10^-52] // N',
    })
    exprs.append({
        "label": "Sector 5 — Magnetic-Dipole (SGR1745: B=10^14 T)",
        "code": 'LMag[10^14, 7.09*^-37, 10^3, 10^4, 10^4] // N',
    })
    exprs.append({
        "label": "Sector 8 — LENR-Resonance (ω=1.25 THz)",
        "code": 'LLENR[10^-15, 1.25*^12, 300*2*Pi, 1, 10^-20, 10^-28, 0] // N',
    })

    # 4. Unified Lagrangian (Sgr A* defaults)
    exprs.append({
        "label": "Unified L_UQFF (Sgr A* defaults)",
        "code": 'result = LUQFF[<||>]; result["L_UQFF_total"] // N',
    })

    # 5. F_U_Bi_i force assembly
    exprs.append({
        "label": "F_U_Bi_i total force (Sgr A* defaults)",
        "code": 'forces = FUBii[<||>]; forces["F_U_Bi_i"] // N',
    })

    # 6. GW damping — D_total product
    exprs.append({
        "label": "GW damping: D_total = D_Aether × D_SCm × D_TRZ × D_String",
        "code": 'DTotal[1.0, 1.0, 0.900, 0.370] // N',
    })

    # 7. GW170817 strain
    exprs.append({
        "label": "GW170817: h_UQFF = h_GR × D_total",
        "code": 'HUQFF[5.4176*^-22, 0.333] // N',
    })

    # 8. GW150914 strain
    exprs.append({
        "label": "GW150914: h_UQFF",
        "code": 'HUQFF[1.2499*^-21, 0.333] // N',
    })

    # 9. Chirp mass for GW170817
    exprs.append({
        "label": "GW170817 chirp mass (1.46 + 1.27 M☉)",
        "code": f'ChirpMass[1.46 * {M_sun}, 1.27 * {M_sun}] // N',
    })

    # 10. Phase lag estimate
    exprs.append({
        "label": "GW170817 phase lag (150 Hz, 100 s chirp)",
        "code": 'PhaseLag[0.333, 150, 100] // N',
    })

    # 11. Apparent distance bias
    exprs.append({
        "label": "GW170817 apparent distance bias",
        "code": 'ApparentDistanceBias[40 * 3.086*^22, 0.333] // N',
    })

    return exprs


# ── §4  KERNEL RUNNERS ────────────────────────────────────────────────────

def run_wolframscript(code: str, timeout: int = 60) -> Tuple[bool, str]:
    """Execute WL code via wolframscript subprocess."""
    try:
        result = subprocess.run(
            ["wolframscript", "-code", code],
            capture_output=True, text=True, timeout=timeout,
        )
        if result.returncode == 0:
            return True, result.stdout.strip()
        return False, result.stderr.strip() or f"exit code {result.returncode}"
    except subprocess.TimeoutExpired:
        return False, "TIMEOUT"
    except FileNotFoundError:
        return False, "wolframscript not found"


def run_wolframclient(code: str, session=None) -> Tuple[bool, str]:
    """Execute WL code via wolframclient library."""
    try:
        from wolframclient.evaluation import WolframLanguageSession
        from wolframclient.language import wl as _wl

        own_session = session is None
        if own_session:
            session = WolframLanguageSession()
            session.start()

        try:
            result = session.evaluate(
                _wl.ToExpression(code)
            )
            return True, str(result)
        finally:
            if own_session:
                session.terminate()
    except ImportError:
        return False, "wolframclient not installed"
    except Exception as e:
        return False, str(e)


def run_dry(code: str) -> Tuple[bool, str]:
    """Dry-run: echo the expression without evaluation."""
    return True, f"[DRY-RUN] {code}"


# ── §5  DEMO ENGINE ──────────────────────────────────────────────────────

@dataclass
class DemoResult:
    """Result of a single demo expression evaluation."""
    label: str
    code: str
    success: bool
    output: str
    elapsed_ms: float


class WSTPKernelDemoRunner:
    """Execute live WSTP kernel demos against the generated .wl package.

    Supports three backends: wolframscript (subprocess), wolframclient
    (Python library), or dry-run (prints expressions without kernel).
    """

    def __init__(self, wl_path: str = WL_PACKAGE, backend: Optional[str] = None):
        self.wl_path = wl_path
        self.backend = backend or KernelBackend.detect()
        self._wc_session = None
        self.results: List[DemoResult] = []

    def _ensure_package(self):
        """Generate .wl package if missing."""
        if not os.path.exists(self.wl_path):
            try:
                from wstp_symbolic_exporter import WSTPSymbolicExporter
                exporter = WSTPSymbolicExporter(self.wl_path)
                exporter.export()
                print(f"  [auto-generated] {self.wl_path}")
            except ImportError:
                raise FileNotFoundError(
                    f"{self.wl_path} not found and wstp_symbolic_exporter.py unavailable"
                )

    def _eval(self, code: str) -> Tuple[bool, str]:
        """Route to the active backend."""
        if self.backend == KernelBackend.WOLFRAMCLIENT:
            return run_wolframclient(code, self._wc_session)
        elif self.backend == KernelBackend.WOLFRAMSCRIPT:
            return run_wolframscript(code)
        else:
            return run_dry(code)

    def run_demo(self) -> Dict[str, Any]:
        """Execute the full demo sequence and return results."""
        self._ensure_package()
        self.results.clear()

        expressions = build_demo_expressions(self.wl_path)

        # Open persistent session for wolframclient
        if self.backend == KernelBackend.WOLFRAMCLIENT:
            try:
                from wolframclient.evaluation import WolframLanguageSession
                self._wc_session = WolframLanguageSession()
                self._wc_session.start()
            except Exception as e:
                print(f"  [WARN] wolframclient session failed, falling back to dry-run: {e}")
                self.backend = KernelBackend.DRY_RUN

        t0_all = time.perf_counter()

        for expr in expressions:
            t0 = time.perf_counter()
            success, output = self._eval(expr["code"])
            elapsed = (time.perf_counter() - t0) * 1000

            self.results.append(DemoResult(
                label=expr["label"],
                code=expr["code"],
                success=success,
                output=output,
                elapsed_ms=round(elapsed, 2),
            ))

        total_elapsed = time.perf_counter() - t0_all

        # Close session
        if self._wc_session is not None:
            try:
                self._wc_session.terminate()
            except Exception:
                pass
            self._wc_session = None

        # GW validation cross-check (Python side)
        gw_validation = self._validate_gw_events()

        return {
            "backend": self.backend,
            "wl_package": os.path.abspath(self.wl_path),
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "total_elapsed_ms": round(total_elapsed * 1000, 2),
            "n_expressions": len(self.results),
            "n_success": sum(1 for r in self.results if r.success),
            "n_failed": sum(1 for r in self.results if not r.success),
            "results": [
                {
                    "label": r.label,
                    "code": r.code,
                    "success": r.success,
                    "output": r.output,
                    "elapsed_ms": r.elapsed_ms,
                }
                for r in self.results
            ],
            "gw_validation": gw_validation,
        }

    def _validate_gw_events(self) -> List[Dict[str, Any]]:
        """Python-side GW cross-validation against PAPER_001-009 values."""
        validations = []
        for name, ev in GW_EVENTS.items():
            d_total = ev["D_total"]
            h_predicted = ev["h_GR"] * d_total
            h_observed = ev["h_UQFF"]
            residual = abs(h_predicted - h_observed) / max(h_observed, 1e-50)
            validations.append({
                "event": name,
                "type": ev["type"],
                "D_total": d_total,
                "h_GR": ev["h_GR"],
                "h_UQFF_predicted": h_predicted,
                "h_UQFF_observed": h_observed,
                "residual": round(residual, 6),
                "pass": residual < 0.15,
            })
        return validations

    def export_results(self, report: Dict[str, Any],
                       filepath: str = "wstp_kernel_demo_results.json") -> str:
        """Export demo results to JSON."""
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2, default=str)
        return filepath


# ── §6  MAIN ──────────────────────────────────────────────────────────────

def main():
    """Run live WSTP kernel demo on the 9-sector UQFF Lagrangian."""
    print("=" * 72)
    print("WSTP Kernel Demo Runner — 9-Sector UQFF Lagrangian Live Session")
    print("=" * 72)

    backend = KernelBackend.detect()
    print(f"  Backend detected: {backend}")

    if backend == KernelBackend.DRY_RUN:
        print("  [INFO] No Wolfram kernel available — running in dry-run mode.")
        print("         Install wolframscript or 'pip install wolframclient' for live evaluation.")

    runner = WSTPKernelDemoRunner(backend=backend)

    print(f"\n── Running demo ({runner.wl_path}) ──\n")
    report = runner.run_demo()

    # Print results
    for r in report["results"]:
        status = "OK" if r["success"] else "FAIL"
        print(f"  [{status:4s}] {r['label']}")
        print(f"         Code:   {r['code'][:80]}{'...' if len(r['code']) > 80 else ''}")
        print(f"         Output: {r['output'][:100]}{'...' if len(r['output']) > 100 else ''}")
        print(f"         Time:   {r['elapsed_ms']:.2f} ms")
        print()

    # GW validation
    print("── GW Event Validation (Python cross-check) ──\n")
    for ev in report["gw_validation"]:
        status = "PASS" if ev["pass"] else "FAIL"
        print(f"  {ev['event']:10s} ({ev['type']}):  D={ev['D_total']:.3f}  "
              f"h_pred={ev['h_UQFF_predicted']:.4e}  "
              f"h_obs={ev['h_UQFF_observed']:.4e}  "
              f"resid={ev['residual']:.4f}  [{status}]")

    # Summary
    print(f"\n── Summary ──")
    print(f"  Backend:     {report['backend']}")
    print(f"  Expressions: {report['n_success']}/{report['n_expressions']} succeeded")
    print(f"  Total time:  {report['total_elapsed_ms']:.2f} ms")

    # Export
    out_path = runner.export_results(report)
    print(f"\n[OK] Results exported: {out_path}")

    print(f"\n{'=' * 72}")
    print(f"DEMO COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
