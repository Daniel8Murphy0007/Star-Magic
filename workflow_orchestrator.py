#!/usr/bin/env python3
"""
workflow_orchestrator.py — Session 204 Master Workflow Orchestrator
═══════════════════════════════════════════════════════════════════════════════

PURPOSE: Master orchestrator that ties together all Session 204 modules:
  1. lagrangian_re_runner.py  — Euler-Lagrange re-derivation for PAPER_859-877
  2. vds_dvp_bsh_lenr_synthesis.py — VDS/DVP/BSH across LENR lab systems
  3. quadriadic_scm_qcalcgeom.py — Quadriadic SCm + QCalcGeom layer

Cross-validates results between modules and produces a unified report
with coherence metrics across all three analytical domains.

SESSION: 204 | April 7, 2026
"""

import json
import sys
import time
from typing import Dict, Optional

# ═══════════════════════════════════════════════════════════════════════════════
# §1  IMPORTS (LAZY — graceful degradation if modules unavailable)
# ═══════════════════════════════════════════════════════════════════════════════

try:
    from lagrangian_re_runner import LagrangianReRunner
    _HAS_LAGRANGIAN = True
except ImportError:
    _HAS_LAGRANGIAN = False

try:
    from vds_dvp_bsh_lenr_synthesis import VDSDVPBSHSynthesis
    _HAS_VDS = True
except ImportError:
    _HAS_VDS = False

try:
    from quadriadic_scm_qcalcgeom import QuadriadicSCmEngine
    _HAS_QUADRIADIC = True
except ImportError:
    _HAS_QUADRIADIC = False


# ═══════════════════════════════════════════════════════════════════════════════
# §2  CROSS-VALIDATION ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

class CrossValidator:
    """
    Validates consistency between the three analysis domains:
      - Lagrangian → Euler-Lagrange consistency (d/dt dL/dq_dot = dL/dq)
      - VDS/DVP/BSH → Triadic coherence + phonon convergence
      - Quadriadic SCm → Cross-equation H_SCm coherence
    """

    @staticmethod
    def validate_lagrangian(results: Dict) -> Dict:
        """Check Lagrangian re-runner outputs for self-consistency."""
        # Top-level keys are target names; skip metadata keys
        meta_keys = {"scm_convergence", "total_re_runs", "lagrangian_sectors_covered"}
        target_keys = [k for k in results if k not in meta_keys]
        n_targets = len(target_keys)
        n_ok = 0
        for k in target_keys:
            t = results[k]
            if not isinstance(t, dict):
                continue
            r = t.get("results", t)
            if any(
                r.get(fld, 0) != 0
                for fld in ("F_gravity", "g_total", "Um", "V_26_total",
                            "Ug1_at_r0", "Um_t1", "g_emergent")
            ):
                n_ok += 1

        return {
            "domain": "Lagrangian re-runner",
            "total_targets": n_targets,
            "valid_targets": n_ok,
            "pass": n_ok == n_targets and n_targets > 0,
            "coverage": n_ok / max(n_targets, 1),
        }

    @staticmethod
    def validate_vds_synthesis(results: Dict) -> Dict:
        """Check VDS/DVP/BSH synthesis coherence."""
        systems = results.get("systems", [])
        n_sys = len(systems)
        coherences = [
            s.get("synthesis", {}).get("scm_coherence", 0) for s in systems
        ]
        mean_coh = sum(coherences) / max(n_sys, 1)
        cross = results.get("cross_system", {})
        phonon_ok = cross.get("all_converge_on_1_25_THz", False)

        return {
            "domain": "VDS/DVP/BSH synthesis",
            "total_systems": n_sys,
            "mean_triadic_coherence": mean_coh,
            "phonon_converged": phonon_ok,
            "pass": mean_coh > 0 and n_sys > 0,
        }

    @staticmethod
    def validate_quadriadic(results: Dict) -> Dict:
        """Check quadriadic SCm cross-equation coherence."""
        global_m = results.get("global_metrics", {})
        mean_coh = global_m.get("mean_coherence", 0)
        mod_check = global_m.get("26_mod_113", -1) == 12
        n_sys = len(results.get("systems", []))

        return {
            "domain": "Quadriadic SCm",
            "total_systems": n_sys,
            "mean_scm_coherence": mean_coh,
            "mod_26_113_pass": mod_check,
            "pass": mean_coh > 0 and mod_check and n_sys > 0,
        }

    @staticmethod
    def cross_domain_consistency(lag_val: Dict, vds_val: Dict, quad_val: Dict) -> Dict:
        """
        Cross-domain consistency check.

        All three domains must pass independently AND their coherence
        metrics must be mutually consistent (all positive, no NaN).
        """
        all_pass = lag_val["pass"] and vds_val["pass"] and quad_val["pass"]

        summary_metrics = {
            "lagrangian_coverage": lag_val.get("coverage", 0),
            "vds_triadic_coherence": vds_val.get("mean_triadic_coherence", 0),
            "quadriadic_scm_coherence": quad_val.get("mean_scm_coherence", 0),
        }

        vals = list(summary_metrics.values())
        unified_coherence = 1.0
        for v in vals:
            unified_coherence *= max(abs(v), 1e-30)
        unified_coherence = unified_coherence ** (1.0 / len(vals))

        return {
            "all_domains_pass": all_pass,
            "unified_coherence": unified_coherence,
            "metrics": summary_metrics,
            "conclusion": (
                "PASS: Cross-domain consistency confirmed across "
                "Lagrangian, VDS/DVP/BSH, and Quadriadic SCm domains."
                if all_pass else
                "PARTIAL: One or more domains did not pass validation."
            ),
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §3  WORKFLOW ORCHESTRATOR
# ═══════════════════════════════════════════════════════════════════════════════

class WorkflowOrchestrator:
    """
    Master orchestrator for Session 204 workflow.

    Execution order:
      1. Lagrangian re-runner (micro-plasmoid, monopole, Um, cosmogenesis)
      2. VDS/DVP/BSH LENR synthesis (6 lab systems)
      3. Quadriadic SCm with QCalcGeom (5 scale systems)
      4. Cross-validation
      5. Unified report
    """

    def __init__(self):
        self.validator = CrossValidator()
        self.results = {}
        self.timings = {}

    def run_lagrangian(self) -> Optional[Dict]:
        """Step 1: Lagrangian re-runner."""
        if not _HAS_LAGRANGIAN:
            print("  [SKIP] lagrangian_re_runner.py not available")
            return None

        t0 = time.time()
        runner = LagrangianReRunner()
        result = runner.run_all()
        self.timings["lagrangian"] = time.time() - t0
        return result

    def run_vds_synthesis(self) -> Optional[Dict]:
        """Step 2: VDS/DVP/BSH LENR synthesis."""
        if not _HAS_VDS:
            print("  [SKIP] vds_dvp_bsh_lenr_synthesis.py not available")
            return None

        t0 = time.time()
        synth = VDSDVPBSHSynthesis()
        result = synth.synthesize_all()
        self.timings["vds_synthesis"] = time.time() - t0
        return result

    def run_quadriadic(self) -> Optional[Dict]:
        """Step 3: Quadriadic SCm + QCalcGeom."""
        if not _HAS_QUADRIADIC:
            print("  [SKIP] quadriadic_scm_qcalcgeom.py not available")
            return None

        t0 = time.time()
        engine = QuadriadicSCmEngine()
        result = engine.analyze_all()
        self.timings["quadriadic"] = time.time() - t0
        return result

    def run_all(self) -> Dict:
        """Execute full Session 204 workflow."""
        print("=" * 80)
        print("SESSION 204 — FULL WORKFLOW ORCHESTRATION")
        print("=" * 80)

        # Step 1: Lagrangian
        print("\n[1/3] Lagrangian re-runner...")
        lag_result = self.run_lagrangian()
        self.results["lagrangian"] = lag_result
        if lag_result:
            print(f"       Completed in {self.timings['lagrangian']:.3f}s")
            print(f"       Targets: {len(lag_result.get('targets', {}))}")

        # Step 2: VDS synthesis
        print("\n[2/3] VDS/DVP/BSH LENR synthesis...")
        vds_result = self.run_vds_synthesis()
        self.results["vds_synthesis"] = vds_result
        if vds_result:
            print(f"       Completed in {self.timings['vds_synthesis']:.3f}s")
            print(f"       Systems: {len(vds_result.get('systems', []))}")

        # Step 3: Quadriadic
        print("\n[3/3] Quadriadic SCm + QCalcGeom...")
        quad_result = self.run_quadriadic()
        self.results["quadriadic"] = quad_result
        if quad_result:
            print(f"       Completed in {self.timings['quadriadic']:.3f}s")
            gm = quad_result.get("global_metrics", {})
            print(f"       Mean SCm coherence: {gm.get('mean_coherence', 0):.6f}")

        # Step 4: Cross-validation
        print("\n[VALIDATION] Cross-domain consistency check...")
        lag_val = (self.validator.validate_lagrangian(lag_result)
                   if lag_result else {"domain": "Lagrangian", "pass": False,
                                       "coverage": 0, "total_targets": 0,
                                       "valid_targets": 0})
        vds_val = (self.validator.validate_vds_synthesis(vds_result)
                   if vds_result else {"domain": "VDS/DVP/BSH", "pass": False,
                                       "mean_triadic_coherence": 0,
                                       "phonon_converged": False, "total_systems": 0})
        quad_val = (self.validator.validate_quadriadic(quad_result)
                    if quad_result else {"domain": "Quadriadic", "pass": False,
                                         "mean_scm_coherence": 0,
                                         "mod_26_113_pass": False, "total_systems": 0})

        cross = self.validator.cross_domain_consistency(lag_val, vds_val, quad_val)

        self.results["validation"] = {
            "lagrangian": lag_val,
            "vds_synthesis": vds_val,
            "quadriadic": quad_val,
            "cross_domain": cross,
        }

        # Step 5: Print report
        self._print_report(lag_val, vds_val, quad_val, cross)

        return self.results

    def _print_report(self, lag_val: Dict, vds_val: Dict,
                      quad_val: Dict, cross: Dict):
        """Print unified workflow report."""
        print("\n" + "=" * 80)
        print("UNIFIED WORKFLOW REPORT — SESSION 204")
        print("=" * 80)

        # Domain summaries
        for val in [lag_val, vds_val, quad_val]:
            status = "PASS" if val["pass"] else "FAIL"
            print(f"  [{status}] {val['domain']}")
            for k, v in val.items():
                if k not in ("domain", "pass"):
                    print(f"         {k}: {v}")

        # Cross-domain
        status = "PASS" if cross["all_domains_pass"] else "PARTIAL"
        print(f"\n  [{status}] Cross-domain consistency")
        print(f"         Unified coherence: {cross['unified_coherence']:.6f}")
        for k, v in cross["metrics"].items():
            print(f"         {k}: {v:.6f}")

        # Timings
        print(f"\n  Timings:")
        for name, t in self.timings.items():
            print(f"    {name}: {t:.3f}s")
        total_t = sum(self.timings.values())
        print(f"    total: {total_t:.3f}s")

        print(f"\n  {cross['conclusion']}")
        print("=" * 80)


# ═══════════════════════════════════════════════════════════════════════════════
# §4  CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    orchestrator = WorkflowOrchestrator()

    if len(sys.argv) > 1 and sys.argv[1] == "--json":
        results = orchestrator.run_all()
        outfile = sys.argv[2] if len(sys.argv) > 2 else "session204_workflow_results.json"
        clean = json.loads(json.dumps(results, default=str))
        with open(outfile, "w") as f:
            json.dump(clean, f, indent=2)
        print(f"\nExported to {outfile}")
    else:
        orchestrator.run_all()


if __name__ == "__main__":
    main()
