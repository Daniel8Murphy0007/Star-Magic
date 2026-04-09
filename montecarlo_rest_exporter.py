"""
Monte Carlo REST Exporter — Cosmogenesis MC → Production REST API Endpoint

Session 204 | Daniel Murphy
PURPOSE: Bridge between cosmogenesis_montecarlo_v2.py (GW-constrained VDS/DVP/BSH
         Monte Carlo) and the production REST API. Provides:
         1. A Python FastAPI/http.server endpoint: POST /api/cosmogenesis/montecarlo
         2. A Node.js route patch for uqff_server.js (port 3141)
         3. JSON export utilities for the MC → REST pipeline

ARCHITECTURE:
  cosmogenesis_montecarlo_v2.py (MC engine)
        ↓ this module wraps
  POST /api/cosmogenesis/montecarlo  (Python HTTP server, port 8443)
        ↓ also generates
  uqff_server_mc_route.js  (drop-in route for uqff_server.js, port 3141)

ENDPOINTS:
  POST /api/cosmogenesis/montecarlo — run MC with optional config override
  GET  /api/cosmogenesis/status     — last run statistics
  GET  /api/cosmogenesis/health     — health check
"""

import json
import os
import sys
import time
import threading
from datetime import datetime, timezone
from http.server import HTTPServer, BaseHTTPRequestHandler
from typing import Any, Dict, Optional
from urllib.parse import urlparse

# Import the MC engine
try:
    from cosmogenesis_montecarlo_v2 import (
        CosmogenesisMonteCarloV2, MCConfig, GW_EVENTS,
        F_COMBINED, D_TRZ, D_STRING, SSQ,
    )
    HAS_MC_ENGINE = True
except ImportError:
    HAS_MC_ENGINE = False


# ── §1  CONFIGURATION ────────────────────────────────────────────────────

MC_PORT = int(os.environ.get("MC_API_PORT", 8443))
MC_HOST = os.environ.get("MC_API_HOST", "127.0.0.1")

# Default MC parameters (can be overridden per request)
DEFAULT_N_SAMPLES = 100_000
DEFAULT_SEED = 42
MAX_N_SAMPLES = 1_000_000  # hard cap for safety


# ── §2  MC RUNNER (thread-safe) ──────────────────────────────────────────

class MCRunner:
    """Thread-safe Monte Carlo runner with result caching."""

    def __init__(self):
        self._lock = threading.Lock()
        self._last_result: Optional[Dict[str, Any]] = None
        self._running = False

    def run(self, config_override: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Execute MC run, return JSON-serializable results."""
        if not HAS_MC_ENGINE:
            return {"error": "cosmogenesis_montecarlo_v2.py not importable"}

        with self._lock:
            if self._running:
                return {"error": "MC run already in progress"}
            self._running = True

        try:
            # Build config from defaults + overrides
            n_samples = DEFAULT_N_SAMPLES
            seed = DEFAULT_SEED
            ssq_range = (0.40, 0.75)
            rho_ratio_range = (1.0, 3.0)
            kappa_range = (1e-10, 1e-7)
            beta_i_range = (0.3, 0.9)
            f_combined_target = F_COMBINED
            f_combined_tol = 0.05

            if config_override:
                n_samples = min(
                    int(config_override.get("n_samples", n_samples)),
                    MAX_N_SAMPLES,
                )
                seed = int(config_override.get("seed", seed))
                if "ssq_range" in config_override:
                    ssq_range = tuple(config_override["ssq_range"])
                if "rho_ratio_range" in config_override:
                    rho_ratio_range = tuple(config_override["rho_ratio_range"])
                if "f_combined_tol" in config_override:
                    f_combined_tol = float(config_override["f_combined_tol"])

            config = MCConfig(
                n_samples=n_samples,
                seed=seed,
                ssq_range=ssq_range,
                rho_ratio_range=rho_ratio_range,
                kappa_range=kappa_range,
                beta_i_range=beta_i_range,
                f_combined_target=f_combined_target,
                f_combined_tol=f_combined_tol,
            )

            engine = CosmogenesisMonteCarloV2(config)
            result = engine.run()

            # Add REST metadata
            result["api_version"] = "1.0.0"
            result["endpoint"] = "/api/cosmogenesis/montecarlo"
            result["timestamp"] = datetime.now(timezone.utc).isoformat()

            with self._lock:
                self._last_result = result

            return result
        finally:
            with self._lock:
                self._running = False

    def get_status(self) -> Dict[str, Any]:
        """Return last run statistics."""
        with self._lock:
            if self._last_result is None:
                return {"status": "no_runs_yet"}
            return {
                "status": "ready",
                "last_run": {
                    "n_accepted": self._last_result.get("n_accepted", 0),
                    "n_rejected": self._last_result.get("n_rejected", 0),
                    "acceptance_rate": self._last_result.get("acceptance_rate", 0),
                    "elapsed_seconds": self._last_result.get("elapsed_seconds", 0),
                    "samples_per_second": self._last_result.get("samples_per_second", 0),
                    "timestamp": self._last_result.get("timestamp", ""),
                },
                "running": self._running,
            }


# Singleton runner
_runner = MCRunner()


# ── §3  HTTP REQUEST HANDLER ─────────────────────────────────────────────

class MCRequestHandler(BaseHTTPRequestHandler):
    """HTTP handler for Monte Carlo REST endpoints."""

    def _send_json(self, status: int, data: Any):
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()
        body = json.dumps(data, indent=2, default=str)
        self.wfile.write(body.encode("utf-8"))

    def _read_body(self) -> Optional[Dict]:
        content_len = int(self.headers.get("Content-Length", 0))
        if content_len == 0:
            return None
        raw = self.rfile.read(content_len)
        try:
            return json.loads(raw.decode("utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError):
            return None

    def do_OPTIONS(self):
        """Handle CORS preflight."""
        self._send_json(204, "")

    def do_GET(self):
        path = urlparse(self.path).path.rstrip("/")

        if path == "/api/cosmogenesis/health":
            self._send_json(200, {
                "status": "ok",
                "service": "montecarlo_rest_exporter",
                "engine_available": HAS_MC_ENGINE,
                "port": MC_PORT,
                "gw_events": list(GW_EVENTS.keys()) if HAS_MC_ENGINE else [],
            })
        elif path == "/api/cosmogenesis/status":
            self._send_json(200, _runner.get_status())
        else:
            self._send_json(404, {"error": f"unknown endpoint: {path}"})

    def do_POST(self):
        path = urlparse(self.path).path.rstrip("/")

        if path == "/api/cosmogenesis/montecarlo":
            body = self._read_body()
            result = _runner.run(config_override=body)
            status = 200 if "error" not in result else 503
            self._send_json(status, result)
        else:
            self._send_json(404, {"error": f"unknown endpoint: {path}"})

    def log_message(self, fmt, *args):
        """Prefix log with timestamp."""
        ts = datetime.now(timezone.utc).strftime("%H:%M:%S")
        print(f"  [{ts}] {fmt % args}")


# ── §4  NODE.JS ROUTE GENERATOR ──────────────────────────────────────────

def generate_node_route(output_path: str = "uqff_server_mc_route.js") -> str:
    """Generate a drop-in Express route for uqff_server.js.

    This JS file spawns cosmogenesis_montecarlo_v2.py as a child process
    and returns JSON results to the REST client.
    """
    js_code = """\
/**
 * Monte Carlo REST route for uqff_server.js
 * Drop-in: add to uqff_server.js routes section
 *
 * Usage:
 *   POST /api/cosmogenesis/montecarlo
 *   Body (optional): {"n_samples": 100000, "seed": 42}
 *
 * Generated by montecarlo_rest_exporter.py
 */

const { execFile } = require('child_process');
const path = require('path');

const MC_SCRIPT = path.join(__dirname, 'cosmogenesis_montecarlo_v2.py');
const MC_RESULTS = path.join(__dirname, 'cosmogenesis_mc_v2_results.json');
const PYTHON = process.env.PYTHON_PATH || 'python';

/**
 * Handle POST /api/cosmogenesis/montecarlo
 * Runs the Python MC engine and returns JSON results.
 */
function handleMonteCarloRequest(req, res, body) {
    // Run MC script
    execFile(PYTHON, [MC_SCRIPT], {
        cwd: __dirname,
        timeout: 120000,  // 2 min max
        maxBuffer: 10 * 1024 * 1024,  // 10 MB
    }, (error, stdout, stderr) => {
        if (error) {
            res.writeHead(500, {'Content-Type': 'application/json'});
            res.end(JSON.stringify({
                error: 'MC execution failed',
                detail: error.message,
            }));
            return;
        }

        // Read the results JSON
        const fs = require('fs');
        try {
            const results = JSON.parse(fs.readFileSync(MC_RESULTS, 'utf8'));
            results.endpoint = '/api/cosmogenesis/montecarlo';
            results.api_version = '1.0.0';
            res.writeHead(200, {'Content-Type': 'application/json'});
            res.end(JSON.stringify(results, null, 2));
        } catch (readErr) {
            res.writeHead(500, {'Content-Type': 'application/json'});
            res.end(JSON.stringify({
                error: 'Failed to read MC results',
                detail: readErr.message,
            }));
        }
    });
}

/**
 * Handle GET /api/cosmogenesis/health
 */
function handleCosmogenesisHealth(req, res) {
    const fs = require('fs');
    const hasScript = fs.existsSync(MC_SCRIPT);
    res.writeHead(200, {'Content-Type': 'application/json'});
    res.end(JSON.stringify({
        status: 'ok',
        service: 'uqff_server_mc_route',
        mc_script_available: hasScript,
        results_file: MC_RESULTS,
    }));
}

module.exports = { handleMonteCarloRequest, handleCosmogenesisHealth };
"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(js_code)
    return output_path


# ── §5  SERVER LAUNCHER ──────────────────────────────────────────────────

def start_server(port: int = MC_PORT, host: str = MC_HOST):
    """Start the Monte Carlo REST API server."""
    server = HTTPServer((host, port), MCRequestHandler)
    print(f"\n  Monte Carlo REST API listening on http://{host}:{port}")
    print(f"  Endpoints:")
    print(f"    POST /api/cosmogenesis/montecarlo  — run MC (body: optional config)")
    print(f"    GET  /api/cosmogenesis/status       — last run statistics")
    print(f"    GET  /api/cosmogenesis/health       — health check")
    print(f"\n  Press Ctrl+C to stop.\n")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n  Server stopped.")
        server.server_close()


# ── §6  MAIN ──────────────────────────────────────────────────────────────

def main():
    """Run MC, export results, optionally start REST server."""
    print("=" * 72)
    print("Monte Carlo REST Exporter — Cosmogenesis MC → REST API")
    print("=" * 72)
    print(f"  MC engine available: {HAS_MC_ENGINE}")
    print(f"  Port: {MC_PORT}")

    # 1. Generate Node.js route file
    js_path = generate_node_route()
    print(f"\n[OK] Node.js route generated: {js_path}")
    print(f"     Add to uqff_server.js for /api/cosmogenesis/* on port 3141")

    # 2. Run one MC batch and export
    if HAS_MC_ENGINE:
        print(f"\n── Running MC batch (n={DEFAULT_N_SAMPLES:,}) ──")
        result = _runner.run()

        if "error" in result:
            print(f"  [ERROR] {result['error']}")
        else:
            print(f"  Accepted:   {result['n_accepted']:,} / {result['config']['n_samples']:,}")
            print(f"  Rate:       {result['acceptance_rate']:.4f}")
            print(f"  Throughput: {result['samples_per_second']:,.0f} samples/s")
            print(f"  Elapsed:    {result['elapsed_seconds']:.3f} s")

            # Export
            out_path = "montecarlo_rest_export.json"
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump(result, f, indent=2, default=str)
            print(f"\n[OK] Results exported: {out_path}")

            # GW validation summary
            gw = result.get("gw_event_validation", [])
            if gw:
                print(f"\n── GW Event Validation ──")
                for ev in gw:
                    status = "PASS" if ev.get("consistent", False) else "FAIL"
                    print(f"  {ev['event']:10s}:  resid={ev.get('residual_fraction', 0):.4f}  [{status}]")

        # Status check
        status = _runner.get_status()
        print(f"\n── API Status ──")
        print(f"  {json.dumps(status, indent=2, default=str)}")
    else:
        print("\n  [WARN] cosmogenesis_montecarlo_v2.py not importable — skipping MC run")

    # 3. Offer to start server
    print(f"\n── Server Mode ──")
    print(f"  To start the REST API server:")
    print(f"    python montecarlo_rest_exporter.py --serve")
    print(f"  Then: curl http://127.0.0.1:{MC_PORT}/api/cosmogenesis/montecarlo -X POST")

    if "--serve" in sys.argv:
        start_server()

    print(f"\n{'=' * 72}")
    print(f"EXPORT COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
