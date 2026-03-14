"""
uqff_ipc.py — Python-side UQFF IPC server
==========================================
Serves the StarMagic_UQFF named pipe (Windows) or a TCP socket (cross-platform),
routing incoming UQFF IPC messages (defined in uqff_ipc.h) to the appropriate
Python calculator:

  MessageType.COMPUTE_SINGLE / GET_PARAMETER → QCalc.py::UnifiedFieldSolver
  MessageType.CP2_TRIGGER                    → CondensedPhysics2.py classes
  MessageType.WOLFRAM_EVAL                   → Deferred (WSTP not available in Python)
  MessageType.GROK_QUERY                     → AI SDK (optional, if configured)

Calibrated constants (GrokThread7b0e961f, Sept 2025):
  kappa   = 0.0005 / day
  [SSq]   = 0.57
  H_SCm   ≈ 0.99
  U_UA    = 0.0001
  F_rel   = 4.31e33 N  (2024 LEP reanalysis)

Usage
-----
  # Start server (blocks, handle connections in threads)
  python uqff_ipc.py --serve [--port 41000]

  # One-shot compute (for subprocess mode from C++)
  python uqff_ipc.py --json-b64 <base64-encoded-json>

Author: Daniel T. Murphy
Date  : March 2026 (v3.2 — Session 60 sync: CP3 112 classes, Aggregator v2.4.0, 243 papers)
"""

from __future__ import annotations

import argparse
import base64
import json
import logging
import os
import socket
import struct
import sys
import time
import threading
from typing import Any, Dict, Optional

# ──────────────────────────────────────────────────────────────
# PROTOCOL CONSTANTS  (mirror of uqff_ipc.h)
# ──────────────────────────────────────────────────────────────
HEADER_MAGIC   = 0x55514646   # 'UQFF'
PROTOCOL_VER   = 0x0301
HEADER_FORMAT  = "<IHHIIII"   # little-endian: magic, version, type, session, seq, payload_len, crc32
HEADER_SIZE    = struct.calcsize(HEADER_FORMAT)  # should be 22

# Message type values
class MsgType:
    COMPUTE_SINGLE   = 0x0101
    COMPUTE_BATCH    = 0x0102
    QUERY_PHYSICS_TERM = 0x0103
    SEND_DATASET     = 0x0202
    SESSION_START    = 0x0301
    SESSION_END      = 0x0302
    PING             = 0x0303
    CANCEL           = 0x0304
    SET_PARAMETER    = 0x0305
    GET_PARAMETER    = 0x0306
    RESPONSE_DATA    = 0x0401
    RESPONSE_JSON    = 0x0402
    RESPONSE_EQUATION= 0x0403
    RESPONSE_EQUATION_SET = 0x0404
    RESPONSE_PONG    = 0x0406
    ERROR_GENERAL    = 0x0501
    ERROR_UNKNOWN_SYSTEM  = 0x0502
    ERROR_TIMEOUT    = 0x0504
    STATUS_BUSY      = 0x0507
    STATUS_READY     = 0x0508
    CP2_TRIGGER      = 0x0A00
    CP2_RESPONSE     = 0x0A01
    GROK_QUERY       = 0x0A04
    GROK_RESPONSE    = 0x0A05

# ──────────────────────────────────────────────────────────────
# CALIBRATED DEFAULTS  (GrokThread7b0e)
# ──────────────────────────────────────────────────────────────
DEFAULTS: Dict[str, Any] = {
    "kappa":   0.0005,     # /day
    "SSq":     0.57,
    "H_SCm":   0.99,
    "U_UA":    0.0001,     # Gaia DR4 calibration
    "F_rel":   4.31e33,    # N — 2024 LEP reanalysis
    "t_n":    -2512.0,
}

# ──────────────────────────────────────────────────────────────
# CRC-32  (pure Python, matches uqff_ipc.h compute_crc32)
# ──────────────────────────────────────────────────────────────
def _crc32(data: bytes) -> int:
    crc = 0xFFFFFFFF
    for byte in data:
        crc ^= byte
        for _ in range(8):
            crc = (crc >> 1) ^ (0xEDB88320 & -(crc & 1))
    return crc ^ 0xFFFFFFFF

# ──────────────────────────────────────────────────────────────
# FRAMING HELPERS
# ──────────────────────────────────────────────────────────────
def build_frame(msg_type: int, session_id: int, sequence: int, payload: bytes) -> bytes:
    crc = _crc32(payload) if payload else 0
    header = struct.pack(
        HEADER_FORMAT,
        HEADER_MAGIC, PROTOCOL_VER, msg_type,
        session_id, sequence, len(payload), crc,
    )
    return header + payload


def parse_header(data: bytes) -> Optional[Dict[str, Any]]:
    if len(data) < HEADER_SIZE:
        return None
    magic, ver, msg_type, sess, seq, plen, crc = struct.unpack(HEADER_FORMAT, data[:HEADER_SIZE])
    if magic != HEADER_MAGIC or ver != PROTOCOL_VER:
        return None
    return {"type": msg_type, "session": sess, "seq": seq,
            "payload_len": plen, "crc": crc}


def validate_payload(payload: bytes, expected_crc: int) -> bool:
    return _crc32(payload) == expected_crc


# ──────────────────────────────────────────────────────────────
# CALCULATOR ROUTING
# ──────────────────────────────────────────────────────────────
CP2_TRIGGERS = {
    "Orb10", "Orb11", "Orb12", "Orb13", "Orb14", "Orb15", "Orb16",
    "Red Mercury", "Plasmoid", "UFEQFET", "Monte Carlo", "Relativistic",
    "ResonanceGravity", "AsymCap", "FractalTime", "VacuumProb",
    "26Layer", "CompressedGravity", "BuoyancyProof",
    "ShapiroWilk", "qwaveNormal", "rotorCS", "H2OH2PES",
    "DPMfreqMUGE", "BECalpha", "complexUi", "SuperCondUI",
    # Thread 7b0e additions
    "AtomicUQFF", "BoseNuclear", "NeutrinoSED", "Z_DPM", "HVP_tau",
}


def _route(dataset: dict) -> str:
    trigger = dataset.get("trigger", "") or dataset.get("system_name", "")
    for kw in CP2_TRIGGERS:
        if kw.lower() in trigger.lower():
            return "cp2"
    return "qcalc"


def _compute_qcalc(dataset: dict) -> dict:
    """Import QCalc.py and call UnifiedFieldSolver."""
    try:
        import importlib.util, pathlib
        spec = importlib.util.spec_from_file_location(
            "QCalc", pathlib.Path(__file__).parent / "QCalc.py")
        qc = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(qc)                     # type: ignore[union-attr]
        solver = qc.UnifiedFieldSolver()
        return solver.compute(dataset)
    except Exception as exc:
        return {"error": str(exc), "backend": "QCalc"}


def _compute_cp2(dataset: dict) -> dict:
    """Import CondensedPhysics2.py and route to the right class."""
    try:
        import importlib.util, pathlib
        spec = importlib.util.spec_from_file_location(
            "CondensedPhysics2",
            pathlib.Path(__file__).parent / "CondensedPhysics2.py")
        cp2 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cp2)                    # type: ignore[union-attr]

        # Try ALL_CALCULATORS dict first (CondensedPhysicsAggregator pattern)
        if hasattr(cp2, "ALL_CALCULATORS"):
            trigger = dataset.get("trigger", "")
            for name, cls in cp2.ALL_CALCULATORS.items():
                if trigger and trigger.lower() in name.lower():
                    calculator = cls()
                    return calculator.compute(dataset)

        # Fallback: first class in module that has a .compute() method
        for attr in dir(cp2):
            obj = getattr(cp2, attr)
            if isinstance(obj, type) and hasattr(obj, "compute") and attr != "object":
                try:
                    return obj().compute(dataset)
                except Exception:
                    continue

        return {"error": "No suitable CP2 calculator found", "backend": "CP2"}
    except Exception as exc:
        return {"error": str(exc), "backend": "CondensedPhysics2"}


def handle_message(msg_type: int, payload_bytes: bytes, session_id: int) -> bytes:
    """
    Process one decoded IPC message and return the response frame bytes.
    """
    seq = 0  # server-to-client sequence (simplified)

    if msg_type == MsgType.PING:
        body = json.dumps({"status": "pong", "version": "3.1"}).encode()
        return build_frame(MsgType.RESPONSE_PONG, session_id, seq, body)

    if msg_type in (MsgType.COMPUTE_SINGLE, MsgType.CP2_TRIGGER,
                    MsgType.SEND_DATASET, MsgType.COMPUTE_BATCH):
        try:
            dataset = json.loads(payload_bytes.decode("utf-8")) if payload_bytes else {}
        except json.JSONDecodeError as e:
            err = json.dumps({"error": f"JSON decode: {e}"}).encode()
            return build_frame(MsgType.ERROR_GENERAL, session_id, seq, err)

        # Apply calibrated defaults where not provided
        for key, val in DEFAULTS.items():
            dataset.setdefault(key, val)

        t0 = time.perf_counter()
        route = "cp2" if msg_type == MsgType.CP2_TRIGGER else _route(dataset)
        result = _compute_cp2(dataset) if route == "cp2" else _compute_qcalc(dataset)
        elapsed = (time.perf_counter() - t0) * 1000.0

        result["elapsed_ms"] = elapsed
        result["backend"] = route
        body = json.dumps(result).encode("utf-8")
        resp_type = MsgType.CP2_RESPONSE if route == "cp2" else MsgType.RESPONSE_JSON
        return build_frame(resp_type, session_id, seq, body)

    if msg_type == MsgType.SESSION_START:
        body = json.dumps({"session_id": session_id, "status": "ok"}).encode()
        return build_frame(MsgType.STATUS_READY, session_id, seq, body)

    if msg_type == MsgType.SESSION_END:
        body = json.dumps({"status": "closed"}).encode()
        return build_frame(MsgType.STATUS_READY, session_id, seq, body)

    err = json.dumps({"error": f"Unhandled message type: 0x{msg_type:04X}"}).encode()
    return build_frame(MsgType.ERROR_GENERAL, session_id, seq, err)


# ──────────────────────────────────────────────────────────────
# TCP SERVER LOOP
# ──────────────────────────────────────────────────────────────
def _handle_connection(conn: socket.socket, addr: tuple) -> None:
    log = logging.getLogger("uqff_ipc")
    log.info("Connection from %s", addr)
    try:
        while True:
            # Read header
            raw = b""
            while len(raw) < HEADER_SIZE:
                chunk = conn.recv(HEADER_SIZE - len(raw))
                if not chunk:
                    return
                raw += chunk

            hdr = parse_header(raw)
            if hdr is None:
                log.warning("Invalid header from %s", addr)
                return

            # Read payload
            payload = b""
            remaining = hdr["payload_len"]
            while remaining > 0:
                chunk = conn.recv(remaining)
                if not chunk:
                    return
                payload += chunk
                remaining -= len(chunk)

            if hdr["payload_len"] > 0 and not validate_payload(payload, hdr["crc"]):
                log.warning("CRC mismatch from %s", addr)
                err = json.dumps({"error": "CRC mismatch"}).encode()
                conn.sendall(build_frame(MsgType.ERROR_GENERAL, hdr["session"], 0, err))
                continue

            response = handle_message(hdr["type"], payload, hdr["session"])
            conn.sendall(response)

    except (ConnectionResetError, BrokenPipeError):
        pass
    except Exception as exc:
        log.exception("Error handling connection: %s", exc)
    finally:
        conn.close()


def serve(host: str = "127.0.0.1", port: int = 41000) -> None:
    logging.basicConfig(
        format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
        level=logging.INFO,
    )
    log = logging.getLogger("uqff_ipc")
    log.info("UQFF IPC server v3.1 starting on %s:%d", host, port)
    log.info("Calibrated constants: F_rel=4.31e33, U_UA=0.0001, kappa=0.0005, SSq=0.57")

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as srv:
        srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        srv.bind((host, port))
        srv.listen(16)
        log.info("Listening …")
        while True:
            try:
                conn, addr = srv.accept()
                t = threading.Thread(target=_handle_connection, args=(conn, addr), daemon=True)
                t.start()
            except KeyboardInterrupt:
                log.info("Shutting down.")
                break


# ──────────────────────────────────────────────────────────────
# ONE-SHOT JSON MODE  (subprocess from C++ PythonBackend)
# ──────────────────────────────────────────────────────────────
def run_oneshot(json_b64: str) -> None:
    try:
        dataset = json.loads(base64.b64decode(json_b64).decode("utf-8"))
    except Exception as exc:
        print(json.dumps({"error": f"Decode failed: {exc}"}))
        sys.exit(1)

    for key, val in DEFAULTS.items():
        dataset.setdefault(key, val)

    route = _route(dataset)
    result = _compute_cp2(dataset) if route == "cp2" else _compute_qcalc(dataset)
    print(json.dumps(result))


# ──────────────────────────────────────────────────────────────
# ENTRY POINT
# ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="UQFF IPC server / one-shot runner")
    parser.add_argument("--serve",      action="store_true",  help="Start TCP IPC server")
    parser.add_argument("--host",       default="127.0.0.1",  help="Bind host (default 127.0.0.1)")
    parser.add_argument("--port",       type=int, default=41000, help="Bind port (default 41000)")
    parser.add_argument("--json-b64",   metavar="B64",        help="Base-64 encoded JSON dataset for one-shot compute")
    args = parser.parse_args()

    if args.serve:
        serve(args.host, args.port)
    elif args.json_b64:
        run_oneshot(args.json_b64)
    else:
        parser.print_help()
        sys.exit(1)
