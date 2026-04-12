---
paper_id: PAPER_988
title: "REST F_U_Bi_i Endpoint — POST /api/fubi/master"
session: 217
date: 2026-04-12
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [REST, API, endpoint, server, HTTP, JSON, UQFF]
crosslinks: [PAPER_979, PAPER_985, PAPER_987]
calibration: {port: 3141, route: "/api/fubi/master", method: "POST"}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_988: REST F_U_Bi_i Endpoint — POST /api/fubi/master

## Abstract

We document the 16th REST API endpoint added to `uqff_server.js` (Port 3141): `POST /api/fubi/master`. This endpoint accepts physical parameters $(M, r, t, \Gamma)$ and returns the complete 6-layer $F_{U,\text{Bi}_i}$ decomposition as JSON. It serves as the HTTP interface for all external consumers of the master buoyancy force.

## 1. Endpoint Specification

| Property | Value |
|----------|-------|
| Method | POST |
| Path | `/api/fubi/master` |
| Port | 3141 |
| Content-Type | application/json |

### Request Body

```json
{
    "M_kg": 1.989e30,
    "r_m": 1.496e11,
    "t_s": 86400,
    "Gamma_THz": 0.1
}
```

### Response

```json
{
    "F_U_Bi_i": -2.405685e-02,
    "layers": {
        "Ug_26layer": 5.930e-03,
        "Ub_26layer": -3.573e-03,
        "Um_magnetism": 1.234e-10,
        "UA_aether": -5.678e-12,
        "Fn_neutron": 1.0e-09,
        "Phi_phonon": 19.56,
        "E_net": -1.234,
        "S26": 19.56
    },
    "metadata": {
        "SSq": 0.57,
        "beta_i": 0.603,
        "kappa_per_day": 0.0005,
        "omega_SCm_THz": 1.25,
        "layer_count": 6,
        "system_count": 1
    }
}
```

## 2. Error Handling

| Status | Condition |
|--------|-----------|
| 200 | Successful computation |
| 400 | Missing required fields or non-numeric values |
| 500 | Internal computation error (NaN/Inf guard) |

## 3. Integration with uqff_server.js

The endpoint is the 16th route in the REST API server, complementing existing routes:
- `/api/gravity/compressed` (POST)
- `/api/gravity/resonance` (POST)
- `/api/systems` (GET)
- `/api/fubi/master` (POST) — **NEW**
- ... (12 others)

## 4. Usage Example

```bash
curl -X POST http://localhost:3141/api/fubi/master \
  -H "Content-Type: application/json" \
  -d '{"M_kg": 1.989e30, "r_m": 1.496e11, "t_s": 86400, "Gamma_THz": 0.1}'
```

## 5. Implementation

Route handler added to `uqff_server.js` — computes all 6 layers server-side using the same constants and equations as `fubi_master_calculator.py`.

## References
- PAPER_979: Complete 6-Layer F_U_Bi_i
- PAPER_985: Production Kernel

---

## §A. Cosmogenesis-Linked Lagrangian

The REST endpoint makes the Lagrangian-derived $F_{U,\text{Bi}_i}$ accessible to any HTTP client, enabling distributed computation across the 6-tier architecture.

## §B. VDS/DVP/BSH Deep Synthesis

- **VDS:** The response `metadata.kappa_per_day` exposes the vacuum density growth rate to API consumers.
- **DVP:** Layer decomposition in the response reveals the DPM contribution via `Um_magnetism`.
- **BSH:** The `S26` field in the response is the buoyancy harmonic sum, central to all force calculations.
