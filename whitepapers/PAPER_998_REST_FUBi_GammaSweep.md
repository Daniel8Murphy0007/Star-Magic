---
paper_id: PAPER_998
title: "REST F_U_Bi Gamma Sweep Endpoints — /api/fubi/inside-outside + /api/fubi/gamma-sweep"
session: 218
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [REST, API, endpoint, F_U_Bi, Gamma, sweep, inside-outside]
crosslinks: [PAPER_988, PAPER_989, PAPER_995]
calibration: {port: 3141, routes: 18, new_routes: 2}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_998: REST F_U_Bi Gamma Sweep Endpoints

## Abstract

We add two new REST API endpoints to `uqff_server.js` (Port 3141), bringing the total to 18 routes:

## 1. POST /api/fubi/inside-outside

Computes the F_U_Bi inside-to-outside buoyancy mass portion.

### Request

```json
{"M_kg": 1.989e30, "r_m": 1.496e11}
```

### Response

```json
{
    "F_U_Bi": 2.33e40, "ratio": 0.606,
    "Ug": 5.93e-03, "Ub": -3.57e-03,
    "rho_SCm": 1.0e-10, "V_region": 1e48, "S26": 19.56
}
```

## 2. POST /api/fubi/gamma-sweep

Computes the aggregate F_U_Bi_i across systems at multiple Γ values.

### Request

```json
{"gamma_THz_list": [0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0]}
```

### Response

```json
{
    "systems": 20, "gamma_count": 7,
    "results": [{"gamma_THz": 0.01, "aggregate": -6.11e13}, ...]
}
```

## 3. Implementation

File: `uqff_server.js`, routes 17–18. CP4 class #582.
