/**
 * @file uqff_server.js
 * @brief HTTP REST API server for UQFF Star-Magic JavaScript engine
 * 
 * Connects the 23,790-line index.js computational engine to external systems.
 * 
 * Usage:
 *   node uqff_server.js
 *   
 *   # Then from curl or any HTTP client:
 *   curl http://localhost:3141/api/compute -X POST -H "Content-Type: application/json" \
 *     -d '{"system":"SgrA*", "M":8.155e36, "r":4.4e19}'
 * 
 * Endpoints:
 *   GET  /api/health          - Health check
 *   GET  /api/constants       - Get all physics constants
 *   GET  /api/systems         - List available systems
 *   POST /api/compute         - Compute UQFF for a system
 *   POST /api/batch           - Batch computation
 *   GET  /api/modules         - List available UQFF modules
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 */

const http = require('http');
const url = require('url');
const fs = require('fs');
const path = require('path');

// Import the main UQFF engine
let uqffEngine;
try {
    uqffEngine = require('./index.js');
    console.log('✅ UQFF Engine loaded successfully');
} catch (err) {
    console.error('❌ Failed to load index.js:', err.message);
    process.exit(1);
}

// Server configuration
const PORT = process.env.UQFF_PORT || 3141;  // π × 1000 ≈ 3141
const HOST = process.env.UQFF_HOST || '127.0.0.1';

// CORS headers for cross-origin requests
const CORS_HEADERS = {
    'Access-Control-Allow-Origin': '*',
    'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
    'Access-Control-Allow-Headers': 'Content-Type, Authorization',
    'Content-Type': 'application/json'
};

/**
 * Physics constants from index.js
 */
const CONSTANTS = {
    SOLAR_MASS: 1.989e30,
    SOLAR_RADIUS: 6.96e8,
    SPEED_OF_LIGHT: 2.998e8,
    GRAVITATIONAL_CONSTANT: 6.674e-11,
    PLANCK_CONSTANT: 1.055e-34,
    RHO_VAC_UA: 7.09e-36,
    RHO_VAC_SCM: 7.09e-37,
    B_CRIT_MAGNETAR: 4.4e13,
    KAPPA: 0.0005,
    SSQ: 0.57,
    BETA_I: 0.603,
    H_SCM: 0.99,
    U_UA: 0.0001
};

/**
 * Pre-defined astrophysical systems
 */
const SYSTEMS = {
    'SgrA*': {
        name: 'Sagittarius A*',
        type: 'SMBH',
        M: 8.155e36,
        r: 4.4e19,
        B0: 1e-4,
        description: 'Milky Way Central Black Hole'
    },
    'SGR1745': {
        name: 'SGR1745 Magnetar',
        type: 'Magnetar',
        M: 2.8e30,
        r: 10000,
        B0: 1e14,
        description: 'Magnetar near Galactic Center'
    },
    'M87*': {
        name: 'M87* Black Hole',
        type: 'SMBH',
        M: 6.5e39,
        r: 1.55e20,
        B0: 1e-3,
        description: 'Virgo A Supermassive Black Hole'
    },
    'Sun': {
        name: 'Sun',
        type: 'Star',
        M: 1.989e30,
        r: 6.96e8,
        B0: 1e-4,
        description: 'G-type Main Sequence Star'
    },
    'NGC3596': {
        name: 'NGC 3596',
        type: 'Galaxy',
        M: 1e41,
        r: 3.086e22,
        B0: 1e-9,
        description: 'Spiral Galaxy'
    }
};

/**
 * Compute Ug1: Magnetic dipole gravity component
 */
function computeUg1(params) {
    const G = CONSTANTS.GRAVITATIONAL_CONSTANT;
    const mu_0 = 1.25663706212e-6;
    const { M, r, B0 = 1e-4 } = params;
    
    if (r <= 0 || M <= 0) return 0;
    
    const mu = B0 * r * r * r;
    return (mu_0 * mu * mu) / (4 * Math.PI * Math.pow(r, 5));
}

/**
 * Compute Ug2: Charge/reactivity gravity component
 */
function computeUg2(params) {
    const G = CONSTANTS.GRAVITATIONAL_CONSTANT;
    const { M, r } = params;
    
    if (r <= 0 || M <= 0) return 0;
    
    const q_ratio = 1e-10;
    return G * M * q_ratio / (r * r);
}

/**
 * Compute Ug3: String rotation gravity component
 */
function computeUg3(params) {
    const { M, r, t = 0 } = params;
    const G = CONSTANTS.GRAVITATIONAL_CONSTANT;
    
    if (r <= 0 || M <= 0) return 0;
    
    const omega = 2 * Math.PI / (24 * 3600);
    return G * M / (r * r) * Math.cos(omega * t) * 1e-6;
}

/**
 * Compute Ug4: Vacuum concentration gravity component
 */
function computeUg4(params) {
    const { r } = params;
    const c = CONSTANTS.SPEED_OF_LIGHT;
    
    if (r <= 0) return 0;
    
    const rho_ratio = CONSTANTS.RHO_VAC_UA / CONSTANTS.RHO_VAC_SCM;
    return (rho_ratio - 1) * c * c / r * 1e-20;
}

/**
 * Compute F_U_Bi_i: Master buoyancy-gravity force
 */
function computeFUBi(params) {
    const { M, t = 0 } = params;
    
    const Ug1 = computeUg1(params);
    const Ug2 = computeUg2(params);
    const Ug3 = computeUg3(params);
    const Ug4 = computeUg4(params);
    
    const Ug_sum = Ug1 + Ug2 + Ug3 + Ug4;
    const beta_term = CONSTANTS.BETA_I * Ug_sum;
    
    const k_LENR = 1e-10;
    const omega_LENR = 1.2e12;
    const lenr_factor = 1.0 + k_LENR * Math.sin(omega_LENR * t);
    
    const scm_factor = CONSTANTS.SSQ * (1.0 - Math.exp(-CONSTANTS.KAPPA * t));
    
    return M * beta_term * lenr_factor * (1.0 + scm_factor);
}

/**
 * Compute compressed gravity (26-layer)
 */
function computeCompressedG(params) {
    const G = CONSTANTS.GRAVITATIONAL_CONSTANT;
    const { M, r } = params;
    
    if (r <= 0 || M <= 0) return 0;
    
    const g_N = G * M / (r * r);
    
    const Ug1 = computeUg1(params);
    const Ug2 = computeUg2(params);
    const Ug3 = computeUg3(params);
    const Ug4 = computeUg4(params);
    
    let layer_sum = 0;
    for (let i = 1; i <= 26; i++) {
        const Q_i = 1.0;
        const UA_i = 1.0 / (1.0 + i * 0.01);
        const SCm_i = 0.99;
        layer_sum += (Ug1 + Ug2 + Ug3 + Ug4) * Q_i * UA_i * SCm_i / 26.0;
    }
    
    return g_N + layer_sum;
}

/**
 * Full UQFF computation
 */
function computeFull(params) {
    const Ug1 = computeUg1(params);
    const Ug2 = computeUg2(params);
    const Ug3 = computeUg3(params);
    const Ug4 = computeUg4(params);
    const Ug_sum = Ug1 + Ug2 + Ug3 + Ug4;
    const F_U_Bi_i = computeFUBi(params);
    const g_compressed = computeCompressedG(params);
    
    return {
        system_name: params.name || 'Custom',
        F_U_Bi_i,
        g_compressed,
        Ug1,
        Ug2,
        Ug3,
        Ug4,
        Ug_sum,
        F_U: F_U_Bi_i,
        status: 'success',
        timestamp: new Date().toISOString(),
        computation_time_ms: 0,
        framework_version: '3.0.0'
    };
}

/**
 * Request handler
 */
function handleRequest(req, res) {
    const parsedUrl = url.parse(req.url, true);
    const pathname = parsedUrl.pathname;
    const method = req.method;
    
    // Handle CORS preflight
    if (method === 'OPTIONS') {
        res.writeHead(204, CORS_HEADERS);
        res.end();
        return;
    }
    
    // Route handlers
    if (pathname === '/api/health' && method === 'GET') {
        res.writeHead(200, CORS_HEADERS);
        res.end(JSON.stringify({
            status: 'healthy',
            version: '3.0.0',
            engine: 'UQFF Star-Magic',
            uptime: process.uptime(),
            memory: process.memoryUsage()
        }));
        return;
    }
    
    if (pathname === '/api/constants' && method === 'GET') {
        res.writeHead(200, CORS_HEADERS);
        res.end(JSON.stringify(CONSTANTS));
        return;
    }
    
    if (pathname === '/api/systems' && method === 'GET') {
        res.writeHead(200, CORS_HEADERS);
        res.end(JSON.stringify(SYSTEMS));
        return;
    }
    
    if (pathname === '/api/modules' && method === 'GET') {
        // List available modules from index.js exports
        const modules = Object.keys(uqffEngine || {}).filter(k => k.includes('Module'));
        res.writeHead(200, CORS_HEADERS);
        res.end(JSON.stringify({
            count: modules.length,
            modules: modules
        }));
        return;
    }
    
    if (pathname === '/api/compute' && method === 'POST') {
        let body = '';
        req.on('data', chunk => body += chunk);
        req.on('end', () => {
            try {
                const params = JSON.parse(body);
                
                // If system name provided, merge with predefined params
                if (params.system && SYSTEMS[params.system]) {
                    Object.assign(params, SYSTEMS[params.system], params);
                }
                
                const startTime = Date.now();
                const result = computeFull(params);
                result.computation_time_ms = Date.now() - startTime;
                
                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify(result));
            } catch (err) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({
                    status: 'error',
                    error: err.message
                }));
            }
        });
        return;
    }
    
    if (pathname === '/api/batch' && method === 'POST') {
        let body = '';
        req.on('data', chunk => body += chunk);
        req.on('end', () => {
            try {
                const { systems } = JSON.parse(body);
                
                const startTime = Date.now();
                const results = systems.map(params => {
                    if (params.system && SYSTEMS[params.system]) {
                        Object.assign(params, SYSTEMS[params.system], params);
                    }
                    return computeFull(params);
                });
                
                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify({
                    count: results.length,
                    total_time_ms: Date.now() - startTime,
                    results
                }));
            } catch (err) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({
                    status: 'error',
                    error: err.message
                }));
            }
        });
        return;
    }
    
    // POST /api/phonon/jet — Phonon-modulated jet power computation
    if (pathname === '/api/phonon/jet' && method === 'POST') {
        let body = '';
        req.on('data', chunk => body += chunk);
        req.on('end', () => {
            try {
                const params = JSON.parse(body);
                const M_bh = params.M_bh || 6.5e9 * CONSTANTS.SOLAR_MASS;
                const a_spin = params.a_spin || 0.9;
                const B_field = params.B_field || 50;
                const Gamma_THz = params.Gamma_THz || 0.1;
                
                const Gamma_rad = Gamma_THz * 2 * Math.PI * 1e12;
                const r_S = 2 * CONSTANTS.GRAVITATIONAL_CONSTANT * M_bh / (CONSTANTS.SPEED_OF_LIGHT ** 2);
                const r_H = r_S / 2 * (1 + Math.sqrt(Math.max(1 - a_spin ** 2, 0)));
                const P_BZ = (B_field ** 2 / (8 * Math.PI)) * (r_H / CONSTANTS.SPEED_OF_LIGHT) ** 2 * a_spin ** 2 * CONSTANTS.SPEED_OF_LIGHT;
                
                const OMEGA_SCM = 2 * Math.PI * 1.25e12;
                const GAMMA_SCM = 2 * Math.PI * 0.1e12;
                const A_JET = 1.5;
                const SIGMA_G = 0.08 * 2 * Math.PI * 1e12;
                
                const delta = Gamma_rad - GAMMA_SCM;
                const M_jet = 1 + A_JET * Math.exp(-delta ** 2 / (2 * SIGMA_G ** 2));
                const P_jet = P_BZ * (1 + M_jet);
                const enhancement = P_jet / P_BZ;
                
                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify({
                    M_bh, a_spin, B_field, Gamma_THz,
                    M_jet,
                    P_BZ,
                    P_jet,
                    enhancement,
                    equation: `P_jet = P_BZ * (1 + M_jet) = ${P_BZ.toExponential(6)} * (1 + ${M_jet.toFixed(4)}) = ${P_jet.toExponential(6)} W`
                }));
            } catch (err) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({
                    status: 'error',
                    error: err.message
                }));
            }
        });
        return;
    }

    // ── POST /api/phonon/jet/cena — Centaurus A jet power (Session 213) ──
    if (method === 'POST' && pathname === '/api/phonon/jet/cena') {
        let body = '';
        req.on('data', c => body += c);
        req.on('end', () => {
            try {
                const d = JSON.parse(body);
                const M = (d.M_Msun || 5.5e7) * 1.989e30;
                const a = d.a_spin || 0.70;
                const B = d.B_field || 3000;
                const gTHz = d.Gamma_THz || 0.10;
                const Ajet = d.A_jet || 0.95;
                const rS = 2 * 6.674e-11 * M / (3e8)**2;
                const rH = rS / 2 * (1 + Math.sqrt(Math.max(1 - a*a, 0)));
                const PBZ = (B*B / (8 * Math.PI)) * (rH / 3e8)**2 * a*a * 3e8;
                const G0 = 2 * Math.PI * 0.1e12;
                const sG = 0.08 * 2 * Math.PI * 1e12;
                const Grad = 2 * Math.PI * gTHz * 1e12;
                const delta = Grad - G0;
                const Mj = 1 + Ajet * Math.exp(-delta*delta / (2 * sG*sG));
                const Pjet = PBZ * (1 + Mj);
                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify({
                    system: 'Centaurus_A', P_BZ: PBZ, P_jet: Pjet,
                    enhancement: Pjet / PBZ, Gamma_THz: gTHz,
                    P_BZ_erg_s: PBZ * 1e7,
                }));
            } catch (e) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({ error: e.message }));
            }
        });
        return;
    }

    // ── POST /api/phonon/jet/txs0506 — TXS 0506+056 jet power (Session 213) ──
    if (method === 'POST' && pathname === '/api/phonon/jet/txs0506') {
        let body = '';
        req.on('data', c => body += c);
        req.on('end', () => {
            try {
                const d = JSON.parse(body);
                const M = (d.M_Msun || 3e8) * 1.989e30;
                const a = d.a_spin || 0.85;
                const B = d.B_field || 8000;
                const gTHz = d.Gamma_THz || 0.10;
                const Ajet = d.A_jet || 1.20;
                const rS = 2 * 6.674e-11 * M / (3e8)**2;
                const rH = rS / 2 * (1 + Math.sqrt(Math.max(1 - a*a, 0)));
                const PBZ = (B*B / (8 * Math.PI)) * (rH / 3e8)**2 * a*a * 3e8;
                const G0 = 2 * Math.PI * 0.1e12;
                const sG = 0.08 * 2 * Math.PI * 1e12;
                const Grad = 2 * Math.PI * gTHz * 1e12;
                const delta = Grad - G0;
                const Mj = 1 + Ajet * Math.exp(-delta*delta / (2 * sG*sG));
                const Pjet = PBZ * (1 + Mj);
                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify({
                    system: 'TXS_0506+056', P_BZ: PBZ, P_jet: Pjet,
                    enhancement: Pjet / PBZ, Gamma_THz: gTHz,
                    icecube: 'IceCube-170922A',
                    P_BZ_erg_s: PBZ * 1e7,
                }));
            } catch (e) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({ error: e.message }));
            }
        });
        return;
    }

    // POST /api/phonon/bcs — BCS gap equation solve
    if (method === 'POST' && url === '/api/phonon/bcs') {
        let body = '';
        req.on('data', chunk => { body += chunk; });
        req.on('end', () => {
            try {
                const params = body ? JSON.parse(body) : {};
                const T = parseFloat(params.T_K || 0);
                const omegaSCm = 2 * Math.PI * 1.25e12;
                const hbar = 1.055e-34;
                const kB = 1.381e-23;
                const ssq = 0.57;
                let S26 = 0;
                for (let k = 1; k <= 26; k++) S26 += Math.exp(-ssq * k / 26);
                const fubi = 0.6;
                const deltaPF = hbar * omegaSCm / 2;
                let delta = deltaPF * S26 * fubi;
                if (T > 0) {
                    for (let i = 0; i < 500; i++) {
                        const arg = Math.min(delta / (2 * kB * T), 500);
                        const newD = deltaPF * Math.tanh(arg) * S26 * fubi;
                        if (Math.abs(newD - delta) < 1e-30) { delta = newD; break; }
                        delta = newD;
                    }
                }
                const Tc = (1.13 * hbar * omegaSCm / kB) * Math.exp(-1 / 0.27);
                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify({
                    status: 'ok', T_K: T,
                    Delta_J: delta, Delta_eV: delta / 1.602e-19,
                    T_c_K: Tc, session: 214
                }));
            } catch (e) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({ status: 'error', error: e.message }));
            }
        });
        return;
    }

    // POST /api/phonon/spectral-ladder — 26-state HRes ladder
    if (method === 'POST' && url === '/api/phonon/spectral-ladder') {
        let body = '';
        req.on('data', chunk => { body += chunk; });
        req.on('end', () => {
            try {
                const params = body ? JSON.parse(body) : {};
                const gammaTHz = parseFloat(params.Gamma_THz || 0.10);
                const gamma = 2 * Math.PI * gammaTHz * 1e12;
                const hbar = 1.055e-34;
                const omegaSCm = 2 * Math.PI * 1.25e12;
                const ssq = 0.57;
                let S26 = 0;
                for (let k = 1; k <= 26; k++) S26 += Math.exp(-ssq * k / 26);
                const E0 = hbar * omegaSCm;
                const levels = [];
                for (let n = 1; n <= 26; n++) {
                    const En = E0 * Math.pow(2 * Math.PI, n / 3) * S26;
                    const omegaN = En / hbar;
                    const Qn = omegaN / (2 * gamma);
                    levels.push({ n, E_eV: En / 1.602e-19, Q: Qn,
                        regime: Qn > 10 ? 'narrow' : (Qn > 3 ? 'optimal' : 'broad') });
                }
                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify({
                    status: 'ok', Gamma_THz: gammaTHz,
                    levels, session: 214
                }));
            } catch (e) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({ status: 'error', error: e.message }));
            }
        });
        return;
    }

    // POST /api/phonon/ns — NS phonon effects for GW190425
    if (method === 'POST' && url === '/api/phonon/ns') {
        let body = '';
        req.on('data', c => body += c);
        req.on('end', () => {
            try {
                const p = JSON.parse(body || '{}');
                const gammaTHz = parseFloat(p.Gamma_THz || 0.10);
                const hGR = parseFloat(p.h_GR || 1e-21);
                const LambdaGR = parseFloat(p.Lambda_GR || 400);

                const OMEGA_SCM = 2 * Math.PI * 1.25e12;
                const S26 = (() => { let s = 0; for (let k = 1; k <= 26; k++) s += Math.exp(-0.57 * k / 26); return s; })();
                const gamma = 2 * Math.PI * gammaTHz * 1e12;
                const phi = Math.exp(0) * S26;  // on-resonance
                const F_UBI = 0.6;
                const hUQFF = hGR * 0.5297;
                const lambdaRatio = 1 - F_UBI * phi;
                const LambdaUQFF = LambdaGR * (1 + F_UBI * phi * 0.1);
                const P_NS = 0.5 - 0.01 * phi;

                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify({
                    status: 'ok', Gamma_THz: gammaTHz,
                    h_UQFF: hUQFF, strain_reduction_pct: 47.03,
                    lambda_ratio: lambdaRatio, Lambda_UQFF: LambdaUQFF,
                    P_NS, P_BH: 1 - P_NS, session: 215
                }));
            } catch (e) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({ status: 'error', error: e.message }));
            }
        });
        return;
    }

    // POST /api/ramanujan/26d — 26D Ramanujan summation S_26(z)
    if (method === 'POST' && url === '/api/ramanujan/26d') {
        let body = '';
        req.on('data', c => body += c);
        req.on('end', () => {
            try {
                const p = JSON.parse(body || '{}');
                const z = parseFloat(p.z || 0.57);
                const N = parseInt(p.N || 50);

                let S = 0;
                for (let n = 1; n <= N; n++) {
                    // Simplified R_n for REST performance
                    let Rn = 1;
                    for (let f = 2; f <= Math.min(n, 20); f++) Rn /= f;
                    S += Math.pow(z, n) / Math.pow(n, 26) * Rn;
                }

                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify({
                    status: 'ok', z, N,
                    S_26: S, session: 215
                }));
            } catch (e) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({ status: 'error', error: e.message }));
            }
        });
        return;
    }

    // POST /api/qgp/density — QGP vacuum density + Yang-Mills mass gap (Session 216)
    if (method === 'POST' && pathname === '/api/qgp/density') {
        let body = '';
        req.on('data', c => body += c);
        req.on('end', () => {
            try {
                const p = JSON.parse(body || '{}');
                const T_K = parseFloat(p.T_K || 2e12);
                const order = parseInt(p.order || 3);

                const SSq = 0.57;
                const Tc = 1.5e12;
                const LamQCD = 0.217e9;
                const rhoSCm = 1e-10;

                // S_26^{(k)}
                let S26k = 0;
                for (let n = 1; n <= 50; n++) {
                    let Rn = 1;
                    for (let f = 2; f <= Math.min(n, 20); f++) Rn /= f;
                    S26k += Math.pow(SSq, n) / Math.pow(n, 26) * Rn;
                }

                const expF = Math.exp(Math.max(Math.min(-(Tc - T_K) / T_K, 500), -500));
                const rho_QGP = rhoSCm * S26k * expF;
                const Delta_YM = T_K < Tc ? LamQCD * (1 - T_K / Tc) * S26k : 0;

                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify({
                    status: 'ok', T_K, order,
                    rho_QGP, Delta_YM, S26k,
                    phase: T_K > Tc ? 'QGP' : 'hadron',
                    session: 216
                }));
            } catch (e) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({ status: 'error', error: e.message }));
            }
        });
        return;
    }

    // POST /api/master/99system — 99-system master equation F_U^{(99)} (Session 216)
    if (method === 'POST' && pathname === '/api/master/99system') {
        let body = '';
        req.on('data', c => body += c);
        req.on('end', () => {
            try {
                const p = JSON.parse(body || '{}');
                const nSystems = parseInt(p.n_systems || 99);
                const SSq = 0.57;
                const G0 = 6.674e-11;
                const Msun = 1.989e30;
                const betaI = 0.603;

                let S26 = 0;
                for (let k = 1; k <= 26; k++) S26 += Math.exp(-SSq * k / 26);

                let totalFU = 0;
                for (let i = 1; i <= nSystems; i++) {
                    const M = (0.1 + i * 2) * Msun;
                    const r = 1e9 * (1 + i * 0.3);
                    let Ug = 0, Ub = 0;
                    for (let j = 1; j <= 26; j++) {
                        Ug += G0 * M / (r * r) * SSq * j / 26;
                        Ub += G0 * M / (r * r) * Math.exp(-SSq * j / 26) * betaI;
                    }
                    const Um = G0 * M / (r * r) * SSq * 0.1;
                    const UA = G0 * M / (r * r) * 1e-10;
                    totalFU += Ug + Um + UA - Ub + 1e-10 * S26 * S26;
                }

                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify({
                    status: 'ok', n_systems: nSystems,
                    F_U_99: totalFU,
                    session: 216
                }));
            } catch (e) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({ status: 'error', error: e.message }));
            }
        });
        return;
    }

    // POST /api/fubi/master — Complete 6-layer F_U_Bi_i master buoyancy (Session 217)
    if (method === 'POST' && pathname === '/api/fubi/master') {
        let body = '';
        req.on('data', c => body += c);
        req.on('end', () => {
            try {
                const p = JSON.parse(body || '{}');
                const M = parseFloat(p.M_kg || 1.989e30);
                const r = Math.max(parseFloat(p.r_m || 1.496e11), 1.0);
                const t = parseFloat(p.t_s || 86400);
                const gTHz = parseFloat(p.Gamma_THz || 0.10);
                const SSq = 0.57;
                const G0 = 6.674e-11;
                const betaI = 0.603;
                const kappa = 0.0005 / 86400;
                const wSCm = 2 * Math.PI * 1.25e12;
                const gammaRad = 2 * Math.PI * gTHz * 1e12;

                let S26 = 0;
                for (let k = 1; k <= 26; k++) S26 += Math.exp(-SSq * k / 26);

                // 26-layer Ug and Ub
                let Ug = 0, Ub = 0;
                const r2 = r * r;
                for (let i = 1; i <= 26; i++) {
                    Ug += G0 * M / r2 * SSq * i / 26;
                    Ub += G0 * M / r2 * Math.exp(-SSq * i / 26) * betaI;
                }
                const Um = G0 * M / r2 * SSq * 0.1;
                const UA = G0 * M / r2 * 1e-10;
                const Fgrav = Ug + Um + UA - Ub;

                // Phonon + neutron
                const Fn = 1e-10 * S26;
                const Phi = Math.exp(-Math.pow(wSCm - wSCm, 2) / (2 * gammaRad * gammaRad)) * S26;
                const Fpn = Fn * S26 * Phi;

                // E_net modulation
                const FUbare = Math.abs(Fgrav) + Math.abs(Fpn) + 1e-300;
                const ratio = Math.abs(Ub) / FUbare;
                const Enet = (2 * ratio - 1) * Math.exp(Math.min(kappa * t, 500)) * S26;

                const FUBii = Fgrav + Fpn * Enet;

                res.writeHead(200, CORS_HEADERS);
                res.end(JSON.stringify({
                    status: 'ok',
                    F_U_Bi_i: FUBii,
                    F_gravity: Fgrav,
                    Ug: Ug, Ub: Ub, Um: Um, UA: UA,
                    F_phonon_neutron: Fpn,
                    E_net: Enet,
                    Phi_phonon: Phi,
                    S26: S26,
                    Gamma_THz: gTHz,
                    session: 217
                }));
            } catch (e) {
                res.writeHead(400, CORS_HEADERS);
                res.end(JSON.stringify({ status: 'error', error: e.message }));
            }
        });
        return;
    }

    // 404 for unknown routes
    res.writeHead(404, CORS_HEADERS);
    res.end(JSON.stringify({
        status: 'error',
        error: 'Not Found',
        available_endpoints: [
            'GET  /api/health',
            'GET  /api/constants',
            'GET  /api/systems',
            'GET  /api/modules',
            'POST /api/compute',
            'POST /api/batch',
            'POST /api/phonon/jet',
            'POST /api/phonon/jet/cena',
            'POST /api/phonon/jet/txs0506',
            'POST /api/phonon/bcs',
            'POST /api/phonon/spectral-ladder',
            'POST /api/phonon/ns',
            'POST /api/ramanujan/26d',
            'POST /api/qgp/density',
            'POST /api/master/99system',
            'POST /api/fubi/master'
        ]
    }));
}

/**
 * Start the server
 */
const server = http.createServer(handleRequest);

server.listen(PORT, HOST, () => {
    console.log('');
    console.log('═══════════════════════════════════════════════════════════════════');
    console.log('  UQFF Star-Magic REST API Server');
    console.log('═══════════════════════════════════════════════════════════════════');
    console.log(`  URL:     http://${HOST}:${PORT}`);
    console.log(`  Version: 3.0.0`);
    console.log(`  Engine:  index.js (23,790 lines, 106 systems)`);
    console.log('');
    console.log('  Endpoints:');
    console.log('    GET  /api/health     - Health check');
    console.log('    GET  /api/constants  - Physics constants');
    console.log('    GET  /api/systems    - Predefined systems');
    console.log('    GET  /api/modules    - Available UQFF modules');
    console.log('    POST /api/compute    - Compute UQFF for system');
    console.log('    POST /api/batch      - Batch computation');
    console.log('═══════════════════════════════════════════════════════════════════');
    console.log('');
});

// Graceful shutdown
process.on('SIGINT', () => {
    console.log('\nShutting down UQFF server...');
    server.close(() => {
        console.log('Server closed.');
        process.exit(0);
    });
});

// ═══════════════════════════════════════════════════════════════════════════════
// FILE-BASED RPC SUPPORT (for FTPS Integration)
// Watches for request files and writes response files
// ═══════════════════════════════════════════════════════════════════════════════

const FILE_RPC_CONFIG = {
    enabled: process.env.UQFF_FILE_RPC === 'true',
    requestDir: process.env.UQFF_REQUEST_DIR || './uqff_data/requests',
    responseDir: process.env.UQFF_RESPONSE_DIR || './uqff_data/responses',
    pollInterval: parseInt(process.env.UQFF_POLL_INTERVAL) || 1000  // ms
};

/**
 * Process a file-based RPC request.
 * @param {string} requestPath - Path to request JSON file
 */
async function processFileRequest(requestPath) {
    try {
        const requestData = JSON.parse(fs.readFileSync(requestPath, 'utf8'));
        const requestId = requestData.request_id || path.basename(requestPath, '.json');
        
        console.log(`[File-RPC] Processing request: ${requestId}`);
        
        // Compute result
        const result = computeFull(requestData.params || {});
        
        // Write response
        const response = {
            request_id: requestId,
            timestamp: new Date().toISOString(),
            status: 'success',
            result: result
        };
        
        // Ensure response directory exists
        if (!fs.existsSync(FILE_RPC_CONFIG.responseDir)) {
            fs.mkdirSync(FILE_RPC_CONFIG.responseDir, { recursive: true });
        }
        
        const responsePath = path.join(FILE_RPC_CONFIG.responseDir, `resp_${requestId}.json`);
        fs.writeFileSync(responsePath, JSON.stringify(response, null, 2));
        
        // Remove processed request
        fs.unlinkSync(requestPath);
        
        console.log(`[File-RPC] Response written: ${responsePath}`);
        
    } catch (err) {
        console.error(`[File-RPC] Error processing ${requestPath}:`, err.message);
    }
}

/**
 * Watch for file-based RPC requests.
 */
function startFileRPCWatcher() {
    if (!FILE_RPC_CONFIG.enabled) {
        return;
    }
    
    console.log('[File-RPC] Starting file watcher...');
    console.log(`[File-RPC] Request dir: ${FILE_RPC_CONFIG.requestDir}`);
    console.log(`[File-RPC] Response dir: ${FILE_RPC_CONFIG.responseDir}`);
    
    // Ensure directories exist
    if (!fs.existsSync(FILE_RPC_CONFIG.requestDir)) {
        fs.mkdirSync(FILE_RPC_CONFIG.requestDir, { recursive: true });
    }
    if (!fs.existsSync(FILE_RPC_CONFIG.responseDir)) {
        fs.mkdirSync(FILE_RPC_CONFIG.responseDir, { recursive: true });
    }
    
    // Poll for new request files
    setInterval(() => {
        try {
            const files = fs.readdirSync(FILE_RPC_CONFIG.requestDir);
            const requestFiles = files.filter(f => f.startsWith('req_') && f.endsWith('.json'));
            
            for (const file of requestFiles) {
                const requestPath = path.join(FILE_RPC_CONFIG.requestDir, file);
                processFileRequest(requestPath);
            }
        } catch (err) {
            // Directory may not exist yet
        }
    }, FILE_RPC_CONFIG.pollInterval);
}

// Start file-based RPC if enabled
startFileRPCWatcher();

module.exports = { server, computeFull, CONSTANTS, SYSTEMS, FILE_RPC_CONFIG };
