/**
 * _session281_js_layer_bridge.js
 * ================================
 * Session 281 — Closes UQFF_CALIBRATION_AUDIT Gap #4.
 *
 * JavaScript port of DPMActiveLayerCounter (Session 278) + dynamic layer
 * bridge for index.js calculateUg1/2/3/4 / calculateCompressedGravity.
 *
 * The canonical JS implementations in index.js hardcode `layers = 26`.
 * This module derives N_active(T) from system temperature using the
 * same law as Python (_session278_dpm_layer_counter.py) and C++:
 *
 *     N_active(T) = max(N_FLOOR, min(N_DECOUPLE, (T/T_SCm)^(3/2)))
 *
 * Additive wrapper — index.js canonical functions are NOT modified.
 * Callers opt-in via:
 *     const { dynamicLayers, withDynamicLayers } = require('./_session281_js_layer_bridge');
 *     const N = dynamicLayers(T_keV);
 *     const Ug1 = calculateUg1(r, t, mass, N);
 *
 * Closes Gap #4 (JS hardcoded N=26 inconsistent with Py/C++).
 */

'use strict';

// ---------------------------------------------------------------------------
// Canonical constants — must match _session278_dpm_layer_counter.py
// ---------------------------------------------------------------------------
const T_SCM_KEV   = 0.086;    // 1 MK Aether floor (keV)
const N_FLOOR     = 26;       // base 13+13 DPM stack
const N_DECOUPLE  = 1.0e4;    // decoupling ceiling (Sgr A*, GW remnants)
const PHI         = (1 + Math.sqrt(5)) / 2;

// ---------------------------------------------------------------------------
// Ring meshes (phi^k) — for cross-validation against Python anchor table
// ---------------------------------------------------------------------------
const RING_MESHES_PHI = {
    1: Math.pow(PHI, 0),    // ~ 1
    2: Math.pow(PHI, 3),    // ~ 4.24
    3: Math.pow(PHI, 6),    // ~ 17.94
    4: Math.pow(PHI, 9),    // ~ 76.01
    5: Math.pow(PHI, 12),   // ~ 322.0
    6: Math.pow(PHI, 15),   // ~ 1364
};

// ---------------------------------------------------------------------------
// dynamicLayers(T_keV) — primary export
// ---------------------------------------------------------------------------
function dynamicLayers(T_keV) {
    if (T_keV == null || !isFinite(T_keV) || T_keV <= 0) {
        return N_FLOOR;
    }
    const ratio = T_keV / T_SCM_KEV;
    const N_raw = Math.pow(ratio, 1.5);
    if (N_raw <= N_FLOOR) return N_FLOOR;
    if (N_raw >= N_DECOUPLE) return N_DECOUPLE;
    return N_raw;
}

// ---------------------------------------------------------------------------
// classifyRegime — debugging / parity with Session 278
// ---------------------------------------------------------------------------
function classifyRegime(T_keV) {
    if (T_keV == null || T_keV <= T_SCM_KEV) return 'cool';
    if (T_keV < 1.0)   return 'warm';
    if (T_keV < 10.0)  return 'hot';
    if (T_keV < 100.0) return 'relativistic';
    return 'super';
}

// ---------------------------------------------------------------------------
// nearestRing — find closest phi-mesh ring to a given N
// ---------------------------------------------------------------------------
function nearestRing(N) {
    let best = null;
    let bestDelta = Infinity;
    for (const k of Object.keys(RING_MESHES_PHI)) {
        const v = RING_MESHES_PHI[k];
        const d = Math.abs(Math.log(N / v));
        if (d < bestDelta) { bestDelta = d; best = parseInt(k, 10); }
    }
    return { ring: best, mesh: RING_MESHES_PHI[best] };
}

// ---------------------------------------------------------------------------
// withDynamicLayers(calcFn, T_keV, arity) -> wrapped fn taking (r, t, ...)
// `arity` = total positional parameter count of calcFn including `layers`.
// Defaults to 4 (matches calculateUg1/2/4 signatures). Use 4 for Ug3 too
// (r, theta, t, layers).
// ---------------------------------------------------------------------------
function withDynamicLayers(calcFn, T_keV, arity = 4) {
    const N = dynamicLayers(T_keV);
    return function (...args) {
        const padded = args.slice(0, arity - 1);
        while (padded.length < arity - 1) padded.push(undefined);
        padded.push(N);
        return calcFn(...padded);
    };
}

// ---------------------------------------------------------------------------
// describe — diagnostic dict matching Session 278 output shape
// ---------------------------------------------------------------------------
function describe(T_keV) {
    const N = dynamicLayers(T_keV);
    const ring = nearestRing(N);
    return {
        T_keV: T_keV,
        T_SCM_keV: T_SCM_KEV,
        N_active: N,
        regime: classifyRegime(T_keV),
        ring: ring.ring,
        ring_mesh: ring.mesh,
        floor_hit: N === N_FLOOR,
        ceiling_hit: N === N_DECOUPLE,
    };
}

// ---------------------------------------------------------------------------
// Smoke tests
// ---------------------------------------------------------------------------
function runTests() {
    let pass = 0, fail = 0;
    const check = (name, cond, detail) => {
        if (cond) { pass++; console.log(`  [PASS] ${name}  ${detail || ''}`); }
        else      { fail++; console.log(`  [FAIL] ${name}  ${detail || ''}`); }
    };

    console.log('Session 281 — JS layer bridge smoke tests');
    console.log('-'.repeat(60));

    // L-1: cool plasma -> floor
    check('L-1 Earth crust (T<<T_SCm) -> floor',
          dynamicLayers(1e-6) === N_FLOOR, `N=${dynamicLayers(1e-6)}`);

    // L-2: T = T_SCm -> floor (ratio=1, raw=1, < floor)
    check('L-2 T=T_SCm -> floor',
          dynamicLayers(T_SCM_KEV) === N_FLOOR);

    // L-3: Solar corona (T~0.1 keV) -> floor (ratio^1.5 = 1.26 < 26)
    check('L-3 Solar corona -> floor',
          dynamicLayers(0.1) === N_FLOOR,
          `N=${dynamicLayers(0.1).toFixed(2)}`);

    // L-4: Perseus core (T=4 keV) -> ~317
    const N_perseus = dynamicLayers(4.0);
    const expected_perseus = Math.pow(4.0 / T_SCM_KEV, 1.5);
    check('L-4 Perseus core ~ 317',
          Math.abs(N_perseus - expected_perseus) < 1e-6,
          `N=${N_perseus.toFixed(2)} expected=${expected_perseus.toFixed(2)}`);

    // L-5: Sgr A* (T=50 keV) -> ceiling-bound
    const N_sgrA = dynamicLayers(50.0);
    check('L-5 Sgr A* near or at ceiling',
          N_sgrA >= 1000, `N=${N_sgrA.toFixed(0)}`);

    // L-6: GW150914 (T=100 keV) -> ceiling
    check('L-6 GW150914 capped at N_DECOUPLE',
          dynamicLayers(1e6) === N_DECOUPLE);

    // L-7: ring matching for Perseus (~Ring 5)
    const ring_perseus = nearestRing(N_perseus).ring;
    check('L-7 Perseus matches Ring 5 (phi^12)',
          ring_perseus === 5, `ring=${ring_perseus}`);

    // L-8: regime classification
    check('L-8 cool regime', classifyRegime(1e-6) === 'cool');
    check('L-9 hot regime',  classifyRegime(4.0)  === 'hot');
    check('L-10 super regime', classifyRegime(1000) === 'super');

    // L-11: withDynamicLayers wrapper sanity
    // Dummy function with 4-arg signature mirroring calculateUg1
    const dummy = (r, t, M = 1, layers = 26) => layers;
    const wrapped = withDynamicLayers(dummy, 4.0);  // Perseus
    const N_via_wrapper = wrapped(1e10, 1.0, 1e30);
    check('L-11 withDynamicLayers passes N correctly',
          Math.abs(N_via_wrapper - expected_perseus) < 1e-6,
          `N_passed=${N_via_wrapper.toFixed(2)}`);

    // L-12: describe() returns full diagnostic dict
    const d = describe(4.0);
    check('L-12 describe() Perseus diagnostic',
          d.regime === 'hot' && d.ring === 5 && !d.floor_hit && !d.ceiling_hit,
          `regime=${d.regime} ring=${d.ring}`);

    // L-13: parity check vs Python S278 anchor table
    // Earth crust T=1e-7 keV (3e-6 K) -> floor
    check('L-13 Earth_crust anchor -> floor',
          dynamicLayers(1e-7) === N_FLOOR);

    // L-14: Crab nebula filament T=0.001 keV (~10 MK) -> floor
    check('L-14 Crab filament -> floor',
          dynamicLayers(0.001) === N_FLOOR);

    // L-15: GW170817 kilonova T=10 keV -> hot, finite N
    const N_kilonova = dynamicLayers(10.0);
    check('L-15 GW170817 kilonova finite N (warm-hot transition)',
          N_kilonova > N_FLOOR && N_kilonova < N_DECOUPLE,
          `N=${N_kilonova.toFixed(0)}`);

    console.log('-'.repeat(60));
    console.log(`Results: ${pass} PASS, ${fail} FAIL`);
    return fail === 0 ? 0 : 1;
}

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
module.exports = {
    // constants
    T_SCM_KEV,
    N_FLOOR,
    N_DECOUPLE,
    PHI,
    RING_MESHES_PHI,
    // functions
    dynamicLayers,
    classifyRegime,
    nearestRing,
    withDynamicLayers,
    describe,
    runTests,
};

// Auto-run tests when invoked as script
if (require.main === module) {
    process.exit(runTests());
}
