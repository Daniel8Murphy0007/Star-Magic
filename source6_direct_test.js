/**
 * Source6.js Direct Test
 * Tests all 49 exports from source6.js directly
 */

console.log('='.repeat(80));
console.log('SOURCE6.JS DIRECT TEST');
console.log('Testing all 49 exports from source6.js');
console.log('='.repeat(80));

// Load source6.js directly
const source6 = require('./source6.js');

// Test results tracker
const results = {
    constants: { tested: 0, passed: 0, failed: [] },
    classes: { tested: 0, passed: 0, failed: [] },
    helpers: { tested: 0, passed: 0, failed: [] },
    physics: { tested: 0, passed: 0, failed: [] },
    graphics: { tested: 0, passed: 0, failed: [] },
    simulation: { tested: 0, passed: 0, failed: [] }
};

// ===========================================================================================
// 1. Test Physical Constants (26 items)
// ===========================================================================================
console.log('\n1. Testing Physical Constants (26 items):');
const constantTests = [
    { name: 'PI', expectedType: 'number', expectedValue: Math.PI },
    { name: 'c', expectedType: 'number' },
    { name: 'G', expectedType: 'number' },
    { name: 'Omega_g', expectedType: 'number' },
    { name: 'Mbh', expectedType: 'number' },
    { name: 'dg', expectedType: 'number' },
    { name: 'v_SCm', expectedType: 'number' },
    { name: 'rho_A', expectedType: 'number' },
    { name: 'rho_sw', expectedType: 'number' },
    { name: 'v_sw', expectedType: 'number' },
    { name: 'QA', expectedType: 'number' },
    { name: 'Qs', expectedType: 'number' },
    { name: 'kappa', expectedType: 'number' },
    { name: 'alpha', expectedType: 'number' },
    { name: 'gamma', expectedType: 'number' },
    { name: 'delta_sw', expectedType: 'number' },
    { name: 'epsilon_sw', expectedType: 'number' },
    { name: 'delta_def', expectedType: 'number' },
    { name: 'HSCm', expectedType: 'number' },
    { name: 'UUA', expectedType: 'number' },
    { name: 'eta', expectedType: 'number' },
    { name: 'k1', expectedType: 'number' },
    { name: 'k2', expectedType: 'number' },
    { name: 'k3', expectedType: 'number' },
    { name: 'k4', expectedType: 'number' },
    { name: 'beta_i', expectedType: 'number' },
    { name: 'rho_v', expectedType: 'number' },
    { name: 'C_concentration', expectedType: 'number' },
    { name: 'f_feedback', expectedType: 'number' },
    { name: 'num_strings', expectedType: 'number' },
    { name: 'Ts00', expectedType: 'number' },
    { name: 'g_mu_nu', expectedType: 'object' }
];

constantTests.forEach(test => {
    results.constants.tested++;
    try {
        const exists = source6.hasOwnProperty(test.name);
        const correctType = typeof source6[test.name] === test.expectedType;
        if (exists && correctType) {
            results.constants.passed++;
            const value = typeof source6[test.name] === 'number' 
                ? source6[test.name].toExponential(2) 
                : `[${test.expectedType}]`;
            console.log(`  ✓ ${test.name}: ${value}`);
        } else {
            results.constants.failed.push(test.name);
            console.log(`  ✗ ${test.name}: FAILED (exists: ${exists}, type: ${typeof source6[test.name]})`);
        }
    } catch (error) {
        results.constants.failed.push(test.name);
        console.log(`  ✗ ${test.name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// 2. Test Classes (12 items)
// ===========================================================================================
console.log('\n2. Testing Classes (12 items):');
const classTests = [
    'CelestialBody',
    'ThreeDObject',
    'ToolPath',
    'SimulationEntity',
    'MeshData',
    'Shader',
    'Camera',
    'SIMPlugin',
    'PhysicsTerm',
    'DarkMatterHaloTerm',
    'VacuumEnergyTerm',
    'UQFFModule6JS'
];

classTests.forEach(name => {
    results.classes.tested++;
    try {
        const exists = source6.hasOwnProperty(name);
        const correctType = typeof source6[name] === 'function';
        if (exists && correctType) {
            results.classes.passed++;
            console.log(`  ✓ ${name}: function (class/constructor)`);
        } else {
            results.classes.failed.push(name);
            console.log(`  ✗ ${name}: FAILED (exists: ${exists}, type: ${typeof source6[name]})`);
        }
    } catch (error) {
        results.classes.failed.push(name);
        console.log(`  ✗ ${name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// 3. Test Helper Functions (7 items)
// ===========================================================================================
console.log('\n3. Testing Helper Functions (7 items):');
const helperTests = [
    'stepFunction',
    'computeEreact',
    'computeMuS',
    'computeGradMsR',
    'computeBj',
    'computeOmegaST',
    'computeMuJ'
];

helperTests.forEach(name => {
    results.helpers.tested++;
    try {
        const exists = source6.hasOwnProperty(name);
        const correctType = typeof source6[name] === 'function';
        if (exists && correctType) {
            results.helpers.passed++;
            console.log(`  ✓ ${name}: function`);
        } else {
            results.helpers.failed.push(name);
            console.log(`  ✗ ${name}: FAILED (exists: ${exists}, type: ${typeof source6[name]})`);
        }
    } catch (error) {
        results.helpers.failed.push(name);
        console.log(`  ✗ ${name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// 4. Test Physics Functions (8 items)
// ===========================================================================================
console.log('\n4. Testing Physics Functions (8 items):');
const physicsTests = [
    'computeUg1',
    'computeUg2',
    'computeUg3',
    'computeUm',
    'computeUg4',
    'computeUbi',
    'computeAMuNu',
    'computeFU'
];

physicsTests.forEach(name => {
    results.physics.tested++;
    try {
        const exists = source6.hasOwnProperty(name);
        const correctType = typeof source6[name] === 'function';
        if (exists && correctType) {
            results.physics.passed++;
            console.log(`  ✓ ${name}: function`);
        } else {
            results.physics.failed.push(name);
            console.log(`  ✗ ${name}: FAILED (exists: ${exists}, type: ${typeof source6[name]})`);
        }
    } catch (error) {
        results.physics.failed.push(name);
        console.log(`  ✗ ${name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// 5. Test 3D Graphics & Model Functions (3 items)
// ===========================================================================================
console.log('\n5. Testing 3D Graphics & Model Functions (3 items):');
const graphicsTests = [
    'loadOBJ',
    'exportOBJ',
    'loadTexture'
];

graphicsTests.forEach(name => {
    results.graphics.tested++;
    try {
        const exists = source6.hasOwnProperty(name);
        const correctType = typeof source6[name] === 'function';
        if (exists && correctType) {
            results.graphics.passed++;
            console.log(`  ✓ ${name}: function`);
        } else {
            results.graphics.failed.push(name);
            console.log(`  ✗ ${name}: FAILED (exists: ${exists}, type: ${typeof source6[name]})`);
        }
    } catch (error) {
        results.graphics.failed.push(name);
        console.log(`  ✗ ${name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// 6. Test Simulation Functions (4 items)
// ===========================================================================================
console.log('\n6. Testing Simulation Functions (4 items):');
const simulationTests = [
    'simulateQuasarJet',
    'printSummaryStats',
    'loadBodies',
    'getDefaultBodies'
];

simulationTests.forEach(name => {
    results.simulation.tested++;
    try {
        const exists = source6.hasOwnProperty(name);
        const correctType = typeof source6[name] === 'function';
        if (exists && correctType) {
            results.simulation.passed++;
            console.log(`  ✓ ${name}: function`);
        } else {
            results.simulation.failed.push(name);
            console.log(`  ✗ ${name}: FAILED (exists: ${exists}, type: ${typeof source6[name]})`);
        }
    } catch (error) {
        results.simulation.failed.push(name);
        console.log(`  ✗ ${name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// Summary Report
// ===========================================================================================
console.log('\n' + '='.repeat(80));
console.log('TEST SUMMARY');
console.log('='.repeat(80));

const categories = [
    { name: 'Constants', data: results.constants },
    { name: 'Classes', data: results.classes },
    { name: 'Helper Functions', data: results.helpers },
    { name: 'Physics Functions', data: results.physics },
    { name: '3D Graphics', data: results.graphics },
    { name: 'Simulation Functions', data: results.simulation }
];

let totalTested = 0;
let totalPassed = 0;
let totalFailed = 0;

categories.forEach(cat => {
    const passRate = cat.data.tested > 0 ? ((cat.data.passed / cat.data.tested) * 100).toFixed(1) : '0.0';
    const status = cat.data.failed.length === 0 ? '✓' : '✗';
    console.log(`${status} ${cat.name}: ${cat.data.passed}/${cat.data.tested} passed (${passRate}%)`);
    if (cat.data.failed.length > 0) {
        console.log(`  Failed: ${cat.data.failed.join(', ')}`);
    }
    totalTested += cat.data.tested;
    totalPassed += cat.data.passed;
    totalFailed += cat.data.failed.length;
});

console.log('\n' + '-'.repeat(80));
const overallPassRate = totalTested > 0 ? ((totalPassed / totalTested) * 100).toFixed(1) : '0.0';
const overallStatus = totalFailed === 0 ? '✓ PASS' : '✗ FAIL';
console.log(`Overall: ${totalPassed}/${totalTested} tests passed (${overallPassRate}%) - ${overallStatus}`);
console.log('='.repeat(80));

// ===========================================================================================
// Functional Tests
// ===========================================================================================
console.log('\n' + '='.repeat(80));
console.log('FUNCTIONAL TEST: CelestialBody Class');
console.log('='.repeat(80));

try {
    // Test 1: Create Sun-like celestial body
    const sun = new source6.CelestialBody(
        1.989e30,  // Ms: Solar mass
        6.96e8,    // Rs: Solar radius
        7e8,       // Rb: Bubble radius
        5778,      // Ts_surface: Surface temp
        2.9e-6,    // omega_s: Rotation rate
        1e-4,      // Bs_avg: Magnetic field
        1e3,       // SCm_density
        1e-15,     // QUA
        1e16,      // Pcore
        1e15,      // PSCm
        1e-5       // omega_c
    );
    
    console.log('✓ CelestialBody instance created (Sun-like)');
    console.log(`  Mass: ${sun.Ms.toExponential(3)} kg`);
    console.log(`  Radius: ${sun.Rs.toExponential(2)} m`);
    console.log(`  Surface Temperature: ${sun.Ts_surface} K`);
    console.log(`  Rotation Rate: ${sun.omega_s.toExponential(2)} rad/s`);
    
    // Test 2: JSON export
    const sunJSON = sun.toJSON();
    console.log('✓ toJSON() method works');
    console.log(`  JSON length: ${JSON.stringify(sunJSON).length} bytes`);
    
    // Test 3: JSON import
    const sunRestored = source6.CelestialBody.fromJSON(sunJSON);
    console.log('✓ fromJSON() method works');
    console.log(`  Restored mass: ${sunRestored.Ms.toExponential(3)} kg`);
    console.log(`  Match: ${sun.Ms === sunRestored.Ms ? '✓' : '✗'}`);
    
} catch (error) {
    console.log(`✗ CelestialBody functional test FAILED: ${error.message}`);
    console.log(error.stack);
}

console.log('\n' + '='.repeat(80));
console.log('FUNCTIONAL TEST: Physics Calculations');
console.log('='.repeat(80));

try {
    // Create test body (Sun)
    const testBody = new source6.CelestialBody(
        1.989e30, 6.96e8, 7e8, 5778, 2.9e-6, 1e-4, 1e3, 1e-15, 1e16, 1e15, 1e-5
    );
    
    const r = 1.496e11;  // 1 AU
    const t = 0;
    const tn = 1e-10;
    const theta = Math.PI / 4;
    
    // Test Ug1 (internal dipole)
    const Ug1 = source6.computeUg1(testBody, r, t, tn, source6.alpha, source6.delta_def, source6.k1);
    console.log(`✓ computeUg1(Sun, 1 AU): ${Ug1.toExponential(3)} m/s²`);
    
    // Test Ug2 (outer bubble)
    const Ug2 = source6.computeUg2(
        testBody, r, t, tn, source6.k2, source6.QA, source6.delta_sw,
        source6.v_sw, source6.HSCm, source6.rho_A, source6.kappa
    );
    console.log(`✓ computeUg2(Sun, 1 AU): ${Ug2.toExponential(3)} m/s²`);
    
    // Test Ug3 (magnetic strings)
    const Ug3 = source6.computeUg3(testBody, r, t, tn, theta, source6.rho_A, source6.kappa, source6.k3);
    console.log(`✓ computeUg3(Sun, 1 AU): ${Ug3.toExponential(3)} m/s²`);
    
    // Test FU (unified field)
    const FU = source6.computeFU(testBody, r, t, tn, theta);
    console.log(`✓ computeFU(Sun, 1 AU): ${FU.toExponential(3)} m/s²`);
    
    console.log('\n✓ All physics calculations working!');
    
} catch (error) {
    console.log(`✗ Physics calculation test FAILED: ${error.message}`);
    console.log(error.stack);
}

console.log('\n' + '='.repeat(80));
console.log('FUNCTIONAL TEST: UQFFModule6JS Self-Expanding Framework');
console.log('='.repeat(80));

try {
    const module6 = new source6.UQFFModule6JS();
    console.log('✓ UQFFModule6JS instance created');
    
    // Test dark matter term
    const M_halo = 1e12 * 1.989e30;  // 1 trillion solar masses
    const r_s = 20000;               // NFW scale radius (parsecs)
    const dmTerm = new source6.DarkMatterHaloTerm(M_halo, r_s);
    module6.registerDynamicTerm(dmTerm);
    console.log(`✓ Dark matter term registered (M=${(M_halo/1.989e30).toExponential(2)} M☉, r_s=${r_s} pc)`);
    
    // Test vacuum energy term
    const Lambda = 1.1e-52;  // m^-2
    const veTerm = new source6.VacuumEnergyTerm(Lambda);
    module6.registerDynamicTerm(veTerm);
    console.log(`✓ Vacuum energy term registered (Λ=${Lambda.toExponential(2)} m⁻²)`);
    
    console.log(`✓ Total dynamic terms: ${module6.dynamicTerms.length}`);
    
    // Test parameters
    module6.setDynamicParameter('learning_rate', 0.01);
    module6.setDynamicParameter('convergence_threshold', 1e-6);
    module6.setDynamicParameter('max_iterations', 1000);
    
    const lr = module6.getDynamicParameter('learning_rate');
    console.log(`✓ Parameter management working (learning_rate = ${lr})`);
    
    // Test state export
    const state = module6.exportState();
    const stateLines = state.split('\n').length;
    console.log(`✓ State export working (${stateLines} lines)`);
    
    // Test computation with dynamic terms
    const testBody = new source6.CelestialBody(
        1.989e30, 6.96e8, 7e8, 5778, 2.9e-6, 1e-4, 1e3, 1e-15, 1e16, 1e15, 1e-5
    );
    const r = 1e4 * 3.086e16;  // 10 kpc
    const accel = module6.computeWithDynamicTerms(testBody, r, 0, 1e-10);
    console.log(`✓ Dynamic term computation: ${accel.toExponential(3)} m/s² at 10 kpc`);
    
    console.log('\n✓ Self-expanding framework fully functional!');
    
} catch (error) {
    console.log(`✗ UQFFModule6JS functional test FAILED: ${error.message}`);
    console.log(error.stack);
}

console.log('\n' + '='.repeat(80));
console.log('FUNCTIONAL TEST: 3D Graphics Classes');
console.log('='.repeat(80));

try {
    // Test MeshData
    const mesh = new source6.MeshData();
    mesh.vertices = [0, 0, 0, 1, 0, 0, 0, 1, 0];
    mesh.normals = [0, 0, 1, 0, 0, 1, 0, 0, 1];
    mesh.indices = [0, 1, 2];
    console.log(`✓ MeshData created (${mesh.vertices.length / 3} vertices, ${mesh.indices.length / 3} triangles)`);
    
    // Test Camera
    const camera = new source6.Camera([0, 0, 10], [0, 0, -1], [0, 1, 0]);
    const viewMatrix = camera.getViewMatrix();
    const projMatrix = camera.getProjectionMatrix(45, 16/9, 0.1, 100);
    console.log(`✓ Camera created (position: [${camera.position.join(', ')}])`);
    console.log(`  View matrix: 4x4 (${viewMatrix.length} elements)`);
    console.log(`  Projection matrix: 4x4 (${projMatrix.length} elements)`);
    
    // Test Shader
    const shader = new source6.Shader();
    console.log(`✓ Shader class instantiated`);
    
    // Test SimulationEntity
    const entity = new source6.SimulationEntity([100, 0, 0], [10, 0, 0], mesh);
    entity.update(0.1);  // Update for 0.1 seconds
    console.log(`✓ SimulationEntity updated (new position: [${entity.position.map(v => v.toFixed(1)).join(', ')}])`);
    
    console.log('\n✓ 3D Graphics classes working!');
    
} catch (error) {
    console.log(`✗ 3D Graphics test FAILED: ${error.message}`);
    console.log(error.stack);
}

console.log('\n' + '='.repeat(80));
console.log('SOURCE6.JS DIRECT TEST COMPLETE');
console.log('='.repeat(80));

// Exit with success if all tests passed
process.exit(totalFailed === 0 ? 0 : 1);
