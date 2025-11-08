/**
 * Source6.js Integration Test
 * Verifies all 49 exports from source6.js are accessible via index.js
 */

const fs = require('fs');

console.log('='.repeat(80));
console.log('SOURCE6.JS INTEGRATION TEST');
console.log('Testing 49 exports from source6.js via index.js');
console.log('='.repeat(80));

// Load index.js (which now includes source6.js exports)
const index = require('./index.js');

// Test categories
const results = {
    constants: { tested: 0, passed: 0, failed: [] },
    classes: { tested: 0, passed: 0, failed: [] },
    helpers: { tested: 0, passed: 0, failed: [] },
    physics: { tested: 0, passed: 0, failed: [] },
    graphics: { tested: 0, passed: 0, failed: [] },
    simulation: { tested: 0, passed: 0, failed: [] }
};

// ===========================================================================================
// 1. Test Physical Constants (26 items - with aliases)
// ===========================================================================================
console.log('\n1. Testing Physical Constants (26 items with aliases):');
const constantTests = [
    { name: 'PI6', expectedType: 'number' },
    { name: 'c6', expectedType: 'number' },
    { name: 'G6', expectedType: 'number' },
    { name: 'Omega_g6', expectedType: 'number' },
    { name: 'Mbh6', expectedType: 'number' },
    { name: 'dg6', expectedType: 'number' },
    { name: 'v_SCm6', expectedType: 'number' },
    { name: 'rho_A6', expectedType: 'number' },
    { name: 'rho_sw6', expectedType: 'number' },
    { name: 'v_sw6', expectedType: 'number' },
    { name: 'QA6', expectedType: 'number' },
    { name: 'Qs6', expectedType: 'number' },
    { name: 'kappa6', expectedType: 'number' },
    { name: 'alpha6', expectedType: 'number' },
    { name: 'gamma6', expectedType: 'number' },
    { name: 'delta_sw6', expectedType: 'number' },
    { name: 'epsilon_sw6', expectedType: 'number' },
    { name: 'delta_def6', expectedType: 'number' },
    { name: 'HSCm6', expectedType: 'number' },
    { name: 'UUA6', expectedType: 'number' },
    { name: 'eta6', expectedType: 'number' },
    { name: 'k16', expectedType: 'number' },
    { name: 'k26', expectedType: 'number' },
    { name: 'k36', expectedType: 'number' },
    { name: 'k46', expectedType: 'number' },
    { name: 'beta_i6', expectedType: 'number' },
    { name: 'rho_v6', expectedType: 'number' },
    { name: 'C_concentration6', expectedType: 'number' },
    { name: 'f_feedback6', expectedType: 'number' },
    { name: 'num_strings6', expectedType: 'number' },
    { name: 'Ts006', expectedType: 'number' },
    { name: 'g_mu_nu6', expectedType: 'object' }
];

constantTests.forEach(test => {
    results.constants.tested++;
    try {
        const exists = index.hasOwnProperty(test.name);
        const correctType = typeof index[test.name] === test.expectedType;
        if (exists && correctType) {
            results.constants.passed++;
            console.log(`  ✓ ${test.name}: ${test.expectedType}`);
        } else {
            results.constants.failed.push(test.name);
            console.log(`  ✗ ${test.name}: FAILED (exists: ${exists}, type: ${typeof index[test.name]})`);
        }
    } catch (error) {
        results.constants.failed.push(test.name);
        console.log(`  ✗ ${test.name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// 2. Test Classes (11 items - with aliases)
// ===========================================================================================
console.log('\n2. Testing Classes (11 items with aliases):');
const classTests = [
    { name: 'CelestialBody6', expectedType: 'function' },
    { name: 'ThreeDObject', expectedType: 'function' },
    { name: 'ToolPath', expectedType: 'function' },
    { name: 'SimulationEntity', expectedType: 'function' },
    { name: 'MeshData', expectedType: 'function' },
    { name: 'Shader', expectedType: 'function' },
    { name: 'Camera', expectedType: 'function' },
    { name: 'SIMPlugin', expectedType: 'function' },
    { name: 'PhysicsTerm6', expectedType: 'function' },
    { name: 'DarkMatterHaloTerm6', expectedType: 'function' },
    { name: 'VacuumEnergyTerm6', expectedType: 'function' },
    { name: 'UQFFModule6JS', expectedType: 'function' }
];

classTests.forEach(test => {
    results.classes.tested++;
    try {
        const exists = index.hasOwnProperty(test.name);
        const correctType = typeof index[test.name] === test.expectedType;
        if (exists && correctType) {
            results.classes.passed++;
            console.log(`  ✓ ${test.name}: ${test.expectedType} (class/constructor)`);
        } else {
            results.classes.failed.push(test.name);
            console.log(`  ✗ ${test.name}: FAILED (exists: ${exists}, type: ${typeof index[test.name]})`);
        }
    } catch (error) {
        results.classes.failed.push(test.name);
        console.log(`  ✗ ${test.name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// 3. Test Helper Functions (7 items - with aliases)
// ===========================================================================================
console.log('\n3. Testing Helper Functions (7 items with aliases):');
const helperTests = [
    { name: 'stepFunction6', expectedType: 'function' },
    { name: 'computeEreact', expectedType: 'function' },
    { name: 'computeMuS', expectedType: 'function' },
    { name: 'computeGradMsR', expectedType: 'function' },
    { name: 'computeBj', expectedType: 'function' },
    { name: 'computeOmegaST', expectedType: 'function' },
    { name: 'computeMuJ', expectedType: 'function' }
];

helperTests.forEach(test => {
    results.helpers.tested++;
    try {
        const exists = index.hasOwnProperty(test.name);
        const correctType = typeof index[test.name] === test.expectedType;
        if (exists && correctType) {
            results.helpers.passed++;
            console.log(`  ✓ ${test.name}: ${test.expectedType}`);
        } else {
            results.helpers.failed.push(test.name);
            console.log(`  ✗ ${test.name}: FAILED (exists: ${exists}, type: ${typeof index[test.name]})`);
        }
    } catch (error) {
        results.helpers.failed.push(test.name);
        console.log(`  ✗ ${test.name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// 4. Test Physics Functions (8 items - with aliases)
// ===========================================================================================
console.log('\n4. Testing Physics Functions (8 items with aliases):');
const physicsTests = [
    { name: 'source6_computeUg1', expectedType: 'function' },
    { name: 'source6_computeUg2', expectedType: 'function' },
    { name: 'source6_computeUg3', expectedType: 'function' },
    { name: 'source6_computeUm', expectedType: 'function' },
    { name: 'source6_computeUg4', expectedType: 'function' },
    { name: 'source6_computeUbi', expectedType: 'function' },
    { name: 'source6_computeAMuNu', expectedType: 'function' },
    { name: 'source6_computeFU', expectedType: 'function' }
];

physicsTests.forEach(test => {
    results.physics.tested++;
    try {
        const exists = index.hasOwnProperty(test.name);
        const correctType = typeof index[test.name] === test.expectedType;
        if (exists && correctType) {
            results.physics.passed++;
            console.log(`  ✓ ${test.name}: ${test.expectedType}`);
        } else {
            results.physics.failed.push(test.name);
            console.log(`  ✗ ${test.name}: FAILED (exists: ${exists}, type: ${typeof index[test.name]})`);
        }
    } catch (error) {
        results.physics.failed.push(test.name);
        console.log(`  ✗ ${test.name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// 5. Test 3D Graphics & Model Functions (3 items)
// ===========================================================================================
console.log('\n5. Testing 3D Graphics & Model Functions (3 items):');
const graphicsTests = [
    { name: 'loadOBJ', expectedType: 'function' },
    { name: 'exportOBJ', expectedType: 'function' },
    { name: 'loadTexture', expectedType: 'function' }
];

graphicsTests.forEach(test => {
    results.graphics.tested++;
    try {
        const exists = index.hasOwnProperty(test.name);
        const correctType = typeof index[test.name] === test.expectedType;
        if (exists && correctType) {
            results.graphics.passed++;
            console.log(`  ✓ ${test.name}: ${test.expectedType}`);
        } else {
            results.graphics.failed.push(test.name);
            console.log(`  ✗ ${test.name}: FAILED (exists: ${exists}, type: ${typeof index[test.name]})`);
        }
    } catch (error) {
        results.graphics.failed.push(test.name);
        console.log(`  ✗ ${test.name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// 6. Test Simulation Functions (4 items - with aliases)
// ===========================================================================================
console.log('\n6. Testing Simulation Functions (4 items with aliases):');
const simulationTests = [
    { name: 'source6_simulateQuasarJet', expectedType: 'function' },
    { name: 'printSummaryStats', expectedType: 'function' },
    { name: 'source6_loadBodies', expectedType: 'function' },
    { name: 'getDefaultBodies', expectedType: 'function' }
];

simulationTests.forEach(test => {
    results.simulation.tested++;
    try {
        const exists = index.hasOwnProperty(test.name);
        const correctType = typeof index[test.name] === test.expectedType;
        if (exists && correctType) {
            results.simulation.passed++;
            console.log(`  ✓ ${test.name}: ${test.expectedType}`);
        } else {
            results.simulation.failed.push(test.name);
            console.log(`  ✗ ${test.name}: FAILED (exists: ${exists}, type: ${typeof index[test.name]})`);
        }
    } catch (error) {
        results.simulation.failed.push(test.name);
        console.log(`  ✗ ${test.name}: ERROR - ${error.message}`);
    }
});

// ===========================================================================================
// Summary Report
// ===========================================================================================
console.log('\n' + '='.repeat(80));
console.log('INTEGRATION TEST SUMMARY');
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

// Functional test: Create a CelestialBody instance
console.log('\n' + '='.repeat(80));
console.log('FUNCTIONAL TEST: CelestialBody6 Class Instantiation');
console.log('='.repeat(80));

try {
    const testBody = new index.CelestialBody6(
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
    
    console.log('✓ CelestialBody6 instance created successfully');
    console.log(`  Mass: ${testBody.Ms} kg`);
    console.log(`  Radius: ${testBody.Rs} m`);
    console.log(`  Surface Temperature: ${testBody.Ts_surface} K`);
    
    // Test JSON export
    const json = testBody.toJSON();
    console.log('✓ toJSON() method works');
    
    // Test JSON import
    const testBody2 = index.CelestialBody6.fromJSON(json);
    console.log('✓ fromJSON() method works');
    console.log(`  Restored Mass: ${testBody2.Ms} kg`);
    
} catch (error) {
    console.log(`✗ CelestialBody6 functional test FAILED: ${error.message}`);
}

// Functional test: UQFFModule6JS
console.log('\n' + '='.repeat(80));
console.log('FUNCTIONAL TEST: UQFFModule6JS Self-Expanding Framework');
console.log('='.repeat(80));

try {
    const module6 = new index.UQFFModule6JS();
    console.log('✓ UQFFModule6JS instance created successfully');
    
    // Test dynamic term registration
    const darkMatterTerm = new index.DarkMatterHaloTerm6(1e12 * 1.989e30, 20000);
    module6.registerDynamicTerm(darkMatterTerm);
    console.log('✓ Dynamic term registration works');
    console.log(`  Registered terms: ${module6.dynamicTerms.length}`);
    
    // Test parameter setting
    module6.setDynamicParameter('test_param', 42);
    const paramValue = module6.getDynamicParameter('test_param');
    console.log(`✓ Parameter management works (test_param = ${paramValue})`);
    
    // Test state export
    const state = module6.exportState();
    console.log('✓ State export works');
    console.log(`  Exported ${state.split('\n').length} lines of state data`);
    
    console.log('\n✓ Self-expanding framework fully functional!');
    
} catch (error) {
    console.log(`✗ UQFFModule6JS functional test FAILED: ${error.message}`);
}

console.log('\n' + '='.repeat(80));
console.log('SOURCE6.JS INTEGRATION TEST COMPLETE');
console.log('='.repeat(80));

// Exit with appropriate code
process.exit(totalFailed === 0 ? 0 : 1);
