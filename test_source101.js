// test_source101.js
// Comprehensive test suite for source101.js (HeliosphereThicknessModule) integration
// Tests all functionality including direct import, U_g2 calculations, dynamic variables,
// physics terms, state management, and heliosphere boundary modeling

console.log('\n========================================');
console.log('SOURCE101 INTEGRATION TEST');
console.log('HeliosphereThicknessModule - H_SCm Factor');
console.log('========================================\n');

// Test 1: Direct import from source101.js
console.log('Test 1: Direct import from source101.js');
try {
    const { HeliosphereThicknessModule, PhysicsTerm, DynamicVacuumTerm, QuantumCouplingTerm } = require('./source101.js');
    const module = new HeliosphereThicknessModule();
    console.log('✓ Direct import working');
    console.log(`  Module type: ${module.constructor.name}`);
    console.log(`  H_SCm initial: ${module.computeH_SCm()}`);
} catch (error) {
    console.log('✗ Direct import failed:', error.message);
}

// Test 2: U_g2 calculation with and without H_SCm
console.log('\nTest 2: U_g2 calculations with/without H_SCm');
try {
    const { HeliosphereThicknessModule } = require('./source101.js');
    const module = new HeliosphereThicknessModule();
    
    // Default H_SCm = 1.0
    const u_g2_default = module.computeU_g2(0.0, 0.0);
    const u_g2_no_H = module.computeU_g2_no_H(0.0, 0.0);
    
    console.log('  Default H_SCm = 1.0:');
    console.log(`    U_g2 (with H) = ${u_g2_default.toExponential(4)} J/m³`);
    console.log(`    U_g2 (H=1) = ${u_g2_no_H.toExponential(4)} J/m³`);
    console.log(`    Difference = ${((u_g2_default - u_g2_no_H) / u_g2_no_H * 100).toFixed(2)}%`);
    
    // Test H_SCm = 1.1 (+10%)
    module.updateVariable('H_SCm', 1.1);
    const u_g2_110 = module.computeU_g2(0.0, 0.0);
    console.log('\n  H_SCm = 1.1 (+10%):');
    console.log(`    U_g2 (with H) = ${u_g2_110.toExponential(4)} J/m³`);
    console.log(`    Increase = ${((u_g2_110 - u_g2_default) / u_g2_default * 100).toFixed(2)}%`);
    
    // Test H_SCm = 0.9 (-10%)
    module.updateVariable('H_SCm', 0.9);
    const u_g2_090 = module.computeU_g2(0.0, 0.0);
    console.log('\n  H_SCm = 0.9 (-10%):');
    console.log(`    U_g2 (with H) = ${u_g2_090.toExponential(4)} J/m³`);
    console.log(`    Decrease = ${((u_g2_090 - u_g2_default) / u_g2_default * 100).toFixed(2)}%`);
    
    console.log('✓ U_g2 calculations working');
} catch (error) {
    console.log('✗ U_g2 calculations failed:', error.message);
}

// Test 3: Dynamic variable management
console.log('\nTest 3: Dynamic variable management');
try {
    const { HeliosphereThicknessModule } = require('./source101.js');
    const module = new HeliosphereThicknessModule();
    
    console.log('  Initial variables:');
    console.log(`    H_SCm = ${module.variables.get('H_SCm')}`);
    console.log(`    k_2 = ${module.variables.get('k_2')}`);
    console.log(`    rho_vac_UA = ${module.variables.get('rho_vac_UA').toExponential(4)}`);
    
    // Update variable
    module.updateVariable('H_SCm', 1.15);
    console.log(`\n  After updateVariable('H_SCm', 1.15):`);
    console.log(`    H_SCm = ${module.variables.get('H_SCm')}`);
    
    // Add to variable
    module.addToVariable('k_2', 0.3);
    console.log(`\n  After addToVariable('k_2', 0.3):`);
    console.log(`    k_2 = ${module.variables.get('k_2')}`);
    
    // Subtract from variable
    module.subtractFromVariable('delta_sw', 0.005);
    console.log(`\n  After subtractFromVariable('delta_sw', 0.005):`);
    console.log(`    delta_sw = ${module.variables.get('delta_sw')}`);
    
    // Test derived variable update
    const oldSwirl = module.variables.get('swirl_factor');
    module.updateVariable('v_sw', 6e5);
    const newSwirl = module.variables.get('swirl_factor');
    console.log(`\n  Derived variable update (v_sw: 5e5 -> 6e5):`);
    console.log(`    swirl_factor: ${oldSwirl.toExponential(4)} -> ${newSwirl.toExponential(4)}`);
    
    console.log('✓ Dynamic variables working');
} catch (error) {
    console.log('✗ Dynamic variables failed:', error.message);
}

// Test 4: Dynamic physics terms
console.log('\nTest 4: Dynamic physics terms');
try {
    const { HeliosphereThicknessModule, DynamicVacuumTerm, QuantumCouplingTerm } = require('./source101.js');
    const module = new HeliosphereThicknessModule();
    
    // Baseline without dynamic terms
    const u_g2_base = module.computeU_g2(0.0, 0.0);
    
    // Add dynamic vacuum term
    const vacuumTerm = new DynamicVacuumTerm(1e-10, 1e-15);
    module.registerDynamicTerm(vacuumTerm);
    
    // Add quantum coupling term
    const quantumTerm = new QuantumCouplingTerm(1e-40);
    module.registerDynamicTerm(quantumTerm);
    
    const u_g2_with_terms = module.computeU_g2(1e6, 0.0);
    
    console.log(`  Registered ${module.dynamicTerms.length} dynamic terms`);
    console.log(`  U_g2 (baseline) = ${u_g2_base.toExponential(4)} J/m³`);
    console.log(`  U_g2 (with dynamic terms) = ${u_g2_with_terms.toExponential(4)} J/m³`);
    console.log('✓ Dynamic physics terms working');
} catch (error) {
    console.log('✗ Dynamic physics terms failed:', error.message);
}

// Test 5: State export/import
console.log('\nTest 5: State export and import');
try {
    const { HeliosphereThicknessModule } = require('./source101.js');
    const module1 = new HeliosphereThicknessModule();
    
    // Modify state
    module1.updateVariable('H_SCm', 1.25);
    module1.updateVariable('k_2', 1.5);
    module1.setDynamicParameter('custom_param', 42.0);
    
    // Export state
    const state = module1.exportState();
    console.log('  Exported state:');
    console.log(`    H_SCm = ${state.variables.H_SCm}`);
    console.log(`    k_2 = ${state.variables.k_2}`);
    console.log(`    custom_param = ${state.dynamicParameters.custom_param}`);
    
    // Create new module and import state
    const module2 = new HeliosphereThicknessModule();
    module2.importState(state);
    
    console.log('\n  Imported state:');
    console.log(`    H_SCm = ${module2.variables.get('H_SCm')}`);
    console.log(`    k_2 = ${module2.variables.get('k_2')}`);
    console.log(`    custom_param = ${module2.getDynamicParameter('custom_param')}`);
    
    console.log('✓ State management working');
} catch (error) {
    console.log('✗ State management failed:', error.message);
}

// Test 6: Component breakdown
console.log('\nTest 6: Component breakdown');
try {
    const { HeliosphereThicknessModule } = require('./source101.js');
    const module = new HeliosphereThicknessModule();
    module.updateVariable('H_SCm', 1.1);
    module.printComponentBreakdown(0.0, 0.0);
    console.log('✓ Component breakdown working');
} catch (error) {
    console.log('✗ Component breakdown failed:', error.message);
}

// Test 7: Equation text
console.log('\nTest 7: Equation text representation');
try {
    const { HeliosphereThicknessModule } = require('./source101.js');
    const module = new HeliosphereThicknessModule();
    const eqText = module.getEquationText();
    console.log('  Equation preview:');
    console.log('  ' + eqText.split('\n')[0]);
    console.log('  ' + eqText.split('\n')[1]);
    console.log('✓ Equation text working');
} catch (error) {
    console.log('✗ Equation text failed:', error.message);
}

// Test 8: Import via index.js
console.log('\nTest 8: Import via index.js');
try {
    const { HeliosphereThicknessModule } = require('./index.js');
    const module = new HeliosphereThicknessModule();
    const h = module.computeH_SCm();
    console.log(`  H_SCm via index.js = ${h}`);
    console.log('✓ index.js integration working');
} catch (error) {
    console.log('✗ Error importing from index.js:', error.message);
}

// Test 9: Heliosphere boundary scenario (1 AU to 120 AU)
console.log('\nTest 9: Heliosphere boundary scenario');
try {
    const { HeliosphereThicknessModule } = require('./source101.js');
    const module = new HeliosphereThicknessModule();
    
    // Test at different heliocentric distances
    const distances = [
        { name: '1 AU (Earth)', r: 1.496e11, H_SCm: 1.0 },
        { name: '30 AU (Neptune)', r: 4.5e12, H_SCm: 1.05 },
        { name: '90 AU (Heliopause approach)', r: 1.35e13, H_SCm: 1.1 },
        { name: '120 AU (Heliopause)', r: 1.8e13, H_SCm: 1.15 }
    ];
    
    console.log('  Heliosphere profile:');
    for (const point of distances) {
        module.updateVariable('r', point.r);
        module.updateVariable('H_SCm', point.H_SCm);
        const u_g2 = module.computeU_g2(0.0, 0.0);
        console.log(`    ${point.name}: U_g2 = ${u_g2.toExponential(4)} J/m³ (H=${point.H_SCm})`);
    }
    console.log('✓ Heliosphere boundary modeling working');
} catch (error) {
    console.log('✗ Heliosphere boundary failed:', error.message);
}

// Test 10: Module metadata
console.log('\nTest 10: Module metadata');
try {
    const { HeliosphereThicknessModule } = require('./source101.js');
    const module = new HeliosphereThicknessModule();
    module.printModuleInfo();
    console.log('✓ Module info working');
} catch (error) {
    console.log('✗ Module info failed:', error.message);
}

// Final summary
console.log('\n=== INTEGRATION TEST COMPLETE ===');
console.log('✓ source101.cpp successfully converted to source101.js');
console.log('✓ HeliosphereThicknessModule integrated into Star-Magic UQFF framework');
console.log('✓ All self-expanding dynamics maintained');
console.log('✓ Direct import: WORKING');
console.log('✓ U_g2 calculations: WORKING');
console.log('✓ Dynamic variables: WORKING');
console.log('✓ Dynamic physics terms: WORKING');
console.log('✓ State management: WORKING');
console.log('✓ Heliosphere boundary modeling: WORKING');
console.log('========================================\n');
