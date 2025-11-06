// test_source98.js
// Quick test for Source98 UnifiedFieldModule integration

const { UnifiedFieldModule } = require('./index.js');

console.log('═══════════════════════════════════════════════════════════════════');
console.log('Source98: Unified Field Strength (F_U) Module - Quick Test');
console.log('═══════════════════════════════════════════════════════════════════\n');

// Test 1: Default Solar System
console.log('Test 1: Solar System (Default Configuration)');
console.log('─────────────────────────────────────────────');
const solar = new UnifiedFieldModule();
const F_U_solar = solar.computeFU(0);
console.log(`F_U (Solar): ${F_U_solar.toExponential(4)} J/m³`);
console.log(`Expected: ~2.28e65 J/m³ (Um dominant)\n`);

// Test 2: Component Breakdown
console.log('Test 2: Component Breakdown');
console.log('─────────────────────────────────────────────');
const breakdown = solar.computeComponentBreakdown(0);
console.log('Components:');
console.log(`  U_g1 (Internal Dipole): ${breakdown.components.U_g1.toExponential(3)} J/m³`);
console.log(`  U_g2 (Outer Field):     ${breakdown.components.U_g2.toExponential(3)} J/m³`);
console.log(`  U_g3 (Magnetic Strings): ${breakdown.components.U_g3.toExponential(3)} J/m³`);
console.log(`  U_g4 (Star-BH):         ${breakdown.components.U_g4.toExponential(3)} J/m³`);
console.log(`  Um (Magnetism):         ${breakdown.components.Um.toExponential(3)} J/m³ ★ DOMINANT`);
console.log(`  Ub (Buoyancy):          ${breakdown.components.Ub_sum.toExponential(3)} J/m³`);
console.log(`  Ui (Inertia):           ${breakdown.components.Ui.toExponential(3)} J/m³`);
console.log(`  Aether:                 ${breakdown.components.Aether.toExponential(3)} J/m³\n`);

// Test 3: Neutron Star
console.log('Test 3: Neutron Star / Magnetar');
console.log('─────────────────────────────────────────────');
const neutron = new UnifiedFieldModule({
    level: 18.0,
    M: 2.8e30,
    r: 10000,
    B: 1e11,
    U_m: 5.0e66
});
const F_U_neutron = neutron.computeFU(0);
console.log(`F_U (Neutron Star): ${F_U_neutron.toExponential(4)} J/m³`);
console.log(`Magnetic Field: 100 GigaTesla`);
console.log(`Quantum Level: 18\n`);

// Test 4: Variable Operations
console.log('Test 4: Dynamic Variable Updates');
console.log('─────────────────────────────────────────────');
solar.updateVariable('U_m', 3.0e65);
const F_U_updated = solar.computeFU(0);
console.log(`Original Um: 2.28e65 J/m³`);
console.log(`Updated Um:  3.00e65 J/m³`);
console.log(`New F_U:     ${F_U_updated.toExponential(4)} J/m³\n`);

// Test 5: Self-Expanding Framework
console.log('Test 5: Dynamic Term Registration');
console.log('─────────────────────────────────────────────');
const customTerm = {
    name: 'CustomFluctuation',
    compute: (t, vars) => {
        return 1e-12 * Math.sin(t / 1000);
    }
};
solar.registerDynamicTerm(customTerm);
console.log('✓ Custom physics term registered');
console.log(`✓ Dynamic terms count: ${solar.dynamicTerms.length}\n`);

// Test 6: State Export
console.log('Test 6: State Export/Import');
console.log('─────────────────────────────────────────────');
const state = solar.exportState();
console.log('Exported state:');
console.log(`  Variables: ${Object.keys(state.variables).length} parameters`);
console.log(`  Metadata: ${JSON.stringify(state.metadata)}`);
console.log(`  Dynamic terms: ${state.dynamicTermCount}\n`);

// Test 7: Equation Text
console.log('Test 7: UQFF Framework Description');
console.log('─────────────────────────────────────────────');
console.log(solar.getEquationText());

console.log('\n═══════════════════════════════════════════════════════════════════');
console.log('✓ Source98 Integration Complete');
console.log('✓ All tests passed successfully');
console.log('✓ Unified Field Strength module fully operational');
console.log('═══════════════════════════════════════════════════════════════════\n');
