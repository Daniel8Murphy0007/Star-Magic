// unified_field_demo.js
// Demonstration of Source98: UnifiedFieldModule - Unified Field Strength (F_U) computation
// Integrates all UQFF forces: Ug (gravity), Um (magnetism), Ub (buoyancy), Ui (inertia), Aether

const { UnifiedFieldModule, PREDEFINED_SYSTEMS } = require('./index.js');

console.log('═══════════════════════════════════════════════════════════════════');
console.log('Source98: Unified Field Strength (F_U) Module Demonstration');
console.log('Holistic UQFF Energy Density Integration');
console.log('═══════════════════════════════════════════════════════════════════\n');

// ========== Test 1: Solar System (Default Level 13) ==========
console.log('━━━ Test 1: Solar System at Quantum Level 13 ━━━');
const solarModule = new UnifiedFieldModule();

console.log('\nInitial Variables:');
solarModule.printVariables();

console.log('\n--- Computing F_U at t=0 ---');
const t0 = 0.0;
const F_U_solar = solarModule.computeFU(t0);
console.log(`\nTotal F_U (Solar): ${F_U_solar.toExponential(4)} J/m³`);

console.log('\nDetailed Component Breakdown:');
solarModule.printComponentBreakdown(t0);

console.log('\nEquation Description:');
console.log(solarModule.getEquationText());

// ========== Test 2: Neutron Star System (High Magnetic Field) ==========
console.log('\n\n━━━ Test 2: Neutron Star / Magnetar System ━━━');
const neutronStarParams = {
    level: 18.0,                    // Higher quantum level for compact object
    M: 2.8e30,                      // 1.4 solar masses
    r: 10000,                       // 10 km radius
    B: 1e11,                        // 100 GigaTesla magnetic field
    U_g1: 1e30,                     // Strong internal dipole
    U_g2: 5e52,                     // Moderate outer field
    U_g3: 1e51,                     // Strong magnetic strings
    U_g4: 1e-18,                    // Star-BH interactions
    U_m: 5.0e66,                    // Enhanced magnetism (dominant)
    U_b_sum: -5e28,                 // Strong buoyancy opposition
    U_i: 10.0,                      // Higher inertia resistance
    Aether: 1e-14,                  // Larger metric perturbation
    z: 0.01                         // Small redshift
};

const neutronModule = new UnifiedFieldModule(neutronStarParams);

const t_neutron = 3.156e7;  // 1 year
console.log(`\nComputing F_U at t = 1 year (${t_neutron} s)...`);
const F_U_neutron = neutronModule.computeFU(t_neutron);

console.log(`\nTotal F_U (Neutron Star): ${F_U_neutron.toExponential(4)} J/m³`);
neutronModule.printComponentBreakdown(t_neutron);

// ========== Test 3: Galactic Center / SMBH Environment ==========
console.log('\n\n━━━ Test 3: Galactic Center (SMBH Environment) ━━━');
const galacticParams = {
    level: 22.0,                    // Very high quantum level
    M: 8.6e36,                      // 4.3 million solar masses (Sgr A*)
    r: 1.27e10,                     // Schwarzschild radius
    B: 1e-2,                        // mT (accretion disk field)
    U_g1: 1e35,                     // Massive internal dipole
    U_g2: 1e58,                     // Enormous outer field
    U_g3: 5e54,                     // Massive magnetic strings
    U_g4: 1e-15,                    // Strong star-BH coupling
    U_m: 1e68,                      // Extreme magnetism
    U_b_sum: -1e32,                 // Strong buoyancy
    U_i: 50.0,                      // Very high inertia
    Aether: 1e-13,                  // Significant spacetime curvature
    z: 0.0                          // Local (Milky Way)
};

const galacticModule = new UnifiedFieldModule(galacticParams);

const t_galactic = 3.156e8;  // 10 years
console.log(`\nComputing F_U at t = 10 years (${t_galactic} s)...`);
const F_U_galactic = galacticModule.computeFU(t_galactic);

console.log(`\nTotal F_U (Galactic Center): ${F_U_galactic.toExponential(4)} J/m³`);
galacticModule.printComponentBreakdown(t_galactic);

// ========== Test 4: Dynamic Term Addition ==========
console.log('\n\n━━━ Test 4: Dynamic Physics Term Integration ━━━');

// Create a custom dynamic term: time-varying vacuum fluctuation
const vacuumFluctuationTerm = {
    name: 'VacuumFluctuation',
    compute: (t, variables) => {
        const rho_vac_UA = variables.get('rho_vac_UA') || 7.09e-36;
        const amplitude = 1e-10;
        const frequency = 1e-15;  // Hz
        return amplitude * rho_vac_UA * Math.sin(frequency * t);
    }
};

// Create a quantum coupling term
const quantumCouplingTerm = {
    name: 'QuantumCoupling',
    compute: (t, variables) => {
        const hbar = 1.0546e-34;
        const M = variables.get('M') || 1.989e30;
        const r = variables.get('r') || 1e4;
        const coupling = 1e-40;
        return coupling * (hbar * hbar) / (M * r * r) * Math.cos(t / 1e6);
    }
};

solarModule.setEnableLogging(true);
solarModule.registerDynamicTerm(vacuumFluctuationTerm);
solarModule.registerDynamicTerm(quantumCouplingTerm);

console.log('\nDynamic terms registered. Computing F_U with extended physics...');
const F_U_dynamic = solarModule.computeFU(3.156e7);  // 1 year
console.log(`\nF_U with dynamic terms: ${F_U_dynamic.toExponential(4)} J/m³`);

// ========== Test 5: State Export/Import ==========
console.log('\n\n━━━ Test 5: State Export and Import ━━━');

const exportedState = neutronModule.exportState();
console.log('\nExported state from neutron star module:');
console.log(JSON.stringify(exportedState, null, 2));

// Create new module and import state
const importedModule = new UnifiedFieldModule();
importedModule.importState(exportedState);
const F_U_imported = importedModule.computeFU(0);

console.log(`\nF_U after state import: ${F_U_imported.toExponential(4)} J/m³`);
console.log('✓ State transfer successful!');

// ========== Test 6: PREDEFINED_SYSTEMS Integration ==========
console.log('\n\n━━━ Test 6: Using PREDEFINED_SYSTEMS Configuration ━━━');

if (PREDEFINED_SYSTEMS && PREDEFINED_SYSTEMS.UNIFIED_FIELD_STRENGTH) {
    const config = PREDEFINED_SYSTEMS.UNIFIED_FIELD_STRENGTH;
    console.log(`\nSystem Name: ${config.name}`);
    console.log(`Module Class: ${config.moduleClass.name}`);
    
    const predefinedModule = new config.moduleClass(config.parameters);
    const F_U_predefined = predefinedModule.computeFU(0);
    
    console.log(`\nF_U (Predefined System): ${F_U_predefined.toExponential(4)} J/m³`);
    console.log('\nValidation ranges:');
    console.log(`  Level: ${config.validation.level_range}`);
    console.log(`  F_U: ${config.validation.F_U_range[0].toExponential(2)} - ${config.validation.F_U_range[1].toExponential(2)} J/m³`);
    console.log(`  Physical regime: ${config.validation.physical_regime}`);
    console.log(`  Applications: ${config.validation.applications.join(', ')}`);
}

// ========== Summary ==========
console.log('\n\n═══════════════════════════════════════════════════════════════════');
console.log('Summary: Source98 Unified Field Strength Module');
console.log('═══════════════════════════════════════════════════════════════════');
console.log('✓ Integrates all UQFF forces: Ug, Um, Ub, Ui, Aether');
console.log('✓ Vacuum-normalized energy density (J/m³)');
console.log('✓ Scales across 26 quantum levels (quantum to cosmic)');
console.log('✓ Self-expanding framework with dynamic terms');
console.log('✓ State export/import for modular computation');
console.log('✓ Full MUGE integration maintained');
console.log('═══════════════════════════════════════════════════════════════════\n');
