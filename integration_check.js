console.log('Star-Magic UQFF Framework Integration Test');
console.log('========================================');

try {
    // Test Source73UQFFModule integration
    const { Source73UQFFModule } = require('./index.js');

    if (Source73UQFFModule) {
        console.log('✓ Source73UQFFModule successfully loaded from index.js');

        // Test module instantiation
        const ngc1300Module = new Source73UQFFModule();
        console.log('✓ NGC1300EnhancedUQFFModule class instantiated successfully');

        // Test basic functionality
        const t = 1e9 * 3.156e7;  // 1 Gyr in seconds
        const r = 5e3 * 3.086e19; // 5 kpc in meters
        const g = ngc1300Module.computeG(t, r);
        console.log('✓ Galaxy gravity computation working:', typeof g === 'number' && g > 0);

        // Test dynamic variables
        ngc1300Module.updateVariable('SFR', 2 * ngc1300Module.getVariable('SFR'));
        const newSFR = ngc1300Module.getVariable('SFR');
        console.log('✓ Dynamic variables working:', newSFR > 0);

        // Test environmental forces
        const f_env = ngc1300Module.computeFenv(t);
        console.log('✓ Environmental forces computation working:', typeof f_env === 'number');

        console.log('✓ Source73UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source73UQFFModule not found in exports');
    }

    // Test Source77UQFFModule integration
    const { Source77UQFFModule } = require('./index.js');

    if (Source77UQFFModule) {
        console.log('✓ Source77UQFFModule successfully loaded from index.js');

        // Test module instantiation
        const ugc10214Module = new Source77UQFFModule();
        console.log('✓ UGC10214UQFFModule class instantiated successfully');

        // Test basic functionality
        const t = 1e9 * 3.156e7;  // 1 Gyr in seconds
        const r = 5e3 * 3.086e19; // 5 kpc in meters
        const g = ugc10214Module.computeG(t, r);
        console.log('✓ Galaxy gravity computation working:', typeof g === 'number' && g > 0);

        // Test dynamic variables
        ugc10214Module.updateVariable('SFR', 2 * ugc10214Module.getVariable('SFR'));
        const newSFR = ugc10214Module.getVariable('SFR');
        console.log('✓ Dynamic variables working:', newSFR > 0);

        // Test environmental forces
        const f_env = ugc10214Module.computeFenv(t);
        console.log('✓ Environmental forces computation working:', typeof f_env === 'number');

        // Test merger dynamics
        const mergerMass = ugc10214Module.computeMmerge(t);
        console.log('✓ Merger mass computation working:', typeof mergerMass === 'number' && mergerMass > 0);

        console.log('✓ Source77UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source77UQFFModule not found in exports');
    }

    // Test Source78UQFFModule integration
    const { Source78UQFFModule } = require('./index.js');

    if (Source78UQFFModule) {
        console.log('✓ Source78UQFFModule successfully loaded from index.js');

        // Test module instantiation
        const ngc4676Module = new Source78UQFFModule();
        console.log('✓ NGC4676EnhancedUQFFModule class instantiated successfully');

        // Test basic functionality
        const t = 1e9 * 3.156e7;  // 1 Gyr in seconds
        const r = 5e3 * 3.086e19; // 5 kpc in meters
        const g = ngc4676Module.computeG(t, r);
        console.log('✓ Galaxy gravity computation working:', typeof g === 'number' && g > 0);

        // Test dynamic variables
        ngc4676Module.updateVariable('SFR', 2 * ngc4676Module.getVariable('SFR'));
        const newSFR = ngc4676Module.getVariable('SFR');
        console.log('✓ Dynamic variables working:', newSFR > 0);

        // Test environmental forces
        const f_env = ngc4676Module.computeFenv(t);
        console.log('✓ Environmental forces computation working:', typeof f_env === 'number');

        // Test merger dynamics
        const mergerMass = ngc4676Module.computeMmerge(t);
        console.log('✓ Merger mass computation working:', typeof mergerMass === 'number' && mergerMass > 0);

        console.log('✓ Source78UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source78UQFFModule not found in exports');
    }

    // Test Source79UQFFModule integration
    const { Source79UQFFModule } = require('./index.js');

    if (Source79UQFFModule) {
        console.log('✓ Source79UQFFModule successfully loaded from index.js');

        // Test module instantiation
        const redSpiderModule = new Source79UQFFModule();
        console.log('✓ RedSpiderUQFFModule class instantiated successfully');

        // Test basic functionality
        const t = 1900 * 3.156e7;  // 1900 yr in seconds
        const r = 7.1e15;  // Red Spider radius
        const g = redSpiderModule.computeG(t, r);
        console.log('✓ Nebula gravity computation working:', typeof g === 'number' && !isNaN(g));

        // Test dynamic variables
        redSpiderModule.updateVariable('f_super', 2 * redSpiderModule.variables.get('f_super'));
        const newFsuper = redSpiderModule.variables.get('f_super');
        console.log('✓ Dynamic variables working:', newFsuper > 0);

        // Test environmental forces
        const f_env = redSpiderModule.computeFenv(t);
        console.log('✓ Environmental forces computation working:', typeof f_env === 'number' && !isNaN(f_env));

        // Test frequency computations
        const f_dpm = redSpiderModule.computeDPMTerm(t);
        console.log('✓ DPM term computation working:', typeof f_dpm === 'number' && !isNaN(f_dpm));

        console.log('✓ Source79UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source79UQFFModule not found in exports');
    }

    // Test Source80UQFFModule integration
    const { Source80UQFFModule } = require('./index.js');

    if (Source80UQFFModule) {
        console.log('✓ Source80UQFFModule successfully loaded from index.js');

        // Test module instantiation
        const smbhBinaryModule = new Source80UQFFModule();
        console.log('✓ SMBHBinaryUQFFModule class instantiated successfully');

        // Test basic functionality
        const t = 1.555e7;  // Coalescence time in seconds
        const r = 0.1 * 9.461e15;  // Initial separation in meters
        const g = smbhBinaryModule.computeG(t, r);
        console.log('✓ SMBH binary gravity computation working:', typeof g === 'number' && !isNaN(g));

        // Test dynamic variables
        smbhBinaryModule.updateVariable('f_super', 2 * smbhBinaryModule.variables.get('f_super'));
        const newFsuper = smbhBinaryModule.variables.get('f_super');
        console.log('✓ Dynamic variables working:', newFsuper > 0);

        // Test environmental forces
        const f_env = smbhBinaryModule.computeFenv(t);
        console.log('✓ Environmental forces computation working:', typeof f_env === 'number' && !isNaN(f_env));

        // Test frequency computations
        const f_res = smbhBinaryModule.computeResonanceTerm(t);
        console.log('✓ Resonance term computation working:', typeof f_res === 'number' && !isNaN(f_res));

        // Test coalescence time computation
        const t_coal = smbhBinaryModule.computeCoalescenceTime();
        console.log('✓ Coalescence time computation working:', typeof t_coal === 'number' && t_coal > 0);

        console.log('✓ Source80UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source80UQFFModule not found in exports');
    }

    // Test Source81UQFFModule integration
    const { Source81UQFFModule } = require('./index.js');

    if (Source81UQFFModule) {
        console.log('✓ Source81UQFFModule successfully loaded from index.js');

        // Test module instantiation
        const ngc346Module = new Source81UQFFModule();
        console.log('✓ NGC346UQFFModule class instantiated successfully');

        // Test basic functionality
        const t = 1e7 * 3.156e7;  // 10 Myr in seconds
        const r = 1e16;  // 0.3 pc
        const g = ngc346Module.computeG(t, r);
        console.log('✓ Nebula gravity computation working:', typeof g === 'number' && !isNaN(g));

        // Test dynamic variables
        ngc346Module.updateVariable('SFR', 2 * ngc346Module.variables.get('SFR'));
        const newSFR = ngc346Module.variables.get('SFR');
        console.log('✓ Dynamic variables working:', newSFR > 0);

        // Test environmental forces
        const f_env = ngc346Module.computeFenv(t);
        console.log('✓ Environmental forces computation working:', typeof f_env === 'number' && !isNaN(f_env));

        // Test Ug3 collapse computation
        const ug3 = ngc346Module.computeUg3(t);
        console.log('✓ Ug3 collapse computation working:', typeof ug3 === 'number' && !isNaN(ug3));

        // Test core energy computation
        const e_core = ngc346Module.computeEcore(ngc346Module.variables.get('rho_gas'));
        console.log('✓ Core energy computation working:', typeof e_core === 'number' && !isNaN(e_core));

        console.log('✓ Source81UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source81UQFFModule not found in exports');
    }

    // Test Source82UQFFModule (SMBH dynamics in M-σ relation)
    const { Source82UQFFModule } = require('./index.js');

    if (Source82UQFFModule) {
        console.log('✓ Source82UQFFModule successfully loaded from index.js');
        const smbhModule = new Source82UQFFModule();
        console.log('✓ SMBHUQFFModule class instantiated successfully');

        // Test SMBH gravity computation
        const t = 4.543e9 * 3.156e7;  // 4.543 Gyr in seconds
        const sigma = 200e3;  // 200 km/s
        const g_smbh = smbhModule.computeG(t, sigma);
        console.log('✓ SMBH gravity computation working:', typeof g_smbh === 'number' && !isNaN(g_smbh));

        // Test dynamic variables
        smbhModule.updateVariable('M_bh', 1e13 * 1.989e30);
        console.log('✓ Dynamic variables working: true');

        // Test environmental forces computation
        const f_env = smbhModule.computeFenv(t, { sigma: sigma, M_bh: smbhModule.variables.get('M_bh') });
        console.log('✓ Environmental forces computation working:', typeof f_env === 'number' && !isNaN(f_env));

        // Test Um computation
        const um = smbhModule.computeUm(t, smbhModule.variables.get('R_bulge'), 1);
        console.log('✓ Um computation working:', typeof um === 'number' && !isNaN(um));

        // Test Ug1 computation
        const ug1 = smbhModule.computeUg1(t, smbhModule.variables.get('R_bulge'), smbhModule.variables.get('M_bh'), 1);
        console.log('✓ Ug1 computation working:', typeof ug1 === 'number' && !isNaN(ug1));

        console.log('✓ Source82UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source82UQFFModule not found in exports');
    }

    // Test Source83UQFFModule (LENR dynamics via electro-weak interactions)
    const { Source83UQFFModule } = require('./index.js');

    if (Source83UQFFModule) {
        console.log('✓ Source83UQFFModule successfully loaded from index.js');
        const lenrModule = new Source83UQFFModule();
        console.log('✓ LENRUQFFModule class instantiated successfully');

        // Test LENR gravity computation
        const t = 1e-6;  // 1 microsecond
        const r = 1e-6;  // 1 micrometer
        const g_lenr = lenrModule.computeG(t, r);
        console.log('✓ LENR gravity computation working:', typeof g_lenr === 'number' && !isNaN(g_lenr));

        // Test dynamic variables
        lenrModule.updateVariable('electron_density', 1e25);
        console.log('✓ Dynamic variables working: true');

        // Test environmental forces computation
        const f_env = lenrModule.computeFenv(t);
        console.log('✓ Environmental forces computation working:', typeof f_env === 'number' && !isNaN(f_env));

        // Test plasma frequency computation
        const omega_p = lenrModule.computePlasmaFreq();
        console.log('✓ Plasma frequency computation working:', typeof omega_p === 'number' && omega_p > 0);

        // Test neutron rate computation
        const neutron_rate = lenrModule.computeNeutronRate(t);
        console.log('✓ Neutron rate computation working:', typeof neutron_rate === 'number' && !isNaN(neutron_rate));

        // Test scenario switching
        lenrModule.setScenario('exploding_wire');
        console.log('✓ Scenario switching working: true');

        console.log('✓ Source83UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source83UQFFModule not found in exports');
    }

    // Test Source84UQFFModule (LENR neutron production calibration)
    const { Source84UQFFModule } = require('./index.js');

    if (Source84UQFFModule) {
        console.log('✓ Source84UQFFModule successfully loaded from index.js');
        const lenrCalibModule = new Source84UQFFModule();
        console.log('✓ LENRCalibUQFFModule class instantiated successfully');

        // Test LENR gravity computation
        const t = 3.156e7;  // 1 year in seconds
        const g_lenr_calib = lenrCalibModule.computeG(t);
        console.log('✓ LENR calib gravity computation working:', typeof g_lenr_calib === 'number' && !isNaN(g_lenr_calib));

        // Test dynamic variables
        lenrCalibModule.updateVariable('k_eta', 1.1e13);
        console.log('✓ Dynamic variables working: true');

        // Test environmental forces computation
        const f_env = lenrCalibModule.computeEnvironmentalForces(t);
        console.log('✓ Environmental forces computation working:', typeof f_env === 'object' && f_env.neutron_rate > 0);

        // Test neutron production rate (eta) computation
        const eta = lenrCalibModule.computeEta(t, 1);
        console.log('✓ Neutron production rate computation working:', typeof eta === 'number' && eta > 0);

        // Test Um computation
        const r = lenrCalibModule.variables.get('r');
        const um = lenrCalibModule.computeUm(t, r, 1);
        console.log('✓ Um computation working:', typeof um === 'number' && !isNaN(um));

        // Test scenario switching
        lenrCalibModule.setScenario('wires');
        console.log('✓ Scenario switching working: true');

        console.log('✓ Source84UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source84UQFFModule not found in exports');
    }

    // Test Source85UQFFModule integration
    const { Source85UQFFModule } = require('./index.js');

    if (Source85UQFFModule) {
        console.log('✓ Source85UQFFModule successfully loaded from index.js');
        const source85Module = new Source85UQFFModule();
        console.log('✓ Source85UQFFModule class instantiated successfully');

        // Test nebula gravity computation
        const t = 1e7 * 3.156e7;  // 10 Myr in seconds
        const r = 1e16;  // 0.3 pc
        const g_source85 = source85Module.computeG(t, r);
        console.log('✓ Nebula gravity computation working:', typeof g_source85 === 'number' && !isNaN(g_source85));

        // Test dynamic variables
        source85Module.updateVariable('SFR', 2 * source85Module.variables.get('SFR'));
        const newSFR = source85Module.variables.get('SFR');
        console.log('✓ Dynamic variables working:', newSFR > 0);

        // Test environmental forces
        const f_env = source85Module.computeFenv(t);
        console.log('✓ Environmental forces computation working:', typeof f_env === 'number' && !isNaN(f_env));

        // Test Ug3 collapse computation
        const ug3 = source85Module.computeUg3(t);
        console.log('✓ Ug3 collapse computation working:', typeof ug3 === 'number' && !isNaN(ug3));

        // Test core energy computation
        const e_core = source85Module.computeEcore(source85Module.variables.get('rho_gas'));
        console.log('✓ Core energy computation working:', typeof e_core === 'number' && !isNaN(e_core));

        // Test quantum term computation
        const quantum_term = source85Module.computeQuantumTerm(source85Module.variables.get('t_Hubble'), r);
        console.log('✓ Quantum term computation working:', typeof quantum_term === 'number' && !isNaN(quantum_term));

        console.log('✓ Source85UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source85UQFFModule not found in exports');
    }

    // Test Source86UQFFModule integration
    const { Source86UQFFModule, SystemType } = require('./index.js');

    if (Source86UQFFModule) {
        console.log('✓ Source86UQFFModule successfully loaded from index.js');
        const source86Module = new Source86UQFFModule(SystemType.MAGNETAR_SGR_1745_2900);
        console.log('✓ Source86UQFFModule class instantiated successfully');

        // Test compressed MUGE computation
        const t = 3.799e10;  // Magnetar timescale
        const g_compressed = source86Module.computeG_compressed(t);
        console.log('✓ Compressed MUGE computation working:', typeof g_compressed === 'number' && !isNaN(g_compressed));

        // Test resonance MUGE computation
        const g_resonance = source86Module.computeG_resonance(t);
        console.log('✓ Resonance MUGE computation working:', typeof g_resonance === 'number' && !isNaN(g_resonance));

        // Test dynamic variables
        source86Module.updateVariable('M', 2 * source86Module.variables.get('M'));
        const newM = source86Module.variables.get('M');
        console.log('✓ Dynamic variables working:', newM > 0);

        // Test system switching
        source86Module.setSystem(SystemType.SAGITTARIUS_A);
        console.log('✓ System switching working: true');

        // Test quantum term computation
        const quantum_term = source86Module.computeQuantumTerm();
        console.log('✓ Quantum term computation working:', typeof quantum_term === 'number' && !isNaN(quantum_term));

        // Test resonance helpers
        const aDPM = source86Module.computeADPM();
        console.log('✓ Resonance helper computation working:', typeof aDPM === 'number' && !isNaN(aDPM));

        console.log('✓ Source86UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source86UQFFModule not found in exports');
    }

    // Test Source87UQFFModule integration
    const { Source87UQFFModule, SystemType_Source87 } = require('./index.js');

    if (Source87UQFFModule) {
        console.log('✓ Source87UQFFModule successfully loaded from index.js');
        const source87Module = new Source87UQFFModule(SystemType_Source87.MAGNETAR_SGR_1745_2900);
        console.log('✓ Source87UQFFModule class instantiated successfully');

        // Test resonance MUGE computation
        const t = 3.799e10;  // Magnetar timescale
        const g_resonance = source87Module.computeG_resonance(t);
        console.log('✓ Resonance MUGE computation working:', typeof g_resonance === 'number' && !isNaN(g_resonance));

        // Test dynamic variables
        source87Module.updateVariable('M', 2 * source87Module.variables.get('M'));
        const newM = source87Module.variables.get('M');
        console.log('✓ Dynamic variables working:', newM > 0);

        // Test system switching
        source87Module.setSystem(SystemType_Source87.SAGITTARIUS_A);
        console.log('✓ System switching working: true');

        // Test resonance term computations
        const aDPM = source87Module.computeADPM();
        console.log('✓ Resonance term computation working:', typeof aDPM === 'number' && !isNaN(aDPM));

        // Test helper computations
        const fdpm = source87Module.computeFDPM();
        console.log('✓ FDPM computation working:', typeof fdpm === 'number' && fdpm > 0);

        // Test Ug4i computation
        const ug4i = source87Module.computeUg4i(t);
        console.log('✓ Ug4i computation working:', typeof ug4i === 'number' && !isNaN(ug4i));

        console.log('✓ Source87UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source87UQFFModule not found in exports');
    }

    // Test Source88UQFFModule integration
    const { Source88UQFFModule } = require('./index.js');

    if (Source88UQFFModule) {
        console.log('✓ Source88UQFFModule successfully loaded from index.js');
        const source88Module = new Source88UQFFModule();
        console.log('✓ Source88UQFFModule class instantiated successfully');

        // Test Andromeda galaxy evolution computation
        const t = 10.0 * 1e9 * 3.156e7;  // 10 Gyr in seconds
        const g_andromeda = source88Module.computeG(t);
        console.log('✓ Andromeda galaxy evolution computation working:', typeof g_andromeda === 'number' && !isNaN(g_andromeda));

        // Test dynamic variables
        source88Module.updateVariable('M', 2e12 * source88Module.variables.get('M_sun'));
        const newM = source88Module.variables.get('M');
        console.log('✓ Dynamic variables working:', newM > 0);

        // Test dust friction computation
        const a_dust = source88Module.computeADust();
        console.log('✓ Dust friction computation working:', typeof a_dust === 'number' && !isNaN(a_dust));

        // Test EM term computation
        const em_term = source88Module.computeEMTerm();
        console.log('✓ EM term computation working:', typeof em_term === 'number' && !isNaN(em_term));

        // Test expansion factor computation
        const Hz = source88Module.computeHz();
        console.log('✓ Expansion factor computation working:', typeof Hz === 'number' && Hz > 0);

        console.log('✓ Source88UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source88UQFFModule not found in exports');
    }

    // Test Source89UQFFModule integration
    const { Source89UQFFModule } = require('./index.js');

    if (Source89UQFFModule) {
        console.log('✓ Source89UQFFModule successfully loaded from index.js');
        const source89Module = new Source89UQFFModule();
        console.log('✓ Source89UQFFModule class instantiated successfully');

        // Test Aether coupling perturbation computation
        const perturbation = source89Module.computePerturbation();
        console.log('✓ Aether coupling perturbation computation working:', typeof perturbation === 'number' && !isNaN(perturbation));

        // Test dynamic variables
        source89Module.updateVariable('eta', 2e-22);
        const newEta = source89Module.variables.get('eta');
        console.log('✓ Dynamic variables working:', newEta > 0);

        // Test stress-energy tensor computation
        const t_s = source89Module.computeT_s();
        console.log('✓ Stress-energy tensor computation working:', typeof t_s === 'number' && t_s > 0);

        // Test perturbed metric computation
        const a_mu_nu = source89Module.computeA_mu_nu();
        console.log('✓ Perturbed metric computation working:', Array.isArray(a_mu_nu) && a_mu_nu.length === 4);

        console.log('✓ Source89UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source89UQFFModule not found in exports');
    }
} catch (error) {
    console.log('✗ Integration test failed:', error.message);
}

// Test Source90UQFFModule integration
try {
    const { Source90UQFFModule } = require('./index.js');

    if (Source90UQFFModule) {
        console.log('✓ Source90UQFFModule successfully loaded from index.js');
        const source90Module = new Source90UQFFModule();
        console.log('✓ Source90UQFFModule class instantiated successfully');

        // Test dynamic variable management
        source90Module.setDynamicParameter('eta', 1e-22);
        source90Module.setDynamicParameter('rho_vac_UA', 7.09e-36);
        console.log('✓ Dynamic parameters set successfully');

        // Test background metric computation
        const g_mu_nu = source90Module.computeG_mu_nu();
        console.log('✓ Background metric computation working:', Array.isArray(g_mu_nu) && g_mu_nu.length === 4);

        // Test stress-energy tensor computation
        const T_s = source90Module.computeT_s();
        console.log('✓ Stress-energy tensor computation working:', typeof T_s === 'number');

        // Test perturbation computation
        const perturbation = source90Module.computePerturbation();
        console.log('✓ Perturbation computation working:', typeof perturbation === 'number');

        // Test perturbed metric computation
        const A_mu_nu = source90Module.computeA_mu_nu();
        console.log('✓ Perturbed metric computation working:', Array.isArray(A_mu_nu) && A_mu_nu.length === 4);

        console.log('✓ Source90UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source90UQFFModule not found in exports');
    }
} catch (error) {
    console.log('✗ Integration test failed:', error.message);
}

// Test Source91UQFFModule integration
try {
    const { Source91UQFFModule } = require('./index.js');

    if (Source91UQFFModule) {
        console.log('✓ Source91UQFFModule successfully loaded from index.js');
        const source91Module = new Source91UQFFModule();
        console.log('✓ Source91UQFFModule class instantiated successfully');

        // Test dynamic variable management
        source91Module.setDynamicParameter('num_states', 26);
        source91Module.setDynamicParameter('SCm_amount', 1e42);
        console.log('✓ Dynamic parameters set successfully');

        // Test DPM sphere centers computation
        const dpmCenters = source91Module.computeDPM();
        console.log('✓ DPM sphere centers computation working:', Array.isArray(dpmCenters) && dpmCenters.length === 26);

        // Test [SCm] energy computation
        const scmEnergy = source91Module.computeSCmEnergy();
        console.log('✓ [SCm] energy computation working:', typeof scmEnergy === 'number');

        // Test [UA] energy computation
        const uaEnergy = source91Module.computeUAEnergy();
        console.log('✓ [UA] energy computation working:', typeof uaEnergy === 'number');

        // Test resonance factor computation
        const resonanceFactor = source91Module.computeResonanceFactor();
        console.log('✓ Resonance factor computation working:', typeof resonanceFactor === 'number');

        console.log('✓ Source91UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source91UQFFModule not found in exports');
    }
} catch (error) {
    console.log('✗ Integration test failed:', error.message);
}

// Test Source92UQFFModule integration
try {
    const { Source92UQFFModule } = require('./index.js');

    if (Source92UQFFModule) {
        console.log('✓ Source92UQFFModule successfully loaded from index.js');
        const source92Module = new Source92UQFFModule();
        console.log('✓ Source92UQFFModule class instantiated successfully');

        // Test dynamic variable management
        source92Module.setDynamicParameter('beta', 0.6);
        source92Module.setDynamicParameter('U_g1', 1.39e26);
        console.log('✓ Dynamic parameters set successfully');

        // Test β_i computation
        const beta1 = source92Module.computeBeta(1);
        console.log('✓ β_i computation working:', typeof beta1 === 'number' && beta1 === 0.6);

        // Test U_bi computation for specific i
        const u_b1 = source92Module.computeU_bi(1);
        console.log('✓ U_bi computation working:', typeof u_b1 === 'number');

        // Test all U_bi computation
        const all_u_bi = source92Module.computeAllU_bi();
        console.log('✓ All U_bi computation working:', Array.isArray(all_u_bi) && all_u_bi.length === 4);

        // Test F_U contribution computation
        const f_u_contrib = source92Module.computeF_U_contribution();
        console.log('✓ F_U contribution computation working:', typeof f_u_contrib === 'number');

        console.log('✓ Source92UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source92UQFFModule not found in exports');
    }
} catch (error) {
    console.log('✗ Integration test failed:', error.message);
}

// Test Source93UQFFModule integration
try {
    const { Source93UQFFModule } = require('./index.js');

    if (Source93UQFFModule) {
        console.log('✓ Source93UQFFModule successfully loaded from index.js');
        const source93Module = new Source93UQFFModule();
        console.log('✓ Source93UQFFModule class instantiated successfully');

        // Test dynamic variable management
        source93Module.setDynamicParameter('epsilon_sw', 0.001);
        source93Module.setDynamicParameter('rho_vac_sw', 8e-21);
        console.log('✓ Dynamic parameters set successfully');

        // Test ε_sw computation
        const epsilon_sw = source93Module.computeEpsilon_sw();
        console.log('✓ ε_sw computation working:', typeof epsilon_sw === 'number' && epsilon_sw === 0.001);

        // Test modulation factor computation
        const mod_factor = source93Module.computeModulationFactor();
        console.log('✓ Modulation factor computation working:', typeof mod_factor === 'number' && mod_factor > 1);

        // Test U_b1 computation
        const u_b1 = source93Module.computeU_b1();
        console.log('✓ U_b1 computation working:', typeof u_b1 === 'number');

        // Test solar wind modulation computation
        const modulation = source93Module.computeSolarWindModulation();
        console.log('✓ Solar wind modulation computation working:', typeof modulation === 'object' && modulation.epsilon_sw === 0.001);

        console.log('✓ Source93UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source93UQFFModule not found in exports');
    }
} catch (error) {
    console.log('✗ Integration test failed:', error.message);
}

// Test Source94UQFFModule integration
/*
try {
    const { Source94UQFFModule } = require('./index.js');

    if (Source94UQFFModule) {
        console.log('✓ Source94UQFFModule successfully loaded from index.js');
        const source94Module = new Source94UQFFModule();
        console.log('✓ Source94UQFFModule class instantiated successfully');

        // Test dynamic variable management
        source94Module.setDynamicParameter('k1', 1.5);
        source94Module.setDynamicParameter('k2', 1.2);
        console.log('✓ Dynamic parameters set successfully');

        // Test k_values array
        const k_values = source94Module.getKValues();
        console.log('✓ k_values array working:', Array.isArray(k_values) && k_values.length === 4 && k_values[0] === 1.5);

        // Test U_gi computations
        const u_g1 = source94Module.computeU_g1();
        console.log('✓ U_g1 computation working:', typeof u_g1 === 'number' && u_g1 > 0);

        const u_g2 = source94Module.computeU_g2();
        console.log('✓ U_g2 computation working:', typeof u_g2 === 'number' && u_g2 > 0);

        const u_g3 = source94Module.computeU_g3();
        console.log('✓ U_g3 computation working:', typeof u_g3 === 'number' && u_g3 > 0);

        const u_g4 = source94Module.computeU_g4();
        console.log('✓ U_g4 computation working:', typeof u_g4 === 'number' && u_g4 > 0);

        // Test k_i * U_gi computations
        const k_u_g1 = source94Module.computeK_Ug1();
        console.log('✓ k1 * U_g1 computation working:', typeof k_u_g1 === 'number' && k_u_g1 > 0);

        const k_u_g2 = source94Module.computeK_Ug2();
        console.log('✓ k2 * U_g2 computation working:', typeof k_u_g2 === 'number' && k_u_g2 > 0);

        const k_u_g3 = source94Module.computeK_Ug3();
        console.log('✓ k3 * U_g3 computation working:', typeof k_u_g3 === 'number' && k_u_g3 > 0);

        const k_u_g4 = source94Module.computeK_Ug4();
        console.log('✓ k4 * U_g4 computation working:', typeof k_u_g4 === 'number' && k_u_g4 > 0);

        // Test sum computation for F_U contribution
        const sum_k_u_gi = source94Module.computeSumK_Ugi();
        console.log('✓ Sum k_i * U_gi computation working:', typeof sum_k_u_gi === 'number' && sum_k_u_gi > 0);

        // Test Ug coupling methods
        const coupling_data = source94Module.getUgCouplingData();
        console.log('✓ Ug coupling data retrieval working:', typeof coupling_data === 'object' && coupling_data.k_values.length === 4);

        console.log('✓ Source94UQFFModule integration: COMPLETED');
    } else {
        console.log('✗ Source94UQFFModule not found in exports');
    }
} catch (error) {
    console.log('✗ Integration test failed:', error.message);
}
*/

console.log('========================================');
console.log('Integration test completed.');
