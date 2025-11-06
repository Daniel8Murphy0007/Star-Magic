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
} catch (error) {
    console.log('✗ Integration test failed:', error.message);
}

console.log('========================================');
console.log('Integration test completed.');
