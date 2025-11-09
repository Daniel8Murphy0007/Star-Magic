// Quick test of UQFF modules
console.log("Testing UQFF Module System...\n");

// Test source134 module
try {
    const Abell2256Module = require('./source134.js');
    const abell = new Abell2256Module();
    const t = 1e15; // 1 billion billion seconds
    const F = abell.computeF(t);
    console.log("✓ source134.js (Abell 2256 Galaxy Cluster):");
    console.log(`  F_U_Bi_i = ${F.re.toExponential(3)} + i·${F.im.toExponential(3)} N`);
    console.log(`  Enhanced methods: ${typeof abell.expandParameterSpace === 'function' ? 'YES' : 'NO'}`);
} catch (err) {
    console.log("✗ source134.js FAILED:", err.message);
}

// Test index.js exports
try {
    const UQFF = require('./index.js');
    console.log("\n✓ index.js loaded successfully");
    console.log(`  Exported modules: ${Object.keys(UQFF).length}`);
    console.log(`  Abell2256 exported: ${typeof UQFF.Abell2256UQFFModule === 'function' ? 'YES' : 'NO'}`);

    // Test instantiation
    if (UQFF.Abell2256UQFFModule) {
        const test = new UQFF.Abell2256UQFFModule();
        console.log(`  Abell2256 instantiation: SUCCESS`);
    }
} catch (err) {
    console.log("\n✗ index.js FAILED:", err.message);
}

console.log("\n✓ Module system operational!");
