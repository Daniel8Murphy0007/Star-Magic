// Verify all Source13-162 integrations in index.js
const fs = require('fs');

console.log("VERIFYING SOURCE13-162 INTEGRATIONS IN INDEX.JS");
console.log("=".repeat(80));

const indexContent = fs.readFileSync('index.js', 'utf-8');
const lines = indexContent.split('\n');

// Find module.exports line
let exportsStart = -1;
for (let i = 0; i < lines.length; i++) {
    if (lines[i].includes('module.exports = {')) {
        exportsStart = i;
        break;
    }
}

console.log(`module.exports starts at line: ${exportsStart + 1}\n`);

// List of modules to check (Source13-162)
const modulesToCheck = [
    { name: 'MagnetarSGR1745_2900', source: 'Source13' },
    { name: 'UQFFBuoyancyCNBModule', source: 'Source156' },
    { name: 'UQFFBuoyancyModule157', source: 'Source157' },
    { name: 'UQFFBuoyancyModule158', source: 'Source158' },
    { name: 'UQFFBuoyancyModule159', source: 'Source159' },
    { name: 'UQFFBuoyancyModule160', source: 'Source160' },
    { name: 'UQFFBuoyancyModule161', source: 'Source161' },
    { name: 'UQFFBuoyancyCNBModule162', source: 'Source162' }
];

console.log("Checking module definitions and exports:\n");

modulesToCheck.forEach(module => {
    // Find class definition
    const classRegex = new RegExp(`class ${module.name}\\b`);
    let classLine = -1;
    for (let i = 0; i < lines.length; i++) {
        if (classRegex.test(lines[i])) {
            classLine = i;
            break;
        }
    }

    // Find export reference
    const exportRegex = new RegExp(`\\b${module.name}\\b`);
    let exportLine = -1;
    for (let i = exportsStart; i < Math.min(exportsStart + 500, lines.length); i++) {
        if (exportRegex.test(lines[i]) && lines[i].includes(',')) {
            exportLine = i;
            break;
        }
    }

    const status = classLine > 0 && exportLine > 0 && classLine < exportsStart ? '✅' : '❌';
    const defStatus = classLine > 0 ? `Line ${classLine + 1}` : 'NOT FOUND';
    const expStatus = exportLine > 0 ? `Line ${exportLine + 1}` : 'NOT FOUND';
    const orderOK = classLine > 0 && classLine < exportsStart ? 'BEFORE exports' : 'AFTER exports ❌';

    console.log(`${status} ${module.source} - ${module.name}`);
    console.log(`   Definition: ${defStatus} (${orderOK})`);
    console.log(`   Export: ${expStatus}`);
    console.log();
});

console.log("=".repeat(80));
console.log("\nSUMMARY:");
const allPassed = modulesToCheck.every(module => {
    const classRegex = new RegExp(`class ${module.name}\\b`);
    let classLine = -1;
    for (let i = 0; i < lines.length; i++) {
        if (classRegex.test(lines[i])) {
            classLine = i;
            break;
        }
    }
    return classLine > 0 && classLine < exportsStart;
});

if (allPassed) {
    console.log("✅ ALL MODULES PROPERLY INTEGRATED (defined before module.exports)");
} else {
    console.log("❌ SOME MODULES ARE DEFINED AFTER module.exports - NEED RELOCATION");
}
