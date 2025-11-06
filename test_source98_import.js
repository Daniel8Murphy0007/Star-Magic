// test_source98_import.js
// Test importing UnifiedFieldModule from index.js

const { UnifiedFieldModule } = require('./index.js');

console.log('Testing source98.js import through index.js...\n');

const module = new UnifiedFieldModule({ level: 13 });
const F_U = module.computeFU(0);

console.log(`✓ UnifiedFieldModule imported successfully`);
console.log(`✓ F_U = ${F_U.toExponential(4)} J/m³`);
console.log(`✓ source98.js integration complete\n`);
