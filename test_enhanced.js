// test_enhanced.js - Quick test of enhanced modules
const files = [131, 132, 133, 134, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162];

console.log('Testing Enhanced Framework Modules\n');
console.log('File      | Clone | Save  | Report | Status');
console.log('----------|-------|-------|--------|--------');

files.forEach(n => {
    try {
        const M = require(`./source${n}.js`);
        const m = new M();
        const clone = typeof m.clone === 'function' ? '✓' : '✗';
        const save = typeof m.saveState === 'function' ? '✓' : '✗';
        const report = typeof m.generateReport === 'function' ? '✓' : '✗';
        const status = (clone === '✓' && save === '✓' && report === '✓') ? 'FULL' : 'PARTIAL';
        console.log(`source${n.toString().padStart(3)} |   ${clone}   |   ${save}   |   ${report}    | ${status}`);
    } catch (e) {
        console.log(`source${n.toString().padStart(3)} | ERROR: ${e.message.substring(0, 40)}`);
    }
});

console.log();
