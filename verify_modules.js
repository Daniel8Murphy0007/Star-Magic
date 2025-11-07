// verify_modules.js - Comprehensive Module Verification Script
// Verifies source13-source125 modules for:
// 1. File existence (.js files)
// 2. Enhanced dynamics framework (25 methods)
// 3. Core functionality (imports, basic operations)

const fs = require('fs');
const path = require('path');

console.log('='.repeat(80));
console.log('STAR-MAGIC MODULE VERIFICATION SCRIPT');
console.log('Verifying source13-source125 integration and enhanced dynamics');
console.log('='.repeat(80));
console.log();

// Enhanced framework methods to check for
const enhancedMethods = [
    'createVariable', 'removeVariable', 'cloneVariable', 'listVariables', 'getSystemName',
    'transformVariableGroup', 'scaleVariableGroup',
    'expandParameterSpace',
    'autoRefineParameters', 'calibrateToObservations', 'optimizeForMetric',
    'generateVariations',
    'mutateParameters', 'evolveSystem',
    'saveState', 'restoreState', 'listSavedStates', 'exportState',
    'sensitivityAnalysis', 'generateReport', 'validateConsistency', 'autoCorrectAnomalies',
    'registerDynamicTerm', 'setDynamicParameter', 'getDynamicParameter'
];

const results = {
    filesExist: [],
    filesMissing: [],
    withEnhancedDynamics: [],
    withoutEnhancedDynamics: [],
    loadable: [],
    notLoadable: [],
    withClone: [],
    withoutClone: []
};

function checkFile(num) {
    const filename = `source${num}.js`;
    const filepath = path.join(__dirname, filename);
    
    if (!fs.existsSync(filepath)) {
        results.filesMissing.push(num);
        return null;
    }
    
    results.filesExist.push(num);
    
    try {
        const content = fs.readFileSync(filepath, 'utf8');
        
        // Check for enhanced dynamics
        const hasEnhanced = enhancedMethods.filter(method => 
            content.includes(method)
        ).length;
        
        if (hasEnhanced >= 20) {  // At least 20 of 25 methods
            results.withEnhancedDynamics.push(num);
        } else {
            results.withoutEnhancedDynamics.push(num);
        }
        
        // Check for clone method
        if (content.includes('clone()')) {
            results.withClone.push(num);
        } else {
            results.withoutClone.push(num);
        }
        
        // Try to load module
        try {
            const Module = require(filepath);
            if (typeof Module === 'function') {
                const instance = new Module();
                results.loadable.push(num);
                return { num, Module, instance, hasEnhanced };
            } else {
                results.loadable.push(num);
                return { num, Module, instance: null, hasEnhanced };
            }
        } catch (loadError) {
            results.notLoadable.push(num);
            return { num, error: loadError.message, hasEnhanced };
        }
        
    } catch (readError) {
        console.error(`Error reading ${filename}:`, readError.message);
        return null;
    }
}

console.log('Phase 1: File Existence Check\n');
for (let i = 13; i <= 125; i++) {
    checkFile(i);
}

console.log(`Files Found: ${results.filesExist.length}/113`);
console.log(`Files Missing: ${results.filesMissing.length}/113`);

if (results.filesMissing.length > 0 && results.filesMissing.length <= 20) {
    console.log('Missing files:', results.filesMissing.join(', '));
}

console.log();
console.log('Phase 2: Enhanced Dynamics Framework Analysis\n');
console.log(`With Enhanced Framework (≥20 methods): ${results.withEnhancedDynamics.length}`);
console.log(`Without Enhanced Framework: ${results.withoutEnhancedDynamics.length}`);

console.log();
console.log('Phase 3: Clone Method Support\n');
console.log(`With clone() method: ${results.withClone.length}`);
console.log(`Without clone() method: ${results.withoutClone.length}`);

console.log();
console.log('Phase 4: Module Loadability\n');
console.log(`Loadable: ${results.loadable.length}`);
console.log(`Not Loadable: ${results.notLoadable.length}`);

if (results.notLoadable.length > 0 && results.notLoadable.length <= 15) {
    console.log('Not loadable:', results.notLoadable.join(', '));
}

console.log();
console.log('='.repeat(80));
console.log('SUMMARY');
console.log('='.repeat(80));

const coverage = {
    files: (results.filesExist.length / 113 * 100).toFixed(1),
    enhanced: (results.withEnhancedDynamics.length / results.filesExist.length * 100).toFixed(1),
    clone: (results.withClone.length / results.filesExist.length * 100).toFixed(1),
    loadable: (results.loadable.length / results.filesExist.length * 100).toFixed(1)
};

console.log();
console.log(`File Coverage:              ${coverage.files}%  (${results.filesExist.length}/113 files exist)`);
console.log(`Enhanced Framework:         ${coverage.enhanced}%  (${results.withEnhancedDynamics.length}/${results.filesExist.length} modules)`);
console.log(`Clone Support:              ${coverage.clone}%  (${results.withClone.length}/${results.filesExist.length} modules)`);
console.log(`Loadable Modules:           ${coverage.loadable}%  (${results.loadable.length}/${results.filesExist.length} modules)`);

console.log();
console.log('='.repeat(80));
console.log('RECENT ADDITIONS');
console.log('='.repeat(80));

// Check recently modified files
const recentFiles = [];
for (const num of results.filesExist.slice(-5)) {
    const filename = `source${num}.js`;
    const filepath = path.join(__dirname, filename);
    const stats = fs.statSync(filepath);
    recentFiles.push({ num, mtime: stats.mtime });
}

recentFiles.sort((a, b) => b.mtime - a.mtime);
console.log('\nMost recently modified (top 5):');
recentFiles.slice(0, 5).forEach((file, i) => {
    const hasEnhanced = results.withEnhancedDynamics.includes(file.num) ? '✓ Enhanced' : '✗ Legacy';
    console.log(`  ${i + 1}. source${file.num}.js - ${file.mtime.toLocaleDateString()} ${file.mtime.toLocaleTimeString()} - ${hasEnhanced}`);
});

console.log();
console.log('='.repeat(80));
console.log('RECOMMENDATIONS');
console.log('='.repeat(80));
console.log();

if (results.filesMissing.length > 0) {
    console.log(`⚠ ${results.filesMissing.length} files missing - consider creating these modules`);
}

if (results.withoutEnhancedDynamics.length > 0) {
    console.log(`⚠ ${results.withoutEnhancedDynamics.length} modules lack enhanced dynamics framework`);
    console.log(`  → Priority: Upgrade source${results.withoutEnhancedDynamics[0]}-source${results.withoutEnhancedDynamics[Math.min(4, results.withoutEnhancedDynamics.length-1)]} first`);
}

if (results.withoutClone.length > 0) {
    console.log(`⚠ ${results.withoutClone.length} modules lack clone() method for parallel processing`);
}

if (results.notLoadable.length > 0) {
    console.log(`⚠ ${results.notLoadable.length} modules have loading errors - requires debugging`);
}

if (parseFloat(coverage.enhanced) === 100) {
    console.log(`✓ All modules have enhanced dynamics framework!`);
}

if (parseFloat(coverage.clone) === 100) {
    console.log(`✓ All modules support parallel processing with clone()!`);
}

console.log();
console.log('Verification complete.');
console.log();

// Export results for programmatic use
module.exports = results;
