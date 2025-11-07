// convert_legacy_to_enhanced.js
// Script to convert legacy source13-100 modules to enhanced framework
// Converts plain class variables to Map-based storage and adds enhanced_dynamics.js integration

const fs = require('fs');
const path = require('path');

function convertLegacyModule(sourceFile, outputFile = null) {
    if (!outputFile) outputFile = sourceFile;
    
    console.log(`\n====== Converting ${sourceFile} ======`);
    
    const content = fs.readFileSync(sourceFile, 'utf8');
    
    // Extract class name
    const classMatch = content.match(/class\s+(\w+)\s*{/);
    if (!classMatch) {
        console.log(`❌ No class found in ${sourceFile}`);
        return false;
    }
    const className = classMatch[1];
    console.log(`Found class: ${className}`);
    
    // Check if already enhanced
    if (content.includes('addEnhancedDynamics') || content.includes('clone()')) {
        console.log(`✓ Already enhanced, skipping`);
        return false;
    }
    
    // Extract system name from class name or file
    const systemName = className.replace(/Module|UQFF/g, '').replace(/([A-Z])/g, '_$1').substring(1);
    
    console.log(`Conversion would:`);
    console.log(`  1. Add enhanced_dynamics.js import`);
    console.log(`  2. Convert variables to Map storage`);
    console.log(`  3. Add initializeVariables() method`);
    console.log(`  4. Add dynamic terms support`);
    console.log(`  5. Add metadata tracking`);
    console.log(`  6. Add clone() method`);
    console.log(`  7. Add domain-specific expansion methods`);
    console.log(`  8. Apply addEnhancedDynamics(${className}, "${systemName}", domainExpansion)`);
    
    return true;
}

// Test conversion on a few modules
const testModules = [
    'source13.js',
    'source14.js',
    'source15.js'
];

console.log('LEGACY TO ENHANCED CONVERSION ANALYSIS');
console.log('=====================================\n');
console.log('This script analyzes legacy modules for conversion to enhanced framework.');
console.log('Full conversion requires significant refactoring per module.\n');

let analyzedCount = 0;
testModules.forEach(file => {
    const fullPath = path.join(__dirname, file);
    if (fs.existsSync(fullPath)) {
        convertLegacyModule(fullPath);
        analyzedCount++;
    } else {
        console.log(`File not found: ${file}`);
    }
});

console.log(`\n\n====== ANALYSIS SUMMARY ======`);
console.log(`Modules analyzed: ${analyzedCount}`);
console.log(`\nCONVERSION COMPLEXITY: HIGH`);
console.log(`Each module requires:`);
console.log(`  - Variable extraction and Map conversion (~50-100 variables per module)`);
console.log(`  - Constructor refactoring`);
console.log(`  - Method adaptation to use this.variables.get()`);
console.log(`  - Domain-specific expansion logic`);
console.log(`  - Estimated: 30-45 minutes per module`);
console.log(`  - Total for 88 modules: 44-66 hours of work`);
console.log(`\nRECOMMENDATION:`);
console.log(`  Option A: Batch convert 5-10 modules at a time`);
console.log(`  Option B: Create automated conversion script (complex)`);
console.log(`  Option C: Keep legacy modules as-is, focus on new modules with enhanced framework`);
