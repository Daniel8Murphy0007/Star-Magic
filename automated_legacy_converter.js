// automated_legacy_converter.js
// Automated conversion of legacy modules to enhanced framework
// Preserves all physics calculations and structure while adding enhanced dynamics

const fs = require('fs');
const path = require('path');

class LegacyConverter {
    constructor() {
        this.complexHelpers = `
// Complex number helpers
function complexAdd(a, b) { return {re: a.re + b.re, im: a.im + b.im}; }
function complexSub(a, b) { return {re: a.re - b.re, im: a.im - b.im}; }
function complexMul(a, b) { return {re: a.re*b.re - a.im*b.im, im: a.re*b.im + a.im*b.re}; }
function complexDiv(a, b) {
  const denom = b.re*b.re + b.im*b.im;
  return {re: (a.re*b.re + a.im*b.im)/denom, im: (a.im*b.re - a.re*b.im)/denom};
}
function complexPow(base, exponent) {
  const r = Math.sqrt(base.re*base.re + base.im*base.im);
  const theta = Math.atan2(base.im, base.re);
  const newR = Math.pow(r, exponent);
  const newTheta = theta * exponent;
  return {re: newR * Math.cos(newTheta), im: newR * Math.sin(newTheta)};
}
function complexScale(a, s) { return {re: a.re*s, im: a.im*s}; }
function complexNeg(a) { return {re: -a.re, im: -a.im}; }
function complexAbs(a) { return Math.sqrt(a.re*a.re + a.im*a.im); }
function toComplex(x) { return typeof x === 'object' ? x : {re: x, im: 0}; }
`;
    }

    extractVariables(content) {
        const variables = new Set();
        const lines = content.split('\n');
        let inConstructor = false;
        let inInitDefaults = false;
        
        for (const line of lines) {
            if (line.includes('constructor()') || line.includes('initializeDefaults()')) {
                inConstructor = true;
                inInitDefaults = true;
            }
            if (inConstructor && line.trim() === '}') {
                inConstructor = false;
                inInitDefaults = false;
            }
            
            if (inInitDefaults) {
                const match = line.match(/this\.(\w+)\s*=/);
                if (match && match[1] !== 'updateCache') {
                    variables.add(match[1]);
                }
            }
        }
        
        return Array.from(variables);
    }

    extractClassName(content) {
        const match = content.match(/class\s+(\w+)\s*{/);
        return match ? match[1] : null;
    }

    extractSystemName(className, content) {
        // Try to extract from comments
        const commentMatch = content.match(/\/\/.*?-\s*(.+?)\s+Module/i);
        if (commentMatch) return commentMatch[1].replace(/[^a-zA-Z0-9]/g, '_');
        
        // Fallback to class name
        return className.replace(/Module|UQFF/g, '').replace(/([A-Z])/g, '_$1').substring(1);
    }

    preservePhysicsMethods(content) {
        // Extract all compute methods to preserve
        const methods = [];
        const methodRegex = /(\s+\w+\([^)]*\)\s*{[\s\S]*?^\s+})/gm;
        let match;
        
        while ((match = methodRegex.exec(content)) !== null) {
            const method = match[1];
            if (method.includes('compute') || method.includes('get') || method.includes('set') || 
                method.includes('add') || method.includes('subtract') || method.includes('print') ||
                method.includes('example')) {
                methods.push(method);
            }
        }
        
        return methods;
    }

    convertToEnhanced(sourceFile, dryRun = false) {
        console.log(`\n${'='.repeat(60)}`);
        console.log(`Converting: ${path.basename(sourceFile)}`);
        console.log('='.repeat(60));
        
        if (!fs.existsSync(sourceFile)) {
            console.log(`❌ File not found`);
            return false;
        }

        let content = fs.readFileSync(sourceFile, 'utf8');
        
        // Check if already enhanced
        if (content.includes('addEnhancedDynamics') || content.includes('const { addEnhancedDynamics }')) {
            console.log(`✓ Already enhanced, skipping`);
            return false;
        }

        const className = this.extractClassName(content);
        if (!className) {
            console.log(`❌ No class found`);
            return false;
        }

        const variables = this.extractVariables(content);
        const systemName = this.extractSystemName(className, content);
        
        console.log(`  Class: ${className}`);
        console.log(`  System: ${systemName}`);
        console.log(`  Variables: ${variables.length} extracted`);

        // SIMPLER APPROACH: Inject enhanced dynamics into existing class
        // Find the class closing brace
        const classMatch = content.match(/(class\s+\w+\s*{[\s\S]*)(}\s*module\.exports)/);
        if (!classMatch) {
            console.log(`❌ Could not parse class structure`);
            return false;
        }

        // Insert enhanced infrastructure before class constructor ends
        const constructorEndMatch = content.match(/(constructor\(\)\s*{[\s\S]*?)(    }\s)/);
        if (constructorEndMatch) {
            const enhancedInit = `\n    // Enhanced dynamics infrastructure
    this.variables = new Map();
    this.dynamicTerms = [];
    this.dynamicParameters = new Map();
    this.metadata = new Map();
    this.metadata.set("enhanced", true);
    this.metadata.set("version", "2.0.0");
    this.metadata.set("system_name", "${systemName}");
    this.enableLogging = false;
    this.learningRate = 0.01;
    this.initializeVariablesMap();
`;
            content = content.replace(constructorEndMatch[2], enhancedInit + constructorEndMatch[2]);
        }

        // Add initializeVariablesMap method after constructor
        const initMethod = `
    initializeVariablesMap() {
        const vars = [${variables.map(v => `"${v}"`).join(', ')}];
        vars.forEach(v => {
            if (this[v] !== undefined) {
                this.variables.set(v, toComplex(this[v]));
            }
        });
    }
`;
        content = content.replace(/(    }\s+)(    \w+\(|$)/, `$1${initMethod}\n$2`);

        // Add enhanced methods before class closing brace
        const enhancedMethods = `
    setEnableLogging(enable) { this.enableLogging = enable; }
    registerDynamicTerm(term) { this.dynamicTerms.push(term); }
    setDynamicParameter(name, value) { this.dynamicParameters.set(name, value); }
    getDynamicParameter(name) { return this.dynamicParameters.get(name); }

    clone() {
        const cloned = new ${className}();
        cloned.variables = new Map(this.variables);
        cloned.dynamicParameters = new Map(this.dynamicParameters);
        cloned.metadata = new Map(this.metadata);
        cloned.enableLogging = this.enableLogging;
        cloned.learningRate = this.learningRate;
        return cloned;
    }
`;
        content = content.replace(/(}\s*module\.exports)/, `${enhancedMethods}}\n\nconst domainExpansion = {
    expandSystemScale(massFactor, radiusFactor) {
        if (this.variables.has("M")) this.variables.get("M").re *= massFactor;
        if (this.variables.has("r")) this.variables.get("r").re *= radiusFactor;
        if (this.enableLogging) console.log(\`Expanded: M×\${massFactor}, r×\${radiusFactor}\`);
    },
    expandPhysicsScale(factor) {
        const keys = ["B0", "B", "L0_W", "rho_fluid", "A_osc"];
        keys.forEach(k => {
            if (this.variables.has(k)) this.variables.get(k).re *= factor;
        });
        if (this.enableLogging) console.log(\`Expanded physics: ×\${factor}\`);
    }
};

addEnhancedDynamics(${className}, "${systemName}", domainExpansion);\n\n$1`);

        // Add imports at top
        const headerEnd = content.search(/class\s+\w+/);
        const header = content.substring(0, headerEnd);
        const body = content.substring(headerEnd);
        content = header + `const { addEnhancedDynamics } = require('./enhanced_dynamics.js');\n${this.complexHelpers}\n` + body;

        if (dryRun) {
            console.log(`  ✓ Conversion prepared (dry run)`);
            return true;
        }

        // Backup original
        const backupFile = sourceFile.replace('.js', '.legacy.js');
        if (!fs.existsSync(backupFile)) {
            fs.copyFileSync(sourceFile, backupFile);
            console.log(`  📦 Backup: ${path.basename(backupFile)}`);
        }

        // Write enhanced version
        fs.writeFileSync(sourceFile, content);
        console.log(`  ✅ Enhanced version written`);
        console.log(`  📊 Variables: ${variables.length}, Original structure preserved`);

        return true;
    }

    batchConvert(startNum, endNum, dryRun = false) {
        console.log(`\n${'='.repeat(60)}`);
        console.log(`BATCH CONVERSION: source${startNum}-source${endNum}`);
        console.log(`Mode: ${dryRun ? 'DRY RUN' : 'LIVE CONVERSION'}`);
        console.log('='.repeat(60));

        let converted = 0;
        let skipped = 0;
        let errors = 0;

        for (let i = startNum; i <= endNum; i++) {
            const fileName = `source${i}.js`;
            const filePath = path.join(__dirname, fileName);
            
            if (!fs.existsSync(filePath)) {
                console.log(`\n⊘ source${i}.js - Not found, skipping`);
                skipped++;
                continue;
            }

            try {
                const result = this.convertToEnhanced(filePath, dryRun);
                if (result) {
                    converted++;
                } else {
                    skipped++;
                }
            } catch (error) {
                console.log(`\n❌ source${i}.js - Error: ${error.message}`);
                errors++;
            }
        }

        console.log(`\n${'='.repeat(60)}`);
        console.log(`BATCH CONVERSION SUMMARY`);
        console.log('='.repeat(60));
        console.log(`  Converted: ${converted}`);
        console.log(`  Skipped: ${skipped}`);
        console.log(`  Errors: ${errors}`);
        console.log(`  Total processed: ${converted + skipped + errors}`);
        console.log('='.repeat(60));

        return { converted, skipped, errors };
    }
}

// Execute conversion
const converter = new LegacyConverter();

// Command line arguments
const args = process.argv.slice(2);
const dryRun = args.includes('--dry-run');
const testMode = args.includes('--test');

if (testMode) {
    console.log('\n🧪 TEST MODE: Converting source13-15 only\n');
    converter.batchConvert(13, 15, dryRun);
} else {
    console.log('\n🚀 FULL CONVERSION: source13-100\n');
    const readline = require('readline').createInterface({
        input: process.stdin,
        output: process.stdout
    });
    
    readline.question('This will convert 88 modules. Continue? (yes/no): ', (answer) => {
        if (answer.toLowerCase() === 'yes') {
            converter.batchConvert(13, 100, dryRun);
        } else {
            console.log('Conversion cancelled.');
        }
        readline.close();
    });
}

module.exports = LegacyConverter;
