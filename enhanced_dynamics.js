// enhanced_dynamics.js
// Reusable 25-method self-expansion framework for UQFF modules
// Implements full dynamic self-update capabilities demonstrated by source147.cpp
// Copyright - Daniel T. Murphy, Nov 2025

/**
 * Adds 25 enhanced dynamic methods to a UQFF module class
 * @param {Object} moduleClass - The module class to enhance
 * @param {string} systemName - Unique system identifier
 * @param {Function} domainExpansion - Domain-specific expansion logic
 */
function addEnhancedDynamics(moduleClass, systemName, domainExpansion) {
  const savedStates = new Map();

  // ========== 1. Variable Management (5 methods) ==========
  
  moduleClass.prototype.createVariable = function(name, value) {
    this.variables.set(name, value);
    if (this.enableLogging) console.log(`Created variable: ${name} = ${value.re} + i ${value.im}`);
  };

  moduleClass.prototype.removeVariable = function(name) {
    this.variables.delete(name);
    if (this.enableLogging) console.log(`Removed variable: ${name}`);
  };

  moduleClass.prototype.cloneVariable = function(source, destination) {
    if (this.variables.has(source)) {
      this.variables.set(destination, {...this.variables.get(source)});
      if (this.enableLogging) console.log(`Cloned ${source} → ${destination}`);
    }
  };

  moduleClass.prototype.listVariables = function() {
    return Array.from(this.variables.keys());
  };

  moduleClass.prototype.getSystemName = function() {
    return systemName;
  };

  // ========== 2. Batch Operations (2 methods) ==========
  
  moduleClass.prototype.transformVariableGroup = function(names, transformFunc) {
    names.forEach(name => {
      if (this.variables.has(name)) {
        this.variables.set(name, transformFunc(this.variables.get(name)));
      }
    });
  };

  moduleClass.prototype.scaleVariableGroup = function(names, scaleFactor) {
    this.transformVariableGroup(names, val => ({
      re: val.re * scaleFactor.re - val.im * scaleFactor.im,
      im: val.re * scaleFactor.im + val.im * scaleFactor.re
    }));
  };

  // ========== 3. Self-Expansion (4 methods) ==========
  
  moduleClass.prototype.expandParameterSpace = function(globalScale) {
    for (const [key, value] of this.variables) {
      this.variables.set(key, {re: value.re * globalScale, im: value.im * globalScale});
    }
    if (this.enableLogging) console.log(`Expanded parameter space by factor ${globalScale}`);
  };

  // Domain-specific expansion methods injected by caller
  if (domainExpansion) {
    Object.assign(moduleClass.prototype, domainExpansion);
  }

  // ========== 4. Self-Refinement (3 methods) ==========
  
  moduleClass.prototype.autoRefineParameters = function(targetMetric) {
    // Monte Carlo parameter refinement
    const paramNames = ['M', 'r', 'DPM_momentum', 'k_LENR'];
    const ranges = {
      M: [1e30, 1e45],
      r: [1e15, 1e23],
      DPM_momentum: [0.01, 10.0],
      k_LENR: [1e-12, 1e-8]
    };
    
    paramNames.forEach(name => {
      if (this.variables.has(name) && ranges[name]) {
        const [min, max] = ranges[name];
        const newValue = min + Math.random() * (max - min);
        this.variables.get(name).re = newValue;
      }
    });
    
    if (this.enableLogging) console.log(`Auto-refined for metric: ${targetMetric}`);
  };

  moduleClass.prototype.calibrateToObservations = function(observedValues) {
    Object.entries(observedValues).forEach(([key, value]) => {
      if (this.variables.has(key)) {
        this.variables.set(key, value);
      }
    });
    if (this.enableLogging) console.log(`Calibrated to ${Object.keys(observedValues).length} observations`);
  };

  moduleClass.prototype.optimizeForMetric = function(metricName) {
    // Subclasses should override with domain-specific optimization presets
    if (this.enableLogging) console.log(`Optimizing for metric: ${metricName}`);
  };

  // ========== 5. Parameter Exploration (1 method) ==========
  
  moduleClass.prototype.generateVariations = function(count, variationPercent) {
    const variations = [];
    const delta = variationPercent / 100.0;
    
    for (let i = 0; i < count; i++) {
      const variant = new Map();
      for (const [key, value] of this.variables) {
        const deltaRe = (Math.random() * 2 - 1) * delta * value.re;
        const deltaIm = (Math.random() * 2 - 1) * delta * value.im;
        variant.set(key, {re: value.re + deltaRe, im: value.im + deltaIm});
      }
      variations.push(variant);
    }
    
    return variations;
  };

  // ========== 6. Adaptive Evolution (2 methods) ==========
  
  moduleClass.prototype.mutateParameters = function(mutationRate) {
    for (const [key, value] of this.variables) {
      const deltaRe = (Math.random() * 2 - 1) * mutationRate * value.re;
      const deltaIm = (Math.random() * 2 - 1) * mutationRate * value.im;
      this.variables.set(key, {re: value.re + deltaRe, im: value.im + deltaIm});
    }
    if (this.enableLogging) console.log(`Mutated parameters with rate ${mutationRate}`);
  };

  moduleClass.prototype.evolveSystem = function(generations, fitnessFunc) {
    let bestFitness = fitnessFunc(this);
    const bestState = new Map(this.variables);
    
    for (let gen = 0; gen < generations; gen++) {
      this.mutateParameters(0.05);
      const currentFitness = fitnessFunc(this);
      
      if (currentFitness > bestFitness) {
        bestFitness = currentFitness;
        bestState.clear();
        for (const [k, v] of this.variables) bestState.set(k, v);
      } else {
        this.variables = new Map(bestState);
      }
    }
    
    if (this.enableLogging) console.log(`Evolved ${generations} generations, best fitness: ${bestFitness}`);
  };

  // ========== 7. State Management (4 methods) ==========
  
  moduleClass.prototype.saveState = function(stateName) {
    savedStates.set(stateName, new Map(this.variables));
    if (this.enableLogging) console.log(`Saved state: ${stateName}`);
  };

  moduleClass.prototype.restoreState = function(stateName) {
    if (savedStates.has(stateName)) {
      this.variables = new Map(savedStates.get(stateName));
      if (this.enableLogging) console.log(`Restored state: ${stateName}`);
    }
  };

  moduleClass.prototype.listSavedStates = function() {
    return Array.from(savedStates.keys());
  };

  moduleClass.prototype.exportState = function() {
    let output = `System: ${systemName}\n`;
    for (const [key, value] of this.variables) {
      output += `${key} = ${value.re.toExponential()} + i*${value.im.toExponential()}\n`;
    }
    return output;
  };

  // ========== 8. System Analysis (4 methods) ==========
  
  moduleClass.prototype.sensitivityAnalysis = function(paramNames, deltaPercent) {
    const results = {};
    const baselineF = this.computeF ? this.computeF(this.variables.get("t").re) : {re: 0, im: 0};
    
    paramNames.forEach(name => {
      if (this.variables.has(name)) {
        const original = {...this.variables.get(name)};
        const delta = original.re * (deltaPercent / 100.0);
        
        this.variables.get(name).re += delta;
        const upperF = this.computeF ? this.computeF(this.variables.get("t").re) : {re: 0, im: 0};
        
        this.variables.set(name, original);
        this.variables.get(name).re -= delta;
        const lowerF = this.computeF ? this.computeF(this.variables.get("t").re) : {re: 0, im: 0};
        
        this.variables.set(name, original);
        
        results[name] = {
          sensitivity: ((upperF.re - lowerF.re) / (2 * delta)),
          baseline: baselineF.re,
          upper: upperF.re,
          lower: lowerF.re
        };
      }
    });
    
    return results;
  };

  moduleClass.prototype.generateReport = function() {
    let report = `=== ${systemName} UQFF Report ===\n`;
    report += `Variables: ${this.variables.size}\n`;
    report += `Dynamic Terms: ${this.dynamicTerms.length}\n`;
    report += `Dynamic Parameters: ${this.dynamicParameters.size}\n`;
    report += `Enhanced: ${this.metadata.get("enhanced")}\n`;
    report += `Version: ${this.metadata.get("version")}\n`;
    
    if (this.computeF) {
      const F = this.computeF(this.variables.get("t").re);
      report += `F(t) = ${F.re.toExponential()} + i*${F.im.toExponential()} N\n`;
    }
    
    return report;
  };

  moduleClass.prototype.validateConsistency = function() {
    // Check for NaN, Infinity, or physically impossible values
    for (const [key, value] of this.variables) {
      if (!isFinite(value.re) || !isFinite(value.im)) {
        console.error(`Invalid value for ${key}: ${value.re} + i ${value.im}`);
        return false;
      }
      
      // Physical constraints
      if (key === "M" && value.re <= 0) return false;
      if (key === "r" && value.re <= 0) return false;
      if (key === "t" && value.re < 0) return false;
    }
    
    return true;
  };

  moduleClass.prototype.autoCorrectAnomalies = function() {
    for (const [key, value] of this.variables) {
      // Replace NaN/Infinity with safe defaults
      if (!isFinite(value.re)) value.re = 0;
      if (!isFinite(value.im)) value.im = 0;
      
      // Enforce positivity for physical quantities
      if ((key === "M" || key === "r" || key === "L_X") && value.re < 0) {
        value.re = Math.abs(value.re);
      }
    }
    
    if (this.enableLogging) console.log("Auto-corrected anomalies");
  };
}

module.exports = { addEnhancedDynamics };
