#include "../include/dimensional_analysis.h"
#include <stdexcept>
#include <cmath>

PhysicalQuantity::PhysicalQuantity(double val, const DimensionVector& dims)
    : value(val), dimensions(dims)
{
}

PhysicalQuantity PhysicalQuantity::operator+(const PhysicalQuantity& other) const {
    if (dimensions != other.dimensions) {
        throw std::runtime_error("Incompatible dimensions for addition");
    }
    return PhysicalQuantity(value + other.value, dimensions);
}

PhysicalQuantity PhysicalQuantity::operator*(const PhysicalQuantity& other) const {
    DimensionVector newDims = dimensions;
    for (size_t i = 0; i < 7; i++) {
        newDims[i] += other.dimensions[i];
    }
    return PhysicalQuantity(value * other.value, newDims);
}

PhysicalQuantity PhysicalQuantity::operator/(const PhysicalQuantity& other) const {
    DimensionVector newDims = dimensions;
    for (size_t i = 0; i < 7; i++) {
        newDims[i] -= other.dimensions[i];
    }
    return PhysicalQuantity(value / other.value, newDims);
}

bool PhysicalQuantity::isCompatibleWith(const PhysicalQuantity& other) const {
    return dimensions == other.dimensions;
}

DimensionalSystem::DimensionalSystem() {
    initializeBaseQuantities();
}

void DimensionalSystem::initializeBaseQuantities() {
    // M=0, L=1, T=2, I=3, Θ=4, N=5, J=6
    baseQuantities_["length"] = PhysicalQuantity(1.0, {0,1,0,0,0,0,0});
    baseQuantities_["mass"] = PhysicalQuantity(1.0, {1,0,0,0,0,0,0});
    baseQuantities_["time"] = PhysicalQuantity(1.0, {0,0,1,0,0,0,0});
    baseQuantities_["current"] = PhysicalQuantity(1.0, {0,0,0,1,0,0,0});
    baseQuantities_["temperature"] = PhysicalQuantity(1.0, {0,0,0,0,1,0,0});
    baseQuantities_["amount"] = PhysicalQuantity(1.0, {0,0,0,0,0,1,0});
    baseQuantities_["luminosity"] = PhysicalQuantity(1.0, {0,0,0,0,0,0,1});
    
    // Derived UQFF quantities
    baseQuantities_["force"] = PhysicalQuantity(1.0, {1,1,-2,0,0,0,0}); // kg·m/s²
    baseQuantities_["energy"] = PhysicalQuantity(1.0, {1,2,-2,0,0,0,0}); // kg·m²/s²
    baseQuantities_["power"] = PhysicalQuantity(1.0, {1,2,-3,0,0,0,0}); // kg·m²/s³
    baseQuantities_["frequency"] = PhysicalQuantity(1.0, {0,0,-1,0,0,0,0}); // 1/s
}

PhysicalQuantity DimensionalSystem::getBaseQuantity(const std::string& name) const {
    auto it = baseQuantities_.find(name);
    if (it == baseQuantities_.end()) {
        throw std::runtime_error("Unknown base quantity: " + name);
    }
    return it->second;
}

bool DimensionalSystem::validateDimensions(const PhysicalQuantity& quantity) const {
    // Check if all dimension exponents are integers (they are by design)
    // Could add more validation here
    return true;
}
