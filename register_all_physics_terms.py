#!/usr/bin/env python3
"""
Extract all PhysicsTerm class names from MAIN_1_CoAnQi.cpp and generate registration code.
"""

import re

def main():
    # Read MAIN_1_CoAnQi.cpp
    with open("MAIN_1_CoAnQi.cpp", "r", encoding="utf-8") as f:
        content = f.read()
    
    # Find all PhysicsTerm class declarations
    pattern = r'^class (\w+Term) : public PhysicsTerm'
    matches = re.findall(pattern, content, re.MULTILINE)
    
    print(f"Found {len(matches)} PhysicsTerm classes")
    
    # Remove duplicates (keep order)
    seen = set()
    unique_classes = []
    for cls in matches:
        if cls not in seen:
            seen.add(cls)
            unique_classes.append(cls)
    
    print(f"Unique classes: {len(unique_classes)}")
    
    # Generate registration function
    output = []
    output.append("// AUTO-GENERATED: Register all physics terms")
    output.append("void registerAllPhysicsTerms(CalculatorCore& core) {")
    
    for cls in unique_classes:
        # Create instance and register
        output.append(f"    core.registerPhysicsTerm(\"{cls}\", std::make_unique<{cls}>(), \"auto-registered\");")
    
    output.append("}")
    output.append("")
    output.append(f"// Total registered: {len(unique_classes)} physics terms")
    
    # Write to file
    with open("physics_registration.cpp", "w", encoding="utf-8") as f:
        f.write("\n".join(output))
    
    print(f"Generated physics_registration.cpp with {len(unique_classes)} registrations")
    print(f"\nInsert this before main() function in MAIN_1_CoAnQi.cpp")
    print(f"Then call: registerAllPhysicsTerms(g_calculatorCore);")

if __name__ == "__main__":
    main()
