#!/usr/bin/env python3
"""
Session 152 build fix script.
Fixes all known compilation errors in wolfram_sources_bridge.cpp and MAIN_1_CoAnQi.cpp.
"""
import re

print("=== Fix Script Start ===")

# ============================================================
# 1. Fix wolfram_sources_bridge.cpp
# ============================================================
print("\n[1] Fixing wolfram_sources_bridge.cpp...")
with open("wolfram_sources_bridge.cpp", "r", encoding="utf-8-sig") as f:
    wsb = f.read()

orig_size = len(wsb)

# 1a. Add file-scope constant c (speed of light) right after #include <complex>
# Only if not already present
if "// speed of light (file-scope for bridge classes)" not in wsb:
    wsb = wsb.replace(
        "#include <complex>",
        '#include <complex>\nstatic constexpr double c   = 2.998e8;   // speed of light (file-scope for bridge classes)\nstatic constexpr double H0_bridge  = 2.27e-18;  // Hubble constant (file-scope)'
    )
    print("  Added file-scope c and H0_bridge constants")

# 1b. Fix static constexpr double PI -> PI_const (conflicts with #define PI in MAIN_1)
wsb = wsb.replace("static constexpr double PI = 3.14159265358979323846;", "static constexpr double PI_const = 3.14159265358979323846;")
wsb = wsb.replace("static constexpr double PI = 3.141592653589793;", "static constexpr double PI_const = 3.141592653589793;")
print(f"  PI_const occurrences: {wsb.count('PI_const')}")

# 1c. Fix P member variable in SGR1745MagnetarSpinEMTerm (WSTP P(s) macro conflict)
wsb = wsb.replace("    double r, P;\r\npublic:", "    double r, spin_period;\r\npublic:")
wsb = wsb.replace("    double r, P;\npublic:", "    double r, spin_period;\npublic:")
# Fix initializer and usage
wsb = wsb.replace("r(1e4), P(3.76) {\r\n        v_spin = (2 * 3.141592653589793 * r) / P;",
                  "r(1e4), spin_period(3.76) {\r\n        v_spin = (2 * 3.141592653589793 * r) / spin_period;")
wsb = wsb.replace("r(1e4), P(3.76) {\n        v_spin = (2 * 3.141592653589793 * r) / P;",
                  "r(1e4), spin_period(3.76) {\n        v_spin = (2 * 3.141592653589793 * r) / spin_period;")

# 1d. Fix P member variable in SGR1745OscillatoryWaveTerm
wsb = wsb.replace("    double A, k, omega, x, pi, P;", "    double A, k, omega, x, pi, spin_period;")
wsb = wsb.replace("pi(3.141592653589793), P(3.76) {\r\n        omega = 2 * pi / P;",
                  "pi(3.141592653589793), spin_period(3.76) {\r\n        omega = 2 * pi / spin_period;")
wsb = wsb.replace("pi(3.141592653589793), P(3.76) {\n        omega = 2 * pi / spin_period;",
                  "pi(3.141592653589793), spin_period(3.76) {\n        omega = 2 * pi / spin_period;")
wsb = wsb.replace("pi(3.141592653589793), P(3.76) {\n        omega = 2 * pi / P;",
                  "pi(3.141592653589793), spin_period(3.76) {\n        omega = 2 * pi / spin_period;")

print(f"  spin_period occurrences: {wsb.count('spin_period')}")

# 1e. Fix compute(double t_time) -> correct 2-arg signature
wsb = wsb.replace(
    "double compute(double t_time) const override {",
    "double compute(double t, const std::map<std::string, double>& /*params*/) const override {"
)
print(f"  compute(t_time) remaining: {wsb.count('compute(double t_time)')}")

# 1f. Fix H0 in default arg -> literal 
wsb = wsb.replace("double Hz = H0)", "double Hz = 2.27e-18)")
print(f"  H0 default arg remaining: {wsb.count('= H0)')}")

with open("wolfram_sources_bridge.cpp", "w", encoding="utf-8") as f:
    f.write(wsb)
print(f"  Written wolfram_sources_bridge.cpp ({orig_size} -> {len(wsb)} chars)")

# ============================================================
# 2. Fix source177_wolfram_field_unity.cpp
# ============================================================
print("\n[2] Fixing source177_wolfram_field_unity.cpp...")
with open("source177_wolfram_field_unity.cpp", "r", encoding="utf-8-sig") as f:
    s177 = f.read()

# Fix complex<double> -> double conversion
s177 = s177.replace(
    "double magnetic_field = decoder.getMagneticField(0, 0.0);",
    "double magnetic_field = decoder.getMagneticField(0, 0.0).real();"
)
print(f"  .real() applied: {'.real()' in s177}")

with open("source177_wolfram_field_unity.cpp", "w", encoding="utf-8") as f:
    f.write(s177)
print("  Written source177_wolfram_field_unity.cpp")

# ============================================================
# 3. Fix MAIN_1_CoAnQi.cpp
# ============================================================
print("\n[3] Fixing MAIN_1_CoAnQi.cpp (large file, please wait)...")
with open("MAIN_1_CoAnQi.cpp", "r", encoding="utf-8-sig") as f:
    main = f.read()

print(f"  File size: {len(main)} chars")

# 3a. Fix validate() override mismatch (3 classes)
# Base class: virtual bool validate(const std::map<std::string, double>& params) const
# Derived: bool validate() const override { return true; }
main = main.replace(
    "    bool validate() const override { return true; }\n};\n\nclass PIInfinityTerm",
    "    bool validate(const std::map<std::string, double>& /*params*/) const override { return true; }\n};\n\nclass PIInfinityTerm"
)
main = main.replace(
    "    bool validate() const override { return true; }\n};\n\nclass ConsciousnessTerm",
    "    bool validate(const std::map<std::string, double>& /*params*/) const override { return true; }\n};\n\nclass ConsciousnessTerm"
)
main = main.replace(
    "    bool validate() const override { return true; }\n};\n\n// ===== END SESSION 129",
    "    bool validate(const std::map<std::string, double>& /*params*/) const override { return true; }\n};\n\n// ===== END SESSION 129"
)
# Also handle CRLF versions
main = main.replace(
    "    bool validate() const override { return true; }\r\n};\r\n\r\nclass PIInfinityTerm",
    "    bool validate(const std::map<std::string, double>& /*params*/) const override { return true; }\r\n};\r\n\r\nclass PIInfinityTerm"
)
main = main.replace(
    "    bool validate() const override { return true; }\r\n};\r\n\r\nclass ConsciousnessTerm",
    "    bool validate(const std::map<std::string, double>& /*params*/) const override { return true; }\r\n};\r\n\r\nclass ConsciousnessTerm"
)
main = main.replace(
    "    bool validate() const override { return true; }\r\n};\r\n\r\n// ===== END SESSION 129",
    "    bool validate(const std::map<std::string, double>& /*params*/) const override { return true; }\r\n};\r\n\r\n// ===== END SESSION 129"
)

print(f"  validate() const override remaining: {main.count('bool validate() const override')}")

# 3b. Fix SOURCE4 namespace misplaced inside function body
# Find the anchor points using unique strings
registration_func_start = "void registerAllPhysicsTerms(CalculatorCore &core)"
source4_comment = "// SOURCE4 PhysicsTerm class definitions (moved before BATCH 20 registration)"
source4_ns_start = "namespace SOURCE4 {"
source4_ns_end = "} // namespace SOURCE4"

# Find positions
pos_func = main.find(registration_func_start)
pos_s4_comment = main.find(source4_comment, pos_func)
pos_s4_ns_end = main.find(source4_ns_end, pos_s4_comment)

if pos_func == -1 or pos_s4_comment == -1 or pos_s4_ns_end == -1:
    print(f"  ERROR: Could not find anchor points: func={pos_func}, s4={pos_s4_comment}, s4end={pos_s4_ns_end}")
else:
    print(f"  Found: func at {pos_func}, SOURCE4 comment at {pos_s4_comment}, SOURCE4 end at {pos_s4_ns_end}")
    
    # Extract the SOURCE4 block (from the comment up to end of namespace incl. newline)
    # Find the exact end: after "} // namespace SOURCE4\n"
    s4_block_end = pos_s4_ns_end + len(source4_ns_end)
    # Include trailing newlines
    while s4_block_end < len(main) and main[s4_block_end] in '\r\n':
        s4_block_end += 1
    
    source4_block = main[pos_s4_comment:s4_block_end]
    print(f"  SOURCE4 block: {len(source4_block)} chars, starts with: {source4_block[:50]!r}")
    
    # Remove SOURCE4 block from current location
    main_without_s4 = main[:pos_s4_comment] + main[s4_block_end:]
    
    # Find insertion point: just before "void registerAllPhysicsTerms"
    # Find the comment block just before the function
    insert_target = "void registerAllPhysicsTerms(CalculatorCore &core)"
    insert_pos = main_without_s4.find(insert_target)
    
    if insert_pos == -1:
        print("  ERROR: Could not find insert target after removing block")
    else:
        # Insert SOURCE4 block (with separator) before the function
        main_fixed = (main_without_s4[:insert_pos] + 
                      source4_block + "\n" +
                      main_without_s4[insert_pos:])
        print(f"  New file size: {len(main_fixed)} chars")
        
        # Verify: source4 namespace is no longer inside the function
        pos_func_new = main_fixed.find(registration_func_start)
        pos_s4_new = main_fixed.find(source4_comment)
        print(f"  Verification: SOURCE4 at {pos_s4_new}, function at {pos_func_new}")
        if pos_s4_new < pos_func_new:
            print("  OK: SOURCE4 is now BEFORE the function")
            main = main_fixed
        else:
            print("  ERROR: Ordering wrong!")

with open("MAIN_1_CoAnQi.cpp", "w", encoding="utf-8") as f:
    f.write(main)
print(f"  Written MAIN_1_CoAnQi.cpp ({len(main)} chars)")

print("\n=== Fix Script Complete ===")
