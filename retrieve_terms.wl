(* Wolfram Language Script: Retrieve comprehensive physics terms database *)
output = OpenWrite["wolfram_physics_terms_5034.txt"];

(* 1. Physical Quantities *)
WriteString[output, "=== PHYSICAL QUANTITIES ===\n"];
WriteString[output, ToString[EntityList["PhysicalQuantity"], InputForm] <> "\n\n"];

(* 2. Physical Constants *)
WriteString[output, "=== PHYSICAL CONSTANTS ===\n"];
WriteString[output, ToString[EntityList["PhysicalConstant"], InputForm] <> "\n\n"];

(* 3. Particles *)
WriteString[output, "=== PARTICLES ===\n"];
WriteString[output, ToString[EntityList["Particle"], InputForm] <> "\n\n"];

(* 4. Isotopes *)
WriteString[output, "=== ISOTOPES ===\n"];
WriteString[output, ToString[Take[EntityList["Isotope"], UpTo[500]], InputForm] <> "\n\n"];

(* 5. Stars *)
WriteString[output, "=== STARS ===\n"];
WriteString[output, ToString[Take[EntityList["Star"], UpTo[500]], InputForm] <> "\n\n"];

(* 6. Galaxies *)
WriteString[output, "=== GALAXIES ===\n"];
WriteString[output, ToString[Take[EntityList["Galaxy"], UpTo[500]], InputForm] <> "\n\n"];

(* 7. Mathematical Physics Functions *)
WriteString[output, "=== SPECIAL FUNCTIONS ===\n"];
funcs = Join[
  Names["System`*Bessel*"],
  Names["System`*Spherical*"],
  Names["System`*Legendre*"],
  Names["System`*Hermite*"],
  Names["System`*Laguerre*"],
  Names["System`*Jacobi*"],
  Names["System`*Hypergeometric*"]
];
WriteString[output, ToString[funcs, InputForm] <> "\n\n"];

(* 8. Quantum Functions *)
WriteString[output, "=== QUANTUM MECHANICS ===\n"];
WriteString[output, ToString[Names["System`*Quantum*"], InputForm] <> "\n\n"];

(* 9. Relativity/Tensor Functions *)
WriteString[output, "=== RELATIVITY & TENSORS ===\n"];
relFuncs = Join[
  Names["System`*Metric*"],
  Names["System`*Tensor*"],
  Names["System`*Riemann*"],
  Names["System`*Ricci*"],
  Names["System`*Einstein*"],
  Names["System`*Christoffel*"]
];
WriteString[output, ToString[relFuncs, InputForm] <> "\n\n"];

(* 10. Electromagnetism *)
WriteString[output, "=== ELECTROMAGNETISM ===\n"];
emFuncs = Join[
  Names["System`*Electric*"],
  Names["System`*Magnetic*"],
  Names["System`*Field*"]
];
WriteString[output, ToString[emFuncs, InputForm] <> "\n\n"];

(* 11. Gravity & Cosmology *)
WriteString[output, "=== GRAVITY & COSMOLOGY ===\n"];
gravFuncs = Join[
  Names["System`*Gravitational*"],
  Names["System`*Cosmological*"]
];
WriteString[output, ToString[gravFuncs, InputForm] <> "\n\n"];

(* 12. Physics-related System symbols *)
WriteString[output, "=== PHYSICS SYSTEM SYMBOLS ===\n"];
physicsSymbols = Select[Names["System`*"], 
  StringContainsQ[#, "Mass" | "Energy" | "Force" | "Charge" | "Momentum" | 
                    "Angular" | "Spin" | "Wavelength" | "Frequency" | "Planck" |
                    "Vacuum" | "Coupling" | "Decay" | "Cross" | "Scatter"] &
];
WriteString[output, ToString[physicsSymbols, InputForm] <> "\n\n"];

(* Count total *)
totalCount = Length[EntityList["PhysicalQuantity"]] + 
             Length[EntityList["PhysicalConstant"]] +
             Length[EntityList["Particle"]] +
             500 + 500 + 500 + 
             Length[funcs] + Length[relFuncs] + Length[emFuncs] + 
             Length[gravFuncs] + Length[physicsSymbols];

WriteString[output, "\n=== TOTAL TERMS RETRIEVED: " <> ToString[totalCount] <> " ===\n"];

Close[output];
Print["Physics terms saved to wolfram_physics_terms_5034.txt"];
Print["Total terms: ", totalCount];
