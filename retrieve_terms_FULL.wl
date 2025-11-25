(* Extended Wolfram Physics Terms Retrieval - Target 5000+ *)
output = OpenWrite["wolfram_physics_terms_FULL.txt"];

(* Core Entity Types *)
WriteString[output, "=== PHYSICAL QUANTITIES ===\n"];
pq = EntityList["PhysicalQuantity"];
WriteString[output, ToString[pq, InputForm] <> "\n\n"];

WriteString[output, "=== PHYSICAL CONSTANTS ===\n"];
pc = EntityList["PhysicalConstant"];
WriteString[output, ToString[pc, InputForm] <> "\n\n"];

WriteString[output, "=== PARTICLES ===\n"];
particles = EntityList["Particle"];
WriteString[output, ToString[particles, InputForm] <> "\n\n"];

WriteString[output, "=== ISOTOPES (First 1000) ===\n"];
isotopes = Take[EntityList["Isotope"], UpTo[1000]];
WriteString[output, ToString[isotopes, InputForm] <> "\n\n"];

WriteString[output, "=== STARS (First 1000) ===\n"];
stars = Take[EntityList["Star"], UpTo[1000]];
WriteString[output, ToString[stars, InputForm] <> "\n\n"];

WriteString[output, "=== GALAXIES (First 1000) ===\n"];
galaxies = Take[EntityList["Galaxy"], UpTo[1000]];
WriteString[output, ToString[galaxies, InputForm] <> "\n\n"];

WriteString[output, "=== EXOPLANETS (First 500) ===\n"];
exoplanets = Take[EntityList["Exoplanet"], UpTo[500]];
WriteString[output, ToString[exoplanets, InputForm] <> "\n\n"];

WriteString[output, "=== PLANETARY MOONS ===\n"];
moons = EntityList["PlanetaryMoon"];
WriteString[output, ToString[moons, InputForm] <> "\n\n"];

WriteString[output, "=== SPECIAL FUNCTIONS ===\n"];
specialFuncs = Join[
  Names["System`*Bessel*"], Names["System`*Airy*"],
  Names["System`*Spherical*"], Names["System`*Legendre*"],
  Names["System`*Hermite*"], Names["System`*Laguerre*"],
  Names["System`*Jacobi*"], Names["System`*Chebyshev*"],
  Names["System`*Hypergeometric*"], Names["System`*Zeta*"],
  Names["System`*Gamma*"], Names["System`*Beta*"],
  Names["System`*Elliptic*"], Names["System`*Mathieu*"],
  Names["System`*Struve*"], Names["System`*Whittaker*"]
];
WriteString[output, ToString[DeleteDuplicates[specialFuncs], InputForm] <> "\n\n"];

WriteString[output, "=== QUANTUM MECHANICS ===\n"];
quantumFuncs = Names["System`*Quantum*"];
WriteString[output, ToString[quantumFuncs, InputForm] <> "\n\n"];

WriteString[output, "=== TENSOR & DIFFERENTIAL GEOMETRY ===\n"];
tensorFuncs = Join[
  Names["System`*Tensor*"], Names["System`*Metric*"],
  Names["System`*Riemann*"], Names["System`*Ricci*"],
  Names["System`*Einstein*"], Names["System`*Christoffel*"],
  Names["System`*Covariant*"], Names["System`*Geodesic*"],
  Names["System`*Curvature*"]
];
WriteString[output, ToString[DeleteDuplicates[tensorFuncs], InputForm] <> "\n\n"];

WriteString[output, "=== ELECTROMAGNETISM ===\n"];
emFuncs = Join[
  Names["System`*Electric*"], Names["System`*Magnetic*"],
  Names["System`*Field*"], Names["System`*Maxwell*"],
  Names["System`*Poynting*"], Names["System`*Lorentz*"]
];
WriteString[output, ToString[DeleteDuplicates[emFuncs], InputForm] <> "\n\n"];

WriteString[output, "=== THERMODYNAMICS & STATISTICAL MECHANICS ===\n"];
thermoFuncs = Join[
  Names["System`*Entropy*"], Names["System`*Temperature*"],
  Names["System`*Heat*"], Names["System`*Boltzmann*"],
  Names["System`*Fermi*"], Names["System`*Bose*"],
  Names["System`*Partition*"]
];
WriteString[output, ToString[DeleteDuplicates[thermoFuncs], InputForm] <> "\n\n"];

WriteString[output, "=== GRAVITY & RELATIVITY ===\n"];
gravFuncs = Join[
  Names["System`*Gravitational*"], Names["System`*Cosmological*"],
  Names["System`*Schwarzschild*"], Names["System`*Kerr*"]
];
WriteString[output, ToString[DeleteDuplicates[gravFuncs], InputForm] <> "\n\n"];

WriteString[output, "=== NUCLEAR PHYSICS ===\n"];
nuclearFuncs = Join[
  Names["System`*Nuclear*"], Names["System`*Atomic*"],
  Names["System`*Decay*"], Names["System`*Fission*"],
  Names["System`*Fusion*"]
];
WriteString[output, ToString[DeleteDuplicates[nuclearFuncs], InputForm] <> "\n\n"];

WriteString[output, "=== PHYSICS SYMBOLS (Mass, Energy, Force, etc.) ===\n"];
physicsSymbols = Select[Names["System`*"], 
  StringContainsQ[#, "Mass" | "Energy" | "Force" | "Charge" | "Momentum" | 
                    "Angular" | "Spin" | "Wavelength" | "Frequency" | "Planck" |
                    "Vacuum" | "Coupling" | "Decay" | "Cross" | "Scatter" |
                    "Potential" | "Kinetic" | "Harmonic" | "Oscillat" | "Wave" |
                    "Photon" | "Electron" | "Proton" | "Neutron" | "Quark"] &
];
WriteString[output, ToString[physicsSymbols, InputForm] <> "\n\n"];

WriteString[output, "=== OPTICS & WAVES ===\n"];
opticsFuncs = Join[
  Names["System`*Refract*"], Names["System`*Reflect*"],
  Names["System`*Diffract*"], Names["System`*Interfere*"],
  Names["System`*Polariz*"], Names["System`*Lens*"]
];
WriteString[output, ToString[DeleteDuplicates[opticsFuncs], InputForm] <> "\n\n"];

WriteString[output, "=== FLUID DYNAMICS ===\n"];
fluidFuncs = Join[
  Names["System`*Flow*"], Names["System`*Viscosity*"],
  Names["System`*Turbulence*"], Names["System`*Reynolds*"],
  Names["System`*Navier*"], Names["System`*Stokes*"]
];
WriteString[output, ToString[DeleteDuplicates[fluidFuncs], InputForm] <> "\n\n"];

(* Count total *)
totalCount = Length[pq] + Length[pc] + Length[particles] + 
             Length[isotopes] + Length[stars] + Length[galaxies] +
             Length[exoplanets] + Length[moons] +
             Length[DeleteDuplicates[specialFuncs]] + 
             Length[quantumFuncs] + Length[DeleteDuplicates[tensorFuncs]] + 
             Length[DeleteDuplicates[emFuncs]] + Length[DeleteDuplicates[thermoFuncs]] +
             Length[DeleteDuplicates[gravFuncs]] + Length[DeleteDuplicates[nuclearFuncs]] +
             Length[physicsSymbols] + Length[DeleteDuplicates[opticsFuncs]] +
             Length[DeleteDuplicates[fluidFuncs]];

WriteString[output, "\n=== TOTAL TERMS RETRIEVED: " <> ToString[totalCount] <> " ===\n"];

Close[output];
Print["Extended physics terms saved to wolfram_physics_terms_FULL.txt"];
Print["Total terms: ", totalCount];
