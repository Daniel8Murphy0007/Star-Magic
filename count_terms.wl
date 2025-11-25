(* Count entities in saved file *)
data = Get["wolfram_physics_terms_FULL.txt"];

(* Count each section *)
counts = <|
  "PhysicalQuantities" -> If[KeyExistsQ[data, "pq"], Length[data["pq"]], 0],
  "PhysicalConstants" -> If[KeyExistsQ[data, "pc"], Length[data["pc"]], 0],
  "Particles" -> If[KeyExistsQ[data, "particles"], Length[data["particles"]], 0],
  "Isotopes" -> If[KeyExistsQ[data, "isotopes"], Length[data["isotopes"]], 0],
  "Stars" -> If[KeyExistsQ[data, "stars"], Length[data["stars"]], 0],
  "Galaxies" -> If[KeyExistsQ[data, "galaxies"], Length[data["galaxies"]], 0],
  "Exoplanets" -> If[KeyExistsQ[data, "exoplanets"], Length[data["exoplanets"]], 0],
  "Moons" -> If[KeyExistsQ[data, "moons"], Length[data["moons"]], 0],
  "SpecialFunctions" -> If[KeyExistsQ[data, "sf"], Length[data["sf"]], 0],
  "QuantumFunctions" -> If[KeyExistsQ[data, "qf"], Length[data["qf"]], 0],
  "TensorFunctions" -> If[KeyExistsQ[data, "tf"], Length[data["tf"]], 0],
  "EMFunctions" -> If[KeyExistsQ[data, "em"], Length[data["em"]], 0],
  "ThermoFunctions" -> If[KeyExistsQ[data, "thermo"], Length[data["thermo"]], 0],
  "GravityFunctions" -> If[KeyExistsQ[data, "grav"], Length[data["grav"]], 0],
  "NuclearFunctions" -> If[KeyExistsQ[data, "nuc"], Length[data["nuc"]], 0],
  "PhysicsSymbols" -> If[KeyExistsQ[data, "phys"], Length[data["phys"]], 0],
  "OpticsFunctions" -> If[KeyExistsQ[data, "opt"], Length[data["opt"]], 0],
  "FluidFunctions" -> If[KeyExistsQ[data, "fluid"], Length[data["fluid"]], 0]
|>;

Total[Values[counts]]
