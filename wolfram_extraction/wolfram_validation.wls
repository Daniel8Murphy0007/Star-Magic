(* =========================================================================== *)
(* WOLFRAM KNOWLEDGEBASE VALIDATION SCRIPT                                     *)
(* Generated: 11/30/2025 15:25:05                                                      *)
(* Purpose: Validate UQFF extracted entities and retrieve canonical values     *)
(* =========================================================================== *)

(* Initialize results list *)
validatedEntities = <||>;

(* Physical Constants Validation *)
Print["Validating Physical Constants..."];
constants = {
  "PI",
};

(* Query Wolfram for each constant *)
constantResults = Association[
  Table[
    name -> <|
      "Exists" -> EntityTypeName[name, "PhysicalConstant"] =!= Missing["NotAvailable"],
      "Value" -> If[EntityTypeName[name, "PhysicalConstant"] =!= Missing["NotAvailable"],
        QuantityMagnitude[Entity["PhysicalConstant", name]["Value"]], 
        Null
      ]
    |>,
    {name, constants}
  ]
];

(* Astrophysical Systems Validation *)
Print["Validating Astrophysical Systems..."];
systems = {
  "Magnetar",
};

systemResults = Association[
  Table[
    name -> <|
      "Exists" -> EntityTypeName[name, "Galaxy"] =!= Missing["NotAvailable"] || 
                  EntityTypeName[name, "Star"] =!= Missing["NotAvailable"],
      "Type" -> Which[
        EntityTypeName[name, "Galaxy"] =!= Missing["NotAvailable"], "Galaxy",
        EntityTypeName[name, "Star"] =!= Missing["NotAvailable"], "Star",
        True, "Unknown"
      ]
    |>,
    {name, systems}
  ]
];

(* Export results *)
Export["c:/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic/wolfram_extraction/wolfram_validation_results.json", 
  <|"Constants" -> constantResults, "Systems" -> systemResults|>
];

Print["Validation complete. Results exported."];
