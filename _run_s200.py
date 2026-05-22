import json
import CondensedPhysics4 as cp4

classes = [
    "WaterReactorBirkelandH2ElectrolysisEfficiencyCalc",     # PAPER_863
    "LRCPseudoMonopoleSparkGapResonanceCalc",                # PAPER_864
    "FieldGeneratorSpookyNonLocalTempDropCalc",              # PAPER_865
    "DCEACEReversalNdFeBCaduceusMotorCalc",                  # PAPER_866
    "MosquitoBioThermalEfficiencyBenchmarkCalc",             # PAPER_867
    "TopoconductorQuantumCoolingComparisonCalc",             # PAPER_868
    "MilkyWay82DayStarTrackingUFTCalc",                      # PAPER_869
]
for c in classes:
    cls = getattr(cp4, c, None)
    if cls is None:
        print(f"=== {c} === NOT FOUND")
        continue
    print(f"\n=== {c} ===")
    try:
        out = cls().compute()
        print(json.dumps(out, indent=2, default=str))
    except Exception as e:
        print(f"ERROR: {e!r}")
