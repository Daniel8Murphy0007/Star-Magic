"""Fix double-backslash LaTeX escaping in PAPER_674-687 and PAPER_702-715 whitepapers."""
import os

ROOT = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
WP_DIR = os.path.join(ROOT, "whitepapers")

cls_map = [
    (674, "UQFFComparedToLIGOData"),
    (675, "UQFFComparedToGW170817"),
    (676, "UQFFComparedToGW190425"),
    (677, "UQFFPredictionsForLISA"),
    (678, "LISAVsLIGOComparisons"),
    (679, "AetherSuperfluidDynamics"),
    (680, "VortexQuantization"),
    (681, "GrossPitaevskiiVortexSimulation"),
    (682, "UQFFStabilityNumericallyForSgrA"),
    (683, "UQFFHawkingTemperatureModulation"),
    (684, "UQFFPrimordialBHEvaporation"),
    (685, "UQFFPBHDarkMatterImplications"),
    (686, "UQFFModulationForM87"),
    (687, "M87MassEvolutionSimulation"),
    # Session 175: PAPER_702-715
    (702, "SaturnRingSystemUQFF"),
    (703, "NGC1275MagneticMonsterUQFF"),
    (704, "HorseheadNebulaBarnard33UQFF"),
    (705, "NGC3603StarCluster2UQFF"),
    (706, "NGC3603StarClusterPrimaryUQFF"),
    (707, "NGC2525BarredSpiral2UQFF"),
    (708, "PillarsOfCreationM16UQFF"),
    (709, "Westerlund2StarClusterUQFF"),
    (710, "NGC2014NGC2020StarformingUQFF"),
    (711, "NGC2014NGC2020Variant2UQFF"),
    (712, "PillarsOfCreationM16v2UQFF"),
    (713, "UQFFKnowledgeBaseKB19"),
    (714, "UQFFKnowledgeBaseKB18"),
    (715, "UQFFKnowledgeBaseKB17"),
    # Session 176: PAPER_716-730
    (716, "UQFFKnowledgeBaseKB1"),
    (717, "UQFFKnowledgeBaseKB2"),
    (718, "UQFFKnowledgeBaseKB3"),
    (719, "UQFFKnowledgeBaseKB4"),
    (720, "UQFFKnowledgeBaseKB5"),
    (721, "UQFFKnowledgeBaseKB6"),
    (722, "UQFFKnowledgeBaseKB8"),
    (723, "UQFFKnowledgeBaseKB9"),
    (724, "UQFFKnowledgeBaseKB10"),
    (725, "UQFFKnowledgeBaseKB11"),
    (726, "UQFFKnowledgeBaseKB12"),
    (727, "UQFFKnowledgeBaseKB13"),
    (728, "UQFFKnowledgeBaseKB14"),
    (729, "UQFFKnowledgeBaseKB15"),
    (730, "UQFFKnowledgeBaseKB16"),
    # Session 177: PAPER_731
    (731, "NGC1316MergerEvolution"),
]

targets = (
    [(n, cls, os.path.join(ROOT, f"PAPER_{n}_{cls}.md")) for n, cls in cls_map]
    + [(n, cls, os.path.join(WP_DIR, f"PAPER_{n}_{cls}.md")) for n, cls in cls_map]
)

BS2 = chr(92) + chr(92)   # two backslash characters
BS1 = chr(92)              # one backslash character

for n, cls, path in targets:
    if not os.path.exists(path):
        print(f"MISSING: {path}")
        continue
    with open(path, "r", encoding="utf-8") as f:
        content = f.read()

    count_before = content.count(BS2)
    new = content.replace(BS2, BS1)
    count_after = new.count(BS2)

    if new == content:
        print(f"  {n} ({os.path.basename(os.path.dirname(path))}): no double-backslashes found")
    else:
        with open(path, "w", encoding="utf-8") as f:
            f.write(new)
        print(f"  {n}: fixed {count_before} double-backslash -> single (remaining pairs: {count_after})")

print("Done.")
