"""Job B Phase B1 v2: Filename-only categorization.

READ-ONLY. Filenames in this repo are intentionally descriptive,
so matching against filename alone gives a clean signal without
frontmatter false positives (e.g. "paper_index:" YAML key).

Bucket priority (first match wins):
  C = SI constants derivations
  F = Index / catalog / master synthesis
  G = Millennium Prize
  D = KK / extra-dim / sub-mm / Yukawa
  E = Falsifiability / predictions / observational tests
  B = Unified field / Lagrangian / EFE / closure DERIVATIONS
  A = Cosmology / Lambda / dark energy / vacuum energy
  I = LENR / nuclear / atomic (framework derivation)
  J = Audits / doctrine / methodology
  H = Specific system or framework application  (catch-all)
"""
import re
import csv
import pathlib
from collections import Counter

ROOT = pathlib.Path(__file__).parent
WPDIR = ROOT / "whitepapers"
OUT = ROOT / "_job_b_categorization.csv"

RULES = [
    ("C", "Add three-anchor SI closure banner (Sess 237-241); STRUCTURAL tier",
     ["_Planck_Constant", "_Fine_Structure_Constant", "_Speed_of_Light",
      "_Gravitational_Constant", "Three_Anchor", "SI_Closure",
      "FineStructureConstant", "Newton_G_"]),

    ("F", "Update tables; add P6-P14, CP4 #254-#264, ref PAPER_1167/1170/1174",
     ["_Index_", "_Catalog_", "_Catalog.md", "Catalog_All",
      "_All29_", "All_29_Systems", "Master_Synthesis", "MasterSynthesis",
      "_Registry_", "_Registry.md", "Summary_of_All",
      "Complete_UQFF_Equations", "AllPapers", "Equation_Catalog",
      "Knowledge_Base", "KnowledgeBase", "PaperSummary",
      "Variable_Reference_Table", "Canonical_Body_Reference"]),

    ("G", "Verify v5.78 closed-Lagrangian/ledger cross-ref present (Sess 225 baseline)",
     ["Millennium", "Yang_Mills", "YangMills", "Yang-Mills",
      "Navier_Stokes", "NavierStokes", "Navier-Stokes",
      "Riemann_Hypothesis", "RiemannHypothesis",
      "P_vs_NP", "PvsNP", "Hodge_Conjecture",
      "Mass_Gap_UQFF", "MassGap_UQFF", "MassGap_SCm"]),

    ("D", "Forward-pointer to PAPER_1171/1172 (xi=13/3 R26+KK lock) + PAPER_1173",
     ["Kaluza", "Kaluza_Klein", "KK_Tower", "KK_Regulator", "KK_Mode",
      "KK_Hbar", "Extra_Dimension", "Extra-Dimension",
      "Sub_mm_Gravity", "Sub-mm", "SubMM",
      "_Yukawa_", "Yukawa_Coupling", "Yukawa_Bound",
      "Compactification", "Calabi_Yau", "Calabi-Yau"]),

    ("E", "Forward-pointer to PAPER_1174 P1-P14 suite + PAPER_1177-1180",
     ["Falsifi", "Falsifiable", "Falsification",
      "Falsifier", "Predictions_Master", "Prediction_Suite",
      "_P06_", "_P07_", "_P08_", "_P09_", "_P10_", "_P11_", "_P12_", "_P13_", "_P14_",
      "_P6_", "_P7_", "_P8_", "_P9_",
      "Joint_Falsifier", "QuadrupleFalsifier", "Triple_Falsifier",
      "DESI_Y5", "Euclid_Sigma", "CMB_S4", "CMB-S4",
      "LIGO_O5", "Ringdown_Spectral", "mu_Distortion",
      "EP10_", "EP11_", "EP12_", "EP13_", "EP14_"]),

    ("B", "Add v5.78 closed-Lagrangian cross-ref (PAPER_1159-1167 G1-G8) + CP4 #254",
     ["Closed_Lagrangian", "ClosedLagrangian",
      "Lagrangian_Full_Closure", "Lagrangian_Derivation",
      "Lagrangian_Re_derivation", "Lagrangian_Unified",
      "Expansion_Lagrangian", "Euler_Lagrange",
      "Master_Equation_Derivation", "MasterEquation_Derivation",
      "Unified_Field_Theory", "UnifiedFieldTheory",
      "Unification_", "Action_Principle",
      "Mexican_Hat", "MexicanHat",
      "Spontaneous_Symmetry_Breaking", "Gauge_Closure",
      "Phi_Res_Codimension", "F_TRZ_SO5", "Pochhammer_Closure",
      "DPM_SO2", "T22_Moduli", "G1_Closure", "G2_Closure",
      "G3_Closure", "G4_Closure", "G5_Closure", "G6_Closure",
      "G7_Closure", "G8_Closure",
      "Wolfram_Field_Unity"]),

    ("A", "Add closing section: v5.78 27-decade vacuum ledger (PAPER_1170) + xi=13/3",
     ["Cosmological_", "_Cosmological",
      "Cosmolog", "Cosmogenesis",
      "Lambda_CDM", "LambdaCDM", "Lambda-CDM",
      "Lambda_Closure", "_Lambda_",
      "Cosmological_Constant", "CosmologicalConstant",
      "Dark_Energy", "DarkEnergy", "Dark-Energy",
      "Vacuum_Energy", "VacuumEnergy",
      "Vacuum_Density_Series", "Vacuum_Density_Evolution",
      "Vacuum_Ledger", "Energy_Ledger",
      "_CMB_", "_CMB.md",
      "Big_Bang", "BigBang", "BBN_",
      "Inflation_", "_Inflation",
      "Reionization", "Recombination",
      "Universe_Diameter", "Hubble_Tension", "H0_Tension",
      "Quintessence", "_w_de", "Equation_of_State_w",
      "Sigma8_", "_Sigma_8",
      "Cosmic_Dawn", "Cosmic_Constant", "Cosmic_Reionization",
      "Universe_Anchor", "Universal_Anchor"]),

    ("I", "NO UPDATE unless ledger/xi enters; verify only",
     ["LENR_", "_LENR", "Cold_Fusion", "ColdFusion",
      "Kozima", "Holmlid", "Widom_Larsen", "Widom-Larsen",
      "Pons_Fleischmann", "PonsFleischmann", "Mizuno_",
      "Parkhomov",
      "Periodic_Table", "Element_Z_",
      "Beta_Decay", "Cold_Neutron",
      "Atomic_Creation", "AtomicCreation"]),

    ("J", "Case-by-case; add v5.78 reference if doctrinal paper",
     ["_Audit_", "_Audit.md", "Doctrine", "Methodology",
      "Coherence_Audit", "CVW_v",
      "G6_SM_Anchor", "G6_Anchor_Gate",
      "Validation_Pipeline", "_Style_Guide",
      "Notation_Guide", "_Doctrinal"]),
]

ID_RE = re.compile(r"^PAPER_(\d+)", re.IGNORECASE)


def classify_by_filename(name: str) -> tuple[str, str, str]:
    nl = name.lower()
    for bucket, action, kws in RULES:
        for kw in kws:
            if kw.lower() in nl:
                return bucket, action, kw
    return "H", "NO UPDATE (specific system / framework application)", ""


def extract_title(path: pathlib.Path) -> str:
    try:
        with path.open(encoding="utf-8") as f:
            head = "".join([next(f, "") for _ in range(30)])
    except Exception:
        return ""
    m = re.search(r'^title:\s*"?([^"\n]+)"?', head, re.MULTILINE)
    if m:
        return m.group(1).strip()
    m = re.search(r"^#\s+(.+)$", head, re.MULTILINE)
    if m:
        return m.group(1).strip()
    return ""


def main():
    rows = []
    counts = Counter()
    paths = sorted(WPDIR.glob("PAPER_*.md"))
    for p in paths:
        name = p.stem
        m = ID_RE.match(name)
        pid = m.group(1) if m else ""
        title = extract_title(p)
        bucket, action, kw = classify_by_filename(name)
        rows.append({
            "paper_id": pid,
            "filename": p.name,
            "title": title,
            "bucket": bucket,
            "matched_keyword": kw,
            "suggested_action": action,
        })
        counts[bucket] += 1

    with OUT.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "paper_id", "filename", "title", "bucket",
            "matched_keyword", "suggested_action",
        ])
        w.writeheader()
        w.writerows(rows)

    print(f"Total papers scanned: {len(rows)}")
    print(f"Output: {OUT.relative_to(ROOT)}")
    print("")
    labels = {
        "A": "Cosmology / Lambda / dark energy / vacuum",
        "B": "Unified field / Lagrangian / EFE / master eq",
        "C": "SI constants derivations",
        "D": "KK / extra-dim / sub-mm / Yukawa",
        "E": "Falsifiability / predictions",
        "F": "Index / catalog / master synthesis",
        "G": "Millennium Prize",
        "H": "Specific system / framework application (NO UPDATE)",
        "I": "LENR / nuclear / atomic",
        "J": "Audits / doctrine / methodology",
    }
    print("Bucket counts:")
    for b in "CFGDEBAIJH":
        print(f"  {b}  {labels[b]:<54}  {counts.get(b, 0):>5}")
    upd = sum(counts.get(b, 0) for b in "ABCDEFGJ")
    no_upd = counts.get("H", 0) + counts.get("I", 0)
    print("")
    print(f"NEEDS UPDATE / verify (A,B,C,D,E,F,G,J): {upd}")
    print(f"NO UPDATE             (H,I):            {no_upd}")


if __name__ == "__main__":
    main()
