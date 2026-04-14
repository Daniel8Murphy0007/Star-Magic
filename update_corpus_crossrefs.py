#!/usr/bin/env python3
"""
update_corpus_crossrefs.py — Batch cross-reference updater for PAPER_001-1081

Session 225 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Scans all whitepapers and appends a "Session 225 Cross-References" section
linking each paper to relevant papers from the 1000-1081 range based on
tag/keyword matching.

Insertion point: BEFORE the "## Appendix: Session 204 Codebase Upgrade
Reference" section if present, otherwise at end of file.

Idempotent: Will not re-add if section already exists.
────────────────────────────────────────────────────────────────────────────────
"""

import os
import re
import glob

REPO = r"c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
WP_DIR = os.path.join(REPO, "whitepapers")

# ── Topic → Paper mapping for PAPER_1000-1081 ─────────────────────────────
# Each entry: (paper_id, short_title, [keywords])
# Keywords are matched case-insensitively against paper tags + title + first 80 lines

CROSSREF_MAP = [
    # Gravitational waves / mergers
    (1000, "NS Merger F_U_Bi Strain Suppression & BCS Gap",
     ["GW190425", "neutron.star.merger", "BCS.gap", "strain.suppress"]),
    (1001, "SMBH Binary Merger F_U_Bi Phonon Damping",
     ["SMBH.merge", "binary.inspiral", "phonon.damp", "BH.merger"]),
    (1011, "GW170817 NS Merger F_U_Bi_i 66.7% Strain Reduction",
     ["GW170817", "strain.reduc", "LIGO", "BNS"]),
    (1012, "GW190425 Upgraded F_U_Bi_i with S26(3)",
     ["GW190425", "Ramanujan.correct"]),
    (1014, "SMBH Merger Inspiral-Coalescence-Ringdown",
     ["SMBH.merge", "inspiral", "coalescence", "ringdown", "QNM"]),
    (1022, "GW Phonon Strain SCm Modulation of h(t)",
     ["gravitational.wave", "strain.modul", "phonon.strain", "h.t"]),

    # AGN / jets / quasars
    (1002, "AGN Buoyancy-Corrected Eddington Luminosity",
     ["AGN.accret", "Eddington", "buoyancy.correct"]),
    (1009, "3C273 AGN F_U_Bi_i Jet Modulation",
     ["3C273", "AGN.jet", "jet.modul"]),
    (1010, "TON618 AGN F_U_Bi_i Jet Modulation",
     ["TON618", "ultramassive", "AGN.jet"]),
    (1037, "AGN Buoyancy Jet Calculator — SCm Jet Launching",
     ["AGN.jet", "Blandford.Znajek", "jet.launch"]),
    (1048, "M-Sigma Phonon-Corrected Relation",
     ["M.sigma", "velocity.dispersion", "SMBH.mass.power"]),

    # QGP / particle physics
    (1004, "QGP Vacuum Density with SCm S26 Phonon Coupling",
     ["QGP", "quark.gluon", "deconfinement", "vacuum.density"]),
    (1005, "Yang-Mills Mass Gap via SCm BCS Phonon Coupling",
     ["Yang.Mills", "mass.gap", "QCD", "running.coupling"]),
    (1006, "ALICE Multiplicity SCm Phonon Scaling",
     ["ALICE", "Pb.Pb", "multiplicity", "centrality"]),
    (1007, "Deconfinement Phase Diagram SCm Phonon Boundary",
     ["QGP.phase", "deconfinement", "hadron"]),
    (1013, "QGP ALICE Centrality F_U_Bi_i dN/deta Scaling",
     ["ALICE.central", "dN.deta"]),
    (1059, "Color Glass Condensate BK Saturation SCm",
     ["CGC", "BK.equation", "gluon.saturat", "HERA"]),
    (1064, "Resummation Effective Coupling BFKL/Sudakov SCm",
     ["BFKL", "Sudakov", "resummation", "QCD"]),

    # Galaxy clusters / ICM
    (1039, "SCm Galaxy Cluster Buoyancy Profile ICM Beta-Model",
     ["galaxy.cluster", "ICM", "beta.model", "cluster.buoyancy"]),
    (1040, "SCm Cluster Merger Shock Mach Number Phonon Damping",
     ["cluster.merger", "shock", "Mach.number"]),
    (1041, "SCm Cool-Core Buoyancy Balance AGN Feedback",
     ["cool.core", "AGN.feedback", "cooling.flow"]),
    (1044, "SCm Cluster Thermal SZ Effect Compton-y Phonon",
     ["SZ.effect", "Compton.y", "CMB.cluster", "Sunyaev"]),
    (1045, "SCm Cluster Radio Relic Polarization",
     ["radio.relic", "polarization", "cluster.shock"]),
    (1046, "SCm Cluster Lensing Mass Phonon Correction",
     ["weak.lensing", "convergence.map", "cluster.mass", "kappa.map"]),
    (1079, "Galaxy Cluster Cooling-Flow Buoyancy Suppression",
     ["cooling.flow", "AGN.feedback", "ICM", "Perseus", "cluster.cool"]),

    # Dark matter
    (1015, "SCm Dark Matter Halos NFW Rotation Curve",
     ["dark.matter", "NFW.halo", "rotation.curve", "halo"]),
    (1019, "Dark Matter Phonon Buoyancy NFW Coupling",
     ["dark.matter.phonon", "NFW.coupling", "DM.buoyancy"]),

    # Dark energy
    (1076, "SCm Dark Energy with Phonon Linewidth Gamma-Modulation",
     ["dark.energy", "Gamma.modul", "LCDM", "equation.state", "w.z"]),

    # SCm phonon phenomenology
    (1020, "Cosmic Ray Phonon Acceleration DSA Spectrum",
     ["cosmic.ray", "DSA", "accelerat", "spectrum"]),
    (1021, "Pulsar Timing Phonon TOA Residual",
     ["pulsar.timing", "PTA", "TOA.residual"]),
    (1023, "Neutrino Oscillation Phonon PMNS Matrix SCm",
     ["neutrino.oscill", "PMNS.matrix", "neutrino.mix"]),
    (1024, "Magnetar Giant Flare SCm Phonon Reservoir",
     ["magnetar", "giant.flare", "SGR.energy"]),
    (1025, "Black Hole Shadow Phonon Photon Ring Correction",
     ["BH.shadow", "EHT", "photon.ring"]),
    (1026, "Reionization Bubble Phonon Stromgren Sphere",
     ["reionization", "Stromgren", "cosmic.dawn"]),
    (1027, "Tidal Disruption Event SCm Fallback",
     ["tidal.disrupt", "TDE", "fallback"]),
    (1028, "Cosmic String Gravitational Lens SCm",
     ["cosmic.string", "deficit.angle", "string.lens"]),
    (1029, "Barocentric Earth Orbital Buoyancy",
     ["Earth.orbit", "solar.system", "barycent"]),
    (1030, "Quantum Gravity Minimum Length GUP-SCm",
     ["quantum.gravity", "GUP", "Planck.length", "minimum.length"]),
    (1031, "Photon Sphere Phonon Orbital SCm",
     ["photon.sphere", "critical.impact", "critical.orbit"]),
    (1032, "ISM Dust Grain Buoyancy SCm",
     ["ISM.dust", "dust.grain", "interstellar.particle"]),
    (1033, "Galactic Bar Resonance SCm Pattern Speed",
     ["galactic.bar", "resonance", "ILR", "pattern.speed"]),
    (1034, "FRB Dispersion Measure Buoyancy",
     ["FRB", "dispersion.measure", "IGM"]),
    (1035, "Kilonova Buoyancy Light Curve r-Process",
     ["kilonova", "r.process", "light.curve"]),
    (1036, "Primordial Nucleosynthesis BBN Phonon",
     ["BBN", "nucleosynthesis", "primordial", "helium"]),
    (1038, "White Dwarf Crystallization Buoyancy",
     ["white.dwarf", "crystalliz", "Gaia", "latent.heat"]),
    (1042, "Mock-Theta Phonon Partition Ramanujan q-Series",
     ["mock.theta", "q.series", "Ramanujan.partition"]),
    (1043, "F_U_Bi_i Multi-System Buoyancy Curve Sweep",
     ["multi.system", "buoyancy.curve", "gamma.sweep"]),
    (1047, "Type Iax Supernova Buoyancy Reversal",
     ["Type.Iax", "supernova", "deflagration", "buoyancy.reversal"]),
    (1072, "SCm Activation Function Phonon Threshold",
     ["activation.function", "Heaviside", "threshold"]),
    (1073, "SCm Phonon-Driven Inflation Vacuum Buoyancy",
     ["inflation", "inflaton", "Hubble", "slow.roll", "vacuum.buoyancy"]),

    # LENR / nuclear
    (1060, "VDS LENR Isotopic Transmutation Chain",
     ["VDS.LENR", "isotopic", "transmutation"]),
    (1061, "Kozima SCm Integration Neutron-Drop",
     ["Kozima", "neutron.drop", "LENR"]),
    (1081, "SCm LENR COP Linewidth Parametric",
     ["LENR.COP", "COP.parametric", "Fleischmann", "micro.plasmoid"]),

    # Formal theory / Lagrangian
    (1003, "Spectral Ladder Merger 26-State Hierarchy",
     ["spectral.ladder", "26.state", "energy.hierarchy"]),
    (1051, "Universal Duality SCm-UA Synthesis",
     ["duality", "SCm.UA", "expansion.collapse", "Sign.E_net"]),
    (1052, "TQFT Anyon Braiding Chern-Simons",
     ["TQFT", "anyon", "Chern.Simons", "topological"]),
    (1053, "Swampland Conjecture SCm Bridge",
     ["Swampland", "WGC", "de.Sitter"]),
    (1054, "SUSY Breaking Soft Terms SCm Mediation",
     ["SUSY", "soft.terms", "supersymmetry"]),
    (1055, "cMERA Entanglement RG Holographic SCm",
     ["cMERA", "entanglement.RG", "holographic"]),
    (1056, "Quantum Error Correction Topological SCm",
     ["QEC", "surface.code", "topological.code", "error.correct"]),
    (1057, "Non-Commutative Geometry Matrix Model SCm",
     ["NCG", "non.commutative", "spectral.triple"]),
    (1058, "LQG Ashtekar Area Spectrum SCm",
     ["LQG", "loop.quantum", "Ashtekar", "area.spectrum", "Barbero"]),
    (1063, "Higher-Curvature Gravity EFT Gauss-Bonnet SCm",
     ["Gauss.Bonnet", "f.R.", "Horndeski", "higher.curvat"]),
    (1065, "Buoyancy Lagrangian EOM Variational Derivation",
     ["Lagrangian", "variational", "Euler.Lagrange", "EOM"]),
    (1066, "UQFF Lagrangian First Principles Field Theory",
     ["Lagrangian.deriv", "first.principles", "field.theory"]),

    # Number systems / bridges
    (1069, "VDS-DVP-BSH Hybrid Calculator Unified",
     ["VDS", "DVP", "BSH", "number.system"]),
    (1070, "Yang-Mills Mass Gap VDS Bridge",
     ["Yang.Mills.VDS", "vacuum.density.gap", "confinement"]),
    (1080, "Ramanujan Binomial Expansion Proof R_n^{(26,3)}",
     ["Ramanujan.binom", "binomial.expans", "R_n.26", "polylog"]),

    # QCalcGeom
    (1067, "QCalc Geometry Bridge Python Solver",
     ["QCalcGeom", "QCalc.geometry", "geometry.bridge"]),
    (1078, "QCalcGeom Master Equation Derivation",
     ["QCalcGeom", "BSFG.crossover", "26D.compactific", "master.equation"]),

    # Wormholes
    (1062, "Wormhole Traversability Morris-Thorne SCm",
     ["wormhole", "Morris.Thorne", "traversab", "exotic.matter"]),

    # ALMA / observational
    (1049, "Source10 GPU DPM Spectral Atlas ALMA Overlay",
     ["DPM.spectral", "ALMA.overlay", "GPU.atlas"]),
    (1074, "GPU-Vectorized DPM S26 Spectral Atlas",
     ["GPU.DPM", "spectral.atlas", "vectorize"]),
    (1077, "ALMA Cycle 12 F_U_Bi_i Line Profile Validation",
     ["ALMA.Cycle", "ALMA.valid", "line.profile.valid"]),
    (1075, "3D Volumetric MUGE Gravitational Field Generator",
     ["3D.MUGE", "volumetric.gravity", "MUGE.field"]),

    # JWST
    (1071, "JWST Synthesis Multi-Instrument UQFF",
     ["JWST", "NIRCam", "MIRI", "NIRSpec"]),

    # Wolfram
    (1068, "Wolfram Physics Bridge WSTP Symbolic Export",
     ["Wolfram.bridge", "WSTP", "symbolic.export", "Mathematica"]),

    # Production / scaling
    (1008, "Production Scaling v14 — 600k calc/s 24 Kernels",
     ["production.scaling", "600k", "benchmark.kernel"]),
    (1018, "Production Scaling v15 — 650k calc/s 30 Kernels",
     ["production.scaling", "650k", "benchmark.kernel"]),
    (1050, "MUGE F_U_Bi_i Unified 9-System Synthesis",
     ["MUGE.unified", "9.system", "multiplier.table"]),
]

# ── Broader tag-based matching rules ──────────────────────────────────────
# Map YAML tags to relevant 1000-1081 papers

TAG_MAP = {
    "GW": [1000, 1001, 1011, 1012, 1014, 1022],
    "gravitational-wave": [1000, 1001, 1011, 1012, 1014, 1022],
    "merger": [1000, 1001, 1011, 1012, 1014, 1035],
    "neutron-star": [1000, 1011, 1012],
    "LIGO": [1011, 1022],
    "AGN": [1002, 1009, 1010, 1037, 1048, 1041, 1079],
    "jet": [1009, 1010, 1037],
    "quasar": [1009, 1010],
    "QGP": [1004, 1005, 1006, 1007, 1013, 1059],
    "QCD": [1005, 1059, 1064],
    "Yang-Mills": [1005, 1070],
    "ALICE": [1006, 1013],
    "cluster": [1039, 1040, 1041, 1044, 1045, 1046, 1079],
    "ICM": [1039, 1040, 1041, 1044, 1079],
    "dark-matter": [1015, 1019],
    "dark-energy": [1076],
    "NFW": [1015, 1019],
    "halo": [1015, 1019],
    "cosmic-ray": [1020],
    "pulsar": [1021],
    "PTA": [1021],
    "neutrino": [1023],
    "magnetar": [1024],
    "EHT": [1025],
    "shadow": [1025],
    "reionization": [1026],
    "TDE": [1027],
    "cosmic-string": [1028],
    "GUP": [1030],
    "FRB": [1034],
    "kilonova": [1035],
    "BBN": [1036],
    "nucleosynthesis": [1036],
    "white-dwarf": [1038],
    "supernova": [1047],
    "Type-Iax": [1047],
    "inflation": [1073],
    "inflaton": [1073],
    "LENR": [1060, 1061, 1081],
    "Kozima": [1061],
    "TQFT": [1052],
    "anyon": [1052],
    "Swampland": [1053],
    "SUSY": [1054],
    "supersymmetry": [1054],
    "cMERA": [1055],
    "holographic": [1055],
    "QEC": [1056],
    "NCG": [1057],
    "LQG": [1058],
    "Ashtekar": [1058],
    "Gauss-Bonnet": [1063],
    "Horndeski": [1063],
    "Lagrangian": [1065, 1066],
    "variational": [1065],
    "wormhole": [1062],
    "Morris-Thorne": [1062],
    "VDS": [1069, 1070, 1080],
    "DVP": [1069, 1080],
    "BSH": [1069, 1080],
    "Ramanujan": [1042, 1080],
    "QCalcGeom": [1067, 1078],
    "DPM": [1049, 1074],
    "ALMA": [1049, 1074, 1077],
    "JWST": [1071],
    "MUGE": [1050, 1075],
    "WSTP": [1068],
    "Wolfram": [1068],
    "phonon": [1020, 1021, 1022, 1023, 1024, 1072, 1073],
    "SCm": [1072, 1073],
    "vacuum": [1004, 1049, 1069],
    "buoyancy": [1037, 1043, 1065, 1079],
    "SZ": [1044],
    "Sunyaev": [1044],
    "lensing": [1021, 1028, 1046],
    "M-sigma": [1048],
}

# Build paper_id → (title) lookup
PAPER_LOOKUP = {p[0]: p[1] for p in CROSSREF_MAP}

S204_MARKER = "## Appendix: Session 204 Codebase Upgrade Reference"
S225_MARKER = "## Appendix: Session 225 Cross-References (PAPER_1000–1081)"


def extract_paper_number(filename):
    """Extract paper number from filename like PAPER_001_... or PAPER_1081_..."""
    m = re.search(r'PAPER_(\d+)', filename)
    return int(m.group(1)) if m else None


def extract_tags(content):
    """Extract YAML tags from frontmatter."""
    m = re.search(r'^tags:\s*\[([^\]]*)\]', content, re.MULTILINE)
    if m:
        return [t.strip().strip('"').strip("'") for t in m.group(1).split(',')]
    return []


def extract_title(content):
    """Extract title from YAML frontmatter."""
    m = re.search(r'^title:\s*["\']?([^"\'\n]*)["\']?', content, re.MULTILINE)
    return m.group(1).strip() if m else ""


def keyword_match(content_lower, keywords):
    """Check if any keyword pattern matches in content (case-insensitive)."""
    for kw in keywords:
        # Convert dot-separated keywords to regex patterns
        pattern = kw.replace('.', r'[\s._-]*')
        if re.search(pattern, content_lower):
            return True
    return False


def find_relevant_crossrefs(paper_num, tags, title, content_head):
    """Find relevant papers from 1000-1081 for this paper."""
    relevant = set()

    # Skip papers that are themselves in the 1000-1081 range (they reference each other)
    # But still process them for cross-references to other 1000-1081 papers

    # 1. Tag-based matching
    for tag in tags:
        tag_clean = tag.strip()
        if tag_clean in TAG_MAP:
            for ref in TAG_MAP[tag_clean]:
                if ref != paper_num:
                    relevant.add(ref)

    # 2. Keyword-based matching on title + first ~80 lines
    search_text = (title + "\n" + content_head).lower()
    for paper_id, _title, keywords in CROSSREF_MAP:
        if paper_id != paper_num and keyword_match(search_text, keywords):
            relevant.add(paper_id)

    return sorted(relevant)


def build_crossref_section(relevant_papers):
    """Build the cross-reference appendix section."""
    lines = [
        "",
        "---",
        "",
        S225_MARKER,
        "",
        "> *Auto-generated cross-reference appendix linking this paper to",
        "> Sessions 204–225 extensions (PAPER_1000–1081). Added by",  
        "> `update_corpus_crossrefs.py` (Session 225, April 2026).*",
        "",
    ]

    # Group by theme
    themes = {
        "Gravitational Waves & Mergers": [1000, 1001, 1011, 1012, 1014, 1022],
        "AGN / Jets / Quasars": [1002, 1009, 1010, 1037, 1048],
        "QGP / Particle Physics": [1004, 1005, 1006, 1007, 1013, 1059, 1064],
        "Galaxy Clusters / ICM": [1039, 1040, 1041, 1044, 1045, 1046, 1079],
        "Dark Matter & Dark Energy": [1015, 1019, 1076],
        "SCm Phonon Phenomenology": [1020, 1021, 1023, 1024, 1025, 1026, 1027,
                                      1028, 1029, 1030, 1031, 1032, 1033, 1034,
                                      1035, 1036, 1038, 1042, 1043, 1047, 1072, 1073],
        "LENR / Nuclear": [1060, 1061, 1081],
        "Formal Theory & Lagrangian": [1003, 1051, 1052, 1053, 1054, 1055,
                                        1056, 1057, 1058, 1063, 1065, 1066],
        "Number Systems & Bridges": [1067, 1068, 1069, 1070, 1078, 1080],
        "Observational / Validation": [1049, 1050, 1071, 1074, 1075, 1077],
        "Production Scaling": [1008, 1018],
    }

    # Only include themes that have relevant papers
    any_listed = False
    for theme_name, theme_papers in themes.items():
        matching = [p for p in relevant_papers if p in theme_papers]
        if matching:
            if not any_listed:
                lines.append("| Paper | Title |")
                lines.append("|-------|-------|")
                any_listed = True
            for p in matching:
                title = PAPER_LOOKUP.get(p, f"PAPER_{p}")
                lines.append(f"| PAPER_{p:04d} | {title} |")

    if not any_listed:
        # No matches — don't add section
        return None

    lines.append("")
    lines.append(f"*{len(relevant_papers)} cross-reference(s) identified.*")
    lines.append("")

    return "\n".join(lines)


def process_paper(filepath):
    """Process a single paper file, adding cross-references if needed."""
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()

    # Skip if already has S225 cross-references
    if S225_MARKER in content:
        return 0  # already updated

    filename = os.path.basename(filepath)
    paper_num = extract_paper_number(filename)
    if paper_num is None:
        return 0

    tags = extract_tags(content)
    title = extract_title(content)

    # Get first ~80 lines for keyword matching
    content_lines = content.split('\n')
    content_head = '\n'.join(content_lines[:min(80, len(content_lines))])

    # Find relevant cross-references
    relevant = find_relevant_crossrefs(paper_num, tags, title, content_head)
    if not relevant:
        return 0

    # Build cross-reference section
    crossref_section = build_crossref_section(relevant)
    if crossref_section is None:
        return 0

    # Insert before S204 appendix, or at end
    if S204_MARKER in content:
        # Find the position — insert before the "---" that precedes S204
        idx = content.index(S204_MARKER)
        # Look for the preceding "---" separator
        preceding = content[:idx].rstrip()
        if preceding.endswith('---'):
            insert_pos = preceding.rfind('---')
            new_content = content[:insert_pos] + crossref_section + "\n" + content[insert_pos:]
        else:
            new_content = content[:idx] + crossref_section + "\n" + content[idx:]
    else:
        # Append at end
        new_content = content.rstrip() + "\n" + crossref_section

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(new_content)

    return len(relevant)


def main():
    print("=" * 72)
    print("update_corpus_crossrefs.py — Session 225 Cross-Reference Update")
    print("=" * 72)

    # Gather all paper files
    paper_files = []
    for f in glob.glob(os.path.join(WP_DIR, "PAPER_*.md")):
        paper_files.append(f)
    for f in glob.glob(os.path.join(REPO, "PAPER_*.md")):
        paper_files.append(f)

    paper_files.sort(key=lambda x: extract_paper_number(os.path.basename(x)) or 0)

    print(f"Found {len(paper_files)} paper files")
    print()

    updated = 0
    skipped = 0
    total_refs = 0

    for filepath in paper_files:
        filename = os.path.basename(filepath)
        paper_num = extract_paper_number(filename)
        refs = process_paper(filepath)
        if refs > 0:
            print(f"  [+] {filename}: {refs} cross-references added")
            updated += 1
            total_refs += refs
        else:
            skipped += 1

    print()
    print(f"  Updated: {updated} papers")
    print(f"  Skipped: {skipped} papers (no matches or already updated)")
    print(f"  Total cross-references added: {total_refs}")
    print("=" * 72)


if __name__ == "__main__":
    main()
