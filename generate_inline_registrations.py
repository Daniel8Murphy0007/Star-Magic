"""
generate_inline_registrations.py
Extracts getName() and getCategory() for the 49 unregistered PhysicsTerm classes
in MAIN_1_CoAnQi.cpp and generates Batch 20 registration calls.
"""
import re

UNREGISTERED = [
    "ReactorEfficiencyTerm", "MuSTerm", "GradMsRTerm", "BjTerm", "OmegaSTTerm",
    "MuJTerm", "UniversalGravity1Term", "UniversalGravity2Term",
    "UniversalGravity3Term", "UniversalGravity4Term", "UniversalBuoyancyTerm",
    "UniversalMagnetismTerm", "UniversalAetherTerm", "UnifiedFieldTerm",
    "MUGECompressedBaseTerm", "MUGEExpansionTerm", "MUGESuperAdjustmentTerm",
    "MUGEEnvelopeTerm", "MUGEUgSumTerm", "MUGECosmologicalTerm",
    "MUGEQuantumTerm", "MUGEFluidTerm", "MUGEPerturbationTerm",
    "MUGEResonanceADPMTerm", "MUGEResonanceATHzTerm", "MUGEResonanceAvacDiffTerm",
    "MUGEResonanceASuperFreqTerm", "MUGEResonanceAAetherResTerm",
    "MUGEResonanceUg4iTerm", "MUGEResonanceAQuantumFreqTerm",
    "MUGEResonanceAAetherFreqTerm", "MUGEResonanceAFluidFreqTerm",
    "MUGEResonanceOscTerm", "MUGEResonanceAExpFreqTerm",
    "MUGEResonanceFTRZTerm", "MUGEResonanceWormholeTerm",
    "TauLeptonDipoleTerm", "CKMVcbTerm", "LFVBranchingTerm",
    "VectorLikeQuarkTerm", "MUGESuperAdjTerm", "MUGECosmTerm",
    "SGR1745MagnetarTerm", "SagittariusAStarTerm", "TapestryStarbirthTerm",
    "Westerlund2ClusterTerm", "PillarsCreationTerm", "RingsRelativityTerm",
    "StudentGuideUniverseTerm",
]

def extract_class_body(content, class_name):
    """Extract full class body using brace counting."""
    pattern = rf'class\s+{re.escape(class_name)}\s*[:\s][^{{]*\{{'
    m = re.search(pattern, content)
    if not m:
        return None
    start = m.start()
    brace_count = 0
    i = m.end() - 1  # The opening brace
    while i < len(content):
        if content[i] == '{':
            brace_count += 1
        elif content[i] == '}':
            brace_count -= 1
            if brace_count == 0:
                return content[start:i+2]  # Include trailing ;
        i += 1
    return None

def extract_return_string(body, method_name):
    """Extract the literal string returned by a method."""
    # Match: methodName() ... return "STRING";
    pattern = rf'{method_name}\s*\(\s*\).*?return\s+"([^"]+)"'
    m = re.search(pattern, body, re.DOTALL)
    if m:
        return m.group(1)
    return None

def infer_name(class_name):
    """Infer a registration name from the class name."""
    name = class_name.replace('Term', '')
    # Add underscores between camelCase words where appropriate
    name = re.sub(r'([a-z])([A-Z])', r'\1_\2', name)
    return name

def infer_category(class_name):
    """Infer a category from the class name."""
    if 'MUGE' in class_name:
        return 'MUGE'
    if 'Resonance' in class_name:
        return 'Resonance'
    if 'Universal' in class_name:
        return 'UQFF-Universal'
    if 'SGR' in class_name or 'Magnetar' in class_name:
        return 'SOURCE4-Magnetar'
    if 'Sagittarius' in class_name:
        return 'SOURCE4-SMBH'
    if 'Tapestry' in class_name or 'Westerlund' in class_name or 'Pillars' in class_name:
        return 'SOURCE4-StarFormation'
    if 'Rings' in class_name or 'Student' in class_name:
        return 'SOURCE4-Cosmology'
    if 'Reactor' in class_name:
        return 'SMBH-Reactor'
    if 'Lepton' in class_name or 'CKM' in class_name or 'LFV' in class_name or 'VectorLike' in class_name:
        return 'BSM-SM-Extension'
    if 'Unified' in class_name or 'Universal' in class_name:
        return 'UQFF-Unified'
    return 'UQFF'

print("Reading MAIN_1_CoAnQi.cpp (this may take a moment)...")
with open('MAIN_1_CoAnQi.cpp', 'r', encoding='utf-8', errors='ignore') as f:
    content = f.read()

print(f"File loaded: {len(content):,} characters")

results = []
for cls in UNREGISTERED:
    body = extract_class_body(content, cls)
    if body is None:
        print(f"WARNING: Could not find class body for {cls}")
        term_name = infer_name(cls)
        category = infer_category(cls)
    else:
        term_name = extract_return_string(body, 'getName')
        if not term_name:
            term_name = infer_name(cls)
        category = extract_return_string(body, 'getCategory')
        if not category:
            category = infer_category(cls)
    results.append((cls, term_name, category))
    print(f"  {cls}: name='{term_name}' category='{category}'")

# Generate registration calls
lines = [
    "    // ========== BATCH 20: PREVIOUSLY UNREGISTERED INLINE CLASSES (49 terms) ==========",
    "    // SOURCE4 astrophysical systems, MUGE compressed/resonance, BSM particle physics",
    "    // Reactor efficiency, universal field components, and unified field terms",
]
for cls, term_name, category in results:
    lines.append(f'    core.registerPhysicsTerm("{term_name}", std::make_unique<{cls}>(), "{category}");')
lines.append("    // ========== END BATCH 20: 49 INLINE TERMS REGISTERED ==========")

output = '\n'.join(lines)

with open('wolfram_extraction/batch20_registrations.txt', 'w', encoding='utf-8') as f:
    f.write(output)

print(f"\nGenerated {len(results)} registration calls.")
print("Output written to wolfram_extraction/batch20_registrations.txt")
print("\nRegistration block preview:")
print(output[:2000])
