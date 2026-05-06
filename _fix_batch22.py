content = open('MAIN_1_CoAnQi.cpp', 'r', encoding='utf-8', errors='replace').read()

# 1. Remove BATCH 22 from registerAllPhysicsTerms_SOURCE4
batch22_block = (
    '\n    // ========== BATCH 22: SESSION 137 grok_share_84a767d3 TERMS (6 terms) ==========\n'
    '    // PI Infinity Decoder, Wolfram Field Unity, Sacred Time Phase, Hypergraph Dimension,\n'
    '    // Buoyant Gravity Hypergraph, WSTP Bridge Validation\n'
    '    // Source: grok_share_84a767d3.txt (4310 lines, Nov 2025 WSTP dev history)\n'
    '    // Papers: PAPER_502 (WSTP), PAPER_506 (PI Decoder), PAPER_507 (Hypergraph), PAPER_508 (Sacred Time)\n'
    '    core.registerPhysicsTerm("PIInfinityDecoder_84A767D3",\n'
    '        std::make_unique<PIInfinityDecoderTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    core.registerPhysicsTerm("WolframFieldUnity_Resonance84A767D3",\n'
    '        std::make_unique<WolframFieldUnityTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    core.registerPhysicsTerm("SacredTimePhase_84A767D3",\n'
    '        std::make_unique<SacredTimePhaseTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    core.registerPhysicsTerm("HypergraphDimension_84A767D3",\n'
    '        std::make_unique<HypergraphDimensionTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    core.registerPhysicsTerm("BuoyantGravityHypergraph_84A767D3",\n'
    '        std::make_unique<BuoyantGravityHypergraphTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    core.registerPhysicsTerm("WSTPBridgeValidation_84A767D3",\n'
    '        std::make_unique<WSTPBridgeValidationTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    g_logger.log("Batch 22 complete: 6 Session137 grok_share_84a767d3 terms registered", 1);\n'
    '    // ========== END BATCH 22: 6 SESSION 137 TERMS REGISTERED ==========\n'
)
if batch22_block in content:
    content = content.replace(batch22_block, '\n    // BATCH 22 moved to registerPhysicsTerms_84A767D3Batch() (after class declarations)\n', 1)
    print('BATCH 22 removed from registerAllPhysicsTerms_SOURCE4 OK')
else:
    print('ERROR: batch22 block not found - checking partial match')
    # Show context for debugging
    idx = content.find('PIInfinityDecoder_84A767D3')
    print(f'First occurrence of PIInfinityDecoder_84A767D3 at char offset {idx}')
    print(repr(content[max(0,idx-200):idx+300]))

# 2. Add new late-registration function after WSTPBridgeValidationTerm_84A767D3 class
new_func = (
    '\nvoid registerPhysicsTerms_84A767D3Batch(CalculatorCore& core) {\n'
    '    // BATCH 22 (late): SESSION 137 grok_share_84a767d3 TERMS (6 terms)\n'
    '    // Moved here because class definitions follow their original call site\n'
    '    core.registerPhysicsTerm("PIInfinityDecoder_84A767D3",\n'
    '        std::make_unique<PIInfinityDecoderTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    core.registerPhysicsTerm("WolframFieldUnity_Resonance84A767D3",\n'
    '        std::make_unique<WolframFieldUnityTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    core.registerPhysicsTerm("SacredTimePhase_84A767D3",\n'
    '        std::make_unique<SacredTimePhaseTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    core.registerPhysicsTerm("HypergraphDimension_84A767D3",\n'
    '        std::make_unique<HypergraphDimensionTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    core.registerPhysicsTerm("BuoyantGravityHypergraph_84A767D3",\n'
    '        std::make_unique<BuoyantGravityHypergraphTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    core.registerPhysicsTerm("WSTPBridgeValidation_84A767D3",\n'
    '        std::make_unique<WSTPBridgeValidationTerm_84A767D3>(), "Session137-grok84A767D3");\n'
    '    g_logger.log("Batch 22 (late) complete: 6 grok_share_84a767d3 terms registered", 1);\n'
    '}\n'
    '\n'
)
marker = 'inline void runSession137PhysicsTerms('
if marker in content:
    content = content.replace(marker, new_func + marker, 1)
    print('registerPhysicsTerms_84A767D3Batch function inserted OK')
else:
    print('ERROR: runSession137PhysicsTerms marker not found')

# 3. Add the call after registerAllPhysicsTerms_SOURCE4
old_call = '        registerAllPhysicsTerms_SOURCE4(g_calculatorCore);\n        registerPhysicsTerms_84A767D3Batch(g_calculatorCore);\n'
if old_call in content:
    print('Call already present - no change needed')
else:
    old_call2 = '        registerAllPhysicsTerms_SOURCE4(g_calculatorCore);\n\n        if (regSpan)'
    new_call2 = '        registerAllPhysicsTerms_SOURCE4(g_calculatorCore);\n        registerPhysicsTerms_84A767D3Batch(g_calculatorCore);\n\n        if (regSpan)'
    if old_call2 in content:
        content = content.replace(old_call2, new_call2, 1)
        print('Call to registerPhysicsTerms_84A767D3Batch added OK')
    else:
        print('ERROR: call site not found')

open('MAIN_1_CoAnQi.cpp', 'w', encoding='utf-8').write(content)
print('File written OK')
