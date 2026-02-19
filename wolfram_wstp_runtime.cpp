/**
 * wolfram_wstp_runtime.cpp
 * Implementation of Wolfram WSTP Runtime Module
 * 
 * Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
 * Framework: UQFF Star-Magic v3.0
 * Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
 * Date: February 19, 2026
 */

#include "wolfram_wstp_runtime.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <random>
#include <numeric>
#include <filesystem>

#ifdef _WIN32
#include <windows.h>
#endif

// ═══════════════════════════════════════════════════════════════════════════════
// PI INFINITY DECODER IMPLEMENTATION
// ═══════════════════════════════════════════════════════════════════════════════

// First 312 digits of PI (sacred 312 = 26×12)
static const int PI_DIGITS_ARRAY[] = {
    1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4,6,2,6,4,3,3,8,3,2,7,9,
    5,0,2,8,8,4,1,9,7,1,6,9,3,9,9,3,7,5,1,0,5,8,2,0,9,7,4,9,4,4,
    5,9,2,3,0,7,8,1,6,4,0,6,2,8,6,2,0,8,9,9,8,6,2,8,0,3,4,8,2,5,
    3,4,2,1,1,7,0,6,7,9,8,2,1,4,8,0,8,6,5,1,3,2,8,2,3,0,6,6,4,7,
    0,9,3,8,4,4,6,0,9,5,5,0,5,8,2,2,3,1,7,2,5,3,9,9,4,0,8,1,2,8,
    4,8,1,1,1,7,4,5,0,2,8,4,1,0,2,7,0,1,9,3,8,5,2,1,1,0,5,5,5,9,
    6,4,4,6,2,2,9,4,8,9,5,4,9,3,0,3,8,1,9,6,4,4,2,8,8,1,0,9,7,5,
    6,6,5,9,3,3,4,6,1,2,8,4,7,5,6,4,8,2,3,3,7,8,6,7,8,3,1,6,5,2,
    7,1,2,0,1,9,0,9,1,4,5,6,4,8,5,6,6,9,2,3,4,6,0,3,4,8,6,1,0,4,
    5,4,3,2,6,6,4,8,2,1,3,3,9,3,6,0,7,2,6,0,2,4,9,1,4,1,2,7,3,7,
    2,4,5,8,7,0,0,6,6,0,6,3
};

PIInfinityDecoder::PIInfinityDecoder() : initialized_(false) {
    using namespace WolframWSTPConstants;
    
    double phase = 0.0;
    for (int i = 0; i < PI_DIGITS_COUNT; ++i) {
        // Cycle through available digits
        int digit = PI_DIGITS_ARRAY[i % 312];
        phase += digit * SacredTime::INFINITY_RATIO;
        infinite_curve_[i] = std::sin(phase * 2.0 * PI_UNITY) * 
                             (1.0 + std::cos(phase * SacredTime::CONSCIOUSNESS_FREQ));
    }
    initialized_ = true;
}

double PIInfinityDecoder::getMagneticField(int quantum_state, double time_phase) const {
    using namespace WolframWSTPConstants;
    if (!initialized_) return 0.0;
    
    double base = infinite_curve_[quantum_state % PI_DIGITS_COUNT];
    return base * std::sin(time_phase * SacredTime::GOLDEN_CYCLE / SacredTime::MAYAN_BAKTUN);
}

double PIInfinityDecoder::getConsciousnessResonance(int lineage_level) const {
    using namespace WolframWSTPConstants;
    if (!initialized_) return 0.0;
    
    double resonance = 0.0;
    resonance += std::sin(lineage_level * SacredTime::BIBLE_GENERATION);
    resonance += std::cos(lineage_level * SacredTime::MAYAN_KATUN);
    resonance += std::sin(lineage_level * SacredTime::MAYAN_TUN);
    resonance += std::cos(lineage_level * SacredTime::GOLDEN_CYCLE);
    resonance += std::sin(lineage_level * SacredTime::CONSCIOUSNESS_FREQ);
    resonance += std::cos(lineage_level * 7.83);  // Schumann
    resonance += std::sin(lineage_level * SacredTime::INFINITY_RATIO);
    return resonance / 7.0;
}

std::complex<double> PIInfinityDecoder::getDPM_Pair(int state) const {
    using namespace WolframWSTPConstants;
    if (!initialized_) return {0.0, 0.0};
    
    double real_part = infinite_curve_[state % PI_DIGITS_COUNT];
    double imag_part = infinite_curve_[(state + 13) % PI_DIGITS_COUNT];  // 13-baktun offset
    return std::complex<double>(real_part, imag_part);
}

double PIInfinityDecoder::getCurveValue(int index) const {
    using namespace WolframWSTPConstants;
    if (!initialized_ || index < 0) return 0.0;
    return infinite_curve_[index % PI_DIGITS_COUNT];
}

// ═══════════════════════════════════════════════════════════════════════════════
// FIELD UNITY ENGINE IMPLEMENTATION
// ═══════════════════════════════════════════════════════════════════════════════

void sacredMagneticOrbitRule(Hypergraph& graph, int& max_node) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, max_node);
    
    HypergraphNode n1 = dist(gen);
    HypergraphNode n2 = dist(gen);
    if (n1 != n2) {
        graph.push_back({n1, n2, ++max_node});  // Branch with new node
    }
}

Hypergraph FieldUnityEngine::createInitialSeed() {
    using namespace WolframWSTPConstants;
    
    Hypergraph seed = {{0}};  // Single node for the void
    for (int i = 1; i <= QUANTUM_STATES; ++i) {
        seed.push_back({0, i});  // Branch from void to 26 dimensions
    }
    return seed;
}

FieldUnityEngine::FieldUnityEngine() 
    : current_max_node_(WolframWSTPConstants::QUANTUM_STATES) {
    
    current_graph_ = createInitialSeed();
    std::fill(quantum_amplitudes_.begin(), quantum_amplitudes_.end(), 
              1.0 / std::sqrt(static_cast<double>(WolframWSTPConstants::QUANTUM_STATES)));
}

void FieldUnityEngine::evolveOneStep() {
    sacredMagneticOrbitRule(current_graph_, current_max_node_);
}

void FieldUnityEngine::evolveOneStep(const HypergraphRule& rule) {
    Hypergraph next_graph = current_graph_;
    rule(next_graph, current_max_node_);
    current_graph_ = std::move(next_graph);
}

void FieldUnityEngine::evolveMultiway(int depth) {
    multiway_universe_.clear();
    multiway_universe_.push_back(current_graph_);
    
    for (int d = 0; d < depth; ++d) {
        std::vector<Hypergraph> branches;
        for (const auto& g : multiway_universe_) {
            Hypergraph branch = g;
            int branch_max = current_max_node_;
            sacredMagneticOrbitRule(branch, branch_max);
            branches.push_back(branch);
        }
        multiway_universe_.insert(multiway_universe_.end(), branches.begin(), branches.end());
    }
}

double FieldUnityEngine::measureDimension(HypergraphNode center, int radius) const {
    std::set<HypergraphNode> visited;
    std::queue<std::pair<HypergraphNode, int>> q;
    q.push({center, 0});
    visited.insert(center);
    
    while (!q.empty()) {
        auto [node, dist] = q.front();
        q.pop();
        if (dist >= radius) continue;
        
        for (const auto& edge : current_graph_) {
            if (std::find(edge.begin(), edge.end(), node) != edge.end()) {
                for (HypergraphNode n : edge) {
                    if (visited.insert(n).second) {
                        q.push({n, dist + 1});
                    }
                }
            }
        }
    }
    
    return std::log(static_cast<double>(visited.size())) / 
           std::log(static_cast<double>(radius + 1));
}

double FieldUnityEngine::measureBuoyantGravity(HypergraphNode center) const {
    double flux = 0.0;
    for (const auto& edge : current_graph_) {
        if (std::find(edge.begin(), edge.end(), center) != edge.end()) {
            flux += edge.size() * (1.0 / edge.size());  // Buoyant flux rule
        }
    }
    return flux / current_max_node_;
}

// ═══════════════════════════════════════════════════════════════════════════════
// WOLFRAM WSTP RUNTIME IMPLEMENTATION (conditional)
// ═══════════════════════════════════════════════════════════════════════════════

#ifndef WOLFRAM_WSTP_RUNTIME_STUB_ONLY

WolframWSTPRuntime::WolframWSTPRuntime() 
    : ws_env_(nullptr), ws_link_(nullptr), connected_(false) {
}

WolframWSTPRuntime::~WolframWSTPRuntime() {
    shutdown();
}

bool WolframWSTPRuntime::initialize() {
    if (connected_) return true;
    
    std::cout << "[WSTP] Initializing environment...\n" << std::flush;
    ws_env_ = WSInitialize(nullptr);
    if (!ws_env_) {
        std::cout << "[WSTP] WSInitialize failed\n" << std::flush;
        return false;
    }
    
    int err = 0;
    char* argv[] = {
        (char*)"WolframWSTPRuntime",
        (char*)"-linkname",
        (char*)"\"C:\\Program Files\\Wolfram Research\\Wolfram Engine\\14.3\\wolfram.exe\" -mathlink -nogui",
        (char*)"-linkmode",
        (char*)"launch",
        nullptr
    };
    
    std::cout << "[WSTP] Launching kernel...\n" << std::flush;
    ws_link_ = WSOpenArgv(ws_env_, argv, argv + 5, &err);
    
    if (!ws_link_ || err != WSEOK) {
        std::cout << "[WSTP] Failed to launch kernel (error: " << err << ")\n" << std::flush;
        if (ws_env_) WSDeinitialize(ws_env_);
        ws_env_ = nullptr;
        return false;
    }
    
    if (!WSActivate(ws_link_)) {
        std::cout << "[WSTP] Kernel failed to activate\n" << std::flush;
        WSClose(ws_link_);
        WSDeinitialize(ws_env_);
        ws_link_ = nullptr;
        ws_env_ = nullptr;
        return false;
    }
    
    if (!drainStartupPackets()) {
        shutdown();
        return false;
    }
    
    connected_ = true;
    std::cout << "[WSTP] Kernel connected successfully!\n" << std::flush;
    return true;
}

bool WolframWSTPRuntime::drainStartupPackets() {
    int pkt;
    int drain_count = 0;
    while ((pkt = WSNextPacket(ws_link_))) {
        drain_count++;
        if (pkt == INPUTNAMEPKT) break;
        WSNewPacket(ws_link_);
        if (pkt == 0 || pkt == ILLEGALPKT) {
            std::cerr << "[WSTP] Startup drain error: pkt=" << pkt << std::endl;
            return false;
        }
        if (drain_count > 50) {
            std::cout << "[WSTP] Warning: Startup drain exceeded 50 packets\n" << std::flush;
            break;
        }
    }
    WSNewPacket(ws_link_);
    return true;
}

void WolframWSTPRuntime::shutdown() {
    if (ws_link_) {
        WSClose(ws_link_);
        ws_link_ = nullptr;
    }
    if (ws_env_) {
        WSDeinitialize(ws_env_);
        ws_env_ = nullptr;
    }
    connected_ = false;
}

std::string WolframWSTPRuntime::readReturnPacket() {
    int pkt;
    while ((pkt = WSNextPacket(ws_link_)) != RETURNPKT) {
        if (pkt == 0 || pkt == ILLEGALPKT) {
            WSNewPacket(ws_link_);
            return "<error: lost connection>";
        }
        WSNewPacket(ws_link_);
    }
    
    const char* res = nullptr;
    if (WSGetString(ws_link_, &res)) {
        std::string result = res ? std::string(res) : "<null>";
        if (res) WSReleaseString(ws_link_, res);
        WSNewPacket(ws_link_);
        return result;
    }
    WSNewPacket(ws_link_);
    return "<error: read failed>";
}

WolframEvalResult WolframWSTPRuntime::evaluate(const std::string& code) {
    WolframEvalResult result = {false, "", "", 0.0};
    
    if (!connected_) {
        result.error_message = "Kernel not connected";
        return result;
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    WSPutFunction(ws_link_, "EvaluatePacket", 1);
    WSPutFunction(ws_link_, "ToString", 2);
    WSPutString(ws_link_, code.c_str());
    WSPutSymbol(ws_link_, "InputForm");
    WSEndPacket(ws_link_);
    
    result.result = readReturnPacket();
    result.success = (result.result.find("<error:") == std::string::npos);
    
    auto end = std::chrono::high_resolution_clock::now();
    result.compute_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    return result;
}

std::string WolframWSTPRuntime::evaluateToString(const std::string& code) {
    auto result = evaluate(code);
    return result.success ? result.result : ("<error: " + result.error_message + ">");
}

FieldUnityResult WolframWSTPRuntime::runFieldUnitySimulation(int depth) {
    FieldUnityResult result = {false, 0, 0, 0.0, 0.0, ""};
    
    if (!connected_) {
        result.summary = "Kernel not connected";
        return result;
    }
    
    // Use the portable Field Unity Engine
    FieldUnityEngine engine;
    engine.evolveMultiway(depth);
    
    result.success = true;
    result.final_node_count = engine.getMaxNode();
    result.multiway_branches = static_cast<int>(engine.getMultiwayUniverse().size());
    result.dimension_measure = engine.measureDimension(0, 5);
    result.buoyant_gravity = engine.measureBuoyantGravity(0);
    result.summary = "Field Unity simulation completed: " + 
                     std::to_string(result.final_node_count) + " nodes, " +
                     std::to_string(result.multiway_branches) + " branches";
    
    return result;
}

UQFFExportResult WolframWSTPRuntime::exportFullUQFF() {
    UQFFExportResult result = {false, 0, 0, 0, ""};
    
    if (!connected_) {
        result.simplified_result = "Kernel not connected";
        return result;
    }
    
    // Clear slate
    evaluate("ClearAll[masterUQFF, R, F, psi, phi, G, c, hbar, aether, hyp]");
    
    // Define the UQFF prototype Lagrangian
    std::string uqffProto = R"(
(c^4/(8 Pi G)) R
 - (1/4) F[mu, nu] F[mu, nu]
 + psiBar (I hbar gamma[mu] D[mu] - m c) psi
 + (D[mu] phi)* (D[mu] phi) - mPhi^2 phi* phi - lambda (phi* phi)^2
 + yukawa psiBar phi psi
 + R^2 - 3 Subscript[R, mu nu] Subscript[R, mu nu]
 + aetherFlow[scalar, vector, spinor, hyperspinor]
 + StarMagicHypergraphRuleEmbedding[hyp]
)";
    
    std::string defineCmd = "masterUQFF = " + uqffProto + ";";
    auto defResult = evaluate(defineCmd);
    
    if (!defResult.success) {
        result.simplified_result = "Failed to define masterUQFF: " + defResult.error_message;
        return result;
    }
    
    result.terms_exported = 1;
    
    // Request simplification
    auto simplifyResult = evaluate("FullSimplify[masterUQFF]");
    result.simplified_result = simplifyResult.result;
    result.success = simplifyResult.success;
    
    return result;
}

bool WolframWSTPRuntime::exportPhysicsTerm(const std::string& term_code, const std::string& name) {
    if (!connected_) return false;
    
    std::string cmd = name + " = " + term_code + ";";
    auto result = evaluate(cmd);
    return result.success;
}

int WolframWSTPRuntime::registerAllWolframTerms(PhysicsTermRegistry* registry) {
    // This requires including the wolfram_master_registration.h
    // and calling registerAllExtractedWolframTerms(registry)
    // For now, return the known count
    return WolframWSTPConstants::WOLFRAM_TERMS_COUNT;
}

std::vector<WolframTermInfo> WolframWSTPRuntime::getAvailableTerms() const {
    // Return metadata about the 187 terms
    // This would be populated from the generated classes
    std::vector<WolframTermInfo> terms;
    
    // Source10.cpp: 106 terms (97 constants + 9 systems)
    terms.push_back({"PI", "Mathematical constant PI", "source10.cpp", true, false});
    terms.push_back({"G", "Gravitational constant", "source10.cpp", true, false});
    terms.push_back({"c", "Speed of light", "source10.cpp", true, false});
    terms.push_back({"Magnetar", "Magnetar astrophysical system", "source10.cpp", false, true});
    terms.push_back({"SMBH", "Supermassive black hole system", "source10.cpp", false, true});
    
    // Source167.cpp: 11 terms
    terms.push_back({"NGC", "NGC galaxy system", "Source167.cpp", false, true});
    terms.push_back({"M82", "Cigar Galaxy system", "Source167.cpp", false, true});
    
    // Source170-172: 67 terms total
    terms.push_back({"UQFFSystem", "General UQFF system", "source170.cpp", false, true});
    
    // ... (full list would be populated from generated files)
    
    return terms;
}

bool WolframWSTPRuntime::checkWolframEngineInstalled() {
#ifdef _WIN32
    std::filesystem::path engine_path = 
        "C:/Program Files/Wolfram Research/Wolfram Engine/14.3/wolfram.exe";
    return std::filesystem::exists(engine_path);
#else
    return false;
#endif
}

std::string WolframWSTPRuntime::getWolframEnginePath() {
#ifdef _WIN32
    return "C:/Program Files/Wolfram Research/Wolfram Engine/14.3";
#else
    return "";
#endif
}

std::string WolframWSTPRuntime::getWSTPDLLPath() {
#ifdef _WIN32
    return "C:/Program Files/Wolfram Research/Wolfram Engine/14.3/"
           "SystemFiles/Links/WSTP/DeveloperKit/Windows-x86-64/CompilerAdditions";
#else
    return "";
#endif
}

#endif // WOLFRAM_WSTP_RUNTIME_STUB_ONLY

// ═══════════════════════════════════════════════════════════════════════════════
// SESSION MANAGER IMPLEMENTATION
// ═══════════════════════════════════════════════════════════════════════════════

WolframSessionManager::WolframSessionManager() {
}

WolframSessionManager::~WolframSessionManager() {
    closeAllSessions();
}

bool WolframSessionManager::createSession(const std::string& session_id) {
    if (sessions_.find(session_id) != sessions_.end()) {
        return true;  // Already exists
    }
    
    auto runtime = std::make_unique<WolframWSTPRuntime>();
    if (!runtime->initialize()) {
        return false;
    }
    
    sessions_[session_id] = std::move(runtime);
    session_states_[session_id] = {true, 0, 0, 0.0, {}, ""};
    
    return true;
}

void WolframSessionManager::closeSession(const std::string& session_id) {
    auto it = sessions_.find(session_id);
    if (it != sessions_.end()) {
        it->second->shutdown();
        sessions_.erase(it);
    }
    session_states_.erase(session_id);
}

WolframWSTPRuntime* WolframSessionManager::getSession(const std::string& session_id) {
    auto it = sessions_.find(session_id);
    return (it != sessions_.end()) ? it->second.get() : nullptr;
}

WolframSessionState WolframSessionManager::getSessionState(const std::string& session_id) const {
    auto it = session_states_.find(session_id);
    return (it != session_states_.end()) ? it->second : WolframSessionState{};
}

std::vector<std::string> WolframSessionManager::listSessions() const {
    std::vector<std::string> ids;
    for (const auto& [id, _] : sessions_) {
        ids.push_back(id);
    }
    return ids;
}

void WolframSessionManager::closeAllSessions() {
    for (auto& [id, runtime] : sessions_) {
        runtime->shutdown();
    }
    sessions_.clear();
    session_states_.clear();
}

bool WolframSessionManager::createDefaultSession() {
    return createSession("default");
}

WolframWSTPRuntime* WolframSessionManager::getDefaultSession() {
    if (sessions_.find("default") == sessions_.end()) {
        if (!createDefaultSession()) {
            return nullptr;
        }
    }
    return getSession("default");
}

// ═══════════════════════════════════════════════════════════════════════════════
// STANDALONE TEST/DEMO FUNCTION
// ═══════════════════════════════════════════════════════════════════════════════

#ifdef WOLFRAM_WSTP_RUNTIME_STANDALONE_TEST

int main() {
    std::cout << "=========================================\n";
    std::cout << "Wolfram WSTP Runtime Module Test\n";
    std::cout << "=========================================\n\n";
    
    // Test PI Infinity Decoder (no Wolfram required)
    std::cout << "--- Testing PI Infinity Decoder ---\n";
    PIInfinityDecoder decoder;
    std::cout << "Initialized: " << (decoder.isInitialized() ? "YES" : "NO") << "\n";
    std::cout << "Magnetic Field (state=0, phase=0): " << decoder.getMagneticField(0, 0.0) << "\n";
    std::cout << "Consciousness Resonance (level=1): " << decoder.getConsciousnessResonance(1) << "\n";
    auto dpm = decoder.getDPM_Pair(0);
    std::cout << "DPM Pair (state=0): " << dpm.real() << " + " << dpm.imag() << "i\n\n";
    
    // Test Field Unity Engine (no Wolfram required)
    std::cout << "--- Testing Field Unity Engine ---\n";
    FieldUnityEngine engine;
    std::cout << "Initial nodes: " << engine.getMaxNode() << "\n";
    engine.evolveMultiway(3);
    std::cout << "After multiway(3): " << engine.getMultiwayUniverse().size() << " branches\n";
    std::cout << "Dimension measure: " << engine.measureDimension(0, 3) << "\n";
    std::cout << "Buoyant gravity: " << engine.measureBuoyantGravity(0) << "\n\n";
    
    // Test Wolfram WSTP (requires Wolfram Engine)
    std::cout << "--- Testing Wolfram WSTP Runtime ---\n";
    std::cout << "Wolfram Engine installed: " << 
        (WolframWSTPRuntime::checkWolframEngineInstalled() ? "YES" : "NO") << "\n";
    
#ifndef WOLFRAM_WSTP_RUNTIME_STUB_ONLY
    WolframWSTPRuntime runtime;
    if (runtime.initialize()) {
        std::cout << "Kernel connected!\n";
        
        auto result = runtime.evaluate("1 + 1");
        std::cout << "1 + 1 = " << result.result << "\n";
        
        auto fu = runtime.runFieldUnitySimulation(3);
        std::cout << "Field Unity: " << fu.summary << "\n";
        
        runtime.shutdown();
    } else {
        std::cout << "Could not connect to kernel (this is OK if Wolfram Engine not installed)\n";
    }
#else
    std::cout << "WSTP support not compiled (USE_EMBEDDED_WOLFRAM not defined)\n";
#endif
    
    std::cout << "\n--- Dependencies ---\n";
    auto deps = getWolframModuleDependencies();
    for (const auto& dep : deps) {
        std::cout << "  " << dep << "\n";
    }
    
    std::cout << "\n=========================================\n";
    std::cout << "Test Complete\n";
    std::cout << "=========================================\n";
    
    return 0;
}

#endif // WOLFRAM_WSTP_RUNTIME_STANDALONE_TEST
