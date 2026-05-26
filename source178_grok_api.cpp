// source178_grok_api.cpp
// Grok AI Integration for UQFF Error Diagnostics and Code Assistance
// xAI API Implementation using Python wrapper (GrokAPI.py)
// Integration Date: 2025-11-22 18:50 (Phase 23: AI Integration)
// Updated: 2026-02-13 - Replaced Qt networking with Python subprocess (console app compatible)

#include <iostream>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <array>

#ifdef _WIN32
#define popen _popen
#define pclose _pclose
#endif

// Helper function to execute command and capture output
std::string exec_command(const std::string& cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
    
    if (!pipe) {
        return "[ERROR] Failed to execute command";
    }
    
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    
    return result;
}

// Function to call xAI Grok API via Python wrapper
// Returns the AI response as string; empty on error
// Improved (Activation Plan): delegates key management to Python (config file or env),
// tries multiple python executables, and gives clearer diagnostics.
std::string callGrokAPI(const std::string& prompt) {
    // Build Python command (escape quotes in prompt)
    std::string escaped_prompt = prompt;
    size_t pos = 0;
    while ((pos = escaped_prompt.find("\"", pos)) != std::string::npos) {
        escaped_prompt.replace(pos, 1, "\\\"");
        pos += 2;
    }

    // Try several python launchers (python, python3, py) + also look next to the exe if possible.
    // This makes the feature far more reliable when the user double-clicks the .exe.
    const std::vector<std::string> pythonCandidates = {
        "python",
        "python3",
        "py -3",
        "py"
    };

    std::string jsonOutput;
    std::string usedPython;

    for (const auto& py : pythonCandidates) {
        std::string cmd = py + " GrokAPI.py \"" + escaped_prompt + "\" \"grok-4-1-fast-reasoning\"";
        std::cout << "[Grok API] Trying: " << cmd << "\n" << std::flush;
        jsonOutput = exec_command(cmd);
        usedPython = py;

        if (!jsonOutput.empty() && jsonOutput != "[ERROR] Failed to execute command") {
            break;
        }
    }

    if (jsonOutput.empty() || jsonOutput == "[ERROR] Failed to execute command") {
        std::cout << "[Grok API] ERROR: Could not execute GrokAPI.py with any Python launcher.\n";
        std::cout << "[Grok API] Troubleshooting:\n";
        std::cout << "  1. Install Python 3.9+ and ensure it is in PATH\n";
        std::cout << "  2. pip install requests\n";
        std::cout << "  3. GrokAPI.py and APIKeyManager.py must be next to MAIN_1_CoAnQi.exe (CMake copies them automatically)\n";
        std::cout << "  4. From a command prompt in the exe directory, test manually:\n";
        std::cout << "        python GrokAPI.py \"Hello\" \"grok-4-1-fast-reasoning\"\n" << std::flush;
        return "";
    }

    // If Python side returned a structured error about missing key, surface it nicely
    if (jsonOutput.find("XAI_API_KEY") != std::string::npos && jsonOutput.find("not set") != std::string::npos) {
        std::cout << "[Grok API] Key not found (neither env var nor grok_api_config.json).\n";
        std::cout << "             Use the menu 'Configure Grok API Key' option first.\n";
    }

    // Parse JSON response (simple manual parsing to avoid Qt dependency)
    // Format: {"success": true/false, "response": "...", "error": "..."}
    
    size_t successPos = jsonOutput.find("\"success\"");
    size_t responsePos = jsonOutput.find("\"response\"");
    size_t errorPos = jsonOutput.find("\"error\"");
    
    bool success = false;
    if (successPos != std::string::npos) {
        size_t truePos = jsonOutput.find("true", successPos);
        size_t falsePos = jsonOutput.find("false", successPos);
        
        if (truePos != std::string::npos && (falsePos == std::string::npos || truePos < falsePos)) {
            size_t commaPos = jsonOutput.find(",", successPos);
            if (commaPos != std::string::npos && truePos < commaPos) {
                success = true;
            }
        }
    }
    
    std::string responseContent;
    if (success && responsePos != std::string::npos) {
        // Extract response content between "response": "..." and next "
        size_t contentStart = jsonOutput.find("\"", responsePos + 10);  // Skip "response"
        if (contentStart != std::string::npos) {
            contentStart++;  // Skip opening quote
            size_t contentEnd = contentStart;
            
            // Find closing quote (handling escaped quotes)
            while (contentEnd < jsonOutput.length()) {
                if (jsonOutput[contentEnd] == '\"' && (contentEnd == 0 || jsonOutput[contentEnd-1] != '\\')) {
                    break;
                }
                contentEnd++;
            }
            
            if (contentEnd < jsonOutput.length()) {
                responseContent = jsonOutput.substr(contentStart, contentEnd - contentStart);
                
                // Unescape common sequences
                size_t pos2 = 0;
                while ((pos2 = responseContent.find("\\n", pos2)) != std::string::npos) {
                    responseContent.replace(pos2, 2, "\n");
                    pos2++;
                }
                pos2 = 0;
                while ((pos2 = responseContent.find("\\\"", pos2)) != std::string::npos) {
                    responseContent.replace(pos2, 2, "\"");
                    pos2++;
                }
            }
        }
        
        std::cout << "[Grok API] Response received successfully.\n" << std::flush;
    } else if (errorPos != std::string::npos) {
        // Extract error message
        size_t errorStart = jsonOutput.find("\"", errorPos + 7);
        if (errorStart != std::string::npos) {
            errorStart++;
            size_t errorEnd = jsonOutput.find("\"", errorStart);
            if (errorEnd != std::string::npos) {
                std::string errorMsg = jsonOutput.substr(errorStart, errorEnd - errorStart);
                std::cout << "[Grok API] ERROR: " << errorMsg << "\n" << std::flush;
            }
        }
    }
    
    return responseContent;
}

// Convenience wrapper for C++ error diagnostics
std::string diagnoseCompilationError(const std::string& errorMessage, const std::string& sourceFile = "", int lineNumber = -1) {
    std::string prompt = "C++ Compilation Error Diagnostic:\n\n";
    prompt += "Error Message: " + errorMessage + "\n";
    
    if (!sourceFile.empty()) {
        prompt += "Source File: " + sourceFile + "\n";
    }
    
    if (lineNumber > 0) {
        prompt += "Line Number: " + std::to_string(lineNumber) + "\n";
    }
    
    prompt += "\nContext: This is from the UQFF (Unified Quantum Field Framework) project, a C++20 codebase with:\n";
    prompt += "- Compiler: MSVC 14.44+ (Visual Studio 2022)\n";
    prompt += "- Dependencies: Qt6, ANTLR4, SymEngine, Wolfram WSTP\n";
    prompt += "- Physics: 894 PhysicsTerm classes for quantum field calculations\n\n";
    prompt += "Please explain the error and suggest a fix.";
    
    return callGrokAPI(prompt);
}

// Wrapper for physics equation assistance
std::string explainPhysicsEquation(const std::string& equationName, const std::string& equationCode = "") {
    std::string prompt = "Physics Equation Explanation:\n\n";
    prompt += "Equation: " + equationName + "\n\n";
    
    if (!equationCode.empty()) {
        prompt += "C++ Implementation:\n" + equationCode + "\n\n";
    }
    
    prompt += "Context: This is from the UQFF (Unified Quantum Field Framework) - a quantum field unity theory.\n";
    prompt += "Please explain the physics behind this equation and its role in field unification.";
    
    return callGrokAPI(prompt);
}

// Wrapper for code review and optimization suggestions
std::string reviewPhysicsCode(const std::string& codeSnippet, const std::string& concernedAspect = "performance") {
    std::string prompt = "C++ Code Review for Physics Simulation:\n\n";
    prompt += "Code:\n" + codeSnippet + "\n\n";
    prompt += "Review Focus: " + concernedAspect + "\n\n";
    prompt += "Context: UQFF quantum field calculations with C++20, MSVC optimization (/Os /GL /LTCG).\n";
    prompt += "Please review for correctness, performance, and suggest improvements.";
    
    return callGrokAPI(prompt);
}

// Test function to verify Grok API connectivity
void testGrokAPI() {
    std::cout << "\n=== Testing Grok API Connectivity ===\n" << std::flush;
    
    std::string testPrompt = "Hello Grok! Please confirm you can assist with C++ physics simulations. "
                            "Respond with a brief message acknowledging your capabilities.";
    
    std::string response = callGrokAPI(testPrompt);
    
    if (!response.empty()) {
        std::cout << "\n[Grok API Test] SUCCESS!\n" << std::flush;
        std::cout << "Grok Response:\n" << response << "\n" << std::flush;
    } else {
        std::cout << "\n[Grok API Test] FAILED - No response received.\n" << std::flush;
        std::cout << "Check:\n";
        std::cout << "1. XAI_API_KEY environment variable is set\n";
        std::cout << "2. Internet connection is active\n";
        std::cout << "3. API key is valid (get free key at https://x.ai/api)\n";
        std::cout << "4. Python 3 is installed with 'requests' library (pip install requests)\n";
        std::cout << "5. GrokAPI.py is in the same directory as MAIN_1_CoAnQi.exe\n" << std::flush;
    }
    
    std::cout << "=====================================\n\n" << std::flush;
}

// ===========================================================================================
// Grok API Configuration & Status Helpers (Activation Plan 2026)
// ===========================================================================================

// Prints a clear one-time status about Grok API key availability at startup or on demand.
// Checks environment variable first (fast path) then falls back to asking Python APIKeyManager.
void printGrokAPIStatus() {
    const char* envKey = std::getenv("XAI_API_KEY");
    bool envHasKey = (envKey && strlen(envKey) > 8); // rough sanity check for real key length

    std::cout << "\n[Grok API Status] ";
    if (envHasKey) {
        std::cout << "ACTIVE (XAI_API_KEY found in environment)\n";
        return;
    }

    // Ask Python side (it checks config file first, then env)
    std::string cmd = "python -c \"from APIKeyManager import get_api_key_status; print(get_api_key_status())\" 2>NUL";
#ifdef _WIN32
    // Windows: suppress stderr noise if module missing
#endif
    std::string status = exec_command(cmd);

    // Clean the output
    // Remove trailing newlines / carriage returns
    while (!status.empty() && (status.back() == '\n' || status.back() == '\r')) {
        status.pop_back();
    }

    if (status.find("Not configured") != std::string::npos || status.empty()) {
        std::cout << "NOT CONFIGURED\n";
        std::cout << "  -> Use menu option 'Configure Grok API Key' (or set XAI_API_KEY env var)\n";
        std::cout << "  -> Get a free key at https://x.ai/api\n";
    } else {
        std::cout << status << "\n";
    }
}

// Interactive configuration helper. Persists to grok_api_config.json via APIKeyManager
// and sets the key for the current process so the session can use Grok immediately.
void configureGrokAPIKey() {
    std::cout << "\n=== Configure Grok API Key ===\n";
    std::cout << "The key will be saved to grok_api_config.json (preferred) and also set for this session.\n";
    std::cout << "Get your free API key at: https://x.ai/api\n";
    std::cout << "\nPaste your xAI key (starts with xai-): ";

    std::string apiKey;
    std::getline(std::cin, apiKey);

    // Trim whitespace
    apiKey.erase(0, apiKey.find_first_not_of(" \t\r\n"));
    apiKey.erase(apiKey.find_last_not_of(" \t\r\n") + 1);

    if (apiKey.empty()) {
        std::cout << "No API key entered. Configuration cancelled.\n";
        return;
    }

    // 1. Set for current process (so this run of the program can use it immediately)
#ifdef _WIN32
    _putenv_s("XAI_API_KEY", apiKey.c_str());
#else
    setenv("XAI_API_KEY", apiKey.c_str(), 1);
#endif

    // 2. Persist via APIKeyManager.py (writes grok_api_config.json)
    // We use a Python one-liner to avoid needing a dedicated --set CLI in GrokAPI.py
    std::string persistCmd = "python -c \"from APIKeyManager import set_xai_api_key; import sys; sys.exit(0 if set_xai_api_key('" + apiKey + "') else 1)\"";
    std::string persistResult = exec_command(persistCmd);

    std::cout << "\n? Grok API key configured for this session.\n";

    if (persistResult.find("Error") == std::string::npos && !persistResult.empty()) {
        std::cout << "? Also saved to grok_api_config.json (will be used by future runs + GrokAPI.py)\n";
    } else {
        std::cout << "  (Note: Could not persist to config file this time - key is still active for this session via env var)\n";
        std::cout << "  Run this in the same directory later: python -c \"from APIKeyManager import set_xai_api_key; set_xai_api_key('your_key')\"\n";
    }

    std::cout << "You can now select 'Test Grok AI Integration' from the menu.\n\n";
}
