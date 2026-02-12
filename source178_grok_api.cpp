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
std::string callGrokAPI(const std::string& prompt) {
    // Check if XAI_API_KEY is set
    const char* apiKey = std::getenv("XAI_API_KEY");
    if (!apiKey || strlen(apiKey) == 0) {
        std::cout << "[Grok API] WARNING: XAI_API_KEY environment variable not found.\n";
        std::cout << "[Grok API] Set with PowerShell: $env:XAI_API_KEY=\"your_key_here\"\n";
        std::cout << "[Grok API] Get your free API key at: https://x.ai/api\n" << std::flush;
        return "";
    }

    // Build Python command (escape quotes in prompt)
    std::string escaped_prompt = prompt;
    size_t pos = 0;
    while ((pos = escaped_prompt.find("\"", pos)) != std::string::npos) {
        escaped_prompt.replace(pos, 1, "\\\"");
        pos += 2;
    }

    std::string cmd = "python GrokAPI.py \"" + escaped_prompt + "\" \"grok-4-1-fast-reasoning\"";
    
    std::cout << "[Grok API] Calling Python wrapper...\n" << std::flush;
    
    // Execute Python script and capture JSON output
    std::string jsonOutput = exec_command(cmd);
    
    if (jsonOutput.empty() || jsonOutput == "[ERROR] Failed to execute command") {
        std::cout << "[Grok API] ERROR: Failed to execute GrokAPI.py\n";
        std::cout << "[Grok API] Check:\n";
        std::cout << "  1. Python 3 is installed and in PATH\n";
        std::cout << "  2. GrokAPI.py exists in current directory\n";
        std::cout << "  3. requests library installed: pip install requests\n" << std::flush;
        return "";
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
