
/// source176_auto_full_uqff.cpp  ←  ENHANCED VERSION (January 26, 2026)
/// Enhanced scanning: More files, more patterns, better discovery for self-expansion
/// Preserves all mathematical methods and variable equation solutions
#ifdef USE_EMBEDDED_WOLFRAM

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <windows.h>
#include <iomanip>
#include <algorithm>
#include <set>

#include "source174_wolfram_bridge_embedded.cpp"

void AutoExportFullUQFF()
{
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "AUTO FULL UQFF EXPORT - ENHANCED (January 26, 2026)" << std::endl;
    std::cout << "Self-Expanding Term Discovery for Modular Physics" << std::endl;
    std::cout << std::string(70, '=') << "\n" << std::flush;

    if (!InitializeWolframKernel())
    {
        std::cout << "[ERROR] InitializeWolframKernel() returned false\n" << std::flush;
        return;
    }

    std::cout << "[OK] Kernel initialized successfully\n" << std::flush;

    // Force project root as search directory
    wchar_t buf[MAX_PATH];
    GetModuleFileNameW(nullptr, buf, MAX_PATH);
    std::filesystem::path exePath = buf;
    std::filesystem::path root = exePath.parent_path();
    
    if (root.filename() == "Release" || root.filename() == "Debug") {
        root = root.parent_path();
        if (root.filename().string().find("build") != std::string::npos) {
            root = root.parent_path();
        }
    }
    
    if (!std::filesystem::exists(root / "source1.cpp")) {
        std::cout << "[WARN] source1.cpp not found in " << root.string() << "\n" << std::flush;
        root = exePath.parent_path();
    }

    std::cout << "Scanning: " << root.string() << "\n" << std::flush;

    // ==== PHASE 1: Scan for WOLFRAM_TERM macros ====
    std::cout << "\n--- PHASE 1: Scanning for WOLFRAM_TERM Macros ---\n" << std::flush;
    
    std::vector<std::string> wolframTerms;
    std::regex term_regex(R"(#define\s+WOLFRAM_TERM\s*\"(.*)\")");
    
    int files_scanned = 0;
    int wolfram_terms_found = 0;

    try
    {
        for (const auto &entry : std::filesystem::directory_iterator(root))
        {
            std::string path_str = entry.path().string();
            if (path_str.length() < 4 || path_str.substr(path_str.length() - 4) != ".cpp")
                continue;
            
            files_scanned++;

            std::ifstream f(entry.path());
            if (!f) continue;

            std::string line;
            int lines_read = 0;
            while (std::getline(f, line) && lines_read < 2000)  // Increased limit for larger files
            {
                lines_read++;
                std::smatch m;
                if (std::regex_search(line, m, term_regex) && m.size() > 1)
                {
                    wolframTerms.push_back(m.str(1));
                    wolfram_terms_found++;
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cout << "[WARN] Phase 1 exception: " << e.what() << "\n" << std::flush;
    }

    std::cout << "Files scanned: " << files_scanned << "\n";
    std::cout << "WOLFRAM_TERM macros found: " << wolfram_terms_found << "\n" << std::flush;

    // ==== PHASE 2: Scan for PhysicsTerm class names ====
    std::cout << "\n--- PHASE 2: Scanning for PhysicsTerm Classes ---\n" << std::flush;
    
    std::set<std::string> physicsTermClasses;
    std::regex class_regex(R"(class\s+(\w+Term)\s*:\s*public\s+PhysicsTerm)");
    std::regex class_regex2(R"(class\s+(\w+Module)\s*:\s*public\s+UQFFModule)");
    
    int physics_classes_found = 0;
    files_scanned = 0;

    try
    {
        for (const auto &entry : std::filesystem::directory_iterator(root))
        {
            std::string path_str = entry.path().string();
            if (path_str.length() < 4 || path_str.substr(path_str.length() - 4) != ".cpp")
                continue;
            
            files_scanned++;

            std::ifstream f(entry.path());
            if (!f) continue;

            std::string line;
            while (std::getline(f, line))
            {
                std::smatch m;
                if (std::regex_search(line, m, class_regex) && m.size() > 1)
                {
                    physicsTermClasses.insert(m.str(1));
                    physics_classes_found++;
                }
                if (std::regex_search(line, m, class_regex2) && m.size() > 1)
                {
                    physicsTermClasses.insert(m.str(1));
                    physics_classes_found++;
                }
            }
        }
        
        // Also scan Core/Modules subdirectory
        std::filesystem::path modulesDir = root / "Core" / "Modules";
        if (std::filesystem::exists(modulesDir))
        {
            for (const auto &entry : std::filesystem::directory_iterator(modulesDir))
            {
                std::string path_str = entry.path().string();
                if (path_str.length() < 4 || path_str.substr(path_str.length() - 4) != ".cpp")
                    continue;

                std::ifstream f(entry.path());
                if (!f) continue;

                std::string line;
                while (std::getline(f, line))
                {
                    std::smatch m;
                    if (std::regex_search(line, m, class_regex) && m.size() > 1)
                    {
                        physicsTermClasses.insert(m.str(1));
                        physics_classes_found++;
                    }
                    if (std::regex_search(line, m, class_regex2) && m.size() > 1)
                    {
                        physicsTermClasses.insert(m.str(1));
                        physics_classes_found++;
                    }
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cout << "[WARN] Phase 2 exception: " << e.what() << "\n" << std::flush;
    }

    std::cout << "PhysicsTerm/UQFFModule classes found: " << physicsTermClasses.size() << "\n" << std::flush;

    // ==== PHASE 3: Build and send to Wolfram ====
    std::cout << "\n--- PHASE 3: Building Wolfram Expressions ---\n" << std::flush;
    
    // Send WOLFRAM_TERM expressions
    if (!wolframTerms.empty())
    {
        std::cout << "Building master expression from " << wolframTerms.size() << " WOLFRAM_TERM macros...\n" << std::flush;
        
        std::ostringstream expr_builder;
        expr_builder << "masterUQFF = " << wolframTerms[0];
        
        for (size_t i = 1; i < wolframTerms.size(); ++i)
        {
            if (i % 50 == 0)
            {
                std::cout << "  Added " << i << " / " << wolframTerms.size() << " terms...\n" << std::flush;
            }
            expr_builder << " + " << wolframTerms[i];
        }
        expr_builder << ";";

        std::string expr = expr_builder.str();
        size_t expr_size_kb = expr.size() / 1024;
        std::cout << "Expression size: " << expr_size_kb << " KB\n" << std::flush;

        std::cout << "Sending to Wolfram kernel...\n" << std::flush;
        WolframEvalToString(expr);
        
        std::cout << "Running FullSimplify (may take time for large expressions)...\n" << std::flush;
        std::string result = WolframEvalToString("ToString[FullSimplify[masterUQFF], InputForm]");
        
        std::cout << "\n=== WOLFRAM SIMPLIFICATION RESULT ===\n" << result << "\n" << std::flush;
        
        if (result.find("0") != std::string::npos && result.length() < 10)
            std::cout << "★★★ PERFECT CANCELLATION ACHIEVED ★★★\n" << std::flush;
    }
    else
    {
        std::cout << "No WOLFRAM_TERM macros found.\n" << std::flush;
    }

    // Send class names list for discovery
    if (!physicsTermClasses.empty())
    {
        std::cout << "\nSending " << physicsTermClasses.size() << " class names to Wolfram...\n" << std::flush;
        
        std::ostringstream classList;
        classList << "uqffClasses = {";
        bool first = true;
        for (const auto& cls : physicsTermClasses)
        {
            if (!first) classList << ", ";
            first = false;
            classList << "\"" << cls << "\"";
        }
        classList << "};";
        
        WolframEvalToString(classList.str());
        
        std::ostringstream countExpr;
        countExpr << "uqffClassCount = " << physicsTermClasses.size() << ";";
        WolframEvalToString(countExpr.str());
    }

    // ==== PHASE 4: Summary ====
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "EXPORT COMPLETE - SELF-EXPANSION DATA AVAILABLE" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::cout << "\nWolfram Variables Created:" << std::endl;
    std::cout << "  masterUQFF      - Sum of " << wolframTerms.size() << " WOLFRAM_TERM expressions" << std::endl;
    std::cout << "  uqffClasses     - List of " << physicsTermClasses.size() << " PhysicsTerm class names" << std::endl;
    std::cout << "  uqffClassCount  - Total class count" << std::endl;
    
    std::cout << "\nFor self-expansion discovery:" << std::endl;
    std::cout << "  - Analyze masterUQFF for new mathematical patterns" << std::endl;
    std::cout << "  - Use uqffClasses to map term relationships" << std::endl;
    std::cout << "  - Run: Variables[masterUQFF] for all variables" << std::endl;
    std::cout << "  - Run: CoefficientRules[masterUQFF] for structure analysis" << std::endl;
    
    std::cout << "\n" << std::string(70, '=') << "\n" << std::flush;
}

#endif
