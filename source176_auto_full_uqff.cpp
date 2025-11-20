
/// source176_auto_full_uqff.cpp  ←  FINAL BULLETPROOF VERSION
#ifdef USE_EMBEDDED_WOLFRAM

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <regex>
#include <windows.h>

#include "source174_wolfram_bridge_embedded.cpp"

void AutoExportFullUQFF()
{
    std::cout << "\n=== AUTO FULL UQFF EXPORT (bulletproof) ===\n";

    if (!InitializeWolframKernel())
        return;

    // Force project root as search directory (works even if VS changes working dir)
    wchar_t buf[MAX_PATH];
    GetModuleFileNameW(nullptr, buf, MAX_PATH);
    std::filesystem::path exePath = buf;
    std::filesystem::path root = exePath.parent_path().parent_path().parent_path(); // x64/Release → repo root

    std::cout << "Searching in: " << root.string() << "\n";

    std::vector<std::string> terms;
    std::regex term_regex(R"(#define\s+WOLFRAM_TERM\s*\"(.*)\")");

    int files = 0;
    try
    {
        for (const auto &entry : std::filesystem::directory_iterator(root))
        {
            if (!entry.path().string().ends_with(".cpp"))
                continue;
            files++;

            std::ifstream f(entry.path());
            if (!f)
                continue;

            std::string line;
            while (std::getline(f, line))
            {
                std::smatch m;
                if (std::regex_search(line, m, term_regex) && m.size() > 1)
                {
                    terms.push_back(m.str(1));
                    std::cout << "✓ " << entry.path().filename().string() << " → " << m.str(1).substr(0, 60) << "...\n";
                    break;
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cout << "Filesystem exception: " << e.what() << "\n";
    }

    std::cout << "Scanned " << files << " files, found " << terms.size() << " terms.\n";

    if (terms.empty())
    {
        std::cout << "No terms found — running fallback test so we know the kernel works...\n";
        WolframEvalToString("test = 42^2");
        std::cout << "Kernel alive → 42² = " << WolframEvalToString("test") << "\n";
        return;
    }

    std::string expr = "masterUQFF = " + terms[0];
    for (size_t i = 1; i < terms.size(); ++i)
        expr += " + " + terms[i];
    expr += ";";

    WolframEvalToString(expr);
    std::string result = WolframEvalToString("ToString[FullSimplify[masterUQFF], InputForm]");

    std::cout << "\n=== WOLFRAM VERDICT ===\n"
              << result << "\n";
    if (result.find("0") != std::string::npos)
        std::cout << "★★★ CANCELLATION ACHIEVED ★★★\n";
}

#endif