
// source176_auto_full_uqff.cpp  ←  FIXED VERSION – copy-paste over the old one
#ifdef USE_EMBEDDED_WOLFRAM

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <regex>
#include <windows.h> // for GetModuleFileNameW

#include "source174_wolfram_bridge_embedded.cpp"

std::string GetExecutableDir()
{
    wchar_t path[MAX_PATH] = {0};
    GetModuleFileNameW(nullptr, path, MAX_PATH);
    std::wstring wpath(path);
    std::filesystem::path p = std::filesystem::path(wpath).parent_path();
    return p.string();
}

void AutoExportFullUQFF()
{
    std::cout << "\n=== AUTO FULL 26D UQFF EXPORT (debug build) ===\n";

    if (!InitializeWolframKernel())
    {
        std::cout << "KERNEL FAILED TO START – cannot continue.\n";
        return;
    }

    std::string root = GetExecutableDir();
    std::cout << "Searching for source files in: " << root << "\n";

    std::vector<std::string> terms;
    std::regex term_regex(R"(#define\s+WOLFRAM_TERM\s*\"(.*)\")", std::regex::icase);

    int files_scanned = 0;

    for (const auto &entry : std::filesystem::directory_iterator(root))
    {
        std::string path = entry.path().string();
        if (!path.ends_with(".cpp"))
            continue;

        std::ifstream file(path);
        if (!file.is_open())
            continue;

        std::string line;
        bool found = false;
        while (std::getline(file, line) && !found)
        {
            std::smatch m;
            if (std::regex_search(line, m, term_regex) && m.size() > 1)
            {
                terms.push_back(m.str(1));
                std::cout << "Found term in " << entry.path().filename().string() << " : " << m.str(1) << "\n";
                found = true;
            }
        }
        files_scanned++;
    }

    std::cout << "Scanned " << files_scanned << " files, found " << terms.size() << " WOLFRAM_TERM definitions.\n";

    if (terms.empty())
    {
        std::cout << "No terms found! Did you run the Python injector? Or are you running from wrong folder?\n";
        std::cout << "Sending a tiny test expression instead so we can prove the kernel is alive...\n";
        WolframEvalToString("testUQFF = 1 + 1;");
        std::cout << "Test result: " << WolframEvalToString("testUQFF") << "\n";
        return;
    }

    std::string fullExpr = "masterUQFF = " + terms[0];
    for (size_t i = 1; i < terms.size(); ++i)
        fullExpr += " + " + terms[i];
    fullExpr += ";";

    std::cout << "Sending " << fullExpr.size() << " bytes to kernel...\n";
    WolframEvalToString(fullExpr);

    std::cout << "Simplifying...\n";
    std::string result = WolframEvalToString("ToString[FullSimplify[masterUQFF], InputForm]");

    std::cout << "\n=== WOLFRAM VERDICT ===\n"
              << result << "\n";

    if (result == "0")
        std::cout << "★★★ UNIFICATION ACHIEVED ★★★\n";
}

#endif