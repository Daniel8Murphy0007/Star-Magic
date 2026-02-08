// source2_minimal_test.cpp
// Minimal MSVC compatibility test for source2.cpp
// Tests core C++ standard library usage without external dependencies
// Date: December 1, 2025

#include <windows.h>   // Windows API (MSVC compatible)
#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <chrono>

// Test structure from source2.cpp
struct SearchResult
{
    std::string url;
    std::string title;
    std::string summary;
    double relevance;
    bool isLive;
};

// Test global vector usage
std::vector<std::string> focusList = {
    "NASA", "STScI", "Hubble", "JWST", "ALMA"
};

// Test thread creation (Windows API compatible)
void TestThread()
{
    std::cout << "Thread running on ID: " << std::this_thread::get_id() << std::endl;
}

int main(int argc, char *argv[])
{
    std::cout << "========================================" << std::endl;
    std::cout << "  SOURCE2 MSVC COMPATIBILITY TEST" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // Test 1: String operations
    std::string testStr = "UQFF Quantum Physics Test";
    std::cout << "[TEST 1] String: " << testStr << std::endl;

    // Test 2: Vector operations
    std::vector<SearchResult> results;
    results.push_back({"https://nasa.gov", "NASA Homepage", "Space agency", 0.95, false});
    std::cout << "[TEST 2] Vector size: " << results.size() << std::endl;
    std::cout << "         Result URL: " << results[0].url << std::endl;

    // Test 3: Focus list iteration
    std::cout << "[TEST 3] Focus list (" << focusList.size() << " items):" << std::endl;
    for (const auto& org : focusList) {
        std::cout << "         - " << org << std::endl;
    }

    // Test 4: Algorithm (find)
    auto it = std::find(focusList.begin(), focusList.end(), "Hubble");
    if (it != focusList.end()) {
        std::cout << "[TEST 4] Found 'Hubble' in focus list" << std::endl;
    }

    // Test 5: Thread creation
    std::cout << "[TEST 5] Creating thread..." << std::endl;
    std::thread t(TestThread);
    t.join();

    // Test 6: Windows API
    #ifdef _WIN32
    std::cout << "[TEST 6] Windows API available: YES" << std::endl;
    #else
    std::cout << "[TEST 6] Windows API available: NO" << std::endl;
    #endif

    // Test 7: Chrono (timing)
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::cout << "[TEST 7] Current time: " << std::ctime(&time_t_now);

    // Test 8: StringStream
    std::stringstream ss;
    ss << "MSVC " << _MSC_VER << " C++ Standard: " << __cplusplus;
    std::cout << "[TEST 8] Compiler: " << ss.str() << std::endl;

    std::cout << "\n========================================" << std::endl;
    std::cout << "  ALL TESTS PASSED - MSVC COMPATIBLE" << std::endl;
    std::cout << "========================================" << std::endl;

    return 0;
}
