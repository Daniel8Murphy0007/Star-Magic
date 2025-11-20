// physics_term_extractor_main.cpp
// Standalone executable to extract all 500+ physics terms using Wolfram Engine 14.3

#include "extract_all_physics_terms.cpp"

int main() {
    std::cout << "========================================\n";
    std::cout << "Physics Term Extractor (Wolfram-Powered)\n";
    std::cout << "========================================\n\n";
    
#ifdef USE_EMBEDDED_WOLFRAM
    ExtractAllPhysicsTerms();
    
    std::cout << "\n\nExtraction complete!\n";
    std::cout << "Output files:\n";
    std::cout << "  - extracted_physics_terms_complete.h\n";
    std::cout << "  - register_all_500plus_terms.cpp\n";
#else
    std::cout << "ERROR: Wolfram integration not enabled!\n";
    std::cout << "Rebuild with: cmake -DUSE_EMBEDDED_WOLFRAM=ON\n";
#endif
    
    return 0;
}
