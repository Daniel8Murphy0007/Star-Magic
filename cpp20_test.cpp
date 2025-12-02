#include <iostream>
#include <version>
int main() {
    std::cout << "C++ Standard: " << __cplusplus << std::endl;
    #ifdef _MSVC_LANG
    std::cout << "MSVC Language: " << _MSVC_LANG << std::endl;
    #endif
    std::cout << "MSVC Version: " << _MSC_VER << std::endl;
    return 0;
}
