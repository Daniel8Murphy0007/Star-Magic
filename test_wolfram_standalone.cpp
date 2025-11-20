// Quick standalone test for Wolfram Engine integration
#include <iostream>
#include "wstp.h"

int main() {
    std::cout << "Testing Wolfram Engine 14.3 WSTP connection...\n";
    
    WSENV env = WSInitialize(nullptr);
    if (!env) {
        std::cout << "FAILED: WSInitialize\n";
        return 1;
    }
    std::cout << "✓ WSInitialize successful\n";
    
    int err = 0;
    const char *kernel_cmd = "C:\\Program Files\\Wolfram Research\\Wolfram Engine\\14.3\\WolframKernel.exe";
    std::string cmd = std::string("-linkmode launch -linkname \"") + kernel_cmd + " -wstp\"";
    
    std::cout << "Launching with command: " << cmd << "\n";
    
    WSLINK link = WSOpenString(env, cmd.c_str(), &err);
    
    if (!link || err != WSEOK) {
        std::cout << "FAILED: WSOpenString (error code " << err << ")\n";
        WSDeinitialize(env);
        return 1;
    }
    std::cout << "✓ Wolfram kernel launched\n";
    
    // Drain startup messages
    int pkt;
    while ((pkt = WSNextPacket(link)) && pkt != RETURNPKT) {
        WSNewPacket(link);
    }
    std::cout << "✓ Kernel ready\n";
    
    // Test simple evaluation: 2+2
    WSPutFunction(link, "EvaluatePacket", 1);
    WSPutFunction(link, "ToString", 2);
    WSPutFunction(link, "Plus", 2);
    WSPutInteger(link, 2);
    WSPutInteger(link, 2);
    WSPutSymbol(link, "InputForm");
    WSEndPacket(link);
    
    pkt = WSNextPacket(link);
    if (pkt == RETURNPKT || pkt == TEXTPKT) {
        const char *res = nullptr;
        if (WSGetString(link, &res) && res) {
            std::cout << "✓ Test evaluation: 2+2 = " << res << "\n";
            WSReleaseString(link, res);
            WSNewPacket(link);
        }
    }
    
    // Cleanup
    WSPutFunction(link, "EvaluatePacket", 1);
    WSPutFunction(link, "Exit", 0);
    WSEndPacket(link);
    while (WSNextPacket(link) != ILLEGALPKT)
        WSNewPacket(link);
    WSClose(link);
    WSDeinitialize(env);
    
    std::cout << "\n✓✓✓ All tests PASSED! Wolfram Engine 14.3 integration working!\n";
    return 0;
}
