
// source174_wolfram_bridge_embedded.cpp
// FULL WSTP implementation with MSVC

#ifdef USE_EMBEDDED_WOLFRAM

#include <iostream>
#include <string>
#include <array> // MSVC requirement
#include "wstp.h"

static WSENV ws_env = nullptr;
static WSLINK ws_link = nullptr;

bool InitializeWolframKernel()
{
    if (ws_link)
        return true;

    std::cout << "[WSTP] Initializing environment...\n" << std::flush;
    ws_env = WSInitialize(nullptr);
    if (!ws_env)
    {
        std::cout << "[WSTP] WSInitialize failed\n" << std::flush;
        return false;
    }
    std::cout << "[WSTP] Environment initialized successfully\n" << std::flush;

    int err = 0;
    // Launch kernel directly with wolfram.exe -mathlink
    char* argv[] = {
        (char*)"CoAnQi",
        (char*)"-linkname",
        (char*)"\"C:\\Program Files\\Wolfram Research\\Wolfram Engine\\14.3\\wolfram.exe\" -mathlink -nogui",
        (char*)"-linkmode",
        (char*)"launch",
        nullptr
    };
    
    std::cout << "[WSTP] About to launch kernel with WSOpenArgv...\n" << std::flush;
    std::cout << "[WSTP] Launch string: " << argv[2] << "\n" << std::flush;
    std::cout << "[DEBUG] Launch command: " << argv[4] << "\n" << std::flush;
    
    ws_link = WSOpenArgv(ws_env, argv, argv + (sizeof(argv)/sizeof(argv[0]) - 1), &err);

    std::cout << "[WSTP] WSOpenArgv returned, ws_link = " << (void*)ws_link << ", err = " << err << "\n" << std::flush;

    if (!ws_link || err != WSEOK)
    {
        std::cout << "[WSTP] Failed to launch kernel (error code: " << err << ")\n" << std::flush;
        if (ws_link) {
            const char* errmsg = WSErrorMessage(ws_link);
            if (errmsg) {
                std::cout << "[WSTP] Error message: " << errmsg << "\n" << std::flush;
            }
        }
        if (ws_env) WSDeinitialize(ws_env);
        ws_env = nullptr;
        return false;
    }

    std::cout << "[WSTP] Kernel link created, activating...\n" << std::flush;
    
    // Activate the link
    if (!WSActivate(ws_link))
    {
        int ws_error = WSError(ws_link);
        std::cout << "[WSTP] Kernel failed to activate (error code: " << ws_error << ")\n" << std::flush;
        const char* errmsg = WSErrorMessage(ws_link);
        if (errmsg) {
            std::cout << "[WSTP] Error message: " << errmsg << "\n" << std::flush;
        }
        WSClose(ws_link);
        WSDeinitialize(ws_env);
        ws_link = nullptr;
        ws_env = nullptr;
        return false;
    }

    std::cout << "[WSTP] Kernel activated, draining startup packets...\n" << std::flush;
    
    // Drain startup packets with safety limit
    int pkt;
    int drain_count = 0;
    while ((pkt = WSNextPacket(ws_link)) && pkt != RETURNPKT && drain_count < 20)  // Safety limit
    {
        WSNewPacket(ws_link);
        drain_count++;
    }

    if (drain_count >= 20)
    {
        std::cout << "[WSTP] Warning: Startup drain took too many packets — possible loop. Continuing...\n" << std::flush;
    }

    std::cout << "[WSTP] Startup packets drained (" << drain_count << " packets).\n" << std::flush;
    std::cout << "[WSTP] Wolfram kernel connected successfully!\n" << std::flush;
    return true;
}

std::string WolframEvalToString(const std::string &code)
{
    if (!InitializeWolframKernel())
        return "KERNEL_NOT_AVAILABLE";

    // Send evaluation packet
    WSPutFunction(ws_link, "EvaluatePacket", 1);
    WSPutFunction(ws_link, "ToString", 2);
    WSPutString(ws_link, code.c_str());
    WSPutSymbol(ws_link, "InputForm");
    WSEndPacket(ws_link);

    // Read result - skip all packets until we get ReturnPacket
    int pkt;
    while ((pkt = WSNextPacket(ws_link)) != RETURNPKT)
    {
        if (pkt == 0 || pkt == ILLEGALPKT)
        {
            WSNewPacket(ws_link);
            return "<error: lost connection>";
        }
        WSNewPacket(ws_link);
    }

    // Now read the actual return value
    const char *res = nullptr;
    if (WSGetString(ws_link, &res))
    {
        std::string result = res ? std::string(res) : "<null>";
        if (res) WSReleaseString(ws_link, res);
        WSNewPacket(ws_link);
        return result;
    }

    int err = WSError(ws_link);
    WSNewPacket(ws_link);
    return "<error: failed to read result, code " + std::to_string(err) + ">";
}

void WolframCleanup()
{
    if (ws_link)
    {
        // Politely ask kernel to exit
        WSPutFunction(ws_link, "EvaluatePacket", 1);
        WSPutFunction(ws_link, "Exit", 0);
        WSEndPacket(ws_link);
        while (WSNextPacket(ws_link) != ILLEGALPKT)
            WSNewPacket(ws_link);
        WSClose(ws_link);
        ws_link = nullptr;
    }
    if (ws_env)
    {
        WSDeinitialize(ws_env);
        ws_env = nullptr;
    }
}

void WolframEmbeddedBridge()
{
    std::cout << "\n=== WSTP Kernel Interface (MSVC + TCPIP) ===\n";

    if (!InitializeWolframKernel())
        return;

    std::string test = "c^4/(8 Pi G) R - F_{\\mu\\nu} F^{\\mu\\nu}/4 + (D_\\mu \\phi)^2";

    std::cout << "Sending   : " << test << "\n";
    std::string result = WolframEvalToString(test);
    std::cout << "Received  : " << result << "\n\n";

    // Demonstrate statefulness
    WolframEvalToString("uqffTest = c^4 R/(8 Pi G) - Fμν Fμν/4 + StarMagicAether;");
    std::cout << "Stateful  : " << WolframEvalToString("uqffTest") << "\n";
    
    std::cout << "\n** Protocol: WSTP via TCPIP (no file dialogs) **\n";
    std::cout << "** Compiler: MSVC **\n";
}

#endif // USE_EMBEDDED_WOLFRAM