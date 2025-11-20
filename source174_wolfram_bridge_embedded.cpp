
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

    ws_env = WSInitialize(nullptr);
    if (!ws_env)
    {
        std::cout << "WSInitialize failed\n";
        return false;
    }

    int err = 0;
    // Use loopback connection - no file dialogs, no PATH issues
    // First create listener, then connect kernel to it
    std::string listen_cmd = "-linkmode listen -linkprotocol TCPIP";
    
    ws_link = WSOpenString(ws_env, listen_cmd.c_str(), &err);

    if (!ws_link || err != WSEOK)
    {
        std::cout << "Failed to create WSTP listener (error code: " << err << ")\n";
        if (ws_env) WSDeinitialize(ws_env);
        ws_env = nullptr;
        return false;
    }

    // Get the link name that the kernel should connect to
    const char *link_name = WSLinkName(ws_link);
    if (!link_name)
    {
        std::cout << "Failed to get listener link name\n";
        WSClose(ws_link);
        WSDeinitialize(ws_env);
        ws_link = nullptr;
        ws_env = nullptr;
        return false;
    }

    std::cout << "Waiting for kernel connection on: " << link_name << "\n";
    
    // Launch kernel in background to connect to our listener
    std::string kernel_launch = "start /B \"\" \"C:\\Program Files\\Wolfram Research\\Wolfram Engine\\14.3\\WolframKernel.exe\" -wstp -linkmode connect -linkprotocol TCPIP -linkname \"" + std::string(link_name) + "\"";
    system(kernel_launch.c_str());

    // Wait for kernel to connect (with timeout)
    if (!WSActivate(ws_link))
    {
        std::cout << "Kernel failed to connect (timeout or error)\n";
        std::cout << "Hint: Ensure Wolfram Engine is activated with 'wolframscript -activate'\n";
        WSClose(ws_link);
        WSDeinitialize(ws_env);
        ws_link = nullptr;
        ws_env = nullptr;
        return false;
    }

    // Drain startup messages
    int pkt;
    while ((pkt = WSNextPacket(ws_link)) && pkt != RETURNPKT)
    {
        WSNewPacket(ws_link);
    }

    std::cout << "Wolfram kernel connected successfully via TCPIP.\n";
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