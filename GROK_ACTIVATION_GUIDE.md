# Grok AI Integration - Activation Guide

**Integration Date:** November 25, 2025 (Commit 33cdfcb)  
**Updated:** April 2026 (Activation Plan execution)  
**Phase:** Phase 23 - AI Integration  
**Status:** ✅ **CODE COMPLETE + ACTIVATION IMPROVEMENTS APPLIED**  
**File:** `source178_grok_api.cpp` (included directly by MAIN_1_CoAnQi.cpp)  
**Python side:** `GrokAPI.py` + `APIKeyManager.py` (persistent config file support)  
**CMake:** Automatic copy of the two .py files next to the executable

**Current Architecture (2026):** Pure Python subprocess (no Qt dependency for Grok). The C++ side shells out to `python GrokAPI.py`. API keys are managed via `grok_api_config.json` (preferred) + `XAI_API_KEY` environment variable (fallback).

---

## Overview

The **Grok listener** is an xAI API integration that provides AI-assisted C++ error diagnostics, physics equation explanations, and code review capabilities specifically tailored for the UQFF (Unified Quantum Field Framework) project.

### What Grok Knows About Your Codebase

Grok is pre-configured with UQFF context:
- **91,384-line C++20 codebase** implementing 894 PhysicsTerm classes
- **Quantum field unity calculations** and 26-dimensional geometric framework
- **MSVC 14.44+ compiler** requirements and C++20 standard specifics
- **Windows threading model** (SimpleMutex, SimpleLockGuard)
- **Physics terminology** (Universal Gravity, SCm, 26D polynomial structures)

### Available Functions (C++ side)

The Grok integration exposes these functions (usable from MAIN_1_CoAnQi menu or directly in code):

1. **`callGrokAPI(const std::string& prompt)`** - Core function (implemented in source178_grok_api.cpp)
   - Delegates to GrokAPI.py (which prefers grok_api_config.json then XAI_API_KEY env var)
   - Returns plain string response (empty on failure)

2. **`diagnoseCompilationError(...)`, `explainPhysicsEquation(...)`, `reviewPhysicsCode(...)`** - Convenience wrappers that build good prompts and call the above.

3. **`testGrokAPI()`** - Menu-driven connectivity test (highly recommended first step after configuration).

4. **`configureGrokAPIKey()`** - Interactive key setup that persists to config file + sets env for the session.

5. **`printGrokAPIStatus()`** - One-line status printed at interactive startup.

### API Configuration

- **Endpoint:** `https://api.x.ai/v1/chat/completions`
- **Model:** `grok-beta`
- **Temperature:** 0.7 (balanced creativity/precision)
- **Max Tokens:** 2048 (response length limit)
- **Authentication:** Bearer token via `XAI_API_KEY` environment variable

---

## Activation Steps

### Step 1: Obtain xAI API Key

1. Visit https://x.ai/api
2. Sign up for a free xAI account (or log in if you already have one)
3. Navigate to API Keys section
4. Generate a new API key (starts with `xai-`)
5. **Copy the key immediately** - you won't be able to see it again

### Step 2: Set Environment Variable

Choose **EITHER** temporary OR permanent setup:

#### Option A: Temporary (Current PowerShell Session Only)
```powershell
$env:XAI_API_KEY = "xai-your-actual-key-here"
```
- **Pros:** Quick testing, no system changes
- **Cons:** Expires when you close PowerShell

#### Option B: Permanent (System-Wide, All Sessions)
```powershell
[System.Environment]::SetEnvironmentVariable('XAI_API_KEY', 'xai-your-actual-key-here', 'User')
```
- **Pros:** Persistent across sessions, no re-entry needed
- **Cons:** Requires restarting terminal/VS Code to take effect

**Verify Environment Variable:**
```powershell
$env:XAI_API_KEY
# Should print: xai-your-actual-key-here
```

If empty, restart your terminal and check again.

### Step 3: Rebuild MAIN_1_CoAnQi (Optional)

Grok is already compiled into your executable, but if you want to ensure the latest version:

```powershell
# Navigate to Star-Magic directory
cd C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic

# Clean rebuild (Visual Studio 2022)
Remove-Item -Recurse -Force build_msvc -ErrorAction SilentlyContinue
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi
```

### Step 4: Test Connectivity

Add test code to `MAIN_1_CoAnQi.cpp` or create a standalone test file:

```cpp
#include "source178_grok_api.cpp"

int main() {
    // Test 1: Basic connectivity
    testGrokAPI();  // Should print "SUCCESS!" with Grok's response
    
    // Test 2: Error diagnostics
    QString errorDiag = diagnoseCompilationError(
        "error C2065: 'undefined_var' : undeclared identifier",
        "MAIN_1_CoAnQi.cpp",
        1234
    );
    std::cout << "\n=== Error Diagnosis ===\n" << errorDiag.toStdString() << std::endl;
    
    // Test 3: Physics explanation
    QString physicsExplain = explainPhysicsEquation(
        "F_U Primary Equation",
        "F_U = c^4/(8πG) * R_μν"
    );
    std::cout << "\n=== Physics Explanation ===\n" << physicsExplain.toStdString() << std::endl;
    
    return 0;
}
```

Compile and run:
```powershell
cmake --build build_msvc --config Release
.\build_msvc\Release\MAIN_1_CoAnQi.exe
```

**Expected Output:**
```
[Grok API] Calling Grok API with prompt: Hello, Grok! This is a test...
[Grok API] Response received (200)
[Grok API] Token usage - Prompt: 42, Completion: 18, Total: 60
SUCCESS! Grok response: [Grok's actual response here]

=== Error Diagnosis ===
[Detailed error explanation and fix suggestions from Grok]

=== Physics Explanation ===
[Physics equation explanation in UQFF context]
```

**If You See Errors (2026 improved diagnostics):**
- Key not found (config or env) → Use the in-app "Configure Grok API Key" menu option first.
- Python execution failure → Install Python 3.9+, run `pip install requests`, ensure the two .py files are next to the .exe (CMake now copies them automatically).
- 401 / rate limit → Check key validity at https://x.ai/api.

---

## Integration in MAIN_1_CoAnQi.cpp

Grok is already integrated in your error handling workflow:

### Error Message Links (Line 34128)
```cpp
// MathErrorListener includes Grok Explain links
<a href=\"grok_explain:line:34128:error_type\">Grok Explain</a>
```

### AI Response Calls (Line 34570)
```cpp
QString aiResponse = callGrokAPI(aiPrompt);
// Grok provides context-aware assistance
```

### Placeholder Code (Line 59925)
```cpp
// Example usage (commented out):
QString apiKey = "your_grok_api_key";  // Replace with environment variable
```

---

## CMake / Runtime Notes (2026)

- No Qt is required for Grok anymore.
- CMake automatically stages `GrokAPI.py` + `APIKeyManager.py` next to the executable.
- Run the program from a shell with Python 3 + requests installed. The menu options will guide you.

---

## Usage Examples

### Example 1: Diagnose MSVC Compilation Error
```cpp
QString errorMsg = "error C2664: 'void PhysicsTerm::calculate(double)' : cannot convert argument 1 from 'int' to 'double'";
QString sourceFile = "source14.cpp";
int lineNumber = 456;

QString diagnosis = diagnoseCompilationError(errorMsg, sourceFile, lineNumber);
std::cout << diagnosis.toStdString() << std::endl;
```

**Grok Response (Example):**
```
The error occurs because you're passing an integer to a function expecting a double.
In C++20 with MSVC, implicit int→double conversion is strict in some contexts.

Fix:
1. Explicit cast: calculate(static_cast<double>(myInt));
2. Or define variable as double: double myValue = 123.0;

In your UQFF context, ensure physics constants use double literals (1.0, not 1).
```

### Example 2: Explain Physics Equation
```cpp
QString eqName = "26D Compressed Gravity (Ug1-Ug4)";
QString eqCode = "g(r,t) = Σ(i=1 to 26) [Ug1_i + Ug2_i + Ug3_i + Ug4_i]";

QString explanation = explainPhysicsEquation(eqName, eqCode);
std::cout << explanation.toStdString() << std::endl;
```

**Grok Response (Example):**
```
This equation represents the 26-layer compressed gravity framework in UQFF:
- Ug1: Internal dipole strength (stellar irregularities)
- Ug2: Spherical outer field bubble (heliospheres)
- Ug3: Disk of magnetic strings (planetary orbit stability)
- Ug4: Observable star/black hole interactions

The summation over 26 dimensions reflects UQFF's geometric folding substrate,
where each layer i has quantum state factors Q_i, [UA]_i, [SCm]_i.

This is fundamentally different from Newtonian G-based gravity - it emerges
from 26D→3D geometric projection, not force laws.
```

### Example 3: Review Code Performance
```cpp
QString code = R"(
for (int i = 0; i < systems.size(); i++) {
    double result = 0.0;
    for (int j = 0; j < 1000000; j++) {
        result += sqrt(systems[i].mass * j);
    }
    systems[i].energyDensity = result;
}
)";

QString review = reviewPhysicsCode(code, "performance");
std::cout << review.toStdString() << std::endl;
```

**Grok Response (Example):**
```
**Performance Issues:**
1. Inner loop executes 1M iterations per system - highly inefficient
2. sqrt() called 1M times unnecessarily (result is constant for each system)
3. No parallelization despite independent system calculations

**Optimizations:**
1. Cache sqrt(systems[i].mass) outside inner loop
2. Use analytical formula: result = sqrt(mass) * (sum from 0 to 999999)
3. Parallelize outer loop with Windows threads (SimpleMutex already available)

**Corrected Code:**
for (int i = 0; i < systems.size(); i++) {
    double sqrtMass = sqrt(systems[i].mass);
    double sum = 1000000.0 * 999999.0 / 2.0;  // Analytical sum
    systems[i].energyDensity = sqrtMass * sum;
}

**UQFF Context:** Ensure mass units are consistent with your SCm framework.
```

---

## Troubleshooting

### "XAI_API_KEY not set in environment"
**Cause:** Environment variable missing or not loaded  
**Fix:**
```powershell
# Set variable
[System.Environment]::SetEnvironmentVariable('XAI_API_KEY', 'xai-your-key', 'User')

# Restart PowerShell/VS Code
# Verify
$env:XAI_API_KEY
```

### "QNetworkReply::NetworkError"
**Cause:** Network connectivity issue or invalid endpoint  
**Fix:**
- Check internet connection
- Verify firewall allows HTTPS to api.x.ai
- Test manually: `curl -X POST https://api.x.ai/v1/chat/completions -H "Authorization: Bearer YOUR_KEY" -H "Content-Type: application/json" -d '{"messages":[{"role":"user","content":"test"}],"model":"grok-beta"}'`

### "401 Unauthorized"
**Cause:** Invalid or expired API key  
**Fix:**
- Regenerate key at https://x.ai/api
- Check for typos in environment variable
- Ensure key starts with `xai-`

### "429 Rate Limit Exceeded"
**Cause:** Too many API calls in short time  
**Fix:**
- Implement exponential backoff in code
- Reduce call frequency
- Upgrade xAI plan if needed

### "Empty Response from Grok"
**Cause:** API key not set, network error, or malformed request  
**Fix:**
- Check `qWarning()` and `qDebug()` console output for error details
- Enable verbose logging: Add `QLoggingCategory::setFilterRules("qt.network.ssl=true");` to main()
- Test with `testGrokAPI()` first before complex calls

---

## Best Practices

### 1. Cache Responses
Don't call Grok repeatedly for the same question:
```cpp
std::map<QString, QString> responseCache;

QString getCachedResponse(const QString& prompt) {
    if (responseCache.find(prompt) != responseCache.end()) {
        return responseCache[prompt];
    }
    QString response = callGrokAPI(prompt);
    responseCache[prompt] = response;
    return response;
}
```

### 2. Handle Empty Responses
```cpp
QString response = callGrokAPI(prompt);
if (response.isEmpty()) {
    qWarning() << "Grok API call failed - check environment variable and network";
    // Fallback to local error handling
}
```

### 3. Respect Rate Limits
Implement delays between rapid calls:
```cpp
#include <chrono>
#include <thread>

void rateLimitedGrokCall(const QString& prompt) {
    static auto lastCallTime = std::chrono::steady_clock::now();
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - lastCallTime).count();
    
    if (elapsed < 1000) {  // Minimum 1 second between calls
        std::this_thread::sleep_for(std::chrono::milliseconds(1000 - elapsed));
    }
    
    callGrokAPI(prompt);
    lastCallTime = std::chrono::steady_clock::now();
}
```

### 4. Sanitize User Input
If prompts include user-provided code:
```cpp
QString sanitizeCode(const QString& code) {
    QString sanitized = code;
    sanitized.replace("\"", "\\\"");  // Escape quotes for JSON
    sanitized.replace("\n", "\\n");   // Preserve newlines
    return sanitized;
}
```

---

## API Cost Considerations

**xAI Grok Beta Pricing (as of Nov 2025):**
- Free tier: Limited requests/month
- Paid tier: Per-token pricing

**Token Usage Tracking:**
The `callGrokAPI()` function automatically logs token usage:
```
[Grok API] Token usage - Prompt: 142, Completion: 87, Total: 229
```

Monitor console output to track API costs.

---

## Security Notes

1. **Never commit API keys to git:**
   ```bash
   # Add to .gitignore (already included):
   .env
   *.key
   ```

2. **Use environment variables only:**
   ```cpp
   // ✅ CORRECT
   QString apiKey = qgetenv("XAI_API_KEY");
   
   // ❌ WRONG (hardcoded key)
   QString apiKey = "xai-abc123...";
   ```

3. **Rotate keys periodically:**
   - Generate new key at x.ai/api
   - Update environment variable
   - Delete old key from xAI dashboard

---

## Next Steps

After activation:
1. ✅ Set `XAI_API_KEY` environment variable
2. ✅ Test connectivity with `testGrokAPI()`
3. ✅ Integrate into error handling workflow
4. ✅ Use for physics equation explanations
5. ✅ Leverage for code review and optimization
6. ⏳ Update workspace documentation (README.md/PLAN.md) ← **YOU ARE HERE**
7. ⏳ Document custom prompts for UQFF-specific queries
8. ⏳ Create issue templates for Grok-assisted debugging

---

## References

- **xAI Documentation:** https://docs.x.ai/
- **API Endpoint:** https://api.x.ai/v1/chat/completions
- **Integration File:** `source178_grok_api.cpp`
- **CMake Config:** `CMakeLists.txt` (lines 182-220)
- **Git Commit:** 33cdfcb (November 25, 2025)
- **Phase Documentation:** Phase 23 - AI Integration

---

**Status Summary:**
- ✅ Code implementation complete
- ✅ CMake integration verified
- ✅ Qt6 Network dependencies linked
- ⚠️ **AWAITING:** User API key setup (Step 1-2 above)
- ⏳ **NEXT:** Test connectivity and integrate into workflow

---

*For questions or issues with Grok integration, contact Daniel T. Murphy at daniel.murphy00@gmail.com*
