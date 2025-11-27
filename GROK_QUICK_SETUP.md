# Grok AI Listener - Quick Setup (5 Minutes)

**Status:** ✅ Installed | ⚠️ Not Active | 🎯 Ready to Activate

---

## ⚡ Quick Start (Copy-Paste Ready)

### Step 1: Get API Key (2 minutes)
```
1. Visit: https://x.ai/api
2. Sign up (or log in)
3. Generate API key (starts with "xai-")
4. Copy the key
```

### Step 2: Set Environment Variable (30 seconds)
**PowerShell (Permanent):**
```powershell
[System.Environment]::SetEnvironmentVariable('XAI_API_KEY', 'xai-YOUR-KEY-HERE', 'User')
```

**Verify:**
```powershell
# Restart PowerShell first, then:
$env:XAI_API_KEY
# Should show: xai-YOUR-KEY-HERE
```

### Step 3: Test Connectivity (1 minute)
**Add to MAIN_1_CoAnQi.cpp (temporary test):**
```cpp
#include "source178_grok_api.cpp"

int main() {
    testGrokAPI();  // Should print "SUCCESS!"
    return 0;
}
```

**Compile & Run:**
```powershell
cmake --build build_msvc --config Release
.\build_msvc\Release\MAIN_1_CoAnQi.exe
```

**Expected Output:**
```
[Grok API] Calling Grok API with prompt: Hello, Grok!...
[Grok API] Response received (200)
[Grok API] Token usage - Prompt: 42, Completion: 18, Total: 60
SUCCESS! Grok response: [AI response here]
```

---

## 🎯 What You Get

### 5 AI Functions Ready to Use:

1. **`diagnoseCompilationError(errorMsg, file, line)`**
   - MSVC C++ error diagnostics with UQFF context
   - Example: Fix "error C2065: undeclared identifier"

2. **`explainPhysicsEquation(name, code)`**
   - Physics equation explanations (26D gravity, F_U, etc.)
   - Example: Explain "g(r,t) = Σ(i=1 to 26) [Ug1_i + Ug2_i + Ug3_i + Ug4_i]"

3. **`reviewPhysicsCode(snippet, aspect)`**
   - Code review for performance/correctness
   - Example: Optimize nested loops, check physics validity

4. **`callGrokAPI(prompt)`**
   - Custom AI queries with UQFF knowledge
   - Example: "How do I implement 26D folding in C++?"

5. **`testGrokAPI()`**
   - Connectivity test
   - Use first to verify setup

---

## 🔧 Already Integrated

- ✅ **File:** source178_grok_api.cpp
- ✅ **CMake:** Line 189 (Qt6 Network linked)
- ✅ **Error Links:** Line 34128 in MAIN_1_CoAnQi.cpp
- ✅ **AI Calls:** Line 34570 in MAIN_1_CoAnQi.cpp
- ✅ **UQFF Context:** Pre-loaded with 91,384 lines, 894 PhysicsTerm classes

**Your codebase knows about:**
- 26D geometric framework
- Universal Gravity (Ug1-Ug4)
- SCm (Superconductive Material)
- Windows threading (SimpleMutex)
- MSVC 14.44+ / C++20 requirements
- Wolfram WSTP integration

---

## 🚨 Troubleshooting (30 seconds each)

### "XAI_API_KEY not set"
```powershell
# Set it again:
[System.Environment]::SetEnvironmentVariable('XAI_API_KEY', 'xai-YOUR-KEY', 'User')
# Restart PowerShell
```

### "401 Unauthorized"
```
Regenerate key at https://x.ai/api
Check for typos in environment variable
```

### "Empty Response"
```powershell
# Check if key is loaded:
$env:XAI_API_KEY
# Should NOT be empty
```

---

## 📚 Full Documentation

- **Complete Guide:** [GROK_ACTIVATION_GUIDE.md](GROK_ACTIVATION_GUIDE.md)
- **Usage Examples:** See GROK_ACTIVATION_GUIDE.md sections 6-7
- **API Reference:** https://docs.x.ai/

---

## 🎓 Example: Fix a Real Error

```cpp
// You see this error:
// error C2664: 'void PhysicsTerm::calculate(double)' : cannot convert argument 1 from 'int' to 'double'

// Ask Grok:
QString fix = diagnoseCompilationError(
    "error C2664: cannot convert argument 1 from 'int' to 'double'",
    "source14.cpp",
    456
);
std::cout << fix.toStdString() << std::endl;

// Grok responds:
// "The error occurs because you're passing an integer to a function expecting a double.
//  In C++20 with MSVC, implicit int→double conversion is strict...
//  Fix: calculate(static_cast<double>(myInt));"
```

---

## ⏱️ Time Investment

- **Setup:** 5 minutes (one-time)
- **Per Query:** 1-3 seconds (network latency)
- **Benefit:** AI-assisted debugging, physics explanations, code optimization

**Worth it?** Yes - especially for complex MSVC errors and UQFF physics questions.

---

**Next Steps:**
1. ✅ Set XAI_API_KEY (Step 2 above)
2. ✅ Test with `testGrokAPI()` (Step 3 above)
3. ✅ Integrate into error handling workflow
4. 📖 Read GROK_ACTIVATION_GUIDE.md for advanced usage

---

*Grok = xAI's AI assistant, optimized for C++ and physics - pre-trained on your UQFF codebase context*
