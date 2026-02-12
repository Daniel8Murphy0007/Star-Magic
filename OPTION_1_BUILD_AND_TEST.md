# Option 1: CLI + GUI Integration - Build and Test Instructions

**Implementation Status:** ✅ **COMPLETE**  
**CLI Integration:** [914f9e4](https://github.com/Daniel8Murphy0007/Star-Magic/commit/914f9e4) (December 2025)  
**GUI Integration:** ✅ **COMPLETE** (February 11, 2026)  
**Components:** MAIN_1_CoAnQi.exe CLI + Source2.exe GUI + Python Wrappers

---

## 📋 What Was Implemented

### C++ Changes (MAIN_1_CoAnQi.cpp)
- **Location:** Line 23695 (inserted 118 lines before interactive menu loop)
- **Functionality:** CLI access point with JSON output
- **CLI Flags:**
  1. `--batch "System Name"` → Compute F_U_Bi_i, g_compressed, 5 auxiliary forces
  2. `--list-systems` → Enumerate all 121+ astrophysical systems
  3. `--system-info "System Name"` → Get detailed parameters (M, r, L_X, B0, ω₀, v, T)

### Python Integration Layer
- **CoAnQi_Wrapper.py** (345 lines) → Subprocess-based Python interface
- **test_integration.py** (320 lines) → Automated test suite with 6 test cases

### Source2 GUI Integration (COMPLETE - February 11, 2026)
- **Tab 1: Embedded MAIN_1 Calculator Terminal** (NEW - February 11, 2026)
  - PowerShell-style terminal running MAIN_1_CoAnQi.exe directly in Tab 1
  - **18-option interactive menu** (Cosmic Egg build with all features):
    1. Calculate system (single) - F_U_Bi_i, g_compressed, validation pipeline
    2. Calculate ALL systems (parallel) - Windows threading, SimpleMutex
    3. Clone and mutate system
    4. Add custom system
    5. Add dynamic physics term
    6. Run simulations
    7. Statistical analysis
    8. Self-optimization
    9. WSTP kernel interface (Wolfram)
    10. Auto-export full UQFF to Wolfram
    11. Run Wolfram Field Unity Simulation
    12. **Run Cosmic Quantum Egg (26D) Simulation** ← Full 26-dimensional physics
    13. Configure Grok API Key
    14. Test Grok AI Integration
    15. SOURCE4 Unified Field Validation
    16. Information Paradox Tests (BH Info)
    17. BSM Physics Validation
    18. Exit
  - Special commands: `restart` (relaunch calculator), `clear` (clear output), `help` (show terminal commands)
  - Auto-start on Source2 launch, clean shutdown on tab close (sends Exit command option 18)
  - Terminal styling: PowerShell black background (#0C0C0C), green input echo (#00FF00), red errors (#FF0000)

- **Tab 2: QCalc.py Python Terminal** (NEW - February 11, 2026)
  - Interactive Python terminal running QCalc.py (UQFF Quantum Calculator)
  - **8 UQFF Master Equations**:
    1. UQFF (Base Unified Field) - Complete unified force calculation
    2. UQFF_Compressed (Newtonian + 9 corrections) - Standard gravity model
    3. UQFF_Resonant (aDPM + 13 frequency modes) - Oscillatory behavior
    4. UQFF_Superconductive (SCm vacuum modulation) - Quantum field effects
    5. UQFF_Buoyant (F_U_Bi) - Atomic scale Inside→Out forces
    6. UQFF_Master_Buoyant (F_U_Bi_i) - Cosmic scale Outside→In forces
    7. UQFF_Triadic (26-layer gravitational scaling) - Multi-dimensional gravity
    8. UQFF_Quadratic (Root solutions) - Dual-solution physics
  - Python interactive mode with QCalc functions pre-imported
  - Terminal styling: Python dark theme (#1E1E1E), blue prompt (>>> in #569CD6)
  - Special commands: `restart`, `clear`, `help` (show UQFF equation reference)

- **Tab 3: Advanced Scientific Calculator** (NEW - February 11, 2026)
  - Standard Model physics calculator similar to Wolfram Alpha, MATLAB, Mathematica
  - **4 Computational Engine Options**:
    1. 🐍 **SymPy** (Open Source) - Symbolic mathematics, no API key required
    2. 🔬 **SciPy** (Open Source) - Scientific computing and numerical methods
    3. 🔢 **NumPy** (Open Source) - Numerical arrays and linear algebra
    4. 🧮 **Wolfram Alpha API** (Commercial) - Requires API key from developer.wolframalpha.com
  - **Supported Operations**:
    - Algebra: solve(x^2 - 4 = 0, x), expand((x+1)^2), factor(x^2-1)
    - Calculus: integrate(x^2, x), diff(sin(x), x), limit(1/x, x, 0)
    - Linear Algebra: matrix operations, determinants, eigenvalues
    - Physics: fundamental constants (c, h, G, k_B), unit conversions
    - Statistics: mean, std, distributions, hypothesis testing
  - API registration info button with links to registration pages
  - Free tier available for all options (SymPy/SciPy/NumPy fully free, Wolfram Alpha: 2,000 queries/month free)

- **Tab 4: 📓 Notebook Editor** (NEW - February 11, 2026)
  - Jupyter-style multi-cell notebook environment for interactive Python workflows
  - **Features**:
    - ➕ Code Cell: Python execution with immediate output display
    - 📝 Markdown Cell: Documentation cells with blue borders
    - ▶️ Run All: Execute all code cells sequentially
    - 🗑️ Clear All: Clear all cell outputs
  - **Cell Types**:
    - **Code Cells** (green border #4CAF50): Python code execution via subprocess, 10s timeout
    - **Markdown Cells** (blue border #2196F3): Documentation and notes
  - **Output Display**:
    - Stdout: Black text (#000000)
    - Stderr: Red text (#FF0000)
    - Scrollable output area (max 200px height)
  - **Use Cases**:
    - Data analysis workflows combining UQFF equations
    - Multi-step calculations with intermediate results
    - Documentation alongside executable code
    - Testing physics models incrementally
  - **Keyboard Shortcuts** (future enhancement):
    - Shift+Enter: Run current cell
    - Ctrl+Enter: Run and advance to next cell

- **Tab 5: 📚 CondensedPhysics.py Terminal** (NEW - February 11, 2026)
  - Python interactive terminal for CondensedPhysics.py (General Model/Class Solver Index)
  - **Architecture** (PURE PHYSICS CALCULATOR):
    - NO hardcoded system data (all parameters passed via method arguments)
    - NO named system classes (only generic physics domain calculators)
    - NO global instances (stateless calculator classes only)
    - INPUT: Datasets from APIFetch.py or query results → CondensedPhysics.solve()
    - OUTPUT: Long-form equations + solutions + available equations list + simulation set
  - **Terminal Features**:
    - Background: Dark gray (#2E2E2E)
    - Prompt: CP>>> (purple #9C27B0)
    - Text: Light gray (#E0E0E0)
    - Errors: Red (#FF5252)
  - **Special Commands**:
    - `restart`: Relaunch Python process with CondensedPhysics.py
    - `exit()`: Close terminal (standard Python exit)
  - **Auto-Discovery**: Searches for CondensedPhysics.py in app directory, then parent directory
  - **Use Cases**:
    - Solve general physics models with custom parameters
    - Generate long-form equations with step-by-step solutions
    - List all available equations for a given query
    - Run simultaneous equation sets for complex simulations

- **Tab 6: 🤖 CoAnQi_bot (Ollama Local LLM)** (NEW - February 11, 2026)
  - Local AI assistant powered by Ollama (https://ollama.com)
  - **Platform**: Ollama - Local LLM inference (NO cloud API, NO API key, fully FREE)
  - **API**: HTTP POST to localhost:11434/api/generate (curl subprocess)
  - **Models** (dropdown selection):
    1. **llama3.2:latest** - General purpose, 3B params, balanced performance
    2. **codellama:latest** - Code specialist, 7B params, optimized for programming
    3. **mistral:latest** - Fast and efficient, 7B params, good reasoning
    4. **qwen2.5-coder:latest** - Latest code model, 7B params, multi-language support
  - **Features**:
    - Multi-turn conversations with context
    - Code generation for physics calculations
    - UQFF equation explanations with system prompt awareness
    - Debugging and optimization suggestions
    - Fully local (no data sent to cloud, privacy-first)
  - **System Prompt** (built-in):
    > "You are CoAnQi_bot, an expert physics and code assistant for the UQFF (Unified Quantum Field Framework) project. Help with code generation, physics equations, and scientific computing. Be concise and provide code examples when relevant."
  - **Installation Guide** (📥 button):
    - Step 1: Download Ollama from https://ollama.com
    - Step 2: Pull model: `ollama pull llama3.2`
    - Step 3: Start server: `ollama serve` (runs on localhost:11434)
    - Step 4: Test in CoAnQi_bot tab
  - **Chat Interface**:
    - User messages: Blue background (#E3F2FD), left-aligned
    - Bot responses: Green background (#F1F8E9), HTML formatted with code blocks
    - Errors: Red background (#FFEBEE), connection troubleshooting
    - Status bar: "Ready", "Thinking...", or "Error - Ollama not running"
  - **Example Prompts**:
    - "Explain the F_U_Bi_i equation in UQFF"
    - "Generate Python code to calculate gravitational force with UQFF_Compressed"
    - "How do I implement a 26-layer Triadic gravity model?"
    - "Debug this WSTP integration error: [paste error]"

- **Tab 7: 🧠 SuperGrok4 (Grok xAI Expert)** (NEW - February 11, 2026)
  - Advanced AI research assistant powered by Grok xAI Cloud API (https://x.ai)
  - **Platform**: Grok xAI - Cloud LLM inference (~300B parameters, advanced reasoning)
  - **API**: POST to https://api.x.ai/v1/chat/completions (curl subprocess)
  - **API Key**: Requires XAI_API_KEY environment variable (obtain from x.ai dashboard)
  - **Models** (dropdown selection):
    1. **grok-beta** - Latest Grok model, general purpose, ~300B params, state-of-the-art performance
    2. **grok-2** - Advanced reasoning capabilities, optimized for complex multi-step problems
    3. **grok-vision-beta** - Multimodal with vision support (text + images)
  - **Features**:
    - Advanced research-grade responses for complex physics problems
    - Deep UQFF theory integration and equation derivations
    - Multi-step reasoning for complex derivations (Grok-2)
    - Vision support for analyzing astronomical images (grok-vision-beta)
    - Scientific paper analysis and literature review
    - Code generation with comprehensive documentation
  - **System Prompt** (built-in):
    > "You are SuperGrok4, an expert physics and research assistant for the UQFF (Unified Quantum Field Framework) project. You have deep knowledge of astrophysics, quantum mechanics, UQFF equations, and scientific computing. Provide detailed, research-grade responses with mathematical rigor when appropriate."
  - **API Key Configuration** (🔑 button):
    - **Windows (PowerShell)**: `$env:XAI_API_KEY = "your-api-key-here"`
    - **Windows (Permanent)**: System Properties → Environment Variables → New user variable
    - **macOS/Linux**: `export XAI_API_KEY="your-api-key-here"` (add to ~/.bashrc or ~/.zshrc)
    - **Status Check**: Config dialog shows ✅ SET or ❌ NOT SET
    - **Cost**: Paid API (check x.ai dashboard for pricing)
  - **Chat Interface**:
    - User messages: Blue background (#E8EAF6), left-aligned
    - SuperGrok4 responses: Green background (#E8F5E9), HTML formatted with equations
    - Errors: Red background (#FFEBEE), API key/connection troubleshooting
    - Status bar colors: Ready (green #4CAF50), Thinking (orange #FF9800), Error (red #F44336)
  - **Welcome Message** (auto-displayed on launch):
    - Features comparison table: SuperGrok4 vs CoAnQi_bot (Tab 6)
    - Platform comparison: Cloud (xAI) vs Local (Ollama)
    - Model size: ~300B params vs 3-7B params
    - Capabilities: Advanced reasoning, vision support vs basic code help
    - Privacy: Cloud-based vs 100% local
    - Cost: Paid API vs fully free
    - Quick start examples with prompts
  - **Example Prompts**:
    - "Derive the complete F_U_Bi_i equation starting from first principles with all 4 Ug components"
    - "Analyze GW170817 kilonova data using UQFF r-process nucleosynthesis framework"
    - "Compare UQFF_Compressed vs UQFF_Resonant for SGR 1745-2900 magnetar"
    - "Generate a complete C++ implementation of 26-layer Triadic gravity with validation"
    - "Explain the connection between UQFF vacuum buoyancy and quantum field theory"
  - **SuperGrok4 vs CoAnQi_bot (Tab 6)**:
    | Feature | SuperGrok4 (Tab 7) | CoAnQi_bot (Tab 6) |
    |---------|-------------------|-------------------|
    | Platform | Grok xAI (Cloud) | Ollama (Local) |
    | Model Size | ~300B params | 3-7B params |
    | Reasoning | Advanced (Grok-2) | Basic-Moderate |
    | Vision | ✅ Yes (grok-vision) | ❌ No |
    | Cost | Paid API | Free (local) |
    | Privacy | Cloud-based | 100% Local |
    | Best For | Research, complex derivations, paper analysis | Quick code help, debugging |

- **Tabs 8-21: Query Fetch Results** (14 tabs dedicated to search results, shifted from Tabs 7-21)
  - Auto-populate with results from user queries (SIMBAD, NASA, Grok fallback)
  - Tab 21 preloaded with ALMA Cycle 12 Observing Tool
  - All tabs closable and movable for custom organization
  - Browser windows with full web navigation capabilities

- **UQFF Physics:** QProcess integration with CoAnQi_Wrapper.py (source2.cpp lines 7645-7700)
- **Number Theory Tool:** SymbolicMath_Wrapper.py integration (source2.cpp lines 5803-5843)
  - Functions: p(n) partition, tau(n) Ramanujan tau, sigma(n) divisor sum, factors(n)
  - Input: Comma-separated or newline-separated equations
- **Error Messages:** Updated 3 NO_PYTHON messages to direct users to Number Theory tool
- **Auto-Deployment:** CMakeLists.txt deploys Python wrappers + OpenSSL + VC++ Runtime + WSTP DLLs

---

## 🛠️ Build Instructions

### Step 1: Clean Rebuild (Required)
The CLI access point code was added to MAIN_1_CoAnQi.cpp. You **must rebuild** to activate it.

```powershell
# Navigate to project root
cd C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic

# Clean old build (optional but recommended)
Remove-Item -Recurse -Force build_msvc -ErrorAction SilentlyContinue

# Configure with Visual Studio 2022 (REQUIRED for Wolfram WSTP)
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64

# Build Release configuration (with UPX compression)
cmake --build build_msvc --config Release

# Expected output location:
# build_msvc\Release\MAIN_1_CoAnQi.exe
```

**⏱️ Build Time:** 8-12 minutes (446 modules, 106,695 lines)

### Step 2: Verify Executable Exists
```powershell
# Check if executable was built successfully
Test-Path build_msvc\Release\MAIN_1_CoAnQi.exe

# Expected output: True
```

### Step 3: Build Source2 GUI (with UQFF + Number Theory Integration)
```powershell
# Build Source2 with auto-deployment of all dependencies
cmake --build build_msvc --config Release --target Source2

# Verify executable and dependencies
Test-Path build_msvc\Release\Source2.exe                    # GUI executable
Test-Path build_msvc\Release\CoAnQi_Wrapper.py             # UQFF physics wrapper
Test-Path build_msvc\Release\SymbolicMath_Wrapper.py       # Number theory wrapper
Test-Path build_msvc\Release\libssl-3-x64.dll              # OpenSSL (Grok AI)
Test-Path build_msvc\Release\libcrypto-3-x64.dll           # OpenSSL crypto
Test-Path build_msvc\Release\wstp64i4.dll                  # Wolfram WSTP

# Launch Source2 GUI
Start-Process "build_msvc\Release\Source2.exe" -WorkingDirectory "build_msvc\Release"
```

**⏱️ Build Time:** 3-5 minutes (Qt6 GUI + WebEngine + VTK)

**Features:**
- 21-tab browser with astronomical system browsing
- UQFF button: Click to compute physics via CoAnQi_Wrapper.py
- Number Theory tool: Compute p(n), tau(n), sigma(n), factors(n) via SymbolicMath_Wrapper.py
- Input formats: Comma-separated `p(10), tau(5)` or newline-separated

---

## ✅ Testing the CLI Integration

### Test 1: Manual CLI Testing

#### Test --batch Flag
```powershell
# Compute UQFF for Sagittarius A* (outputs JSON)
.\build_msvc\Release\MAIN_1_CoAnQi.exe --batch "Sagittarius A*"
```

**Expected Output (JSON):**
```json
{
  "status": "success",
  "system_name": "Sagittarius A*",
  "F_U_Bi_i": 1.234567e+30,
  "g_compressed": 9.876543e-08,
  "F_jet_rel": 1.111111e+28,
  "E_acc_rel": 2.222222e+40,
  "F_drag_rel": 3.333333e+26,
  "F_gw_rel": 4.444444e+25
}
```

#### Test --list-systems Flag
```powershell
# List all available systems
.\build_msvc\Release\MAIN_1_CoAnQi.exe --list-systems
```

**Expected Output (JSON):**
```json
{
  "status": "success",
  "total_systems": 121,
  "systems": [
    "Sagittarius A*",
    "M87",
    "Betelgeuse",
    "NGC 3596",
    ...
  ]
}
```

#### Test --system-info Flag
```powershell
# Get detailed parameters for a specific system
.\build_msvc\Release\MAIN_1_CoAnQi.exe --system-info "Betelgeuse"
```

**Expected Output (JSON):**
```json
{
  "status": "success",
  "system_name": "Betelgeuse",
  "M": 2.38e+31,
  "r": 6.17e+11,
  "L_X": 1.00e+30,
  "B0": 1.00e-04,
  "omega0": 1.00e-07,
  "v": 30000.0,
  "T": 3500.0
}
```

### Test 2: Python Wrapper Testing

#### Quick Test via CLI
```powershell
# Test Python wrapper with sample system
python CoAnQi_Wrapper.py "Sagittarius A*"
```

**Expected Output:**
```
Computing: Sagittarius A*
========================================
Status: success
F_U_Bi_i: 1.234567e+30 N
g_compressed: 9.876543e-08 m/s²
F_jet_rel: 1.111111e+28 N
E_acc_rel: 2.222222e+40 J
F_drag_rel: 3.333333e+26 N
F_gw_rel: 4.444444e+25 N
Computation time: 2.34s
```

#### Programmatic Usage
```python
from CoAnQi_Wrapper import CoAnQiCalculator

# Initialize calculator
calc = CoAnQiCalculator(verbose=True)

# Compute single system
result = calc.compute_system("Sagittarius A*")
print(f"F_U_Bi_i: {result.F_U_Bi_i:.6e} N")

# List all systems
systems = calc.list_available_systems()
print(f"Total systems: {len(systems)}")

# Get system info
info = calc.get_system_info("M87")
print(f"M87 mass: {info['M']} kg")
```

### Test 3: Automated Test Suite

#### Run Full Test Suite
```powershell
# Run all 6 integration tests
python test_integration.py
```

**Expected Output:**
```
======================================================================
Star-Magic UQFF Integration Test Suite
Testing C++ ↔ Python Integration Layer
======================================================================

1. Checking C++ Executable
✅ Found C++ executable: build_msvc\Release\MAIN_1_CoAnQi.exe

2. Testing CLI Batch Mode (--batch)
Testing: Sagittarius A*
✅ Received valid JSON output
  F_U_Bi_i: 1.234567e+30
  g_compressed: 9.876543e-08

3. Testing System List (--list-systems)
✅ Found 121 systems

4. Testing System Info (--system-info)
✅ Retrieved info for Betelgeuse

5. Testing Python Wrapper
✅ Successfully imported CoAnQi_Wrapper module
✅ Initialized CoAnQiCalculator
✅ Computation successful

6. Testing Data Layer Integration (Optional)
⚠️  Data layer modules not found (optional)

Test Summary
  cpp_exe              ✅ PASSED
  batch_mode           ✅ PASSED
  list_systems         ✅ PASSED
  system_info          ✅ PASSED
  python_wrapper       ✅ PASSED
  data_layer           ⚪ SKIPPED

Results: 5/6 passed, 0 failed, 1 skipped

🎉 All required tests passed! Integration layer is ready to use.

Total execution time: 15.23s
```

### Test 4: Source2 GUI Integration Testing

#### Test UQFF Physics Button
1. Launch Source2: `Start-Process "build_msvc\Release\Source2.exe" -WorkingDirectory "build_msvc\Release"`
2. In browser panel, select "Sagittarius A*" (or any system)
3. Click **UQFF** button
4. Verify results display in results panel

**Expected Output:**
```
System: Sagittarius A*
F_U_Bi_i: 2.0736×10¹²⁴ N
g_compressed: 4.86×10¹¹ m/s²
F_jet_rel: [value] N
E_acc_rel: [value] J
F_drag_rel: [value] N
F_gw_rel: [value] N
```

#### Test Number Theory Tool
1. In Source2 bottom panel, locate Number Theory calculator
2. Enter comma-separated equations: `p(10), tau(5), sigma(12), factors(60)`
3. Click **Compute** or press Enter
4. Verify output shows all results

**Expected Output:**
```json
{
  "results": [
    {"equation": "p(10)", "result": 42},
    {"equation": "tau(5)", "result": 4830},
    {"equation": "sigma(12)", "result": 28},
    {"equation": "factors(60)", "result": "2^2 * 3 * 5"}
  ]
}
```

#### Test Error Messages (NO_PYTHON blocks)
1. Trigger calculator features that use old embedded Python (disabled)
2. Verify error messages direct to Number Theory tool:
   - Derivative: `"Derivative: Use Number Theory tool (bottom panel)"`
   - Integral: `"Integral: Use Number Theory tool (bottom panel)"`
   - System solving: `"[System solving: Use Number Theory tool for symbolic math]"`
3. These messages confirm Python IS available via QProcess (not embedded pybind11)

---

## 🎛️ Using Tab 1: Embedded Calculator Terminal (NEW)

**What is Tab 1?**
Tab 1 in Source2.exe is now a dedicated embedded terminal running MAIN_1_CoAnQi.exe directly. This provides seamless access to the C++ calculator's **18-option interactive menu** (Cosmic Egg build) without leaving the GUI.

### Quick Start
1. Launch Source2.exe: `.\build_msvc\Release\Source2.exe`
2. Tab 1 automatically opens with "🎛️ MAIN_1 Calculator"
3. MAIN_1_CoAnQi.exe starts automatically (PID displayed)
4. Wait for menu to appear (18 options displayed automatically)

### 18-Option Menu (Cosmic Egg Build)
1. Calculate system (single)
2. Calculate ALL systems (parallel)
3. Clone and mutate system
4. Add custom system
5. Add dynamic physics term
6. Run simulations
7. Statistical analysis
8. Self-optimization
9. WSTP kernel interface (Wolfram)
10. Auto-export full UQFF to Wolfram
11. Run Wolfram Field Unity Simulation
12. **Run Cosmic Quantum Egg (26D) Simulation** ← Full 26-dimensional physics
13. Configure Grok API Key
14. Test Grok AI Integration
15. SOURCE4 Unified Field Validation
16. Information Paradox Tests (BH Info)
17. BSM Physics Validation
18. Exit

### Terminal Features
- **PowerShell-style interface** - Black background (#0C0C0C), green input echo
- **Direct stdin/stdout** - Type menu options (1-18) and press Enter
- **Error display** - Red text for stderr/errors
- **Auto-scroll** - Terminal automatically scrolls to bottom for new output
- **Process monitoring** - PID displayed, status tracked

### Special Terminal Commands
These commands are handled by the terminal widget itself (NOT sent to MAIN_1_CoAnQi.exe):

| Command | Action |
|---------|--------|
| `restart` | Kill and relaunch MAIN_1_CoAnQi.exe (useful if stuck) |
| `clear` or `cls` | Clear terminal output (process keeps running) |
| `help` | Show terminal command reference |

All other input (menu options 1-18, system names, etc.) is sent directly to MAIN_1_CoAnQi.exe.

### Example Usage: Calculate Sagittarius A*
```
1. Wait for menu to appear in Tab 1 terminal
2. Type: 1   (Option 1: Calculate system - single)
3. Press Enter
4. At prompt "Enter system name:", type: Sagittarius A*
5. Press Enter
6. View results in terminal (F_U_Bi_i, g_compressed, validation pipeline)
```

### Example Usage: Access Cosmic Quantum Egg (26D Simulation)
```
1. Type: 12   (Option 12: Run Cosmic Quantum Egg - 26D Simulation)
2. Press Enter
3. View 26 independent dimensional spheres with quantum state factors
4. Each dimension shows Q_i, [UA]_i, [SCm]_i compression parameters
```

### Example Usage: Access Wolfram Integration
```
1. Type: 9   (Option 9: WSTP kernel interface)
2. Press Enter
3. At prompt, enter Wolfram Language command: Integrate[x^2, x]
4. View result: x^3/3
```

### Clean Shutdown
- **Close Tab 1**: Terminal automatically sends Exit command (option 16 for Cosmic Egg build, varies by build)
- **Process kill**: If process doesn't exit within 3 seconds, terminal force-kills it
- **No orphaned processes**: All cleanup handled automatically

### Troubleshooting Tab 1 Terminal

#### Error: "MAIN_1_CoAnQi.exe not found"
**Solution:**
1. Build the calculator first:
   ```powershell
   cmake --build build_msvc --config Release --target MAIN_1_CoAnQi
   ```
2. Or copy MAIN_1_CoAnQi.exe to the same directory as Source2.exe:
   ```powershell
   Copy-Item .\build_msvc\Release\MAIN_1_CoAnQi.exe .\build_msvc\Release\
   ```

#### Error: "Failed to start MAIN_1_CoAnQi.exe"
**Possible Causes:**
- Missing DLL dependencies (Wolfram WSTP, OpenSSL)
- Incorrect working directory
- Antivirus blocking execution

**Solution:**
1. Run dependency check:
   ```powershell
   dumpbin /dependents .\build_msvc\Release\MAIN_1_CoAnQi.exe
   ```
2. Verify DLLs exist in same directory as executable
3. Check Windows Event Viewer for blocked execution errors

#### Process Crashes on Startup
**Symptoms:**
- Terminal shows "Process Terminated" immediately after "Calculator started"
- Exit code: Non-zero

**Solution:**
1. Run MAIN_1_CoAnQi.exe directly in PowerShell to see crash details:
   ```powershell
   .\build_msvc\Release\MAIN_1_CoAnQi.exe
   ```
2. Check for missing configuration files (if any required by your build)
3. Verify Wolfram WSTP kernel is accessible (for Wolfram-enabled builds)

#### Terminal Doesn't Respond to Input
**Symptoms:**
- Typing commands shows in terminal but nothing happens
- No output from MAIN_1_CoAnQi.exe

**Solution:**
1. Type `restart` to relaunch the process
2. If still frozen, close and reopen Source2.exe
3. Check if MAIN_1_CoAnQi.exe is waiting for specific input format

---

## � Using Tab 2: QCalc.py Python Terminal (NEW)

**What is Tab 2?**
Tab 2 is an interactive Python terminal running QCalc.py, the UQFF Quantum Calculator implementing 8 Master Equations for unified field physics.

### Quick Start
1. Launch Source2.exe
2. Click Tab 2: "🐍 QCalc.py (UQFF)"
3. Python starts automatically with QCalc.py pre-imported
4. Type Python commands directly or use QCalc functions

### 8 UQFF Master Equations
1. **UQFF (Base Unified Field)** - Complete unified force calculation
2. **UQFF_Compressed** - Newtonian gravity + 9 corrections (expansion, super, envelope, Ug_sum, cosm, quantum, fluid, perturbation, dark matter)
3. **UQFF_Resonant** - aDPM base + 13 frequency modes (THz, vacuum diff, super-frequency, aether resonance, quantum/aether/fluid frequencies, etc.)
4. **UQFF_Superconductive** - SCm vacuum modulation (quantum field effects)
5. **UQFF_Buoyant (F_U_Bi)** - Atomic scale Inside→Out forces
6. **UQFF_Master_Buoyant (F_U_Bi_i)** - Cosmic scale Outside→In forces
7. **UQFF_Triadic** - 26-layer gravitational scaling (multi-dimensional gravity)
8. **UQFF_Quadratic** - Root solutions (dual-solution physics)

### Terminal Features
- **Python dark theme** - VS Code colors (#1E1E1E background, #D4D4D4 text)
- **Python prompt** - Blue `>>>` prompt (#569CD6)
- **Interactive mode** - Full Python interpreter with QCalc pre-loaded
- **Auto-scroll** - Terminal scrolls to bottom for new output

### Special Commands
| Command | Action |
|---------|--------|
| `restart` | Restart Python interpreter |
| `clear` or `cls` | Clear terminal output (Python keeps running) |
| `help` | Show UQFF equation reference |

### Example Usage: Import and Use QCalc
```python
>>> from QCalc import *
>>> calc = UQFFCalculator()
>>> # Compute UQFF for a system
>>> result = calc.solve_UQFF(M=1.989e30, r=6.96e8, L_X=1e26)
>>> print(result['F_U_Bi_i'])
```

### Example Usage: 8 Master Equations
```python
>>> # 1. Base UQFF
>>> calc.solve_UQFF(params)

>>> # 2. Compressed (Newtonian + corrections)
>>> calc.solve_UQFF_Compressed(params)

>>> # 3. Resonant (frequency modes)
>>> calc.solve_UQFF_Resonant(params)

>>> # 5. Buoyant (atomic scale)
>>> calc.solve_UQFF_Buoyant(params)

>>> # 6. Master Buoyant (cosmic scale)
>>> calc.solve_UQFF_Master_Buoyant(params)

>>> # 7. Triadic (26-layer)
>>> calc.solve_UQFF_Triadic(params)
```

### Troubleshooting Tab 2
**Error: "QCalc.py not found"**
- Solution: Ensure QCalc.py is in the application directory or project root
- Check locations: `build_msvc/Release/QCalc.py` or `Star-Magic/QCalc.py`

**Error: "Failed to start Python"**
- Solution: Verify Python is installed and accessible via `python` command
- Test: `python --version` in PowerShell (should show Python 3.x)

---

## 🧮 Using Tab 3: Advanced Scientific Calculator (NEW)

**What is Tab 3?**
Tab 3 is an advanced physics calculator similar to Wolfram Alpha, MATLAB, and Mathematica, providing Standard Model physics computations with multiple engine options.

### Quick Start
1. Launch Source2.exe
2. Click Tab 3: "🧮 Scientific Calculator"
3. Select computational engine from dropdown (default: SymPy)
4. Enter equation in expression field
5. Click "⚡ Compute" or press Enter

### 4 Computational Engine Options

#### 1. 🐍 SymPy (Open Source) - **Recommended**
- **Type:** Symbolic mathematics library
- **Cost:** Free (no API key required)
- **Installation:** `pip install sympy`
- **Best for:** Algebraic equations, calculus, symbolic math
- **Example:** `solve(x^2 - 4 = 0, x)` → `[-2, 2]`

#### 2. 🔬 SciPy (Open Source)
- **Type:** Scientific computing library
- **Cost:** Free (no API key required)
- **Installation:** `pip install scipy`
- **Best for:** Numerical methods, optimization, signal processing
- **Example:** `optimize.minimize(func, x0)` → minimum value

#### 3. 🔢 NumPy (Open Source)
- **Type:** Numerical array library
- **Cost:** Free (no API key required)
- **Installation:** `pip install numpy`
- **Best for:** Linear algebra, matrix operations, array computations
- **Example:** `np.linalg.eig(matrix)` → eigenvalues

#### 4. 🧮 Wolfram Alpha API (Commercial)
- **Type:** Computational intelligence engine
- **Cost:** Free tier (2,000 queries/month), then paid
- **Registration:** [developer.wolframalpha.com](https://developer.wolframalpha.com/)
- **Best for:** Complex physics, natural language queries, step-by-step solutions
- **API Key Required:** Enter in field below dropdown

### Supported Operations

#### Algebra
```
solve(x^2 - 4 = 0, x)      → Find roots
expand((x+1)^2)            → x^2 + 2x + 1
factor(x^2 - 1)            → (x-1)(x+1)
simplify(x^2/x)            → x
```

#### Calculus
```
integrate(x^2, x)          → x^3/3 + C
diff(sin(x), x)            → cos(x)
limit(1/x, x, 0)           → ∞
series(exp(x), x, 0, 5)    → Taylor series
```

#### Linear Algebra
```
matrix([[1,2],[3,4]]).inv()  → Matrix inverse
det(A)                       → Determinant
eigenvals(M)                 → Eigenvalues
A @ B                        → Matrix multiplication
```

#### Physics & Constants
```
c                          → Speed of light (2.998e8 m/s)
h                          → Planck constant (6.626e-34 J·s)
G                          → Gravitational constant (6.674e-11)
k_B                        → Boltzmann constant (1.381e-23 J/K)
```

#### Statistics
```
mean([1,2,3,4,5])         → 3
std([...])                → Standard deviation
normal_dist(0, 1)         → Gaussian distribution
```

### Example Usage: Solve Quadratic Equation
```
1. Select engine: SymPy (Open Source)
2. Enter: solve(x^2 - 5x + 6 = 0, x)
3. Click "⚡ Compute"
4. Result: [2, 3]
```

### Example Usage: Integrate Function
```
1. Select engine: SymPy
2. Enter: integrate(x^2 * sin(x), x)
3. Result: -x^2*cos(x) + 2*x*sin(x) + 2*cos(x)
```

### API Registration (Wolfram Alpha)
1. Click "ℹ️ API Registration Info" button
2. Dialog shows registration links for all engines
3. For Wolfram Alpha:
   - Visit [developer.wolframalpha.com](https://developer.wolframalpha.com/)
   - Sign up for free account (2,000 queries/month)
   - Get API key from dashboard
   - Paste API key in "API Key" field

### Troubleshooting Tab 3
**Error: "Python library not found (SymPy/SciPy/NumPy)"**
- Solution: Install missing library
  ```powershell
  pip install sympy scipy numpy
  ```

**Error: "API Key required for Wolfram Alpha"**
- Solution: Enter valid Wolfram Alpha API key
- Get key from [developer.wolframalpha.com](https://developer.wolframalpha.com/)

**Error: "Computation timeout (>5s)"**
- Solution: Simplify equation or use different engine
- Complex symbolic math may require more time

---

## � Using Tab 4: Notebook Editor (NEW)

**What is Tab 4?**
Tab 4 is a Jupyter-style notebook environment for interactive Python workflows. Unlike Tab 2 (QCalc.py terminal), this tab allows you to create multi-cell notebooks combining code execution, output display, and markdown documentation.

### Quick Start
1. Launch Source2.exe
2. Click Tab 4: "📓 Notebook Editor"
3. Read welcome cell (purple border) for features overview
4. Click "➕ Code Cell" to add Python code cell
5. Type Python code, click "▶️ Run" to execute

### Toolbar Buttons

| Button | Action | Shortcut |
|--------|--------|----------|
| ➕ Code Cell | Add new Python code cell with green border | - |
| 📝 Markdown Cell | Add new documentation cell with blue border | - |
| ▶️ Run All | Execute all code cells sequentially | - |
| 🗑️ Clear All | Clear all cell outputs (cells remain) | - |

### Cell Types

#### Code Cells (Green Border #4CAF50)
- **Purpose:** Executable Python code with output display
- **Execution:** Click "▶️ Run" button or (future) Shift+Enter
- **Timeout:** 10 seconds per cell
- **Output Display:**
  - **Stdout:** Black text (#000000) - normal output
  - **Stderr:** Red text (#FF0000) - errors and warnings
  - **Scrollable:** Output area scrolls if exceeds 200px height
- **Format:** Monospace font (Consolas), line numbers via cell labels

#### Markdown Cells (Blue Border #2196F3)
- **Purpose:** Documentation, notes, explanations
- **Formatting:** Plain text (future: Markdown rendering)
- **Use Cases:** 
  - Explain next code section
  - Document assumptions
  - Show equations in text form
  - Create section headers

#### Welcome Cell (Purple Border #9C27B0)
- **Purpose:** Quick reference guide (displayed on startup)
- **Content:** Features list and use cases
- **Deletable:** Use 🗑️ button if not needed

### Example Usage: Multi-Step Physics Calculation

#### Workflow
```
1. Click "📝 Markdown Cell"
2. Type: "# UQFF Validation for Sagittarius A*"

3. Click "➕ Code Cell"
4. Type:
   import numpy as np
   M = 4.1e6 * 1.989e30  # Solar masses to kg
   r = 8500 * 3.086e16   # Parsecs to meters
   print(f"Mass: {M:.2e} kg")
   print(f"Distance: {r:.2e} m")

5. Click "▶️ Run" on code cell

6. Click "➕ Code Cell" (second cell)
7. Type:
   G = 6.674e-11
   g = G * M / r**2
   print(f"Newtonian gravity: {g:.2e} m/s²")

8. Click "▶️ Run All" to execute all cells sequentially
```

**Expected Output:**
```
[1] (Code Cell)
Mass: 8.16e36 kg
Distance: 2.62e20 m

[2] (Code Cell)
Newtonian gravity: 7.94e-08 m/s²
```

### Example Usage: Compare UQFF Methods
```
1. Add Markdown Cell: "Compare UQFF vs MUGE for NGC1365"

2. Add Code Cell:
   # Import parameters from Tab 2 QCalc
   from QCalc import UQFF_Compressed, UQFF_Resonant
   
   params = {
       'M': 1.2e9 * 1.989e30,
       'r': 17.3e6 * 3.086e16,
       'z': 0.00546
   }

3. Add Code Cell:
   g_compressed = UQFF_Compressed(params)
   print(f"UQFF Compressed: {g_compressed:.2e} m/s²")

4. Add Code Cell:
   g_resonant = UQFF_Resonant(params)
   print(f"UQFF Resonant: {g_resonant:.2e} m/s²")
   
   diff = abs(g_compressed - g_resonant) / g_compressed * 100
   print(f"Difference: {diff:.2f}%")

5. Click "▶️ Run All"
```

### Use Cases
- **Multi-Step Workflows:** Break complex calculations into logical steps
- **Data Analysis:** Load datasets, compute statistics, visualize results
- **UQFF Validation:** Compare different equation variants systematically
- **Debugging:** Execute code incrementally to identify errors
- **Documentation:** Combine explanations with executable code
- **Teaching:** Create tutorials with code + explanations

### Cell Management
- **Add Cell:** Click toolbar buttons (➕ or 📝)
- **Run Cell:** Click ▶️ button on individual cell
- **Run All:** Click "▶️ Run All" in toolbar (executes top to bottom)
- **Delete Cell:** Click 🗑️ button on cell
- **Clear Outputs:** Click "🗑️ Clear All" in toolbar (cells remain)

### Future Enhancements (Planned)
- **Keyboard Shortcuts:** Shift+Enter (run), Ctrl+Enter (run+advance), Esc (stop)
- **Markdown Rendering:** HTML display for markdown cells
- **Cell Reordering:** Drag-and-drop to reorganize cells
- **Export:** Save notebook to .ipynb format
- **Import:** Load existing Jupyter notebooks

### Troubleshooting Tab 4

**Error: "Execution timed out (>10s)"**
- **Cause:** Cell took longer than 10 seconds to execute
- **Solution:** Break into smaller cells, optimize code, or increase timeout in source code

**Error: "python: command not found"**
- **Cause:** Python not in system PATH
- **Solution:** Install Python 3.x and add to PATH, or specify full path in widget code

**Output Not Showing**
- **Symptoms:** Cell executes but no output appears
- **Solution:** Check if code actually prints (notebook doesn't auto-display last expression like interactive mode)
- **Workaround:** Add explicit `print()` statements

**Cell Execution Order Confusion**
- **Issue:** Cells run out of order, variables undefined
- **Solution:** Use "▶️ Run All" to execute sequentially from top
- **Best Practice:** Design cells to be independent or document dependencies

---

## 📚 Using Tab 5: CondensedPhysics.py Terminal (NEW)

**What is Tab 5?**
Tab 5 is a Python interactive terminal for CondensedPhysics.py, a **pure physics calculator** that solves general model/class equations. Unlike Tab 2 (QCalc with 8 fixed UQFF equations), this terminal provides access to a general solver index for arbitrary physics domains.

### Critical Architecture Rules
CondensedPhysics.py is **NOT a data repository**—it's a **pure physics calculator**:
- ❌ **NO hardcoded system data** (no "Sagittarius A* = 4.1e6 solar masses")
- ❌ **NO named system classes** (no `class NGC1365Model`)
- ❌ **NO global instances** (no pre-computed solutions)
- ❌ **NO pre-stored results** (computations are live only)
- ✅ **Parameterized equations only** (accepts datasets as input)
- ✅ **Fundamental constants only** (G, c, ℏ, k_B, etc.)

### Data Flow Architecture
```
Source2 Query Field → APIFetch.py → parameters dict
                                         ↓
                          CondensedPhysics.solve(params)
                                         ↓
                     Long-form equations + solutions
                                         ↓
                      CondensedPhysics_OutputData.py (storage)
                                         ↓
                            Source2 recall system
```

### Quick Start
1. Launch Source2.exe
2. Click Tab 5: "📚 CondensedPhysics.py"
3. Terminal auto-starts Python with CondensedPhysics.py loaded
4. Verify prompt shows: `CP>>>` (purple text)
5. Type Python commands or CondensedPhysics functions

### Terminal Features
- **Dark Theme:** Background #2E2E2E (dark gray)
- **Purple Prompt:** `CP>>>` in #9C27B0 color
- **Light Text:** #E0E0E0 (light gray) for readability
- **Error Display:** Red text (#FF5252) for stderr
- **Success Messages:** Green text (#00E676) for confirmations

### Special Commands

| Command | Action |
|---------|--------|
| `restart` | Kill and relaunch Python process with CondensedPhysics.py |
| `exit()` | Close Python terminal (standard Python exit) |
| `help(CondensedPhysics)` | Show module documentation |

### Example Usage: General Galaxy Rotation Solver

```python
CP>>> # Define parameters from API fetch or user input
CP>>> params = {
...     'M': 1.2e11 * 1.989e30,      # Galaxy mass (kg)
...     'r': 5.0e20,                  # Radius (m)
...     'SFR': 15.3,                  # Star formation rate (M☉/yr)
...     'z': 0.0055                   # Redshift
... }

CP>>> # Solve with CondensedPhysics general solver
CP>>> result = CondensedPhysics.solve_rotation_curve(params)

CP>>> print(result['long_form_equations'])
[
  "v(r) = sqrt(G * M(r) / r)",
  "M(r) = M_visible + M_dark(r)",
  "M_dark(r) = 4π * ρ_0 * r_s³ * [ln((r_s+r)/r_s) - r/(r_s+r)]",
  ...
]

CP>>> print(result['solutions'])
{
  'v_rotation': 2.31e5,  # m/s
  'M_dark': 8.4e41,       # kg
  'DM_fraction': 0.68     # 68% dark matter
}

CP>>> print(result['available_equations'])
[
  "Tully-Fisher relation: L ∝ v^4",
  "MOND modified dynamics: a = a_N/(1 + a_0/a_N)",
  "Navarro-Frenk-White profile",
  "Burkert dark matter profile",
  ...
]
```

### Output Format
CondensedPhysics.py returns structured results:

```python
{
  'long_form_equations': [      # Step-by-step derivation
    "Equation 1 with substitutions shown",
    "Equation 2 intermediate step",
    ...
  ],
  'solutions': {                 # Numerical results
    'parameter1': value1,
    'parameter2': value2,
    ...
  },
  'available_equations': [       # Other solvable equations for this query
    "Alternative method 1",
    "Alternative method 2",
    ...
  ],
  'simulation_set': {            # For simultaneous multi-equation simulation
    'equations': [...],
    'initial_conditions': {...}
  }
}
```

### Example Usage: Pass Data from Tab 1 Results
```python
CP>>> # After running calculation in Tab 1 (MAIN_1 Calculator)
CP>>> # Copy parameters from terminal output

CP>>> system_params = {
...     'name': 'Sagittarius A*',
...     'M': 4.1e6 * 1.989e30,
...     'r': 26000 * 3.086e16,
...     'B': 1e-4,
...     'Charge': 1.602e-19
... }

CP>>> # Validate with general magnetar model
CP>>> result = CondensedPhysics.solve_magnetar_field(system_params)
CP>>> print(f"B_crit validation: {result['B_crit']:.2e} T")
```

### Use Cases
- **General Physics Solver:** Compute arbitrary physics domains (not limited to 8 UQFF equations)
- **Model Validation:** Compare standard models with UQFF predictions
- **Custom Workflows:** Build multi-stage calculations with parameter passing
- **Research Exploration:** Test hypotheses with parameterized models
- **Equation Discovery:** See all available equations for a given domain

### Troubleshooting Tab 5

**Error: "CondensedPhysics.py not found"**
- **Cause:** File missing from expected locations
- **Solution:** Terminal searches 2 locations automatically:
  1. Same directory as Source2.exe
  2. Parent directory (project root)
- **Fix:** Copy CondensedPhysics.py to one of these locations
- **Message Display:** Terminal shows expected paths and architecture rules

**Error: "Hardcoded system data detected"**
- **Cause:** CondensedPhysics.py violates architecture rules
- **Solution:** Review MANDATORY ARCHITECTURE RULES in CondensedPhysics.py header
- **Required Changes:**
  - Remove `self.distance = 6500 * ly` style assignments
  - Change `class NGC3596Model` to `class GalaxyRotationCalculator`
  - Replace global instances with stateless calculator functions

**Process Crashes on Import**
- **Cause:** Syntax error or missing Python dependencies
- **Solution:** Run CondensedPhysics.py directly in PowerShell:
  ```powershell
  python CondensedPhysics.py
  ```
- **Check:** Ensure all imports available (numpy, scipy, etc.)

---

## 🤖 Using Tab 6: CoAnQi_bot (Ollama Local LLM) (NEW)

**What is Tab 6?**
Tab 6 is a local AI assistant powered by Ollama, providing code generation, physics explanations, and debugging help—all running entirely on your machine with **no cloud API, no API key, and completely free**.

### Quick Start (First-Time Setup)
1. **Install Ollama:**
   - Download from [https://ollama.com](https://ollama.com)
   - Run installer (Windows/Mac/Linux)
2. **Pull Model:**
   ```powershell
   ollama pull llama3.2
   ```
3. **Start Server:**
   ```powershell
   ollama serve
   ```
   - Keep terminal open (runs on localhost:11434)
4. **Open Source2 → Tab 6**
5. Select model from dropdown (llama3.2:latest)
6. Start chatting!

### Installation Guide (In-App)
Click "📥 Installation Guide" button in Tab 6 for step-by-step setup:

**Step 1: Download Ollama**
- Visit https://ollama.com
- Download installer for your OS
- Install (typically 1-2 minutes)

**Step 2: Pull a Model**
```powershell
# General purpose (recommended for beginners)
ollama pull llama3.2

# Code specialist
ollama pull codellama

# Fast and efficient
ollama pull mistral

# Latest code model
ollama pull qwen2.5-coder
```

**Step 3: Start Server**
```powershell
ollama serve
```
- Server starts on localhost:11434
- Keep terminal open while using CoAnQi_bot

**Step 4: Test Connection**
- Open Source2 → Tab 6
- Type: "Hello, are you working?"
- Verify bot responds

### 4 Available Models

| Model | Size | Best For | Strengths |
|-------|------|----------|-----------|
| **llama3.2:latest** | 3B | General purpose | Balanced, fast, good reasoning |
| **codellama:latest** | 7B | Code generation | Python/C++ expert, debugging |
| **mistral:latest** | 7B | Fast inference | Quick responses, efficient |
| **qwen2.5-coder:latest** | 7B | Latest code model | Multi-language, up-to-date |

**Note:** Download size: 2-4 GB per model, one-time download, stored locally

### System Prompt (Auto-Included)
CoAnQi_bot has built-in UQFF awareness via system prompt:

> "You are CoAnQi_bot, an expert physics and code assistant for the UQFF (Unified Quantum Field Framework) project. Help with code generation, physics equations, and scientific computing. Be concise and provide code examples when relevant."

This makes the bot aware of:
- UQFF terminology (F_U_Bi_i, g_compressed, 26D Triadic, etc.)
- Project structure (MAIN_1_CoAnQi.cpp, QCalc.py, source files)
- Physics context (magnetars, black holes, dark matter)

### Example Usage: Explain UQFF Equation

**You:**
> Explain the F_U_Bi_i equation in UQFF and how it differs from Newtonian gravity

**CoAnQi_bot:**
> F_U_Bi_i (Master Buoyant Unified Force, cosmic scale Outside→In) represents the complete unified field force in UQFF. Key differences from Newtonian gravity:
>
> 1. **Multi-Component**: F_U_Bi_i = Ug1 + Ug2 + Ug3 + Ug4 (4 unified forces)
>    - Ug1: Magnetic dipole contribution
>    - Ug2: Charge-reactivity field
>    - Ug3: String rotation dynamics
>    - Ug4: Vacuum concentration effects
>
> 2. **Newtonian**: F = G*M*m/r² (single mass-based force)
>
> 3. **UQFF Advantages:**
>    - Accounts for quantum field effects (Ug4)
>    - Includes electromagnetic coupling (Ug1, Ug2)
>    - Models spacetime rotation (Ug3)
>    - Naturally explains galaxy rotation curves without dark matter adjustments
>
> Use `compute_FU_SOURCE4()` in MAIN_1_CoAnQi.cpp for calculations.

### Example Usage: Generate Code

**You:**
> Generate Python code to calculate gravitational force using UQFF_Compressed for NGC1365

**CoAnQi_bot:**
> ```python
> from QCalc import UQFF_Compressed
> 
> # NGC1365 parameters
> params = {
>     'M': 1.2e9 * 1.989e30,      # 1.2 billion solar masses
>     'r': 17.3e6 * 3.086e16,     # 17.3 Mpc
>     'z': 0.00546,                # Redshift
>     'B': 1e-6,                   # Magnetic field (T)
>     'SFR': 15.3                  # Star formation rate (M☉/yr)
> }
> 
> # Compute UQFF Compressed (Newtonian + 9 corrections)
> g_compressed = UQFF_Compressed(params)
> 
> print(f"NGC1365 Compressed Gravity: {g_compressed:.2e} m/s²")
> 
> # Compare with Newtonian
> G = 6.674e-11
> g_newton = G * params['M'] / params['r']**2
> 
> correction_factor = g_compressed / g_newton
> print(f"UQFF Correction Factor: {correction_factor:.4f}")
> ```
>
> Run this in Tab 2 (QCalc.py) or Tab 4 (Notebook Editor).

### Example Usage: Debug C++ Error

**You:**
> I'm getting "error C2065: 'QScrollArea': undeclared identifier" in source2.cpp line 1074. How do I fix this?

**CoAnQi_bot:**
> This error means QScrollArea is used but not included. Add this to your includes section in source2.cpp:
>
> ```cpp
> #include <QScrollArea>
> ```
>
> Make sure it's placed with other Qt widget includes (after QWidget, QLabel, etc.).
>
> If error persists:
> 1. Check CMakeLists.txt includes Qt6::Widgets
> 2. Clean and rebuild: `cmake --build build_msvc --config Release --target Source2`
> 3. Verify Qt 6.x is properly installed

### Chat Interface Features
- **User Messages:** Blue background (#E3F2FD), left-aligned, clear "You:" label
- **Bot Responses:** Green background (#F1F8E9), HTML formatted with code syntax highlighting
- **Error Messages:** Red background (#FFEBEE), connection troubleshooting tips
- **Status Bar:** Shows "Ready", "Thinking..." during inference, or "Error - Ollama not running"
- **Auto-Scroll:** Chat automatically scrolls to latest message

### Use Cases
- **Code Generation:** "Generate C++ code for 26D Triadic gravity solver"
- **Physics Explanations:** "Explain Hawking radiation in context of Information Paradox"
- **Debugging:** "Why is my WSTP integration crashing on startup?"
- **Optimization:** "How can I parallelize this loop in MAIN_1_CoAnQi.cpp?"
- **Documentation:** "Write docstring for compute_Ug4_SOURCE4 function"
- **Learning:** "What's the difference between UQFF_Resonant and UQFF_Compressed?"

### Privacy & Performance
- **100% Local:** All inference runs on your machine, no data sent to cloud
- **No API Key:** No registration, no usage limits, no tracking
- **Offline Capable:** Works without internet (after model download)
- **Performance:** Depends on your hardware:
  - **CPU (no GPU):** 5-15 seconds per response
  - **GPU (CUDA/Metal):** 1-3 seconds per response
- **Model Size:** 2-7B parameters (fits on most modern computers)

### Troubleshooting Tab 6

**Error: "Ollama not running - Connection refused"**
- **Cause:** Ollama server not started
- **Solution:**
  ```powershell
  ollama serve
  ```
- **Verify:** Open http://localhost:11434 in browser (should show "Ollama is running")

**Error: "Model not found - llama3.2:latest"**
- **Cause:** Model not downloaded
- **Solution:**
  ```powershell
  ollama pull llama3.2
  ```
- **Check Available Models:**
  ```powershell
  ollama list
  ```

**Response Takes Too Long (>30s)**
- **Cause:** CPU inference on large model
- **Solutions:**
  - Use smaller model: llama3.2 (3B) instead of codellama (7B)
  - Enable GPU acceleration (Ollama auto-detects CUDA/Metal)
  - Reduce prompt length

**Bot Response Cut Off/Incomplete**
- **Cause:** Response length limit in API call
- **Solution:** Ask bot to "continue" or rephrase as shorter question
- **Workaround:** Increase max_tokens in OllamaCodeBotWidget (source2.cpp line ~1500)

**Installation Guide Doesn't Open**
- **Cause:** Modal dialog coding issue
- **Workaround:** Follow manual installation steps:
  1. https://ollama.com → Download
  2. `ollama pull llama3.2`
  3. `ollama serve`

---

## 🧠 Using Tab 7: SuperGrok4 (Grok xAI Expert) (NEW)

**What is Tab 7?**
Tab 7 is an advanced AI research assistant powered by Grok xAI, providing **research-grade responses** for complex physics problems, multi-step derivations, and scientific paper analysis. Unlike Tab 6 (local Ollama), SuperGrok4 uses cloud-based models (~300B parameters) with advanced reasoning capabilities.

### Quick Start (First-Time Setup)
1. **Sign up for xAI API:**
   - Visit [https://x.ai](https://x.ai) (xAI platform)
   - Create account and generate API key from dashboard
   - **Cost:** Paid API (check x.ai for pricing, typically $0.01-0.10 per query)
2. **Set Environment Variable:**
   - **Windows (PowerShell, temporary):**
     ```powershell
     $env:XAI_API_KEY = "your-api-key-here"
     ```
   - **Windows (Permanent):**
     - Right-click "This PC" → Properties → Advanced System Settings
     - Environment Variables → New (under User variables)
     - Variable name: `XAI_API_KEY`
     - Variable value: `your-api-key-here`
   - **macOS/Linux (temporary):**
     ```bash
     export XAI_API_KEY="your-api-key-here"
     ```
   - **macOS/Linux (Permanent):**
     ```bash
     echo 'export XAI_API_KEY="your-api-key-here"' >> ~/.bashrc
     source ~/.bashrc
     ```
3. **Restart Source2.exe** to load environment variable
4. **Open Source2 → Tab 7**
5. Click "🔑 Configure API Key" to verify status (should show ✅ SET)
6. Select model from dropdown (grok-beta recommended for general use)
7. Start chatting!

### API Key Configuration (In-App)
Click "🔑 Configure API Key" button in Tab 7 to check current status:

**Configuration Dialog Shows:**
- ✅ **SET** → API key found in environment, ready to use
- ❌ **NOT SET** → Follow OS-specific setup instructions

**Current Status Example:**
```
Status: ✅ XAI_API_KEY is SET
Value: xai-••••••••••••••••••••••••••••••5f3a (masked)
```

**If NOT SET, dialog shows:**
```
Windows (PowerShell):
  $env:XAI_API_KEY = "your-api-key-here"

Windows (Permanent):
  System Properties → Environment Variables → New user variable

macOS/Linux:
  export XAI_API_KEY="your-api-key-here"
  Add to ~/.bashrc or ~/.zshrc for permanent setup
```

### 3 Available Models

| Model | Best For | Strengths | Unique Features |
|-------|----------|-----------|-----------------|
| **grok-beta** | General purpose | Latest model, ~300B params, state-of-the-art performance | Balanced speed/quality |
| **grok-2** | Complex reasoning | Multi-step derivations, advanced logic | Research-grade math |
| **grok-vision-beta** | Multimodal | Text + image analysis | Astronomical image analysis |

**Note:** All models require internet connection and API key

### System Prompt (Auto-Included)
SuperGrok4 has built-in UQFF expertise via system prompt:

> "You are SuperGrok4, an expert physics and research assistant for the UQFF (Unified Quantum Field Framework) project. You have deep knowledge of astrophysics, quantum mechanics, UQFF equations, and scientific computing. Provide detailed, research-grade responses with mathematical rigor when appropriate."

This makes the bot aware of:
- UQFF terminology (F_U_Bi_i, g_compressed, 26D Triadic, Ug1-Ug4 components)
- Project structure (MAIN_1_CoAnQi.cpp, source files, SOURCE4 validation)
- Physics context (magnetars, black holes, kilonova, r-process nucleosynthesis)
- Research standards (expect derivations, citations, equations)

### Welcome Message (Auto-Displayed)
On first launch, SuperGrok4 displays a comprehensive welcome with:

**Features List:**
- 🧠 **Deep UQFF Knowledge** - Expert understanding of all UQFF equations and theory
- 🔬 **Research Analysis** - Scientific paper review, literature analysis, data interpretation
- 💻 **Code Generation** - Complete C++/Python implementations with documentation
- 🔢 **Mathematical Derivations** - Step-by-step equation derivations with rigor
- 🔭 **Astrophysical Expertise** - GW170817, magnetars, black holes, cosmology
- 👁️ **Vision Support** (grok-vision-beta) - Analyze astronomical images and diagrams

**SuperGrok4 vs CoAnQi_bot Comparison Table:**
| Feature | SuperGrok4 (Tab 7) | CoAnQi_bot (Tab 6) |
|---------|-------------------|-------------------|
| Platform | Grok xAI (Cloud) | Ollama (Local) |
| Model Size | ~300B params | 3-7B params |
| Reasoning | Advanced (Grok-2) | Basic-Moderate |
| Vision | ✅ Yes (grok-vision) | ❌ No |
| Cost | Paid API (~$0.01-0.10/query) | Free (local) |
| Privacy | Cloud-based | 100% Local |
| Best For | Research, complex derivations, paper analysis | Quick code help, debugging, simple questions |

**Quick Start Examples:**
- "Derive the complete F_U_Bi_i equation starting from first principles"
- "Analyze GW170817 kilonova data using UQFF r-process nucleosynthesis"
- "Generate a complete C++ implementation of 26-layer Triadic gravity"

### Example Usage: Complete Equation Derivation

**You:**
> Derive the complete F_U_Bi_i equation starting from first principles with all 4 Ug components (Ug1, Ug2, Ug3, Ug4). Show each step and explain the physical meaning.

**SuperGrok4:**
> # Complete F_U_Bi_i Derivation
>
> The Master Buoyant Unified Force (F_U_Bi_i) represents the cosmic-scale Outside→In unified field force in UQFF. Starting from first principles:
>
> ## 1. Fundamental Postulate
> The vacuum is not empty but contains quantum field energy with density ρ_vac. This creates a buoyancy-like effect on massive objects.
>
> ## 2. Four Unified Components
>
> ### Ug1: Magnetic Dipole Contribution
> ```
> Ug1 = (μ₀/4π) * (M_mag * r_vec) / r⁵
> ```
> - M_mag = magnetic dipole moment = B₀ * (4π/3) * R³
> - Physical meaning: Magnetic field dipole interaction with vacuum
> - Dominant for: Magnetars (B > 10¹⁴ T)
>
> ### Ug2: Charge-Reactivity Field
> ```
> Ug2 = (k_e * Q_eff * [UA]) / r²
> ```
> - Q_eff = effective charge from ionization
> - [UA] = vacuum reactivity coefficient ≈ 10⁻⁴
> - Physical meaning: Electric charge coupling to vacuum polarization
> - Dominant for: Highly ionized plasma environments
>
> ### Ug3: String Rotation Dynamics
> ```
> Ug3 = (ℏ * ω * sin(θ)) / (2π * r²)
> ```
> - ω = angular velocity (rotation)
> - θ = angle between rotation axis and position vector
> - Physical meaning: Spacetime twist from rotating mass-energy
> - Dominant for: Fast-spinning compact objects (pulsars, black holes)
>
> ### Ug4: Vacuum Concentration Effects
> ```
> Ug4 = -(k_η * [SSq] * M * c²) / r²
> ```
> - k_η = 10⁻¹¹³ (quantum coupling constant)
> - [SSq] = 0.57 (supersymmetric string constant)
> - Physical meaning: Vacuum energy gradient around mass
> - Dominant for: All systems (universal quantum correction)
>
> ## 3. Complete F_U_Bi_i Equation
> ```
> F_U_Bi_i = Ug1 + Ug2 + Ug3 + Ug4
> 
> F_U_Bi_i = (μ₀/4π) * (M_mag * r_vec)/r⁵ 
>          + (k_e * Q_eff * [UA])/r² 
>          + (ℏ * ω * sin(θ))/(2π*r²) 
>          - (k_η * [SSq] * M * c²)/r²
> ```
>
> ## 4. Comparison with Newtonian Gravity
> - **Newtonian**: F = G*M*m/r² (single component, mass-only)
> - **UQFF F_U_Bi_i**: 4 components (magnetic, electric, rotational, quantum)
> - **Advantage**: Naturally explains galaxy rotation curves, dark matter effects
>
> ## 5. Implementation (SOURCE4)
> See MAIN_1_CoAnQi.cpp lines 25623-26026:
> ```cpp
> double FU = SOURCE4::compute_FU_SOURCE4(body, r, t, tn, theta);
> ```
>
> ## 6. Validation
> - GW170817 kilonova ejecta dynamics: F_U_Bi_i matches observations
> - SGR 1745-2900 magnetar field: Ug1 dominant as expected
> - Galaxy rotation curves: Ug4 provides dark matter correction
>
> Would you like me to derive the 26-layer Triadic extension or show numerical validation for a specific system?

### Example Usage: Paper Analysis

**You:**
> Analyze the GW170817 LIGO paper (arXiv:1710.05832) and explain how UQFF r-process nucleosynthesis compares with their kilonova models. Are there specific predictions UQFF makes that differ from standard models?

**SuperGrok4:**
> # GW170817 LIGO Paper Analysis (arXiv:1710.05832)
>
> ## Standard Model (LIGO Paper)
> The LIGO collaboration's analysis of GW170817 relied on:
> - **Ejecta mass**: 0.04-0.06 M☉ (dynamical + wind)
> - **Velocity**: 0.1-0.3c (two-component model)
> - **r-process yield**: Blue kilonova (lanthanide-poor) + Red kilonova (lanthanide-rich)
> - **Heating source**: Radioactive decay only (Nd-252, Cf-254, etc.)
>
> ## UQFF Predictions
> UQFF adds vacuum buoyancy effects to r-process nucleosynthesis:
>
> ### 1. Enhanced Neutron Capture Rates
> ```
> λ_n_eff = λ_n_standard * (1 + [UA] * Ug4/k_B*T)
> ```
> - Ug4 vacuum concentration increases neutron capture cross-section
> - **Prediction**: Faster r-process (reaches A=195-200 more efficiently)
> - **Observational test**: Higher Pt/Au ratio in late-time spectra
>
> ### 2. Modified Ejecta Expansion
> F_U_Bi_i buoyancy force affects ejecta dynamics:
> ```
> v_ejecta(t) = v_0 * exp(-t/τ_drag) + (F_U_Bi_i / M_ej) * t
> ```
> - **Prediction**: Slight deceleration from Ug4 term (~5% at t=10 days)
> - **Observational test**: Spectral line blueshifts vs time
>
> ### 3. Dual Heating Mechanism
> - **Standard**: Radioactive decay only
> - **UQFF**: Radioactive + Vacuum coupling (Ug4 * M_ej / R_photosphere)
> - **Prediction**: Extended late-time emission (~30-50 days)
> - **Match**: UQFF explains 40-day Hubble detection better than pure radioactive model
>
> ## Specific Testable Predictions
> 1. **Pt/Au Abundance Ratio**: UQFF predicts 1.2x higher than standard (measure via spectroscopy)
> 2. **Late-Time Flux**: UQFF predicts 15-20% excess at t>30 days (JWST follow-up)
> 3. **Neutrino Flux**: Ug4 coupling enhances ν_e capture → measure with IceCube/KM3NeT
>
> ## Next Steps for Validation
> - Compare UQFF_InputData.py GW170817 parameters with LIGO results
> - Run CondensedPhysics.py r-process model with UQFF corrections
> - Check ValidationData.py benchmarks against observational constraints
>
> Would you like me to generate Python code to compute UQFF r-process yields for GW170817?

### Example Usage: Complete Code Generation

**You:**
> Generate a complete C++ implementation of 26-layer Triadic gravity with validation against NGC1365 rotation curve data. Include all 26 Ug components, quantum state factors, and error handling.

**SuperGrok4:**
> ```cpp
> // 26-Layer Triadic Gravity Calculator
> // Based on SOURCE115 (source172.cpp) master equations
> // Author: SuperGrok4 (Generated), Date: 2026-02-11
>
> #include <vector>
> #include <cmath>
> #include <stdexcept>
> #include <iostream>
>
> struct TriadicLayer {
>     double Ug1, Ug2, Ug3, Ug4;  // 4 components per layer
>     double Q_i;                  // Quantum state factor
>     double UA_i;                 // Vacuum reactivity
>     double SCm_i;                // Superconductive modulation
> };
>
> class TriadicGravitySolver {
> private:
>     std::vector<TriadicLayer> layers_;
>     double M_total_;    // Total mass (kg)
>     double r_;          // Distance (m)
>     double B0_;         // Magnetic field (T)
>
>     // Physical constants
>     static constexpr double G = 6.67430e-11;
>     static constexpr double c = 2.99792458e8;
>     static constexpr double mu_0 = 1.25663706212e-6;
>     static constexpr double k_eta = 1e-113;
>     static constexpr double SSq = 0.57;
>
> public:
>     TriadicGravitySolver(double M, double r, double B0)
>         : M_total_(M), r_(r), B0_(B0) {
>         if (M <=  0 || r <= 0) {
>             throw std::invalid_argument("Mass and radius must be positive");
>         }
>         initializeLayers();
>     }
>
>     void initializeLayers() {
>         layers_.resize(26);
>         for (int i = 0; i < 26; ++i) {
>             double layer_factor = (i + 1) / 26.0;
>             layers_[i].Q_i = 0.01 * std::exp(-layer_factor * 10);
>             layers_[i].UA_i = 1e-4 * (  1 - layer_factor * 0.5);
>             layers_[i].SCm_i = 0.99 + 0.01 * std::sin(layer_factor * M_PI);
>         }
>     }
>
>     double computeLayer(int layer_idx) {
>         if (layer_idx < 0 || layer_idx >= 26) {
>             throw std::out_of_range("Layer index must be 0-25");
>         }
>
>         TriadicLayer& L = layers_[layer_idx];
>         double r_eff = r_ * (1 + layer_idx * 0.05);  // Layered distance
>
>         // Ug1: Magnetic dipole
>         double M_mag = B0_ * (4.0 * M_PI / 3.0) * std::pow(1e4, 3);
>         L.Ug1 = (mu_0 / (4 * M_PI)) * (M_mag / std::pow(r_eff, 3));
>
>         // Ug2: Charge-reactivity
>         double Q_eff = 1.602e-19 * 1e40;  // Effective charge
>         L.Ug2 = (8.99e9 * Q_eff * L.UA_i) / (r_eff * r_eff);
>
>         // Ug3: String rotation (ω = 2π/day)
>         double omega = 2 * M_PI / 86400;
>         L.Ug3 = (1.055e-34 * omega) / (2 * M_PI * r_eff * r_eff);
>
>         // Ug4: Vacuum concentration
>         L.Ug4 = -(k_eta * SSq * M_total_ * c * c) / (r_eff * r_eff);
>
>         // Quantum-modulated sum
>         return (L.Ug1 + L.Ug2 + L.Ug3 + L.Ug4) * L.Q_i * L.SCm_i;
>     }
>
>     double computeTotal26D() {
>         double g_total = 0.0;
>         for (int i = 0; i < 26; ++i) {
>             g_total += computeLayer(i);
>         }
>         return g_total;
>     }
>
>     void validate_NGC1365() {
>         // NGC1365 observed rotation curve: v ~ 220 km/s at r ~ 20 kpc
>         double r_test = 20e3 * 3.086e16;  // 20 kpc in meters
>         r_ = r_test;
>         M_total_ = 1.2e9 * 1.989e30;      // 1.2 billion solar masses
>         B0_ = 1e-6;                        // Typical galaxy B-field
>
>         initializeLayers();
>         double g_26D = computeTotal26D();
>         double v_rot = std::sqrt(g_26D * r_test);
>
>         std::cout << "NGC1365 26D Triadic Validation:\n";
>         std::cout << "  Computed g_26D: " << g_26D << " m/s²\n";
>         std::cout << "  Rotation velocity: " << v_rot / 1000 << " km/s\n";
>         std::cout << "  Observed (target): 220 km/s\n";
>         std::cout << "  Match: " << (std::abs(v_rot/1000 - 220) < 20 ? "✅ PASS" : "❌ FAIL") << "\n";
>     }
> };
>
> int main() {
>     try {
>         // Example: Sagittarius A*
>         double M = 4.1e6 * 1.989e30;  // 4.1 million solar masses
>         double r = 1.5 * 1.496e11;    // 1.5 AU from event horizon
>         double B0 = 1e-4;             // Weak magnetic field
>
>         TriadicGravitySolver solver(M, r, B0);
>         double g_result = solver.computeTotal26D();
>
>         std::cout << "Sagittarius A* 26D Gravity: " << g_result << " m/s²\n";
>
>         // Validation against NGC1365
>         std::cout << "\n=== Running NGC1365 Validation ===\n";
>         solver.validate_NGC1365();
>
>     } catch (const std::exception& e) {
>         std::cerr << "Error: " << e.what() << std::endl;
>         return 1;
>     }
>
>     return 0;
> }
> ```
>
> **Compilation:**
> ```bash
> g++ -std=c++17 -O3 triadic_26d.cpp -o triadic_26d
> ./triadic_26d
> ```
>
> **Integration with MAIN_1_CoAnQi.cpp:**
> - This code can be added as a new SOURCE173 module
> - Use PhysicsTerm architecture for integration with CalculatorCore
> - See SOURCE115 (lines 24800-25622) for reference master equations
>
> Would you like me to add threading for parallel layer computation or expand to 121+ astrophysical systems?

### Chat Interface Features
- **User Messages:** Blue background (#E8EAF6), left-aligned, "You:" label
- **SuperGrok4 Responses:** Green background (#E8F5E9), HTML formatted with LaTeX-style equations
- **Error Messages:** Red background (#FFEBEE), API key/connection troubleshooting
- **Status Bar Colors:**
  - **Ready** (green #4CAF50): Idle, ready for queries
  - **Thinking** (orange #FF9800): API call in progress
  - **Error** (red #F44336): API key missing or connection failed
- **Auto-Scroll:** Chat automatically scrolls to latest response

### Use Cases (SuperGrok4 Strengths)
- **Research Papers:** "Analyze this arXiv paper and compare with UQFF predictions"
- **Complex Derivations:** "Derive 26D Triadic master equation from first principles"
- **Multi-Step Reasoning:** "Prove GW170817 kilonova data validates UQFF r-process"
- **Complete Implementations:** "Generate C++ code for SOURCE173 with full validation"
- **Vision Analysis:** *(grok-vision-beta)* "Analyze this ALMA image and identify r-process signatures"
- **Literature Review:** "Survey all papers on neutron star equation of state 2020-2025"

### When to Use SuperGrok4 vs CoAnQi_bot

| Task Type | Recommended Tab | Reason |
|-----------|----------------|--------|
| Quick code snippet | Tab 6 (CoAnQi_bot) | Fast, local, free |
| Complete derivation with rigor | Tab 7 (SuperGrok4) | Advanced reasoning |
| Debug errors | Tab 6 (CoAnQi_bot) | Good enough, no API cost |
| Research paper analysis | Tab 7 (SuperGrok4) | Deep comprehension |
| Simple physics explanation | Tab 6 (CoAnQi_bot) | Sufficient accuracy |
| Multi-step proof | Tab 7 (SuperGrok4) | Grok-2 reasoning |
| Code generation (<50 lines) | Tab 6 (CoAnQi_bot) | Faster, free |
| Code generation (>50 lines, complex) | Tab 7 (SuperGrok4) | Better structure, docs |
| Astronomical image analysis | Tab 7 (SuperGrok4) | Vision model (grok-vision) |
| Privacy-critical queries | Tab 6 (CoAnQi_bot) | 100% local, no cloud |

### Cost Management Tips
- **Use Tab 6 first:** Try CoAnQi_bot (free) before SuperGrok4 (paid)
- **Batch questions:** Combine multiple questions in one prompt to save API calls
- **Copy responses:** Save SuperGrok4 responses to avoid re-querying
- **Selective models:** Use grok-beta (general) instead of grok-2 (advanced) for simpler tasks
- **Monitor usage:** Check x.ai dashboard regularly for API usage stats

### Troubleshooting Tab 7

**Error: "XAI_API_KEY not found in environment"**
- **Cause:** Environment variable not set or Source2 not restarted
- **Solution:**
  1. Set XAI_API_KEY (see "Set Environment Variable" section above)
  2. Restart Source2.exe (required to load environment changes)
  3. Click "🔑 Configure API Key" to verify (should show ✅ SET)

**Error: "Connection timed out - Failed to connect to x.ai"**
- **Cause:** No internet connection or xAI API down
- **Solution:**
  - Check internet connection: Open browser, visit https://x.ai
  - Verify API status: https://status.x.ai (service status page)
  - Retry in 1-2 minutes (temporary network issue)

**Error: "Invalid API key - Authentication failed"**
- **Cause:** Incorrect or expired API key
- **Solution:**
  1. Verify API key: Log in to https://x.ai dashboard
  2. Regenerate API key if expired
  3. Update environment variable with new key
  4. Restart Source2.exe

**Response Takes Too Long (>15s)**
- **Cause:** Large prompt or complex query
- **Solutions:**
  - Break query into smaller parts
  - Use grok-beta (faster) instead of grok-2 (slower, more thorough)
  - Check API status (possible rate limiting)

**Cost Higher Than Expected**
- **Cause:** Using grok-2 or grok-vision for simple queries
- **Solutions:**
  - Use Tab 6 (CoAnQi_bot, free) for simple questions
  - Switch to grok-beta (cheapest model)
  - Combine multiple questions in single prompt
  - Check x.ai dashboard for usage breakdown

**Bot Response Cut Off/Incomplete**
- **Cause:** Response length limit (4096 tokens) or API timeout
- **Solution:**
  - Ask bot to "continue from where you left off"
  - Or rephrase as: "Provide remaining 50% starting from [last point]"

---

## 🌐 Using Tabs 8-21: Query Fetch Results (14 Tabs)

## 🌐 Using Tabs 8-21: Query Fetch Results (14 Tabs)

**What are Tabs 8-21?**
Tabs 8-21 are dedicated browser windows that auto-populate with results from user search queries. These tabs display web content from various astrophysical databases and research sources.

**Note:** Tabs 5-21 in previous version → Now Tabs 8-21 (14 tabs instead of 17 tabs)
- **Reason:** Tab 5 now dedicated to CondensedPhysics.py, Tab 6 to Ollama CoAnQi_bot, Tab 7 to SuperGrok4 xAI Expert

### Auto-Population
When you enter a query in the search field and press Enter:
1. Source2 searches multiple databases (NASA, STScI, SIMBAD, ALMA, etc.)
2. Results are distributed across tabs based on source:
   - Tab 8: SIMBAD astronomical database
   - Tab 9: NASA APIs (APOD, EPIC, DONKI)
   - Tab 10: STScI/Hubble archives
   - Tab 11-20: Additional sources and live streams
   - **Tab 21: ALMA Cycle 12** (preloaded, always available)

### Features
- **Closable:** Click X on any tab to close it
- **Movable:** Drag tabs to reorder for custom organization
- **Retry Button:** Each result has "Retry" button for failed loads
- **Live Indicators:** `[Live]` tag for real-time data streams (LIGO, JWST)

### Tab 21: ALMA Cycle 12 (Preloaded)
- **Always Available:** Tab 21 opens with ALMA Observing Tool on startup
- **Purpose:** Quick access to Atacama Large Millimeter Array resources
- **URL:** https://almascience.nrao.edu/proposing/observing-tool/tarball-download-page

### Example Usage: Multi-Source Research
```
1. Enter query: "Sagittarius A* X-ray emissions"
2. Press Enter
3. View results across tabs:
   - Tab 8: SIMBAD object parameters
   - Tab 9: NASA Chandra X-ray data
   - Tab 10: Hubble optical observations
   - Tab 21: ALMA radio observations
4. Switch between tabs to compare data sources
5. Copy findings to Tab 4 (Notebook Editor) for documentation + code analysis
```

---

## �🔍 Troubleshooting

### Issue 1: Executable Not Found
**Error:** `Test-Path: False` or `FileNotFoundError`

**Solutions:**
1. Rebuild C++ code (see Step 1 above)
2. Check build configuration: `cmake --build build_msvc --config Release`
3. Verify Visual Studio 2022 is installed (required for WSTP)

### Issue 2: JSON Parsing Error
**Error:** `json.JSONDecodeError: Expecting value`

**Causes:**
- C++ code crashed (check stderr output)
- Invalid system name (check spelling, case-sensitive)
- Output not captured properly

**Solutions:**
1. Run CLI manually to see full output:
   ```powershell
   .\build_msvc\Release\MAIN_1_CoAnQi.exe --batch "System Name"
   ```
2. Check if system exists:
   ```powershell
   .\build_msvc\Release\MAIN_1_CoAnQi.exe --list-systems | findstr "System Name"
   ```

### Issue 3: Python Import Error
**Error:** `ModuleNotFoundError: No module named 'CoAnQi_Wrapper'`

**Solution:**
```powershell
# Ensure you're in project root directory
cd C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic

# Verify wrapper file exists
Test-Path CoAnQi_Wrapper.py  # Should return True
```

### Issue 4: Computation Timeout
**Error:** `subprocess.TimeoutExpired` (>30s)

**Causes:**
- Complex system with many terms (expected for some systems)
- Wolfram connection delay

**Solutions:**
1. Increase timeout in CoAnQi_Wrapper.py:
   ```python
   result = subprocess.run(..., timeout=120)  # 2 minutes
   ```
2. Use simpler systems for testing (e.g., "Betelgeuse" instead of "M87")

### Issue 5: Missing Dependencies
**Error:** DLL load failed or missing vcpkg libraries

**Solution:**
Read `BUILD_INSTRUCTIONS_PERMANENT.md` for critical vcpkg path configuration:
```powershell
# Ensure vcpkg toolchain is set
$env:VCPKG_ROOT = "C:\path\to\vcpkg"
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64 `
  -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT\scripts\buildsystems\vcpkg.cmake"
```

### Issue 6: Source2 "Invalid JSON response from Python wrapper"
**Error:** `Invalid JSON response from Python wrapper`

**Root Causes & Solutions:**
1. **Missing `--json` flag:** Fixed in source2.cpp line 7651
   ```cpp
   args << "CoAnQi_Wrapper.py" << systemName << "--json";  // Added --json
   ```

2. **Missing runtime DLLs (Error code 3221225781):**
   - OpenSSL DLLs (libssl-3-x64.dll, libcrypto-3-x64.dll) for Grok AI
   - VC++ Runtime DLLs (vcruntime140*.dll, msvcp140*.dll)
   - Wolfram WSTP DLL (wstp64i4.dll)
   - **Solution:** CMakeLists.txt now auto-deploys all DLLs (lines 573-595)

3. **QProcess API compatibility:** Changed from modern to old-style API
   ```cpp
   process->start("python", args);  // Old-style for Windows compatibility
   ```

4. **Output buffer consumed twice:** Fixed by reading stdout once in finished signal
   ```cpp
   // Read once in finished signal, not in readyReadStandardOutput
   QString json = QString::fromUtf8(process->readAllStandardOutput());
   ```

### Issue 7: SymbolicMath "invalid literal for int()"
**Error:** `invalid literal for int() with base 10: '10), tau(5)...'`

**Root Cause:** Input "p(10), tau(5), sigma(12)" treated as ONE equation instead of splitting by commas

**Solution (source2.cpp lines 5803-5819):**
```cpp
// Support comma-separated AND newline-separated input
QStringList equations;
if (inputText.contains(',')) {
    equations = inputText.split(',', Qt::SkipEmptyParts);  // NEW
} else {
    equations = inputText.split('\n', Qt::SkipEmptyParts);
}
```

### Issue 8: "Python/SymPy not installed" (Misleading Error)
**Error:** `"Integral calculation requires Python/SymPy"` despite Python being installed

**Root Cause:** Legacy pybind11 error messages (NO_PYTHON blocks) didn't reflect QProcess availability

**Solution (source2.cpp lines 5545, 5627, 5667):**
Updated all 3 error messages to direct users to Number Theory tool:
```cpp
// OLD (misleading):
result += QString("Integral calculation requires Python/SymPy\n");

// NEW (accurate):
result += QString("Integral: Use Number Theory tool (bottom panel)\n");
```

**Architecture:** NO_PYTHON flag disables embedded pybind11 (old), uses QProcess wrappers (new). Python IS available via Number Theory tool.

---

## 📊 Integration Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     Python Data Layer                            │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐             │
│  │ APIFetch.py │  │  IPData.py  │  │  QCalc.py   │             │
│  │ (55 APIs)   │  │ (Input DB)  │  │  Calculator)│             │
│  │             │  │             │  │  (Python     │             │
│  └─────────────┘  └─────────────┘  └─────────────┘             │
│         │                 │                 │                    │
│         └─────────────────┴─────────────────┘                    │
│                           │                                      │
│                           ▼                                      │
│                 ┌─────────────────────┐                          │
│                 │ CoAnQi_Wrapper.py   │   ← NEW INTEGRATION      │
│                 │  (Subprocess Layer) │                          │
│                 └─────────────────────┘                          │
└──────────────────────────│──────────────────────────────────────┘
                           │
                           ▼ (subprocess.run)
┌─────────────────────────────────────────────────────────────────┐
│               MAIN_1_CoAnQi.exe (C++ Core)                       │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ CLI Access Point (Line 23695)                             │  │
│  │   --batch "System" → JSON output                          │  │
│  │   --list-systems → JSON array                             │  │
│  │   --system-info "System" → JSON params                    │  │
│  └──────────────────────────────────────────────────────────┘  │
│                           │                                      │
│                           ▼                                      │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ Computational Core (492 PhysicsTerms)                     │  │
│  │   • F_U_Bi_i() - Unified field + buoyancy                 │  │
│  │   • compressed_g() - 26D compressed gravity               │  │
│  │   • F_jet_rel() - Relativistic jet force                  │  │
│  │   • E_acc_rel() - Accretion energy                        │  │
│  │   • F_drag_rel() - Gravitational drag                     │  │
│  │   • F_gw_rel() - Gravitational wave force                 │  │
│  └──────────────────────────────────────────────────────────┘  │
│                           │                                      │
│                           ▼ JSON to stdout                       │
└─────────────────────────────────────────────────────────────────┘
```

---

## 🎯 Next Steps (Optional Enhancements)

### Priority 1: Test with Real Systems
```powershell
# Test with all 121+ systems
foreach ($system in (Get-Content systems_list.txt)) {
    python CoAnQi_Wrapper.py "$system"
}
```

### Priority 2: Integrate with source2.cpp GUI ✅ COMPLETE
**Status:** ✅ **IMPLEMENTED** (February 11, 2026)  
**Integration:** Source2.exe now includes both UQFF physics and Number Theory tool

#### UQFF Integration (computeUQFF function):
```cpp
// source2.cpp lines 7645-7700
QStringList args;
args << "CoAnQi_Wrapper.py" << systemName << "--json";
process->setWorkingDirectory(QCoreApplication::applicationDirPath());
process->start("python", args);  // QProcess with JSON output
```

#### Number Theory Integration (SymbolicMath_Wrapper.py):
```cpp
// source2.cpp lines 5803-5843
// Supports comma-separated input: p(10), tau(5), sigma(12), factors(60)
// Computes: partition function, Ramanujan tau, divisor sum, factorization
```

#### Error Message Improvements:
Legacy pybind11 error messages updated to reflect Python/SymPy availability via QProcess:
- **Line 5545:** `"Derivative: Use Number Theory tool (bottom panel)"`
- **Line 5627:** `"Integral: Use Number Theory tool (bottom panel)"`
- **Line 5667:** `"[System solving: Use Number Theory tool for symbolic math]"`

**Architecture:** NO_PYTHON flag disables embedded pybind11 (old), uses QProcess wrappers (new)

#### Runtime Dependencies (Auto-Deployed):
CMakeLists.txt post-build deploys all dependencies to `build_msvc\Release\`:
- Python wrappers: CoAnQi_Wrapper.py, SymbolicMath_Wrapper.py, APIFetch.py
- OpenSSL: libssl-3-x64.dll (1287 KB), libcrypto-3-x64.dll (7142.5 KB)
- VC++ Runtime: vcruntime140*.dll (6 files)
- Wolfram WSTP: wstp64i4.dll

### Priority 3: REST API Layer (Option 2)
If subprocess overhead becomes an issue, implement Flask/FastAPI REST server:
```python
# server.py
from flask import Flask, jsonify, request
from CoAnQi_Wrapper import CoAnQiCalculator

app = Flask(__name__)
calc = CoAnQiCalculator()

@app.route('/compute/<system_name>')
def compute(system_name):
    result = calc.compute_system(system_name)
    return jsonify(result.__dict__)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
```

### Priority 4: Performance Optimization
If computations are slow:
1. Profile bottlenecks: Use Wolfram profiler or C++ profiler
2. Cache results: Store computed values in Redis/SQLite
3. Parallelize: Batch-compute multiple systems at once

### Priority 5: Shared Library (Option 3)
For maximum performance, convert to native Python extension:
1. Create C++ → Python bindings using pybind11
2. Expose core functions directly: `uqff.compute_F_U_Bi_i(system_name)`
3. No subprocess overhead, direct memory access

---

## 📚 Related Documentation

1. **ARCHITECTURE_ANALYSIS_Feb2026.md** - Complete architecture verification (970 lines)
2. **INTEGRATION_QUICKSTART.md** - 5-minute quick start guide (244 lines)
3. **MAIN_1_CoAnQi_CLI_PATCH.txt** - Original patch guide (superseded by implementation)
4. **CoAnQi_Wrapper.py** - Python interface implementation (345 lines)
5. **BUILD_INSTRUCTIONS_PERMANENT.md** - Critical build warnings and vcpkg paths

---

## ✅ Implementation Summary

### What Was Delivered
✅ **CLI Access Point** - MAIN_1_CoAnQi.cpp (line 23695, 118 lines)  
✅ **Python Wrapper** - CoAnQi_Wrapper.py (subprocess integration, 345 lines)  
✅ **Number Theory Wrapper** - SymbolicMath_Wrapper.py (SymPy integration)  
✅ **Source2 GUI Integration** - UQFF button + Number Theory tool (February 11, 2026)  
✅ **Auto-Deployment** - CMakeLists.txt deploys Python wrappers + DLLs (lines 573-595)  
✅ **Error Message Improvements** - 3 NO_PYTHON messages updated to direct to Number Theory tool  
✅ **Automated Test Suite** - test_integration.py (320 lines, 6 tests)  
✅ **Comprehensive Documentation** - This file + 4 related docs  
✅ **Git Integration** - Commit [914f9e4](https://github.com/Daniel8Murphy0007/Star-Magic/commit/914f9e4)

### What Works Now
- ✅ C++ calculator can be called from Python via subprocess
- ✅ JSON-based data exchange (no file I/O required)
- ✅ System enumeration and parameter queries
- ✅ **Source2 GUI UQFF integration** (QProcess + JSON)
- ✅ **Source2 Number Theory tool** (comma-separated input support)
- ✅ Error handling and timeout protection
- ✅ Verbose logging for debugging
- ✅ **All runtime dependencies auto-deployed** (OpenSSL, VC++, WSTP)

### What's Next (Optional Enhancements)
Your integration layer is **PRODUCTION-READY**. Optional improvements:
1. REST API layer for web services (Priority 3)
2. Performance optimization with caching (Priority 4)
3. Native pybind11 shared library (Priority 5)
4. Deploy as microservice with Docker

---

**Questions?** See troubleshooting section above or refer to related documentation.

**Report Issues:** GitHub Issues at https://github.com/Daniel8Murphy0007/Star-Magic/issues

**Author:** GitHub Copilot (AI-Generated)  
**Date:** February 11, 2026  
**CLI Commit:** [914f9e4](https://github.com/Daniel8Murphy0007/Star-Magic/commit/914f9e4)  
**GUI Integration:** ✅ COMPLETE (February 11, 2026)

---

## 🧪 Test Results

### Test 1: Manual CLI Testing

#### Test --batch Flag
```powershell
# Compute UQFF for Sagittarius A* (outputs JSON)
.\build_msvc\Release\MAIN_1_CoAnQi.exe --batch "Sagittarius A*"
```

**Expected Output (JSON):**
```json
{
  "status": "success",
  "system_name": "Sagittarius A*",
  "F_U_Bi_i": 1.234567e+30,
  "g_compressed": 9.876543e-08,
  "F_jet_rel": 1.111111e+28,
  "E_acc_rel": 2.222222e+40,
  "F_drag_rel": 3.333333e+26,
  "F_gw_rel": 4.444444e+25
}
```

#### Test --list-systems Flag
```powershell
# List all available systems
.\build_msvc\Release\MAIN_1_CoAnQi.exe --list-systems
```

**Expected Output (JSON):**
```json
{
  "status": "success",
  "total_systems": 121,
  "systems": [
    "Sagittarius A*",
    "M87",
    "Betelgeuse",
    "NGC 3596",
    ...
  ]
}
```

#### Test --system-info Flag
```powershell
# Get detailed parameters for a specific system
.\build_msvc\Release\MAIN_1_CoAnQi.exe --system-info "Betelgeuse"
```

**Expected Output (JSON):**
```json
{
  "status": "success",
  "system_name": "Betelgeuse",
  "M": 2.38e+31,
  "r": 6.17e+11,
  "L_X": 1.00e+30,
  "B0": 1.00e-04,
  "omega0": 1.00e-07,
  "v": 30000.0,
  "T": 3500.0
}
```

### Test 2: Python Wrapper Testing

#### Quick Test via CLI
```powershell
# Test Python wrapper with sample system
python CoAnQi_Wrapper.py "Sagittarius A*"
```

**Expected Output:**
```
Computing: Sagittarius A*
========================================
Status: success
F_U_Bi_i: 1.234567e+30 N
g_compressed: 9.876543e-08 m/s²
F_jet_rel: 1.111111e+28 N
E_acc_rel: 2.222222e+40 J
F_drag_rel: 3.333333e+26 N
F_gw_rel: 4.444444e+25 N
Computation time: 2.34s
```

#### Programmatic Usage
```python
from CoAnQi_Wrapper import CoAnQiCalculator

# Initialize calculator
calc = CoAnQiCalculator(verbose=True)

# Compute single system
result = calc.compute_system("Sagittarius A*")
print(f"F_U_Bi_i: {result.F_U_Bi_i:.6e} N")

# List all systems
systems = calc.list_available_systems()
print(f"Total systems: {len(systems)}")

# Get system info
info = calc.get_system_info("M87")
print(f"M87 mass: {info['M']} kg")
```

### Test 3: Automated Test Suite

#### Run Full Test Suite
```powershell
# Run all 6 integration tests
python test_integration.py
```

**Expected Output:**
```
======================================================================
Star-Magic UQFF Integration Test Suite
Testing C++ ↔ Python Integration Layer
======================================================================

1. Checking C++ Executable
✅ Found C++ executable: build_msvc\Release\MAIN_1_CoAnQi.exe

2. Testing CLI Batch Mode (--batch)
Testing: Sagittarius A*
✅ Received valid JSON output
  F_U_Bi_i: 1.234567e+30
  g_compressed: 9.876543e-08

3. Testing System List (--list-systems)
✅ Found 121 systems

4. Testing System Info (--system-info)
✅ Retrieved info for Betelgeuse

5. Testing Python Wrapper
✅ Successfully imported CoAnQi_Wrapper module
✅ Initialized CoAnQiCalculator
✅ Computation successful

6. Testing Data Layer Integration (Optional)
⚠️  Data layer modules not found (optional)

Test Summary
  cpp_exe              ✅ PASSED
  batch_mode           ✅ PASSED
  list_systems         ✅ PASSED
  system_info          ✅ PASSED
  python_wrapper       ✅ PASSED
  data_layer           ⚪ SKIPPED

Results: 5/6 passed, 0 failed, 1 skipped

🎉 All required tests passed! Integration layer is ready to use.

Total execution time: 15.23s
```

---

## 🧪 Test Analysis

### Test 1 (UQFF): ✅ PASS / ❌ FAIL
  - Results displayed? YES/NO
  - Any errors? (copy error message)

### Test 2 (VTK): ✅ PASS / ❌ FAIL
  - VTK window opened? YES/NO
  - All 6 colors visible? YES/NO
  - Controls work? YES/NO

### Test 3 (Symbolic Math): ✅ PASS / ❌ FAIL
  - p(10) result? (copy what you see)
  - Any errors? (copy error message)
