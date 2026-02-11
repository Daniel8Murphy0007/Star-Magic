# source2.cpp - Comprehensive Analysis & Integration Plan
**Analysis Date:** February 11, 2026  
**File:** source2.cpp (7,642 lines)  
**Build Status:** ✅ **COMPILES SUCCESSFULLY** (Release configuration)  
**Integration Status:** 🟡 **STANDALONE** - Not connected to MAIN_1_CoAnQi.cpp physics engine

---

## 🎯 Executive Summary

**source2.cpp** is the **PRIMARY USER INTERFACE** (HEAD PROGRAM) for the Star-Magic UQFF platform. It's a **Qt6-based scientific search and analysis application** designed to query 55+ astronomical data sources (NASA, STScI, JWST, ALMA, etc.) and display results across 21 simultaneous browser windows.

###Current Reality:
- ✅ **Builds cleanly** with Qt6, VTK, CURL, SQLite, OpenCV (NO_AWS flag bypasses AWS SDK errors)
- ✅ **Full GUI implementation** - 21-tab browser, scientific calculators, voice/video input stubs
- ❌ **NOT INTEGRATED** with MAIN_1_CoAnQi.cpp (the 446-module physics calculator)
- ⚠️ **Missing:** Connection between GUI queries and UQFF physics computations

---

## 📊 Architecture Breakdown

### **1. Main Classes**

| Class | Purpose | Lines | Status |
|-------|---------|-------|--------|
| **MainWindow** | Primary Qt6 GUI window | 7300-7642 | ✅ Functional |
| **BrowserWindow** | Individual tab browser (21x) | 6034-6085 | ✅ Functional |
| **EnhancedMainWindow** | Plugin system extension | 2960-3100 | ✅ Functional |
| **ScientificCalculatorDialog** | Qalculate-based solver | 5287-5500 | ✅ Functional |
| **RamanujanCalculatorDialog** | Ramanujan equations solver | 5739-5900 | ✅ Functional |
| **CalculusButtonField** | Derivative/integral toolbar | 6500-6700 | ✅ Functional |
| **PluginRepositoryManager** | Cloud plugin loader | 3070-3200 | ✅ Functional |
| **PluginManager** | Dynamic plugin system | 2500-2700 | ✅ Functional |
| **TestFramework** | Unit testing infrastructure | 2700-2900 | ✅ Functional |
| **DistributedComputing** | Worker pool system | 2200-2300 | ✅ Functional |
| **MLIntegration** | Machine learning loader | 2300-2400 | ✅ Functional |

### **2. Core Functionality**

#### **A. Multi-Source Scientific Search**
```cpp
void PerformSearch(const std::string &query, 
                   std::vector<std::string> &focus, 
                   bool online, 
                   const std::string &oauth_token)
```
- **Purpose:** Orchestrates queries across multiple astronomical data sources
- **Supported APIs:** (55+ endpoints)
  - NASA APOD, EPIC, DONKI (Space Weather)
  - MAST (Hubble, JWST, Chandra archives)
  - JPL Horizons (Ephemeris data)
  - JPL Periodic Orbits (Halo/DRO trajectories)
  - JPL JD-Cal (Julian Date conversion)
  - Event Horizon Telescope (WebSocket stream)
  - LIGO (Gravitational wave alerts)
  - ALMA Cycle 12 (Radio telescope)
- **Results Distribution:** 21 browser windows (tabs #0-20)
- **Query Routing:** Keyword detection (e.g., "ephemeris" → JPL Horizons)

#### **B. Input Modalities**
1. **Text Input** (QLineEdit) - Primary search field (max 6,000 chars)
2. **Voice Input** (PocketSphinx stub) - Microphone button (🎤)
   - **Status:** ⚠️ Function stub exists (ProcessVoiceInput), no implementation
3. **Video Gesture** (OpenCV stub) - Camera button (📹)
   - **Status:** ⚠️ Function stub exists (ProcessVideoInput), no implementation
4. **Drag-and-Drop** - File upload via QMimeData
   - **Status:** ✅ Implemented in QTextEdit

#### **C. Computational Tools**
1. **Scientific Calculator** (Qalculate library)
   - **Status:** ✅ Full symbolic math solver (algebra, calculus, chemistry)
2. **Ramanujan Calculator** (Custom equations)
   - **Status:** ✅ Specialized solver for Ramanujan's work
3. **Calculus Toolbar** (d/dx, ∫)
   - **Status:** ✅ Derivative/integral buttons with Qalculate integration

#### **D. Visualization**
- **VTK Scatter Plots** (RenderScatterPlot function)
  - **Status:** ⚠️ Stub function exists, placeholder label in GUI
  - **Purpose:** Visualize orbital data from JPL Horizons
- **Left Sidebar Dock** (Qt::LeftDockWidgetArea)
  - **Status:** ✅ Dock widget created, shows "Dataset Graph Placeholder"

#### **E. Data Persistence**
1. **SQLite Cache** (coanqi_cache.db)
   - **Status:** ✅ Functional database with `cache` table
   - **Schema:** `(url TEXT, title TEXT, summary TEXT, isLive INTEGER)`
   - **Purpose:** Offline search (OfflineSearch function)
2. **AWS S3 Sync** (cloud backup)
   - **Status:** ❌ Disabled via NO_AWS flag (linker errors)
3. **AWS Cognito** (authentication)
   - **Status:** ❌ Disabled via NO_AWS flag

### **3. Dependency Status**

| Dependency | Status | Purpose | Compile Flag |
|------------|--------|---------|--------------|
| **Qt6** | ✅ ACTIVE | GUI framework | (required) |
| **Qt6::WebEngineWidgets** | ✅ ACTIVE | Browser tabs | (required) |
| **VTK** | ✅ ACTIVE | 3D visualizations | NO_VTK |
| **CURL** | ✅ ACTIVE | HTTP/API requests | NO_CURL |
| **SQLite3** | ✅ ACTIVE | Local caching | NO_SQLITE |
| **nlohmann/json** | ✅ ACTIVE | JSON parsing | (required) |
| **OpenCV** | ⚠️ OPTIONAL | Video gesture input | NO_OPENCV |
| **PocketSphinx** | ⚠️ OPTIONAL | Voice recognition | NO_POCKETSPHINX |
| **Qalculate** | ⚠️ OPTIONAL | Calculator backend | NO_QALCULATE |
| **pybind11** | ⚠️ OPTIONAL | Python integration | NO_PYTHON |
| **AWS SDK** | ❌ DISABLED | Cloud sync | NO_AWS |

---

## 🔗 Integration Gap: source2.cpp ↔ MAIN_1_CoAnQi.cpp

### **Current Problem**
source2.cpp and MAIN_1_CoAnQi.cpp are **completely separate executables** with NO communication:

```
┌─────────────────────────────────────┐    ❌ NO CONNECTION   ┌──────────────────────────────────┐
│       source2.cpp (GUI)             │                       │  MAIN_1_CoAnQi.cpp (Physics)     │
│                                     │                       │                                  │
│  - User searches "Sagittarius A*"  │  ──X──────────────X── │  - 492 PhysicsTerms             │
│  - Queries NASA/MAST/JPL APIs      │                       │  - SOURCE1-116 modules          │
│  - Displays in 21 browser tabs     │                       │  - F_U_Bi_i calculations        │
│  - NO PHYSICS COMPUTATIONS         │                       │  - 6,643 registered terms       │
└─────────────────────────────────────┘                       └──────────────────────────────────┘
```

### **What's Missing**

#### **1. API Fetch → Physics Computation Pipeline**
When user searches for "Sagittarius A*" in source2.cpp:
- ❌ **Does NOT call** MAIN_1_CoAnQi.cpp to compute UQFF predictions
- ❌ **Does NOT display** F_U_Bi_i, Ug1-4, Ubi forces
- ❌ **Does NOT visualize** 26D compressed gravity fields
- ❌ **Does NOT compare** NASA data vs UQFF theoretical values

#### **2. No CoAnQi_Wrapper.py Integration**
The Python wrapper we built (commit f27247b) exists but is **not called** from source2.cpp:
- ✅ **EXISTS:** CoAnQi_Wrapper.py (Python subprocess wrapper)
- ✅ **EXISTS:** CLI flags in MAIN_1_CoAnQi.cpp (--batch, --list-systems, --system-info)
- ❌ **MISSING:** QProcess calls from source2.cpp to execute wrapper
- ❌ **MISSING:** JSON parsing of wrapper output in Qt GUI
- ❌ **MISSING:** GUI widgets to display physics results

#### **3. Visualization Gap**
- ⚠️ VTK scatter plots **stubbed** (line 7417: "Dataset Graph Placeholder")
- ❌ NO visualization of UQFF force fields
- ❌ NO comparison charts (NASA data vs UQFF predictions)
- ❌ NO 26D hypergraph rendering (SOURCE200 Cosmic Quantum Egg)

#### **4. Input Data Gap**
- ✅ source2.cpp **CAN fetch** astronomical parameters (M, r, L_X, B0, etc.) from APIs
- ❌ source2.cpp **DOES NOT format** this data for MAIN_1_CoAnQi.cpp input
- ❌ NO mapping from API JSON → SystemParams struct
- ❌ NO automatic CSV generation (bodies_YYYYMMDD_HHMMSS.csv format)

---

## 🛠️ Implementation Plan: Full GUI-Physics Integration

### **Phase 1: Basic Connection (2-4 hours)**
**Goal:** Enable source2.cpp to call MAIN_1_CoAnQi.cpp and display results

#### **Step 1.1: Add QProcess Wrapper in source2.cpp**
```cpp
// In MainWindow class (add new method)
void MainWindow::computeUQFF(const QString& systemName) {
    // Create QProcess to run CoAnQi_Wrapper.py
    QProcess* process = new QProcess(this);
    
    // Setup process
    process->setProgram("python");
    process->setArguments(QStringList() << "CoAnQi_Wrapper.py" << systemName);
    process->setWorkingDirectory(QCoreApplication::applicationDirPath());
    
    // Connect signals
    connect(process, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
            [this, process](int exitCode, QProcess::ExitStatus exitStatus) {
        if (exitCode == 0) {
            // Parse JSON output from wrapper
            QString jsonOutput = process->readAllStandardOutput();
            parseAndDisplayUQFFResults(jsonOutput);
        } else {
            QMessageBox::warning(this, "UQFF Error", 
                                "Failed to compute physics: " + process->readAllStandardError());
        }
        process->deleteLater();
    });
    
    // Start computation
    process->start();
}
```

#### **Step 1.2: Add JSON Parsing**
```cpp
void MainWindow::parseAndDisplayUQFFResults(const QString& jsonStr) {
    try {
        json data = json::parse(jsonStr.toStdString());
        
        // Extract values
        double FU = data["F_U_Bi_i"];
        double g_compressed = data["g_compressed"];
        std::string system = data["system_name"];
        
        // Display in new dialog window
        QString message = QString("UQFF Results for %1:\n\n"
                                 "F_U_Bi_i: %2\n"
                                 "g_compressed: %3\n"
                                 "Ug1: %4\n"
                                 "Ug2: %5\n"
                                 "Ug3: %6\n"
                                 "Ug4: %7\n"
                                 "Ubi: %8")
                         .arg(QString::fromStdString(system))
                         .arg(FU, 0, 'e', 6)
                         .arg(g_compressed, 0, 'e', 6)
                         .arg(data["Ug1"].get<double>(), 0, 'e', 6)
                         .arg(data["Ug2"].get<double>(), 0, 'e', 6)
                         .arg(data["Ug3"].get<double>(), 0, 'e', 6)
                         .arg(data["Ug4"].get<double>(), 0, 'e', 6)
                         .arg(data["Ubi"].get<double>(), 0, 'e', 6);
        
        QMessageBox::information(this, "UQFF Physics Results", message);
        
    } catch (const std::exception& e) {
        QMessageBox::critical(this, "JSON Parse Error", e.what());
    }
}
```

#### **Step 1.3: Wire Up Button in Search Handler**
```cpp
// In MainWindow constructor, modify search handler:
connect(queryField, &QLineEdit::returnPressed, [=]() {
    std::string query = queryField->text().toStdString();
    
    // ... existing API search code ...
    PerformSearch(query, focusList, online, oauth_token);
    
    // NEW: Trigger UQFF computation
    QPushButton* uqffBtn = new QPushButton("Compute UQFF Physics");
    connect(uqffBtn, &QPushButton::clicked, [this, query]() {
        computeUQFF(QString::fromStdString(query));
    });
    // Add button to top bar or results area
});
```

---

### **Phase 2: Advanced Visualization (4-8 hours)**
**Goal:** Display UQFF force fields in VTK sidebar

#### **Step 2.1: Create UQFF Results Widget**
```cpp
class UQFFResultsWidget : public QWidget {
    Q_OBJECT
private:
    QTextEdit* rawDataDisplay;
    QPushButton* exportBtn;
    QPushButton* visualizeBtn;
    
public:
    UQFFResultsWidget(QWidget* parent = nullptr) : QWidget(parent) {
        QVBoxLayout* layout = new QVBoxLayout(this);
        
        // Raw data display
        rawDataDisplay = new QTextEdit();
        rawDataDisplay->setReadOnly(true);
        layout->addWidget(new QLabel("UQFF Computation Results:"));
        layout->addWidget(rawDataDisplay);
        
        // Action buttons
        QHBoxLayout* btnLayout = new QHBoxLayout();
        exportBtn = new QPushButton("Export to CSV");
        visualizeBtn = new QPushButton("Visualize Forces");
        btnLayout->addWidget(exportBtn);
        btnLayout->addWidget(visualizeBtn);
        layout->addLayout(btnLayout);
        
        setLayout(layout);
    }
    
    void setResults(const json& data) {
        QString text;
        text += QString("System: %1\n").arg(QString::fromStdString(data["system_name"]));
        text += QString("F_U_Bi_i: %1\n").arg(data["F_U_Bi_i"].get<double>(), 0, 'e', 6);
        text += QString("g_compressed: %1\n").arg(data["g_compressed"].get<double>(), 0, 'e', 6);
        // ... add all fields ...
        rawDataDisplay->setText(text);
    }
};
```

#### **Step 2.2: Replace VTK Placeholder**
```cpp
// In MainWindow constructor, replace line 7417:
// OLD: visLayout->addWidget(new QLabel("Dataset Graph Placeholder"));
// NEW:
UQFFResultsWidget* uqffWidget = new UQFFResultsWidget();
visLayout->addWidget(uqffWidget);

// Store as member variable for updates
this->uqffResultsWidget = uqffWidget;
```

#### **Step 2.3: Add 3D Force Visualization**
```cpp
void MainWindow::visualizeUQFFForces(const json& data) {
    #ifndef NO_VTK
    // Create VTK scatter plot
    std::vector<double> radii, forces;
    
    // Generate sample data (r = 1AU to 100AU)
    for (double r = AU_TO_METERS; r < 100*AU_TO_METERS; r *= 1.1) {
        radii.push_back(r / AU_TO_METERS);  // Convert to AU
        
        // Call MAIN_1_CoAnQi physics (via wrapper)
        // ... compute FU at distance r ...
        forces.push_back(/* computed FU value */);
    }
    
    // Render in VTK widget
    RenderScatterPlot(sidebar, radii, forces);
    #endif
}
```

---

### **Phase 3: API Data → UQFF Parameters (6-12 hours)**
**Goal:** Automatically convert NASA/MAST data to SystemParams

#### **Step 3.1: Create Parameter Mapper**
```cpp
class APItoUQFFMapper {
public:
    static json mapSIMBADtoSystemParams(const json& simbadData) {
        json params;
        
        // Extract from SIMBAD JSON
        params["M"] = extractMass(simbadData);           // Solar masses
        params["r"] = extractDistance(simbadData);       // Parsecs
        params["L_X"] = extractXrayLuminosity(simbadData); // erg/s
        params["B0"] = extractMagneticField(simbadData); // Gauss
        params["omega_0"] = extractRotationPeriod(simbadData); // rad/s
        params["v"] = extractVelocity(simbadData);       // km/s
        params["T"] = extractTemperature(simbadData);    // Kelvin
        
        return params;
    }
    
    static json mapNASAtoSystemParams(const json& nasaData) {
        // Similar for NASA API responses
        // ...
    }
    
    static json mapMASTtoSystemParams(const json& mastData) {
        // Similar for MAST archive
        // ...
    }
};
```

#### **Step 3.2: Integrate in Search Pipeline**
```cpp
// In PerformSearch function (line 7093), add after API calls:
void PerformSearch(const std::string &query, ...) {
    // ... existing NASA/MAST queries ...
    
    // NEW: Auto-compute UQFF for first result
    if (!nasaResults.empty()) {
        json nasaData = /* parse first NASA result */;
        json systemParams = APItoUQFFMapper::mapNASAtoSystemParams(nasaData);
        
        // Write to bodies.csv
        writeSystemParamsToCSV(systemParams, "bodies_auto.csv");
        
        // Trigger UQFF computation
        computeUQFF(QString::fromStdString(query));
    }
}
```

---

### **Phase 4: Real-Time Streaming (8-16 hours)**
**Goal:** Live UQFF updates for LIGO/EHT WebSocket streams

#### **Step 4.1: WebSocket Integration**
```cpp
class LiveDataProcessor : public QObject {
    Q_OBJECT
private:
    QWebSocket* ligoSocket;
    QWebSocket* ehtSocket;
    
public:
    LiveDataProcessor(QObject* parent = nullptr) : QObject(parent) {
        ligoSocket = new QWebSocket();
        ehtSocket = new QWebSocket();
        
        connect(ligoSocket, &QWebSocket::textMessageReceived,
                this, &LiveDataProcessor::processLIGOAlert);
        connect(ehtSocket, &QWebSocket::textMessageReceived,
                this, &LiveDataProcessor::processEHTData);
        
        // Connect to live streams
        ligoSocket->open(QUrl("wss://ligo.org/alerts"));
        ehtSocket->open(QUrl("wss://eventhorizontelescope.org/data"));
    }
    
signals:
    void newLIGOEvent(const json& gwData);
    void newEHTData(const json& bhData);
    
private slots:
    void processLIGOAlert(const QString& message) {
        json gwData = json::parse(message.toStdString());
        
        // Extract gravitational wave parameters
        double chirp_mass = gwData["chirp_mass"];
        double distance_Mpc = gwData["distance"];
        
        // Convert to UQFF parameters
        json params = convertGWtoUQFF(gwData);
        
        // Trigger real-time computation
        emit newLIGOEvent(params);
    }
};
```

---

## 🎯 Priority Tasks (Immediate Implementation)

### **Critical Path (Hours 1-4):**
1. ✅ **Add computeUQFF() method** to MainWindow (30 min)
2. ✅ **Add parseAndDisplayUQFFResults()** JSON parser (30 min)
3. ✅ **Create "Compute UQFF" button** in search results (15 min)
4. ✅ **Test basic integration** with Sagittarius A* (15 min)
5. ✅ **Create UQFFResultsWidget** for sidebar (1 hour)
6. ✅ **Replace VTK placeholder** with actual widget (15 min)

**Expected Outcome:** User searches "Sagittarius A*" → NASA APIs fetch data → Click "Compute UQFF" → Python wrapper calls MAIN_1_CoAnQi → JSON results displayed in Qt dialog

### **High Priority (Hours 5-12):**
1. ⚠️ **APItoUQFFMapper** class for automatic parameter extraction
2. ⚠️ **VTK force visualization** (scatter plots of Ug1-4 vs radius)
3. ⚠️ **Auto-generate bodies.csv** from API responses
4. ⚠️ **Add "Compare NASA vs UQFF"** toggle button

### **Medium Priority (Hours 13-24):**
1. 🔵 **LiveDataProcessor** for LIGO/EHT WebSocket streams
2. 🔵 **Batch computation** mode (process all 21 tab results)
3. 🔵 **Export results** to JSON/CSV
4. 🔵 **Historical comparison** (store results in SQLite cache)

### **Low Priority (Future Work):**
1. 🟢 **Voice input** actual implementation (PocketSphinx)
2. 🟢 **Video gesture** recognition (OpenCV)
3. 🟢 **AWS S3 sync** fix (resolve iostream LNK2001 errors)
4. 🟢 **Plugin system** for custom UQFF modules

---

## 📝 Code Modifications Checklist

### **Files to Modify:**

#### **source2.cpp** (PRIMARY)
- [ ] Add `#include "CoAnQi_Wrapper.py"` reference
- [ ] Add `computeUQFF(QString)` method to MainWindow class
- [ ] Add `parseAndDisplayUQFFResults(QString)` method
- [ ] Add `UQFFResultsWidget` class definition
- [ ] Modify search handler to include UQFF button
- [ ] Replace VTK placeholder (line 7417)
- [ ] Add APItoUQFFMapper class
- [ ] Add LiveDataProcessor class (optional)

#### **source2_mainwindow.h** (MODIFY)
- [ ] Add forward declaration: `class UQFFResultsWidget;`
- [ ] Add member variable: `UQFFResultsWidget* uqffResultsWidget;`
- [ ] Add method declarations:
  ```cpp
  void computeUQFF(const QString& systemName);
  void parseAndDisplayUQFFResults(const QString& jsonStr);
  void visualizeUQFFForces(const json& data);
  ```

#### **CMakeLists.txt** (NO CHANGES NEEDED)
- ✅ Source2 already links Qt6::Core, Qt6::Network (QProcess available)
- ✅ nlohmann/json already included
- ✅ NO_AWS flag set (bypasses linker errors)

#### **CoAnQi_Wrapper.py** (NO CHANGES NEEDED)
- ✅ Already validates JSON output (commit f27247b)
- ✅ Already handles system not found errors
- ✅ Already provides correct JSON keys (system_name, F_U_Bi_i, etc.)

---

## 🐛 Known Bugs & TODOs

### **From source2.cpp Code Comments:**
1. **Line 6299:** `// TODO: Parse JSON response to extract actual access_token`
   - GetOAuthToken() returns placeholder token
   - **Impact:** Authenticated APIs (MAST, some NASA) may fail
   
2. **Line 6426:** `// TODO: Replace with actual audio capture and processing`
   - ProcessVoiceInput() returns empty string
   - **Impact:** Voice button non-functional
   
3. **Line 6469:** `// TODO: Replace with actual gesture recognition using OpenCV`
   - ProcessVideoInput() returns empty string
   - **Impact:** Video button non-functional
   
4. **Line 6508:** `// TODO: Add actual implementation:`
   - RenderScatterPlot() is empty stub
   - **Impact:** VTK visualizations not working
   
5. **Line 7417:** `// TODO: Add actual VTK plots`
   - Visualization sidebar shows placeholder
   - **Impact:** No data visualization
   
6. **Line 7478:** `// TODO: add actual connectivity check via ping or curl`
   - Online mode hardcoded to `true`
   - **Impact:** No offline mode detection

### **CMakeLists.txt Configuration:**
7. **Line 513:** `set(NO_AWS ON)` - AWS SDK disabled
   - AWS S3 sync non-functional
   - AWS Cognito authentication non-functional
   - **Fix:** Resolve iostream LNK2001 errors in AWS SDK

---

## 🚀 Deployment Notes

### **Current Build Configuration:**
```powershell
# Configure (already done)
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64

# Build Source2 (GUI)
cmake --build build_msvc --config Release --target Source2

# Run GUI
.\build_msvc\Release\Source2.exe
```

### **Runtime Dependencies:**
- ✅ Qt6 DLLs (Qt6Core, Qt6Widgets, Qt6Network, Qt6WebEngineWidgets)
- ✅ VTK libraries (libvtk*.dll)
- ✅ CURL libraries (libcurl.dll)
- ✅ SQLite3 (sqlite3.dll)
- ✅ OpenSSL (libssl-3-x64.dll, libcrypto-3-x64.dll) - for HTTPS APIs
- ⚠️ Python 3.10+ (for CoAnQi_Wrapper.py)
- ⚠️ MAIN_1_CoAnQi.exe (must be in same directory or PATH)

### **Packaging Requirements:**
```
Release/
├── Source2.exe
├── MAIN_1_CoAnQi.exe        # Physics engine
├── CoAnQi_Wrapper.py         # Python integration layer
├── bodies_*.csv              # System parameter databases
├── Qt6Core.dll
├── Qt6Widgets.dll
├── Qt6Network.dll
├── Qt6WebEngineWidgets.dll
├── libssl-3-x64.dll          # OpenSSL (for Qt HTTPS)
├── libcrypto-3-x64.dll
├── libcurl.dll
├── sqlite3.dll
├── vtkCommonCore-9.2.dll     # VTK libraries (50+ DLLs)
└── ... (other VTK dependencies)
```

---

## 📚 Additional Documentation

### **Related Files:**
- [MAIN_1_CoAnQi_integration_status.json](MAIN_1_CoAnQi_integration_status.json) - Physics term inventory
- [CoAnQi_Wrapper.py](CoAnQi_Wrapper.py) - Python subprocess wrapper
- [test_integration.py](test_integration.py) - Integration test suite
- [INTEGRATION_QUICKSTART.md](INTEGRATION_QUICKSTART.md) - Quick setup guide
- [BUILD_INSTRUCTIONS_PERMANENT.md](BUILD_INSTRUCTIONS_PERMANENT.md) - Build system guide

### **Copilot Instructions:**
See [.github/copilot-instructions.md](.github/copilot-instructions.md) for full architecture documentation:
- MAIN_1_CoAnQi.cpp: 106,892 lines, 492 PhysicsTerms, SOURCE1-116
- source2.cpp: 7,642 lines, Qt6 GUI, 21-window browser, 55+ APIs

---

**Analysis Complete** ✅  
**Next Step:** Implement Phase 1 (Basic Connection) to enable GUI-Physics integration

