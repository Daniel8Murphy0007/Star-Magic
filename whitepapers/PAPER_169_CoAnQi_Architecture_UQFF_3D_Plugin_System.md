# PAPER_169: CoAnQi Architecture — Multi-Tier UQFF+3D+Plugin System
## Whitepaper §2.4-A | Thread 381a8fe7 | Session 48

### Abstract
CoAnQi (CoAnQi Labs) is a multi-tier simulation framework that unifies the
Unified Quantum Field Framework (UQFF) physics engine with a full 3D graphics
pipeline and a cross-platform runtime plugin system. This paper documents the
canonical six-tier architecture extracted from the CoAnQi codebase shared in
Grok thread 381a8fe78e1a4ecbaf32a88aa386df30.

---

### 1. System Tiers

| Tier | Component | Role |
|------|-----------|------|
| 1 | source2.cpp (Qt6 GUI, 15,753 L) | Principal user interface — all workflows start here |
| 2 | MAIN_1_CoAnQi.cpp (107,019 L) | C++ physics calculator (446 modules) |
| 3 | QCalc.py, CondensedPhysics.py/2.py | Python parallel calculators |
| 4 | uqff_server.js (imports index.js) | REST API (port 3141) |
| 5 | source2(HEAD PROGRAM).cpp | VR/VM GPU backend |
| 6 | physics_backend.cpp | CPU-bound headless server |

Additionally the codebase exposes:
- **IPC**: `uqff_ipc.h`, `python_bridge.h`, `physics_service.h`
- **Storage**: `bodies_*.csv`, `uqff_results.json`, `CondensedPhysics_OutputData.py`

---

### 2. CoAnQi File Structure (from thread)

```
CelestialBody.h/cpp   — 12-field body descriptor + Ug1/Ug2/Ug3/Um helpers
MUGE.h/cpp            — ResonanceParams, MUGESystem, compressed/resonance functions
main.cpp              — global constants, Ug4, Ubi, A_μν, FU, quasar jet sim
FluidSolver.h/cpp     — 2D Navier-Stokes (N=32, dt=0.1)
UnitTests.cpp         — 26 validated unit tests
ModelLoader.h/cpp     — OBJ import/export
Texture.h/cpp         — stb_image OpenGL loader
Shader.h/cpp          — GLSL compile/link
Camera.h/cpp          — lookAt, multi-viewport render
Animation.h/cpp       — Bone/BoneInfo, SLERP skeletal animation
Landscape.h/cpp       — procedural terrain generation
Modeling.h/cpp        — extrudeMesh, booleanUnion stubs
LaTeXRenderer.h/cpp   — MicroTeX integration
PluginModule.h/cpp    — SIMPlugin cross-platform dynamic loader
CoAnQiNode.py         — Python mirror (OpenGL/Vulkan/PyQt5/vtk/OpenCV)
```

---

### 3. Data Flow

```
User query (source2.cpp)
        ↓
APIFetch.py → bodies_*.csv
        ↓
5 parallel calculators (MAIN_1 + QCalc + CP + CP2 + uqff_server)
        ↓
OPData.py → uqff_results.json
        ↓
3D Simulation loop → renderMultiViewports → glfwSwapBuffers
```

---

### 4. Physics ↔ Visual Loop Integration

The function `populate_simulation_entities()` maps each **MUGESystem** instance
directly to a **SimulationEntity** containing:
- `position[3]`: initialized from system distances
- `velocity[3]`: from system vexp
- `model`: a 3DObject (MeshData) loaded from OBJ or procedurally generated
- `archive`: media files stored at time of simulation

The function `simulate_fluids_for_muge()` runs a **FluidSolver** step coupled
to `compute_resonance_MUGE()` for each system, so the Navier-Stokes solver
receives the UQFF gravity as an external body force.

---

### 5. Plugin Architecture — SIMPlugin

```cpp
struct SIMPlugin {
    void* handle;                     // dlopen / LoadLibrary
    void (*playAPI)(const char*);     // exported function pointer
    std::string name, path;
};
```

Load/unload uses `dlopen`+`dlsym` (POSIX) or `LoadLibrary`+`GetProcAddress`
(Windows) wrapped behind a uniform interface. Plugins inject into the
simulation loop via the `playAPI` callback.

---

### 6. Build Configuration

```cmake
Generator: Visual Studio 17 2022 / x64
Standard: C++20
Flags: /O2 /GL /DUSE_COSMIC_QUANTUM_EGG /DUSE_EMBEDDED_WOLFRAM
Post-build: UPX 5.0.2 compression (1.43 MB final)
Dependencies: Qt6, OpenGL, GLFW, GLM, stb_image, MicroTeX, Wolfram WSTP
```

---

### 7. References
- Thread 381a8fe78e1a4ecbaf32a88aa386df30 (grok_share_381a8f.txt)
- ARCHITECTURE_FLOW_DIAGRAM.md v4.4.0 (this repository)
- BUILD_INSTRUCTIONS_PERMANENT.md (this repository)
- copilot-instructions.md §Big Picture Architecture
