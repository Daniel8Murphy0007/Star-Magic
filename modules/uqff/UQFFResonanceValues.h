// UQFFResonanceValues.h
// Hardcoded resonance solution table for UQFF validation.
// Values derived from Grok analysis Oct 09, 2025 (grok_share_5fa36e4e035.txt).
// Watermark: Copyright - Daniel T. Murphy, analyzed Oct 09, 2025.

#ifndef UQFF_RESONANCE_VALUES_H
#define UQFF_RESONANCE_VALUES_H

#include <map>
#include <string>

namespace UQFF {

// Expected g_resonance values (m/s²) per system.
// Used for validation: |computeG(t) - resonance| / resonance < tolerance
inline std::map<std::string, double> getResonanceTable() {
    return {
        // Doc 34a systems
        {"UniverseDiameter",         7.579e53},
        {"HydrogenAtom",             1.975e-7},
        {"LagoonNebula",             1.667e29},
        {"SpiralsSupernovae",        4.353e35},
        {"NGC6302",                  4.113e20},
        {"OrionNebula",              3.458e26},
        {"UniverseGuide",            3.958e14},
        // Doc 34b systems
        {"GalaxiesGalore",           4.353e35},
        {"StellarForge",             1.001e27},
        {"NewStars",                 1.001e27},
        {"SombreroGalaxy",           1.000e36},
        {"Saturn",                   7.401e3},
        {"CrabNebula",               8.343e24},
        // Doc 39b (magnetar, mode-specific)
        {"MagnetarSGR1745_compressed", 1.782e39},
        {"MagnetarSGR1745_frequency",  1.773e-9},
    };
}

} // namespace UQFF

#endif // UQFF_RESONANCE_VALUES_H
