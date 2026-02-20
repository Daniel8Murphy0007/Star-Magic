/**
 * @file astro_graphics.h
 * @brief Astronomical Graphics Engine for UQFF Star-Magic VR Runtime
 * 
 * Integration point for third-party astronomical graphics programs:
 * - Load external executable or library
 * - GPU tasking for astronomical visualization
 * - Data streaming from physics backend
 * - Real-time UQFF field overlay
 * 
 * Supported Graphics Backends:
 * - Custom OpenGL/Vulkan programs
 * - VTK scientific visualization
 * - Stellarium-style star catalogs
 * - Point cloud renderers
 * 
 * Data Sources:
 * - Physics Backend (IPC) for field data
 * - SIMBAD/NASA API for object catalogs
 * - Local CSV/JSON files for system parameters
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 3 - VR Runtime Scaffold
 */

#ifndef ASTRO_GRAPHICS_H
#define ASTRO_GRAPHICS_H

#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <map>

namespace VR {

// Forward declarations
class VulkanCompositor;

/**
 * @struct CatalogEntry
 * @brief Single astronomical object from catalog
 */
struct CatalogEntry {
    std::string id;                     // SIMBAD identifier
    std::string name;                   // Common name
    std::string type;                   // Object type (star, galaxy, nebula, etc.)
    
    // Position
    double ra = 0;                      // Right ascension (degrees)
    double dec = 0;                     // Declination (degrees)
    double distance_pc = 0;             // Distance (parsecs)
    
    // Physical properties
    double mass_solar = 0;              // Mass in solar masses
    double radius_solar = 0;            // Radius in solar radii
    double luminosity_solar = 0;        // Luminosity in solar luminosities
    double temperature_k = 0;           // Surface temperature (K)
    double magnetic_field_t = 0;        // Magnetic field (Tesla)
    double spin_period_s = 0;           // Rotation period (seconds)
    
    // Visual
    double apparent_magnitude = 0;
    double absolute_magnitude = 0;
    float color_bv = 0;                 // B-V color index
    
    // UQFF properties (calculated)
    double F_U = 0;                     // Local unified field
    double field_gradient = 0;          // Field gradient magnitude
};

/**
 * @struct FieldOverlayConfig
 * @brief Configuration for UQFF field overlay on astronomical scene
 */
struct FieldOverlayConfig {
    bool enabled = true;
    double sample_radius_pc = 100.0;    // Radius to sample field (parsecs)
    int grid_resolution = 64;           // Points per dimension
    
    // Field display options
    bool show_field_magnitude = true;
    bool show_field_lines = false;
    bool show_equipotentials = false;
    
    // Color mapping
    enum class ColorScheme {
        HeatMap,        // Blue → Red
        Viridis,        // Perceptually uniform
        Plasma,         // High contrast
        Custom
    } color_scheme = ColorScheme::HeatMap;
    
    // Transparency
    float opacity = 0.5f;
    bool depth_fade = true;             // Fade with distance
};

/**
 * @struct SceneConfig
 * @brief Overall scene configuration
 */
struct SceneConfig {
    // Camera
    double camera_position[3] = {0, 0, 0};  // In parsecs from galactic center
    double camera_target[3] = {0, 0, 0};
    double fov_degrees = 60.0;
    
    // Display
    double min_magnitude = -30.0;       // Brightest to display
    double max_magnitude = 20.0;        // Faintest to display
    double scale_factor = 1.0;          // Visual scale
    
    // Time
    double epoch_jd = 2460000.0;        // Julian date
    bool animate = false;
    double time_scale = 1.0;            // Real-time multiplier
    
    // Labels
    bool show_object_labels = true;
    bool show_coordinate_grid = false;
    double label_min_size = 5.0;        // Pixels
    
    // Performance
    int max_visible_objects = 100000;
    double lod_distance_pc = 10.0;      // Distance for LOD switching
};

/**
 * @class AstroGraphics
 * @brief Astronomical graphics engine wrapper
 */
class AstroGraphics {
public:
    AstroGraphics();
    ~AstroGraphics();
    
    // Initialization
    bool initialize(VulkanCompositor* compositor);
    void shutdown();
    
    // External program loading
    bool loadExternalProgram(const std::string& path);
    bool loadVTK(const std::string& data_path);
    bool loadCatalog(const std::string& catalog_file);  // CSV or JSON
    
    // Catalog management
    void addCatalogEntry(const CatalogEntry& entry);
    CatalogEntry* findEntry(const std::string& id);
    std::vector<CatalogEntry*> searchNear(double ra, double dec, double radius_deg);
    std::vector<CatalogEntry*> searchByName(const std::string& pattern);
    int getCatalogSize() const { return static_cast<int>(catalog_.size()); }
    
    // UQFF integration
    void setFieldOverlayConfig(const FieldOverlayConfig& config);
    void updateFieldData(const double* grid_data, int grid_dims[3], double grid_origin[3], double cell_size);
    void calculateFieldAtObjects();      // Batch calculate F_U for all objects
    
    // Scene control
    void setSceneConfig(const SceneConfig& config);
    void setCameraPosition(const double* position);
    void setCameraTarget(const double* target);
    void flyTo(const std::string& object_id, double duration_s);
    void orbit(const std::string& object_id, double radius_pc, double speed);
    
    // Rendering
    void render(int eye_index, const float* view_matrix, const float* proj_matrix);
    void renderOverlay();                // UI elements, labels
    
    // Selection
    std::string pickObject(int screen_x, int screen_y);  // Ray pick
    void selectObject(const std::string& id);
    void clearSelection();
    CatalogEntry* getSelectedObject() { return selected_object_; }
    
    // Data streaming
    using DataCallback = std::function<void(const std::string& object_id, const CatalogEntry& data)>;
    void setDataCallback(DataCallback callback) { data_callback_ = callback; }
    void requestObjectData(const std::string& id);  // Async fetch from API
    
    // State
    struct RenderStats {
        int visible_objects = 0;
        int triangles_rendered = 0;
        double render_time_ms = 0;
    };
    RenderStats getStats() const { return stats_; }

private:
    // Vulkan integration
    VulkanCompositor* compositor_ = nullptr;
    
    // External program (optional)
    void* external_program_ = nullptr;
    std::string external_path_;
    
    // Catalog
    std::vector<CatalogEntry> catalog_;
    std::map<std::string, size_t> catalog_index_;  // ID → index
    
    // Field data
    std::vector<double> field_grid_;
    int field_dims_[3] = {0, 0, 0};
    double field_origin_[3] = {0, 0, 0};
    double field_cell_size_ = 1.0;
    
    // Configuration
    FieldOverlayConfig field_config_;
    SceneConfig scene_config_;
    
    // Selection
    CatalogEntry* selected_object_ = nullptr;
    
    // Callbacks
    DataCallback data_callback_;
    
    // Stats
    RenderStats stats_;
    
    // Rendering resources
    struct RenderResources;
    std::unique_ptr<RenderResources> render_resources_;
    
    // Helper methods
    void buildSpatialIndex();
    void updateVisibleSet();
    void renderStars();
    void renderGalaxies();
    void renderNebulae();
    void renderFieldOverlay();
    void renderLabels();
    double calculateApparentMagnitude(const CatalogEntry& entry);
    float colorIndexToRGB(float bv, float* rgb);
};

/**
 * @namespace StarCatalog
 * @brief Pre-defined star catalog constants
 */
namespace StarCatalog {
    // Catalog file formats
    constexpr const char* FORMAT_CSV = "csv";
    constexpr const char* FORMAT_JSON = "json";
    constexpr const char* FORMAT_FITS = "fits";
    
    // Standard catalog names
    constexpr const char* HIPPARCOS = "hipparcos";
    constexpr const char* GAIA_DR3 = "gaia_dr3";
    constexpr const char* TYCHO2 = "tycho2";
    constexpr const char* SIMBAD = "simbad";
    
    // Object types
    constexpr const char* TYPE_STAR = "Star";
    constexpr const char* TYPE_GALAXY = "Galaxy";
    constexpr const char* TYPE_NEBULA = "Nebula";
    constexpr const char* TYPE_CLUSTER = "Cluster";
    constexpr const char* TYPE_MAGNETAR = "Magnetar";
    constexpr const char* TYPE_PULSAR = "Pulsar";
    constexpr const char* TYPE_BH = "BlackHole";
    constexpr const char* TYPE_SMBH = "SMBH";
}

} // namespace VR

#endif // ASTRO_GRAPHICS_H
