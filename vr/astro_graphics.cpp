/**
 * @file astro_graphics.cpp
 * @brief Astronomical Graphics Engine Implementation for UQFF Star-Magic VR Runtime
 * 
 * Integration point for astronomical visualization
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 4 - Astro Graphics IPC Integration
 */

#include "astro_graphics.h"
#include "vulkan_compositor.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstring>  // Phase 4: for memset, strncpy

// Define M_PI for Windows compatibility
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <regex>

// JSON parsing (using nlohmann/json if available)
#ifdef HAVE_NLOHMANN_JSON
#include <nlohmann/json.hpp>
#endif

namespace VR {

// ============================================================================
// RenderResources internal structure
// ============================================================================

struct AstroGraphics::RenderResources {
    // Placeholder for VTK or custom rendering resources
    bool initialized = false;
};

// ============================================================================
// Constructor / Destructor
// ============================================================================

AstroGraphics::AstroGraphics() 
    : render_resources_(std::make_unique<RenderResources>()) {
}

AstroGraphics::~AstroGraphics() {
    shutdown();
}

// ============================================================================
// Initialization
// ============================================================================

bool AstroGraphics::initialize(VulkanCompositor* compositor) {
    compositor_ = compositor;
    
    if (!compositor_) {
        std::cerr << "AstroGraphics: No compositor provided" << std::endl;
        return false;
    }
    
    render_resources_->initialized = true;
    std::cout << "AstroGraphics: Initialized" << std::endl;
    
    return true;
}

void AstroGraphics::shutdown() {
    catalog_.clear();
    catalog_index_.clear();
    field_grid_.clear();
    
    if (external_program_) {
        // Cleanup external program handle
        external_program_ = nullptr;
    }
    
    render_resources_->initialized = false;
    compositor_ = nullptr;
}

// ============================================================================
// External Program Loading
// ============================================================================

bool AstroGraphics::loadExternalProgram(const std::string& path) {
    external_path_ = path;
    
    // Check if file exists
    std::ifstream file(path);
    if (!file.good()) {
        std::cerr << "AstroGraphics: External program not found: " << path << std::endl;
        return false;
    }
    
    // TODO: Load and validate external program
    // Could be a DLL/SO, executable, or scripted pipeline
    
    std::cout << "AstroGraphics: External program stub loaded: " << path << std::endl;
    return true;
}

bool AstroGraphics::loadVTK(const std::string& data_path) {
#ifdef HAVE_VTK
    // Load VTK data files (e.g., .vtk, .vtp, .vti)
    std::cout << "AstroGraphics: Loading VTK data from " << data_path << std::endl;
    return true;
#else
    std::cout << "AstroGraphics: VTK not available" << std::endl;
    (void)data_path;
    return false;
#endif
}

bool AstroGraphics::loadCatalog(const std::string& catalog_file) {
    std::ifstream file(catalog_file);
    if (!file.good()) {
        std::cerr << "AstroGraphics: Catalog file not found: " << catalog_file << std::endl;
        return false;
    }
    
    // Determine format from extension
    std::string ext = catalog_file.substr(catalog_file.find_last_of('.') + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    
    if (ext == "csv") {
        // Parse CSV
        std::string line;
        std::getline(file, line);  // Skip header
        
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            CatalogEntry entry;
            
            std::string field;
            int col = 0;
            while (std::getline(iss, field, ',')) {
                switch (col++) {
                    case 0: entry.id = field; break;
                    case 1: entry.name = field; break;
                    case 2: entry.type = field; break;
                    case 3: entry.ra = std::stod(field); break;
                    case 4: entry.dec = std::stod(field); break;
                    case 5: entry.distance_pc = std::stod(field); break;
                    case 6: entry.mass_solar = std::stod(field); break;
                    case 7: entry.radius_solar = std::stod(field); break;
                    case 8: entry.apparent_magnitude = std::stod(field); break;
                    case 9: entry.color_bv = std::stof(field); break;
                }
            }
            
            addCatalogEntry(entry);
        }
    }
#ifdef HAVE_NLOHMANN_JSON
    else if (ext == "json") {
        // Parse JSON
        nlohmann::json j;
        file >> j;
        
        for (const auto& item : j["objects"]) {
            CatalogEntry entry;
            entry.id = item.value("id", "");
            entry.name = item.value("name", "");
            entry.type = item.value("type", "star");
            entry.ra = item.value("ra", 0.0);
            entry.dec = item.value("dec", 0.0);
            entry.distance_pc = item.value("distance_pc", 0.0);
            entry.mass_solar = item.value("mass_solar", 1.0);
            entry.radius_solar = item.value("radius_solar", 1.0);
            entry.apparent_magnitude = item.value("apparent_magnitude", 0.0);
            entry.color_bv = item.value("color_bv", 0.6f);
            
            addCatalogEntry(entry);
        }
    }
#endif
    else {
        std::cerr << "AstroGraphics: Unknown catalog format: " << ext << std::endl;
        return false;
    }
    
    std::cout << "AstroGraphics: Loaded " << catalog_.size() << " objects from catalog" << std::endl;
    buildSpatialIndex();
    
    return true;
}

// ============================================================================
// Catalog Management
// ============================================================================

void AstroGraphics::addCatalogEntry(const CatalogEntry& entry) {
    size_t index = catalog_.size();
    catalog_.push_back(entry);
    catalog_index_[entry.id] = index;
}

CatalogEntry* AstroGraphics::findEntry(const std::string& id) {
    auto it = catalog_index_.find(id);
    if (it != catalog_index_.end()) {
        return &catalog_[it->second];
    }
    return nullptr;
}

std::vector<CatalogEntry*> AstroGraphics::searchNear(double ra, double dec, double radius_deg) {
    std::vector<CatalogEntry*> results;
    
    for (auto& entry : catalog_) {
        // Simple angular distance check
        double dra = entry.ra - ra;
        double ddec = entry.dec - dec;
        double dist = std::sqrt(dra*dra + ddec*ddec);
        
        if (dist <= radius_deg) {
            results.push_back(&entry);
        }
    }
    
    return results;
}

std::vector<CatalogEntry*> AstroGraphics::searchByName(const std::string& pattern) {
    std::vector<CatalogEntry*> results;
    
    try {
        std::regex regex(pattern, std::regex_constants::icase);
        
        for (auto& entry : catalog_) {
            if (std::regex_search(entry.name, regex) || 
                std::regex_search(entry.id, regex)) {
                results.push_back(&entry);
            }
        }
    } catch (const std::regex_error&) {
        // Fall back to simple substring match
        std::string pattern_lower = pattern;
        std::transform(pattern_lower.begin(), pattern_lower.end(), pattern_lower.begin(), ::tolower);
        
        for (auto& entry : catalog_) {
            std::string name_lower = entry.name;
            std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), ::tolower);
            
            if (name_lower.find(pattern_lower) != std::string::npos) {
                results.push_back(&entry);
            }
        }
    }
    
    return results;
}

// ============================================================================
// UQFF Integration
// ============================================================================

void AstroGraphics::setFieldOverlayConfig(const FieldOverlayConfig& config) {
    field_config_ = config;
}

void AstroGraphics::updateFieldData(const double* grid_data, int grid_dims[3], 
                                    double grid_origin[3], double cell_size) {
    // Store field grid data
    size_t total_cells = grid_dims[0] * grid_dims[1] * grid_dims[2];
    field_grid_.resize(total_cells);
    std::copy(grid_data, grid_data + total_cells, field_grid_.begin());
    
    std::copy(grid_dims, grid_dims + 3, field_dims_);
    std::copy(grid_origin, grid_origin + 3, field_origin_);
    field_cell_size_ = cell_size;
}

void AstroGraphics::calculateFieldAtObjects() {
    // Batch calculate F_U for all catalog objects
    // Would normally send to physics backend via IPC
    
    for (auto& entry : catalog_) {
        // Convert RA/Dec/Distance to Cartesian
        double ra_rad = entry.ra * M_PI / 180.0;
        double dec_rad = entry.dec * M_PI / 180.0;
        double d = entry.distance_pc * 3.086e16;  // pc to meters
        
        double x = d * std::cos(dec_rad) * std::cos(ra_rad);
        double y = d * std::cos(dec_rad) * std::sin(ra_rad);
        double z = d * std::sin(dec_rad);
        
        // Sample field at position (if grid available)
        if (!field_grid_.empty()) {
            // Convert to grid coordinates
            int ix = static_cast<int>((x - field_origin_[0]) / field_cell_size_);
            int iy = static_cast<int>((y - field_origin_[1]) / field_cell_size_);
            int iz = static_cast<int>((z - field_origin_[2]) / field_cell_size_);
            
            if (ix >= 0 && ix < field_dims_[0] &&
                iy >= 0 && iy < field_dims_[1] &&
                iz >= 0 && iz < field_dims_[2]) {
                size_t idx = iz * field_dims_[0] * field_dims_[1] + iy * field_dims_[0] + ix;
                entry.F_U = field_grid_[idx];
            }
        }
    }
}

// ============================================================================
// IPC-Based Field Calculation (Phase 4)
// ============================================================================

bool AstroGraphics::calculateFieldViaIPC(const std::string& system_name, double r, double t,
                                         double& F_U_out, double& compressed_g_out) {
    if (!physics_channel_ || !physics_channel_->is_connected()) {
        std::cerr << "AstroGraphics: No IPC channel to physics backend" << std::endl;
        return false;
    }
    
    // Build request
    UQFF::IPC::CalculateFieldRequest request;
    std::memset(&request, 0, sizeof(request));
    std::strncpy(request.system_name, system_name.c_str(), 
                 std::min(system_name.size(), size_t(63)));
    request.r = r;
    request.t = t;
    request.tn = 1.0;
    request.theta = 0.0;
    request.flags = 0;
    
    // Send request
    UQFF::IPC::MessageHeader header(UQFF::IPC::MessageType::CALCULATE_FIELD, 
                                    sizeof(UQFF::IPC::CalculateFieldRequest));
    
    if (!physics_channel_->send(header, &request)) {
        std::cerr << "AstroGraphics: Failed to send IPC request" << std::endl;
        return false;
    }
    
    // Receive response
    UQFF::IPC::MessageHeader response_header;
    std::vector<uint8_t> payload;
    
    if (!physics_channel_->receive(response_header, payload, 5000)) {  // 5s timeout
        std::cerr << "AstroGraphics: IPC response timeout" << std::endl;
        return false;
    }
    
    if (response_header.type == UQFF::IPC::MessageType::RESPONSE_DATA &&
        payload.size() >= sizeof(UQFF::IPC::CalculateFieldResponse)) {
        const auto* response = reinterpret_cast<const UQFF::IPC::CalculateFieldResponse*>(payload.data());
        F_U_out = response->F_U;
        compressed_g_out = response->compressed_g;
        return true;
    }
    
    return false;
}

void AstroGraphics::calculateAllFieldsViaIPC() {
    if (!physics_channel_ || !physics_channel_->is_connected()) {
        std::cerr << "AstroGraphics: No IPC channel - using local field grid" << std::endl;
        calculateFieldAtObjects();
        return;
    }
    
    std::cout << "AstroGraphics: Calculating fields for " << catalog_.size() 
              << " objects via IPC..." << std::endl;
    
    int success_count = 0;
    for (auto& entry : catalog_) {
        // Convert position to radial distance from galactic center
        double ra_rad = entry.ra * M_PI / 180.0;
        double dec_rad = entry.dec * M_PI / 180.0;
        double d = entry.distance_pc * 3.086e16;  // pc to meters
        
        double F_U = 0, compressed_g = 0;
        if (calculateFieldViaIPC(entry.id, d, 0.0, F_U, compressed_g)) {
            entry.F_U = F_U;
            entry.field_gradient = compressed_g;  // Store compressed_g in gradient field
            success_count++;
        }
    }
    
    std::cout << "AstroGraphics: Successfully calculated " << success_count 
              << " / " << catalog_.size() << " field values via IPC" << std::endl;
}

// ============================================================================
// Scene Control
// ============================================================================

void AstroGraphics::setSceneConfig(const SceneConfig& config) {
    scene_config_ = config;
}

void AstroGraphics::setCameraPosition(const double* position) {
    std::copy(position, position + 3, scene_config_.camera_position);
}

void AstroGraphics::setCameraTarget(const double* target) {
    std::copy(target, target + 3, scene_config_.camera_target);
}

void AstroGraphics::flyTo(const std::string& object_id, double duration_s) {
    CatalogEntry* entry = findEntry(object_id);
    if (!entry) {
        // Try name search
        auto results = searchByName(object_id);
        if (!results.empty()) {
            entry = results[0];
        }
    }
    
    if (entry) {
        // Convert RA/Dec to position
        double ra_rad = entry->ra * M_PI / 180.0;
        double dec_rad = entry->dec * M_PI / 180.0;
        double d = entry->distance_pc;
        
        double target[3] = {
            d * std::cos(dec_rad) * std::cos(ra_rad),
            d * std::cos(dec_rad) * std::sin(ra_rad),
            d * std::sin(dec_rad)
        };
        
        setCameraTarget(target);
        std::cout << "AstroGraphics: Flying to " << entry->name << " over " << duration_s << "s" << std::endl;
    } else {
        std::cerr << "AstroGraphics: Object not found: " << object_id << std::endl;
    }
}

void AstroGraphics::orbit(const std::string& object_id, double radius_pc, double speed) {
    (void)object_id;
    (void)radius_pc;
    (void)speed;
    // TODO: Implement orbital camera animation
}

// ============================================================================
// Rendering
// ============================================================================

void AstroGraphics::render(int eye_index, const float* view_matrix, const float* proj_matrix) {
    (void)eye_index;
    (void)view_matrix;
    (void)proj_matrix;
    
    if (!compositor_ || !render_resources_->initialized) {
        return;
    }
    
    // Update visible set based on camera
    updateVisibleSet();
    
    // Render objects by type
    renderStars();
    renderGalaxies();
    renderNebulae();
    
    // Render field overlay if enabled
    if (field_config_.enabled && !field_grid_.empty()) {
        renderFieldOverlay();
    }
    
    stats_.visible_objects = static_cast<int>(std::min(
        catalog_.size(), 
        static_cast<size_t>(scene_config_.max_visible_objects)
    ));
}

void AstroGraphics::renderOverlay() {
    if (scene_config_.show_object_labels) {
        renderLabels();
    }
}

// ============================================================================
// Selection
// ============================================================================

std::string AstroGraphics::pickObject(int screen_x, int screen_y) {
    (void)screen_x;
    (void)screen_y;
    // TODO: Implement ray-sphere intersection for picking
    return "";
}

void AstroGraphics::selectObject(const std::string& id) {
    selected_object_ = findEntry(id);
    if (selected_object_) {
        std::cout << "AstroGraphics: Selected " << selected_object_->name << std::endl;
    }
}

void AstroGraphics::clearSelection() {
    selected_object_ = nullptr;
}

// ============================================================================
// Data Streaming
// ============================================================================

void AstroGraphics::requestObjectData(const std::string& id) {
    // Async fetch from SIMBAD/NASA API
    // Would normally use CURL or similar
    
    if (data_callback_) {
        CatalogEntry* entry = findEntry(id);
        if (entry) {
            data_callback_(id, *entry);
        }
    }
}

// ============================================================================
// Helper Methods
// ============================================================================

void AstroGraphics::buildSpatialIndex() {
    // TODO: Build octree or BVH for efficient spatial queries
}

void AstroGraphics::updateVisibleSet() {
    // Frustum culling and LOD selection
    // TODO: Implement based on camera position and view frustum
}

void AstroGraphics::renderStars() {
    if (!compositor_) return;
    
    // Render stars as points or billboards
    size_t star_count = 0;
    for (const auto& entry : catalog_) {
        if (entry.type == StarCatalog::TYPE_STAR) {
            star_count++;
        }
    }
    
    stats_.triangles_rendered += static_cast<int>(star_count * 2);  // Placeholder
}

void AstroGraphics::renderGalaxies() {
    if (!compositor_) return;
    
    // Render galaxies as textured sprites or meshes
    for (const auto& entry : catalog_) {
        if (entry.type == StarCatalog::TYPE_GALAXY) {
            stats_.triangles_rendered += 100;  // Placeholder
        }
    }
}

void AstroGraphics::renderNebulae() {
    if (!compositor_) return;
    
    // Render nebulae as volumetric or billboard
    for (const auto& entry : catalog_) {
        if (entry.type == StarCatalog::TYPE_NEBULA) {
            stats_.triangles_rendered += 500;  // Placeholder
        }
    }
}

void AstroGraphics::renderFieldOverlay() {
    // Pass field data to compositor for visualization
    if (compositor_ && !field_grid_.empty()) {
        compositor_->updateFieldData(field_grid_.data(), field_dims_, field_cell_size_);
    }
}

void AstroGraphics::renderLabels() {
    if (!compositor_) return;
    
    // Render object labels as 2D text
    for (const auto& entry : catalog_) {
        if (entry.apparent_magnitude < scene_config_.max_magnitude) {
            // Would calculate screen position and render text
        }
    }
}

double AstroGraphics::calculateApparentMagnitude(const CatalogEntry& entry) {
    // m = M + 5 * log10(d/10)  where d is in parsecs
    if (entry.distance_pc <= 0) return entry.apparent_magnitude;
    
    double distance_modulus = 5.0 * std::log10(entry.distance_pc / 10.0);
    return entry.absolute_magnitude + distance_modulus;
}

float AstroGraphics::colorIndexToRGB(float bv, float* rgb) {
    // Convert B-V color index to RGB
    // Simplified - actual conversion would use blackbody spectrum
    
    float t;  // Temperature-like parameter
    
    if (bv < 0) {
        // Hot blue stars
        rgb[0] = 0.6f;
        rgb[1] = 0.7f;
        rgb[2] = 1.0f;
    } else if (bv < 0.4) {
        // White stars
        t = bv / 0.4f;
        rgb[0] = 0.6f + 0.4f * t;
        rgb[1] = 0.7f + 0.3f * t;
        rgb[2] = 1.0f;
    } else if (bv < 0.8) {
        // Yellow stars
        t = (bv - 0.4f) / 0.4f;
        rgb[0] = 1.0f;
        rgb[1] = 1.0f - 0.2f * t;
        rgb[2] = 0.9f - 0.4f * t;
    } else if (bv < 1.5) {
        // Orange stars
        t = (bv - 0.8f) / 0.7f;
        rgb[0] = 1.0f;
        rgb[1] = 0.8f - 0.3f * t;
        rgb[2] = 0.5f - 0.3f * t;
    } else {
        // Red stars
        rgb[0] = 1.0f;
        rgb[1] = 0.5f;
        rgb[2] = 0.2f;
    }
    
    return bv;  // Return original value
}

} // namespace VR
