/**
 * @file vulkan_compositor.h
 * @brief Vulkan GPU Compositor for UQFF Star-Magic VR Runtime
 * 
 * GPU rendering pipeline for VR compositor:
 * - Vulkan instance/device management
 * - Stereo rendering to OpenXR swapchains
 * - Physics field visualization (color-coded gradients)
 * - Point cloud rendering for astronomical data
 * - Compute shaders for real-time field calculations
 * 
 * Supported Render Modes:
 * - STANDARD: Solid rendering with lighting
 * - WIREFRAME: Field lines and vector visualization
 * - VOLUME: Volumetric field density
 * - POINT_CLOUD: Raw astronomical data points
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 * Phase: 3 - VR Runtime Scaffold
 */

#ifndef VULKAN_COMPOSITOR_H
#define VULKAN_COMPOSITOR_H

// Conditional compilation for Vulkan availability
#ifdef HAVE_VULKAN
#define VK_USE_PLATFORM_WIN32_KHR
#include <vulkan/vulkan.h>
#endif

#include <string>
#include <vector>
#include <memory>
#include <array>
#include <functional>

namespace VR {

// Forward declarations
class OpenXRSession;
struct SwapchainInfo;

/**
 * @enum RenderMode
 * @brief Rendering style for field visualization
 */
enum class RenderMode {
    Standard,       // PBR with lighting
    Wireframe,      // Field lines
    Volume,         // Volumetric density
    PointCloud,     // Raw data points
    Hybrid          // Combination
};

/**
 * @struct FieldVisualizationConfig
 * @brief Configuration for physics field rendering
 */
struct FieldVisualizationConfig {
    // Field type to visualize
    enum class FieldType {
        F_U_Total,      // Total unified field
        Ug1_Magnetic,   // Magnetic dipole
        Ug2_Charge,     // Charge-reactivity
        Ug3_String,     // String rotation
        Ug4_Vacuum,     // Vacuum concentration
        Um_Magnetism,   // Magnetism
        Ubi_Buoyancy    // Buoyancy opposition
    } field_type = FieldType::F_U_Total;
    
    // Color mapping
    float min_field_value = 0.0f;
    float max_field_value = 1e-8f;
    std::array<float, 4> low_color = {0.0f, 0.0f, 1.0f, 1.0f};   // Blue
    std::array<float, 4> high_color = {1.0f, 0.0f, 0.0f, 1.0f};  // Red
    bool use_log_scale = true;
    
    // Field line rendering
    bool show_field_lines = true;
    float line_density = 100.0f;        // Lines per unit volume
    float line_spacing = 0.1f;          // Meters
    
    // Vector arrows
    bool show_vectors = false;
    float vector_scale = 1.0f;
    float vector_spacing = 0.5f;        // Meters
    
    // Isosurfaces
    bool show_isosurface = false;
    std::vector<float> isovalues = {1e-10, 1e-9, 1e-8};
};

/**
 * @struct AstronomicalObject
 * @brief Renderable astronomical object
 */
struct AstronomicalObject {
    std::string name;
    
    // Position in local space (meters from origin)
    double position[3] = {0, 0, 0};
    double scale = 1.0;                 // Visual scale factor
    
    // Physical properties (for field calculations)
    double mass = 0;                    // kg
    double radius = 0;                  // m
    double magnetic_field = 0;          // T
    double spin_period = 0;             // s
    
    // Visual properties
    std::string texture_path;
    std::array<float, 4> base_color = {1, 1, 1, 1};
    float emission = 0.0f;              // For stars
    
    // Calculated field at object
    double local_F_U = 0;
};

/**
 * @struct RenderStats
 * @brief Per-frame rendering statistics
 */
struct RenderStats {
    double gpu_time_ms = 0;
    uint32_t draw_calls = 0;
    uint32_t triangles = 0;
    uint32_t points = 0;
    double memory_used_mb = 0;
};

/**
 * @class VulkanCompositor
 * @brief GPU rendering engine using Vulkan
 */
class VulkanCompositor {
public:
    VulkanCompositor();
    ~VulkanCompositor();
    
    // Initialization
    bool initialize(OpenXRSession* xr_session);
    void shutdown();
    
    // Frame rendering
    void beginFrame(int eye_index);
    void endFrame(int eye_index);
    
    // Scene composition
    void clearBackground(const float* color);
    void setViewProjection(const float* view_matrix, const float* proj_matrix);
    
    // Field visualization
    void setFieldConfig(const FieldVisualizationConfig& config);
    void updateFieldData(const double* field_grid, int grid_size[3], double cell_size);
    void renderFieldVisualization();
    
    // Astronomical objects
    void addAstronomicalObject(const AstronomicalObject& obj);
    void updateObjectPosition(const std::string& name, const double* position);
    void removeAstronomicalObject(const std::string& name);
    void renderAstronomicalObjects();
    
    // Point cloud (raw data)
    void setPointCloud(const float* positions, const float* colors, uint32_t count);
    void renderPointCloud();
    
    // UI overlay
    void beginUIOverlay();
    void renderText(const std::string& text, float x, float y, float size);
    void renderRect(float x, float y, float width, float height, const float* color);
    void endUIOverlay();
    
    // Render mode
    void setRenderMode(RenderMode mode) { render_mode_ = mode; }
    RenderMode getRenderMode() const { return render_mode_; }
    
    // Statistics
    RenderStats getStats() const { return stats_; }
    
    // Vulkan accessors (for advanced use)
#ifdef HAVE_VULKAN
    VkDevice getDevice() const { return device_; }
    VkQueue getQueue() const { return graphics_queue_; }
    VkCommandBuffer getCurrentCommandBuffer() const { return current_command_buffer_; }
#else
    void* getDevice() const { return device_; }
    void* getQueue() const { return graphics_queue_; }
    void* getCurrentCommandBuffer() const { return current_command_buffer_; }
#endif

private:
    // Vulkan resources
#ifdef HAVE_VULKAN
    VkInstance instance_ = VK_NULL_HANDLE;
    VkPhysicalDevice physical_device_ = VK_NULL_HANDLE;
    VkDevice device_ = VK_NULL_HANDLE;
    VkQueue graphics_queue_ = VK_NULL_HANDLE;
    VkCommandPool command_pool_ = VK_NULL_HANDLE;
    VkCommandBuffer current_command_buffer_ = VK_NULL_HANDLE;
    
    // Pipelines
    VkPipelineLayout pipeline_layout_ = VK_NULL_HANDLE;
    VkPipeline standard_pipeline_ = VK_NULL_HANDLE;
    VkPipeline wireframe_pipeline_ = VK_NULL_HANDLE;
    VkPipeline point_cloud_pipeline_ = VK_NULL_HANDLE;
    VkPipeline field_viz_pipeline_ = VK_NULL_HANDLE;
    
    // Descriptor sets
    VkDescriptorSetLayout desc_set_layout_ = VK_NULL_HANDLE;
    VkDescriptorPool desc_pool_ = VK_NULL_HANDLE;
    
    // Render pass
    VkRenderPass render_pass_ = VK_NULL_HANDLE;
    std::vector<VkFramebuffer> framebuffers_;
    
    // Buffers
    VkBuffer uniform_buffer_ = VK_NULL_HANDLE;
    VkDeviceMemory uniform_memory_ = VK_NULL_HANDLE;
    VkBuffer vertex_buffer_ = VK_NULL_HANDLE;
    VkDeviceMemory vertex_memory_ = VK_NULL_HANDLE;
    VkBuffer field_buffer_ = VK_NULL_HANDLE;
    VkDeviceMemory field_memory_ = VK_NULL_HANDLE;
#else
    // Stubs when Vulkan not available
    void* instance_ = nullptr;
    void* physical_device_ = nullptr;
    void* device_ = nullptr;
    void* graphics_queue_ = nullptr;
    void* command_pool_ = nullptr;
    void* current_command_buffer_ = nullptr;
    void* pipeline_layout_ = nullptr;
    void* standard_pipeline_ = nullptr;
    void* wireframe_pipeline_ = nullptr;
    void* point_cloud_pipeline_ = nullptr;
    void* field_viz_pipeline_ = nullptr;
    void* desc_set_layout_ = nullptr;
    void* desc_pool_ = nullptr;
    void* render_pass_ = nullptr;
    std::vector<void*> framebuffers_;
    void* uniform_buffer_ = nullptr;
    void* uniform_memory_ = nullptr;
    void* vertex_buffer_ = nullptr;
    void* vertex_memory_ = nullptr;
    void* field_buffer_ = nullptr;
    void* field_memory_ = nullptr;
#endif
    
    // OpenXR integration
    OpenXRSession* xr_session_ = nullptr;
    
    // State
    RenderMode render_mode_ = RenderMode::Standard;
    FieldVisualizationConfig field_config_;
    std::vector<AstronomicalObject> objects_;
    RenderStats stats_;
    
    // Current frame state
    int current_eye_ = 0;
    uint32_t current_image_index_ = 0;
    
    // Helper methods
    bool createDevice();
    bool createPipelines();
    bool createBuffers();
    bool createFramebuffers(const SwapchainInfo& swapchain);
    void updateUniforms(const float* view, const float* proj);
    void recordFieldVizCommands();
    uint32_t findMemoryType(uint32_t type_filter, uint32_t properties);
};

} // namespace VR

#endif // VULKAN_COMPOSITOR_H
