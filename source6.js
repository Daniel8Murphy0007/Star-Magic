// source6.js - UQFF Module 6 with 3D Graphics, Model Loading, and Advanced Rendering
// Converted from Source6.cpp
// © 2025 Daniel T. Murphy - Star-Magic UQFF Project

const fs = require('fs');
const path = require('path');

// ============================================================================
// CONSTANTS - Physics and Material Properties
// ============================================================================

const PI = 3.14159265358979323846;
const c = 3e8;           // Speed of light (m/s)
const G = 6.67430e-11;   // Gravitational constant
const Omega_g = 7.3e-16; // Galactic rotation
const Mbh = 8.15e36;     // Black hole mass (kg)
const dg = 2.55e20;      // Galactic distance (m)

// Material properties
const v_SCm = 0.99 * c;
const rho_A = 1e-23;
const rho_sw = 8e-21;
const v_sw = 5e5;
const QA = 1e-10;
const Qs = 0.0;

// Coupling constants
const kappa = 0.0005;
const alpha = 0.001;
const gamma = 0.00005;
const delta_sw = 0.01;
const epsilon_sw = 0.001;
const delta_def = 0.01;
const HSCm = 1.0;
const UUA = 1.0;
const eta = 1e-22;

// Advanced constants
const k1 = 1.5;
const k2 = 1.2;
const k3 = 1.8;
const k4 = 2.0;
const beta_i = 0.6;
const rho_v = 6e-27;
const C_concentration = 1.0;
const f_feedback = 0.1;
const num_strings = 1e9;
const Ts00 = 1.27e3 + 1.11e7;

// Metric tensor (Minkowski)
const g_mu_nu = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, -1.0, 0.0, 0.0],
    [0.0, 0.0, -1.0, 0.0],
    [0.0, 0.0, 0.0, -1.0]
];

// ============================================================================
// CELESTIAL BODY CLASS
// ============================================================================

class CelestialBody {
    constructor(params = {}) {
        this.name = params.name || "Unnamed";
        this.Ms = params.Ms || 0.0;           // Star mass (kg)
        this.Rs = params.Rs || 0.0;           // Star radius (m)
        this.Rb = params.Rb || 0.0;           // Bubble radius (m)
        this.Ts_surface = params.Ts_surface || 0.0;  // Surface temperature (K)
        this.omega_s = params.omega_s || 0.0; // Rotation rate (rad/s)
        this.Bs_avg = params.Bs_avg || 0.0;   // Magnetic field (T)
        this.SCm_density = params.SCm_density || 0.0;
        this.QUA = params.QUA || 0.0;         // Aether charge
        this.Pcore = params.Pcore || 0.0;     // Core pressure (Pa)
        this.PSCm = params.PSCm || 0.0;       // SCm pressure (Pa)
        this.omega_c = params.omega_c || 0.0; // Cycle frequency (rad/s)
    }

    toJSON() {
        return {
            name: this.name,
            Ms: this.Ms,
            Rs: this.Rs,
            Rb: this.Rb,
            Ts_surface: this.Ts_surface,
            omega_s: this.omega_s,
            Bs_avg: this.Bs_avg,
            SCm_density: this.SCm_density,
            QUA: this.QUA,
            Pcore: this.Pcore,
            PSCm: this.PSCm,
            omega_c: this.omega_c
        };
    }

    static fromJSON(json) {
        return new CelestialBody(json);
    }
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

function stepFunction(r, Rb) {
    return r <= Rb ? 1.0 : 0.0;
}

function computeEreact(t, rho_SCm, v_SCm, rho_A, kappa) {
    const decay = Math.exp(-kappa * t);
    return 0.5 * rho_SCm * v_SCm * v_SCm * decay + rho_A * kappa * t;
}

function computeMuS(t, Bs, omega_c, Rs) {
    const Bj = Math.cos(omega_c * t);
    return (4.0 * PI * Rs * Rs * Rs / 3.0) * Bs * Bj;
}

function computeGradMsR(Ms, Rs) {
    if (Rs <= 0.0) throw new Error("Invalid Rs value");
    return G * Ms / (Rs * Rs);
}

function computeBj(t, omega_c) {
    return Math.cos(omega_c * t);
}

function computeOmegaST(t, omega_s, omega_c) {
    return omega_s * (1.0 + 0.1 * Math.cos(omega_c * t));
}

function computeMuJ(t, omega_c, Rs) {
    const Bj = computeBj(t, omega_c);
    return (4.0 * PI * Rs * Rs * Rs / 3.0) * Bj;
}

// ============================================================================
// PHYSICS COMPUTATION FUNCTIONS
// ============================================================================

/**
 * Compute Ug1 - Internal dipole field with defects
 */
function computeUg1(body, r, t, tn, alpha, delta_def, k1) {
    const mu_s = computeMuS(t, body.Bs_avg, body.omega_c, body.Rs);
    const grad_Ms = computeGradMsR(body.Ms, body.Rs);
    const step = stepFunction(r, body.Rb);
    const defect_term = delta_def * body.Pcore * Math.cos(PI * tn);
    const decay = Math.exp(-alpha * t);
    
    if (r <= 0.0) throw new Error("Invalid r value in computeUg1");
    
    return k1 * mu_s * grad_Ms * step * (1.0 + defect_term) * decay / (r * r);
}

/**
 * Compute Ug2 - Outer field bubble with stellar wind
 */
function computeUg2(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa) {
    const Ereact = computeEreact(t, body.SCm_density, v_SCm, rho_A, kappa);
    const omega_s_t = computeOmegaST(t, body.omega_s, body.omega_c);
    const wind_mod = 1.0 + delta_sw * v_sw;
    const step = stepFunction(r, body.Rb);
    
    if (r <= 0.0) throw new Error("Invalid r value in computeUg2");
    
    return k2 * QA * body.QUA * Ereact * omega_s_t * wind_mod * HSCm * step / (r * r);
}

/**
 * Compute Ug3 - Magnetic strings with planetary core
 */
function computeUg3(body, r, t, tn, theta, rho_A, kappa, k3) {
    const mu_j = computeMuJ(t, body.omega_c, body.Rs);
    const core_term = body.Pcore + body.PSCm * Math.cos(PI * tn);
    const theta_dep = Math.sin(theta);
    
    if (r <= 0.0) throw new Error("Invalid r value in computeUg3");
    
    return k3 * mu_j * core_term * theta_dep / (r * r * r);
}

/**
 * Compute Um - Magnetism with multiple strings
 */
function computeUm(body, t, tn, rj, gamma, rho_A, kappa, num_strings) {
    const mu_j = computeMuJ(t, body.omega_c, body.Rs);
    const cycle = Math.cos(PI * tn);
    
    if (rj <= 0.0) throw new Error("Invalid rj value in computeUm");
    
    return gamma * num_strings * mu_j * cycle / (rj * rj);
}

/**
 * Compute Ug4 - Star-black hole interactions
 */
function computeUg4(t, tn, rho_v, C_concentration, Mbh, dg, alpha, f_feedback, k4) {
    if (dg <= 0.0) throw new Error("Invalid dg value");
    
    const decay = Math.exp(-alpha * t);
    const cycle = Math.cos(PI * tn);
    
    return k4 * rho_v * C_concentration * Mbh / dg * decay * cycle * (1 + f_feedback);
}

/**
 * Compute Ubi - Universal buoyancy with wind modulation
 */
function computeUbi(Ugi, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn) {
    if (dg <= 0.0) throw new Error("Invalid dg value");
    
    const wind_mod = 1.0 + epsilon_sw * rho_sw;
    return -beta_i * Ugi * Omega_g * Mbh / dg * wind_mod * UUA * Math.cos(PI * tn);
}

/**
 * Compute A_mu_nu - Metric tensor with temperature modulation
 */
function computeAMuNu(tn, eta, Ts00) {
    const A = g_mu_nu.map(row => [...row]); // Deep copy
    const mod = eta * Ts00 * Math.cos(PI * tn);
    
    for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
            A[i][j] += mod;
        }
    }
    
    return A;
}

/**
 * Compute FU - Unified Field Strength
 */
function computeFU(body, r, t, tn, theta) {
    try {
        // Individual Ug components
        const Ug1 = computeUg1(body, r, t, tn, alpha, delta_def, k1);
        const Ug2 = computeUg2(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa);
        const Ug3 = computeUg3(body, r, t, tn, theta, rho_A, kappa, k3);
        const Ug4 = computeUg4(t, tn, rho_v, C_concentration, Mbh, dg, alpha, f_feedback, k4);
        const sumUgi = Ug1 + Ug2 + Ug3 + Ug4;
        
        // Buoyancy components
        const Ubi1 = computeUbi(Ug1, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        const Ubi2 = computeUbi(Ug2, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        const Ubi3 = computeUbi(Ug3, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        const Ubi4 = computeUbi(Ug4, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
        const sumUbi = Ubi1 + Ubi2 + Ubi3 + Ubi4;
        
        // Magnetism
        const Um = computeUm(body, t, tn, body.Rb, gamma, rho_A, kappa, num_strings);
        
        // Metric tensor contribution
        const A = computeAMuNu(tn, eta, Ts00);
        let A_scalar = 0.0;
        for (let i = 0; i < 4; i++) {
            A_scalar += A[i][i];
        }
        
        return sumUgi + sumUbi + Um + A_scalar;
    } catch (error) {
        console.error(`Error in computeFU for ${body.name}:`, error.message);
        return 0.0;
    }
}

// ============================================================================
// 3D GRAPHICS CLASSES
// ============================================================================

/**
 * 3DObject - Represents a 3D mesh with vertices, normals, and indices
 * Provides abstraction for different rendering backends (WebGL, etc.)
 */
class ThreeDObject {
    constructor() {
        this.vertices = [];  // Flat array of x,y,z coordinates
        this.normals = [];   // Flat array of nx,ny,nz normals
        this.indices = [];   // Array of vertex indices
        this.textureId = null;
        
        // WebGL buffers (if using WebGL)
        this.VAO = null;
        this.VBO = null;
        this.EBO = null;
    }

    /**
     * Setup WebGL buffers (browser environment)
     */
    setupWebGL(gl) {
        if (typeof WebGLRenderingContext === 'undefined') {
            console.warn("WebGL not available in this environment");
            return;
        }

        // Create VAO
        this.VAO = gl.createVertexArray();
        gl.bindVertexArray(this.VAO);

        // Create and bind VBO
        this.VBO = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.VBO);
        
        // Interleave vertices and normals
        const vertexData = [];
        for (let i = 0; i < this.vertices.length / 3; i++) {
            vertexData.push(
                this.vertices[i * 3],
                this.vertices[i * 3 + 1],
                this.vertices[i * 3 + 2],
                this.normals[i * 3],
                this.normals[i * 3 + 1],
                this.normals[i * 3 + 2]
            );
        }
        
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertexData), gl.STATIC_DRAW);

        // Position attribute
        gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 6 * 4, 0);
        gl.enableVertexAttribArray(0);

        // Normal attribute
        gl.vertexAttribPointer(1, 3, gl.FLOAT, false, 6 * 4, 3 * 4);
        gl.enableVertexAttribArray(1);

        // Create and bind EBO
        this.EBO = gl.createBuffer();
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.EBO);
        gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(this.indices), gl.STATIC_DRAW);

        gl.bindVertexArray(null);
    }

    /**
     * Render using WebGL
     */
    renderWebGL(gl) {
        if (!this.VAO) {
            console.warn("WebGL buffers not initialized");
            return;
        }

        gl.bindVertexArray(this.VAO);
        gl.drawElements(gl.TRIANGLES, this.indices.length, gl.UNSIGNED_INT, 0);
        gl.bindVertexArray(null);
    }

    /**
     * Generic render method (can be extended for different backends)
     */
    render(backend = 'webgl', context = null) {
        switch (backend) {
            case 'webgl':
                if (context) this.renderWebGL(context);
                break;
            default:
                console.warn(`Rendering backend '${backend}' not implemented`);
        }
    }
}

/**
 * ToolPath - CNC-like tool path with positions and speeds
 */
class ToolPath {
    constructor() {
        this.points = [];  // Flat array of x,y,z coordinates
        this.speeds = [];  // Array of speeds for each point
    }

    /**
     * Import tool path from CSV file
     */
    importFromCSV(filename) {
        try {
            const content = fs.readFileSync(filename, 'utf-8');
            const lines = content.split('\n');
            
            this.points = [];
            this.speeds = [];
            
            for (const line of lines) {
                const trimmed = line.trim();
                if (!trimmed) continue;
                
                const parts = trimmed.split(',');
                if (parts.length >= 4) {
                    this.points.push(
                        parseFloat(parts[0]),
                        parseFloat(parts[1]),
                        parseFloat(parts[2])
                    );
                    this.speeds.push(parseFloat(parts[3]));
                }
            }
            
            console.log(`Loaded tool path with ${this.speeds.length} points from ${filename}`);
        } catch (error) {
            console.error(`Error importing tool path from ${filename}:`, error.message);
        }
    }

    /**
     * Export tool path to binary format
     */
    exportToBinary(filename) {
        try {
            const data = {
                numPoints: this.speeds.length,
                points: this.points,
                speeds: this.speeds
            };
            
            fs.writeFileSync(filename, JSON.stringify(data));
            console.log(`Exported tool path to ${filename}`);
        } catch (error) {
            console.error(`Error exporting tool path to ${filename}:`, error.message);
        }
    }
}

/**
 * SimulationEntity - Entity in a physics simulation
 */
class SimulationEntity {
    constructor() {
        this.position = [0, 0, 0];
        this.velocity = [0, 0, 0];
        this.model = new ThreeDObject();
    }

    /**
     * Update entity position based on velocity
     */
    update(dt) {
        this.position[0] += this.velocity[0] * dt;
        this.position[1] += this.velocity[1] * dt;
        this.position[2] += this.velocity[2] * dt;
    }
}

// ============================================================================
// MODEL LOADING (OBJ FORMAT)
// ============================================================================

/**
 * MeshData - Container for mesh data
 */
class MeshData {
    constructor() {
        this.vertices = [];   // Array of [x, y, z] or flat array
        this.normals = [];    // Array of [x, y, z] or flat array
        this.texCoords = [];  // Array of [u, v] or flat array
        this.indices = [];    // Array of vertex indices
    }
}

/**
 * Load OBJ file into MeshData
 */
function loadOBJ(filepath, mesh) {
    try {
        const content = fs.readFileSync(filepath, 'utf-8');
        const lines = content.split('\n');
        
        const tempVertices = [];
        const tempNormals = [];
        const tempTexCoords = [];
        
        for (const line of lines) {
            const trimmed = line.trim();
            if (!trimmed || trimmed.startsWith('#')) continue;
            
            const parts = trimmed.split(/\s+/);
            const token = parts[0];
            
            if (token === 'v' && parts.length >= 4) {
                // Vertex
                tempVertices.push([
                    parseFloat(parts[1]),
                    parseFloat(parts[2]),
                    parseFloat(parts[3])
                ]);
            } else if (token === 'vn' && parts.length >= 4) {
                // Normal
                tempNormals.push([
                    parseFloat(parts[1]),
                    parseFloat(parts[2]),
                    parseFloat(parts[3])
                ]);
            } else if (token === 'vt' && parts.length >= 3) {
                // Texture coordinate
                tempTexCoords.push([
                    parseFloat(parts[1]),
                    parseFloat(parts[2])
                ]);
            } else if (token === 'f' && parts.length >= 4) {
                // Face (assuming triangles)
                for (let i = 1; i <= 3; i++) {
                    const face = parts[i].replace(/\//g, ' ');
                    const indices = face.split(/\s+/).map(x => parseInt(x));
                    
                    const v = indices[0];
                    const t = indices[1] || 0;
                    const n = indices[2] || 0;
                    
                    if (v > 0 && v <= tempVertices.length) {
                        mesh.vertices.push(...tempVertices[v - 1]);
                    }
                    if (t > 0 && t <= tempTexCoords.length) {
                        mesh.texCoords.push(...tempTexCoords[t - 1]);
                    }
                    if (n > 0 && n <= tempNormals.length) {
                        mesh.normals.push(...tempNormals[n - 1]);
                    }
                    
                    mesh.indices.push(mesh.indices.length);
                }
            }
        }
        
        console.log(`Loaded OBJ file: ${filepath} (${mesh.indices.length} indices)`);
        return true;
    } catch (error) {
        console.error(`Error loading OBJ file ${filepath}:`, error.message);
        return false;
    }
}

/**
 * Export mesh to OBJ file
 */
function exportOBJ(filepath, mesh) {
    try {
        let content = '';
        
        // Write vertices
        for (let i = 0; i < mesh.vertices.length; i += 3) {
            content += `v ${mesh.vertices[i]} ${mesh.vertices[i + 1]} ${mesh.vertices[i + 2]}\n`;
        }
        
        // Write texture coordinates
        for (let i = 0; i < mesh.texCoords.length; i += 2) {
            content += `vt ${mesh.texCoords[i]} ${mesh.texCoords[i + 1]}\n`;
        }
        
        // Write normals
        for (let i = 0; i < mesh.normals.length; i += 3) {
            content += `vn ${mesh.normals[i]} ${mesh.normals[i + 1]} ${mesh.normals[i + 2]}\n`;
        }
        
        // Write faces
        for (let i = 0; i < mesh.indices.length; i += 3) {
            const i1 = mesh.indices[i] + 1;
            const i2 = mesh.indices[i + 1] + 1;
            const i3 = mesh.indices[i + 2] + 1;
            content += `f ${i1}/${i1}/${i1} ${i2}/${i2}/${i2} ${i3}/${i3}/${i3}\n`;
        }
        
        fs.writeFileSync(filepath, content);
        console.log(`Exported OBJ file: ${filepath}`);
    } catch (error) {
        console.error(`Error exporting OBJ file ${filepath}:`, error.message);
    }
}

// ============================================================================
// TEXTURE LOADING
// ============================================================================

/**
 * Load texture (browser environment with Image API)
 */
function loadTexture(filepath, gl = null) {
    if (typeof Image === 'undefined' || !gl) {
        console.warn("Texture loading requires browser environment with WebGL context");
        return null;
    }

    const texture = gl.createTexture();
    const image = new Image();
    
    image.onload = () => {
        gl.bindTexture(gl.TEXTURE_2D, texture);
        
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image);
        gl.generateMipmap(gl.TEXTURE_2D);
        
        console.log(`Loaded texture: ${filepath}`);
    };
    
    image.onerror = () => {
        console.error(`Failed to load texture: ${filepath}`);
    };
    
    image.src = filepath;
    return texture;
}

// ============================================================================
// SHADER CLASS
// ============================================================================

/**
 * Shader - GLSL shader program wrapper
 */
class Shader {
    constructor(vertexPath = null, fragmentPath = null, gl = null) {
        this.ID = null;
        this.gl = gl;
        
        if (vertexPath && fragmentPath && gl) {
            this.loadFromFiles(vertexPath, fragmentPath);
        }
    }

    /**
     * Load shaders from files
     */
    loadFromFiles(vertexPath, fragmentPath) {
        try {
            const vertexCode = fs.readFileSync(vertexPath, 'utf-8');
            const fragmentCode = fs.readFileSync(fragmentPath, 'utf-8');
            this.compileShaders(vertexCode, fragmentCode);
        } catch (error) {
            console.error("Error loading shader files:", error.message);
        }
    }

    /**
     * Compile shader from source code
     */
    compileShaders(vertexCode, fragmentCode) {
        if (!this.gl) {
            console.error("WebGL context not available");
            return;
        }

        const gl = this.gl;

        // Compile vertex shader
        const vertexShader = gl.createShader(gl.VERTEX_SHADER);
        gl.shaderSource(vertexShader, vertexCode);
        gl.compileShader(vertexShader);
        
        if (!gl.getShaderParameter(vertexShader, gl.COMPILE_STATUS)) {
            console.error("Vertex shader compilation error:", gl.getShaderInfoLog(vertexShader));
            return;
        }

        // Compile fragment shader
        const fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
        gl.shaderSource(fragmentShader, fragmentCode);
        gl.compileShader(fragmentShader);
        
        if (!gl.getShaderParameter(fragmentShader, gl.COMPILE_STATUS)) {
            console.error("Fragment shader compilation error:", gl.getShaderInfoLog(fragmentShader));
            return;
        }

        // Link program
        this.ID = gl.createProgram();
        gl.attachShader(this.ID, vertexShader);
        gl.attachShader(this.ID, fragmentShader);
        gl.linkProgram(this.ID);
        
        if (!gl.getProgramParameter(this.ID, gl.LINK_STATUS)) {
            console.error("Shader program linking error:", gl.getProgramInfoLog(this.ID));
            return;
        }

        // Cleanup
        gl.deleteShader(vertexShader);
        gl.deleteShader(fragmentShader);
        
        console.log("Shader compiled successfully");
    }

    /**
     * Use this shader program
     */
    use() {
        if (this.gl && this.ID) {
            this.gl.useProgram(this.ID);
        }
    }

    /**
     * Set mat4 uniform
     */
    setMat4(name, mat) {
        if (this.gl && this.ID) {
            const location = this.gl.getUniformLocation(this.ID, name);
            this.gl.uniformMatrix4fv(location, false, mat);
        }
    }

    /**
     * Set vec3 uniform
     */
    setVec3(name, x, y, z) {
        if (this.gl && this.ID) {
            const location = this.gl.getUniformLocation(this.ID, name);
            this.gl.uniform3f(location, x, y, z);
        }
    }

    /**
     * Set float uniform
     */
    setFloat(name, value) {
        if (this.gl && this.ID) {
            const location = this.gl.getUniformLocation(this.ID, name);
            this.gl.uniform1f(location, value);
        }
    }
}

// ============================================================================
// CAMERA CLASS
// ============================================================================

/**
 * Camera - FPS-style camera with view/projection matrices
 */
class Camera {
    constructor() {
        this.position = [0, 0, 3];
        this.front = [0, 0, -1];
        this.up = [0, 1, 0];
        this.right = [1, 0, 0];
        this.worldUp = [0, 1, 0];
        
        this.yaw = -90.0;
        this.pitch = 0.0;
        this.fov = 45.0;
        this.movementSpeed = 2.5;
        this.mouseSensitivity = 0.1;
        
        this.updateCameraVectors();
    }

    /**
     * Get view matrix
     */
    getViewMatrix() {
        // lookAt implementation
        const target = [
            this.position[0] + this.front[0],
            this.position[1] + this.front[1],
            this.position[2] + this.front[2]
        ];
        
        return this.lookAt(this.position, target, this.up);
    }

    /**
     * Get projection matrix
     */
    getProjectionMatrix(aspect, near = 0.1, far = 100.0) {
        return this.perspective(this.fov * Math.PI / 180.0, aspect, near, far);
    }

    /**
     * Process keyboard input
     */
    processKeyboard(direction, deltaTime) {
        const velocity = this.movementSpeed * deltaTime;
        
        if (direction === 'FORWARD') {
            this.position[0] += this.front[0] * velocity;
            this.position[1] += this.front[1] * velocity;
            this.position[2] += this.front[2] * velocity;
        }
        if (direction === 'BACKWARD') {
            this.position[0] -= this.front[0] * velocity;
            this.position[1] -= this.front[1] * velocity;
            this.position[2] -= this.front[2] * velocity;
        }
        if (direction === 'LEFT') {
            this.position[0] -= this.right[0] * velocity;
            this.position[1] -= this.right[1] * velocity;
            this.position[2] -= this.right[2] * velocity;
        }
        if (direction === 'RIGHT') {
            this.position[0] += this.right[0] * velocity;
            this.position[1] += this.right[1] * velocity;
            this.position[2] += this.right[2] * velocity;
        }
    }

    /**
     * Process mouse movement
     */
    processMouseMovement(xoffset, yoffset, constrainPitch = true) {
        xoffset *= this.mouseSensitivity;
        yoffset *= this.mouseSensitivity;
        
        this.yaw += xoffset;
        this.pitch += yoffset;
        
        if (constrainPitch) {
            if (this.pitch > 89.0) this.pitch = 89.0;
            if (this.pitch < -89.0) this.pitch = -89.0;
        }
        
        this.updateCameraVectors();
    }

    /**
     * Update camera vectors from yaw/pitch
     */
    updateCameraVectors() {
        const yawRad = this.yaw * Math.PI / 180.0;
        const pitchRad = this.pitch * Math.PI / 180.0;
        
        this.front[0] = Math.cos(yawRad) * Math.cos(pitchRad);
        this.front[1] = Math.sin(pitchRad);
        this.front[2] = Math.sin(yawRad) * Math.cos(pitchRad);
        
        this.normalize(this.front);
        
        this.cross(this.right, this.front, this.worldUp);
        this.normalize(this.right);
        
        this.cross(this.up, this.right, this.front);
        this.normalize(this.up);
    }

    // Vector math helpers
    normalize(v) {
        const len = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (len > 0) {
            v[0] /= len;
            v[1] /= len;
            v[2] /= len;
        }
    }

    cross(out, a, b) {
        out[0] = a[1] * b[2] - a[2] * b[1];
        out[1] = a[2] * b[0] - a[0] * b[2];
        out[2] = a[0] * b[1] - a[1] * b[0];
    }

    lookAt(eye, center, up) {
        const f = [
            center[0] - eye[0],
            center[1] - eye[1],
            center[2] - eye[2]
        ];
        this.normalize(f);
        
        const s = [0, 0, 0];
        this.cross(s, f, up);
        this.normalize(s);
        
        const u = [0, 0, 0];
        this.cross(u, s, f);
        
        return [
            s[0], u[0], -f[0], 0,
            s[1], u[1], -f[1], 0,
            s[2], u[2], -f[2], 0,
            -s[0] * eye[0] - s[1] * eye[1] - s[2] * eye[2],
            -u[0] * eye[0] - u[1] * eye[1] - u[2] * eye[2],
            f[0] * eye[0] + f[1] * eye[1] + f[2] * eye[2],
            1
        ];
    }

    perspective(fovy, aspect, near, far) {
        const f = 1.0 / Math.tan(fovy / 2);
        const nf = 1.0 / (near - far);
        
        return [
            f / aspect, 0, 0, 0,
            0, f, 0, 0,
            0, 0, (far + near) * nf, -1,
            0, 0, 2 * far * near * nf, 0
        ];
    }
}

// ============================================================================
// PLUGIN SYSTEM
// ============================================================================

/**
 * SIMPlugin - Dynamic plugin loader using ES6 dynamic imports
 */
class SIMPlugin {
    constructor(modulePath = null) {
        this.module = null;
        this.modulePath = modulePath;
        
        if (modulePath) {
            this.load(modulePath);
        }
    }

    /**
     * Load plugin module
     */
    async load(modulePath) {
        try {
            this.modulePath = modulePath;
            this.module = await import(modulePath);
            
            if (!this.module.playAPI) {
                throw new Error("Plugin missing 'playAPI' function");
            }
            
            console.log(`Loaded plugin: ${modulePath}`);
        } catch (error) {
            console.error(`Error loading plugin ${modulePath}:`, error.message);
            throw error;
        }
    }

    /**
     * Call plugin's playAPI function
     */
    async playAPI(...args) {
        if (!this.module) {
            throw new Error("Plugin not loaded");
        }
        
        return await this.module.playAPI(...args);
    }
}

// ============================================================================
// SIMULATION FUNCTIONS
// ============================================================================

/**
 * Simulate quasar jet with Navier-Stokes and MUGE gravity
 * Note: Requires FluidSolver and MUGESystem from source5.js
 */
function simulateQuasarJet(initialVelocity, outputFile = "", FluidSolver = null, computeResonanceMUGE = null) {
    try {
        if (!FluidSolver) {
            console.warn("FluidSolver not available (requires source5.js)");
            return;
        }

        const solver = new FluidSolver();
        solver.addJetForce(initialVelocity / 10.0);

        // Create Sagittarius A* MUGE system
        const sagA = {
            name: "Sagittarius A*",
            I: 1e23,
            A: 2.813e30,
            omega1: 1e-5,
            omega2: -1e-5,
            Vsys: 3.552e45,
            vexp: 5e6,
            t: 3.786e14,
            ffluid: 3.465e-8,
            r: 1e12
        };

        let uqffG = 0;
        if (computeResonanceMUGE) {
            uqffG = computeResonanceMUGE(sagA, {});
        }

        console.log(`Simulating quasar jet with Navier-Stokes (10 steps) using UQFF g=${uqffG}...`);
        
        for (let step = 0; step < 10; step++) {
            solver.step(uqffG / 1e30);
        }

        solver.printVelocityField();

        if (outputFile) {
            // Write velocity to CSV if needed
            console.log(`Would write velocity field to ${outputFile}`);
        }
    } catch (error) {
        console.error("Error in simulateQuasarJet:", error.message);
    }
}

/**
 * Print summary statistics
 */
function printSummaryStats(values, name) {
    if (values.length === 0) return;
    
    let min = values[0];
    let max = values[0];
    let sum = 0.0;
    
    for (const val of values) {
        if (val < min) min = val;
        if (val > max) max = val;
        sum += val;
    }
    
    const mean = sum / values.length;
    console.log(`${name} summary - Min: ${min}, Max: ${max}, Mean: ${mean}`);
}

/**
 * Load celestial bodies from JSON or CSV file
 */
function loadBodies(filename) {
    try {
        const ext = path.extname(filename).toLowerCase();
        
        if (ext === '.json') {
            const data = JSON.parse(fs.readFileSync(filename, 'utf-8'));
            return data.map(item => CelestialBody.fromJSON(item));
        } else if (ext === '.csv') {
            const content = fs.readFileSync(filename, 'utf-8');
            const lines = content.split('\n');
            const bodies = [];
            
            for (let i = 1; i < lines.length; i++) { // Skip header
                const line = lines[i].trim();
                if (!line) continue;
                
                const parts = line.split(',');
                if (parts.length >= 12) {
                    bodies.push(new CelestialBody({
                        name: parts[0],
                        Ms: parseFloat(parts[1]),
                        Rs: parseFloat(parts[2]),
                        Rb: parseFloat(parts[3]),
                        Ts_surface: parseFloat(parts[4]),
                        omega_s: parseFloat(parts[5]),
                        Bs_avg: parseFloat(parts[6]),
                        SCm_density: parseFloat(parts[7]),
                        QUA: parseFloat(parts[8]),
                        Pcore: parseFloat(parts[9]),
                        PSCm: parseFloat(parts[10]),
                        omega_c: parseFloat(parts[11])
                    }));
                }
            }
            
            return bodies;
        }
        
        return [];
    } catch (error) {
        console.error(`Error loading bodies from ${filename}:`, error.message);
        return [];
    }
}

/**
 * Get default celestial bodies
 */
function getDefaultBodies() {
    return [
        new CelestialBody({
            name: "Sun",
            Ms: 1.989e30,
            Rs: 6.96e8,
            Rb: 1.5e11,
            Ts_surface: 5778.0,
            omega_s: 2.87e-6,
            Bs_avg: 1e-4,
            SCm_density: 1e-20,
            QUA: 1e-10,
            Pcore: 2.65e16,
            PSCm: 1e10,
            omega_c: 1e-8
        }),
        new CelestialBody({
            name: "Earth",
            Ms: 5.972e24,
            Rs: 6.371e6,
            Rb: 1.5e11,
            Ts_surface: 288.0,
            omega_s: 7.29e-5,
            Bs_avg: 3.1e-5,
            SCm_density: 1e-23,
            QUA: 5e-11,
            Pcore: 3.6e11,
            PSCm: 1e8,
            omega_c: 1e-9
        }),
        new CelestialBody({
            name: "Jupiter",
            Ms: 1.898e27,
            Rs: 6.9911e7,
            Rb: 7.78e11,
            Ts_surface: 165.0,
            omega_s: 1.76e-4,
            Bs_avg: 4.3e-4,
            SCm_density: 5e-22,
            QUA: 8e-11,
            Pcore: 4e12,
            PSCm: 5e9,
            omega_c: 5e-9
        }),
        new CelestialBody({
            name: "Neptune",
            Ms: 1.024e26,
            Rs: 2.4622e7,
            Rb: 4.5e12,
            Ts_surface: 72.0,
            omega_s: 1.08e-4,
            Bs_avg: 2.7e-5,
            SCm_density: 3e-23,
            QUA: 4e-11,
            Pcore: 7e11,
            PSCm: 2e8,
            omega_c: 2e-9
        })
    ];
}

// ============================================================================
// SELF-EXPANDING PHYSICS FRAMEWORK (2.0-Enhanced)
// ============================================================================

/**
 * Base class for dynamic physics terms
 */
class PhysicsTerm {
    constructor(name, description) {
        this.name = name;
        this.description = description;
        this.parameters = new Map();
        this.enabled = true;
    }

    compute(state) {
        throw new Error("compute() must be implemented by subclass");
    }

    setParameter(key, value) {
        this.parameters.set(key, value);
    }

    getParameter(key, defaultValue = 0) {
        return this.parameters.has(key) ? this.parameters.get(key) : defaultValue;
    }

    exportState() {
        return {
            name: this.name,
            description: this.description,
            parameters: Object.fromEntries(this.parameters),
            enabled: this.enabled
        };
    }
}

/**
 * Dark Matter Halo Term
 */
class DarkMatterHaloTerm extends PhysicsTerm {
    constructor(haloMass, haloRadius) {
        super("DarkMatterHalo", "NFW dark matter halo contribution");
        this.setParameter("halo_mass", haloMass);
        this.setParameter("halo_radius", haloRadius);
    }

    compute(state) {
        const M_halo = this.getParameter("halo_mass");
        const r_s = this.getParameter("halo_radius");
        const r = state.r || 1e13;
        
        if (r <= 0 || r_s <= 0) return 0.0;
        
        const x = r / r_s;
        const rho_0 = M_halo / (4 * PI * r_s * r_s * r_s * (Math.log(2) - 0.5));
        
        return -4 * PI * G * rho_0 * r_s * r_s * Math.log(1 + x) / r;
    }
}

/**
 * Vacuum Energy Term
 */
class VacuumEnergyTerm extends PhysicsTerm {
    constructor(lambda) {
        super("VacuumEnergy", "Cosmological vacuum energy contribution");
        this.setParameter("lambda", lambda);
    }

    compute(state) {
        const lambda = this.getParameter("lambda");
        const r = state.r || 1e13;
        
        return lambda * r / 3.0;
    }
}

/**
 * Self-Expanding UQFF Module 6
 */
class UQFFModule6JS {
    constructor() {
        this.metadata = {
            version: "2.0-Enhanced",
            created: new Date().toISOString(),
            author: "Star-Magic UQFF",
            description: "3D Graphics, Model Loading, Advanced Rendering"
        };
        
        this.dynamicTerms = [];
        this.dynamicParameters = new Map();
        this.learningRate = 0.01;
        this.enableLogging = false;
        this.computationHistory = [];
    }

    /**
     * Register a new dynamic physics term
     */
    registerDynamicTerm(term) {
        if (!(term instanceof PhysicsTerm)) {
            throw new Error("Term must be an instance of PhysicsTerm");
        }
        
        this.dynamicTerms.push(term);
        
        if (this.enableLogging) {
            console.log(`Registered dynamic term: ${term.name}`);
        }
    }

    /**
     * Set dynamic parameter
     */
    setDynamicParameter(key, value) {
        this.dynamicParameters.set(key, value);
        
        if (this.enableLogging) {
            console.log(`Set dynamic parameter: ${key} = ${value}`);
        }
    }

    /**
     * Get dynamic parameter
     */
    getDynamicParameter(key, defaultValue = 0) {
        return this.dynamicParameters.has(key) ? this.dynamicParameters.get(key) : defaultValue;
    }

    /**
     * Compute FU with dynamic terms
     */
    computeFUWithDynamics(body, r, t, tn, theta) {
        // Core validated calculation
        const coreFU = computeFU(body, r, t, tn, theta);
        
        // Add dynamic terms
        let dynamicContribution = 0.0;
        
        if (this.dynamicTerms.length > 0) {
            const state = { body, r, t, tn, theta };
            
            for (const term of this.dynamicTerms) {
                if (term.enabled) {
                    dynamicContribution += term.compute(state);
                }
            }
        }
        
        const total = coreFU + dynamicContribution;
        
        if (this.enableLogging) {
            this.computationHistory.push({
                timestamp: Date.now(),
                body: body.name,
                coreFU,
                dynamicContribution,
                total
            });
        }
        
        return total;
    }

    /**
     * Export module state
     */
    exportState(filename = null) {
        const state = {
            metadata: this.metadata,
            dynamicTerms: this.dynamicTerms.map(term => term.exportState()),
            dynamicParameters: Object.fromEntries(this.dynamicParameters),
            learningRate: this.learningRate,
            historyLength: this.computationHistory.length
        };
        
        if (filename) {
            fs.writeFileSync(filename, JSON.stringify(state, null, 2));
            console.log(`Exported module state to ${filename}`);
        }
        
        return state;
    }

    /**
     * Import module state
     */
    importState(filename) {
        try {
            const data = JSON.parse(fs.readFileSync(filename, 'utf-8'));
            
            this.metadata = data.metadata;
            this.learningRate = data.learningRate;
            this.dynamicParameters = new Map(Object.entries(data.dynamicParameters || {}));
            
            console.log(`Imported module state from ${filename}`);
        } catch (error) {
            console.error(`Error importing state from ${filename}:`, error.message);
        }
    }

    /**
     * Set learning rate for optimization
     */
    setLearningRate(rate) {
        this.learningRate = rate;
    }

    /**
     * Enable/disable logging
     */
    setEnableLogging(enable) {
        this.enableLogging = enable;
    }

    /**
     * Get computation history
     */
    getComputationHistory() {
        return this.computationHistory;
    }

    /**
     * Clear computation history
     */
    clearHistory() {
        this.computationHistory = [];
    }
}

// ============================================================================
// EXPORTS
// ============================================================================

module.exports = {
    // Constants
    PI, c, G, Omega_g, Mbh, dg,
    v_SCm, rho_A, rho_sw, v_sw, QA, Qs,
    kappa, alpha, gamma, delta_sw, epsilon_sw, delta_def, HSCm, UUA, eta,
    k1, k2, k3, k4, beta_i, rho_v, C_concentration, f_feedback, num_strings,
    Ts00, g_mu_nu,
    
    // Classes
    CelestialBody,
    ThreeDObject,
    ToolPath,
    SimulationEntity,
    MeshData,
    Shader,
    Camera,
    SIMPlugin,
    PhysicsTerm,
    DarkMatterHaloTerm,
    VacuumEnergyTerm,
    UQFFModule6JS,
    
    // Helper functions
    stepFunction,
    computeEreact,
    computeMuS,
    computeGradMsR,
    computeBj,
    computeOmegaST,
    computeMuJ,
    
    // Physics functions
    computeUg1,
    computeUg2,
    computeUg3,
    computeUm,
    computeUg4,
    computeUbi,
    computeAMuNu,
    computeFU,
    
    // 3D & Model functions
    loadOBJ,
    exportOBJ,
    loadTexture,
    
    // Simulation functions
    simulateQuasarJet,
    printSummaryStats,
    loadBodies,
    getDefaultBodies
};

console.log("source6.js loaded - UQFF Module 6 with 3D Graphics & Advanced Rendering");
