/**
 * @file csv_body_reader.h
 * @brief CSV Reader for bodies_*.csv files from APIFetch.py
 * 
 * Parses astronomical body data exported by APIFetch.py in the format:
 *   name,mass,radius,distance,B_field,SFR,z,v_rot,...
 * 
 * Usage:
 *   #include "csv_body_reader.h"
 *   
 *   auto bodies = CSVBodyReader::read_latest(".");
 *   for (const auto& body : bodies) {
 *       SystemParams sys = body.to_system_params();
 *       // Process...
 *   }
 * 
 * Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
 * Framework: UQFF Star-Magic v3.0
 * Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
 */

#ifndef CSV_BODY_READER_H
#define CSV_BODY_READER_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <regex>
#include <stdexcept>
#include <cmath>
#include <optional>

namespace UQFF {

/**
 * Represents a celestial body parsed from CSV.
 */
struct CelestialBodyCSV {
    std::string name;
    double mass = 0.0;           // kg
    double radius = 0.0;         // m
    double distance = 0.0;       // m (from observer)
    double B_field = 0.0;        // T (magnetic field)
    double SFR = 0.0;            // M_sun/yr (star formation rate)
    double z = 0.0;              // redshift
    double v_rot = 0.0;          // m/s (rotation velocity)
    double luminosity = 0.0;     // W
    double temperature = 0.0;    // K
    double metallicity = 0.0;    // [Fe/H]
    double age = 0.0;            // s
    double spin_period = 0.0;    // s
    double inclination = 0.0;    // radians
    
    // Coordinates
    double ra = 0.0;             // degrees (right ascension)
    double dec = 0.0;            // degrees (declination)
    double galactic_l = 0.0;     // degrees
    double galactic_b = 0.0;     // degrees
    
    // Source metadata
    std::string source_api;      // "SIMBAD", "NASA", "Grok", etc.
    std::string object_type;     // "star", "galaxy", "magnetar", etc.
    std::string timestamp;       // ISO timestamp of fetch
    
    // Additional parameters (key-value)
    std::map<std::string, double> extra_params;
    
    /**
     * Convert to SystemParams struct for UQFF calculations.
     * Requires SystemParams to be defined (forward declaration).
     */
    template<typename SystemParams>
    SystemParams to_system_params() const {
        SystemParams sys;
        sys.name = name;
        sys.M = mass;
        sys.r = (radius > 0) ? radius : distance;  // Use radius or distance
        sys.v = v_rot;
        sys.B0 = B_field;
        sys.t = 0;  // Initial time
        sys.SFR = SFR;
        sys.z = z;
        return sys;
    }
    
    /**
     * Check if essential fields are populated.
     */
    bool is_valid() const {
        return !name.empty() && (mass > 0 || radius > 0 || distance > 0);
    }
};


/**
 * CSV Reader for bodies_*.csv files.
 */
class CSVBodyReader {
public:
    /**
     * Read a specific CSV file.
     * 
     * @param filepath Path to CSV file
     * @return Vector of parsed bodies
     * @throws std::runtime_error if file cannot be opened
     */
    static std::vector<CelestialBodyCSV> read(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open CSV file: " + filepath);
        }
        
        std::vector<CelestialBodyCSV> bodies;
        std::string line;
        std::vector<std::string> headers;
        
        // Read header line
        if (std::getline(file, line)) {
            headers = parse_csv_line(line);
            // Normalize headers to lowercase
            for (auto& h : headers) {
                std::transform(h.begin(), h.end(), h.begin(), ::tolower);
            }
        }
        
        // Read data lines
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            auto values = parse_csv_line(line);
            if (values.empty()) continue;
            
            CelestialBodyCSV body = parse_body(headers, values);
            if (body.is_valid()) {
                bodies.push_back(body);
            }
        }
        
        return bodies;
    }
    
    /**
     * Find and read the most recent bodies_*.csv file in a directory.
     * 
     * @param directory Directory to search (default: current)
     * @return Vector of parsed bodies
     * @throws std::runtime_error if no matching file found
     */
    static std::vector<CelestialBodyCSV> read_latest(const std::string& directory = ".") {
        std::filesystem::path dir(directory);
        std::string latest_file;
        std::filesystem::file_time_type latest_time;
        
        // Pattern: bodies_YYYYMMDD_HHMMSS.csv
        std::regex pattern(R"(bodies_\d{8}_\d{6}\.csv)");
        
        for (const auto& entry : std::filesystem::directory_iterator(dir)) {
            if (!entry.is_regular_file()) continue;
            
            std::string filename = entry.path().filename().string();
            if (std::regex_match(filename, pattern)) {
                auto file_time = entry.last_write_time();
                if (latest_file.empty() || file_time > latest_time) {
                    latest_file = entry.path().string();
                    latest_time = file_time;
                }
            }
        }
        
        if (latest_file.empty()) {
            // Try generic bodies.csv
            std::filesystem::path generic = dir / "bodies.csv";
            if (std::filesystem::exists(generic)) {
                latest_file = generic.string();
            } else {
                throw std::runtime_error("No bodies_*.csv file found in: " + directory);
            }
        }
        
        return read(latest_file);
    }
    
    /**
     * List all bodies_*.csv files in a directory.
     */
    static std::vector<std::string> list_csv_files(const std::string& directory = ".") {
        std::vector<std::string> files;
        std::filesystem::path dir(directory);
        std::regex pattern(R"(bodies.*\.csv)");
        
        for (const auto& entry : std::filesystem::directory_iterator(dir)) {
            if (!entry.is_regular_file()) continue;
            
            std::string filename = entry.path().filename().string();
            if (std::regex_search(filename, pattern)) {
                files.push_back(entry.path().string());
            }
        }
        
        std::sort(files.begin(), files.end());
        return files;
    }

private:
    /**
     * Parse a CSV line handling quoted fields.
     */
    static std::vector<std::string> parse_csv_line(const std::string& line) {
        std::vector<std::string> result;
        std::string field;
        bool in_quotes = false;
        
        for (size_t i = 0; i < line.length(); ++i) {
            char c = line[i];
            
            if (c == '"') {
                in_quotes = !in_quotes;
            } else if (c == ',' && !in_quotes) {
                // Trim whitespace
                size_t start = field.find_first_not_of(" \t");
                size_t end = field.find_last_not_of(" \t");
                if (start != std::string::npos) {
                    result.push_back(field.substr(start, end - start + 1));
                } else {
                    result.push_back("");
                }
                field.clear();
            } else {
                field += c;
            }
        }
        
        // Add last field
        size_t start = field.find_first_not_of(" \t\r\n");
        size_t end = field.find_last_not_of(" \t\r\n");
        if (start != std::string::npos) {
            result.push_back(field.substr(start, end - start + 1));
        } else {
            result.push_back("");
        }
        
        return result;
    }
    
    /**
     * Parse a body from header-value pairs.
     * Supports both APIFetch.py format (mass, radius, distance) and 
     * bodies.csv format (Ms, Rs, Rb, Ts_surface, omega_s, Bs_avg, etc.)
     */
    static CelestialBodyCSV parse_body(const std::vector<std::string>& headers,
                                        const std::vector<std::string>& values) {
        CelestialBodyCSV body;
        
        for (size_t i = 0; i < headers.size() && i < values.size(); ++i) {
            const std::string& key = headers[i];
            const std::string& val = values[i];
            
            if (val.empty() || val == "nan" || val == "NaN" || val == "None") {
                continue;
            }
            
            try {
                // === Name/Identifier ===
                if (key == "name" || key == "object" || key == "identifier") {
                    body.name = val;
                }
                // === Mass (APIFetch or bodies.csv format) ===
                else if (key == "mass" || key == "m" || key == "mass_kg" || key == "ms") {
                    body.mass = parse_number(val);
                }
                // === Radius (APIFetch or bodies.csv format) ===
                else if (key == "radius" || key == "r" || key == "radius_m" || key == "rs") {
                    body.radius = parse_number(val);
                }
                // === Distance/Boundary (APIFetch or bodies.csv Rb) ===
                else if (key == "distance" || key == "d" || key == "distance_m" || key == "rb") {
                    body.distance = parse_number(val);
                }
                // === Magnetic field (bodies.csv: Bs_avg) ===
                else if (key == "b_field" || key == "b" || key == "magnetic_field" || key == "bs_avg") {
                    body.B_field = parse_number(val);
                }
                // === Temperature (bodies.csv: Ts_surface) ===
                else if (key == "temperature" || key == "t" || key == "temp" || key == "t_eff" || key == "ts_surface") {
                    body.temperature = parse_number(val);
                }
                // === Spin rate (bodies.csv: omega_s maps to spin via period) ===
                else if (key == "omega_s") {
                    double omega = parse_number(val);
                    if (omega > 0) {
                        body.spin_period = 2.0 * 3.14159265358979 / omega;  // T = 2π/ω
                    }
                    body.extra_params["omega_s"] = omega;
                }
                // === SCm density (bodies.csv unique) ===
                else if (key == "scm_density") {
                    body.extra_params["scm_density"] = parse_number(val);
                }
                // === QUA - Universal Aether charge (bodies.csv unique) ===
                else if (key == "qua") {
                    body.extra_params["qua"] = parse_number(val);
                }
                // === Pcore - Core pressure ratio (bodies.csv unique) ===
                else if (key == "pcore") {
                    body.extra_params["pcore"] = parse_number(val);
                }
                // === PSCm - SCm pressure ratio (bodies.csv unique) ===
                else if (key == "pscm") {
                    body.extra_params["pscm"] = parse_number(val);
                }
                // === omega_c - Core spin (bodies.csv unique) ===
                else if (key == "omega_c") {
                    body.extra_params["omega_c"] = parse_number(val);
                }
                // === Standard fields ===
                else if (key == "sfr" || key == "star_formation_rate") {
                    body.SFR = parse_number(val);
                } else if (key == "z" || key == "redshift") {
                    body.z = parse_number(val);
                } else if (key == "v_rot" || key == "v" || key == "velocity" || key == "rotation_velocity") {
                    body.v_rot = parse_number(val);
                } else if (key == "luminosity" || key == "l" || key == "luminosity_w") {
                    body.luminosity = parse_number(val);
                } else if (key == "metallicity" || key == "feh" || key == "[fe/h]") {
                    body.metallicity = parse_number(val);
                } else if (key == "age" || key == "age_s") {
                    body.age = parse_number(val);
                } else if (key == "spin_period" || key == "period" || key == "p_rot") {
                    body.spin_period = parse_number(val);
                } else if (key == "inclination" || key == "inc" || key == "i") {
                    body.inclination = parse_number(val);
                } else if (key == "ra") {
                    body.ra = parse_number(val);
                } else if (key == "dec") {
                    body.dec = parse_number(val);
                } else if (key == "l" || key == "galactic_l" || key == "glon") {
                    body.galactic_l = parse_number(val);
                } else if (key == "galactic_b" || key == "glat") {
                    body.galactic_b = parse_number(val);
                } else if (key == "source" || key == "api" || key == "source_api") {
                    body.source_api = val;
                } else if (key == "type" || key == "object_type" || key == "otype") {
                    body.object_type = val;
                } else if (key == "timestamp" || key == "fetch_time") {
                    body.timestamp = val;
                } else {
                    // Store as extra parameter
                    double num_val = parse_number(val);
                    if (!std::isnan(num_val)) {
                        body.extra_params[key] = num_val;
                    }
                }
            } catch (...) {
                // Skip unparseable values
                continue;
            }
        }
        
        return body;
    }
    
    /**
     * Parse a number that may be in scientific notation.
     */
    static double parse_number(const std::string& s) {
        if (s.empty()) return std::nan("");
        
        try {
            // Handle common formats: 1.23e+45, 1.23E-45, 1.23, etc.
            return std::stod(s);
        } catch (...) {
            return std::nan("");
        }
    }
};

} // namespace UQFF

#endif // CSV_BODY_READER_H
