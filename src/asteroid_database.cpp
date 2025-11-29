/**
 * @file asteroid_database.cpp
 * @brief Implementazione database asteroidi locale
 */

#include "ioccultcalc/asteroid_database.h"
#include "ioccultcalc/data_manager.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <ctime>
#include <iomanip>

namespace ioccultcalc {

// ============================================================================
// AsteroidDatabase implementation
// ============================================================================

AsteroidDatabase::AsteroidDatabase(const std::string& dbPath) 
    : dbPath_(dbPath) {
    
    if (dbPath_.empty()) {
        dbPath_ = DataManager::instance().getAsteroidDatabasePath();
    }
    
    lastUpdateDate_ = "";
}

std::string AsteroidDatabase::getDefaultPath() {
    return DataManager::instance().getAsteroidDatabasePath();
}

bool AsteroidDatabase::loadFromFile() {
    std::ifstream file(dbPath_);
    if (!file.is_open()) {
        return false;
    }
    
    try {
        nlohmann::json j;
        file >> j;
        
        // Load metadata
        if (j.contains("metadata")) {
            auto meta = j["metadata"];
            if (meta.contains("last_update")) {
                lastUpdateDate_ = meta["last_update"];
            }
            if (meta.contains("source")) {
                stats_.source = meta["source"];
            }
        }
        
        // Load asteroids
        if (j.contains("asteroids")) {
            for (auto& [key, value] : j["asteroids"].items()) {
                int number = std::stoi(key);
                AsteroidProperties props = parseMPCExtendedJson(value);
                props.number = number;
                asteroids_[number] = props;
            }
        }
        
        updateStats();
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error loading database: " << e.what() << std::endl;
        return false;
    }
}

bool AsteroidDatabase::loadFromFile(const std::string& path) {
    std::string oldPath = dbPath_;
    dbPath_ = path;
    bool result = loadFromFile();
    if (!result) {
        dbPath_ = oldPath;  // Restore on failure
    }
    return result;
}

void AsteroidDatabase::addAsteroid(const AsteroidProperties& props) {
    asteroids_[props.number] = props;
    updateStats();
}

bool AsteroidDatabase::saveToFile() {
    try {
        nlohmann::json j;
        
        // Metadata
        j["metadata"] = {
            {"version", "1.0"},
            {"source", stats_.source},
            {"last_update", getCurrentDate()},
            {"total_asteroids", stats_.total_asteroids},
            {"with_diameter", stats_.with_diameter},
            {"with_albedo", stats_.with_albedo}
        };
        
        // Asteroids
        nlohmann::json asteroidsJson;
        for (const auto& pair : asteroids_) {
            int num = pair.first;
            const auto& props = pair.second;
            
            asteroidsJson[std::to_string(num)] = {
                {"number", props.number},
                {"name", props.name},
                {"designation", props.designation},
                {"diameter", props.diameter},
                {"H", props.H},
                {"albedo", props.albedo},
                {"a", props.a},
                {"e", props.e},
                {"i", props.i},
                {"rotation_period", props.rotation_period},
                {"orbit_class", props.orbit_class},
                {"spectral_type", props.spectral_type}
            };
        }
        j["asteroids"] = asteroidsJson;
        
        // Write to file
        std::ofstream file(dbPath_);
        if (!file.is_open()) {
            return false;
        }
        
        file << std::setw(2) << j << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error saving database: " << e.what() << std::endl;
        return false;
    }
}

bool AsteroidDatabase::saveToFile(const std::string& path) {
    std::string oldPath = dbPath_;
    dbPath_ = path;
    bool result = saveToFile();
    dbPath_ = oldPath;  // Always restore
    return result;
}

AsteroidProperties AsteroidDatabase::parseMPCExtendedJson(const nlohmann::json& j) {
    AsteroidProperties props;
    
    // MPC Extended Format Fields
    // Number - Extract from "(123)" format
    if (j.contains("Number")) {
        std::string numStr = j["Number"];
        // Remove parentheses "(1)" -> "1"
        numStr.erase(std::remove(numStr.begin(), numStr.end(), '('), numStr.end());
        numStr.erase(std::remove(numStr.begin(), numStr.end(), ')'), numStr.end());
        try {
            props.number = std::stoi(numStr);
        } catch (...) {}
    }
    
    // Name
    if (j.contains("Name")) {
        props.name = j["Name"];
    }
    
    // Designation
    if (j.contains("Principal_desig")) {
        props.designation = j["Principal_desig"];
    }
    
    // H magnitude (absolute magnitude)
    if (j.contains("H")) {
        props.H = j["H"];
    }
    
    // Orbital elements (for filtering)
    if (j.contains("a")) props.a = j["a"];
    if (j.contains("e")) props.e = j["e"];
    if (j.contains("i")) props.i = j["i"];
    
    // Orbit type (MBA, NEA, etc.)
    if (j.contains("Orbit_type")) {
        props.orbit_class = j["Orbit_type"];
    }
    
    // Additional useful fields
    if (j.contains("G")) {
        // Slope parameter (for magnitude calculation)
        // Not directly used but good to preserve
    }
    
    // Estimate diameter from H magnitude if not provided
    // Using standard formula: D = 1329 * 10^(-H/5) / sqrt(albedo)
    // Assuming typical albedo 0.15 for C-type asteroids
    if (!props.has_diameter && props.H > 0) {
        double albedo = 0.15; // Default albedo
        props.diameter = 1329.0 * std::pow(10.0, -props.H / 5.0) / std::sqrt(albedo);
        props.has_diameter = true;
    }
    
    return props;
}

std::vector<AsteroidProperties> AsteroidDatabase::query(const AsteroidRange& range) const {
    std::vector<AsteroidProperties> results;
    
    // Get list of candidate numbers
    std::vector<int> numbers = range.getAsteroidList();
    
    for (int num : numbers) {
        auto it = asteroids_.find(num);
        if (it != asteroids_.end()) {
            if (range.matches(it->second)) {
                results.push_back(it->second);
            }
        }
    }
    
    return results;
}

std::vector<int> AsteroidDatabase::queryNumbers(const AsteroidRange& range) const {
    std::vector<int> results;
    
    std::vector<int> numbers = range.getAsteroidList();
    
    for (int num : numbers) {
        auto it = asteroids_.find(num);
        if (it != asteroids_.end()) {
            if (range.matches(it->second)) {
                results.push_back(num);
            }
        }
    }
    
    return results;
}

AsteroidProperties AsteroidDatabase::getProperties(int number) const {
    auto it = asteroids_.find(number);
    if (it != asteroids_.end()) {
        return it->second;
    }
    
    // Return default
    AsteroidProperties props;
    props.number = number;
    return props;
}

bool AsteroidDatabase::hasAsteroid(int number) const {
    return asteroids_.find(number) != asteroids_.end();
}

int AsteroidDatabase::count(const AsteroidRange& range) const {
    return queryNumbers(range).size();
}

DatabaseStats AsteroidDatabase::getStats() const {
    return stats_;
}

void AsteroidDatabase::updateStats() {
    stats_.total_asteroids = asteroids_.size();
    stats_.with_diameter = 0;
    stats_.with_albedo = 0;
    stats_.with_spectral_type = 0;
    stats_.with_rotation_period = 0;
    stats_.numbered = 0;
    stats_.unnumbered = 0;
    
    for (const auto& pair : asteroids_) {
        const auto& props = pair.second;
        
        if (props.has_diameter) stats_.with_diameter++;
        if (props.has_albedo) stats_.with_albedo++;
        if (!props.spectral_type.empty()) stats_.with_spectral_type++;
        if (props.rotation_period > 0) stats_.with_rotation_period++;
        if (props.number > 0) stats_.numbered++;
        else stats_.unnumbered++;
    }
    
    stats_.last_update = lastUpdateDate_;
}

std::string AsteroidDatabase::getCurrentDate() const {
    time_t now = time(nullptr);
    struct tm* tm_now = gmtime(&now);
    
    std::ostringstream oss;
    oss << std::put_time(tm_now, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

bool AsteroidDatabase::needsUpdate(int maxAgeDays) const {
    if (lastUpdateDate_.empty()) {
        return true;
    }
    
    // Simple check: compare with current date
    // For now, just check if database file exists and its age
    return DataManager::instance().getFileAgeDays(dbPath_) > maxAgeDays;
}

size_t AsteroidDatabase::getMemoryUsage() const {
    // Rough estimate: ~1KB per asteroid
    return asteroids_.size() * 1024;
}

void AsteroidDatabase::clear() {
    asteroids_.clear();
    updateStats();
}

bool AsteroidDatabase::exportToJson(const AsteroidRange& range, const std::string& outputPath) const {
    try {
        auto results = query(range);
        
        nlohmann::json j;
        j["metadata"] = {
            {"version", "1.0"},
            {"source", stats_.source},
            {"export_date", getCurrentDate()},
            {"total_asteroids", results.size()}
        };
        
        nlohmann::json asteroidsJson;
        for (const auto& props : results) {
            asteroidsJson[std::to_string(props.number)] = {
                {"number", props.number},
                {"name", props.name},
                {"diameter", props.diameter},
                {"H", props.H},
                {"albedo", props.albedo},
                {"a", props.a},
                {"e", props.e},
                {"i", props.i},
                {"orbit_class", props.orbit_class},
                {"spectral_type", props.spectral_type}
            };
        }
        j["asteroids"] = asteroidsJson;
        
        std::ofstream file(outputPath);
        if (!file.is_open()) {
            return false;
        }
        
        file << std::setw(2) << j << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error exporting: " << e.what() << std::endl;
        return false;
    }
}

int AsteroidDatabase::importFromJson(const std::string& inputPath, bool merge) {
    try {
        std::ifstream file(inputPath);
        if (!file.is_open()) {
            std::cerr << "Cannot open file: " << inputPath << std::endl;
            return 0;
        }
        
        nlohmann::json j;
        file >> j;
        
        if (!merge) {
            asteroids_.clear();
        }
        
        int count = 0;
        
        // Check if MPC Extended Format (array of objects)
        if (j.is_array()) {
            std::cout << "Parsing MPC Extended Format (array)..." << std::endl;
            for (const auto& item : j) {
                AsteroidProperties props = parseMPCExtendedJson(item);
                if (props.number > 0) {
                    asteroids_[props.number] = props;
                    count++;
                }
                
                // Progress every 10k
                if (count % 10000 == 0) {
                    std::cout << "  Loaded " << count << " asteroids..." << std::endl;
                }
            }
        }
        // Legacy format: {"asteroids": {...}}
        else if (j.contains("asteroids")) {
            std::cout << "Parsing legacy format (object)..." << std::endl;
            for (auto& [key, value] : j["asteroids"].items()) {
                int number = std::stoi(key);
                AsteroidProperties props = parseMPCExtendedJson(value);
                props.number = number;
                asteroids_[number] = props;
                count++;
            }
        }
        else {
            std::cerr << "Unknown JSON format" << std::endl;
            return 0;
        }
        
        std::cout << "✓ Loaded " << count << " asteroids from MPC data" << std::endl;
        
        updateStats();
        return count;
        
    } catch (const std::exception& e) {
        std::cerr << "Error importing: " << e.what() << std::endl;
        return 0;
    }
}

// Stub implementations for download functions
bool AsteroidDatabase::downloadDatabase(DataSource source, int maxAsteroids) {
    std::cerr << "Download not yet implemented. Use external tool to download MPC data.\n";
    return false;
}

bool AsteroidDatabase::updateAsteroid(int number, DataSource source) {
    std::cerr << "Update not yet implemented.\n";
    return false;
}

int AsteroidDatabase::updateRange(int from, int to, DataSource source) {
    std::cerr << "Update range not yet implemented.\n";
    return 0;
}

} // namespace ioccultcalc
