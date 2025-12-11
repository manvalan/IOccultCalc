/**
 * @file config_manager.h
 * @brief Configuration management system for IOccultCalc
 * 
 * Inspired by OrbFit's .oop file format, this system provides a flexible
 * configuration mechanism with support for JSON import/export.
 */

#ifndef IOCCULTCALC_CONFIG_MANAGER_H
#define IOCCULTCALC_CONFIG_MANAGER_H

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <optional>
#include <nlohmann/json.hpp>
#ifdef IOC_CONFIG_AVAILABLE
#include "ioc_config/oop_parser.h"
namespace ioc_config {
    class OopParser;
}
#endif

namespace ioccultcalc {

/**
 * @brief Configuration section types (inspired by OrbFit .oop sections)
 */
enum class ConfigSection {
    OBJECT,        // Target asteroid configuration
    PROPAGATION,   // Orbit propagation settings
    EPHEMERIS,     // JPL ephemerides configuration
    OUTPUT,        // Output format and options
    ERROR_MODEL,   // Error model settings
    OPERATIONS,    // Operations to perform
    PERTURBATIONS, // Perturbation model settings
    SEARCH,        // Occultation search parameters
    STAR,          // Star catalog configuration
    IERS,          // Earth orientation parameters
    OBSERVER,      // Observer location constraints
    FILTERING,     // Quality and observability filters
    SCORING,       // Priority scoring system
    PERFORMANCE,   // Performance and optimization
    DATABASE,      // Asteroid database filters
    GAIA,          // Gaia catalog settings
    ASTDYS,        // AstDyS local data directories
    ORBIT_FITTING, // Orbit fitting with observations
    VALIDATION,    // Validation and quality control
    CHEBYSHEV,     // Chebyshev approximation settings (LinOccult method)
    ASTEROIDS,     // Asteroid list files and selection
    CUSTOM         // User-defined sections
};

/**
 * @brief Single configuration parameter
 */
struct ConfigParameter {
    std::string name;
    std::string value;
    std::string type;  // "string", "double", "int", "bool", "array"
    std::string comment;
    
    // Type-safe accessors
    std::string asString() const;
    double asDouble() const;
    int asInt() const;
    bool asBool() const;
    std::vector<std::string> asArray() const;
    
    nlohmann::json toJson() const;
    static ConfigParameter fromJson(const nlohmann::json& j);
};

/**
 * @brief Configuration section containing parameters
 */
class ConfigSectionData {
public:
    ConfigSectionData() : type_(ConfigSection::CUSTOM), name_("") {}
    ConfigSectionData(ConfigSection type, const std::string& name = "");
    
    // Parameter management
    void setParameter(const std::string& name, const std::string& value, 
                     const std::string& comment = "");
    void setParameter(const std::string& name, double value, 
                     const std::string& comment = "");
    void setParameter(const std::string& name, int value, 
                     const std::string& comment = "");
    void setParameter(const std::string& name, bool value, 
                     const std::string& comment = "");
    void setParameter(const std::string& name, const std::vector<std::string>& value, 
                     const std::string& comment = "");
    
    std::optional<ConfigParameter> getParameter(const std::string& name) const;
    bool hasParameter(const std::string& name) const;
    void removeParameter(const std::string& name);
    
    // Accessors
    ConfigSection getType() const { return type_; }
    std::string getName() const { return name_; }
    std::map<std::string, ConfigParameter> getAllParameters() const { return parameters_; }
    
    // Serialization
    nlohmann::json toJson() const;
    static ConfigSectionData fromJson(const nlohmann::json& j);
    std::string toOopFormat() const;
    
    // Helper methods (public for ConfigManager access)
    std::string sectionTypeToString() const;
    static ConfigSection stringToSectionType(const std::string& str);
    
private:
    ConfigSection type_;
    std::string name_;
    std::map<std::string, ConfigParameter> parameters_;
};

/**
 * @brief Main configuration manager (similar to OrbFit .oop file)
 */
class ConfigManager {
public:
    ConfigManager();
    ConfigManager(const ConfigManager& other);
    ConfigManager& operator=(const ConfigManager& other);
    
    // Section management
    void addSection(const ConfigSectionData& section);
    std::optional<ConfigSectionData> getSection(ConfigSection type) const;
    std::optional<ConfigSectionData> getSectionByName(const std::string& name) const;
    std::vector<ConfigSectionData> getAllSections() const;
    void removeSection(ConfigSection type);
    
    // Quick parameter access (across all sections)
    std::optional<ConfigParameter> findParameter(const std::string& name) const;
    
    // File I/O
    void loadFromJson(const std::string& filepath);
    void saveToJson(const std::string& filepath) const;
    void loadFromOop(const std::string& filepath);
    void saveToOop(const std::string& filepath) const;
    
    // JSON serialization
    nlohmann::json toJson() const;
    static ConfigManager fromJson(const nlohmann::json& j);
    
    // Preset configurations
    static ConfigManager createDefault();
    static ConfigManager createHighPrecision();
    static ConfigManager createFastSearch();
    
    // Validation
    bool validate(std::vector<std::string>& errors) const;
    
    // Metadata
    void setMetadata(const std::string& key, const std::string& value);
    std::optional<std::string> getMetadata(const std::string& key) const;
    
private:
    std::map<ConfigSection, ConfigSectionData> sections_;
    std::map<std::string, std::string> metadata_;
    
#ifdef IOC_CONFIG_AVAILABLE
    mutable std::unique_ptr<ioc_config::OopParser> parser_;
    void syncFromIOCConfig();
    void syncToIOCConfig();
#endif
    
    void parseOopLine(const std::string& line, ConfigSection& currentSection, 
                      std::string& currentSectionName);
};

/**
 * @brief Builder for creating configurations programmatically
 */
class ConfigBuilder {
public:
    ConfigBuilder();
    
    // Object configuration
    ConfigBuilder& setObject(const std::string& name, const std::string& id);
    ConfigBuilder& setObjectElements(const std::string& elementFile);
    
    // Propagation settings
    ConfigBuilder& setPropagator(const std::string& type); // "RK4", "RA15", "ORBFIT"
    ConfigBuilder& setStepSize(double days);
    ConfigBuilder& setTimeSpan(double startJD, double endJD);
    
    // Perturbations
    ConfigBuilder& enablePlanets(bool enable = true);
    ConfigBuilder& enableAsteroids(int count = 17); // 0 = none, 17 = AST17, etc.
    ConfigBuilder& enableRelativity(bool enable = true);
    
    // Ephemerides
    ConfigBuilder& setJplEphemeris(const std::string& version); // "DE441", "DE440"
    ConfigBuilder& setAst17File(const std::string& filepath);
    
    // Search parameters
    ConfigBuilder& setSearchRegion(double raMin, double raMax, 
                                   double decMin, double decMax);
    ConfigBuilder& setMagnitudeLimit(double magLimit);
    ConfigBuilder& setSearchInterval(double startJD, double endJD, double stepDays);
    
    // Output settings
    ConfigBuilder& setOutputFormat(const std::string& format); // "JSON", "KML", "TEXT"
    ConfigBuilder& setOutputFile(const std::string& filepath);
    ConfigBuilder& setVerbosity(int level); // 0 = quiet, 1 = normal, 2 = verbose
    
    // Operations
    ConfigBuilder& enableDownloadFromAstDyS(bool enable = true);
    ConfigBuilder& enablePropagation(bool enable = true);
    ConfigBuilder& enableOccultationSearch(bool enable = true);
    
    // Build final configuration
    ConfigManager build() const;
    
private:
    ConfigManager config_;
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_CONFIG_MANAGER_H
