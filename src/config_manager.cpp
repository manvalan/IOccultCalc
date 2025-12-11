/**
 * @file config_manager.cpp
 * @brief Implementation of configuration management system
 * 
 * Uses IOC_Config library internally when available
 */

#include "ioccultcalc/config_manager.h"
#ifdef IOC_CONFIG_AVAILABLE
#include "ioc_config/oop_parser.h"
#endif
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace ioccultcalc {

// ============================================================================
// ConfigParameter Implementation
// ============================================================================

std::string ConfigParameter::asString() const {
    return value;
}

double ConfigParameter::asDouble() const {
    try {
        return std::stod(value);
    } catch (const std::exception& e) {
        throw std::runtime_error("Cannot convert '" + value + "' to double: " + e.what());
    }
}

int ConfigParameter::asInt() const {
    try {
        return std::stoi(value);
    } catch (const std::exception& e) {
        throw std::runtime_error("Cannot convert '" + value + "' to int: " + e.what());
    }
}

bool ConfigParameter::asBool() const {
    std::string lower = value;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "true" || lower == "1" || lower == "t" || lower == ".true.") {
        return true;
    } else if (lower == "false" || lower == "0" || lower == "f" || lower == ".false.") {
        return false;
    }
    throw std::runtime_error("Cannot convert '" + value + "' to bool");
}

std::vector<std::string> ConfigParameter::asArray() const {
    std::vector<std::string> result;
    std::stringstream ss(value);
    std::string item;
    
    while (std::getline(ss, item, ',')) {
        // Trim whitespace
        item.erase(0, item.find_first_not_of(" \t"));
        item.erase(item.find_last_not_of(" \t") + 1);
        if (!item.empty()) {
            result.push_back(item);
        }
    }
    
    return result;
}

nlohmann::json ConfigParameter::toJson() const {
    nlohmann::json j;
    j["name"] = name;
    j["value"] = value;
    j["type"] = type;
    if (!comment.empty()) {
        j["comment"] = comment;
    }
    return j;
}

ConfigParameter ConfigParameter::fromJson(const nlohmann::json& j) {
    ConfigParameter param;
    param.name = j["name"];
    param.value = j["value"];
    param.type = j["type"];
    if (j.contains("comment")) {
        param.comment = j["comment"];
    }
    return param;
}

// ============================================================================
// ConfigSectionData Implementation
// ============================================================================

ConfigSectionData::ConfigSectionData(ConfigSection type, const std::string& name)
    : type_(type), name_(name) {
}

void ConfigSectionData::setParameter(const std::string& name, const std::string& value, 
                                      const std::string& comment) {
    ConfigParameter param;
    param.name = name;
    param.value = value;
    param.type = "string";
    param.comment = comment;
    parameters_[name] = param;
}

void ConfigSectionData::setParameter(const std::string& name, double value, 
                                      const std::string& comment) {
    ConfigParameter param;
    param.name = name;
    param.value = std::to_string(value);
    param.type = "double";
    param.comment = comment;
    parameters_[name] = param;
}

void ConfigSectionData::setParameter(const std::string& name, int value, 
                                      const std::string& comment) {
    ConfigParameter param;
    param.name = name;
    param.value = std::to_string(value);
    param.type = "int";
    param.comment = comment;
    parameters_[name] = param;
}

void ConfigSectionData::setParameter(const std::string& name, bool value, 
                                      const std::string& comment) {
    ConfigParameter param;
    param.name = name;
    param.value = value ? "true" : "false";
    param.type = "bool";
    param.comment = comment;
    parameters_[name] = param;
}

void ConfigSectionData::setParameter(const std::string& name, 
                                      const std::vector<std::string>& value, 
                                      const std::string& comment) {
    ConfigParameter param;
    param.name = name;
    std::stringstream ss;
    for (size_t i = 0; i < value.size(); ++i) {
        if (i > 0) ss << ",";
        ss << value[i];
    }
    param.value = ss.str();
    param.type = "array";
    param.comment = comment;
    parameters_[name] = param;
}

std::optional<ConfigParameter> ConfigSectionData::getParameter(const std::string& name) const {
    auto it = parameters_.find(name);
    if (it != parameters_.end()) {
        return it->second;
    }
    return std::nullopt;
}

bool ConfigSectionData::hasParameter(const std::string& name) const {
    return parameters_.find(name) != parameters_.end();
}

void ConfigSectionData::removeParameter(const std::string& name) {
    parameters_.erase(name);
}

nlohmann::json ConfigSectionData::toJson() const {
    nlohmann::json j;
    j["type"] = sectionTypeToString();
    if (!name_.empty()) {
        j["name"] = name_;
    }
    
    nlohmann::json params = nlohmann::json::array();
    for (const auto& [key, param] : parameters_) {
        params.push_back(param.toJson());
    }
    j["parameters"] = params;
    
    return j;
}

ConfigSectionData ConfigSectionData::fromJson(const nlohmann::json& j) {
    ConfigSection type = stringToSectionType(j["type"]);
    std::string name = j.contains("name") ? j["name"].get<std::string>() : "";
    
    ConfigSectionData section(type, name);
    
    if (j.contains("parameters")) {
        for (const auto& paramJson : j["parameters"]) {
            ConfigParameter param = ConfigParameter::fromJson(paramJson);
            section.parameters_[param.name] = param;
        }
    }
    
    return section;
}

std::string ConfigSectionData::toOopFormat() const {
    std::stringstream ss;
    
    // Section header
    ss << sectionTypeToString();
    if (!name_.empty()) {
        ss << " '" << name_ << "'";
    }
    ss << ".\n";
    
    // Parameters
    for (const auto& [key, param] : parameters_) {
        ss << "        ." << param.name << " = ";
        
        if (param.type == "string") {
            ss << "'" << param.value << "'";
        } else if (param.type == "bool") {
            ss << (param.asBool() ? ".TRUE." : ".FALSE.");
        } else {
            ss << param.value;
        }
        
        if (!param.comment.empty()) {
            ss << "  ! " << param.comment;
        }
        ss << "\n";
    }
    
    ss << "\n";
    return ss.str();
}

std::string ConfigSectionData::sectionTypeToString() const {
    switch (type_) {
        case ConfigSection::OBJECT: return "object";
        case ConfigSection::PROPAGATION: return "propag";
        case ConfigSection::EPHEMERIS: return "ephemeris";
        case ConfigSection::OUTPUT: return "output";
        case ConfigSection::ERROR_MODEL: return "error_model";
        case ConfigSection::OPERATIONS: return "operations";
        case ConfigSection::PERTURBATIONS: return "perturbations";
        case ConfigSection::SEARCH: return "search";
        case ConfigSection::STAR: return "star";
        case ConfigSection::IERS: return "IERS";
        case ConfigSection::OBSERVER: return "observer";
        case ConfigSection::FILTERING: return "filtering";
        case ConfigSection::SCORING: return "scoring";
        case ConfigSection::PERFORMANCE: return "performance";
        case ConfigSection::DATABASE: return "database";
        case ConfigSection::GAIA: return "gaia";
        case ConfigSection::ASTDYS: return "astdys";
        case ConfigSection::ORBIT_FITTING: return "orbit_fitting";
        case ConfigSection::VALIDATION: return "validation";
        case ConfigSection::CHEBYSHEV: return "chebyshev";
        case ConfigSection::CUSTOM: return "custom";
        default: return "unknown";
    }
}

ConfigSection ConfigSectionData::stringToSectionType(const std::string& str) {
    if (str == "object") return ConfigSection::OBJECT;
    if (str == "propag" || str == "propagation") return ConfigSection::PROPAGATION;
    if (str == "ephemeris") return ConfigSection::EPHEMERIS;
    if (str == "output") return ConfigSection::OUTPUT;
    if (str == "error_model") return ConfigSection::ERROR_MODEL;
    if (str == "operations") return ConfigSection::OPERATIONS;
    if (str == "perturbations") return ConfigSection::PERTURBATIONS;
    if (str == "search") return ConfigSection::SEARCH;
    if (str == "star") return ConfigSection::STAR;
    if (str == "IERS" || str == "iers") return ConfigSection::IERS;
    if (str == "observer") return ConfigSection::OBSERVER;
    if (str == "filtering") return ConfigSection::FILTERING;
    if (str == "scoring") return ConfigSection::SCORING;
    if (str == "performance") return ConfigSection::PERFORMANCE;
    if (str == "database") return ConfigSection::DATABASE;
    if (str == "gaia") return ConfigSection::GAIA;
    if (str == "astdys") return ConfigSection::ASTDYS;
    if (str == "orbit_fitting") return ConfigSection::ORBIT_FITTING;
    if (str == "validation") return ConfigSection::VALIDATION;
    if (str == "chebyshev") return ConfigSection::CHEBYSHEV;
    if (str == "asteroids") return ConfigSection::ASTEROIDS;
    return ConfigSection::CUSTOM;
}

// ============================================================================
// ConfigManager Implementation
// ============================================================================

ConfigManager::ConfigManager() {
// Temporarily disabled IOC_Config integration due to compilation issues
// #ifdef IOC_CONFIG_AVAILABLE
//     parser_ = std::make_unique<ioc_config::OopParser>();
// #endif
}

ConfigManager::ConfigManager(const ConfigManager& other)
    : sections_(other.sections_), metadata_(other.metadata_) {
    // parser_ is not copied - will be created on demand if needed
}

ConfigManager& ConfigManager::operator=(const ConfigManager& other) {
    if (this != &other) {
        sections_ = other.sections_;
        metadata_ = other.metadata_;
        // parser_ is not copied - will be created on demand if needed
    }
    return *this;
}

void ConfigManager::addSection(const ConfigSectionData& section) {
    sections_[section.getType()] = section;
}

std::optional<ConfigSectionData> ConfigManager::getSection(ConfigSection type) const {
    auto it = sections_.find(type);
    if (it != sections_.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<ConfigSectionData> ConfigManager::getSectionByName(const std::string& name) const {
    for (const auto& [type, section] : sections_) {
        if (section.getName() == name) {
            return section;
        }
    }
    return std::nullopt;
}

std::vector<ConfigSectionData> ConfigManager::getAllSections() const {
    std::vector<ConfigSectionData> result;
    for (const auto& [type, section] : sections_) {
        result.push_back(section);
    }
    return result;
}

void ConfigManager::removeSection(ConfigSection type) {
    sections_.erase(type);
}

std::optional<ConfigParameter> ConfigManager::findParameter(const std::string& name) const {
    for (const auto& [type, section] : sections_) {
        auto param = section.getParameter(name);
        if (param) {
            return param;
        }
    }
    return std::nullopt;
}

void ConfigManager::loadFromJson(const std::string& filepath) {
// Temporarily using native JSON parser
// #ifdef IOC_CONFIG_AVAILABLE
//     if (!parser_) {
//         parser_ = std::make_unique<ioc_config::OopParser>();
//     }
//     
//     if (!parser_->loadFromJson(filepath)) {
//         throw std::runtime_error("Failed to load JSON: " + parser_->getLastError());
//     }
//     
//     syncFromIOCConfig();
// #else
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    nlohmann::json j;
    file >> j;
    
    *this = fromJson(j);
// #endif
}

void ConfigManager::saveToJson(const std::string& filepath) const {
    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot create file: " + filepath);
    }
    
    nlohmann::json j = toJson();
    file << j.dump(2);
}

void ConfigManager::loadFromOop(const std::string& filepath) {
// Temporarily using native OOP parser
// #ifdef IOC_CONFIG_AVAILABLE
//     if (!parser_) {
//         parser_ = std::make_unique<ioc_config::OopParser>();
//     }
//     
//     if (!parser_->loadFromOop(filepath)) {
//         throw std::runtime_error("Failed to load OOP: " + parser_->getLastError());
//     }
//     
//     syncFromIOCConfig();
// #else
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    ConfigSection currentSection = ConfigSection::CUSTOM;
    std::string currentSectionName;
    std::string line;
    
    while (std::getline(file, line)) {
        parseOopLine(line, currentSection, currentSectionName);
    }
// #endif
}

void ConfigManager::saveToOop(const std::string& filepath) const {
    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot create file: " + filepath);
    }
    
    file << "! IOccultCalc Configuration File\n";
    file << "! Generated automatically\n\n";
    
    // Write metadata
    for (const auto& [key, value] : metadata_) {
        file << "! " << key << ": " << value << "\n";
    }
    file << "\n";
    
    // Write sections
    for (const auto& [type, section] : sections_) {
        file << section.toOopFormat();
    }
}

nlohmann::json ConfigManager::toJson() const {
    nlohmann::json j;
    
    // Metadata
    if (!metadata_.empty()) {
        j["metadata"] = metadata_;
    }
    
    // Sections
    nlohmann::json sections = nlohmann::json::array();
    for (const auto& [type, section] : sections_) {
        sections.push_back(section.toJson());
    }
    j["sections"] = sections;
    
    return j;
}

ConfigManager ConfigManager::fromJson(const nlohmann::json& j) {
    ConfigManager config;
    
    // Load metadata
    if (j.contains("metadata")) {
        config.metadata_ = j["metadata"].get<std::map<std::string, std::string>>();
    }
    
    // Load sections
    if (j.contains("sections")) {
        for (const auto& sectionJson : j["sections"]) {
            ConfigSectionData section = ConfigSectionData::fromJson(sectionJson);
            config.addSection(section);
        }
    }
    
    return config;
}

ConfigManager ConfigManager::createDefault() {
    ConfigBuilder builder;
    return builder
        .setObject("Default", "0")
        .setPropagator("RK4")
        .setStepSize(0.05)
        .enablePlanets(true)
        .enableAsteroids(17)
        .setJplEphemeris("DE441")
        .setOutputFormat("JSON")
        .setVerbosity(1)
        .enableDownloadFromAstDyS(true)   // DEFAULT: scarica da AstDyS
        .enablePropagation(true)
        .enableOccultationSearch(true)
        .build();
}

ConfigManager ConfigManager::createHighPrecision() {
    ConfigBuilder builder;
    return builder
        .setObject("HighPrecision", "0")
        .setPropagator("RA15")
        .setStepSize(0.01)
        .enablePlanets(true)
        .enableAsteroids(17)
        .enableRelativity(true)
        .setJplEphemeris("DE441")
        .setOutputFormat("JSON")
        .setVerbosity(2)
        .enableDownloadFromAstDyS(true)   // DEFAULT: scarica da AstDyS
        .enablePropagation(true)
        .enableOccultationSearch(true)
        .build();
}

ConfigManager ConfigManager::createFastSearch() {
    ConfigBuilder builder;
    return builder
        .setObject("FastSearch", "0")
        .setPropagator("RK4")
        .setStepSize(0.1)
        .enablePlanets(true)
        .enableAsteroids(0)
        .setJplEphemeris("DE441")
        .setOutputFormat("JSON")
        .setVerbosity(0)
        .enableDownloadFromAstDyS(true)   // DEFAULT: scarica da AstDyS
        .enablePropagation(true)
        .enableOccultationSearch(true)
        .build();
}

bool ConfigManager::validate(std::vector<std::string>& errors) const {
    errors.clear();
    
    // Check required sections
    if (!getSection(ConfigSection::OBJECT)) {
        errors.push_back("Missing required section: object");
    }
    
    if (!getSection(ConfigSection::PROPAGATION)) {
        errors.push_back("Missing required section: propagation");
    }
    
    // Validate propagation parameters
    auto propag = getSection(ConfigSection::PROPAGATION);
    if (propag) {
        auto stepSize = propag->getParameter("step_size");
        if (stepSize && stepSize->asDouble() <= 0) {
            errors.push_back("Invalid step_size: must be positive");
        }
    }
    
    return errors.empty();
}

void ConfigManager::setMetadata(const std::string& key, const std::string& value) {
    metadata_[key] = value;
}

std::optional<std::string> ConfigManager::getMetadata(const std::string& key) const {
    auto it = metadata_.find(key);
    if (it != metadata_.end()) {
        return it->second;
    }
    return std::nullopt;
}

void ConfigManager::parseOopLine(const std::string& line, ConfigSection& currentSection, 
                                  std::string& currentSectionName) {
    // Skip comments and empty lines
    std::string trimmed = line;
    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
    if (trimmed.empty() || trimmed[0] == '!') {
        return;
    }
    
    // Section header detection (ends with '.')
    if (trimmed.back() == '.') {
        std::string sectionStr = trimmed.substr(0, trimmed.length() - 1);
        currentSection = ConfigSectionData::stringToSectionType(sectionStr);
        currentSectionName = sectionStr;
        
        ConfigSectionData newSection(currentSection, currentSectionName);
        addSection(newSection);
        return;
    }
    
    // Parameter parsing (format: .name = value ! comment)
    if (trimmed[0] == '.') {
        size_t eqPos = trimmed.find('=');
        if (eqPos == std::string::npos) return;
        
        std::string paramName = trimmed.substr(1, eqPos - 1);
        paramName.erase(paramName.find_last_not_of(" \t") + 1);
        
        std::string valueStr = trimmed.substr(eqPos + 1);
        valueStr.erase(0, valueStr.find_first_not_of(" \t"));
        
        // Extract comment
        std::string comment;
        size_t commentPos = valueStr.find('!');
        if (commentPos != std::string::npos) {
            comment = valueStr.substr(commentPos + 1);
            comment.erase(0, comment.find_first_not_of(" \t"));
            valueStr = valueStr.substr(0, commentPos);
        }
        
        // Remove trailing whitespace from value
        valueStr.erase(valueStr.find_last_not_of(" \t") + 1);
        
        // Remove quotes if present
        if (!valueStr.empty() && valueStr.front() == '\'' && valueStr.back() == '\'') {
            valueStr = valueStr.substr(1, valueStr.length() - 2);
        }
        
        // Add parameter to current section
        auto section = getSection(currentSection);
        if (section) {
            ConfigSectionData updatedSection = *section;
            updatedSection.setParameter(paramName, valueStr, comment);
            addSection(updatedSection);
        }
    }
}

// ============================================================================
// ConfigBuilder Implementation
// ============================================================================

ConfigBuilder::ConfigBuilder() {
}

ConfigBuilder& ConfigBuilder::setObject(const std::string& name, const std::string& id) {
    ConfigSectionData section(ConfigSection::OBJECT);
    section.setParameter("name", name, "Object name");
    section.setParameter("id", id, "Object ID");
    config_.addSection(section);
    return *this;
}

ConfigBuilder& ConfigBuilder::setObjectElements(const std::string& elementFile) {
    auto section = config_.getSection(ConfigSection::OBJECT);
    if (!section) {
        section = ConfigSectionData(ConfigSection::OBJECT);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("element_file", elementFile, "Orbital elements file");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::setPropagator(const std::string& type) {
    ConfigSectionData section(ConfigSection::PROPAGATION);
    section.setParameter("type", type, "Propagator type (RK4, RA15, ORBFIT)");
    config_.addSection(section);
    return *this;
}

ConfigBuilder& ConfigBuilder::setStepSize(double days) {
    auto section = config_.getSection(ConfigSection::PROPAGATION);
    if (!section) {
        section = ConfigSectionData(ConfigSection::PROPAGATION);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("step_size", days, "Integration step size (days)");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::setTimeSpan(double startJD, double endJD) {
    auto section = config_.getSection(ConfigSection::PROPAGATION);
    if (!section) {
        section = ConfigSectionData(ConfigSection::PROPAGATION);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("start_jd", startJD, "Start epoch (Julian Date)");
    updatedSection.setParameter("end_jd", endJD, "End epoch (Julian Date)");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::enablePlanets(bool enable) {
    ConfigSectionData section(ConfigSection::PERTURBATIONS);
    section.setParameter("planets", enable, "Include planetary perturbations");
    config_.addSection(section);
    return *this;
}

ConfigBuilder& ConfigBuilder::enableAsteroids(int count) {
    auto section = config_.getSection(ConfigSection::PERTURBATIONS);
    if (!section) {
        section = ConfigSectionData(ConfigSection::PERTURBATIONS);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("asteroid_count", count, "Number of massive asteroids (0-17)");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::enableRelativity(bool enable) {
    auto section = config_.getSection(ConfigSection::PERTURBATIONS);
    if (!section) {
        section = ConfigSectionData(ConfigSection::PERTURBATIONS);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("relativity", enable, "Include relativistic corrections");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::setJplEphemeris(const std::string& version) {
    ConfigSectionData section(ConfigSection::EPHEMERIS);
    section.setParameter("jpl_version", version, "JPL ephemeris version (DE441, DE440)");
    config_.addSection(section);
    return *this;
}

ConfigBuilder& ConfigBuilder::setAst17File(const std::string& filepath) {
    auto section = config_.getSection(ConfigSection::EPHEMERIS);
    if (!section) {
        section = ConfigSectionData(ConfigSection::EPHEMERIS);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("ast17_file", filepath, "AST17 SPK file path");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::setSearchRegion(double raMin, double raMax, 
                                              double decMin, double decMax) {
    ConfigSectionData section(ConfigSection::SEARCH);
    section.setParameter("ra_min", raMin, "Minimum right ascension (degrees)");
    section.setParameter("ra_max", raMax, "Maximum right ascension (degrees)");
    section.setParameter("dec_min", decMin, "Minimum declination (degrees)");
    section.setParameter("dec_max", decMax, "Maximum declination (degrees)");
    config_.addSection(section);
    return *this;
}

ConfigBuilder& ConfigBuilder::setMagnitudeLimit(double magLimit) {
    auto section = config_.getSection(ConfigSection::SEARCH);
    if (!section) {
        section = ConfigSectionData(ConfigSection::SEARCH);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("mag_limit", magLimit, "Magnitude limit for star search");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::setSearchInterval(double startJD, double endJD, double stepDays) {
    auto section = config_.getSection(ConfigSection::SEARCH);
    if (!section) {
        section = ConfigSectionData(ConfigSection::SEARCH);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("start_jd", startJD, "Search start epoch (JD)");
    updatedSection.setParameter("end_jd", endJD, "Search end epoch (JD)");
    updatedSection.setParameter("step_days", stepDays, "Search sampling interval (days)");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::setOutputFormat(const std::string& format) {
    ConfigSectionData section(ConfigSection::OUTPUT);
    section.setParameter("format", format, "Output format (JSON, KML, TEXT)");
    config_.addSection(section);
    return *this;
}

ConfigBuilder& ConfigBuilder::setOutputFile(const std::string& filepath) {
    auto section = config_.getSection(ConfigSection::OUTPUT);
    if (!section) {
        section = ConfigSectionData(ConfigSection::OUTPUT);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("file", filepath, "Output file path");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::setVerbosity(int level) {
    auto section = config_.getSection(ConfigSection::OUTPUT);
    if (!section) {
        section = ConfigSectionData(ConfigSection::OUTPUT);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("verbosity", level, "Verbosity level (0=quiet, 1=normal, 2=verbose)");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::enableDownloadFromAstDyS(bool enable) {
    auto section = config_.getSection(ConfigSection::OPERATIONS);
    if (!section) {
        section = ConfigSectionData(ConfigSection::OPERATIONS);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("download_elements", enable, "Download orbital elements from AstDyS");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::enablePropagation(bool enable) {
    auto section = config_.getSection(ConfigSection::OPERATIONS);
    if (!section) {
        section = ConfigSectionData(ConfigSection::OPERATIONS);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("propagate", enable, "Propagate orbit");
    config_.addSection(updatedSection);
    return *this;
}

ConfigBuilder& ConfigBuilder::enableOccultationSearch(bool enable) {
    auto section = config_.getSection(ConfigSection::OPERATIONS);
    if (!section) {
        section = ConfigSectionData(ConfigSection::OPERATIONS);
    }
    ConfigSectionData updatedSection = *section;
    updatedSection.setParameter("search_occultations", enable, "Search for occultations");
    config_.addSection(updatedSection);
    return *this;
}

ConfigManager ConfigBuilder::build() const {
    return config_;
}

#ifdef IOC_CONFIG_AVAILABLE
// ============================================================================
// IOC_Config Integration - Sync methods
// ============================================================================

void ConfigManager::syncFromIOCConfig() {
    if (!parser_) return;
    
    sections_.clear();
    
    // Get all sections from IOC_Config
    auto iocSections = parser_->getAllSections();
    for (const auto& iocSection : iocSections) {
        std::string sectionName = iocSection.name;
        
        // Convert section name to ConfigSection enum
        ConfigSection sectionType = ConfigSectionData::stringToSectionType(sectionName);
        
        // Create ConfigSectionData
        ConfigSectionData section(sectionType, sectionName);
        
        // Copy parameters
        for (const auto& [key, iocParam] : iocSection.parameters) {
            // Remove leading dot from key if present
            std::string cleanKey = key;
            if (!cleanKey.empty() && cleanKey[0] == '.') {
                cleanKey = cleanKey.substr(1);
            }
            section.setParameter(cleanKey, iocParam.value);
        }
        
        sections_[sectionType] = section;
    }
}

void ConfigManager::syncToIOCConfig() {
    if (!parser_) {
        parser_ = std::make_unique<ioc_config::OopParser>();
    }
    
    // Clear parser
    parser_->clear();
    
    // Add all sections to IOC_Config
    for (const auto& [type, section] : sections_) {
        std::string sectionName = section.getName();
        if (sectionName.empty()) {
            sectionName = section.sectionTypeToString();
        }
        
        // Add parameters to IOC_Config
        for (const auto& [key, param] : section.getAllParameters()) {
            parser_->setParameter(sectionName, key, param.value);
        }
    }
}
#endif

} // namespace ioccultcalc
