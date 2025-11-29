// Adapter layer: ioccultcalc API → IOC_GaiaLib API (v2.0)
// Uses new UnifiedGaiaCatalog API
#include "ioccultcalc/gaia_cache.h"
#include "ioccultcalc/gaia_client.h"
#include "ioc_gaialib/unified_gaia_catalog.h"
#include "ioc_gaialib/types.h"
#include <memory>
#include <map>
#include <mutex>
#include <iostream>
#include <cstdlib>

namespace ioccultcalc {

// Global mutex for thread safety
static std::mutex g_mutex;
static bool g_catalog_initialized = false;

// Dummy Impl to satisfy header declarations (not used in adapter)
class GaiaCache::Impl {};
class GaiaClient::Impl {};

// Convert ioc::gaia::GaiaStar → ioccultcalc::GaiaStar
static GaiaStar convertStar(const ioc::gaia::GaiaStar& src) {
    GaiaStar dst;
    dst.sourceId = std::to_string(src.source_id);
    dst.pos.ra = src.ra;
    dst.pos.dec = src.dec;
    dst.parallax = src.parallax;
    dst.pmra = src.pmra;
    dst.pmdec = src.pmdec;
    dst.phot_g_mean_mag = src.phot_g_mean_mag;
    dst.phot_bp_mean_mag = src.phot_bp_mean_mag;
    dst.phot_rp_mean_mag = src.phot_rp_mean_mag;
    return dst;
}

// Initialize UnifiedGaiaCatalog (called once)
static bool initializeCatalog() {
    if (g_catalog_initialized) return true;
    
    // Get catalog path from environment or use default
    std::string home = std::getenv("HOME") ? std::getenv("HOME") : "/tmp";
    std::string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
    
    // Check environment override
    const char* env_path = std::getenv("IOC_GAIA_CATALOG_PATH");
    if (env_path) {
        catalog_path = env_path;
    }
    
    // Build JSON config
    std::string json_config = R"({
        "catalog_type": "multifile_v2",
        "multifile_directory": ")" + catalog_path + R"(",
        "max_cached_chunks": 100,
        "log_level": "info"
    })";
    
    std::cout << "[INFO] Initializing UnifiedGaiaCatalog..." << std::endl;
    std::cout << "[INFO] Catalog path: " << catalog_path << std::endl;
    
    if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) {
        std::cerr << "[ERROR] Failed to initialize UnifiedGaiaCatalog" << std::endl;
        return false;
    }
    
    auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    auto info = catalog.getCatalogInfo();
    std::cout << "[INFO] Catalog initialized: " << info.catalog_name 
              << " (" << info.total_stars << " stars, mag≤" << info.magnitude_limit << ")" << std::endl;
    
    g_catalog_initialized = true;
    return true;
}

// GaiaCache implementation
GaiaCache::GaiaCache(const std::string& cacheDir) {
    std::lock_guard<std::mutex> lock(g_mutex);
    initializeCatalog();
}

GaiaCache::~GaiaCache() {
    // Singleton pattern - no cleanup needed
}

void GaiaCache::setGaiaVersion(GaiaVersion version) {
    // UnifiedGaiaCatalog is always DR3
    std::cout << "[INFO] Gaia version set (using DR3 catalog)" << std::endl;
}

bool GaiaCache::loadIndex() {
    std::lock_guard<std::mutex> lock(g_mutex);
    return initializeCatalog();
}

std::vector<GaiaStar> GaiaCache::queryRegion(
    double raDeg, double decDeg, double radiusDeg, double maxMag, bool useCache) {
    
    std::lock_guard<std::mutex> lock(g_mutex);
    
    if (!g_catalog_initialized) {
        if (!initializeCatalog()) {
            std::cerr << "[ERROR] Catalog not initialized" << std::endl;
            return {};
        }
    }
    
    try {
        auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
        
        // Build query parameters
        ioc::gaia::QueryParams params;
        params.ra_center = raDeg;
        params.dec_center = decDeg;
        params.radius = radiusDeg;
        params.max_magnitude = maxMag;
        params.min_parallax = -1.0;  // No parallax filter
        
        // Execute query
        auto ioStars = catalog.queryCone(params);
        
        // Log periodically
        static size_t query_count = 0;
        if (++query_count % 100 == 0) {
            auto stats = catalog.getStatistics();
            std::cout << "[INFO] UnifiedGaiaCatalog stats: " 
                      << stats.total_queries << " queries, "
                      << stats.total_stars_returned << " stars returned, "
                      << "cache hit rate: " << (stats.cache_hit_rate * 100) << "%" << std::endl;
        }
        
        // Convert to ioccultcalc format
        std::vector<GaiaStar> result;
        result.reserve(ioStars.size());
        for (const auto& star : ioStars) {
            result.push_back(convertStar(star));
        }
        
        return result;
        
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Catalog query failed: " << e.what() << std::endl;
        return {};
    }
}

GaiaCacheStats GaiaCache::getStats() const {
    std::lock_guard<std::mutex> lock(g_mutex);
    GaiaCacheStats stats;
    
    if (g_catalog_initialized) {
        try {
            auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
            auto catStats = catalog.getStatistics();
            auto info = catalog.getCatalogInfo();
            
            stats.total_tiles = 0;  // N/A for unified catalog
            stats.total_stars = info.total_stars;
            stats.total_size_mb = catStats.memory_used_mb + catStats.disk_cache_used_mb;
            stats.sky_coverage = 100.0;  // Full sky coverage
        } catch (...) {}
    }
    
    return stats;
}

// ==================== GaiaClient Adapter ====================

GaiaClient::GaiaClient() {
    std::lock_guard<std::mutex> lock(g_mutex);
    initializeCatalog();
}

GaiaClient::~GaiaClient() {
    // Singleton pattern - no cleanup needed
}

std::vector<GaiaStar> GaiaClient::queryRegion(
    double raCenterDeg, double decCenterDeg, double radiusDeg, double maxMagnitude) {
    
    std::lock_guard<std::mutex> lock(g_mutex);
    
    if (!g_catalog_initialized) {
        if (!initializeCatalog()) {
            return {};
        }
    }
    
    try {
        auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
        
        std::cout << "[DEBUG] queryCone RA=" << raCenterDeg 
                  << " Dec=" << decCenterDeg << " r=" << radiusDeg << "°..." << std::flush;
        
        ioc::gaia::QueryParams params;
        params.ra_center = raCenterDeg;
        params.dec_center = decCenterDeg;
        params.radius = radiusDeg;
        params.max_magnitude = maxMagnitude;
        params.min_parallax = -1.0;
        
        auto ioStars = catalog.queryCone(params);
        std::cout << " OK (" << ioStars.size() << " stelle)" << std::endl;
        
        std::vector<GaiaStar> result;
        result.reserve(ioStars.size());
        for (const auto& star : ioStars) {
            result.push_back(convertStar(star));
        }
        return result;
        
    } catch (const std::exception& e) {
        std::cerr << " ERROR: " << e.what() << std::endl;
        return {};
    }
}

std::vector<GaiaStar> GaiaClient::queryCone(
    double raCenterDeg, double decCenterDeg, double radiusDeg, double maxMagnitude) {
    return queryRegion(raCenterDeg, decCenterDeg, radiusDeg, maxMagnitude);
}

void GaiaClient::setGaiaVersion(GaiaVersion version) {
    // UnifiedGaiaCatalog is always DR3
    std::cout << "[INFO] Gaia version set (using DR3 catalog)" << std::endl;
}

// ==================== queryAlongPath using Corridor API ====================

std::vector<GaiaStar> GaiaClient::queryAlongPath(
    const std::vector<EquatorialCoordinates>& path,
    double widthDeg, double maxMagnitude) {
    
    std::lock_guard<std::mutex> lock(g_mutex);
    
    if (!g_catalog_initialized) {
        if (!initializeCatalog()) {
            return {};
        }
    }
    
    if (path.size() < 2) {
        std::cerr << "[ERROR] queryAlongPath requires at least 2 points" << std::endl;
        return {};
    }
    
    try {
        auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
        
        // Build corridor query params
        ioc::gaia::CorridorQueryParams corridor;
        for (const auto& pt : path) {
            corridor.path.push_back({pt.ra, pt.dec});
        }
        corridor.width = widthDeg;
        corridor.max_magnitude = maxMagnitude;
        corridor.min_parallax = -1.0;  // No parallax filter
        corridor.max_results = 0;      // No limit
        
        std::cout << "[DEBUG] queryCorridor: " << path.size() << " points, width=" 
                  << widthDeg << "°, mag<" << maxMagnitude << "..." << std::flush;
        
        auto ioStars = catalog.queryCorridor(corridor);
        std::cout << " OK (" << ioStars.size() << " stelle)" << std::endl;
        
        // Convert to ioccultcalc format
        std::vector<GaiaStar> result;
        result.reserve(ioStars.size());
        for (const auto& star : ioStars) {
            result.push_back(convertStar(star));
        }
        return result;
        
    } catch (const std::exception& e) {
        std::cerr << " ERROR: " << e.what() << std::endl;
        return {};
    }
}

// ==================== GaiaStar Methods ====================

EquatorialCoordinates GaiaStar::propagateToEpoch(const JulianDate& epoch) const {
    // Propagate from J2016.0 to target epoch using proper motion
    double dt_years = (epoch.jd - 2457389.0) / 365.25;  // Years since J2016.0
    
    EquatorialCoordinates result;
    result.ra = pos.ra + (pmra / 3600000.0) * dt_years / std::cos(pos.dec * M_PI / 180.0);
    result.dec = pos.dec + (pmdec / 3600000.0) * dt_years;
    
    return result;
}

} // namespace ioccultcalc
