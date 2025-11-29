// Test UnifiedGaiaCatalog API for Sierks star query
// Target: Gaia DR3 3411546266140512128 / UCAC4 552-011427
// RA=73.416106°, Dec=+20.331662°, G=12.131

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "ioc_gaialib/unified_gaia_catalog.h"
#include "ioc_gaialib/types.h"

int main() {
    std::cout << "=== Test UnifiedGaiaCatalog API ===\n\n";
    
    // Expected star data for (17030) Sierks occultation
    const uint64_t TARGET_SOURCE_ID = 3411546266140512128ULL;
    const double TARGET_RA = 73.416106;
    const double TARGET_DEC = 20.331662;
    const double TARGET_MAG = 12.131;
    
    // Get catalog path
    std::string home = std::getenv("HOME") ? std::getenv("HOME") : "/tmp";
    std::string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
    
    // Check environment override
    const char* env_path = std::getenv("IOC_GAIA_CATALOG_PATH");
    if (env_path) {
        catalog_path = env_path;
    }
    
    std::cout << "Catalog path: " << catalog_path << "\n\n";
    
    // Initialize catalog
    std::string json_config = R"({
        "catalog_type": "multifile_v2",
        "multifile_directory": ")" + catalog_path + R"(",
        "max_cached_chunks": 100,
        "log_level": "info"
    })";
    
    std::cout << "Initializing UnifiedGaiaCatalog...\n";
    if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) {
        std::cerr << "ERROR: Failed to initialize catalog\n";
        return 1;
    }
    
    auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    
    // Get catalog info
    auto info = catalog.getCatalogInfo();
    std::cout << "Catalog: " << info.catalog_name << "\n";
    std::cout << "Total stars: " << info.total_stars << "\n";
    std::cout << "Magnitude limit: " << info.magnitude_limit << "\n\n";
    
    // Test 1: Query by source_id
    std::cout << "=== Test 1: Query by source_id ===\n";
    std::cout << "Looking for Gaia DR3 " << TARGET_SOURCE_ID << "...\n";
    
    auto star_opt = catalog.queryBySourceId(TARGET_SOURCE_ID);
    if (star_opt) {
        auto& star = *star_opt;
        std::cout << "FOUND!\n";
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  RA:       " << star.ra << "° (expected " << TARGET_RA << "°)\n";
        std::cout << "  Dec:      " << star.dec << "° (expected " << TARGET_DEC << "°)\n";
        std::cout << std::setprecision(3);
        std::cout << "  G mag:    " << star.phot_g_mean_mag << " (expected " << TARGET_MAG << ")\n";
        std::cout << "  BP mag:   " << star.phot_bp_mean_mag << "\n";
        std::cout << "  RP mag:   " << star.phot_rp_mean_mag << "\n";
        std::cout << "  Parallax: " << star.parallax << " ± " << star.parallax_error << " mas\n";
        std::cout << "  PM RA:    " << star.pmra << " ± " << star.pmra_error << " mas/yr\n";
        std::cout << "  PM Dec:   " << star.pmdec << " ± " << star.pmdec_error << " mas/yr\n";
        std::cout << "  Desig:    " << star.getDesignation() << "\n";
    } else {
        std::cout << "NOT FOUND by source_id\n";
    }
    
    // Test 2: Cone search around target position
    std::cout << "\n=== Test 2: Cone search ===\n";
    std::cout << "Searching RA=" << TARGET_RA << "°, Dec=" << TARGET_DEC << "°, r=0.1°, mag<14\n";
    
    ioc::gaia::QueryParams params;
    params.ra_center = TARGET_RA;
    params.dec_center = TARGET_DEC;
    params.radius = 0.1;
    params.max_magnitude = 14.0;
    params.min_parallax = -1.0;
    
    auto stars = catalog.queryCone(params);
    std::cout << "Found " << stars.size() << " stars\n";
    
    // Find target in results
    bool found_target = false;
    for (const auto& star : stars) {
        if (star.source_id == static_cast<int64_t>(TARGET_SOURCE_ID)) {
            found_target = true;
            std::cout << "\nTarget star found in cone search:\n";
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "  RA:  " << star.ra << "°\n";
            std::cout << "  Dec: " << star.dec << "°\n";
            std::cout << std::setprecision(3);
            std::cout << "  G:   " << star.phot_g_mean_mag << " mag\n";
            break;
        }
    }
    
    if (!found_target && stars.size() > 0) {
        std::cout << "\nFirst 5 stars found:\n";
        for (size_t i = 0; i < std::min(size_t(5), stars.size()); i++) {
            const auto& star = stars[i];
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "  [" << i+1 << "] RA=" << star.ra << "°, Dec=" << star.dec;
            std::cout << std::setprecision(3);
            std::cout << ", G=" << star.phot_g_mean_mag << ", ID=" << star.source_id << "\n";
        }
    }
    
    // Test 3: Corridor search (for occultation path)
    std::cout << "\n=== Test 3: Corridor search ===\n";
    
    ioc::gaia::CorridorQueryParams corridor;
    corridor.path.push_back({73.0, 20.0});
    corridor.path.push_back({74.0, 21.0});
    corridor.width = 0.5;
    corridor.max_magnitude = 14.0;
    
    std::cout << "Corridor: (73°,20°) → (74°,21°), width=0.5°, mag<14\n";
    
    auto corridor_stars = catalog.queryCorridor(corridor);
    std::cout << "Found " << corridor_stars.size() << " stars in corridor\n";
    
    // Check if target is in corridor
    for (const auto& star : corridor_stars) {
        if (star.source_id == static_cast<int64_t>(TARGET_SOURCE_ID)) {
            std::cout << "Target star FOUND in corridor!\n";
            break;
        }
    }
    
    // Statistics
    std::cout << "\n=== Statistics ===\n";
    auto stats = catalog.getStatistics();
    std::cout << "Total queries: " << stats.total_queries << "\n";
    std::cout << "Avg query time: " << stats.average_query_time_ms << " ms\n";
    std::cout << "Total stars returned: " << stats.total_stars_returned << "\n";
    std::cout << "Cache hit rate: " << (stats.cache_hit_rate * 100) << "%\n";
    
    // Cleanup
    ioc::gaia::UnifiedGaiaCatalog::shutdown();
    
    std::cout << "\n=== Test Complete ===\n";
    return 0;
}
