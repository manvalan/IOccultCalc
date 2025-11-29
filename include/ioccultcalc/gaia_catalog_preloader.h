/**
 * @file gaia_catalog_preloader.h
 * @brief Smart preloader for Gaia Mag18 catalog - Phase 2 optimization
 * 
 * Strategy:
 * 1. Pre-analyze asteroid trajectories (rough estimation)
 * 2. Identify HEALPix tiles needed (10-50 out of 12288)
 * 3. Load only necessary tiles into RAM (~500 MB - 2 GB)
 * 4. Query stars in-memory (instant lookups)
 * 
 * Performance target: 8-12 min → 4-6 min (2× speedup)
 * Memory usage: 2-4 GB RAM (vs 12 GB full preload)
 * 
 * @author Michele Bigi
 * @date 2025-11-27
 */

#ifndef IOCCULTCALC_GAIA_CATALOG_PRELOADER_H
#define IOCCULTCALC_GAIA_CATALOG_PRELOADER_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <memory>
#include "star_catalog.h"
#include "orbital_elements.h"

namespace ioccultcalc {

/**
 * @brief Sky region for trajectory coverage analysis
 */
struct SkyRegion {
    double ra_min, ra_max;      // RA bounds (degrees)
    double dec_min, dec_max;    // Dec bounds (degrees)
    double ra_center, dec_center; // Center (for HEALPix lookup)
};

/**
 * @brief HEALPix tile with cached stars
 */
struct HEALPixTile {
    uint64_t healpix_id;                // HEALPix index (NSIDE=64)
    std::vector<StarData> stars;        // All stars in tile
    double ra_center, dec_center;       // Tile center
    size_t star_count;                  // Number of stars
    bool loaded;                        // Load status
    
    HEALPixTile() : healpix_id(0), ra_center(0), dec_center(0), 
                    star_count(0), loaded(false) {}
};

/**
 * @brief Smart catalog preloader - Phase 2 optimization
 * 
 * Loads only necessary HEALPix tiles based on asteroid trajectories.
 * Much faster than on-demand queries, uses less memory than full preload.
 */
class GaiaCatalogPreloader {
public:
    /**
     * @brief Constructor
     * @param catalogPath Path to Mag18 catalog file
     * @param healpixNSIDE HEALPix NSIDE (default 64 for Mag18 v2)
     */
    explicit GaiaCatalogPreloader(const std::string& catalogPath, int healpixNSIDE = 64);
    
    ~GaiaCatalogPreloader();
    
    // ==================== Pre-Analysis Phase ====================
    
    /**
     * @brief Analyze asteroid trajectories and identify needed HEALPix tiles
     * @param asteroids List of asteroids to analyze
     * @param searchRadiusDeg Search radius around each position
     * @return Number of unique tiles identified
     */
    size_t analyzeTrajectories(const std::vector<OrbitalElements>& asteroids,
                               double searchRadiusDeg = 7.0);
    
    /**
     * @brief Get sky regions covered by asteroid trajectories
     * @param asteroids Asteroid orbital elements
     * @return List of sky regions to cover
     */
    std::vector<SkyRegion> computeSkyRegions(
        const std::vector<OrbitalElements>& asteroids);
    
    /**
     * @brief Convert sky regions to HEALPix tile IDs
     * @param regions Sky regions to cover
     * @param searchRadiusDeg Padding radius
     * @return Set of unique HEALPix IDs
     */
    std::unordered_set<uint64_t> regionsToHEALPix(
        const std::vector<SkyRegion>& regions,
        double searchRadiusDeg = 7.0);
    
    // ==================== Preloading Phase ====================
    
    /**
     * @brief Preload identified HEALPix tiles into RAM
     * @param parallel Use parallel loading (OpenMP)
     * @return Number of tiles loaded, total stars
     */
    std::pair<size_t, size_t> preloadTiles(bool parallel = true);
    
    /**
     * @brief Load specific HEALPix tile
     * @param healpixId HEALPix tile ID
     * @return Number of stars loaded
     */
    size_t loadTile(uint64_t healpixId);
    
    // ==================== Query Interface ====================
    
    /**
     * @brief Query stars in cone (in-memory lookup)
     * @param raDeg RA center (degrees)
     * @param decDeg Dec center (degrees)
     * @param radiusDeg Cone radius (degrees)
     * @param maxMag Maximum magnitude
     * @return List of stars in cone
     */
    std::vector<StarData> queryCone(double raDeg, double decDeg, 
                                   double radiusDeg, double maxMag = 18.0);
    
    /**
     * @brief Check if region is covered by preloaded tiles
     * @param raDeg RA center
     * @param decDeg Dec center
     * @param radiusDeg Radius
     * @return true if fully covered
     */
    bool isCovered(double raDeg, double decDeg, double radiusDeg);
    
    // ==================== Statistics ====================
    
    /**
     * @brief Get preloader statistics
     */
    struct PreloaderStats {
        size_t tiles_identified;    // HEALPix tiles identified
        size_t tiles_loaded;        // Tiles actually loaded
        size_t total_stars;         // Stars in RAM
        double memory_mb;           // Estimated RAM usage
        double preload_time_ms;     // Preload duration
        size_t query_count;         // Number of queries
        double avg_query_time_ms;   // Average query time
    };
    
    PreloaderStats getStats() const { return stats_; }
    void printStats() const;
    
private:
    std::string catalog_path_;
    int healpix_nside_;
    
    // Identified tiles (phase 1)
    std::unordered_set<uint64_t> needed_tiles_;
    
    // Loaded tiles (phase 2)
    std::unordered_map<uint64_t, HEALPixTile> tile_cache_;
    
    // Statistics
    PreloaderStats stats_;
    
    // HEALPix utilities
    uint64_t ang2pix(double ra, double dec);
    void pix2ang(uint64_t healpix_id, double& ra, double& dec);
    std::vector<uint64_t> query_disc(double ra, double dec, double radius);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_GAIA_CATALOG_PRELOADER_H
