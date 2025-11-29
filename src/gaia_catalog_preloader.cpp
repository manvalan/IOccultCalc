/**
 * @file gaia_catalog_preloader.cpp
 * @brief Implementation of smart Gaia catalog preloader - Phase 2
 */

#include "ioccultcalc/gaia_catalog_preloader.h"
#include "ioccultcalc/gaia_cache.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ioccultcalc {

namespace {
    constexpr double DEG2RAD = M_PI / 180.0;
    constexpr double RAD2DEG = 180.0 / M_PI;
}

GaiaCatalogPreloader::GaiaCatalogPreloader(const std::string& catalogPath, int healpixNSIDE)
    : catalog_path_(catalogPath), healpix_nside_(healpixNSIDE) {
    
    stats_ = PreloaderStats();
    stats_.tiles_identified = 0;
    stats_.tiles_loaded = 0;
    stats_.total_stars = 0;
    stats_.memory_mb = 0.0;
    stats_.preload_time_ms = 0.0;
    stats_.query_count = 0;
    stats_.avg_query_time_ms = 0.0;
    
    std::cout << "[GaiaCatalogPreloader] Initialized with catalog: " << catalogPath << "\n";
    std::cout << "[GaiaCatalogPreloader] HEALPix NSIDE=" << healpixNSIDE << "\n";
}

GaiaCatalogPreloader::~GaiaCatalogPreloader() {
    // Cleanup
}

// ============================================================================
// Pre-Analysis Phase
// ============================================================================

size_t GaiaCatalogPreloader::analyzeTrajectories(
    const std::vector<OrbitalElements>& asteroids,
    double searchRadiusDeg) {
    
    std::cout << "\n[Phase 2.1] Analyzing " << asteroids.size() 
             << " asteroid trajectories...\n";
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Step 1: Compute rough sky regions
    auto regions = computeSkyRegions(asteroids);
    std::cout << "  → Identified " << regions.size() << " sky regions\n";
    
    // Step 2: Convert to HEALPix tiles
    needed_tiles_ = regionsToHEALPix(regions, searchRadiusDeg);
    stats_.tiles_identified = needed_tiles_.size();
    
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    std::cout << "  → HEALPix tiles needed: " << needed_tiles_.size() 
             << " out of " << (12 * healpix_nside_ * healpix_nside_) << " total\n";
    std::cout << "  → Coverage: " << std::fixed << std::setprecision(2)
             << (100.0 * needed_tiles_.size() / (12 * healpix_nside_ * healpix_nside_))
             << "% of sky\n";
    std::cout << "  → Analysis time: " << elapsed_ms << " ms\n";
    
    return needed_tiles_.size();
}

std::vector<SkyRegion> GaiaCatalogPreloader::computeSkyRegions(
    const std::vector<OrbitalElements>& asteroids) {
    
    std::vector<SkyRegion> regions;
    const double merge_distance_deg = 15.0; // Merge regions within 15°
    
    for (const auto& ast : asteroids) {
        // Rough position estimate from orbital elements
        // (same as batching logic)
        double ra = ast.Omega * RAD2DEG;
        double dec = (ast.i * RAD2DEG) - 20.0;
        
        // Try to merge with existing region
        bool merged = false;
        for (auto& region : regions) {
            double dra = std::abs(region.ra_center - ra);
            double ddec = std::abs(region.dec_center - dec);
            
            if (dra < merge_distance_deg && ddec < merge_distance_deg) {
                // Expand region bounds
                region.ra_min = std::min(region.ra_min, ra - 7.0);
                region.ra_max = std::max(region.ra_max, ra + 7.0);
                region.dec_min = std::min(region.dec_min, dec - 7.0);
                region.dec_max = std::max(region.dec_max, dec + 7.0);
                
                // Update center
                region.ra_center = (region.ra_min + region.ra_max) / 2.0;
                region.dec_center = (region.dec_min + region.dec_max) / 2.0;
                
                merged = true;
                break;
            }
        }
        
        if (!merged) {
            // New region
            SkyRegion region;
            region.ra_min = ra - 7.0;
            region.ra_max = ra + 7.0;
            region.dec_min = dec - 7.0;
            region.dec_max = dec + 7.0;
            region.ra_center = ra;
            region.dec_center = dec;
            regions.push_back(region);
        }
    }
    
    return regions;
}

std::unordered_set<uint64_t> GaiaCatalogPreloader::regionsToHEALPix(
    const std::vector<SkyRegion>& regions,
    double searchRadiusDeg) {
    
    std::unordered_set<uint64_t> healpix_ids;
    
    for (const auto& region : regions) {
        // Query HEALPix disc around region center
        double radius_deg = std::max(
            std::abs(region.ra_max - region.ra_min),
            std::abs(region.dec_max - region.dec_min)
        ) / 2.0 + searchRadiusDeg;
        
        auto disc_pixels = query_disc(region.ra_center, region.dec_center, radius_deg);
        healpix_ids.insert(disc_pixels.begin(), disc_pixels.end());
    }
    
    return healpix_ids;
}

// ============================================================================
// Preloading Phase
// ============================================================================

std::pair<size_t, size_t> GaiaCatalogPreloader::preloadTiles(bool parallel) {
    std::cout << "\n[Phase 2.2] Preloading " << needed_tiles_.size() 
             << " HEALPix tiles";
    if (parallel) std::cout << " (parallel)";
    std::cout << "...\n";
    
    auto start = std::chrono::high_resolution_clock::now();
    
    size_t tiles_loaded = 0;
    size_t total_stars = 0;
    
    // Convert set to vector for OpenMP
    std::vector<uint64_t> tile_ids(needed_tiles_.begin(), needed_tiles_.end());
    
#ifdef _OPENMP
    if (parallel) {
        #pragma omp parallel for reduction(+:tiles_loaded,total_stars) schedule(dynamic)
        for (size_t i = 0; i < tile_ids.size(); i++) {
            size_t stars = loadTile(tile_ids[i]);
            tiles_loaded++;
            total_stars += stars;
            
            if (i % 10 == 0) {
                #pragma omp critical
                {
                    std::cout << "  [" << (i+1) << "/" << tile_ids.size() << "] "
                             << "Loaded " << tiles_loaded << " tiles, "
                             << total_stars << " stars\r" << std::flush;
                }
            }
        }
    } else
#endif
    {
        // Sequential loading
        for (size_t i = 0; i < tile_ids.size(); i++) {
            size_t stars = loadTile(tile_ids[i]);
            tiles_loaded++;
            total_stars += stars;
            
            if (i % 10 == 0 || i == tile_ids.size() - 1) {
                std::cout << "  [" << (i+1) << "/" << tile_ids.size() << "] "
                         << "Loaded " << tiles_loaded << " tiles, "
                         << total_stars << " stars\r" << std::flush;
            }
        }
    }
    
    std::cout << "\n";
    
    auto end = std::chrono::high_resolution_clock::now();
    stats_.preload_time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    stats_.tiles_loaded = tiles_loaded;
    stats_.total_stars = total_stars;
    
    // Estimate memory (rough: 200 bytes/star)
    stats_.memory_mb = (total_stars * 200.0) / (1024.0 * 1024.0);
    
    std::cout << "  ✓ Preloaded " << tiles_loaded << " tiles in " 
             << stats_.preload_time_ms << " ms\n";
    std::cout << "  ✓ Total stars in RAM: " << total_stars << "\n";
    std::cout << "  ✓ Estimated RAM: " << std::fixed << std::setprecision(1)
             << stats_.memory_mb << " MB\n\n";
    
    return {tiles_loaded, total_stars};
}

size_t GaiaCatalogPreloader::loadTile(uint64_t healpixId) {
    // TODO: Implement actual HEALPix tile loading via IOC_GaiaLib
    // For now, placeholder that returns 0
    
    // This will use GaiaCatalog::queryCone() with HEALPix bounds
    // and cache the results
    
    HEALPixTile tile;
    tile.healpix_id = healpixId;
    tile.loaded = false;
    tile.star_count = 0;
    
    // Placeholder: mark as loaded
    #pragma omp critical
    {
        tile_cache_[healpixId] = tile;
    }
    
    return 0; // TODO: return actual star count
}

// ============================================================================
// Query Interface
// ============================================================================

std::vector<StarData> GaiaCatalogPreloader::queryCone(
    double raDeg, double decDeg, double radiusDeg, double maxMag) {
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<StarData> results;
    
    // Get HEALPix tiles intersecting cone
    auto tiles = query_disc(raDeg, decDeg, radiusDeg);
    
    // Query each tile
    for (auto tile_id : tiles) {
        if (tile_cache_.count(tile_id) && tile_cache_[tile_id].loaded) {
            // Filter stars in cone
            for (const auto& star : tile_cache_[tile_id].stars) {
                double star_ra = star.position.ra * RAD2DEG;
                double star_dec = star.position.dec * RAD2DEG;
                
                // Angular distance
                double dra = (star_ra - raDeg) * std::cos(decDeg * DEG2RAD);
                double ddec = star_dec - decDeg;
                double dist = std::sqrt(dra*dra + ddec*ddec);
                
                if (dist <= radiusDeg && star.G_mag <= maxMag) {
                    results.push_back(star);
                }
            }
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    double query_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    stats_.query_count++;
    stats_.avg_query_time_ms = 
        (stats_.avg_query_time_ms * (stats_.query_count - 1) + query_ms) / stats_.query_count;
    
    return results;
}

bool GaiaCatalogPreloader::isCovered(double raDeg, double decDeg, double radiusDeg) {
    auto tiles = query_disc(raDeg, decDeg, radiusDeg);
    
    for (auto tile_id : tiles) {
        if (!tile_cache_.count(tile_id) || !tile_cache_[tile_id].loaded) {
            return false;
        }
    }
    
    return true;
}

// ============================================================================
// HEALPix Utilities (simplified - use proper library in production)
// ============================================================================

uint64_t GaiaCatalogPreloader::ang2pix(double ra, double dec) {
    // Simplified HEALPix conversion
    // TODO: Use proper HEALPix library (healpy or C++ equivalent)
    // For now, rough approximation
    
    int npix = 12 * healpix_nside_ * healpix_nside_;
    
    // Very rough: divide sky into pixels
    int ra_bin = static_cast<int>((ra / 360.0) * healpix_nside_ * 4);
    int dec_bin = static_cast<int>(((dec + 90.0) / 180.0) * healpix_nside_ * 3);
    
    return (dec_bin * healpix_nside_ * 4 + ra_bin) % npix;
}

void GaiaCatalogPreloader::pix2ang(uint64_t healpix_id, double& ra, double& dec) {
    // Simplified inverse
    // TODO: Use proper HEALPix library
    
    int npix = 12 * healpix_nside_ * healpix_nside_;
    int width = healpix_nside_ * 4;
    
    int dec_bin = healpix_id / width;
    int ra_bin = healpix_id % width;
    
    ra = (ra_bin * 360.0) / width;
    dec = (dec_bin * 180.0) / (healpix_nside_ * 3) - 90.0;
}

std::vector<uint64_t> GaiaCatalogPreloader::query_disc(
    double ra, double dec, double radius) {
    
    std::vector<uint64_t> pixels;
    
    // TODO: Proper HEALPix disc query
    // For now, simple grid search around center
    
    int n_sample = static_cast<int>(radius * 2) + 1;
    for (int i = -n_sample; i <= n_sample; i++) {
        for (int j = -n_sample; j <= n_sample; j++) {
            double test_ra = ra + i * (radius / n_sample);
            double test_dec = dec + j * (radius / n_sample);
            
            // Check if in disc
            double dra = (test_ra - ra) * std::cos(dec * DEG2RAD);
            double ddec = test_dec - dec;
            double dist = std::sqrt(dra*dra + ddec*ddec);
            
            if (dist <= radius) {
                uint64_t pix = ang2pix(test_ra, test_dec);
                pixels.push_back(pix);
            }
        }
    }
    
    // Remove duplicates
    std::sort(pixels.begin(), pixels.end());
    pixels.erase(std::unique(pixels.begin(), pixels.end()), pixels.end());
    
    return pixels;
}

// ============================================================================
// Statistics
// ============================================================================

void GaiaCatalogPreloader::printStats() const {
    std::cout << "\n╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║     GaiaCatalogPreloader Statistics                   ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n\n";
    std::cout << "  HEALPix tiles identified: " << stats_.tiles_identified << "\n";
    std::cout << "  HEALPix tiles loaded:     " << stats_.tiles_loaded << "\n";
    std::cout << "  Total stars in RAM:       " << stats_.total_stars << "\n";
    std::cout << "  Estimated RAM usage:      " << std::fixed << std::setprecision(1)
             << stats_.memory_mb << " MB\n";
    std::cout << "  Preload time:             " << stats_.preload_time_ms << " ms\n";
    std::cout << "  Queries performed:        " << stats_.query_count << "\n";
    std::cout << "  Average query time:       " << std::fixed << std::setprecision(3)
             << stats_.avg_query_time_ms << " ms (in-memory!)\n\n";
}

} // namespace ioccultcalc
