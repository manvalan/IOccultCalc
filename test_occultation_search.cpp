#include "asteroid_occultation_search.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>

int main(int argc, char* argv[]) {
    // Default parameters (can be overridden via command line)
    std::string eq1_file = "17030_astdys.eq1";
    std::string designation = "17030";
    double start_mjd = 61007.0;
    double end_mjd = 61008.0;
    double diameter_km = 7.332;  // Default: 17030 Sierks
    
    // Parse command line arguments
    // Usage: ./test_occultation_search [eq1_file] [designation] [start_mjd] [end_mjd] [diameter_km]
    if (argc > 1) eq1_file = argv[1];
    if (argc > 2) designation = argv[2];
    if (argc > 3) start_mjd = std::atof(argv[3]);
    if (argc > 4) end_mjd = std::atof(argv[4]);
    if (argc > 5) diameter_km = std::atof(argv[5]);
    
    std::cout << "\n╔═══════════════════════════════════════════════════════════╗\n";
    std::cout << "║         ASTEROID OCCULTATION SEARCH - TEST              ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════╝\n\n";
    
    std::cout << "Configuration:\n";
    std::cout << "  Asteroid file:  " << eq1_file << "\n";
    std::cout << "  Designation:    " << designation << "\n";
    std::cout << "  Search window:  MJD " << std::fixed << std::setprecision(1) 
              << start_mjd << " - " << end_mjd << "\n";
    std::cout << "  Diameter:       " << std::setprecision(3) << diameter_km << " km\n\n";
    
    ioccult::AsteroidOccultationSearch search;
    
    // Initialize Gaia catalog
    const char* home = std::getenv("HOME");
    std::string gaia_dir = std::string(home) + "/.catalog/gaia_mag18_v2_multifile";
    if (!search.initializeGaiaCatalog(gaia_dir)) {
        std::cerr << "❌ Failed to initialize Gaia catalog\n";
        return 1;
    }
    
    // Load asteroid
    if (!search.loadAsteroidFromFile(eq1_file)) {
        std::cerr << "❌ Failed to load asteroid from " << eq1_file << "\n";
        return 1;
    }
    
    // Configure search parameters
    ioccult::OccultationSearchConfig config;
    config.search_start_mjd = start_mjd;
    config.search_end_mjd = end_mjd;
    config.corridor_width_deg = 0.005;  // ±18" search corridor
    config.max_approach_arcsec = 10.0;  // Filter candidates beyond this
    config.max_star_magnitude = 14.0;
    config.time_step_hours = 1.0;
    config.refinement_window_minutes = 10.0;
    config.refinement_step_seconds = 1.0;
    config.asteroid_diameter_km = diameter_km;
    config.enable_orbital_fitting = false;  // Can be enabled if needed
    config.asteroid_designation = designation;
    config.max_observations_for_fit = 25;
    
    // Execute search
    auto results = search.search(config);
    
    if (!results.success) {
        std::cerr << "❌ Search failed: " << results.error_message << "\n";
        return 1;
    }
    
    // Display results
    std::cout << "✓ Search completed successfully!\n\n";
    std::cout << "Found " << results.events.size() << " event(s)\n";
    
    if (results.events.empty()) {
        std::cout << "No close approaches or occultations found in the search window.\n";
        return 0;
    }
    
    std::cout << "\n";
    for (size_t i = 0; i < results.events.size(); ++i) {
        const auto& evt = results.events[i];
        std::cout << "─────────────────────────────────────────────────────────────\n";
        std::cout << "Event #" << (i+1);
        if (evt.is_occultation) {
            std::cout << " ★ OCCULTATION ★";
        }
        std::cout << "\n";
        std::cout << "  Gaia DR3 " << evt.source_id << "\n";
        std::cout << "  Star magnitude: " << std::fixed << std::setprecision(2) << evt.star_mag << "\n";
        std::cout << "  Closest approach: " << std::setprecision(3) << evt.refined_ca_arcsec 
                  << "\" @ MJD " << std::setprecision(6) << evt.refined_ca_mjd << "\n";
        std::cout << "  Angular velocity: " << std::setprecision(2) << evt.angular_velocity_arcsec_per_sec 
                  << " \"/sec\n";
        std::cout << "  Position angle: " << std::setprecision(1) << evt.position_angle_deg << "°\n";
        
        if (evt.is_occultation) {
            std::cout << "  Max duration: " << std::setprecision(2) << evt.max_duration_sec << " seconds\n";
        }
    }
    std::cout << "─────────────────────────────────────────────────────────────\n";
    
    // Summary statistics
    std::cout << "\nStatistics:\n";
    std::cout << "  Phase 1 candidates: " << results.total_candidates_found << "\n";
    std::cout << "  Phase 1 time: " << std::setprecision(1) << results.phase1_processing_time_ms << " ms\n";
    std::cout << "  Phase 2 time: " << results.phase2_processing_time_ms << " ms\n";
    std::cout << "  Total occultations: " << results.actual_occultations << "\n";
    
    if (results.orbit_fitting_attempted) {
        std::cout << "\nOrbital fitting:\n";
        std::cout << "  Status: " << (results.orbit_fitting_successful ? "SUCCESS" : "FAILED") << "\n";
        if (results.orbit_fitting_successful) {
            std::cout << "  Observations used: " << results.observations_used << "\n";
            std::cout << "  RMS residuals: " << std::setprecision(2) << results.fit_rms_arcsec << " arcsec\n";
            std::cout << "  Chi²: " << std::setprecision(2) << results.fit_chi_squared << "\n";
        }
    }
    
    std::cout << "\n✓ Test completed!\n\n";
    return 0;
}
