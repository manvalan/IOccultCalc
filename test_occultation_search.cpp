#include "asteroid_occultation_search.h"
#include <iostream>
#include <iomanip>

int main() {
    std::cout << "\n═══ TEST: AsteroidOccultationSearch ═══\n\n";
    
    ioccult::AsteroidOccultationSearch search;
    
    // Percorso Gaia multifile_v2
    const char* home = std::getenv("HOME");
    std::string gaia_dir = std::string(home) + "/.catalog/gaia_mag18_v2_multifile";
    if (!search.initializeGaiaCatalog(gaia_dir)) {
        return 1;
    }
    
    if (!search.loadAsteroidFromFile("17030_astdys.eq1")) {
        return 1;
    }
    
    //search.setAsteroidDiameter(10.0);
    search.setAsteroidDiameter(1000.0);  // 1000 km invece di 10
    
    ioccult::OccultationSearchConfig config;
    config.search_start_mjd = 61000.0;  // Esempio
    config.search_end_mjd = 61010.0;
    config.max_approach_arcsec = 2.0;
    config.max_star_magnitude = 14.0;
    config.time_step_hours = 1.0;
    config.refinement_window_minutes = 10.0;
    config.refinement_step_seconds = 1.0;
    config.enable_orbital_fitting = false; // Per test iniziale
    config.asteroid_designation = "17030";
    config.max_observations_for_fit = 25;
    config.asteroid_diameter_km = 1000.0; // 1000 km per test occultazione
    
    auto results = search.search(config);
    
    if (!results.success) {
        std::cerr << "❌ Search failed: " << results.error_message << "\n";
        return 1;
    }
    
    std::cout << "✓ Found " << results.events.size() << " event(s)\n";
    
    for (size_t i = 0; i < results.events.size(); ++i) {
        const auto& evt = results.events[i];
        std::cout << "\nEvent #" << (i+1) << ":\n";
        std::cout << "  Gaia " << evt.source_id << "\n";
        std::cout << "  CA: " << std::fixed << std::setprecision(3) 
                  << evt.refined_ca_arcsec << "\" @ MJD " << evt.refined_ca_mjd << "\n";
        std::cout << "  Occultation: " << (evt.is_occultation ? "YES" : "NO") << "\n";
    }
    
    std::cout << "\n✓ Test completed!\n";
    return 0;
}
