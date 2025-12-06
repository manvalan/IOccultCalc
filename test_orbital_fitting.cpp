/**
 * Test Orbital Fitting con RWO da AstDyS
 * Asteroide 17030 - Download osservazioni e fit orbitale
 */
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "phase1_candidate_screening.h"
#include "phase2_geocentric.h"
#include "ioc_gaialib/unified_gaia_catalog.h"

using namespace ioccultcalc;

int main() {
    std::cout << "\n═══ Test Orbital Fitting ═══\n";
    std::cout << "Asteroide 17030 - Download RWO + Fit\n\n";
    
    // Inizializza Gaia
    const char* home = std::getenv("HOME");
    std::string json_config = R"({"catalog_type": "multifile_v2", "multifile_directory": ")" 
        + std::string(home) + R"(/.catalog/gaia_mag18_v2_multifile"})";
    
    if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) {
        std::cerr << "✗ Errore Gaia\n";
        return 1;
    }
    std::cout << "✓ Gaia initialized\n\n";
    
    // PHASE 1: Trova candidate
    std::cout << "=== PHASE 1: SCREENING ===\n";
    Phase1CandidateScreening phase1;
    phase1.loadAsteroidFromEQ1("17030_astdys.eq1");
    
    Phase1Config config1;
    config1.start_mjd_tdb = 61007.0;  // 28 Nov 2025
    config1.end_mjd_tdb = 61008.0;
    config1.path_interval_seconds = 3600;
    config1.corridor_width_deg = 0.005;  // ±18"
    config1.max_magnitude = 14.0;
    config1.closest_approach_threshold_arcsec = 10.0;
    
    auto results1 = phase1.screenCandidates(config1);
    
    std::cout << "Candidate trovate: " << results1.candidates.size() << "\n\n";
    
    if (results1.candidates.empty()) {
        std::cout << "Nessuna candidata - termino.\n";
        return 0;
    }
    
    // Mostra stella target
    for (const auto& c : results1.candidates) {
        std::cout << "Gaia " << c.source_id 
                  << " CA=" << std::fixed << std::setprecision(2) << c.closest_approach_arcsec << "\"\n";
    }
    
    // PHASE 2: CON ORBITAL FITTING
    std::cout << "\n=== PHASE 2: CON ORBITAL FITTING ===\n";
    ioccult::Phase2Geocentric phase2;
    phase2.setOrbitalElements(phase1.getOrbitalElements());
    
    ioccult::Phase2GeocentricConfig config2;
    config2.fit_orbit_from_rwo = true;  // ← ATTIVA FITTING
    config2.asteroid_designation = "17030";
    config2.max_observations_for_fit = 25;
    config2.use_all_available_observations = false;
    config2.time_window_minutes = 10.0;
    config2.time_step_seconds = 1.0;
    config2.asteroid_diameter_km = 10.0;
    
    auto results2 = phase2.refineAllCandidates(results1.candidates, config2);
    
    std::cout << "\n=== RISULTATI ===\n";
    std::cout << "Orbital fitting attempted: " << (results2.orbit_refinement_attempted ? "YES" : "NO") << "\n";
    std::cout << "Orbital fitting successful: " << (results2.orbit_refinement_successful ? "YES" : "NO") << "\n";
    
    if (results2.orbit_refinement_successful) {
        std::cout << "\n✓✓✓ ORBIT REFINED ✓✓✓\n";
        std::cout << "Observations downloaded: " << results2.total_observations_downloaded << "\n";
        std::cout << "Observations used: " << results2.observations_used_in_fit << "\n";
        std::cout << "RMS residuals: " << std::fixed << std::setprecision(3) 
                  << results2.fit_rms_arcsec << " arcsec\n";
        std::cout << "Chi²: " << std::setprecision(2) << results2.fit_chi_squared << "\n";
        std::cout << "Notes: " << results2.fit_notes << "\n";
    } else {
        std::cout << "\n⚠ Fitting failed\n";
        std::cout << "Notes: " << results2.fit_notes << "\n";
    }
    
    // Mostra eventi
    std::cout << "\n=== EVENTI (con elementi " 
              << (results2.orbit_refinement_successful ? "RAFFINATI" : "NOMINALI") << ") ===\n";
    
    for (size_t i = 0; i < results2.results.size(); ++i) {
        const auto& e = results2.results[i];
        std::cout << "\n" << (i+1) << ") Gaia " << e.source_id;
        if (e.orbit_was_refined) std::cout << " [REFINED ORBIT]";
        std::cout << "\n";
        std::cout << "   CA: " << std::fixed << std::setprecision(3) << e.refined_ca_arcsec 
                  << "\" @ MJD " << std::setprecision(6) << e.refined_ca_mjd << "\n";
        std::cout << "   Velocità: " << std::setprecision(2) << e.angular_velocity_arcsec_per_sec << " \"/sec\n";
        std::cout << "   PA: " << std::setprecision(1) << e.position_angle_deg << "°\n";
        if (e.is_occultation) {
            std::cout << "   ★ OCCULTAZIONE! Durata: " << std::setprecision(2) << e.max_duration_sec << " sec\n";
        }
    }
    
    std::cout << "\n✓ Test completato!\n\n";
    return 0;
}
