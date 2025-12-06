/**
 * @file asteroid_occultation_search.cpp
 * @brief Implementation of high-level occultation search API
 */

#include "asteroid_occultation_search.h"
#include <iostream>
#include <chrono>
#include "ioc_gaialib/unified_gaia_catalog.h"

namespace ioccult {

struct AsteroidOccultationSearch::Impl {
    ioccultcalc::Phase1CandidateScreening phase1;
    Phase2Geocentric phase2;
    bool gaia_initialized = false;
    bool asteroid_loaded = false;
    double asteroid_diameter_km = 10.0;
};

AsteroidOccultationSearch::AsteroidOccultationSearch() 
    : pimpl_(std::make_unique<Impl>()) {
}

AsteroidOccultationSearch::~AsteroidOccultationSearch() = default;

bool AsteroidOccultationSearch::initializeGaiaCatalog(const std::string& gaia_cache_dir) {
    try {
        // Crea stringa JSON per UnifiedGaiaCatalog
        std::string json_config = R"({"catalog_type": "multifile_v2", "multifile_directory": ")" 
            + gaia_cache_dir + R"("})";
        
        // Inizializza il singleton
        if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) {
            throw std::runtime_error("Failed to initialize UnifiedGaiaCatalog");
        }
        
        // Ottieni il singleton e collegalo a Phase1
        auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
        pimpl_->phase1.setCatalog(&catalog);
        pimpl_->gaia_initialized = true;
        std::cout << "✓ Gaia catalog initialized\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "❌ Failed to initialize Gaia catalog: " << e.what() << "\n";
        return false;
    }
}

bool AsteroidOccultationSearch::loadAsteroidFromFile(const std::string& eq1_path) {
    bool p1_ok = pimpl_->phase1.loadAsteroidFromEQ1(eq1_path);
    bool p2_ok = pimpl_->phase2.loadAsteroidFromEQ1(eq1_path);
    
    if (p1_ok && p2_ok) {
        pimpl_->asteroid_loaded = true;
        std::cout << "✓ Asteroid orbital elements loaded\n";
        return true;
    } else {
        std::cerr << "❌ Failed to load asteroid elements\n";
        return false;
    }
}

void AsteroidOccultationSearch::setAsteroidDiameter(double diameter_km) {
    pimpl_->asteroid_diameter_km = diameter_km;
}

ioccultcalc::Phase1CandidateScreening& AsteroidOccultationSearch::getPhase1() {
    return pimpl_->phase1;
}

ioccult::Phase2Geocentric& AsteroidOccultationSearch::getPhase2() {
    return pimpl_->phase2;
}

OccultationSearchResult AsteroidOccultationSearch::search(const OccultationSearchConfig& config) {
    OccultationSearchResult result;
    
    if (!pimpl_->gaia_initialized) {
        result.error_message = "Gaia catalog not initialized";
        return result;
    }
    
    if (!pimpl_->asteroid_loaded && config.asteroid_eq1_file.empty() && config.asteroid_designation.empty()) {
        result.error_message = "No asteroid data provided";
        return result;
    }
    
    if (config.search_start_mjd >= config.search_end_mjd) {
        result.error_message = "Invalid time window";
        return result;
    }
    
    if (!config.asteroid_eq1_file.empty() && !pimpl_->asteroid_loaded) {
        if (!loadAsteroidFromFile(config.asteroid_eq1_file)) {
            result.error_message = "Failed to load asteroid from file";
            return result;
        }
    }
    
    if (config.asteroid_diameter_km > 0) {
        setAsteroidDiameter(config.asteroid_diameter_km);
    }
    
    std::cout << "\n╔═══════════════════════════════════════════════════════════╗\n";
    std::cout << "║     ASTEROID OCCULTATION SEARCH - TWO PHASE PIPELINE     ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════╝\n\n";
    
    // PHASE 1
    std::cout << "┌─────────────────────────────────────────────────────────┐\n";
    std::cout << "│ PHASE 1: CANDIDATE SCREENING                            │\n";
    std::cout << "└─────────────────────────────────────────────────────────┘\n\n";
    
    auto t1_start = std::chrono::high_resolution_clock::now();
    
    ioccultcalc::Phase1Config phase1_config;
    phase1_config.start_mjd_tdb = config.search_start_mjd;
    phase1_config.end_mjd_tdb = config.search_end_mjd;
    phase1_config.path_interval_seconds = (int)(config.time_step_hours * 3600);
    phase1_config.closest_approach_threshold_arcsec = config.max_approach_arcsec;
    // min mag not available
    phase1_config.max_magnitude = config.max_star_magnitude;
    
    auto phase1_results = pimpl_->phase1.screenCandidates(phase1_config);
    
    auto t1_end = std::chrono::high_resolution_clock::now();
    result.phase1_processing_time_ms = std::chrono::duration<double, std::milli>(t1_end - t1_start).count();
    result.total_candidates_found = static_cast<int>(phase1_results.candidates.size());
    
    std::cout << "📊 Phase 1 Results:\n";
    std::cout << "   Candidates found: " << result.total_candidates_found << "\n";
    for (const auto& c : phase1_results.candidates) {
        std::cout << "   • Gaia " << c.source_id 
                  << " CA=" << std::fixed << std::setprecision(2) << c.closest_approach_arcsec 
                  << "\" @ MJD " << std::setprecision(6) << c.closest_approach_mjd << "\n";
    }
    std::cout << "   Processing time: " << result.phase1_processing_time_ms << " ms\n\n";
    
    if (phase1_results.candidates.empty()) {
        std::cout << "ℹ No candidates found. Search complete.\n";
        result.success = true;
        return result;
    }
    
    // PHASE 2
    std::cout << "┌─────────────────────────────────────────────────────────┐\n";
    std::cout << "│ PHASE 2: GEOCENTRIC REFINEMENT (RKF78)                 │\n";
    std::cout << "└─────────────────────────────────────────────────────────┘\n\n";
    
    auto t2_start = std::chrono::high_resolution_clock::now();
    
    // Passa elementi orbitali da Phase1 a Phase2
    pimpl_->phase2.setOrbitalElements(pimpl_->phase1.getOrbitalElements());
    
    Phase2GeocentricConfig phase2_config;
    phase2_config.time_window_minutes = config.refinement_window_minutes;
    phase2_config.time_step_seconds = config.refinement_step_seconds;
    phase2_config.asteroid_diameter_km = pimpl_->asteroid_diameter_km;
    phase2_config.fit_orbit_from_rwo = false;  // Default
    
    if (config.enable_orbital_fitting && !config.asteroid_designation.empty()) {
        phase2_config.fit_orbit_from_rwo = true;
        phase2_config.asteroid_designation = config.asteroid_designation;
        phase2_config.max_observations_for_fit = config.max_observations_for_fit;
        phase2_config.use_all_available_observations = config.use_all_observations;
    }
    
    auto phase2_results = pimpl_->phase2.refineAllCandidates(phase1_results.candidates, phase2_config);
    
    auto t2_end = std::chrono::high_resolution_clock::now();
    result.phase2_processing_time_ms = std::chrono::duration<double, std::milli>(t2_end - t2_start).count();
    
    result.orbit_fitting_attempted = phase2_results.orbit_refinement_attempted;
    result.orbit_fitting_successful = phase2_results.orbit_refinement_successful;
    result.observations_downloaded = phase2_results.total_observations_downloaded;
    result.observations_used = phase2_results.observations_used_in_fit;
    result.fit_rms_arcsec = phase2_results.fit_rms_arcsec;
    result.fit_chi_squared = phase2_results.fit_chi_squared;
    result.fit_notes = phase2_results.fit_notes;
    
    std::cout << "📊 Phase 2 Results:\n";
    std::cout << "   Events refined: " << phase2_results.results.size() << "\n";
    std::cout << "   Occultations found: " << phase2_results.occultations_found << "\n";
    std::cout << "   Processing time: " << result.phase2_processing_time_ms << " ms\n\n";
    
    if (config.only_occultations) {
        for (const auto& evt : phase2_results.results) {
            if (evt.is_occultation) {
                result.events.push_back(evt);
            }
        }
        result.actual_occultations = static_cast<int>(result.events.size());
    } else {
        result.events = phase2_results.results;
        result.actual_occultations = phase2_results.occultations_found;
    }
    
    result.total_events_refined = static_cast<int>(result.events.size());
    result.success = true;
    
    std::cout << "╔═══════════════════════════════════════════════════════════╗\n";
    std::cout << "║                   SEARCH COMPLETE                         ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════╝\n\n";
    
    return result;
}

} // namespace ioccult
