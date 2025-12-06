/**
 * @file asteroid_occultation_search.h
 * @brief High-level API for two-phase asteroid occultation search
 */

#ifndef ASTEROID_OCCULTATION_SEARCH_H
#define ASTEROID_OCCULTATION_SEARCH_H

#include "phase1_candidate_screening.h"
#include "phase2_geocentric.h"
#include <string>
#include <vector>
#include <memory>

namespace ioccult {

struct OccultationSearchConfig {
    // Asteroid
    std::string asteroid_eq1_file;
    std::string asteroid_designation = "";
    double asteroid_diameter_km = 10.0;
    
    // Time window
    double search_start_mjd = 0.0;
    double search_end_mjd = 0.0;
    
    // Phase 1
    double max_approach_arcsec = 2.0;
    double min_star_magnitude = 0.0;
    double max_star_magnitude = 16.0;
    double time_step_hours = 1.0;
    
    // Phase 2
    bool enable_orbital_fitting = false;
    int max_observations_for_fit = 25;
    bool use_all_observations = false;
    double refinement_window_minutes = 10.0;
    double refinement_step_seconds = 1.0;
    
    // Filters
    bool only_occultations = false;
};

struct OccultationSearchResult {
    int total_candidates_found = 0;
    double phase1_processing_time_ms = 0.0;
    int total_events_refined = 0;
    int actual_occultations = 0;
    double phase2_processing_time_ms = 0.0;
    
    bool orbit_fitting_attempted = false;
    bool orbit_fitting_successful = false;
    int observations_downloaded = 0;
    int observations_used = 0;
    double fit_rms_arcsec = 0.0;
    double fit_chi_squared = 0.0;
    std::string fit_notes;
    
    std::vector<Phase2GeocentricResult> events;
    
    bool success = false;
    std::string error_message;
};

class AsteroidOccultationSearch {
public:
    AsteroidOccultationSearch();
    ~AsteroidOccultationSearch();
    
    bool initializeGaiaCatalog(const std::string& gaia_cache_dir);
    bool loadAsteroidFromFile(const std::string& eq1_path);
    void setAsteroidDiameter(double diameter_km);
    
    OccultationSearchResult search(const OccultationSearchConfig& config);
    
    ioccultcalc::Phase1CandidateScreening& getPhase1();
    Phase2Geocentric& getPhase2();

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

} // namespace ioccult

#endif
