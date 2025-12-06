#pragma once

#ifndef PHASE2_GEOCENTRIC_H
#define PHASE2_GEOCENTRIC_H

#include "phase1_candidate_screening.h"
#include <string>
#include <vector>
#include <memory>

// Forward declarations
namespace astdyn {
    namespace propagation {
        struct KeplerianElements;
    }
}

namespace ioccult {

struct Phase2GeocentricConfig {
    bool fit_orbit_from_rwo = false;
    std::string asteroid_designation = "";
    int max_observations_for_fit = 25;
    bool use_all_available_observations = false;
    double observation_arc_days = 365.0;
    double time_window_minutes = 10.0;
    double time_step_seconds = 1.0;
    double asteroid_diameter_km = 10.0;
    bool apply_light_time_correction = true;
    bool apply_stellar_aberration = true;
    bool apply_planetary_aberration = true;
};

struct Phase2GeocentricResult {
    uint64_t source_id = 0;
    double refined_ca_mjd = 0.0;
    double refined_ca_arcsec = 0.0;
    double position_angle_deg = 0.0;
    double angular_velocity_arcsec_per_sec = 0.0;
    double asteroid_distance_au = 0.0;
    double asteroid_angular_diameter_mas = 0.0;
    bool is_occultation = false;
    double max_duration_sec = 0.0;
    double miss_distance_arcsec = 0.0;
    double star_ra_deg = 0.0;
    double star_dec_deg = 0.0;
    double star_mag = 0.0;
    double asteroid_ra_deg = 0.0;
    double asteroid_dec_deg = 0.0;
    std::string error_message;
    bool success = false;
    bool orbit_was_refined = false;
};

struct Phase2GeocentricResults {
    std::vector<Phase2GeocentricResult> results;
    int total_candidates = 0;
    int occultations_found = 0;
    double processing_time_ms = 0.0;
    bool orbit_refinement_attempted = false;
    bool orbit_refinement_successful = false;
    int total_observations_downloaded = 0;
    int observations_used_in_fit = 0;
    double fit_rms_arcsec = 0.0;
    double fit_chi_squared = 0.0;
    std::string fit_notes;
};

class Phase2Geocentric {
public:
    Phase2Geocentric();
    ~Phase2Geocentric();
    
    bool loadAsteroidFromEQ1(const std::string& eq1_path);
    void setOrbitalElements(const astdyn::propagation::KeplerianElements& elements);
    void setAsteroidDiameter(double diameter_km);
    
    Phase2GeocentricResult refineCandidate(
        const ioccultcalc::CandidateStar& candidate,
        const Phase2GeocentricConfig& config = Phase2GeocentricConfig());
    
    Phase2GeocentricResults refineAllCandidates(
        const std::vector<ioccultcalc::CandidateStar>& candidates,
        const Phase2GeocentricConfig& config = Phase2GeocentricConfig());
    
private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

} // namespace ioccult

#endif
