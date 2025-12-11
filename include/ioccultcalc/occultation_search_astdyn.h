/**
 * @file occultation_search_astdyn.h
 * @brief Wrapper unificato per ricerca occultazioni usando ITALOccultLibrary/AstDyn
 */

#ifndef IOCCULTCALC_OCCULTATION_SEARCH_ASTDYN_H
#define IOCCULTCALC_OCCULTATION_SEARCH_ASTDYN_H

#include "astdyn_interface.h"
#include "phase1_candidate_screening.h"
#include "phase2_occultation_geometry.h"
#include "occultation_predictor.h"  // ← AGGIUNGI per OccultationEvent
#include "orbital_elements.h"
#include "types.h"
#include <vector>
#include <string>
#include <memory>

namespace ioccultcalc {

struct OccultationSearchConfig {
    double start_mjd_tdb;
    double end_mjd_tdb;
    Phase1Config phase1_config;
    Phase2Config phase2_config;
    double asteroid_diameter_km;
    double orbital_uncertainty_km;
    double min_probability;
    double max_magnitude;
    
    OccultationSearchConfig()
        : start_mjd_tdb(0.0)
        , end_mjd_tdb(0.0)
        , phase1_config(Phase1Config::conservative())
        , asteroid_diameter_km(0.0)
        , orbital_uncertainty_km(100.0)
        , min_probability(0.01)
        , max_magnitude(18.0)
    {}
};

struct OccultationSearchResults {
    std::vector<::ioccultcalc::OccultationEvent> events;
    int num_candidates_found;
    int num_stars_in_corridor;
    int num_events_calculated;
    int num_events_failed;
    double phase1_time_ms;
    double phase2_time_ms;
    double total_time_ms;
    std::string asteroid_name;
    int asteroid_number;
    
    OccultationSearchResults()
        : num_candidates_found(0)
        , num_stars_in_corridor(0)
        , num_events_calculated(0)
        , num_events_failed(0)
        , phase1_time_ms(0.0)
        , phase2_time_ms(0.0)
        , total_time_ms(0.0)
        , asteroid_number(0)
    {}
};

class OccultationSearchAstDyn {
public:
    explicit OccultationSearchAstDyn(double tolerance = 1e-12);
    ~OccultationSearchAstDyn();
    
    OccultationSearchAstDyn(const OccultationSearchAstDyn&) = delete;
    OccultationSearchAstDyn& operator=(const OccultationSearchAstDyn&) = delete;
    
    bool loadAsteroid(int asteroid_number);
    bool loadAsteroid(const std::string& designation);
    bool loadAsteroidFromEQ1(const std::string& eq1_path);
    bool loadAsteroidFromJSON(int asteroid_number, const std::string& json_path = "");
    void setAsteroidElements(const AstDySElements& elements);
    
    OccultationSearchResults search(const OccultationSearchConfig& config);
    OccultationSearchResults search(
        double start_mjd_tdb,
        double end_mjd_tdb,
        double max_magnitude = 18.0);
    
    void setAsteroidDiameter(double diameter_km);
    void setOrbitalUncertainty(double sigma_km);
    void setGaiaCatalog(ioc::gaia::UnifiedGaiaCatalog* catalog);
    void addObserverSite(const std::string& name,
                        double lat_deg, double lon_deg, double elev_m);
    
    const AstDySElements& getAsteroidElements() const;
    const astdyn::propagation::KeplerianElements& getKeplerianElements() const;
    bool hasAsteroid() const;
    std::string getAsteroidName() const;
    int getAsteroidNumber() const;
    
private:
    class Impl;
    std::unique_ptr<Impl> pimpl_;
    void syncOrbitalElements();
    astdyn::propagation::KeplerianElements convertToKeplerian(
        const AstDySElements& astdys) const;
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_OCCULTATION_SEARCH_ASTDYN_H