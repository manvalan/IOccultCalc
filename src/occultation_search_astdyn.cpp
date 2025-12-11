/**
 * @file occultation_search_astdyn.cpp
 * @brief Implementazione OccultationSearchAstDyn
 */

 #include "ioccultcalc/occultation_search_astdyn.h"
 #include "ioccultcalc/astdyn_interface.h"
 #include "ioccultcalc/time_utils.h"
#include "ioccultcalc/types.h"
#include "ioccultcalc/orbital_elements.h"
 #include <chrono>
 #include <stdexcept>
 #include <cmath>
#include <iostream>
 
 #include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include "phase2_occultation_geometry.h"
#include "ioccultcalc/astdys_client.h"
 
 namespace ioccultcalc {
 
 constexpr double MJD_TO_JD = 2400000.5;
 
 class OccultationSearchAstDyn::Impl {
 public:
     std::unique_ptr<AstDynPropagator> propagator;
     AstDySElements astdys_elements;
     astdyn::propagation::KeplerianElements keplerian_elements;
     bool has_elements;
     Phase1CandidateScreening phase1;
     Phase2OccultationGeometry phase2;
     double asteroid_diameter_km;
     double orbital_uncertainty_km;
     ioc::gaia::UnifiedGaiaCatalog* gaia_catalog;
     
     Impl(double tolerance)
         : propagator(std::make_unique<AstDynPropagator>(tolerance))
         , has_elements(false)
         , asteroid_diameter_km(0.0)
         , orbital_uncertainty_km(100.0)
         , gaia_catalog(nullptr)
     {
         keplerian_elements.semi_major_axis = 0.0;
         keplerian_elements.eccentricity = 0.0;
         keplerian_elements.inclination = 0.0;
         keplerian_elements.longitude_ascending_node = 0.0;
         keplerian_elements.argument_perihelion = 0.0;
         keplerian_elements.mean_anomaly = 0.0;
         keplerian_elements.epoch_mjd_tdb = 0.0;
         keplerian_elements.gravitational_parameter = 
             1.32712440018e20 / std::pow(1.495978707e11, 3) * std::pow(86400.0, 2);
     }
 };
 
 OccultationSearchAstDyn::OccultationSearchAstDyn(double tolerance)
     : pimpl_(std::make_unique<Impl>(tolerance))
 {
 }
 
 OccultationSearchAstDyn::~OccultationSearchAstDyn() = default;
 
 astdyn::propagation::KeplerianElements 
 OccultationSearchAstDyn::convertToKeplerian(const AstDySElements& astdys) const {
     astdyn::propagation::KeplerianElements kep;
     kep.semi_major_axis = astdys.a;
     kep.eccentricity = astdys.e;
     kep.inclination = astdys.i * DEG_TO_RAD;
     kep.longitude_ascending_node = astdys.Omega * DEG_TO_RAD;
     kep.argument_perihelion = astdys.omega * DEG_TO_RAD;
     kep.mean_anomaly = astdys.M * DEG_TO_RAD;
     kep.epoch_mjd_tdb = astdys.epoch_mjd;
     kep.gravitational_parameter = 
         1.32712440018e20 / std::pow(1.495978707e11, 3) * std::pow(86400.0, 2);
     return kep;
 }
 
 void OccultationSearchAstDyn::syncOrbitalElements() {
     if (!pimpl_->has_elements) return;
     pimpl_->keplerian_elements = convertToKeplerian(pimpl_->astdys_elements);
     pimpl_->phase1.setOrbitalElements(pimpl_->keplerian_elements);
     pimpl_->phase2.setOrbitalElements(pimpl_->keplerian_elements);
 }
 
 bool OccultationSearchAstDyn::loadAsteroid(int asteroid_number) {
     try {
         pimpl_->astdys_elements = AstDySClient::downloadElements(asteroid_number);
         pimpl_->has_elements = true;
         syncOrbitalElements();
         return true;
     } catch (...) {
         return false;
     }
 }
 
 bool OccultationSearchAstDyn::loadAsteroid(const std::string& designation) {
     try {
         pimpl_->astdys_elements = AstDySClient::downloadElements(designation);
         pimpl_->has_elements = true;
         syncOrbitalElements();
         return true;
     } catch (...) {
         return false;
     }
 }
 
 bool OccultationSearchAstDyn::loadAsteroidFromEQ1(const std::string& eq1_path) {
     if (pimpl_->phase1.loadAsteroidFromEQ1(eq1_path)) {
         pimpl_->keplerian_elements = pimpl_->phase1.getOrbitalElements();
         pimpl_->astdys_elements.a = pimpl_->keplerian_elements.semi_major_axis;
         pimpl_->astdys_elements.e = pimpl_->keplerian_elements.eccentricity;
         pimpl_->astdys_elements.i = pimpl_->keplerian_elements.inclination * RAD_TO_DEG;
         pimpl_->astdys_elements.Omega = pimpl_->keplerian_elements.longitude_ascending_node * RAD_TO_DEG;
         pimpl_->astdys_elements.omega = pimpl_->keplerian_elements.argument_perihelion * RAD_TO_DEG;
         pimpl_->astdys_elements.M = pimpl_->keplerian_elements.mean_anomaly * RAD_TO_DEG;
         pimpl_->astdys_elements.epoch_mjd = pimpl_->keplerian_elements.epoch_mjd_tdb;
         pimpl_->has_elements = true;
         syncOrbitalElements();
         return true;
     }
     return false;
 }
 
 bool OccultationSearchAstDyn::loadAsteroidFromJSON(int asteroid_number, const std::string& json_path) {
     if (pimpl_->phase1.loadAsteroidFromJSON(asteroid_number, json_path)) {
         pimpl_->keplerian_elements = pimpl_->phase1.getOrbitalElements();
         pimpl_->astdys_elements.a = pimpl_->keplerian_elements.semi_major_axis;
         pimpl_->astdys_elements.e = pimpl_->keplerian_elements.eccentricity;
         pimpl_->astdys_elements.i = pimpl_->keplerian_elements.inclination * RAD_TO_DEG;
         pimpl_->astdys_elements.Omega = pimpl_->keplerian_elements.longitude_ascending_node * RAD_TO_DEG;
         pimpl_->astdys_elements.omega = pimpl_->keplerian_elements.argument_perihelion * RAD_TO_DEG;
         pimpl_->astdys_elements.M = pimpl_->keplerian_elements.mean_anomaly * RAD_TO_DEG;
         pimpl_->astdys_elements.epoch_mjd = pimpl_->keplerian_elements.epoch_mjd_tdb;
         pimpl_->astdys_elements.number = asteroid_number;
         pimpl_->has_elements = true;
         syncOrbitalElements();
         return true;
     }
     return false;
 }
 
 void OccultationSearchAstDyn::setAsteroidElements(const AstDySElements& elements) {
     pimpl_->astdys_elements = elements;
     pimpl_->has_elements = true;
     syncOrbitalElements();
 }
 
 OccultationSearchResults OccultationSearchAstDyn::search(
     const OccultationSearchConfig& config) {
     if (!pimpl_->has_elements) {
         throw std::runtime_error("Asteroide non caricato. Chiamare loadAsteroid() prima.");
     }
     
     auto t_start = std::chrono::high_resolution_clock::now();
     OccultationSearchResults results;
     
     if (pimpl_->gaia_catalog) {
         pimpl_->phase1.setCatalog(pimpl_->gaia_catalog);
     }
     
     results.asteroid_name = pimpl_->astdys_elements.name;
     results.asteroid_number = pimpl_->astdys_elements.number;
     
     Phase1Config phase1_cfg = config.phase1_config;
     phase1_cfg.start_mjd_tdb = config.start_mjd_tdb;
     phase1_cfg.end_mjd_tdb = config.end_mjd_tdb;
     phase1_cfg.max_magnitude = config.max_magnitude;
     
    std::cerr << "[OccultationSearch] Starting Phase1 (propagation + Gaia query)...\n";
    std::cerr.flush();
     auto phase1_start = std::chrono::high_resolution_clock::now();
     Phase1Results phase1_results = pimpl_->phase1.screenCandidates(phase1_cfg);
     auto phase1_end = std::chrono::high_resolution_clock::now();
    std::cerr << "[OccultationSearch] Phase1 completed: " 
              << phase1_results.num_candidates_filtered << " candidates found\n";
    std::cerr.flush();
     
     results.phase1_time_ms = std::chrono::duration<double, std::milli>(
         phase1_end - phase1_start).count();
     
     results.num_candidates_found = phase1_results.num_candidates_filtered;
     results.num_stars_in_corridor = phase1_results.num_stars_in_corridor;
     
     Phase2Config phase2_cfg = config.phase2_config;
     
     if (config.asteroid_diameter_km > 0) {
         pimpl_->asteroid_diameter_km = config.asteroid_diameter_km;
     }
     if (config.orbital_uncertainty_km > 0) {
         pimpl_->orbital_uncertainty_km = config.orbital_uncertainty_km;
     }
     
     auto phase2_start = std::chrono::high_resolution_clock::now();
     Phase2Results phase2_results = pimpl_->phase2.calculateGeometry(
         phase1_results.candidates, phase2_cfg);
     auto phase2_end = std::chrono::high_resolution_clock::now();
     
     results.phase2_time_ms = std::chrono::duration<double, std::milli>(
         phase2_end - phase2_start).count();
     
    // Converti gli eventi da Phase2Results a ioccultcalc::OccultationEvent
    results.events.clear();
    results.events.reserve(phase2_results.events.size());
    for (const auto& phase2_event : phase2_results.events) {
        ioccultcalc::OccultationEvent event;
        // Conversione completa dei campi da Phase2 a OccultationEvent
        event.eventId = std::to_string(phase2_event.star_source_id) + "_" + 
                       std::to_string(phase2_event.asteroid_number);
        event.timeCA = JulianDate(phase2_event.time_ca_mjd_utc + MJD_TO_JD);
        event.closeApproachDistance = phase2_event.closest_approach_mas / 1000.0; // mas -> arcsec
        event.positionAngle = phase2_event.position_angle_deg;
        event.maxDuration = phase2_event.max_duration_sec;
        event.probability = phase2_event.occultation_probability;
        
        // Converti asteroid da Phase2 (ha solo nome/numero) a EquinoctialElements
        // Usa elementi già caricati (AstDySElements contiene kepleriani)
        if (pimpl_->has_elements) {
            // Converti AstDySElements (kepleriani) a OrbitalElements, poi a EquinoctialElements
            OrbitalElements orbElem;
            orbElem.a = pimpl_->astdys_elements.a;
            orbElem.e = pimpl_->astdys_elements.e;
            orbElem.i = pimpl_->astdys_elements.i * DEG_TO_RAD;
            orbElem.Omega = pimpl_->astdys_elements.Omega * DEG_TO_RAD;
            orbElem.omega = pimpl_->astdys_elements.omega * DEG_TO_RAD;
            orbElem.M = pimpl_->astdys_elements.M * DEG_TO_RAD;
            orbElem.epoch.jd = pimpl_->astdys_elements.epoch_mjd + MJD_TO_JD;
            orbElem.designation = pimpl_->astdys_elements.name;
            orbElem.H = pimpl_->astdys_elements.H;
            orbElem.G = pimpl_->astdys_elements.G;
            
            // Converti OrbitalElements a EquinoctialElements
            EquinoctialElements astElem = EquinoctialElements::fromKeplerian(orbElem);
            astElem.name = phase2_event.asteroid_name;
            event.asteroid = astElem;
        }
        
        // Shadow path (già nel formato corretto)
        event.shadowPath.clear();
        for (const auto& pt : phase2_event.shadow_path) {
            ShadowPathPoint shadowPt;
            shadowPt.time.jd = pt.time_mjd_utc + MJD_TO_JD;
            shadowPt.location.latitude = pt.latitude_deg * DEG_TO_RAD;
            shadowPt.location.longitude = pt.longitude_deg * DEG_TO_RAD;
            shadowPt.location.altitude = 0.0;
            shadowPt.duration = phase2_event.max_duration_sec;
            shadowPt.centerlineDistance = 0.0;  // TODO: calcolare
            event.shadowPath.push_back(shadowPt);
        }
        
        // Incertezze (converti da arcsec a km)
        // Usa AU da types.h (AU = 149597870.7 km)
        double uncertainty_rad = phase2_event.total_uncertainty_arcsec / 3600.0 * DEG_TO_RAD;
        event.uncertaintyNorth = uncertainty_rad * phase2_event.asteroid_distance_au * AU;
        event.uncertaintySouth = event.uncertaintyNorth;
        
        results.events.push_back(event);
    }
    
     results.num_events_calculated = phase2_results.successful_calculations;
     results.num_events_failed = phase2_results.failed_calculations;
     
     auto t_end = std::chrono::high_resolution_clock::now();
     results.total_time_ms = std::chrono::duration<double, std::milli>(
         t_end - t_start).count();
     
     return results;
 }
 
 OccultationSearchResults OccultationSearchAstDyn::search(
     double start_mjd_tdb,
     double end_mjd_tdb,
     double max_magnitude) {
     OccultationSearchConfig config;
     config.start_mjd_tdb = start_mjd_tdb;
     config.end_mjd_tdb = end_mjd_tdb;
     config.max_magnitude = max_magnitude;
     return search(config);
 }
 
 void OccultationSearchAstDyn::setAsteroidDiameter(double diameter_km) {
     pimpl_->asteroid_diameter_km = diameter_km;
 }
 
 void OccultationSearchAstDyn::setOrbitalUncertainty(double sigma_km) {
     pimpl_->orbital_uncertainty_km = sigma_km;
 }
 
 void OccultationSearchAstDyn::setGaiaCatalog(ioc::gaia::UnifiedGaiaCatalog* catalog) {
     pimpl_->gaia_catalog = catalog;
     if (catalog) {
         pimpl_->phase1.setCatalog(catalog);
     }
 }
 
 void OccultationSearchAstDyn::addObserverSite(const std::string& name,
                                               double lat_deg, 
                                               double lon_deg, 
                                               double elev_m) {
     pimpl_->phase2.addObserverSite(name, lat_deg, lon_deg, elev_m);
 }
 
 const AstDySElements& OccultationSearchAstDyn::getAsteroidElements() const {
     if (!pimpl_->has_elements) {
         throw std::runtime_error("Asteroide non caricato");
     }
     return pimpl_->astdys_elements;
 }
 
 const astdyn::propagation::KeplerianElements& 
 OccultationSearchAstDyn::getKeplerianElements() const {
     if (!pimpl_->has_elements) {
         throw std::runtime_error("Asteroide non caricato");
     }
     return pimpl_->keplerian_elements;
 }
 
 bool OccultationSearchAstDyn::hasAsteroid() const {
     return pimpl_->has_elements;
 }
 
 std::string OccultationSearchAstDyn::getAsteroidName() const {
     if (!pimpl_->has_elements) {
         return "";
     }
     return pimpl_->astdys_elements.name;
 }
 
 int OccultationSearchAstDyn::getAsteroidNumber() const {
     if (!pimpl_->has_elements) {
         return 0;
     }
     return pimpl_->astdys_elements.number;
 }
 
 } // namespace ioccultcalc