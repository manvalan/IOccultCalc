/**
 * @file chebyshev_occultation_manager.cpp
 * @brief Implementazione ChebyshevOccultationManager
 * @author IOccultCalc Development Team
 * @date 4 Dicembre 2025
 */

#include "chebyshev_occultation_manager.h"
// Includes AstDyn per propagazione diretta
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/io/parsers/OrbFitEQ1Parser.hpp>
// Includes IOC_GaiaLib per ricerca stelle con UnifiedGaiaCatalog
#include "ioc_gaialib/unified_gaia_catalog.h"
#include "ioc_gaialib/types.h"
// Includes IOccultCalc
#include "ioccultcalc/coordinates.h"
#include "ioccultcalc/time_utils.h"
#include <chrono>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>

namespace ioccultcalc {

// Using declarations per tipi AstDyn
using astdyn::io::IOrbitParser;
using astdyn::propagation::Propagator;
using astdyn::propagation::PropagatorSettings;
using astdyn::propagation::RKF78Integrator;
using astdyn::propagation::KeplerianElements;
using astdyn::propagation::keplerian_to_cartesian;
using astdyn::ephemeris::PlanetaryEphemeris;
using astdyn::io::parsers::OrbFitEQ1Parser;

// ============================================================================
// OccultationSearchConfig Implementation
// ============================================================================

OccultationSearchConfig::OccultationSearchConfig()
    : start_mjd(61000.0)
    , end_mjd(61014.0)
    , num_propagation_points(100)
    , rkf78_tolerance(1e-12)
    , num_chebyshev_coeffs(8)
    , corridor_width_km(100.0)
    , search_step_days(0.1)
    , max_magnitude(16.0)
    , min_altitude_deg(15.0)
    , threshold_arcsec(2.0)
    , observer_lat_deg(40.0)
    , observer_lon_deg(0.0)
    , observer_alt_m(0.0)
{
}

OccultationSearchConfig OccultationSearchConfig::fastSurvey() {
    OccultationSearchConfig cfg;
    cfg.num_propagation_points = 50;
    cfg.num_chebyshev_coeffs = 6;
    cfg.corridor_width_km = 150.0;
    cfg.search_step_days = 0.25;
    cfg.max_magnitude = 14.0;
    return cfg;
}

OccultationSearchConfig OccultationSearchConfig::highPrecision() {
    OccultationSearchConfig cfg;
    cfg.num_propagation_points = 150;
    cfg.rkf78_tolerance = 1e-14;
    cfg.num_chebyshev_coeffs = 10;
    cfg.corridor_width_km = 80.0;
    cfg.search_step_days = 0.05;
    cfg.max_magnitude = 17.0;
    return cfg;
}

// ============================================================================
// ChebyshevOccultationManager::Impl (PIMPL)
// ============================================================================

class ChebyshevOccultationManager::Impl {
public:
    // Orbital elements caricati (formato AstDyn)
    astdyn::io::IOrbitParser::OrbitalElements elements;
    astdyn::propagation::KeplerianElements kep_elements;
    std::string asteroid_name;
    double asteroid_diameter_km = 0.0;
    
    // Propagator e ephemeris
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephemeris;
    
    // Posizioni propagate (cache) - barycentric ICRF J2000
    std::vector<Eigen::Vector3d> propagated_positions;
    std::vector<double> propagated_mjds;
    
    // Statistiche performance
    double propagation_time_ms = 0.0;
    double fitting_time_ms = 0.0;
    double corridor_query_time_ms = 0.0;
    double closest_approach_time_ms = 0.0;
    
    // Accuratezza Chebyshev
    double chebyshev_rms_error_au = 0.0;
    double chebyshev_max_error_au = 0.0;
    
    // Coefficienti Chebyshev (3 array per x, y, z)
    std::vector<double> chebyshev_coeffs_x;
    std::vector<double> chebyshev_coeffs_y;
    std::vector<double> chebyshev_coeffs_z;
    double chebyshev_t0 = 0.0;  // MJD inizio
    double chebyshev_t1 = 0.0;  // MJD fine
};

// ============================================================================
// ChebyshevOccultationManager Implementation
// ============================================================================

ChebyshevOccultationManager::ChebyshevOccultationManager(
    const OccultationSearchConfig& config)
    : pimpl_(std::make_unique<Impl>())
    , config_(config)
    , asteroid_loaded_(false)
    , propagation_done_(false)
    , fitting_done_(false)
    , corridor_searched_(false)
{
}

// Distruttore deve essere nel .cpp per gestire unique_ptr di tipi PIMPL
ChebyshevOccultationManager::~ChebyshevOccultationManager() {
    // Distruzione automatica di pimpl_
}

bool ChebyshevOccultationManager::loadAsteroidFromEQ1(const std::string& eq1_file) {
    try {
        std::cout << "Caricamento asteroide da: " << eq1_file << "\n";
        
        // Parsa file .eq1 con AstDyn parser
        astdyn::io::parsers::OrbFitEQ1Parser parser;
        pimpl_->elements = parser.parse(eq1_file);
        
        // Estrai nome asteroide dal path (formato: XXXXX_name.eq1)
        size_t slash_pos = eq1_file.find_last_of("/\\");
        size_t underscore_pos = eq1_file.find('_', slash_pos != std::string::npos ? slash_pos : 0);
        size_t dot_pos = eq1_file.find(".eq1");
        
        if (underscore_pos != std::string::npos && dot_pos != std::string::npos) {
            pimpl_->asteroid_name = eq1_file.substr(underscore_pos + 1, dot_pos - underscore_pos - 1);
        } else {
            pimpl_->asteroid_name = "Unknown";
        }
        
        // Converti a KeplerianElements per propagatore
        pimpl_->kep_elements.semi_major_axis = pimpl_->elements.semi_major_axis;
        pimpl_->kep_elements.eccentricity = pimpl_->elements.eccentricity;
        pimpl_->kep_elements.inclination = pimpl_->elements.inclination;
        pimpl_->kep_elements.longitude_ascending_node = pimpl_->elements.longitude_asc_node;
        pimpl_->kep_elements.argument_perihelion = pimpl_->elements.argument_perihelion;
        pimpl_->kep_elements.mean_anomaly = pimpl_->elements.mean_anomaly;
        pimpl_->kep_elements.epoch_mjd_tdb = pimpl_->elements.epoch_mjd_tdb;
        // GM Sole in unità AU³/day²
        pimpl_->kep_elements.gravitational_parameter = 1.32712440018e20 / 
            (1.495978707e11 * 1.495978707e11 * 1.495978707e11) * (86400.0 * 86400.0);
        
        // Crea ephemeris planetarie
        pimpl_->ephemeris = std::make_shared<astdyn::ephemeris::PlanetaryEphemeris>();
        
        std::cout << "  ✓ Asteroide: " << pimpl_->asteroid_name << "\n";
        std::cout << "  ✓ Epoca: " << std::fixed << std::setprecision(6) 
                  << pimpl_->elements.epoch_mjd_tdb << " MJD TDB\n";
        std::cout << "  ✓ a=" << pimpl_->elements.semi_major_axis << " AU, "
                  << "e=" << pimpl_->elements.eccentricity << "\n";
        std::cout << "  ✓ Elementi convertiti a formato Kepleriano\n";
        
        asteroid_loaded_ = true;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE nel caricamento asteroide: " << e.what() << "\n";
        return false;
    }
}

bool ChebyshevOccultationManager::propagateAndFit() {
    if (!asteroid_loaded_) {
        std::cerr << "ERRORE: Asteroide non caricato. Chiamare loadAsteroidFromEQ1() prima.\n";
        return false;
    }
    
    try {
        std::cout << "\n=== PROPAGAZIONE RKF78 (AstDyn Propagator) ===\n";
        std::cout << "Intervallo: MJD " << std::fixed << std::setprecision(1)
                  << config_.start_mjd << " - " << config_.end_mjd
                  << " (" << (config_.end_mjd - config_.start_mjd) << " giorni)\n";
        std::cout << "Punti di campionamento: " << config_.num_propagation_points << "\n";
        std::cout << "RKF78 tolerance: " << std::scientific << config_.rkf78_tolerance << " AU\n";
        
        // PASSO 1: Crea propagatore AstDyn con RKF78
        auto integrator = std::make_unique<astdyn::propagation::RKF78Integrator>(
            0.1, config_.rkf78_tolerance);
        
        // Configura tutte le perturbazioni
        astdyn::propagation::PropagatorSettings settings;
        settings.include_planets = true;
        settings.include_asteroids = true;
        settings.include_relativity = true;
        settings.perturb_mercury = true;
        settings.perturb_venus = true;
        settings.perturb_earth = true;
        settings.perturb_mars = true;
        settings.perturb_jupiter = true;
        settings.perturb_saturn = true;
        settings.perturb_uranus = true;
        settings.perturb_neptune = true;
        
        astdyn::propagation::Propagator propagator(
            std::move(integrator), pimpl_->ephemeris, settings);
        
        std::cout << "✓ Propagatore creato con tutte le perturbazioni (8 pianeti + relatività)\n";
        
        // PASSO 2: Propaga per tutti i punti richiesti
        auto start_time = std::chrono::high_resolution_clock::now();
        
        pimpl_->propagated_positions.clear();
        pimpl_->propagated_mjds.clear();
        
        constexpr double EPSILON_J2000 = 23.4392911 * M_PI / 180.0;
        double cos_eps = std::cos(EPSILON_J2000);
        double sin_eps = std::sin(EPSILON_J2000);
        
        for (size_t i = 0; i < config_.num_propagation_points; ++i) {
            double mjd = config_.start_mjd + 
                (i * (config_.end_mjd - config_.start_mjd) / (config_.num_propagation_points - 1));
            
            // Propaga a questa epoca
            auto kep_prop = propagator.propagate_keplerian(pimpl_->kep_elements, mjd);
            auto cart_ecl = astdyn::propagation::keplerian_to_cartesian(kep_prop);
            
            // Converti da eclittico J2000 a ICRF equatoriale (rotazione obliquità)
            Eigen::Vector3d pos_icrf;
            pos_icrf.x() = cart_ecl.position.x();
            pos_icrf.y() = cart_ecl.position.y() * cos_eps - cart_ecl.position.z() * sin_eps;
            pos_icrf.z() = cart_ecl.position.y() * sin_eps + cart_ecl.position.z() * cos_eps;
            
            pimpl_->propagated_positions.push_back(pos_icrf);
            pimpl_->propagated_mjds.push_back(mjd);
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        pimpl_->propagation_time_ms = std::chrono::duration<double, std::milli>(
            end_time - start_time).count();
        
        std::cout << "✓ Propagazione completata in " << std::fixed << std::setprecision(1)
                  << pimpl_->propagation_time_ms << " ms\n";
        std::cout << "  Tempo medio per punto: " << std::setprecision(2)
                  << (pimpl_->propagation_time_ms / config_.num_propagation_points)
                  << " ms\n";
        std::cout << "  Frame: Barycentric ICRF J2000.0\n";
        
        propagation_done_ = true;
        
        // PASSO 3: Fitta polinomi di Chebyshev (implementazione semplificata)
        std::cout << "\n=== FITTING CHEBYSHEV ===\n";
        std::cout << "Coefficienti per asse: " << config_.num_chebyshev_coeffs << "\n";
        std::cout << "Totale coefficienti: " << (config_.num_chebyshev_coeffs * 3) << "\n";
        
        start_time = std::chrono::high_resolution_clock::now();
        
        // Salva parametri temporali per fitting
        pimpl_->chebyshev_t0 = config_.start_mjd;
        pimpl_->chebyshev_t1 = config_.end_mjd;
        
        // TODO: Implementare fitting Chebyshev qui
        // Per ora usiamo interpolazione lineare come placeholder
        std::cout << "⚠ Fitting Chebyshev TODO - usando cache posizioni propagate\n";
        
        end_time = std::chrono::high_resolution_clock::now();
        pimpl_->fitting_time_ms = std::chrono::duration<double, std::milli>(
            end_time - start_time).count();
        
        std::cout << "✓ Preparazione completata in " << std::setprecision(1)
                  << pimpl_->fitting_time_ms << " ms\n";
        
        fitting_done_ = true;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE nella propagazione/fitting: " << e.what() << "\n";
        return false;
    }
}

Eigen::Vector3d ChebyshevOccultationManager::getPositionAtEpoch(double epoch_mjd) const {
    if (!propagation_done_) {
        throw std::runtime_error("Propagazione non eseguita. Chiamare propagateAndFit().");
    }
    
    // Interpolazione lineare dalla cache
    if (epoch_mjd < pimpl_->propagated_mjds.front() || 
        epoch_mjd > pimpl_->propagated_mjds.back()) {
        throw std::runtime_error("Epoca fuori range propagato");
    }
    
    // Trova indici per interpolazione
    size_t idx = 0;
    for (size_t i = 0; i < pimpl_->propagated_mjds.size() - 1; ++i) {
        if (epoch_mjd >= pimpl_->propagated_mjds[i] && 
            epoch_mjd <= pimpl_->propagated_mjds[i + 1]) {
            idx = i;
            break;
        }
    }
    
    // Interpolazione lineare
    double t0 = pimpl_->propagated_mjds[idx];
    double t1 = pimpl_->propagated_mjds[idx + 1];
    double alpha = (epoch_mjd - t0) / (t1 - t0);
    
    return pimpl_->propagated_positions[idx] * (1.0 - alpha) + 
           pimpl_->propagated_positions[idx + 1] * alpha;
}

// TODO: Definire CartesianStateICRF e implementare questo metodo
// CartesianStateICRF ChebyshevOccultationManager::getStateAtEpoch(double epoch_mjd) const {
//     if (!propagation_done_) {
//         throw std::runtime_error("Propagazione non eseguita. Chiamare propagateAndFit().");
//     }
//     
//     CartesianStateICRF state;
//     state.position = getPositionAtEpoch(epoch_mjd);
//     
//     // Velocità approssimata con differenze finite
//     double dt = 0.001; // 1.44 minuti
//     Eigen::Vector3d pos_before = getPositionAtEpoch(epoch_mjd - dt);
//     Eigen::Vector3d pos_after = getPositionAtEpoch(epoch_mjd + dt);
//     state.velocity = (pos_after - pos_before) / (2.0 * dt);
//     
//     return state;
// }

size_t ChebyshevOccultationManager::searchStarsInCorridor() {
    if (!fitting_done_) {
        std::cerr << "ERRORE: Fitting non eseguito. Chiamare propagateAndFit() prima.\n";
        return 0;
    }
    
    std::cout << "\n=== RICERCA STELLE NEL CORRIDOR ===\n";
    std::cout << "Larghezza corridor: " << config_.corridor_width_km << " km\n";
    std::cout << "Passo campionamento: " << config_.search_step_days << " giorni\n";
    std::cout << "Magnitudine max: " << config_.max_magnitude << "\n";
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    candidates_.clear();
    
    try {
        // 1. Genera path asteroide campionando con step temporale
        std::vector<ioc::gaia::CelestialPoint> corridor_path;
        double num_steps = (config_.end_mjd - config_.start_mjd) / config_.search_step_days;
        double total_distance = 0.0;
        
        for (size_t i = 0; i <= static_cast<size_t>(num_steps); ++i) {
            double epoch = config_.start_mjd + i * config_.search_step_days;
            
            // Ottieni posizione baricentrica ICRF dall'interpolazione cache
            Eigen::Vector3d pos_icrf = getPositionAtEpoch(epoch);
            
            // Converti da baricentrico ICRF a RA/Dec
            double ra_rad = std::atan2(pos_icrf.y(), pos_icrf.x());
            if (ra_rad < 0) ra_rad += 2.0 * M_PI;
            
            double r_xy = std::sqrt(pos_icrf.x()*pos_icrf.x() + pos_icrf.y()*pos_icrf.y());
            double dec_rad = std::atan2(pos_icrf.z(), r_xy);
            
            // Aggiungi punto al corridor path
            ioc::gaia::CelestialPoint point;
            point.ra = ra_rad * 180.0 / M_PI;   // gradi
            point.dec = dec_rad * 180.0 / M_PI;  // gradi
            corridor_path.push_back(point);
            
            // Accumula distanza per calcolo medio
            total_distance += pos_icrf.norm();
        }
        
        std::cout << "  Path asteroide: " << corridor_path.size() << " punti campionati\n";
        
        // 2. Converti larghezza corridor da km a gradi
        double avg_distance_au = total_distance / corridor_path.size();
        double avg_distance_km = avg_distance_au * 149597870.7;
        
        double corridor_width_au = config_.corridor_width_km / 149597870.7; // km to AU
        double corridor_width_rad = std::atan(corridor_width_au / avg_distance_au);
        double corridor_width_deg = corridor_width_rad * 180.0 / M_PI;
        
        std::cout << "  Distanza media asteroide: " << std::fixed << std::setprecision(3)
                  << avg_distance_au << " AU (" << std::setprecision(0) 
                  << avg_distance_km << " km)\n";
        std::cout << "  Larghezza corridor angolare: " << std::fixed << std::setprecision(6)
                  << corridor_width_deg << " gradi (" << std::setprecision(1)
                  << corridor_width_deg * 3600.0 << " arcsec)\n";
        
        // 3. Prepara parametri query corridor
        ioc::gaia::CorridorQueryParams params;
        params.path = corridor_path;
        params.width = corridor_width_deg;  // half-width in degrees
        params.max_magnitude = config_.max_magnitude;
        params.min_parallax = -1.0;  // no limit
        params.max_results = 0;      // no limit
        
        // 4. Inizializza UnifiedGaiaCatalog (STATIC method, must call BEFORE getInstance)
        std::cout << "  Inizializzazione UnifiedGaiaCatalog...\n";
        
        // Usa il catalogo locale in ~/.catalog/gaia_mag18_v2_multifile
        std::string home = std::getenv("HOME");
        std::string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
        std::string json_config = R"({
            "catalog_type": "multifile_v2",
            "multifile_directory": ")" + catalog_path + R"(",
            "cache_size_mb": 512,
            "prefetch_enabled": true,
            "log_level": "info"
        })";
        std::cout << "  Path catalogo: " << catalog_path << "\n";
        
        // IMPORTANTE: Inizializza PRIMA di chiamare getInstance()
        if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) {
            throw std::runtime_error("Impossibile inizializzare UnifiedGaiaCatalog. Verificare path: " + catalog_path);
        }
        std::cout << "  ✓ Catalogo inizializzato con successo\n";
        
        // Ora possiamo ottenere l'istanza
        auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
        
        // 5. Query corridor
        std::cout << "  Query corridor con " << params.path.size() << " punti...\n";
        std::vector<ioc::gaia::GaiaStar> stars = catalog.queryCorridor(params);
        
        std::cout << "  Stelle Gaia DR3 trovate: " << stars.size() << "\n";
        
        // 5. Converti stelle Gaia in StarCandidate
        for (const auto& star : stars) {
            StarCandidate candidate;
            candidate.catalog_id = "Gaia DR3 " + std::to_string(star.source_id);
            candidate.ra_deg = star.ra;
            candidate.dec_deg = star.dec;
            candidate.magnitude = star.phot_g_mean_mag;
            candidate.pmra_mas_yr = star.pmra;
            candidate.pmdec_mas_yr = star.pmdec;
            candidate.closest_distance_arcsec = 999999.0;  // Da calcolare
            candidate.closest_epoch_mjd = 0.0;
            candidate.position_angle_deg = 0.0;
            candidate.relative_velocity_km_s = 0.0;
            candidate.is_occultation = false;
            
            candidates_.push_back(candidate);
        }
        
    } catch (const std::exception& e) {
        std::string home = std::getenv("HOME");
        std::string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
        std::cerr << "ERRORE nella ricerca stelle: " << e.what() << "\n";
        std::cerr << "  Dettaglio: Probabilmente il catalogo Gaia non è disponibile\n";
        std::cerr << "  Path catalogo: " << catalog_path << "\n";
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    pimpl_->corridor_query_time_ms = std::chrono::duration<double, std::milli>(
        end_time - start_time).count();
    
    corridor_searched_ = true;
    
    std::cout << "✓ Ricerca completata in " << std::fixed << std::setprecision(1)
              << pimpl_->corridor_query_time_ms << " ms\n";
    std::cout << "  Stelle candidate: " << candidates_.size() << "\n";
    
    return candidates_.size();
}

size_t ChebyshevOccultationManager::computeClosestApproaches() {
    if (!corridor_searched_) {
        std::cerr << "ERRORE: Corridor non cercato. Chiamare searchStarsInCorridor() prima.\n";
        return 0;
    }
    
    std::cout << "\n=== CALCOLO CLOSEST APPROACHES ===\n";
    std::cout << "Stelle da analizzare: " << candidates_.size() << "\n";
    std::cout << "Threshold occultazione: " << config_.threshold_arcsec << " arcsec\n";
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    size_t num_occultations = 0;
    
    // Per ogni stella candidata, trova il closest approach
    for (auto& candidate : candidates_) {
        double min_distance = 999999.0;
        double best_epoch = config_.start_mjd;
        
        // Ricerca a grana fine nell'intervallo temporale
        // Step più piccolo per ricerca closest approach (1 ora = 0.04167 giorni)
        double fine_step = 0.04167;
        size_t num_fine_steps = static_cast<size_t>((config_.end_mjd - config_.start_mjd) / fine_step);
        
        for (size_t i = 0; i <= num_fine_steps; ++i) {
            double epoch = config_.start_mjd + i * fine_step;
            
            // Posizione asteroide baricentrica ICRF
            Eigen::Vector3d ast_pos_icrf = getPositionAtEpoch(epoch);
            
            // Converti a RA/Dec asteroide
            double ast_ra = std::atan2(ast_pos_icrf.y(), ast_pos_icrf.x()) * 180.0 / M_PI;
            if (ast_ra < 0) ast_ra += 360.0;
            
            double r_xy = std::sqrt(ast_pos_icrf.x()*ast_pos_icrf.x() + 
                                   ast_pos_icrf.y()*ast_pos_icrf.y());
            double ast_dec = std::atan2(ast_pos_icrf.z(), r_xy) * 180.0 / M_PI;
            
            // Calcola distanza angolare stella-asteroide (sfera celeste)
            // Formula haversine per great circle distance
            double delta_ra = (ast_ra - candidate.ra_deg) * M_PI / 180.0;
            double delta_dec = (ast_dec - candidate.dec_deg) * M_PI / 180.0;
            
            double a = std::sin(delta_dec/2.0) * std::sin(delta_dec/2.0) +
                      std::cos(candidate.dec_deg * M_PI / 180.0) * 
                      std::cos(ast_dec * M_PI / 180.0) *
                      std::sin(delta_ra/2.0) * std::sin(delta_ra/2.0);
            double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
            double distance_deg = c * 180.0 / M_PI;
            double distance_arcsec = distance_deg * 3600.0;
            
            // Traccia il minimo
            if (distance_arcsec < min_distance) {
                min_distance = distance_arcsec;
                best_epoch = epoch;
            }
        }
        
        // Aggiorna candidate con closest approach trovato
        candidate.closest_distance_arcsec = min_distance;
        candidate.closest_epoch_mjd = best_epoch;
        
        // Determina se è un'occultazione (entro threshold)
        if (min_distance < config_.threshold_arcsec) {
            candidate.is_occultation = true;
            num_occultations++;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    pimpl_->closest_approach_time_ms = std::chrono::duration<double, std::milli>(
        end_time - start_time).count();
    
    std::cout << "✓ Closest approach completato in " << std::fixed << std::setprecision(1)
              << pimpl_->closest_approach_time_ms << " ms\n";
    std::cout << "  Occultazioni trovate: " << num_occultations << "\n";
    
    // Stampa le occultazioni trovate
    if (num_occultations > 0) {
        std::cout << "\nOCCULTAZIONI RILEVATE:\n";
        std::cout << std::string(70, '=') << "\n";
        for (const auto& candidate : candidates_) {
            if (candidate.is_occultation) {
                std::cout << candidate.catalog_id << "\n";
                std::cout << "  Posizione: RA=" << std::fixed << std::setprecision(5)
                          << candidate.ra_deg << "°, Dec=" << candidate.dec_deg << "°\n";
                std::cout << "  Magnitudine: " << std::setprecision(2) << candidate.magnitude << "\n";
                std::cout << "  Closest approach: " << std::setprecision(3)
                          << candidate.closest_distance_arcsec << " arcsec\n";
                std::cout << "  Epoca: MJD " << std::setprecision(5)
                          << candidate.closest_epoch_mjd << "\n";
                std::cout << std::string(70, '-') << "\n";
            }
        }
    }
    
    return num_occultations;
}

OccultationSearchResults ChebyshevOccultationManager::performFullSearch() {
    OccultationSearchResults results;
    
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Step 1: Propaga e fitta
    if (!propagateAndFit()) {
        throw std::runtime_error("Propagazione/fitting fallita");
    }
    
    // Step 2: Cerca stelle nel corridor
    searchStarsInCorridor();
    
    // Step 3: Calcola closest approaches
    computeClosestApproaches();
    
    auto total_end = std::chrono::high_resolution_clock::now();
    
    // Popola risultati
    results.asteroid_name = pimpl_->asteroid_name;
    results.asteroid_diameter_km = pimpl_->asteroid_diameter_km;
    // TODO: Convertire OrbitalElements AstDyn a EquinoctialElements
    // results.elements = pimpl_->elements;
    results.candidates = candidates_;
    results.occultations = getOccultations();
    results.stars_in_corridor = candidates_.size();
    results.propagation_time_ms = pimpl_->propagation_time_ms;
    results.fitting_time_ms = pimpl_->fitting_time_ms;
    results.corridor_query_time_ms = pimpl_->corridor_query_time_ms;
    results.closest_approach_time_ms = pimpl_->closest_approach_time_ms;
    results.search_duration_seconds = std::chrono::duration<double>(
        total_end - total_start).count();
    results.chebyshev_rms_error_au = pimpl_->chebyshev_rms_error_au;
    results.chebyshev_max_error_au = pimpl_->chebyshev_max_error_au;
    
    return results;
}

std::vector<StarCandidate> ChebyshevOccultationManager::getOccultations() const {
    std::vector<StarCandidate> occultations;
    for (const auto& candidate : candidates_) {
        if (candidate.is_occultation) {
            occultations.push_back(candidate);
        }
    }
    return occultations;
}

std::string ChebyshevOccultationManager::getAsteroidName() const {
    return pimpl_->asteroid_name;
}

// TODO: Convertire OrbitalElements AstDyn a EquinoctialElements
// EquinoctialElements ChebyshevOccultationManager::getOrbitalElements() const {
//     return pimpl_->elements;
// }

double ChebyshevOccultationManager::getChebyshevRMSError() const {
    return pimpl_->chebyshev_rms_error_au;
}

void ChebyshevOccultationManager::printSummary() const {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "RIEPILOGO RICERCA OCCULTAZIONI\n";
    std::cout << std::string(70, '=') << "\n\n";
    
    std::cout << "Asteroide: " << pimpl_->asteroid_name << "\n";
    std::cout << "Intervallo: MJD " << std::fixed << std::setprecision(1)
              << config_.start_mjd << " - " << config_.end_mjd
              << " (" << (config_.end_mjd - config_.start_mjd) << " giorni)\n\n";
    
    std::cout << "Performance:\n";
    std::cout << "  Propagazione RKF78:     " << std::setw(8) << std::setprecision(1)
              << pimpl_->propagation_time_ms << " ms\n";
    std::cout << "  Fitting Chebyshev:      " << std::setw(8)
              << pimpl_->fitting_time_ms << " ms\n";
    std::cout << "  Query corridor:         " << std::setw(8)
              << pimpl_->corridor_query_time_ms << " ms\n";
    std::cout << "  Closest approaches:     " << std::setw(8)
              << pimpl_->closest_approach_time_ms << " ms\n\n";
    
    std::cout << "Accuratezza Chebyshev:\n";
    std::cout << "  RMS error: " << std::scientific << std::setprecision(3)
              << pimpl_->chebyshev_rms_error_au << " AU\n";
    std::cout << "  Max error: " << pimpl_->chebyshev_max_error_au << " AU\n\n";
    
    std::cout << "Risultati:\n";
    std::cout << "  Stelle nel corridor:    " << candidates_.size() << "\n";
    std::cout << "  Occultazioni trovate:   " << getOccultations().size() << "\n";
    
    std::cout << "\n" << std::string(70, '=') << "\n";
}

bool ChebyshevOccultationManager::exportResultsJSON(const std::string& output_file) const {
    // TODO: Implementare export JSON
    return false;
}

bool ChebyshevOccultationManager::exportResultsOOP(const std::string& output_file) const {
    // TODO: Implementare export OOP
    return false;
}

} // namespace ioccultcalc
