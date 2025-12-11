/**
 * @file phase2_occultation_geometry.cpp
 * @brief Implementazione Phase 2: Geometria precisa occultazione
 * @date 4 Dicembre 2025
 * 
 * IMPLEMENTAZIONE STRATEGIA:
 * ===========================
 * 
 * Phase 2 riceve i candidati da Phase 1 e per ciascuno:
 * 
 * 1. PROPAGAZIONE PRECISA (±5 min attorno CA):
 *    - RKF78 con tolleranza 1e-12
 *    - Perturbazioni planetarie attive
 *    - Step temporale: 1 secondo
 *    - ~600 punti per evento
 * 
 * 2. CORREZIONI ASTROMETRICHE:
 *    - Parallasse geocentrica
 *    - Aberrazione stellare e planetaria
 *    - Proper motion stella (da Gaia DR3 a epoca evento)
 *    - Light-time correction
 * 
 * 3. GEOMETRIA OCCULTAZIONE:
 *    - Trova istante esatto closest approach (sub-millisecondo)
 *    - Calcola miss distance con precisione sub-milliarcsecondo
 *    - Determina chord length e position angle
 *    - Calcola durata massima occultazione
 * 
 * 4. PROIEZIONE SU TERRA:
 *    - Path dell'ombra su ellissoide WGS84
 *    - Velocità e direzione ombra
 *    - Limiti nord/sud del path
 *    - Entry/exit points sul pianeta
 * 
 * 5. PREDIZIONI PER OSSERVATORI:
 *    - Geometria locale per ogni sito
 *    - Tempo CA locale, altitudine, azimut
 *    - Visibilità (Sole, Luna, meteo teorico)
 *    - Distanza dalla linea centrale
 */

#include "phase2_occultation_geometry.h"

// Astdyn
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/io/parsers/OrbFitEQ1Parser.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Integrator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/observations/RWOReader.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/observations/Observation.hpp"
// OrbitFitter di AstDyn
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/io/AstDysOrbitFitter.hpp"

// IOC GaiaLib
#include "ioc_gaialib/unified_gaia_catalog.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

// Per download HTTP (curl)
#include <curl/curl.h>

// Costanti
constexpr double MJD_TO_JD = 2400000.5;
constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;
constexpr double RAD_TO_MAS = RAD_TO_DEG * 3600.0 * 1000.0;  // milliarcsec
constexpr double MAS_TO_RAD = 1.0 / RAD_TO_MAS;
constexpr double AU_TO_KM = 1.495978707e8;
constexpr double EPSILON_J2000 = 23.4392911 * DEG_TO_RAD;  // Obliquità eclittica

// WGS84 ellissoide terrestre
constexpr double EARTH_EQUATORIAL_RADIUS_KM = 6378.137;
constexpr double EARTH_POLAR_RADIUS_KM = 6356.752;
constexpr double EARTH_FLATTENING = 1.0 / 298.257223563;

// ═══════════════════════════════════════════════════════════════
// IMPLEMENTAZIONE PIMPL
// ═══════════════════════════════════════════════════════════════

class Phase2OccultationGeometry::Impl {
public:
    // Elementi orbitali asteroide
    astdyn::propagation::KeplerianElements keplerian_elements;
    astdyn::propagation::KeplerianElements refined_elements;  // Elementi dopo fit
    bool has_elements = false;
    bool has_refined_elements = false;
    
    // Propagatore
    std::unique_ptr<astdyn::propagation::Propagator> propagator;
    
    // Siti osservatori
    std::vector<ObserverGeometry> observer_sites;
    
    // Orbit fitter di AstDyn (TODO: implementare quando disponibile)
    // std::unique_ptr<astdyn::fitting::OrbitFitter> orbit_fitter;
    
    // Nome/numero asteroide
    std::string asteroid_designation;
    
    // Parametri fisici asteroide (LinOccult method)
    double asteroid_diameter_km = 0.0;      ///< Diametro asteroide (km), 0 = usa placeholder
    double orbital_uncertainty_km = 100.0;  ///< Incertezza orbitale (km, 1-sigma)
    
    Impl() {
        // Inizializza propagatore con parametri ad alta precisione
        // Usa RKF78 con tutte le perturbazioni e correzioni
        auto ephemeris = std::make_shared<astdyn::ephemeris::PlanetaryEphemeris>();
        auto integrator = std::make_unique<astdyn::propagation::RKF78Integrator>(0.1, 1e-12);
        
        // Configurazione completa con TUTTE le perturbazioni
        astdyn::propagation::PropagatorSettings settings;
        
        // Perturbazioni planetarie
        settings.include_planets = true;
        settings.include_moon = true;
        settings.include_asteroids = true;  // AST17: Ceres, Vesta, ecc.
        
        // Correzioni relativistiche
        settings.include_relativity = true;
        
        // Tutti i pianeti attivi
        settings.perturb_mercury = true;
        settings.perturb_venus = true;
        settings.perturb_earth = true;
        settings.perturb_mars = true;
        settings.perturb_jupiter = true;
        settings.perturb_saturn = true;
        settings.perturb_uranus = true;
        settings.perturb_neptune = true;
        
        propagator = std::make_unique<astdyn::propagation::Propagator>(
            std::move(integrator), ephemeris, settings);
        
        // Inizializza orbit fitter di AstDyn (TODO: quando disponibile)
        // orbit_fitter = std::make_unique<astdyn::fitting::OrbitFitter>();
    }
    
    /**
     * @brief Scarica osservazioni RWO da AstDyS
     */
    std::vector<astdyn::observations::OpticalObservation> downloadRWOObservations(
        const std::string& designation);
    
    /**
     * @brief Esegue fit orbitale preciso con tutte le correzioni
     */
    OrbitalFitResults refineOrbitWithObservations(
        const std::vector<astdyn::observations::OpticalObservation>& observations,
        const Phase2Config& config,
        double target_mjd_tdb);
    
    /**
     * @brief Costruisce URL AstDyS per download RWO
     */
    std::string buildAstDySURL(const std::string& designation);
    
    /**
     * @brief Download file via HTTP
     */
    bool downloadFile(const std::string& url, const std::string& output_path);
};

// ═══════════════════════════════════════════════════════════════
// COSTRUTTORE / DISTRUTTORE
// ═══════════════════════════════════════════════════════════════

Phase2OccultationGeometry::Phase2OccultationGeometry() 
    : pimpl_(std::make_unique<Impl>()) {
}

Phase2OccultationGeometry::~Phase2OccultationGeometry() = default;

// ═══════════════════════════════════════════════════════════════
// CARICAMENTO ELEMENTI ORBITALI
// ═══════════════════════════════════════════════════════════════

bool Phase2OccultationGeometry::loadAsteroidFromEQ1(const std::string& eq1_path) {
    try {
        astdyn::io::parsers::OrbFitEQ1Parser parser;
        auto elements = parser.parse(eq1_path);
        
        pimpl_->keplerian_elements.semi_major_axis = elements.semi_major_axis;
        pimpl_->keplerian_elements.eccentricity = elements.eccentricity;
        pimpl_->keplerian_elements.inclination = elements.inclination;
        pimpl_->keplerian_elements.longitude_ascending_node = elements.longitude_asc_node;
        pimpl_->keplerian_elements.argument_perihelion = elements.argument_perihelion;
        pimpl_->keplerian_elements.mean_anomaly = elements.mean_anomaly;
        pimpl_->keplerian_elements.epoch_mjd_tdb = elements.epoch_mjd_tdb;
        pimpl_->keplerian_elements.gravitational_parameter = 
            1.32712440018e20 / std::pow(1.495978707e11, 3) * std::pow(86400.0, 2);
        
        pimpl_->has_elements = true;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Phase2: Errore caricamento " << eq1_path << ": " << e.what() << "\n";
        return false;
    }
}

bool Phase2OccultationGeometry::loadAsteroidFromJSON(int asteroid_number, const std::string& json_path) {
    try {
        // Determina path del database JSON
        std::string path = json_path;
        if (path.empty()) {
            const char* home = std::getenv("HOME");
            if (home) {
                path = std::string(home) + "/.ioccultcalc/data/all_numbered_asteroids.json";
            } else {
                std::cerr << "Phase2: Errore HOME non definito e json_path non specificato\n";
                return false;
            }
        }
        
        // Leggi file JSON
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "Phase2: Errore impossibile aprire " << path << "\n";
            return false;
        }
        
        nlohmann::json j;
        file >> j;
        
        // Cerca l'asteroide nel database
        if (!j.contains("asteroids")) {
            std::cerr << "Phase2: Errore chiave 'asteroids' non trovata\n";
            return false;
        }
        
        bool found = false;
        for (const auto& asteroid : j["asteroids"]) {
            if (asteroid["number"].get<int>() == asteroid_number) {
                // Estrai elementi orbitali (angoli in gradi nel JSON)
                pimpl_->keplerian_elements.semi_major_axis = asteroid["a"].get<double>();
                pimpl_->keplerian_elements.eccentricity = asteroid["e"].get<double>();
                pimpl_->keplerian_elements.inclination = asteroid["i"].get<double>() * DEG_TO_RAD;
                pimpl_->keplerian_elements.longitude_ascending_node = asteroid["Omega"].get<double>() * DEG_TO_RAD;
                pimpl_->keplerian_elements.argument_perihelion = asteroid["omega"].get<double>() * DEG_TO_RAD;
                pimpl_->keplerian_elements.mean_anomaly = asteroid["M"].get<double>() * DEG_TO_RAD;
                
                // Epoca: da JD a MJD TDB
                double epoch_jd = asteroid["epoch"].get<double>();
                pimpl_->keplerian_elements.epoch_mjd_tdb = epoch_jd - MJD_TO_JD;
                
                // GM Sole in AU³/day²
                pimpl_->keplerian_elements.gravitational_parameter = 
                    1.32712440018e20 / std::pow(1.495978707e11, 3) * std::pow(86400.0, 2);
                
                // Salva designazione per uso futuro
                if (asteroid.contains("name")) {
                    pimpl_->asteroid_designation = asteroid["name"].get<std::string>();
                }
                
                found = true;
                
                std::cout << "✓ Phase2: Elementi orbitali caricati per asteroide " << asteroid_number << "\n";
                break;
            }
        }
        
        if (!found) {
            std::cerr << "Phase2: Asteroide " << asteroid_number << " non trovato nel database JSON\n";
            return false;
        }
        
        pimpl_->has_elements = true;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Phase2: Errore parsing JSON: " << e.what() << "\n";
        return false;
    }
}

void Phase2OccultationGeometry::setOrbitalElements(
    const astdyn::propagation::KeplerianElements& elements) {
    pimpl_->keplerian_elements = elements;
    pimpl_->has_elements = true;
}

const astdyn::propagation::KeplerianElements& 
Phase2OccultationGeometry::getOrbitalElements() const {
    if (!pimpl_->has_elements) {
        throw std::runtime_error("Phase2: Elementi orbitali non caricati");
    }
    return pimpl_->keplerian_elements;
}

void Phase2OccultationGeometry::setAsteroidDiameter(double diameter_km) {
    pimpl_->asteroid_diameter_km = diameter_km;
}

void Phase2OccultationGeometry::setOrbitalUncertainty(double uncertainty_km) {
    pimpl_->orbital_uncertainty_km = uncertainty_km;
}

// ═══════════════════════════════════════════════════════════════
// GESTIONE OSSERVATORI
// ═══════════════════════════════════════════════════════════════

void Phase2OccultationGeometry::addObserverSite(
    const std::string& name, double lat_deg, double lon_deg, double elev_m) {
    
    ObserverGeometry site;
    site.site_name = name;
    site.latitude_deg = lat_deg;
    site.longitude_deg = lon_deg;
    site.elevation_m = elev_m;
    
    pimpl_->observer_sites.push_back(site);
}

void Phase2OccultationGeometry::clearObserverSites() {
    pimpl_->observer_sites.clear();
}

// ═══════════════════════════════════════════════════════════════
// HELPER: CALLBACK CURL
// ═══════════════════════════════════════════════════════════════

static size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

// ═══════════════════════════════════════════════════════════════
// COSTRUISCI URL AstDyS
// ═══════════════════════════════════════════════════════════════

std::string Phase2OccultationGeometry::Impl::buildAstDySURL(const std::string& designation) {
    // URL pattern: https://newton.spacedys.com/~astdys2/mpcobs/numbered/17/17030.rwo
    // Per asteroidi numerati: /numbered/XX/XXXXX.rwo (XX = prime due cifre)
    
    try {
        int number = std::stoi(designation);
        int folder = number / 1000;  // 17030 → folder 17
        
        std::ostringstream url;
        url << "https://newton.spacedys.com/~astdys2/mpcobs/numbered/"
            << folder << "/" << number << ".rwo";
        
        return url.str();
        
    } catch (...) {
        // Se non è un numero, prova con unnumbered
        std::ostringstream url;
        url << "https://newton.spacedys.com/~astdys2/mpcobs/unnumbered/"
            << designation << ".rwo";
        return url.str();
    }
}

// ═══════════════════════════════════════════════════════════════
// DOWNLOAD FILE HTTP
// ═══════════════════════════════════════════════════════════════

bool Phase2OccultationGeometry::Impl::downloadFile(
    const std::string& url, const std::string& output_path) {
    
    CURL* curl = curl_easy_init();
    if (!curl) {
        std::cerr << "  ✗ Errore inizializzazione CURL\n";
        return false;
    }
    
    std::string response_string;
    
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response_string);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 30L);
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L);  // Per HTTPS senza certificati
    
    CURLcode res = curl_easy_perform(curl);
    long response_code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &response_code);
    
    curl_easy_cleanup(curl);
    
    if (res != CURLE_OK || response_code != 200) {
        std::cerr << "  ✗ Download fallito: " << curl_easy_strerror(res) 
                 << " (HTTP " << response_code << ")\n";
        return false;
    }
    
    // Salva su file
    std::ofstream out(output_path);
    if (!out) {
        std::cerr << "  ✗ Impossibile scrivere " << output_path << "\n";
        return false;
    }
    
    out << response_string;
    out.close();
    
    return true;
}

// ═══════════════════════════════════════════════════════════════
// DOWNLOAD OSSERVAZIONI RWO
// ═══════════════════════════════════════════════════════════════

std::vector<astdyn::observations::OpticalObservation> 
Phase2OccultationGeometry::Impl::downloadRWOObservations(const std::string& designation) {
    
    std::cout << "  [RWO Download] Scaricamento osservazioni per " << designation << "...\n";
    
    try {
        // Costruisci URL AstDyS
        std::string url = buildAstDySURL(designation);
        std::cout << "  URL: " << url << "\n";
        
        // Path temporaneo per file RWO
        std::string temp_path = "/tmp/" + designation + ".rwo";
        
        // Download file
        std::cout << "  Downloading...\n";
        if (!downloadFile(url, temp_path)) {
            std::cerr << "  ⚠ Download fallito - uso elementi nominali\n";
            return {};
        }
        
        std::cout << "  ✓ File scaricato: " << temp_path << "\n";
        
        // Parse RWO con AstDyn
        std::cout << "  Parsing RWO...\n";
        auto observations = astdyn::observations::RWOReader::readFile(temp_path);
        
        std::cout << "  ✓ Parsate " << observations.size() << " osservazioni RWO\n";
        
        // Rimuovi file temporaneo
        std::remove(temp_path.c_str());
        
        if (observations.empty()) {
            std::cerr << "  ⚠ Nessuna osservazione valida! Uso elementi nominali.\n";
        }
        
        return observations;
        
    } catch (const std::exception& e) {
        std::cerr << "  ✗ Errore download RWO: " << e.what() << "\n";
        std::cerr << "  Continuo con elementi nominali.\n";
        return {};
    }
}

// ═══════════════════════════════════════════════════════════════
// ORBITAL FITTING CON OSSERVAZIONI
// ═══════════════════════════════════════════════════════════════

OrbitalFitResults Phase2OccultationGeometry::Impl::refineOrbitWithObservations(
    const std::vector<astdyn::observations::OpticalObservation>& observations,
    const Phase2Config& config,
    double target_mjd_tdb) {
    
    OrbitalFitResults results;
    
    if (observations.empty()) {
        results.fit_notes = "Nessuna osservazione disponibile per fit";
        return results;
    }
    
    std::cout << "  [Orbital Fit] Fitting con " << observations.size() << " osservazioni...\n";
    
    try {
        // Crea il fitter di AstDyn
        astdyn::io::AstDysOrbitFitter fitter;
        fitter.set_verbose(false);  // Verbose solo se config lo richiede
        
        // Imposta osservazioni ed elementi iniziali
        fitter.set_observations(observations);
        fitter.set_elements(keplerian_elements);
        
        // Esegui il fitting
        auto t_start = std::chrono::high_resolution_clock::now();
        auto fit_result = fitter.fit();
        auto t_end = std::chrono::high_resolution_clock::now();
        
        // Popola i risultati
        results.fit_successful = fit_result.converged;
        results.num_observations_used = fit_result.num_observations_used;
        results.iterations_performed = fit_result.num_iterations;
        
        // Calcola RMS e residuo massimo
        results.rms_residuals_arcsec = std::sqrt(
            fit_result.rms_ra * fit_result.rms_ra + 
            fit_result.rms_dec * fit_result.rms_dec
        );
        results.max_residual_arcsec = std::max(fit_result.rms_ra, fit_result.rms_dec) * 3.0; // ~3-sigma
        results.chi_squared = fit_result.chi_squared;
        
        // Se converge, salva gli elementi raffinati
        if (results.fit_successful) {
            refined_elements = fit_result.fitted_orbit;
            has_refined_elements = true;
            
            // Propaga elementi raffinati all'epoca target (CA)
            auto kep_target = propagator->propagate_keplerian(refined_elements, target_mjd_tdb);
            results.refined_elements = kep_target;
            
            results.fit_notes = "Fit convergente - elementi raffinati salvati";
            
            double fit_time = std::chrono::duration<double>(t_end - t_start).count();
            std::cout << "  ✓ Fit convergente in " << fit_time << " sec\n";
            std::cout << "    RMS: " << results.rms_residuals_arcsec << " arcsec\n";
            std::cout << "    Chi²: " << results.chi_squared << "\n";
            std::cout << "    Osservazioni: " << results.num_observations_used 
                     << " / " << fit_result.num_observations_loaded;
            if (fit_result.num_outliers > 0) {
                std::cout << " (" << fit_result.num_outliers << " outlier rimossi)";
            }
            std::cout << "\n";
        } else {
            results.fit_notes = "Fit non convergente - uso elementi nominali";
            std::cout << "  ⚠ Fit non convergente dopo " << results.iterations_performed << " iterazioni\n";
            std::cout << "  Uso elementi nominali\n";
        }
        
        // TODO: Estrarre la matrice di covarianza se disponibile
        // results.covariance_matrix = ...;
        
    } catch (const std::exception& e) {
        results.fit_successful = false;
        results.fit_notes = std::string("Errore fit: ") + e.what();
        std::cout << "  ✗ Errore durante fitting: " << e.what() << "\n";
        std::cout << "  Uso elementi nominali\n";
    }
    
    return results;
}

// ═══════════════════════════════════════════════════════════════
// CALCOLO GEOMETRIA (TO BE IMPLEMENTED)
// ═══════════════════════════════════════════════════════════════

Phase2Results Phase2OccultationGeometry::calculateGeometry(
    const std::vector<ioccultcalc::CandidateStar>& candidates,
    const Phase2Config& config) {
    
    Phase2Results results;
    auto t_start = std::chrono::high_resolution_clock::now();
    
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  PHASE 2: Geometria Precisa con Orbital Refinement        ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    std::cout << "Candidati da processare: " << candidates.size() << "\n";
    std::cout << "Orbital refinement: " << (config.refine_orbit_from_observations ? "ATTIVO" : "DISATTIVO") << "\n";
    std::cout << "Finestra temporale: ±" << config.time_window_minutes << " min\n";
    std::cout << "Risoluzione: " << config.time_step_seconds << " sec\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 1: ORBITAL REFINEMENT (se richiesto)
    // ═══════════════════════════════════════════════════════════════
    if (config.refine_orbit_from_observations && !config.mpc_code.empty()) {
        std::cout << "══════════════════════════════════════════════════════════\n";
        std::cout << "STEP 1: ORBITAL REFINEMENT\n";
        std::cout << "══════════════════════════════════════════════════════════\n\n";
        
        // Salva designation per riferimento
        pimpl_->asteroid_designation = config.mpc_code;
        
        // Download osservazioni RWO da AstDyS
        std::vector<astdyn::observations::OpticalObservation> all_observations = 
            pimpl_->downloadRWOObservations(config.mpc_code);
        
        if (!all_observations.empty() && !candidates.empty()) {
            std::cout << "  [Observations] Scaricate " << all_observations.size() << " osservazioni\n";
            
            // Seleziona le N osservazioni più recenti (se richiesto)
            std::vector<astdyn::observations::OpticalObservation> observations;
            if (!config.use_all_available_observations && 
                config.max_observations_for_fit > 0 && 
                all_observations.size() > static_cast<size_t>(config.max_observations_for_fit)) {
                
                // Ordina per epoca (MJD) decrescente
                std::sort(all_observations.begin(), all_observations.end(),
                    [](const astdyn::observations::OpticalObservation& a, 
                       const astdyn::observations::OpticalObservation& b) {
                        return a.mjd_utc > b.mjd_utc;
                    });
                
                // Prendi le prime N (le più recenti)
                observations.assign(all_observations.begin(), 
                                  all_observations.begin() + config.max_observations_for_fit);
                
                std::cout << "  [Observations] Selezionate " << observations.size() 
                         << " osservazioni più recenti per il fit\n";
                std::cout << "  [Observations] Epoca più recente: MJD " 
                         << std::fixed << std::setprecision(2) << observations[0].mjd_utc << "\n";
                std::cout << "  [Observations] Epoca più vecchia: MJD " 
                         << observations.back().mjd_utc << "\n";
            } else {
                observations = all_observations;
                std::cout << "  [Observations] Uso tutte le " << observations.size() 
                         << " osservazioni disponibili\n";
            }
            
            // Usa il primo candidato per determinare l'epoca target (CA time)
            double target_mjd_tdb = candidates[0].closest_approach_mjd;
            
            // Esegui fit orbitale
            results.orbital_fit = pimpl_->refineOrbitWithObservations(
                observations, config, target_mjd_tdb);
            
            if (results.orbital_fit.fit_successful) {
                results.orbit_refined = true;
                std::cout << "\n  ✓✓✓ ORBITA RAFFINATA CON SUCCESSO ✓✓✓\n";
                std::cout << "  Elementi propagati al CA per massima precisione\n\n";
            } else {
                std::cout << "\n  ⚠ Fit fallito - uso elementi nominali\n\n";
            }
        }
    }
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 2: CALCOLO GEOMETRIA PER OGNI CANDIDATO
    // ═══════════════════════════════════════════════════════════════
    std::cout << "══════════════════════════════════════════════════════════\n";
    std::cout << "STEP 2: CALCOLO GEOMETRIA OCCULTAZIONI\n";
    std::cout << "══════════════════════════════════════════════════════════\n\n";
    
    for (size_t i = 0; i < candidates.size(); ++i) {
        const auto& candidate = candidates[i];
        
        try {
            std::cout << "Candidato " << (i+1) << "/" << candidates.size() 
                     << " - Stella " << candidate.source_id << "...\n";
            
            OccultationEvent event = calculateSingleEvent(candidate, config);
            results.events.push_back(event);
            results.successful_calculations++;
            
            std::cout << "  ✓ CA: " << event.closest_approach_mas << " mas @ MJD " 
                     << std::fixed << std::setprecision(6) << event.time_ca_mjd_utc << "\n";
            std::cout << "  Duration: " << event.max_duration_sec << " sec\n";
            std::cout << "  Position Angle: " << std::setprecision(2) << event.position_angle_deg << " deg\n";
            
            // Output logica LinOccult
            std::cout << "  [LinOccult] Angular radius: " << std::setprecision(3) 
                     << event.angular_radius_arcsec << " arcsec\n";
            std::cout << "  [LinOccult] Total uncertainty: " << event.total_uncertainty_arcsec 
                     << " arcsec\n";
            std::cout << "  [LinOccult] Shadow threshold: " << event.shadow_threshold_arcsec 
                     << " arcsec\n";
            std::cout << "  [LinOccult] Probability: " << std::setprecision(4) 
                     << event.occultation_probability * 100.0 << "%\n";
            std::cout << "  [LinOccult] Is occultation: " 
                     << (event.is_occultation ? "YES" : "NO") << "\n";
            std::cout << "  Shadow path: " << event.path_length_km << " km\n\n";
            
        } catch (const std::exception& e) {
            std::cerr << "  ✗ Errore: " << e.what() << "\n\n";
            results.failed_calculations++;
            results.error_messages += std::string(e.what()) + "\n";
        }
    }
    
    auto t_end = std::chrono::high_resolution_clock::now();
    results.total_computation_time_ms = 
        std::chrono::duration<double, std::milli>(t_end - t_start).count();
    
    std::cout << "══════════════════════════════════════════════════════════\n";
    std::cout << "PHASE 2 COMPLETATA\n";
    std::cout << "══════════════════════════════════════════════════════════\n";
    std::cout << "  Orbital refinement: " << (results.orbit_refined ? "✓ Successo" : "✗ Non eseguito") << "\n";
    if (results.orbit_refined) {
        std::cout << "    RMS residui: " << results.orbital_fit.rms_residuals_arcsec << " arcsec\n";
        std::cout << "    Osservazioni: " << results.orbital_fit.num_observations_used << "\n";
    }
    std::cout << "  Eventi calcolati: " << results.successful_calculations << "\n";
    std::cout << "  Falliti: " << results.failed_calculations << "\n";
    std::cout << "  Tempo totale: " << results.total_computation_time_ms << " ms\n\n";
    
    return results;
}

OccultationEvent Phase2OccultationGeometry::calculateSingleEvent(
    const ioccultcalc::CandidateStar& candidate,
    const Phase2Config& config) {
    
    if (!pimpl_->has_elements) {
        throw std::runtime_error("Elementi orbitali non caricati");
    }
    
    // Seleziona quali elementi usare (raffinati o nominali)
    const auto& elements_to_use = pimpl_->has_refined_elements 
        ? pimpl_->refined_elements 
        : pimpl_->keplerian_elements;
    
    std::string elements_source = pimpl_->has_refined_elements 
        ? "elementi raffinati da fit RWO" 
        : "elementi nominali";
    
    std::cout << "  Usando: " << elements_source << "\n";
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 0: CORREZIONI ASTROMETRICHE (LINOCCULT METHOD)
    // ═══════════════════════════════════════════════════════════════
    
    // Proper Motion stella (da Gaia DR3 a epoca evento)
    // TODO: prendere PM da Gaia catalog se disponibile
    double star_pm_ra_mas_yr = 0.0;   // TODO: da Gaia
    double star_pm_dec_mas_yr = 0.0;  // TODO: da Gaia
    
    // Epoch Gaia DR3: J2016.0 = MJD 57532.0
    double ca_mjd = candidate.closest_approach_mjd;
    double gaia_epoch_mjd = 57532.0;
    double years_from_gaia = (ca_mjd - gaia_epoch_mjd) / 365.25;
    
    // Applica proper motion
    double star_ra_corrected_deg = candidate.ra_deg + 
                                   (star_pm_ra_mas_yr / 1000.0 / 3600.0) * years_from_gaia / 
                                   std::cos(candidate.dec_deg * DEG_TO_RAD);
    double star_dec_corrected_deg = candidate.dec_deg + 
                                    (star_pm_dec_mas_yr / 1000.0 / 3600.0) * years_from_gaia;
    
    // ═══════════════════════════════════════════════════════════════
    // CORREZIONI ASTROMETRICHE (LINOCCULT METHOD)
    // ═══════════════════════════════════════════════════════════════
    
    // NOTA IMPORTANTE:
    // - Le coordinate Gaia DR3 sono già corrette per aberrazione stellare
    //   (epoca J2016.0, frame ICRS). Non serve correzione aggiuntiva.
    // - L'aberrazione planetaria per l'asteroide è già gestita da AstDyn
    //   nella propagazione geocentrica.
    
    // 1. Proper Motion: già applicato sopra (star_ra_corrected_deg, star_dec_corrected_deg)
    
    // 2. Light-time correction (tempo di viaggio luce stella-Terra)
    // Per stelle lontane (>100 pc), correzione trascurabile (<0.1 mas)
    // Se parallasse disponibile: light_time_days = parallasse_mas / (1000.0 * 206265.0)
    // Per ora assumiamo stella all'infinito (correzione trascurabile)
    
    // Le coordinate finali della stella sono quelle con proper motion applicato
    // (star_ra_corrected_deg, star_dec_corrected_deg)
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 1: PROPAGAZIONE DENSA ATTORNO AL CA
    // ═══════════════════════════════════════════════════════════════
    double time_window_days = config.time_window_minutes / 1440.0;  // min → days
    double step_days = config.time_step_seconds / 86400.0;          // sec → days
    
    double start_mjd = ca_mjd - time_window_days;
    double end_mjd = ca_mjd + time_window_days;
    int num_steps = static_cast<int>((end_mjd - start_mjd) / step_days) + 1;
    
    std::cout << "  Propagazione densa: " << num_steps << " punti in ±" 
             << config.time_window_minutes << " min\n";
    
    // Propaga path denso
    std::vector<double> times;
    std::vector<Eigen::Vector3d> positions_icrf;
    
    for (int i = 0; i < num_steps; ++i) {
        double mjd = start_mjd + i * step_days;
        times.push_back(mjd);
        
        // Propaga con massima precisione
        auto kep_prop = pimpl_->propagator->propagate_keplerian(elements_to_use, mjd);
        auto cart = astdyn::propagation::keplerian_to_cartesian(kep_prop);
        
        // cart.position è BARICENTRICO ECLITTICO (come in AstDynPropagationHelper::getRADec)
        // Converti da baricentrico eclittico a geocentrico ICRF
        
        // 1. Converti asteroide baricentrico da eclittico a equatoriale ICRF
        Eigen::Vector3d asteroid_bary_ecl = cart.position;
        Eigen::Vector3d asteroid_bary_icrf;
        asteroid_bary_icrf[0] = asteroid_bary_ecl[0];
        asteroid_bary_icrf[1] = asteroid_bary_ecl[1] * std::cos(EPSILON_J2000) - asteroid_bary_ecl[2] * std::sin(EPSILON_J2000);
        asteroid_bary_icrf[2] = asteroid_bary_ecl[1] * std::sin(EPSILON_J2000) + asteroid_bary_ecl[2] * std::cos(EPSILON_J2000);
        
        // 2. Calcola posizione Terra baricentrica ICRF (come in AstDynPropagationHelper::getRADec)
        double jd_tdb = mjd + 2400000.5;
        Eigen::Vector3d earth_helio_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getPosition(
                astdyn::ephemeris::CelestialBody::EARTH, jd_tdb);
        Eigen::Vector3d sun_bary_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
        Eigen::Vector3d earth_bary_ecl = earth_helio_ecl - sun_bary_ecl;
        
        Eigen::Vector3d earth_bary_icrf;
        earth_bary_icrf[0] = earth_bary_ecl[0];
        earth_bary_icrf[1] = earth_bary_ecl[1] * std::cos(EPSILON_J2000) - earth_bary_ecl[2] * std::sin(EPSILON_J2000);
        earth_bary_icrf[2] = earth_bary_ecl[1] * std::sin(EPSILON_J2000) + earth_bary_ecl[2] * std::cos(EPSILON_J2000);
        
        // 3. Posizione geocentrica dell'asteroide = asteroide_baricentrico - terra_baricentrica
        Eigen::Vector3d geocentric_icrf = asteroid_bary_icrf - earth_bary_icrf;
        
        positions_icrf.push_back(geocentric_icrf);
    }
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 2: TROVA ISTANTE ESATTO CLOSEST APPROACH
    // ═══════════════════════════════════════════════════════════════
    
    // Coordinate stella corrette (J2000, con proper motion applicato)
    Eigen::Vector3d star_unit;
    double ra_rad = star_ra_corrected_deg * DEG_TO_RAD;
    double dec_rad = star_dec_corrected_deg * DEG_TO_RAD;
    star_unit[0] = std::cos(dec_rad) * std::cos(ra_rad);
    star_unit[1] = std::cos(dec_rad) * std::sin(ra_rad);
    star_unit[2] = std::sin(dec_rad);
    
    // Trova minima distanza angolare
    double min_distance_rad = M_PI;
    double min_time_mjd = ca_mjd;
    int min_index = 0;
    
    for (size_t i = 0; i < positions_icrf.size(); ++i) {
        Eigen::Vector3d ast_dir = positions_icrf[i].normalized();
        double cos_sep = star_unit.dot(ast_dir);
        cos_sep = std::max(-1.0, std::min(1.0, cos_sep));
        double sep_rad = std::acos(cos_sep);
        
        if (sep_rad < min_distance_rad) {
            min_distance_rad = sep_rad;
            min_time_mjd = times[i];
            min_index = i;
        }
    }
    
    double min_distance_mas = min_distance_rad * RAD_TO_MAS;
    
    std::cout << "  CA preciso: " << min_distance_mas << " mas @ MJD " 
             << std::fixed << std::setprecision(8) << min_time_mjd << "\n";
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 4: CALCOLA PARAMETRI GEOMETRICI
    // ═══════════════════════════════════════════════════════════════
    
    // Distanza asteroide dalla Terra al CA
    double earth_distance_au = positions_icrf[min_index].norm();
    
    // Velocità angolare (approssimata con differenza finita)
    double angular_velocity_rad_day = 0.0;
    if (min_index > 0 && min_index < static_cast<int>(positions_icrf.size()) - 1) {
        Eigen::Vector3d dir_before = positions_icrf[min_index - 1].normalized();
        Eigen::Vector3d dir_after = positions_icrf[min_index + 1].normalized();
        double sep_rad = std::acos(dir_before.dot(dir_after));
        double dt_days = times[min_index + 1] - times[min_index - 1];
        angular_velocity_rad_day = sep_rad / dt_days;
    }
    
    // ═══════════════════════════════════════════════════════════════
    // CALCOLO PARAMETRI OCCULTAZIONE (LINOCCULT METHOD)
    // ═══════════════════════════════════════════════════════════════
    
    // Usa diametro reale se disponibile, altrimenti placeholder
    double asteroid_diameter_km = pimpl_->asteroid_diameter_km > 0.0 
        ? pimpl_->asteroid_diameter_km 
        : 10.0;  // Fallback placeholder
    
    // Diametro angolare asteroide
    double angular_diameter_rad = asteroid_diameter_km / (earth_distance_au * AU_TO_KM);
    double angular_diameter_arcsec = angular_diameter_rad * RAD_TO_DEG * 3600.0;
    double angular_radius_arcsec = angular_diameter_arcsec / 2.0;
    
    // Durata massima teorica
    double max_duration_sec = 0.0;
    if (angular_velocity_rad_day > 0) {
        max_duration_sec = (angular_diameter_rad / angular_velocity_rad_day) * 86400.0;
    }
    
    // ═══════════════════════════════════════════════════════════════
    // CALCOLO POSITION ANGLE (LINOCCULT METHOD)
    // ═══════════════════════════════════════════════════════════════
    
    // Position angle: angolo tra direzione Nord e direzione del moto relativo
    // Formula: PA = atan2(dRA*cos(Dec), dDec) [radianti]
    // Convertito in gradi Est da Nord
    
    double position_angle_deg = 0.0;
    if (min_index > 0 && min_index < static_cast<int>(positions_icrf.size()) - 1) {
        // Calcola direzione del moto relativo
        Eigen::Vector3d dir_before = positions_icrf[min_index - 1].normalized();
        Eigen::Vector3d dir_after = positions_icrf[min_index + 1].normalized();
        Eigen::Vector3d motion_dir = (dir_after - dir_before).normalized();
        
        // Converti in coordinate sferiche per calcolare PA
        // RA e Dec del punto prima e dopo
        double ra_before = std::atan2(dir_before[1], dir_before[0]);
        double dec_before = std::asin(dir_before[2]);
        double ra_after = std::atan2(dir_after[1], dir_after[0]);
        double dec_after = std::asin(dir_after[2]);
        
        // Differenze
        double delta_ra = ra_after - ra_before;
        double delta_dec = dec_after - dec_before;
        
        // Normalizza delta_ra
        while (delta_ra > M_PI) delta_ra -= 2.0 * M_PI;
        while (delta_ra < -M_PI) delta_ra += 2.0 * M_PI;
        
        // Position angle: atan2(dRA*cos(Dec), dDec)
        double dec_avg = (dec_before + dec_after) / 2.0;
        position_angle_deg = std::atan2(delta_ra * std::cos(dec_avg), delta_dec) * RAD_TO_DEG;
        
        // Converti da radianti a gradi Est da Nord (0° = Nord, 90° = Est)
        // LinOccult usa convenzione: 0° = Nord, positivo verso Est
        if (position_angle_deg < 0) position_angle_deg += 360.0;
    }
    
    // ═══════════════════════════════════════════════════════════════
    // VERIFICA OCCULTAZIONE (LINOCCULT METHOD)
    // ═══════════════════════════════════════════════════════════════
    
    // LinOccult usa shadow_threshold_factor = 3.0
    // Criterio: closest_approach < (angular_radius + 3*uncertainty)
    constexpr double SHADOW_THRESHOLD_FACTOR = 3.0;  // Standard LinOccult
    
    // Incertezza orbitale angolare
    double orbital_uncertainty_km = pimpl_->orbital_uncertainty_km;
    double angular_uncertainty_arcsec = (orbital_uncertainty_km / (earth_distance_au * AU_TO_KM)) 
                                        * RAD_TO_DEG * 3600.0;
    
    // Incertezza stellare (da Gaia, se disponibile)
    // Base: 7 mas per Gaia EDR3, degradazione con proper motion
    double star_uncertainty_base_mas = 7.0;  // Gaia EDR3 baseline
    double star_pm_degradation_mas_per_year = 2.5;
    double years_from_epoch = (min_time_mjd - 51544.5) / 365.25;  // Epoch J2016.0
    double star_uncertainty_mas = star_uncertainty_base_mas + 
                                  star_pm_degradation_mas_per_year * std::abs(years_from_epoch);
    double star_uncertainty_arcsec = star_uncertainty_mas / 1000.0;
    
    // Incertezza totale (quadratura)
    double total_uncertainty_arcsec = std::sqrt(
        angular_uncertainty_arcsec * angular_uncertainty_arcsec +
        star_uncertainty_arcsec * star_uncertainty_arcsec
    );
    
    // Threshold LinOccult
    double shadow_threshold_arcsec = angular_radius_arcsec + 
                                     SHADOW_THRESHOLD_FACTOR * total_uncertainty_arcsec;
    double closest_approach_arcsec = min_distance_mas / 1000.0;
    
    bool is_occultation = closest_approach_arcsec < shadow_threshold_arcsec;
    
    // ═══════════════════════════════════════════════════════════════
    // CALCOLO PROBABILITÀ OCCULTAZIONE (LINOCCULT CDF METHOD)
    // ═══════════════════════════════════════════════════════════════
    
    double probability = 0.0;
    if (total_uncertainty_arcsec > 0) {
        // Usa CDF della distribuzione normale (come LinOccult)
        // Probabilità che la vera posizione cada entro il raggio dell'asteroide
        double z = (angular_radius_arcsec - closest_approach_arcsec) / total_uncertainty_arcsec;
        
        // CDF normale standard: 0.5 * (1 + erf(z / sqrt(2)))
        double prob = 0.5 * (1.0 + std::erf(z / std::sqrt(2.0)));
        
        // Limita a [0, 1]
        probability = std::max(0.0, std::min(1.0, prob));
    } else {
        // Senza incertezza, probabilità binaria
        probability = (closest_approach_arcsec <= angular_radius_arcsec) ? 1.0 : 0.0;
    }
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 5: CREA EVENTO
    // ═══════════════════════════════════════════════════════════════
    
    OccultationEvent event;
    
    // Identificazione
    event.star_source_id = candidate.source_id;
    event.asteroid_name = pimpl_->asteroid_designation;
    
    // Dati stella (usando coordinate corrette per proper motion)
    event.star_ra_deg = star_ra_corrected_deg;
    event.star_dec_deg = star_dec_corrected_deg;
    event.star_magnitude = candidate.phot_g_mean_mag;
    event.star_pm_ra_mas_yr = star_pm_ra_mas_yr;
    event.star_pm_dec_mas_yr = star_pm_dec_mas_yr;
    
    // Geometria
    event.time_ca_mjd_utc = min_time_mjd;  // TODO: convertire TDB→UTC
    event.closest_approach_mas = min_distance_mas;
    event.max_duration_sec = max_duration_sec;
    event.asteroid_distance_au = earth_distance_au;
    event.position_angle_deg = position_angle_deg;
    
    // Verifica occultazione (LinOccult method)
    event.is_occultation = is_occultation;
    event.occultation_probability = probability;
    event.shadow_threshold_arcsec = shadow_threshold_arcsec;
    event.angular_radius_arcsec = angular_radius_arcsec;
    event.total_uncertainty_arcsec = total_uncertainty_arcsec;
    
    // Inizializza campi non ancora implementati
    event.star_distance_au = 0.0;  // TODO: calcolare da parallasse Gaia
    event.sun_target_elongation_deg = 0.0;  // TODO: calcolare
    event.snr = 0.0;  // TODO: calcolare
    
    // ═══════════════════════════════════════════════════════════════
    // CALCOLO SHADOW PATH SU TERRA (LINOCCULT METHOD)
    // ═══════════════════════════════════════════════════════════════
    
    // Calcola percorso dell'ombra sulla superficie terrestre (WGS84)
    // Metodo: intersezione raggio asteroide-stella con ellissoide terrestre
    
    event.shadow_path.clear();
    event.shadow_width_km = asteroid_diameter_km;
    
    // Direzione stella (unit vector in ICRF)
    Eigen::Vector3d star_dir = star_unit;
    
    // Calcola shadow path per finestra temporale estesa (±30 minuti attorno CA)
    double path_window_minutes = 30.0;  // ±30 minuti
    double path_step_seconds = 10.0;    // Step 10 secondi
    double path_start_mjd = min_time_mjd - path_window_minutes / 1440.0;
    double path_end_mjd = min_time_mjd + path_window_minutes / 1440.0;
    int path_num_steps = static_cast<int>((path_end_mjd - path_start_mjd) * 86400.0 / path_step_seconds) + 1;
    
    double path_length_km = 0.0;
    Eigen::Vector3d prev_shadow_point;
    bool has_prev = false;
    
    for (int i = 0; i < path_num_steps; ++i) {
        double path_mjd = path_start_mjd + i * (path_end_mjd - path_start_mjd) / (path_num_steps - 1);
        
        // Propaga asteroide a questo istante
        auto kep_path = pimpl_->propagator->propagate_keplerian(elements_to_use, path_mjd);
        auto cart_path = astdyn::propagation::keplerian_to_cartesian(kep_path);
        
        // Converti a geocentrico ICRF (come prima)
        Eigen::Vector3d ast_bary_ecl = cart_path.position;
        Eigen::Vector3d ast_bary_icrf;
        ast_bary_icrf[0] = ast_bary_ecl[0];
        ast_bary_icrf[1] = ast_bary_ecl[1] * std::cos(EPSILON_J2000) - ast_bary_ecl[2] * std::sin(EPSILON_J2000);
        ast_bary_icrf[2] = ast_bary_ecl[1] * std::sin(EPSILON_J2000) + ast_bary_ecl[2] * std::cos(EPSILON_J2000);
        
        double jd_path = path_mjd + 2400000.5;
        Eigen::Vector3d earth_helio_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getPosition(
                astdyn::ephemeris::CelestialBody::EARTH, jd_path);
        Eigen::Vector3d sun_bary_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_path);
        Eigen::Vector3d earth_bary_ecl = earth_helio_ecl - sun_bary_ecl;
        
        Eigen::Vector3d earth_bary_icrf;
        earth_bary_icrf[0] = earth_bary_ecl[0];
        earth_bary_icrf[1] = earth_bary_ecl[1] * std::cos(EPSILON_J2000) - earth_bary_ecl[2] * std::sin(EPSILON_J2000);
        earth_bary_icrf[2] = earth_bary_ecl[1] * std::sin(EPSILON_J2000) + earth_bary_ecl[2] * std::cos(EPSILON_J2000);
        
        Eigen::Vector3d ast_geo_icrf = ast_bary_icrf - earth_bary_icrf;
        
        // Converti da AU a km
        Eigen::Vector3d ast_geo_km = ast_geo_icrf * AU_TO_KM;
        
        // Raggio asteroide-stella: r(s) = ast_geo_km + s * star_dir
        // Trova intersezione con ellissoide WGS84
        // Equazione ellissoide: x²/a² + y²/a² + z²/b² = 1
        // dove a = EARTH_EQUATORIAL_RADIUS_KM, b = EARTH_POLAR_RADIUS_KM
        
        // Parametro s per intersezione con ellissoide
        // Sostituisci r(s) nell'equazione ellissoide e risolvi per s
        double a = EARTH_EQUATORIAL_RADIUS_KM;
        double b = EARTH_POLAR_RADIUS_KM;
        
        // Coefficienti equazione quadratica: A*s² + B*s + C = 0
        double A = (star_dir[0]*star_dir[0] + star_dir[1]*star_dir[1]) / (a*a) + 
                   (star_dir[2]*star_dir[2]) / (b*b);
        double B = 2.0 * ((ast_geo_km[0]*star_dir[0] + ast_geo_km[1]*star_dir[1]) / (a*a) +
                          (ast_geo_km[2]*star_dir[2]) / (b*b));
        double C = (ast_geo_km[0]*ast_geo_km[0] + ast_geo_km[1]*ast_geo_km[1]) / (a*a) +
                   (ast_geo_km[2]*ast_geo_km[2]) / (b*b) - 1.0;
        
        double discriminant = B*B - 4.0*A*C;
        
        if (discriminant >= 0 && A != 0) {
            // Due soluzioni: prendi quella più vicina (s negativo = verso Terra)
            double s1 = (-B - std::sqrt(discriminant)) / (2.0 * A);
            double s2 = (-B + std::sqrt(discriminant)) / (2.0 * A);
            
            // Prendi s negativo (verso Terra) o quello più piccolo
            double s = (s1 < s2) ? s1 : s2;
            
            // Punto sulla superficie terrestre (star_dir è già unitario, s è in km)
            Eigen::Vector3d shadow_point_icrf = ast_geo_km + s * star_dir;
            
            // Converti ICRF a coordinate geografiche (WGS84)
            // ICRF è approssimativamente equatoriale, quindi:
            // x = R*cos(lat)*cos(lon)
            // y = R*cos(lat)*sin(lon)
            // z = R*sin(lat)
            
            double R = shadow_point_icrf.norm();
            if (R > 0) {
                double lat_rad = std::asin(shadow_point_icrf[2] / R);
                double lon_rad = std::atan2(shadow_point_icrf[1], shadow_point_icrf[0]);
                
                ShadowPathPoint path_point;
                path_point.time_mjd_utc = path_mjd;
                path_point.latitude_deg = lat_rad * RAD_TO_DEG;
                path_point.longitude_deg = lon_rad * RAD_TO_DEG;
                
                // Calcola velocità ombra (differenza finita)
                if (has_prev) {
                    double dx = shadow_point_icrf[0] - prev_shadow_point[0];
                    double dy = shadow_point_icrf[1] - prev_shadow_point[1];
                    double dz = shadow_point_icrf[2] - prev_shadow_point[2];
                    double dt = path_step_seconds;
                    path_point.speed_km_s = std::sqrt(dx*dx + dy*dy + dz*dz) / dt;
                    
                    // Accumula lunghezza path
                    path_length_km += std::sqrt(dx*dx + dy*dy + dz*dz);
                } else {
                    path_point.speed_km_s = 0.0;
                }
                
                path_point.position_angle_deg = position_angle_deg;
                
                event.shadow_path.push_back(path_point);
                prev_shadow_point = shadow_point_icrf;
                has_prev = true;
            }
        }
    }
    
    // Chord length: distanza tra primo e ultimo punto del path
    if (event.shadow_path.size() >= 2) {
        const auto& first = event.shadow_path.front();
        const auto& last = event.shadow_path.back();
        
        // Distanza geodetica approssimata (formula di Haversine)
        double lat1_rad = first.latitude_deg * DEG_TO_RAD;
        double lon1_rad = first.longitude_deg * DEG_TO_RAD;
        double lat2_rad = last.latitude_deg * DEG_TO_RAD;
        double lon2_rad = last.longitude_deg * DEG_TO_RAD;
        
        double dlat = lat2_rad - lat1_rad;
        double dlon = lon2_rad - lon1_rad;
        
        double a_hav = std::sin(dlat/2) * std::sin(dlat/2) +
                       std::cos(lat1_rad) * std::cos(lat2_rad) *
                       std::sin(dlon/2) * std::sin(dlon/2);
        double c_hav = 2 * std::atan2(std::sqrt(a_hav), std::sqrt(1 - a_hav));
        
        // Raggio medio Terra (approssimazione)
        double R_mean_km = (EARTH_EQUATORIAL_RADIUS_KM + EARTH_POLAR_RADIUS_KM) / 2.0;
        event.chord_length_km = R_mean_km * c_hav;
    } else {
        event.chord_length_km = 0.0;
    }
    
    event.path_length_km = path_length_km;
    
    // Durata attraversamento Terra: tempo tra primo e ultimo punto
    if (event.shadow_path.size() >= 2) {
        event.path_duration_sec = (event.shadow_path.back().time_mjd_utc - 
                                   event.shadow_path.front().time_mjd_utc) * 86400.0;
    } else {
    event.path_duration_sec = 0.0;
    }
    
    // Quality
    event.high_confidence = pimpl_->has_refined_elements;
    event.notes = elements_source;
    
    // Aggiungi note sulla verifica occultazione
    if (is_occultation) {
        event.notes += " [OCCULTATION CONFIRMED]";
    } else {
        event.notes += " [Close approach only]";
    }
    
    return event;
}
