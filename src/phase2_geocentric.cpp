/**
 * @file phase2_geocentric.cpp
 * @brief Phase 2 con calcolo GEOCENTRICO e RKF78 + Orbital Fitting
 * 
 * Usa RKF78 con perturbazioni planetarie per massima precisione.
 * Include orbital fitting con download da AstDyS + RWOReader + AstDysOrbitFitter.
 */

#include "phase2_geocentric.h"

// Header astdyn
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/io/parsers/OrbFitEQ1Parser.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Integrator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/io/AstDysOrbitFitter.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/observations/RWOReader.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <curl/curl.h>

namespace ioccult {

constexpr double MJD_TO_JD = 2400000.5;
constexpr double RAD_TO_DEG = 180.0 / M_PI;
constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double AU_TO_KM = 1.495978707e8;
constexpr double RAD_TO_ARCSEC = 206264.806247;

static Eigen::Vector3d eclipticToEquatorial(const Eigen::Vector3d& ecl) {
    constexpr double obliquity = 23.4392911 * DEG_TO_RAD;
    const double cos_obl = std::cos(obliquity);
    const double sin_obl = std::sin(obliquity);
    return Eigen::Vector3d(ecl.x(), 
                           ecl.y() * cos_obl - ecl.z() * sin_obl,
                           ecl.y() * sin_obl + ecl.z() * cos_obl);
}

static void cartesianToRaDec(const Eigen::Vector3d& pos, double& ra_rad, double& dec_rad) {
    double r = pos.norm();
    dec_rad = std::asin(pos.z() / r);
    ra_rad = std::atan2(pos.y(), pos.x());
    if (ra_rad < 0) ra_rad += 2 * M_PI;
}

static double angularSeparation(double ra1_deg, double dec1_deg, 
                                double ra2_deg, double dec_ast_deg) {
    double ra1 = ra1_deg * DEG_TO_RAD, dec1 = dec1_deg * DEG_TO_RAD;
    double ra2 = ra2_deg * DEG_TO_RAD, dec2 = dec_ast_deg * DEG_TO_RAD;
    double cos_sep = std::sin(dec1) * std::sin(dec2) +
                     std::cos(dec1) * std::cos(dec2) * std::cos(ra1 - ra2);
    return std::acos(std::max(-1.0, std::min(1.0, cos_sep))) * RAD_TO_DEG * 3600.0;
}

static double positionAngle(double ra_star_deg, double dec_star_deg,
                            double ra_ast_deg, double dec_ast_deg) {
    double ra1 = ra_star_deg * DEG_TO_RAD, dec1 = dec_star_deg * DEG_TO_RAD;
    double ra2 = ra_ast_deg * DEG_TO_RAD, dec2 = dec_ast_deg * DEG_TO_RAD;
    double dra = ra2 - ra1;
    double y = std::sin(dra) * std::cos(dec2);
    double x = std::cos(dec1) * std::sin(dec2) - std::sin(dec1) * std::cos(dec2) * std::cos(dra);
    double pa = std::atan2(y, x) * RAD_TO_DEG;
    if (pa < 0) pa += 360.0;
    return pa;
}

// Helper CURL per download da AstDyS
static size_t write_callback(void* contents, size_t size, size_t nmemb, std::string* output) {
    size_t totalSize = size * nmemb;
    output->append((char*)contents, totalSize);
    return totalSize;
}

static bool downloadFileFromURL(const std::string& url, std::string& output) {
    CURL* curl = curl_easy_init();
    if (!curl) return false;
    
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &output);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 30L);
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L); // Per https
    
    CURLcode res = curl_easy_perform(curl);
    long http_code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
    curl_easy_cleanup(curl);
    
    return (res == CURLE_OK && http_code == 200 && !output.empty());
}

static std::string buildAstDySEQ1Url(const std::string& designation) {
    int number = std::stoi(designation);
    int folder = (number / 1000);  // 17030 -> 17
    return "https://newton.spacedys.com/~astdys2/epoch/numbered/" + 
           std::to_string(folder) + "/" + designation + ".eq1";
}


static std::string buildAstDySRWOUrl(const std::string& designation) {
    int number = std::stoi(designation);
    int folder = (number / 1000);  // 17030 -> 17
    return "https://newton.spacedys.com/~astdys2/mpcobs/numbered/" + 
           std::to_string(folder) + "/" + designation + ".rwo";
}

struct Phase2Geocentric::Impl {
    astdyn::propagation::KeplerianElements keplerian_elements;
    astdyn::propagation::KeplerianElements refined_elements;
    std::unique_ptr<astdyn::propagation::Propagator> propagator;
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephemeris;
    double asteroid_diameter_km = 10.0;
    bool has_elements = false;
    bool has_refined_elements = false;
    
    Impl() {
        ephemeris = std::make_shared<astdyn::ephemeris::PlanetaryEphemeris>();
        
        // RKF78 ad alta precisione (tolleranza 1e-12)
        auto integrator = std::make_unique<astdyn::propagation::RKF78Integrator>(0.1, 1e-12);
        
        // Perturbazioni planetarie attive
        astdyn::propagation::PropagatorSettings settings;
        settings.include_planets = true;
        settings.include_moon = true;
        settings.include_asteroids = false;
        
        propagator = std::make_unique<astdyn::propagation::Propagator>(
            std::move(integrator), ephemeris, settings);
    }
    
    // CALCOLO GEOCENTRICO con RKF78
    bool computeGeocentricPosition(double mjd_tdb, double& ra_deg, double& dec_deg, double& distance_au) {
        // Usa elementi raffinati se disponibili, altrimenti nominali
        const auto& elements = has_refined_elements ? refined_elements : keplerian_elements;
        
        // Propaga con RKF78 (kepleriano -> cartesiano)
        auto kep_prop = propagator->propagate_keplerian(elements, mjd_tdb);
        auto cart_ecl = astdyn::propagation::keplerian_to_cartesian(kep_prop);
        Eigen::Vector3d ast_bary_icrf = eclipticToEquatorial(cart_ecl.position);
        
        double jd_tdb = mjd_tdb + MJD_TO_JD;
        Eigen::Vector3d earth_helio_ecl = astdyn::ephemeris::PlanetaryEphemeris::getPosition(
            astdyn::ephemeris::CelestialBody::EARTH, jd_tdb);
        Eigen::Vector3d sun_bary_ecl = astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
        Eigen::Vector3d earth_bary_icrf = eclipticToEquatorial(earth_helio_ecl - sun_bary_ecl);
        
        // GEOCENTRICO!
        Eigen::Vector3d ast_geo_icrf = ast_bary_icrf - earth_bary_icrf;
        
        double ra_rad, dec_rad;
        cartesianToRaDec(ast_geo_icrf, ra_rad, dec_rad);
        ra_deg = ra_rad * RAD_TO_DEG;
        dec_deg = dec_rad * RAD_TO_DEG;
        distance_au = ast_geo_icrf.norm();
        return true;
    }
    
    // Fit orbitale scaricando RWO da AstDyS e usando RWOReader + AstDysOrbitFitter
    bool fitOrbitFromAstDyS(const std::string& designation, const Phase2GeocentricConfig& config, Phase2GeocentricResults& results) {
        try {
            // STEP 1: Download elementi orbitali aggiornati da AstDyS
            std::cout << "📥 Downloading orbital elements (.eq1) from AstDyS...\n";
            std::string eq1_url = buildAstDySEQ1Url(designation);
            std::cout << "   URL: " << eq1_url << "\n";
            
            std::string eq1_content;
            if (!downloadFileFromURL(eq1_url, eq1_content)) {
                results.fit_notes = "Failed to download .eq1 file from AstDyS";
                std::cerr << "❌ EQ1 download failed\n";
                return false;
            }
            
            std::cout << "   Downloaded " << eq1_content.size() << " bytes\n";
            
            // Salva e parsalo
            std::string eq1_temp_file = "/tmp/asteroid_" + designation + ".eq1";
            std::ofstream eq1_ofs(eq1_temp_file);
            eq1_ofs << eq1_content;
            eq1_ofs.close();
            
            // Parsa elementi orbitali
            std::cout << "📊 Parsing orbital elements...\n";
            astdyn::io::parsers::OrbFitEQ1Parser eq1_parser;
            auto astdys_elements = eq1_parser.parse(eq1_temp_file);
            
            // Converti in elementi kepleriani
            astdyn::propagation::KeplerianElements initial_elements;
            initial_elements.semi_major_axis = astdys_elements.semi_major_axis;
            initial_elements.eccentricity = astdys_elements.eccentricity;
            initial_elements.inclination = astdys_elements.inclination;
            initial_elements.longitude_ascending_node = astdys_elements.longitude_asc_node;
            initial_elements.argument_perihelion = astdys_elements.argument_perihelion;
            initial_elements.mean_anomaly = astdys_elements.mean_anomaly;
            initial_elements.epoch_mjd_tdb = astdys_elements.epoch_mjd_tdb;
            initial_elements.gravitational_parameter = 
                1.32712440018e20 / std::pow(1.495978707e11, 3) * std::pow(86400.0, 2);
            
            std::cout << "   Epoch: MJD " << initial_elements.epoch_mjd_tdb << "\n";
            
            // STEP 2: Download osservazioni RWO
            std::cout << "🌐 Downloading RWO for asteroid " << designation << " from AstDyS...\n";
            
            // Costruisci URL e scarica file RWO
            std::string url = buildAstDySRWOUrl(designation);
            std::cout << "   URL: " << url << "\n";
            
            std::string rwo_content;
            if (!downloadFileFromURL(url, rwo_content)) {
                results.fit_notes = "Failed to download RWO file from AstDyS";
                std::cerr << "❌ Download failed\n";
                return false;
            }
            
            std::cout << "   Downloaded " << rwo_content.size() << " bytes\n";
            
            // Salva temporaneamente su disco
            std::string temp_file = "/tmp/asteroid_" + designation + ".rwo";
            std::ofstream ofs(temp_file);
            ofs << rwo_content;
            ofs.close();
            
            // Leggi osservazioni con RWOReader di AstDyn
            std::cout << "📊 Parsing observations with RWOReader...\n";
            auto observations = astdyn::observations::RWOReader::readFile(temp_file);
            results.total_observations_downloaded = observations.size();
            
            if (observations.empty()) {
                results.fit_notes = "No observations found in RWO file";
                std::cerr << "❌ No observations\n";
                return false;
            }
            
            std::cout << "   Found " << observations.size() << " observations\n";
            
            // Limita osservazioni se richiesto (usa le più recenti)
            std::vector<astdyn::observations::OpticalObservation> obs_to_use = observations;
            if (!config.use_all_available_observations && 
                config.max_observations_for_fit > 0 && 
                observations.size() > static_cast<size_t>(config.max_observations_for_fit)) {
                
                std::sort(obs_to_use.begin(), obs_to_use.end(), 
                    [](const auto& a, const auto& b) { return a.mjd_utc > b.mjd_utc; });
                obs_to_use.resize(config.max_observations_for_fit);
                std::cout << "   ⚡ Using " << obs_to_use.size() << " most recent observations (limit: " 
                         << config.max_observations_for_fit << ")\n";
            }
            
            
            // Configura fitter con elementi nominali e osservazioni
            astdyn::io::AstDysOrbitFitter fitter;
            fitter.set_elements(initial_elements);
            fitter.set_observations(obs_to_use);
            fitter.set_verbose(false);
            
            // Esegui fit
            auto fit_result = fitter.fit();
            
            if (fit_result.converged) {
                refined_elements = fit_result.fitted_orbit;
                has_refined_elements = true;
                results.observations_used_in_fit = fit_result.num_observations_used;
                results.fit_rms_arcsec = fit_result.rms_ra;  // RMS in RA
                results.fit_chi_squared = fit_result.chi_squared;
                results.fit_notes = "Fit converged successfully using " + 
                                   std::to_string(fit_result.num_observations_used) + " observations from AstDyS";
                std::cout << "✓ Fit successful!\n";
                std::cout << "   Observations used: " << fit_result.num_observations_used << "\n";
                std::cout << "   RMS RA:  " << fit_result.rms_ra << " arcsec\n";
                std::cout << "   RMS Dec: " << fit_result.rms_dec << " arcsec\n";
                std::cout << "   χ²: " << fit_result.chi_squared << "\n";
                return true;
            } else {
                results.fit_notes = "Fit did not converge after " + 
                                   std::to_string(fit_result.num_iterations) + " iterations";
                std::cout << "⚠ Fit did not converge\n";
                return false;
            }
        } catch (const std::exception& e) {
            results.fit_notes = std::string("Exception during fitting: ") + e.what();
            std::cerr << "❌ Error: " << e.what() << "\n";
            return false;
        }
    }
};

Phase2Geocentric::Phase2Geocentric() : pimpl_(std::make_unique<Impl>()) {}
Phase2Geocentric::~Phase2Geocentric() = default;

bool Phase2Geocentric::loadAsteroidFromEQ1(const std::string& eq1_path) {
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
        std::cerr << "Errore caricamento: " << e.what() << "\n";
        return false;
    }
}

void Phase2Geocentric::setOrbitalElements(const astdyn::propagation::KeplerianElements& elements) {
    pimpl_->keplerian_elements = elements;
    pimpl_->has_elements = true;
}

void Phase2Geocentric::setAsteroidDiameter(double diameter_km) {
    pimpl_->asteroid_diameter_km = diameter_km;
}

Phase2GeocentricResult Phase2Geocentric::refineCandidate(
    const ioccultcalc::CandidateStar& candidate, const Phase2GeocentricConfig& config) {
    
    Phase2GeocentricResult result;
    result.source_id = candidate.source_id;
    result.star_ra_deg = candidate.ra_deg;
    result.star_dec_deg = candidate.dec_deg;
    result.star_mag = candidate.phot_g_mean_mag;
    result.orbit_was_refined = pimpl_->has_refined_elements;
    
    if (!pimpl_->has_elements) {
        result.error_message = "Elementi orbitali non caricati";
        return result;
    }
    
    double star_ra = candidate.ra_deg;
    double star_dec = candidate.dec_deg;
    
    double mjd_start = candidate.closest_approach_mjd - config.time_window_minutes / 1440.0;
    double mjd_end = candidate.closest_approach_mjd + config.time_window_minutes / 1440.0;
    double step_days = config.time_step_seconds / 86400.0;
    
    double best_mjd = candidate.closest_approach_mjd;
    double best_sep = 1e9, best_ra = 0.0, best_dec = 0.0, best_dist = 0.0;
    
    for (double mjd = mjd_start; mjd <= mjd_end; mjd += step_days) {
        double ast_ra, ast_dec, dist_au;
        if (!pimpl_->computeGeocentricPosition(mjd, ast_ra, ast_dec, dist_au)) continue;
        double sep = angularSeparation(star_ra, star_dec, ast_ra, ast_dec);
        if (sep < best_sep) {
            best_sep = sep; best_mjd = mjd; best_ra = ast_ra; best_dec = ast_dec; best_dist = dist_au;
        }
    }
    
    // Golden section search
    double a = best_mjd - step_days, b = best_mjd + step_days;
    const double golden = 0.381966011250105;
    for (int iter = 0; iter < 20; ++iter) {
        double c = a + golden * (b - a), d = b - golden * (b - a);
        double ra_c, dec_c, dist_c, ra_d, dec_d, dist_d;
        pimpl_->computeGeocentricPosition(c, ra_c, dec_c, dist_c);
        pimpl_->computeGeocentricPosition(d, ra_d, dec_d, dist_d);
        double sep_c = angularSeparation(star_ra, star_dec, ra_c, dec_c);
        double sep_d = angularSeparation(star_ra, star_dec, ra_d, dec_d);
        if (sep_c < sep_d) {
            b = d;
            if (sep_c < best_sep) { best_sep = sep_c; best_mjd = c; best_ra = ra_c; best_dec = dec_c; best_dist = dist_c; }
        } else {
            a = c;
            if (sep_d < best_sep) { best_sep = sep_d; best_mjd = d; best_ra = ra_d; best_dec = dec_d; best_dist = dist_d; }
        }
        if (b - a < 1e-8) break;
    }
    
    result.refined_ca_mjd = best_mjd;
    result.refined_ca_arcsec = best_sep;
    result.asteroid_ra_deg = best_ra;
    result.asteroid_dec_deg = best_dec;
    result.asteroid_distance_au = best_dist;
    result.position_angle_deg = positionAngle(star_ra, star_dec, best_ra, best_dec);
    
    double diameter_rad = pimpl_->asteroid_diameter_km / (best_dist * AU_TO_KM);
    result.asteroid_angular_diameter_mas = diameter_rad * RAD_TO_ARCSEC * 1000.0;
    
    double dt = 1.0 / 86400.0;
    double ra1, dec1, dist1, ra2, dec2, dist2;
    pimpl_->computeGeocentricPosition(best_mjd - dt/2, ra1, dec1, dist1);
    pimpl_->computeGeocentricPosition(best_mjd + dt/2, ra2, dec2, dist2);
    result.angular_velocity_arcsec_per_sec = angularSeparation(ra1, dec1, ra2, dec2);
    
    double shadow_radius_arcsec = result.asteroid_angular_diameter_mas / 2000.0;
    result.miss_distance_arcsec = best_sep - shadow_radius_arcsec;
    result.is_occultation = (result.miss_distance_arcsec < 0);
    
    if (result.angular_velocity_arcsec_per_sec > 0) {
        double chord_length = result.asteroid_angular_diameter_mas / 1000.0;
        result.max_duration_sec = chord_length / result.angular_velocity_arcsec_per_sec;
    }
    result.success = true;
    return result;
}

Phase2GeocentricResults Phase2Geocentric::refineAllCandidates(
    const std::vector<ioccultcalc::CandidateStar>& candidates, const Phase2GeocentricConfig& config) {
    
    auto t_start = std::chrono::high_resolution_clock::now();
    Phase2GeocentricResults results;
    results.total_candidates = static_cast<int>(candidates.size());
    results.results.reserve(candidates.size());
    
    // Se richiesto, esegui orbital fitting PRIMA di processare le candidate
    if (config.fit_orbit_from_rwo && !config.asteroid_designation.empty()) {
        std::cout << "\n=== ORBITAL REFINEMENT ===\n";
        results.orbit_refinement_attempted = true;
        bool fit_success = pimpl_->fitOrbitFromAstDyS(config.asteroid_designation, config, results);
        results.orbit_refinement_successful = fit_success;
        
        if (fit_success) {
            std::cout << "✓ Using refined orbital elements for all candidates\n\n";
        } else {
            std::cout << "⚠ Fitting failed, using nominal elements\n\n";
        }
    }
    
    // Processa tutte le candidate (con elementi raffinati se disponibili)
    for (const auto& cand : candidates) {
        auto result = refineCandidate(cand, config);
        if (result.is_occultation) results.occultations_found++;
        results.results.push_back(std::move(result));
    }
    
    auto t_end = std::chrono::high_resolution_clock::now();
    results.processing_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    return results;
}

} // namespace ioccult
