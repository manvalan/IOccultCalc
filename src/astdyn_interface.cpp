/**
 * @file astdyn_interface.cpp
 * @brief Implementazione AstDynPropagator e AstDynOrbitFitter
 */

#include "ioccultcalc/astdyn_interface.h"
#include "ioccultcalc/types.h"
#include "ioccultcalc/time_utils.h"
#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/astdys_client.h"  // ← AGGIUNGI
#include "ioccultcalc/astdyn_propagation_helper.h"  // ← AGGIUNGI per convertFromEquinoctial

// AstDyn headers
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Integrator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/core/Constants.hpp"
#include <Eigen/Dense>  // ← AGGIUNGI per Vector3d

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace ioccultcalc {

// Costanti
constexpr double MJD_TO_JD = 2400000.5;
constexpr double AU_TO_KM = 149597870.7;
constexpr double KM_TO_AU = 1.0 / AU_TO_KM;
constexpr double DEG_TO_RAD_LOCAL = M_PI / 180.0;
constexpr double RAD_TO_DEG_LOCAL = 180.0 / M_PI;

// ============================================================================
// AstDynPropagator Implementation
// ============================================================================

class AstDynPropagator::Impl {
public:
    std::unique_ptr<astdyn::propagation::Propagator> propagator;
    std::shared_ptr<astdyn::ephemeris::PlanetaryEphemeris> ephemeris;
    double tolerance;
    double h_min_days;
    double h_max_days;
    bool planets_enabled;
    bool asteroids_enabled;
    bool relativity_enabled;
    
    // Statistiche ultima propagazione
    int steps_accepted;
    int steps_rejected;
    double min_step;
    double max_step;
    
    Impl(double tol)
        : tolerance(tol)
        , h_min_days(1e-6)
        , h_max_days(1.0)
        , planets_enabled(true)
        , asteroids_enabled(false)
        , relativity_enabled(true)
        , steps_accepted(0)
        , steps_rejected(0)
        , min_step(0.0)
        , max_step(0.0)
    {
        initializePropagator();
    }
    
    void initializePropagator() {
        // Crea ephemeris planetaria
        ephemeris = std::make_shared<astdyn::ephemeris::PlanetaryEphemeris>();
        
        // Crea integrator RKF78
        auto integrator = std::make_unique<astdyn::propagation::RKF78Integrator>(
            h_min_days, tolerance);
        
        // Configura settings con TUTTE le perturbazioni (come Phase2)
        astdyn::propagation::PropagatorSettings settings;
        settings.include_planets = planets_enabled;
        settings.include_moon = planets_enabled;
        settings.include_asteroids = asteroids_enabled;
        settings.include_relativity = relativity_enabled;
        
        // Attiva TUTTI i pianeti individualmente (come Phase2)
        if (planets_enabled) {
            settings.perturb_mercury = true;
            settings.perturb_venus = true;
            settings.perturb_earth = true;
            settings.perturb_mars = true;
            settings.perturb_jupiter = true;
            settings.perturb_saturn = true;
            settings.perturb_uranus = true;
            settings.perturb_neptune = true;
        }
        
        // Crea propagatore
        propagator = std::make_unique<astdyn::propagation::Propagator>(
            std::move(integrator), ephemeris, settings);
    }
    
    // Converti AstDySElements -> astdyn::propagation::KeplerianElements
    astdyn::propagation::KeplerianElements toKeplerian(const AstDySElements& astdys) const {
        astdyn::propagation::KeplerianElements kep;
        kep.semi_major_axis = astdys.a;
        kep.eccentricity = astdys.e;
        kep.inclination = astdys.i * DEG_TO_RAD_LOCAL;
        kep.longitude_ascending_node = astdys.Omega * DEG_TO_RAD_LOCAL;
        kep.argument_perihelion = astdys.omega * DEG_TO_RAD_LOCAL;
        kep.mean_anomaly = astdys.M * DEG_TO_RAD_LOCAL;
        kep.epoch_mjd_tdb = astdys.epoch_mjd;
        // Parametro gravitazionale del Sole (AU^3/day^2)
        kep.gravitational_parameter = 
            1.32712440018e20 / std::pow(1.495978707e11, 3) * std::pow(86400.0, 2);
        return kep;
    }
    
    // Converti astdyn::propagation::KeplerianElements -> AstDySElements
    AstDySElements fromKeplerian(const astdyn::propagation::KeplerianElements& kep) const {
        AstDySElements astdys;
        astdys.a = kep.semi_major_axis;
        astdys.e = kep.eccentricity;
        astdys.i = kep.inclination * RAD_TO_DEG_LOCAL;
        astdys.Omega = kep.longitude_ascending_node * RAD_TO_DEG_LOCAL;
        astdys.omega = kep.argument_perihelion * RAD_TO_DEG_LOCAL;
        astdys.M = kep.mean_anomaly * RAD_TO_DEG_LOCAL;
        astdys.epoch_mjd = kep.epoch_mjd_tdb;
        astdys.H = 0.0;
        astdys.G = 0.15;
        astdys.has_covariance = false;
        return astdys;
    }
    
    // Calcola RA/Dec da stato cartesiano
    std::pair<double, double> computeRADec(
        const astdyn::propagation::CartesianElements& cart,  // ← Cambia in CartesianElements
        double mjd_tdb) const {
        
        // ═══════════════════════════════════════════════════════════════
        // UNITA' DI MISURA: Tutte le posizioni sono in AU (Unità Astronomiche)
        // ═══════════════════════════════════════════════════════════════
        
        // ═══════════════════════════════════════════════════════════════
        // VERIFICA FRAME: keplerian_to_cartesian potrebbe restituire
        // coordinate ELIOCENTRICHE invece di baricentriche!
        // ═══════════════════════════════════════════════════════════════
        // NOTA: Il commento dice "baricentriche" ma potrebbe essere errato
        // Se gli elementi kepleriani sono eliocentrici, allora:
        //   asteroide_geocentrico = asteroide_eliocentrico - terra_eliocentrico
        // Se sono baricentrici, allora:
        //   asteroide_geocentrico = asteroide_baricentrico - terra_baricentrico
        // ═══════════════════════════════════════════════════════════════
        
        // DEBUG: Verifica per asteroide 17030 il 28 nov 2025
        constexpr double TARGET_MJD = 60311.0;
        bool is_debug = (std::abs(mjd_tdb - TARGET_MJD) < 0.1);
        
        Eigen::Vector3d asteroid_pos_ecl = cart.position;  // [AU] eclittico J2000
        
        if (is_debug) {
            std::cerr << "[DEBUG computeRADec] MJD " << mjd_tdb << "\n";
            std::cerr << "[DEBUG] cart.position (eclittico): [" 
                      << asteroid_pos_ecl[0] << ", " 
                      << asteroid_pos_ecl[1] << ", " 
                      << asteroid_pos_ecl[2] << "] AU\n";
            std::cerr << "[DEBUG] Norm asteroid_pos_ecl: " << asteroid_pos_ecl.norm() << " AU\n";
        }
        
        // PROVA: Assumiamo che siano ELIOCENTRICHE invece di baricentriche
        // Se il frame è eliocentrico, dobbiamo usare terra_helio_ecl direttamente
        // Se il frame è baricentrico, dobbiamo usare terra_bary_ecl
        
        // Per ora manteniamo l'assunzione baricentrica, ma aggiungiamo logging
        Eigen::Vector3d asteroid_bary_ecl = asteroid_pos_ecl;  // [AU] eclittico J2000
        
        // Converti da eclittico a equatoriale ICRF J2000
        constexpr double EPSILON_J2000 = 23.4392911 * DEG_TO_RAD_LOCAL;  // Obliquità eclittica J2000
        double cos_eps = std::cos(EPSILON_J2000);
        double sin_eps = std::sin(EPSILON_J2000);
        
        Eigen::Vector3d asteroid_bary_icrf;  // [AU] ICRF J2000
        asteroid_bary_icrf[0] = asteroid_bary_ecl[0];
        asteroid_bary_icrf[1] = asteroid_bary_ecl[1] * cos_eps - asteroid_bary_ecl[2] * sin_eps;
        asteroid_bary_icrf[2] = asteroid_bary_ecl[1] * sin_eps + asteroid_bary_ecl[2] * cos_eps;
        
        // Converti MJD TDB -> JD TDB
        // MJD = JD - 2400000.5, quindi JD = MJD + 2400000.5
        double jd_tdb = mjd_tdb + 2400000.5;  // [JD TDB]
        
        // Usa PlanetaryEphemeris di AstDyn per coerenza (come in phase1_candidate_screening.cpp)
        // UNITA': AU (Unità Astronomiche)
        // SISTEMA: Eclittico J2000, eliocentrico
        Eigen::Vector3d earth_helio_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getPosition(
                astdyn::ephemeris::CelestialBody::EARTH, jd_tdb);  // [AU] eclittico J2000
        
        // UNITA': AU (Unità Astronomiche)
        // SISTEMA: Eclittico J2000, baricentrico del Sistema Solare
        Eigen::Vector3d sun_bary_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);  // [AU] eclittico J2000
        
        // Posizione baricentrica della Terra (eclittico)
        // UNITA': AU (Unità Astronomiche)
        // Formula: Terra_baricentrica = Terra_eliocentrica - Sole_baricentrico
        Eigen::Vector3d earth_bary_ecl = earth_helio_ecl - sun_bary_ecl;  // [AU] eclittico J2000
        
        // Converti asteroide da eclittico a ICRF (sia eliocentrico che baricentrico)
        // Assumiamo che asteroid_pos_ecl possa essere sia eliocentrico che baricentrico
        Eigen::Vector3d asteroid_helio_ecl = asteroid_pos_ecl;  // Prova: assumiamo eliocentrico
        Eigen::Vector3d asteroid_helio_icrf;  // [AU] ICRF J2000
        asteroid_helio_icrf[0] = asteroid_helio_ecl[0];
        asteroid_helio_icrf[1] = asteroid_helio_ecl[1] * cos_eps - asteroid_helio_ecl[2] * sin_eps;
        asteroid_helio_icrf[2] = asteroid_helio_ecl[1] * sin_eps + asteroid_helio_ecl[2] * cos_eps;
        
        // Rimuovi la ridichiarazione - usa quella già definita sopra
        // Eigen::Vector3d asteroid_bary_icrf;  // [AU] ICRF J2000 (già definito sopra)
        // Aggiorna asteroid_bary_icrf (già definito a riga 179)
        asteroid_bary_icrf[0] = asteroid_bary_ecl[0];
        asteroid_bary_icrf[1] = asteroid_bary_ecl[1] * cos_eps - asteroid_bary_ecl[2] * sin_eps;
        asteroid_bary_icrf[2] = asteroid_bary_ecl[1] * sin_eps + asteroid_bary_ecl[2] * cos_eps;
        
        // Converti Terra da eclittico a ICRF
        Eigen::Vector3d earth_helio_icrf;  // [AU] ICRF J2000
        earth_helio_icrf[0] = earth_helio_ecl[0];
        earth_helio_icrf[1] = earth_helio_ecl[1] * cos_eps - earth_helio_ecl[2] * sin_eps;
        earth_helio_icrf[2] = earth_helio_ecl[1] * sin_eps + earth_helio_ecl[2] * cos_eps;
        
        Eigen::Vector3d earth_bary_icrf;  // [AU] ICRF J2000
        earth_bary_icrf[0] = earth_bary_ecl[0];
        earth_bary_icrf[1] = earth_bary_ecl[1] * cos_eps - earth_bary_ecl[2] * sin_eps;
        earth_bary_icrf[2] = earth_bary_ecl[1] * sin_eps + earth_bary_ecl[2] * cos_eps;
        
        // ═══════════════════════════════════════════════════════════════
        // PROVA ENTRAMBE LE ASSUNZIONI: ELIOCENTRICO vs BARICENTRICO
        // ═══════════════════════════════════════════════════════════════
        
        // PROVA 1: Assumiamo ELIOCENTRICHE
        Eigen::Vector3d geocentric_helio;  // [AU] ICRF J2000
        geocentric_helio[0] = asteroid_helio_icrf[0] - earth_helio_icrf[0];
        geocentric_helio[1] = asteroid_helio_icrf[1] - earth_helio_icrf[1];
        geocentric_helio[2] = asteroid_helio_icrf[2] - earth_helio_icrf[2];
        
        // PROVA 2: Assumiamo BARICENTRICHE (metodo originale)
        Eigen::Vector3d geocentric_bary;  // [AU] ICRF J2000
        geocentric_bary[0] = asteroid_bary_icrf[0] - earth_bary_icrf[0];
        geocentric_bary[1] = asteroid_bary_icrf[1] - earth_bary_icrf[1];
        geocentric_bary[2] = asteroid_bary_icrf[2] - earth_bary_icrf[2];
        
        // ═══════════════════════════════════════════════════════════════
        // PROVA 3: Forse cart.position è già GEOCENTRICO?
        // ═══════════════════════════════════════════════════════════════
        Eigen::Vector3d geocentric_direct_icrf;  // [AU] ICRF J2000
        geocentric_direct_icrf[0] = asteroid_pos_ecl[0];
        geocentric_direct_icrf[1] = asteroid_pos_ecl[1] * cos_eps - asteroid_pos_ecl[2] * sin_eps;
        geocentric_direct_icrf[2] = asteroid_pos_ecl[1] * sin_eps + asteroid_pos_ecl[2] * cos_eps;
        
        // Usa ELIOCENTRICO per ora (da testare)
        Eigen::Vector3d geocentric = geocentric_helio;
        
        if (is_debug) {
            double r_helio = geocentric_helio.norm();
            double r_bary = geocentric_bary.norm();
            double r_direct = geocentric_direct_icrf.norm();
            
            double ra_helio_rad = std::atan2(geocentric_helio[1], geocentric_helio[0]);
            double dec_helio_rad = std::asin(geocentric_helio[2] / r_helio);
            double ra_bary_rad = std::atan2(geocentric_bary[1], geocentric_bary[0]);
            double dec_bary_rad = std::asin(geocentric_bary[2] / r_bary);
            double ra_direct_rad = std::atan2(geocentric_direct_icrf[1], geocentric_direct_icrf[0]);
            double dec_direct_rad = std::asin(geocentric_direct_icrf[2] / r_direct);
            
            double ra_helio = ra_helio_rad * 180.0 / M_PI;
            double dec_helio = dec_helio_rad * 180.0 / M_PI;
            double ra_bary = ra_bary_rad * 180.0 / M_PI;
            double dec_bary = dec_bary_rad * 180.0 / M_PI;
            double ra_direct = ra_direct_rad * 180.0 / M_PI;
            double dec_direct = dec_direct_rad * 180.0 / M_PI;
            
            if (ra_helio < 0) ra_helio += 360.0;
            if (ra_bary < 0) ra_bary += 360.0;
            if (ra_direct < 0) ra_direct += 360.0;
            
            std::cerr << "[DEBUG] HELIOCENTRIC assumption: RA=" << ra_helio 
                      << "°, Dec=" << dec_helio << "°, r=" << r_helio << " AU\n";
            std::cerr << "[DEBUG] BARYCENTRIC assumption: RA=" << ra_bary 
                      << "°, Dec=" << dec_bary << "°, r=" << r_bary << " AU\n";
            std::cerr << "[DEBUG] DIRECT (assume geocentric): RA=" << ra_direct 
                      << "°, Dec=" << dec_direct << "°, r=" << r_direct << " AU\n";
            std::cerr << "[DEBUG] Expected: RA=168.10°, Dec=-21.31°\n";
            std::cerr.flush();
        }
        
        // Distanza geocentrica
        double r = geocentric.norm();  // [AU]
        
        if (r < 1e-10) {
            // Asteroide troppo vicino alla Terra (improbabile ma gestiamo)
            return {0.0, 0.0};
        }
        
        // RA e Dec (radianti) dalla posizione geocentrica
        // NOTA: La distanza r è in AU, ma per RA/Dec serve solo la direzione (vettore unitario)
        double ra_rad = std::atan2(geocentric[1], geocentric[0]);  // [rad]
        double dec_rad = std::asin(geocentric[2] / r);  // [rad]
        
        // Converti in gradi
        double ra = ra_rad * RAD_TO_DEG_LOCAL;  // [deg]
        double dec = dec_rad * RAD_TO_DEG_LOCAL;  // [deg]
        
        // Normalizza RA a [0, 360)
        if (ra < 0) ra += 360.0;
        
        return {ra, dec};  // [deg, deg]
    }
};

AstDynPropagator::AstDynPropagator(double tolerance)
    : pimpl_(std::make_unique<Impl>(tolerance))
{
}

AstDynPropagator::~AstDynPropagator() = default;

void AstDynPropagator::setTolerance(double tol) {
    pimpl_->tolerance = tol;
    pimpl_->initializePropagator();
}

void AstDynPropagator::setStepLimits(double h_min_days, double h_max_days) {
    pimpl_->h_min_days = h_min_days;
    pimpl_->h_max_days = h_max_days;
    pimpl_->initializePropagator();
}

void AstDynPropagator::usePlanetPerturbations(bool enable) {
    pimpl_->planets_enabled = enable;
    pimpl_->initializePropagator();
}

void AstDynPropagator::useAsteroidPerturbations(bool enable) {
    pimpl_->asteroids_enabled = enable;
    pimpl_->initializePropagator();
}

void AstDynPropagator::useRelativisticCorrections(bool enable) {
    pimpl_->relativity_enabled = enable;
    pimpl_->initializePropagator();
}

AstDySElements AstDynPropagator::propagate(
    const AstDySElements& elements,
    double target_mjd) {
    
    // Converti in KeplerianElements
    auto kep = pimpl_->toKeplerian(elements);
    
    // Propaga usando propagate_keplerian (non propagate)
    auto propagated = pimpl_->propagator->propagate_keplerian(kep, target_mjd);
    
    // Converti risultato
    return pimpl_->fromKeplerian(propagated);
}

std::pair<double, double> AstDynPropagator::getRADec(
    const AstDySElements& elements,
    double mjd_utc) {
    
    // Converti in KeplerianElements
    auto kep = pimpl_->toKeplerian(elements);
    
    // Propaga e converti in cartesiano
    auto kep_prop = pimpl_->propagator->propagate_keplerian(kep, mjd_utc);
    auto cart = astdyn::propagation::keplerian_to_cartesian(kep_prop);  // ← Restituisce CartesianElements
    
    // Calcola RA/Dec
    return pimpl_->computeRADec(cart, mjd_utc);
}

std::vector<RWOObservation> AstDynPropagator::computeResiduals(
    const AstDySElements& elements,
    const std::vector<RWOObservation>& observations) {
    
    std::vector<RWOObservation> results = observations;
    
    // Converti elementi
    auto kep = pimpl_->toKeplerian(elements);
    
    // Calcola residui per ogni osservazione
    for (auto& obs : results) {
        try {
            // Propaga all'epoca osservazione
            auto radec = getRADec(elements, obs.mjd_utc);
            
            // Calcola residui (O-C)
            obs.ra_residual_arcsec = (radec.first - obs.ra_deg) * 3600.0 * 
                                     std::cos(obs.dec_deg * DEG_TO_RAD_LOCAL);
            obs.dec_residual_arcsec = (radec.second - obs.dec_deg) * 3600.0;
            
            // Chi² normalizzato
            double ra_err_sq = obs.ra_sigma_arcsec * obs.ra_sigma_arcsec;
            double dec_err_sq = obs.dec_sigma_arcsec * obs.dec_sigma_arcsec;
            obs.chi_squared = (obs.ra_residual_arcsec * obs.ra_residual_arcsec / ra_err_sq +
                              obs.dec_residual_arcsec * obs.dec_residual_arcsec / dec_err_sq) / 2.0;
            
            // Outlier detection (3σ)
            obs.is_outlier = (std::abs(obs.ra_residual_arcsec) > 3.0 * obs.ra_sigma_arcsec ||
                             std::abs(obs.dec_residual_arcsec) > 3.0 * obs.dec_sigma_arcsec);
        } catch (...) {
            obs.is_outlier = true;
            obs.chi_squared = 1e10;
        }
    }
    
    return results;
}

int AstDynPropagator::getLastStepsAccepted() const {
    // TODO: Implementare statistiche reali dal propagatore
    return pimpl_->steps_accepted;
}

int AstDynPropagator::getLastStepsRejected() const {
    // TODO: Implementare statistiche reali dal propagatore
    return pimpl_->steps_rejected;
}

double AstDynPropagator::getLastMinStep() const {
    // TODO: Implementare statistiche reali dal propagatore
    return pimpl_->min_step;
}

double AstDynPropagator::getLastMaxStep() const {
    // TODO: Implementare statistiche reali dal propagatore
    return pimpl_->max_step;
}

// ============================================================================
// AstDynOrbitFitter Implementation (stub per ora)
// ============================================================================

class AstDynOrbitFitter::Impl {
public:
    std::unique_ptr<AstDynPropagator> propagator;
    double outlier_threshold;
    int max_iterations;
    double convergence_tolerance;
    bool verbose;
    
    Impl(double tolerance)
        : propagator(std::make_unique<AstDynPropagator>(tolerance))
        , outlier_threshold(3.0)
        , max_iterations(20)
        , convergence_tolerance(1e-6)
        , verbose(false)
    {
    }
};

AstDynOrbitFitter::AstDynOrbitFitter(double tolerance)
    : pimpl_(std::make_unique<Impl>(tolerance))
{
}

AstDynOrbitFitter::~AstDynOrbitFitter() = default;

void AstDynOrbitFitter::setOutlierThreshold(double sigma) {
    pimpl_->outlier_threshold = sigma;
}

void AstDynOrbitFitter::setMaxIterations(int max_iter) {
    pimpl_->max_iterations = max_iter;
}

void AstDynOrbitFitter::setConvergenceTolerance(double tol_au) {
    pimpl_->convergence_tolerance = tol_au;
}

void AstDynOrbitFitter::setVerbose(bool verbose) {
    pimpl_->verbose = verbose;
}

AstDynOrbitFitResult AstDynOrbitFitter::fit(
    const AstDySElements& initial_elements,
    const std::vector<RWOObservation>& observations) {
    
    // TODO: Implementare differential correction
    AstDynOrbitFitResult result;
    result.fitted_elements = initial_elements;
    result.n_observations = observations.size();
    result.n_used = observations.size();
    result.n_outliers = 0;
    result.method = "AstDyn-RKF78 (stub)";
    
    // Calcola residui
    auto residuals = pimpl_->propagator->computeResiduals(initial_elements, observations);
    result.observations = residuals;
    
    // Statistiche base
    double sum_ra_sq = 0.0, sum_dec_sq = 0.0;
    int n_valid = 0;
    for (const auto& obs : residuals) {
        if (!obs.is_outlier) {
            sum_ra_sq += obs.ra_residual_arcsec * obs.ra_residual_arcsec;
            sum_dec_sq += obs.dec_residual_arcsec * obs.dec_residual_arcsec;
            n_valid++;
        }
    }
    
    if (n_valid > 0) {
        result.rms_ra_arcsec = std::sqrt(sum_ra_sq / n_valid);
        result.rms_dec_arcsec = std::sqrt(sum_dec_sq / n_valid);
        result.rms_total_arcsec = std::sqrt((sum_ra_sq + sum_dec_sq) / n_valid);
    }
    
    return result;
}

AstDynOrbitFitResult AstDynOrbitFitter::computeResidualsOnly(
    const AstDySElements& elements,
    const std::vector<RWOObservation>& observations) {
    
    AstDynOrbitFitResult result;
    result.fitted_elements = elements;
    result.n_observations = observations.size();
    result.n_used = observations.size();
    result.n_outliers = 0;
    result.method = "AstDyn-RKF78 (residuals only)";
    
    auto residuals = pimpl_->propagator->computeResiduals(elements, observations);
    result.observations = residuals;
    
    // Conta outliers
    for (const auto& obs : residuals) {
        if (obs.is_outlier) {
            result.n_outliers++;
            result.n_used--;
        }
    }
    
    // Statistiche
    double sum_ra_sq = 0.0, sum_dec_sq = 0.0;
    int n_valid = 0;
    for (const auto& obs : residuals) {
        if (!obs.is_outlier) {
            sum_ra_sq += obs.ra_residual_arcsec * obs.ra_residual_arcsec;
            sum_dec_sq += obs.dec_residual_arcsec * obs.dec_residual_arcsec;
            n_valid++;
        }
    }
    
    if (n_valid > 0) {
        result.rms_ra_arcsec = std::sqrt(sum_ra_sq / n_valid);
        result.rms_dec_arcsec = std::sqrt(sum_dec_sq / n_valid);
        result.rms_total_arcsec = std::sqrt((sum_ra_sq + sum_dec_sq) / n_valid);
    }
    
    return result;
}

} // namespace ioccultcalc
