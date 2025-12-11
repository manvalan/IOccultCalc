/**
 * @file test_propagation_debug.cpp
 * @brief Test di debug per la propagazione asteroide 17030
 * 
 * Verifica la propagazione e confronta con phase1_candidate_screening
 */

#include "ioccultcalc/astdyn_propagation_helper.h"
#include "ioccultcalc/astdyn_interface.h"
#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/time_utils.h"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

// Colori per output
#define RESET   "\033[0m"
#define GREEN   "\033[32m"
#define RED     "\033[31m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define CYAN    "\033[36m"
#define MAGENTA "\033[35m"

void printHeader(const std::string& title) {
    std::cout << "\n" << BLUE << "═══════════════════════════════════════════════════════\n";
    std::cout << title << "\n";
    std::cout << "═══════════════════════════════════════════════════════" << RESET << "\n\n";
}

void printValue(const std::string& label, double value, int precision = 6) {
    std::cout << "  " << CYAN << label << ": " << RESET 
              << std::fixed << std::setprecision(precision) << value << "\n";
}

void printVector(const std::string& label, const Eigen::Vector3d& v, int precision = 6) {
    std::cout << "  " << CYAN << label << ": " << RESET 
              << "[" << std::fixed << std::setprecision(precision) 
              << v[0] << ", " << v[1] << ", " << v[2] << "] AU\n";
}

int main() {
    std::cout << GREEN << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST PROPAGAZIONE DEBUG - ASTEROIDE 17030          ║\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n" << RESET;
    
    try {
        // Data di test: 28 novembre 2025 (metà periodo)
        JulianDate testJD = TimeUtils::calendarToJD(2025, 11, 28, 12, 0, 0);
        double test_mjd = testJD.jd - 2400000.5;
        
        printHeader("Configurazione Test");
        printValue("Data test", test_mjd, 5);
        std::cout << "  " << CYAN << "Data UTC: " << RESET << "2025-11-28 12:00:00\n";
        std::cout << "  " << CYAN << "JD TDB: " << RESET << std::fixed << std::setprecision(8) << testJD.jd << "\n";
        std::cout << "  " << CYAN << "MJD TDB: " << RESET << std::fixed << std::setprecision(8) << test_mjd << "\n\n";
        
        // Carica elementi orbitali asteroide 17030
        printHeader("Caricamento Elementi Orbitali");
        AstDySElements elements = AstDySClient::downloadElements(17030);
        if (elements.a <= 0) {
            std::cerr << RED << "✗ Errore caricamento elementi 17030" << RESET << "\n";
            return 1;
        }
        std::cout << GREEN << "✓ Elementi caricati\n" << RESET;
        printValue("Epoch (MJD TDB)", elements.epoch_mjd, 5);
        printValue("Semi-major axis (AU)", elements.a, 8);
        printValue("Eccentricity", elements.e, 8);
        printValue("Inclination (deg)", elements.i, 6);
        printValue("Longitude of node (deg)", elements.Omega, 6);
        printValue("Argument of perihelion (deg)", elements.omega, 6);
        printValue("Mean anomaly (deg)", elements.M, 6);
        
        // Verifica epoca
        double days_from_epoch = test_mjd - elements.epoch_mjd;
        std::cout << "\n  " << CYAN << "Giorni dall'epoca: " << RESET 
                  << std::fixed << std::setprecision(2) << days_from_epoch << " giorni\n";
        
        // Test 1: Propagazione usando AstDynPropagationHelper
        printHeader("Test 1: AstDynPropagationHelper::getRADec");
        auto& helper = AstDynPropagationHelper::getInstance();
        auto radec_helper = helper.getRADec(elements, test_mjd);
        printValue("RA (deg)", radec_helper.first, 8);
        printValue("Dec (deg)", radec_helper.second, 8);
        
        // Test 2: Propagazione usando AstDynPropagator direttamente
        printHeader("Test 2: AstDynPropagator::getRADec");
        AstDynPropagator propagator(1e-12);
        auto radec_propagator = propagator.getRADec(elements, test_mjd);
        printValue("RA (deg)", radec_propagator.first, 8);
        printValue("Dec (deg)", radec_propagator.second, 8);
        
        // Test 3: Propagazione step-by-step per vedere i valori intermedi
        printHeader("Test 3: Step-by-Step Propagation");
        
        // Converti in KeplerianElements
        auto kep = propagator.propagate(elements, elements.epoch_mjd);  // All'epoca
        std::cout << "  " << CYAN << "Elementi all'epoca:\n" << RESET;
        printValue("  a (AU)", kep.a, 8);
        printValue("  e", kep.e, 8);
        printValue("  i (deg)", kep.i, 6);
        printValue("  M (deg)", kep.M, 6);
        
        // Propaga all'epoca di test
        auto kep_prop = propagator.propagate(elements, test_mjd);
        std::cout << "\n  " << CYAN << "Elementi propagati a MJD " << test_mjd << ":\n" << RESET;
        printValue("  a (AU)", kep_prop.a, 8);
        printValue("  e", kep_prop.e, 8);
        printValue("  i (deg)", kep_prop.i, 6);
        printValue("  M (deg)", kep_prop.M, 6);
        
        // Test 4: Verifica posizione baricentrica eclittica
        printHeader("Test 4: Posizione Baricentrica Eclittica");
        
        // Usa AstDyn direttamente per ottenere posizione cartesiana
        // (gli include sono già presenti tramite astdyn_interface.h)
        
        astdyn::propagation::KeplerianElements kep_astdyn;
        kep_astdyn.semi_major_axis = elements.a;
        kep_astdyn.eccentricity = elements.e;
        kep_astdyn.inclination = elements.i * M_PI / 180.0;
        kep_astdyn.longitude_ascending_node = elements.Omega * M_PI / 180.0;
        kep_astdyn.argument_perihelion = elements.omega * M_PI / 180.0;
        kep_astdyn.mean_anomaly = elements.M * M_PI / 180.0;
        kep_astdyn.epoch_mjd_tdb = elements.epoch_mjd;
        kep_astdyn.gravitational_parameter = 
            1.32712440018e20 / std::pow(1.495978707e11, 3) * std::pow(86400.0, 2);
        
        // Propaga
        auto ephemeris = std::make_shared<astdyn::ephemeris::PlanetaryEphemeris>();
        auto integrator = std::make_unique<astdyn::propagation::RKF78Integrator>(0.1, 1e-12);
        astdyn::propagation::PropagatorSettings settings;
        settings.include_planets = true;
        settings.include_moon = true;
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
        
        astdyn::propagation::Propagator astdyn_prop(std::move(integrator), ephemeris, settings);
        auto kep_prop_astdyn = astdyn_prop.propagate_keplerian(kep_astdyn, test_mjd);
        auto cart_ecl = astdyn::propagation::keplerian_to_cartesian(kep_prop_astdyn);
        
        std::cout << "  " << CYAN << "Posizione baricentrica eclittica (AstDyn):\n" << RESET;
        printVector("  position [AU]", cart_ecl.position, 10);
        printVector("  velocity [AU/day]", cart_ecl.velocity, 10);
        
        // Test 5: Posizione Terra
        printHeader("Test 5: Posizione Terra");
        double jd_tdb = test_mjd + 2400000.5;
        Eigen::Vector3d earth_helio_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getPosition(
                astdyn::ephemeris::CelestialBody::EARTH, jd_tdb);
        Eigen::Vector3d sun_bary_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
        Eigen::Vector3d earth_bary_ecl = earth_helio_ecl - sun_bary_ecl;
        
        std::cout << "  " << CYAN << "Terra eliocentrica eclittica:\n" << RESET;
        printVector("  position [AU]", earth_helio_ecl, 10);
        std::cout << "  " << CYAN << "Sole baricentrico eclittica:\n" << RESET;
        printVector("  position [AU]", sun_bary_ecl, 10);
        std::cout << "  " << CYAN << "Terra baricentrica eclittica:\n" << RESET;
        printVector("  position [AU]", earth_bary_ecl, 10);
        
        // Test 6: Conversione a ICRF e posizione geocentrica
        printHeader("Test 6: Conversione ICRF e Posizione Geocentrica");
        constexpr double EPSILON_J2000 = 23.4392911 * M_PI / 180.0;
        double cos_eps = std::cos(EPSILON_J2000);
        double sin_eps = std::sin(EPSILON_J2000);
        
        Eigen::Vector3d asteroid_bary_icrf;
        asteroid_bary_icrf[0] = cart_ecl.position[0];
        asteroid_bary_icrf[1] = cart_ecl.position[1] * cos_eps - cart_ecl.position[2] * sin_eps;
        asteroid_bary_icrf[2] = cart_ecl.position[1] * sin_eps + cart_ecl.position[2] * cos_eps;
        
        Eigen::Vector3d earth_bary_icrf;
        earth_bary_icrf[0] = earth_bary_ecl[0];
        earth_bary_icrf[1] = earth_bary_ecl[1] * cos_eps - earth_bary_ecl[2] * sin_eps;
        earth_bary_icrf[2] = earth_bary_ecl[1] * sin_eps + earth_bary_ecl[2] * cos_eps;
        
        Eigen::Vector3d geocentric = asteroid_bary_icrf - earth_bary_icrf;
        
        std::cout << "  " << CYAN << "Asteroide baricentrico ICRF:\n" << RESET;
        printVector("  position [AU]", asteroid_bary_icrf, 10);
        std::cout << "  " << CYAN << "Terra baricentrica ICRF:\n" << RESET;
        printVector("  position [AU]", earth_bary_icrf, 10);
        std::cout << "  " << CYAN << "Asteroide geocentrico ICRF:\n" << RESET;
        printVector("  position [AU]", geocentric, 10);
        printValue("  distance [AU]", geocentric.norm(), 10);
        
        // Calcola RA/Dec
        double r = geocentric.norm();
        double ra_rad = std::atan2(geocentric[1], geocentric[0]);
        double dec_rad = std::asin(geocentric[2] / r);
        double ra_deg = ra_rad * 180.0 / M_PI;
        double dec_deg = dec_rad * 180.0 / M_PI;
        if (ra_deg < 0) ra_deg += 360.0;
        
        std::cout << "\n  " << CYAN << "RA/Dec calcolati:\n" << RESET;
        printValue("  RA (deg)", ra_deg, 8);
        printValue("  Dec (deg)", dec_deg, 8);
        
        // Confronto
        printHeader("Confronto Risultati");
        std::cout << "  " << CYAN << "AstDynPropagationHelper: " << RESET 
                  << "RA=" << std::fixed << std::setprecision(6) << radec_helper.first 
                  << "°, Dec=" << radec_helper.second << "°\n";
        std::cout << "  " << CYAN << "AstDynPropagator: " << RESET 
                  << "RA=" << radec_propagator.first 
                  << "°, Dec=" << radec_propagator.second << "°\n";
        std::cout << "  " << CYAN << "Step-by-step: " << RESET 
                  << "RA=" << ra_deg 
                  << "°, Dec=" << dec_deg << "°\n";
        
        double diff_ra_helper = std::abs(radec_helper.first - ra_deg);
        double diff_dec_helper = std::abs(radec_helper.second - dec_deg);
        double diff_ra_prop = std::abs(radec_propagator.first - ra_deg);
        double diff_dec_prop = std::abs(radec_propagator.second - dec_deg);
        
        std::cout << "\n  " << CYAN << "Differenze:\n" << RESET;
        printValue("  Helper vs Step-by-step (RA arcsec)", diff_ra_helper * 3600.0, 3);
        printValue("  Helper vs Step-by-step (Dec arcsec)", diff_dec_helper * 3600.0, 3);
        printValue("  Propagator vs Step-by-step (RA arcsec)", diff_ra_prop * 3600.0, 3);
        printValue("  Propagator vs Step-by-step (Dec arcsec)", diff_dec_prop * 3600.0, 3);
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Exception: " << e.what() << RESET << "\n";
        return 1;
    }
    
    std::cout << "\n" << GREEN << "✓ Test completato\n" << RESET;
    return 0;
}

