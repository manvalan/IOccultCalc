/**
 * @file test_phase1_vs_horizons.cpp
 * @brief Confronto propagazione Phase1 vs JPL Horizons
 */

#include "phase1_candidate_screening.h"
#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/jpl_horizons_client.h"
#include "ioccultcalc/astdyn_propagation_helper.h"
#include "ioccultcalc/time_utils.h"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Integrator.hpp"
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

void printVector(const std::string& label, const Eigen::Vector3d& v, int precision = 8) {
    std::cout << "  " << CYAN << label << ": " << RESET 
              << "[" << std::fixed << std::setprecision(precision) 
              << v[0] << ", " << v[1] << ", " << v[2] << "] AU\n";
}

// Converti RA/Dec a vettore unitario ICRF
Eigen::Vector3d radecToVector(double ra_deg, double dec_deg) {
    double ra_rad = ra_deg * M_PI / 180.0;
    double dec_rad = dec_deg * M_PI / 180.0;
    return Eigen::Vector3d(
        std::cos(dec_rad) * std::cos(ra_rad),
        std::cos(dec_rad) * std::sin(ra_rad),
        std::sin(dec_rad)
    );
}

// Calcola distanza angolare tra due direzioni (arcsec)
double angularDistance(double ra1_deg, double dec1_deg, double ra2_deg, double dec2_deg) {
    Eigen::Vector3d v1 = radecToVector(ra1_deg, dec1_deg);
    Eigen::Vector3d v2 = radecToVector(ra2_deg, dec2_deg);
    double dot = v1.dot(v2);
    dot = std::max(-1.0, std::min(1.0, dot));
    double angle_rad = std::acos(dot);
    return angle_rad * 180.0 / M_PI * 3600.0; // arcsec
}

int main() {
    std::cout << GREEN << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  CONFRONTO PHASE1 vs JPL HORIZONS - ASTEROIDE 17030 ║\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n" << RESET;
    
    try {
        // Data di test: 28 novembre 2025 12:00 UTC
        JulianDate testJD = TimeUtils::calendarToJD(2025, 11, 28, 12, 0, 0);
        double test_mjd = testJD.jd - 2400000.5;
        
        printHeader("Configurazione Test");
        printValue("Data test (MJD TDB)", test_mjd, 5);
        std::cout << "  " << CYAN << "Data UTC: " << RESET << "2025-11-28 12:00:00\n";
        std::cout << "  " << CYAN << "JD TDB: " << RESET << std::fixed << std::setprecision(8) << testJD.jd << "\n\n";
        
        // Carica elementi
        printHeader("Caricamento Elementi Orbitali");
        AstDySElements elements = AstDySClient::downloadElements(17030);
        if (elements.a <= 0) {
            std::cerr << RED << "✗ Errore caricamento elementi" << RESET << "\n";
            return 1;
        }
        std::cout << GREEN << "✓ Elementi caricati\n" << RESET;
        printValue("Epoch (MJD TDB)", elements.epoch_mjd, 5);
        printValue("a (AU)", elements.a, 8);
        printValue("e", elements.e, 8);
        printValue("Giorni dall'epoca", test_mjd - elements.epoch_mjd, 2);
        
        // Test 1: Propagazione Phase1
        printHeader("Test 1: Propagazione Phase1");
        Phase1CandidateScreening phase1;
        
        // Converti AstDySElements in KeplerianElements
        astdyn::propagation::KeplerianElements kep;
        kep.semi_major_axis = elements.a;
        kep.eccentricity = elements.e;
        kep.inclination = elements.i * M_PI / 180.0;
        kep.longitude_ascending_node = elements.Omega * M_PI / 180.0;
        kep.argument_perihelion = elements.omega * M_PI / 180.0;
        kep.mean_anomaly = elements.M * M_PI / 180.0;
        kep.epoch_mjd_tdb = elements.epoch_mjd;
        kep.gravitational_parameter = 
            1.32712440018e20 / std::pow(1.495978707e11, 3) * std::pow(86400.0, 2);
        
        phase1.setOrbitalElements(kep);
        
        // Usa AstDynPropagationHelper per ottenere RA/Dec (stesso metodo usato da Phase1 internamente)
        // Phase1 usa il propagator della libreria, ma con include_planets=false
        // Per ora usiamo Helper che usa tutte le perturbazioni per confronto
        
        // Per ottenere la posizione Phase1, dobbiamo replicare il calcolo
        // Phase1 usa: propagator->propagate_keplerian() con include_planets=false
        // Quindi creiamo un propagator identico a quello di Phase1
        auto ephemeris = std::make_shared<astdyn::ephemeris::PlanetaryEphemeris>();
        auto integrator = std::make_unique<astdyn::propagation::RKF78Integrator>(0.5, 1e-9);
        astdyn::propagation::PropagatorSettings settings;
        settings.include_planets = false;  // Come Phase1
        settings.include_asteroids = false;
        settings.include_relativity = false;
        astdyn::propagation::Propagator phase1_propagator(std::move(integrator), ephemeris, settings);
        
        // Propaga come fa Phase1
        auto kep_prop = phase1_propagator.propagate_keplerian(kep, test_mjd);
        auto cart_ecl = astdyn::propagation::keplerian_to_cartesian(kep_prop);
        
        // Converti da eclittico a equatoriale (come Phase1)
        constexpr double EPSILON_J2000 = 23.4392911 * M_PI / 180.0;
        double cos_eps = std::cos(EPSILON_J2000);
        double sin_eps = std::sin(EPSILON_J2000);
        Eigen::Vector3d ast_bary_icrf;
        ast_bary_icrf[0] = cart_ecl.position[0];
        ast_bary_icrf[1] = cart_ecl.position[1] * cos_eps - cart_ecl.position[2] * sin_eps;
        ast_bary_icrf[2] = cart_ecl.position[1] * sin_eps + cart_ecl.position[2] * cos_eps;
        
        // Posizione Terra (come Phase1)
        double jd_tdb = test_mjd + 2400000.5;
        Eigen::Vector3d earth_helio_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getPosition(
                astdyn::ephemeris::CelestialBody::EARTH, jd_tdb);
        Eigen::Vector3d sun_bary_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
        Eigen::Vector3d earth_bary_ecl = earth_helio_ecl - sun_bary_ecl;
        Eigen::Vector3d earth_bary_icrf;
        earth_bary_icrf[0] = earth_bary_ecl[0];
        earth_bary_icrf[1] = earth_bary_ecl[1] * cos_eps - earth_bary_ecl[2] * sin_eps;
        earth_bary_icrf[2] = earth_bary_ecl[1] * sin_eps + earth_bary_ecl[2] * cos_eps;
        
        // Posizione geocentrica
        Eigen::Vector3d ast_geo_icrf = ast_bary_icrf - earth_bary_icrf;
        double r = ast_geo_icrf.norm();
        double ra_phase1_rad = std::atan2(ast_geo_icrf[1], ast_geo_icrf[0]);
        double dec_phase1_rad = std::asin(ast_geo_icrf[2] / r);
        double ra_phase1_deg = ra_phase1_rad * 180.0 / M_PI;
        double dec_phase1_deg = dec_phase1_rad * 180.0 / M_PI;
        if (ra_phase1_deg < 0) ra_phase1_deg += 360.0;
        
        double ra_phase1_deg_final = ra_phase1_deg;
        double dec_phase1_deg_final = dec_phase1_deg;
        
        std::cout << GREEN << "✓ Posizione Phase1 calcolata\n" << RESET;
        printValue("RA Phase1 (deg)", ra_phase1_deg_final, 8);
        printValue("Dec Phase1 (deg)", dec_phase1_deg_final, 8);
        
        // Test 2: Propagazione diretta con AstDynPropagationHelper
        printHeader("Test 2: AstDynPropagationHelper");
        auto& helper = AstDynPropagationHelper::getInstance();
        auto radec_helper = helper.getRADec(elements, test_mjd);
        printValue("RA Helper (deg)", radec_helper.first, 8);
        printValue("Dec Helper (deg)", radec_helper.second, 8);
        
        // Test 3: JPL Horizons
        printHeader("Test 3: JPL Horizons");
        JPLHorizonsClient horizons;
        
        // Scarica posizione eliocentrica da Horizons
        auto [pos_helio_horizons, vel_helio_horizons] = 
            horizons.getStateVectors("17030", testJD, "@sun");
        
        std::cout << "  " << CYAN << "Posizione eliocentrica (Horizons):\n" << RESET;
        std::cout << "    [" << std::fixed << std::setprecision(10)
                  << pos_helio_horizons.x << ", "
                  << pos_helio_horizons.y << ", "
                  << pos_helio_horizons.z << "] AU\n";
        
        // Scarica posizione Terra
        auto [pos_earth_helio, vel_earth_helio] = 
            horizons.getStateVectors("399", testJD, "@sun");
        
        // Posizione geocentrica asteroide (Horizons)
        Vector3D pos_geo_horizons;
        pos_geo_horizons.x = pos_helio_horizons.x - pos_earth_helio.x;
        pos_geo_horizons.y = pos_helio_horizons.y - pos_earth_helio.y;
        pos_geo_horizons.z = pos_helio_horizons.z - pos_earth_helio.z;
        
        // Converti in RA/Dec
        double r_horizons = std::sqrt(pos_geo_horizons.x*pos_geo_horizons.x + 
                            pos_geo_horizons.y*pos_geo_horizons.y + 
                            pos_geo_horizons.z*pos_geo_horizons.z);
        double ra_horizons_rad = std::atan2(pos_geo_horizons.y, pos_geo_horizons.x);
        double dec_horizons_rad = std::asin(pos_geo_horizons.z / r_horizons);
        double ra_horizons_deg = ra_horizons_rad * 180.0 / M_PI;
        double dec_horizons_deg = dec_horizons_rad * 180.0 / M_PI;
        if (ra_horizons_deg < 0) ra_horizons_deg += 360.0;
        
        printValue("RA Horizons (deg)", ra_horizons_deg, 8);
        printValue("Dec Horizons (deg)", dec_horizons_deg, 8);
        printValue("Distance (AU)", r_horizons, 8);
        
        // Confronto
        printHeader("Confronto Risultati");
        std::cout << "  " << CYAN << "Phase1:\n" << RESET;
        printValue("    RA (deg)", ra_phase1_deg_final, 8);
        printValue("    Dec (deg)", dec_phase1_deg_final, 8);
        
        std::cout << "  " << CYAN << "Helper:\n" << RESET;
        printValue("    RA (deg)", radec_helper.first, 8);
        printValue("    Dec (deg)", radec_helper.second, 8);
        
        std::cout << "  " << CYAN << "Horizons:\n" << RESET;
        printValue("    RA (deg)", ra_horizons_deg, 8);
        printValue("    Dec (deg)", dec_horizons_deg, 8);
        
        // Calcola differenze
        double diff_ra_phase1 = std::abs(ra_phase1_deg_final - ra_horizons_deg);
        double diff_dec_phase1 = std::abs(dec_phase1_deg_final - dec_horizons_deg);
        if (diff_ra_phase1 > 180.0) diff_ra_phase1 = 360.0 - diff_ra_phase1;
        
        double diff_ra_helper = std::abs(radec_helper.first - ra_horizons_deg);
        double diff_dec_helper = std::abs(radec_helper.second - dec_horizons_deg);
        if (diff_ra_helper > 180.0) diff_ra_helper = 360.0 - diff_ra_helper;
        
        double dist_phase1 = angularDistance(ra_phase1_deg_final, dec_phase1_deg_final, 
                                            ra_horizons_deg, dec_horizons_deg);
        double dist_helper = angularDistance(radec_helper.first, radec_helper.second,
                                            ra_horizons_deg, dec_horizons_deg);
        
        std::cout << "\n  " << CYAN << "Differenze vs Horizons:\n" << RESET;
        printValue("    Phase1 RA (arcsec)", diff_ra_phase1 * 3600.0, 3);
        printValue("    Phase1 Dec (arcsec)", diff_dec_phase1 * 3600.0, 3);
        printValue("    Phase1 Distanza angolare (arcsec)", dist_phase1, 3);
        
        printValue("    Helper RA (arcsec)", diff_ra_helper * 3600.0, 3);
        printValue("    Helper Dec (arcsec)", diff_dec_helper * 3600.0, 3);
        printValue("    Helper Distanza angolare (arcsec)", dist_helper, 3);
        
        // Verifica quale propagator usa Phase1
        printHeader("Configurazione Propagator Phase1");
        std::cout << "  " << YELLOW << "NOTA: Phase1 usa propagator della libreria AstDyn\n" << RESET;
        std::cout << "  " << YELLOW << "ma con include_planets = false (solo Keplero)\n" << RESET;
        std::cout << "  " << YELLOW << "per velocità massima nello screening iniziale\n" << RESET;
        
        if (dist_phase1 > 10.0 || dist_helper > 10.0) {
            std::cout << "\n" << RED << "✗ DISCREPANZA SIGNIFICATIVA!\n" << RESET;
            std::cout << "  Le posizioni differiscono da Horizons di più di 10 arcsec\n";
            std::cout << "  Questo potrebbe indicare un problema nella propagazione\n";
        } else {
            std::cout << "\n" << GREEN << "✓ Posizioni coerenti con Horizons\n" << RESET;
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Exception: " << e.what() << RESET << "\n";
        return 1;
    }
    
    std::cout << "\n" << GREEN << "✓ Test completato\n" << RESET;
    return 0;
}

