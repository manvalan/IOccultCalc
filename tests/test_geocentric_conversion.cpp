/**
 * @file test_geocentric_conversion.cpp
 * @brief Test conversione coordinate geocentriche - confronto con Horizons
 */

#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/astdyn_propagation_helper.h"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

using namespace ioccultcalc;

int main() {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST CONVERSIONE COORDINATE GEOCENTRICHE                ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    try {
        // Carica elementi
        AstDySElements elements = AstDySClient::downloadElements(17030);
        
        // Test date: 28 Nov 2025 12:00 UTC = MJD 61007.5 TDB (circa)
        double test_mjd = 61007.5;
        
        std::cout << "Test date: MJD " << std::fixed << std::setprecision(5) << test_mjd 
                  << " (28 Nov 2025 12:00 UTC)\n\n";
        
        // 1. Propaga con AstDyn
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
        
        astdyn::propagation::Propagator prop(std::move(integrator), ephemeris, settings);
        auto kep_prop = prop.propagate_keplerian(kep, test_mjd);
        auto cart = astdyn::propagation::keplerian_to_cartesian(kep_prop);
        
        std::cout << "1. POSIZIONE ASTEROIDE BARICENTRICA ECLITTICA:\n";
        std::cout << "   [" << std::setprecision(10) << cart.position[0] << ", " 
                  << cart.position[1] << ", " << cart.position[2] << "] AU\n\n";
        
        // 2. Posizione Terra
        double jd_tdb = test_mjd + 2400000.5;
        
        // Terra eliocentrica eclittica
        Eigen::Vector3d earth_helio_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getPosition(
                astdyn::ephemeris::CelestialBody::EARTH, jd_tdb);
        
        // Sole baricentrico eclittico
        Eigen::Vector3d sun_bary_ecl = 
            astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
        
        // Terra baricentrica eclittica
        Eigen::Vector3d earth_bary_ecl = earth_helio_ecl - sun_bary_ecl;
        
        std::cout << "2. POSIZIONE TERRA:\n";
        std::cout << "   Terra eliocentrica eclittica: [" 
                  << earth_helio_ecl[0] << ", " << earth_helio_ecl[1] << ", " 
                  << earth_helio_ecl[2] << "] AU\n";
        std::cout << "   Sole baricentrico eclittico: [" 
                  << sun_bary_ecl[0] << ", " << sun_bary_ecl[1] << ", " 
                  << sun_bary_ecl[2] << "] AU\n";
        std::cout << "   Terra baricentrica eclittica: [" 
                  << earth_bary_ecl[0] << ", " << earth_bary_ecl[1] << ", " 
                  << earth_bary_ecl[2] << "] AU\n\n";
        
        // 3. Asteroide geocentrico eclittico
        Eigen::Vector3d asteroid_geoc_ecl = cart.position - earth_bary_ecl;
        
        std::cout << "3. POSIZIONE ASTEROIDE GEOCENTRICA ECLITTICA:\n";
        std::cout << "   [" << asteroid_geoc_ecl[0] << ", " 
                  << asteroid_geoc_ecl[1] << ", " << asteroid_geoc_ecl[2] << "] AU\n";
        double r_ecl = asteroid_geoc_ecl.norm();
        std::cout << "   Distanza: " << r_ecl << " AU\n\n";
        
        // 4. Conversione eclittico -> equatoriale ICRF
        constexpr double EPSILON_J2000 = 23.4392911 * M_PI / 180.0;
        double cos_eps = std::cos(EPSILON_J2000);
        double sin_eps = std::sin(EPSILON_J2000);
        
        Eigen::Vector3d asteroid_geoc_icrf;
        asteroid_geoc_icrf[0] = asteroid_geoc_ecl[0];
        asteroid_geoc_icrf[1] = asteroid_geoc_ecl[1] * cos_eps - asteroid_geoc_ecl[2] * sin_eps;
        asteroid_geoc_icrf[2] = asteroid_geoc_ecl[1] * sin_eps + asteroid_geoc_ecl[2] * cos_eps;
        
        std::cout << "4. POSIZIONE ASTEROIDE GEOCENTRICA ICRF:\n";
        std::cout << "   [" << asteroid_geoc_icrf[0] << ", " 
                  << asteroid_geoc_icrf[1] << ", " << asteroid_geoc_icrf[2] << "] AU\n";
        double r_icrf = asteroid_geoc_icrf.norm();
        std::cout << "   Distanza: " << r_icrf << " AU\n\n";
        
        // 5. RA/Dec
        double ra_rad = std::atan2(asteroid_geoc_icrf[1], asteroid_geoc_icrf[0]);
        double dec_rad = std::asin(asteroid_geoc_icrf[2] / r_icrf);
        double ra_deg = ra_rad * 180.0 / M_PI;
        double dec_deg = dec_rad * 180.0 / M_PI;
        if (ra_deg < 0) ra_deg += 360.0;
        
        std::cout << "5. RA/DEC CALCOLATI:\n";
        std::cout << "   RA: " << std::setprecision(8) << ra_deg << " deg\n";
        std::cout << "   Dec: " << dec_deg << " deg\n\n";
        
        // 6. Confronto con Horizons
        std::cout << "6. CONFRONTO CON HORIZONS (28 Nov 2025 12:00 UTC):\n";
        std::cout << "   Horizons RA: 73.316067 deg (04h 53m 15.856s)\n";
        std::cout << "   Horizons Dec: 20.325164 deg (+20° 19' 30.59\")\n";
        std::cout << "   Calcolo RA: " << ra_deg << " deg\n";
        std::cout << "   Calcolo Dec: " << dec_deg << " deg\n";
        
        double ra_diff = std::abs(ra_deg - 73.316067);
        if (ra_diff > 180.0) ra_diff = 360.0 - ra_diff;
        double dec_diff = std::abs(dec_deg - 20.325164);
        
        std::cout << "\n   Differenza RA: " << ra_diff << " deg = " 
                  << ra_diff * 3600.0 << " arcsec\n";
        std::cout << "   Differenza Dec: " << dec_diff << " deg = " 
                  << dec_diff * 3600.0 << " arcsec\n";
        
        if (ra_diff < 0.1 && dec_diff < 0.1) {
            std::cout << "\n   ✓ COERENTE CON HORIZONS\n";
        } else {
            std::cout << "\n   ⚠ DIFFERENZA SIGNIFICATIVA\n";
            std::cout << "   Verificare:\n";
            std::cout << "     1. Posizione Terra baricentrica\n";
            std::cout << "     2. Conversione eclittico -> equatoriale\n";
            std::cout << "     3. Calcolo RA/Dec\n";
        }
        
        // 7. Verifica usando AstDynPropagationHelper
        std::cout << "\n7. VERIFICA CON AstDynPropagationHelper::getRADec:\n";
        auto& helper = AstDynPropagationHelper::getInstance();
        auto radec_helper = helper.getRADec(elements, test_mjd);
        std::cout << "   RA: " << radec_helper.first << " deg\n";
        std::cout << "   Dec: " << radec_helper.second << " deg\n";
        
        double ra_diff_helper = std::abs(radec_helper.first - ra_deg);
        if (ra_diff_helper > 180.0) ra_diff_helper = 360.0 - ra_diff_helper;
        double dec_diff_helper = std::abs(radec_helper.second - dec_deg);
        
        std::cout << "   Differenza con calcolo manuale:\n";
        std::cout << "     RA: " << ra_diff_helper << " deg\n";
        std::cout << "     Dec: " << dec_diff_helper << " deg\n";
        
        if (ra_diff_helper < 1e-6 && dec_diff_helper < 1e-6) {
            std::cout << "   ✓ COERENTE\n";
        } else {
            std::cout << "   ⚠ DIFFERENZA\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

