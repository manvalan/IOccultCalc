/**
 * @file test_horizons_direct.cpp
 * @brief Confronto diretto con Horizons - verifica step-by-step
 */

#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/astdyn_propagation_helper.h"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

int main() {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  CONFRONTO DIRETTO CON HORIZONS - STEP BY STEP          ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    // Carica elementi
    AstDySElements elements = AstDySClient::downloadElements(17030);
    std::cout << "Elementi iniziali:\n";
    std::cout << "  Epoch: " << std::fixed << std::setprecision(5) << elements.epoch_mjd << " MJD TDB\n";
    std::cout << "  a: " << std::setprecision(10) << elements.a << " AU\n";
    std::cout << "  e: " << elements.e << "\n";
    std::cout << "  i: " << elements.i << " deg\n";
    std::cout << "  Omega: " << elements.Omega << " deg\n";
    std::cout << "  omega: " << elements.omega << " deg\n";
    std::cout << "  M: " << elements.M << " deg\n\n";
    
    // Test date: 28 Nov 2025 12:00 UTC = MJD 61007.5 TDB (circa)
    double test_mjd = 61007.5;
    std::cout << "Propagazione a MJD " << test_mjd << " (28 Nov 2025 12:00 UTC)\n\n";
    
    // Converti in KeplerianElements per AstDyn
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
    
    std::cout << "Elementi kepleriani per AstDyn:\n";
    std::cout << "  epoch_mjd_tdb: " << kep.epoch_mjd_tdb << "\n";
    std::cout << "  a: " << kep.semi_major_axis << " AU\n";
    std::cout << "  e: " << kep.eccentricity << "\n";
    std::cout << "  i: " << kep.inclination * 180.0 / M_PI << " deg\n";
    std::cout << "  Omega: " << kep.longitude_ascending_node * 180.0 / M_PI << " deg\n";
    std::cout << "  omega: " << kep.argument_perihelion * 180.0 / M_PI << " deg\n";
    std::cout << "  M: " << kep.mean_anomaly * 180.0 / M_PI << " deg\n\n";
    
    // Propaga con AstDyn
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
    
    std::cout << "Propagazione in corso...\n";
    auto kep_prop = prop.propagate_keplerian(kep, test_mjd);
    
    std::cout << "\nElementi propagati:\n";
    std::cout << "  epoch_mjd_tdb: " << kep_prop.epoch_mjd_tdb << "\n";
    std::cout << "  a: " << kep_prop.semi_major_axis << " AU\n";
    std::cout << "  e: " << kep_prop.eccentricity << "\n";
    std::cout << "  i: " << kep_prop.inclination * 180.0 / M_PI << " deg\n";
    std::cout << "  Omega: " << kep_prop.longitude_ascending_node * 180.0 / M_PI << " deg\n";
    std::cout << "  omega: " << kep_prop.argument_perihelion * 180.0 / M_PI << " deg\n";
    std::cout << "  M: " << kep_prop.mean_anomaly * 180.0 / M_PI << " deg\n\n";
    
    // Converti in cartesiano
    auto cart = astdyn::propagation::keplerian_to_cartesian(kep_prop);
    
    std::cout << "Posizione baricentrica eclittica J2000:\n";
    std::cout << "  [" << std::setprecision(10) << cart.position[0] 
              << ", " << cart.position[1] 
              << ", " << cart.position[2] << "] AU\n\n";
    
    // Converti a RA/Dec usando il nostro metodo
    auto& helper = AstDynPropagationHelper::getInstance();
    auto radec = helper.getRADec(elements, test_mjd);
    
    std::cout << "RA/Dec calcolati:\n";
    std::cout << "  RA: " << std::setprecision(8) << radec.first << " deg\n";
    std::cout << "  Dec: " << radec.second << " deg\n\n";
    
    // Confronto con Horizons
    std::cout << "Confronto con Horizons (28 Nov 2025 12:00 UTC):\n";
    std::cout << "  Horizons RA: 73.316067 deg (04h 53m 15.856s)\n";
    std::cout << "  Horizons Dec: 20.325164 deg (+20° 19' 30.59\")\n";
    std::cout << "  Calcolo RA: " << radec.first << " deg\n";
    std::cout << "  Calcolo Dec: " << radec.second << " deg\n";
    
    double ra_diff = std::abs(radec.first - 73.316067);
    if (ra_diff > 180.0) ra_diff = 360.0 - ra_diff;
    double dec_diff = std::abs(radec.second - 20.325164);
    
    std::cout << "\n  Differenza RA: " << ra_diff << " deg = " << ra_diff * 3600.0 << " arcsec\n";
    std::cout << "  Differenza Dec: " << dec_diff << " deg = " << dec_diff * 3600.0 << " arcsec\n";
    
    // Verifica posizione Terra
    double jd_tdb = test_mjd + 2400000.5;
    Eigen::Vector3d earth_helio_ecl = 
        astdyn::ephemeris::PlanetaryEphemeris::getPosition(
            astdyn::ephemeris::CelestialBody::EARTH, jd_tdb);
    Eigen::Vector3d sun_bary_ecl = 
        astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
    Eigen::Vector3d earth_bary_ecl = earth_helio_ecl - sun_bary_ecl;
    
    std::cout << "\nPosizione Terra baricentrica eclittica:\n";
    std::cout << "  [" << earth_bary_ecl[0] << ", " << earth_bary_ecl[1] << ", " << earth_bary_ecl[2] << "] AU\n";
    
    std::cout << "\nPosizione asteroide baricentrica eclittica:\n";
    std::cout << "  [" << cart.position[0] << ", " << cart.position[1] << ", " << cart.position[2] << "] AU\n";
    
    Eigen::Vector3d asteroid_geoc_ecl = cart.position - earth_bary_ecl;
    std::cout << "\nPosizione asteroide geocentrica eclittica:\n";
    std::cout << "  [" << asteroid_geoc_ecl[0] << ", " << asteroid_geoc_ecl[1] << ", " << asteroid_geoc_ecl[2] << "] AU\n";
    
    double r_ecl = asteroid_geoc_ecl.norm();
    std::cout << "  Distanza: " << r_ecl << " AU\n";
    
    // Converti a ICRF
    constexpr double EPSILON_J2000 = 23.4392911 * M_PI / 180.0;
    double cos_eps = std::cos(EPSILON_J2000);
    double sin_eps = std::sin(EPSILON_J2000);
    
    Eigen::Vector3d asteroid_geoc_icrf;
    asteroid_geoc_icrf[0] = asteroid_geoc_ecl[0];
    asteroid_geoc_icrf[1] = asteroid_geoc_ecl[1] * cos_eps - asteroid_geoc_ecl[2] * sin_eps;
    asteroid_geoc_icrf[2] = asteroid_geoc_ecl[1] * sin_eps + asteroid_geoc_ecl[2] * cos_eps;
    
    std::cout << "\nPosizione asteroide geocentrica ICRF:\n";
    std::cout << "  [" << asteroid_geoc_icrf[0] << ", " << asteroid_geoc_icrf[1] << ", " << asteroid_geoc_icrf[2] << "] AU\n";
    
    double r_icrf = asteroid_geoc_icrf.norm();
    double ra_manual = std::atan2(asteroid_geoc_icrf[1], asteroid_geoc_icrf[0]) * 180.0 / M_PI;
    double dec_manual = std::asin(asteroid_geoc_icrf[2] / r_icrf) * 180.0 / M_PI;
    if (ra_manual < 0) ra_manual += 360.0;
    
    std::cout << "\nRA/Dec manuali:\n";
    std::cout << "  RA: " << ra_manual << " deg\n";
    std::cout << "  Dec: " << dec_manual << " deg\n";
    
    return 0;
}

