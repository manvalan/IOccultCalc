/**
 * @file test_epoch_verification.cpp
 * @brief Verifica che la propagazione usi correttamente l'epoca
 */

#include "ioccultcalc/astdyn_propagation_helper.h"
#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/time_utils.h"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

int main() {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  VERIFICA EPOCA PROPAGAZIONE - ASTEROIDE 17030          ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    // Carica elementi
    AstDySElements elements = AstDySClient::downloadElements(17030);
    std::cout << "Elementi caricati:\n";
    std::cout << "  Epoch (MJD TDB): " << std::fixed << std::setprecision(5) << elements.epoch_mjd << "\n";
    std::cout << "  a (AU): " << elements.a << "\n";
    std::cout << "  e: " << elements.e << "\n";
    std::cout << "  M (deg): " << elements.M << "\n\n";
    
    // Test: propaga all'epoca stessa (dovrebbe restituire gli stessi elementi)
    std::cout << "Test 1: Propagazione all'epoca stessa (MJD " << elements.epoch_mjd << "):\n";
    auto& helper = AstDynPropagationHelper::getInstance();
    auto radec_at_epoch = helper.getRADec(elements, elements.epoch_mjd);
    std::cout << "  RA: " << std::fixed << std::setprecision(8) << radec_at_epoch.first << " deg\n";
    std::cout << "  Dec: " << radec_at_epoch.second << " deg\n\n";
    
    // Test: propaga a epoca diversa
    double test_mjd = 61007.5;  // 28 Nov 2025 12:00
    std::cout << "Test 2: Propagazione a MJD " << test_mjd << " (28 Nov 2025 12:00):\n";
    auto radec_test = helper.getRADec(elements, test_mjd);
    std::cout << "  RA: " << radec_test.first << " deg\n";
    std::cout << "  Dec: " << radec_test.second << " deg\n\n";
    
    // Confronto con Horizons
    std::cout << "Confronto con Horizons (28 Nov 2025 12:00 UTC):\n";
    std::cout << "  Horizons RA: 73.316067 deg\n";
    std::cout << "  Horizons Dec: 20.325164 deg\n";
    std::cout << "  Calcolo RA: " << radec_test.first << " deg\n";
    std::cout << "  Calcolo Dec: " << radec_test.second << " deg\n";
    
    double ra_diff = std::abs(radec_test.first - 73.316067);
    double dec_diff = std::abs(radec_test.second - 20.325164);
    if (ra_diff > 180.0) ra_diff = 360.0 - ra_diff;
    
    std::cout << "\n  Differenza RA: " << ra_diff << " deg = " << ra_diff * 3600.0 << " arcsec\n";
    std::cout << "  Differenza Dec: " << dec_diff << " deg = " << dec_diff * 3600.0 << " arcsec\n";
    
    // Verifica step-by-step con AstDyn direttamente
    std::cout << "\n\nTest 3: Propagazione step-by-step con AstDyn:\n";
    
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
    
    std::cout << "  Elementi all'epoca:\n";
    std::cout << "    epoch_mjd_tdb: " << kep.epoch_mjd_tdb << "\n";
    std::cout << "    M: " << kep.mean_anomaly * 180.0 / M_PI << " deg\n";
    
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
    
    astdyn::propagation::Propagator prop(std::move(integrator), ephemeris, settings);
    auto kep_prop = prop.propagate_keplerian(kep, test_mjd);
    
    std::cout << "\n  Elementi propagati a MJD " << test_mjd << ":\n";
    std::cout << "    epoch_mjd_tdb: " << kep_prop.epoch_mjd_tdb << "\n";
    std::cout << "    M: " << kep_prop.mean_anomaly * 180.0 / M_PI << " deg\n";
    std::cout << "    a: " << kep_prop.semi_major_axis << " AU\n";
    std::cout << "    e: " << kep_prop.eccentricity << "\n";
    
    // Converti in cartesiano
    auto cart = astdyn::propagation::keplerian_to_cartesian(kep_prop);
    std::cout << "\n  Posizione baricentrica eclittica:\n";
    std::cout << "    [" << cart.position[0] << ", " << cart.position[1] << ", " << cart.position[2] << "] AU\n";
    
    return 0;
}

