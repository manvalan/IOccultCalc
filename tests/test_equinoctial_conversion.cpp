/**
 * @file test_equinoctial_conversion.cpp
 * @brief Verifica conversione elementi equinoctiali - confronto con codice vecchio
 */

#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/astdyn_propagation_helper.h"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

using namespace ioccultcalc;

// Funzione dal codice vecchio (ephemeris.cpp) per convertire elementi equinoctiali
void equinoctialToState_OLD(const EquinoctialElements& elements, double epoch_mjd,
                           Eigen::Vector3d& pos_ecl, Eigen::Vector3d& vel_ecl) {
    constexpr double GAUSS_K = 0.01720209895;  // Gauss constant
    constexpr double OBLIQUITY_J2000 = 23.4392911 * M_PI / 180.0;
    
    // Mean motion
    double n = GAUSS_K / elements.a / sqrt(elements.a);
    
    // Mean longitude at epoch
    double lambda = elements.lambda + n * (epoch_mjd - (elements.epoch.jd - 2400000.5));
    
    // Eccentricity
    double e = sqrt(elements.h * elements.h + elements.k * elements.k);
    
    // Longitude of pericenter
    double omega_plus_Omega = atan2(elements.h, elements.k);
    
    // Mean anomaly
    double M = lambda - omega_plus_Omega;
    while (M < 0) M += 2.0 * M_PI;
    while (M >= 2.0 * M_PI) M -= 2.0 * M_PI;
    
    // Solve Kepler equation
    double E = M;
    for (int i = 0; i < 100; i++) {
        double f = E - e * sin(E) - M;
        double fp = 1.0 - e * cos(E);
        double delta = f / fp;
        E -= delta;
        if (std::abs(delta) < 1e-12) break;
    }
    
    // True anomaly
    double nu = 2.0 * atan2(sqrt(1.0 + e) * sin(E/2.0), sqrt(1.0 - e) * cos(E/2.0));
    
    // Radius
    double r = elements.a * (1.0 - e * cos(E));
    
    // Coordinates in orbital plane
    double x_orb = r * cos(nu);
    double y_orb = r * sin(nu);
    
    // Velocity in orbital plane
    double v_factor = GAUSS_K * sqrt(elements.a) / r;
    double vx_orb = -v_factor * sin(E);
    double vy_orb = v_factor * sqrt(1.0 - e * e) * cos(E);
    
    // Transformation matrix (equinoctial -> reference frame)
    double f = 1.0 + elements.p * elements.p + elements.q * elements.q;
    
    double m11 = (1.0 - elements.p * elements.p + elements.q * elements.q) / f;
    double m12 = 2.0 * elements.p * elements.q / f;
    double m13 = -2.0 * elements.p / f;
    
    double m21 = 2.0 * elements.p * elements.q / f;
    double m22 = (1.0 + elements.p * elements.p - elements.q * elements.q) / f;
    double m23 = 2.0 * elements.q / f;
    
    double m31 = 2.0 * elements.p / f;
    double m32 = -2.0 * elements.q / f;
    double m33 = (1.0 - elements.p * elements.p - elements.q * elements.q) / f;
    
    // Rotate from orbital plane to reference frame
    double cos_wp = cos(omega_plus_Omega);
    double sin_wp = sin(omega_plus_Omega);
    
    double x_ref = x_orb * cos_wp - y_orb * sin_wp;
    double y_ref = x_orb * sin_wp + y_orb * cos_wp;
    double z_ref = 0;
    
    double vx_ref = vx_orb * cos_wp - vy_orb * sin_wp;
    double vy_ref = vx_orb * sin_wp + vy_orb * cos_wp;
    double vz_ref = 0;
    
    // Apply transformation to equinoctial frame
    // NOTE: Equinoctial elements from AstDyS are in ECLIPTIC J2000 (ECLM)
    pos_ecl[0] = m11 * x_ref + m12 * y_ref + m13 * z_ref;
    pos_ecl[1] = m21 * x_ref + m22 * y_ref + m23 * z_ref;
    pos_ecl[2] = m31 * x_ref + m32 * y_ref + m33 * z_ref;
    
    vel_ecl[0] = m11 * vx_ref + m12 * vy_ref + m13 * vz_ref;
    vel_ecl[1] = m21 * vx_ref + m22 * vy_ref + m23 * vz_ref;
    vel_ecl[2] = m31 * vx_ref + m32 * vy_ref + m33 * vz_ref;
}

int main() {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  CONFRONTO CONVERSIONE ELEMENTI EQUINOCTIALI             ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    // Carica elementi
    EquinoctialElements eq = AstDysClient().getElements("17030");
    
    std::cout << "Elementi equinoctiali AstDyS:\n";
    std::cout << "  a: " << std::setprecision(10) << eq.a << " AU\n";
    std::cout << "  h: " << eq.h << "\n";
    std::cout << "  k: " << eq.k << "\n";
    std::cout << "  p: " << eq.p << "\n";
    std::cout << "  q: " << eq.q << "\n";
    std::cout << "  lambda: " << eq.lambda * 180.0 / M_PI << " deg\n";
    std::cout << "  epoch: " << std::fixed << std::setprecision(5) << (eq.epoch.jd - 2400000.5) << " MJD\n\n";
    
    // Test all'epoca
    double test_mjd = eq.epoch.jd - 2400000.5;
    
    // Metodo VECCHIO (ephemeris.cpp)
    Eigen::Vector3d pos_ecl_old, vel_ecl_old;
    equinoctialToState_OLD(eq, test_mjd, pos_ecl_old, vel_ecl_old);
    
    std::cout << "Metodo VECCHIO (ephemeris.cpp):\n";
    std::cout << "  Posizione eclittica: [" << std::setprecision(10) 
              << pos_ecl_old[0] << ", " << pos_ecl_old[1] << ", " << pos_ecl_old[2] << "] AU\n\n";
    
    // Metodo NUOVO (AstDyn)
    AstDySElements astdys = AstDynPropagationHelper::convertFromEquinoctial(eq);
    
    astdyn::propagation::KeplerianElements kep;
    kep.semi_major_axis = astdys.a;
    kep.eccentricity = astdys.e;
    kep.inclination = astdys.i * M_PI / 180.0;
    kep.longitude_ascending_node = astdys.Omega * M_PI / 180.0;
    kep.argument_perihelion = astdys.omega * M_PI / 180.0;
    kep.mean_anomaly = astdys.M * M_PI / 180.0;
    kep.epoch_mjd_tdb = astdys.epoch_mjd;
    kep.gravitational_parameter = 
        1.32712440018e20 / std::pow(1.495978707e11, 3) * std::pow(86400.0, 2);
    
    auto cart = astdyn::propagation::keplerian_to_cartesian(kep);
    
    std::cout << "Metodo NUOVO (AstDyn keplerian_to_cartesian):\n";
    std::cout << "  Posizione eclittica: [" << std::setprecision(10)
              << cart.position[0] << ", " << cart.position[1] << ", " << cart.position[2] << "] AU\n\n";
    
    // Confronto
    Eigen::Vector3d diff = pos_ecl_old - cart.position;
    double diff_mag = diff.norm();
    
    std::cout << "Differenza:\n";
    std::cout << "  [" << diff[0] << ", " << diff[1] << ", " << diff[2] << "] AU\n";
    std::cout << "  Magnitudine: " << diff_mag << " AU = " << diff_mag * 1.495978707e8 << " km\n";
    
    if (diff_mag > 1e-6) {
        std::cout << "\n  ⚠ ATTENZIONE: Differenza significativa!\n";
        std::cout << "  Potrebbe indicare un problema nella conversione degli elementi\n";
    } else {
        std::cout << "\n  ✓ Le conversioni sono coerenti\n";
    }
    
    return 0;
}

