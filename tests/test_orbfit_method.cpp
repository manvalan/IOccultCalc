/**
 * @file test_orbfit_method.cpp
 * @brief Test seguendo metodo OrbFit prop2b per conversione elementi equinoctiali
 */

#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/astdyn_propagation_helper.h"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

using namespace ioccultcalc;

int main() {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST METODO ORBFIT prop2b - ASTEROIDE 17030            ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    try {
        // Carica elementi equinoctiali
        EquinoctialElements eq = AstDysClient().getElements("17030");
        
        std::cout << "Elementi equinoctiali AstDyS:\n";
        std::cout << "  Epoch: " << std::fixed << std::setprecision(5) << (eq.epoch.jd - 2400000.5) << " MJD\n";
        std::cout << "  a: " << std::setprecision(10) << eq.a << " AU\n";
        std::cout << "  h: " << eq.h << "\n";
        std::cout << "  k: " << eq.k << "\n";
        std::cout << "  p: " << eq.p << "\n";
        std::cout << "  q: " << eq.q << "\n";
        std::cout << "  lambda: " << eq.lambda * 180.0 / M_PI << " deg\n\n";
        
        // Test all'epoca
        double test_mjd = eq.epoch.jd - 2400000.5;
        
        // 1. METODO ORBFIT (orbit_propagator.cpp - elementsToState)
        // Questo implementa prop2b di OrbFit e restituisce coordinate EQUATORIALI
        OrbitPropagator prop;
        OrbitState state_orbfit = prop.elementsToState(eq);
        
        std::cout << "1. METODO ORBFIT (elementsToState - prop2b):\n";
        std::cout << "   Coordinate EQUATORIALI (come restituisce OrbFit):\n";
        std::cout << "     [" << std::setprecision(10)
                  << state_orbfit.position.x << ", " 
                  << state_orbfit.position.y << ", " 
                  << state_orbfit.position.z << "] AU\n\n";
        
        // 2. Converti coordinate equatoriali OrbFit a eclittiche per confronto
        constexpr double OBLIQUITY_J2000 = 23.4392911 * M_PI / 180.0;
        double cos_eps = std::cos(OBLIQUITY_J2000);
        double sin_eps = std::sin(OBLIQUITY_J2000);
        
        // Conversione inversa: equatoriale -> eclittico
        Eigen::Vector3d pos_ecl_from_orbfit;
        pos_ecl_from_orbfit[0] = state_orbfit.position.x;
        pos_ecl_from_orbfit[1] = state_orbfit.position.y * cos_eps + state_orbfit.position.z * sin_eps;
        pos_ecl_from_orbfit[2] = -state_orbfit.position.y * sin_eps + state_orbfit.position.z * cos_eps;
        
        std::cout << "   Coordinate ECLITTICHE (convertite da equatoriali OrbFit):\n";
        std::cout << "     [" << pos_ecl_from_orbfit[0] << ", " 
                  << pos_ecl_from_orbfit[1] << ", " 
                  << pos_ecl_from_orbfit[2] << "] AU\n\n";
        
        // 3. METODO ASTDYN (keplerian_to_cartesian)
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
        
        std::cout << "2. METODO ASTDYN (keplerian_to_cartesian):\n";
        std::cout << "   Coordinate ECLITTICHE:\n";
        std::cout << "     [" << std::setprecision(10)
                  << cart.position[0] << ", " 
                  << cart.position[1] << ", " 
                  << cart.position[2] << "] AU\n\n";
        
        // Confronto coordinate eclittiche
        Eigen::Vector3d diff_ecl;
        diff_ecl[0] = pos_ecl_from_orbfit[0] - cart.position[0];
        diff_ecl[1] = pos_ecl_from_orbfit[1] - cart.position[1];
        diff_ecl[2] = pos_ecl_from_orbfit[2] - cart.position[2];
        double diff_mag_ecl = diff_ecl.norm();
        
        std::cout << "CONFRONTO COORDINATE ECLITTICHE:\n";
        std::cout << "  OrbFit (convertito) vs AstDyn:\n";
        std::cout << "    Differenza: [" << diff_ecl[0] << ", " 
                  << diff_ecl[1] << ", " << diff_ecl[2] << "] AU\n";
        std::cout << "    Magnitudine: " << diff_mag_ecl << " AU = " 
                  << diff_mag_ecl * 1.495978707e8 << " km\n";
        
        if (diff_mag_ecl < 1e-6) {
            std::cout << "    ✓ COERENTI\n";
        } else {
            std::cout << "    ⚠ DIFFERENZA SIGNIFICATIVA\n";
            if (std::abs(diff_ecl[2]) > 0.1) {
                std::cout << "    ⚠ Componente Z molto diversa: " << diff_ecl[2] << " AU\n";
                std::cout << "    Questo suggerisce un problema nella conversione elementi equinoctiali -> kepleriani\n";
            }
        }
        
        // 4. Verifica elementi kepleriani calcolati
        std::cout << "\n3. ELEMENTI KEPLERIANI CALCOLATI:\n";
        std::cout << "   a: " << astdys.a << " AU\n";
        std::cout << "   e: " << astdys.e << "\n";
        std::cout << "   i: " << astdys.i << " deg\n";
        std::cout << "   Omega: " << astdys.Omega << " deg\n";
        std::cout << "   omega: " << astdys.omega << " deg\n";
        std::cout << "   M: " << astdys.M << " deg\n\n";
        
        // 5. Verifica passo per passo la conversione OrbFit
        std::cout << "4. VERIFICA PASSO PER PASSO ORBFIT prop2b:\n";
        
        // Mean motion
        constexpr double GM_SUN_AU = 0.0002959122082855911;  // AU^3/day^2
        double enne = std::sqrt(GM_SUN_AU / (eq.a * eq.a * eq.a));
        std::cout << "   Mean motion (n): " << enne << " rad/day\n";
        
        // Mean longitude
        double t0 = eq.epoch.jd - 2400000.5;
        double t1 = t0;  // No propagation
        double pml = eq.lambda + enne * (t1 - t0);
        std::cout << "   Mean longitude (pml): " << pml * 180.0 / M_PI << " deg\n";
        
        // Longitude of pericenter
        double ecc2 = eq.h * eq.h + eq.k * eq.k;
        double pol = std::atan2(eq.h, eq.k);
        std::cout << "   Longitude of pericenter (pol): " << pol * 180.0 / M_PI << " deg\n";
        
        // Broucke vectors
        double upq = 1.0 + eq.p * eq.p + eq.q * eq.q;
        double f_z = -2.0 * eq.p / upq;
        double g_z = 2.0 * eq.q / upq;
        std::cout << "   Broucke f.z: " << f_z << "\n";
        std::cout << "   Broucke g.z: " << g_z << "\n";
        
        // Confronta con matrice vecchia
        double m31 = 2.0 * eq.p / upq;
        double m32 = -2.0 * eq.q / upq;
        std::cout << "   Matrice vecchia m31: " << m31 << " (dovrebbe essere opposto di f.z)\n";
        std::cout << "   Matrice vecchia m32: " << m32 << " (dovrebbe essere opposto di g.z)\n";
        std::cout << "   Differenza segno: f.z = " << f_z << " vs m31 = " << m31 << "\n";
        std::cout << "   Differenza segno: g.z = " << g_z << " vs m32 = " << m32 << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

