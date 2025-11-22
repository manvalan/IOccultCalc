/**
 * @file test_initial_state_eros.cpp
 * @brief Verifica stato iniziale @ MJD 55918.807
 * 
 * Test: Confronto stato cartesiano a T=0
 *   • IOccultCalc: elementi AstDyS → stato
 *   • Horizons: stato @ MJD 55918.807
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <ioccultcalc/astdys_client.h>
#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/orbit_propagator.h>

using namespace ioccultcalc;

int main() {
    std::cout << "==============================================================\n";
    std::cout << "TEST: STATO INIZIALE @ MJD 55918.807 (Dec 23, 2011)\n";
    std::cout << "==============================================================\n\n";
    
    try {
        // Download elementi AstDyS
        AstDysClient astdys;
        auto elements_eq = astdys.getElements("433");
        
        double epoch_mjd = elements_eq.epoch.toMJD();
        std::cout << "Elementi AstDyS @ MJD " << std::fixed << std::setprecision(6) 
                  << epoch_mjd << ":\n";
        std::cout << "  a = " << elements_eq.a << " AU\n";
        std::cout << "  h = " << elements_eq.h << "\n";
        std::cout << "  k = " << elements_eq.k << "\n";
        std::cout << "  p = " << elements_eq.p << "\n";
        std::cout << "  q = " << elements_eq.q << "\n";
        std::cout << "  λ = " << std::setprecision(9) << elements_eq.lambda * 180.0 / M_PI << "°\n\n";
        
        // Converti a stato cartesiano
        PropagatorOptions opts;
        OrbitPropagator prop(opts);
        auto state_ioc = prop.elementsToState(elements_eq);
        
        std::cout << "IOccultCalc stato @ T=0:\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << state_ioc.position.x << ", "
                  << state_ioc.position.y << ", "
                  << state_ioc.position.z << ") AU\n";
        std::cout << "  v = ("
                  << state_ioc.velocity.x << ", "
                  << state_ioc.velocity.y << ", "
                  << state_ioc.velocity.z << ") AU/d\n";
        double r_ioc = std::sqrt(state_ioc.position.x * state_ioc.position.x + 
                                  state_ioc.position.y * state_ioc.position.y +
                                  state_ioc.position.z * state_ioc.position.z);
        std::cout << "  |r| = " << std::setprecision(9) << r_ioc << " AU\n\n";
        
        // Horizons stato @ stessa epoca
        JPLHorizonsClient horizons;
        auto [pos_hor, vel_hor] = horizons.getStateVectors("433", elements_eq.epoch, "@sun");
        
        std::cout << "Horizons stato @ T=0 (JPL#659):\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << pos_hor.x << ", "
                  << pos_hor.y << ", "
                  << pos_hor.z << ") AU\n";
        std::cout << "  v = ("
                  << vel_hor.x << ", "
                  << vel_hor.y << ", "
                  << vel_hor.z << ") AU/d\n";
        double r_hor_mag = std::sqrt(pos_hor.x * pos_hor.x + pos_hor.y * pos_hor.y + pos_hor.z * pos_hor.z);
        std::cout << "  |r| = " << std::setprecision(9) << r_hor_mag << " AU\n\n";
        
        // DIFFERENZA
        double dx = (state_ioc.position.x - pos_hor.x) * 149597870.7;
        double dy = (state_ioc.position.y - pos_hor.y) * 149597870.7;
        double dz = (state_ioc.position.z - pos_hor.z) * 149597870.7;
        double diff = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        std::cout << "==============================================================\n";
        std::cout << "DIFFERENZA @ T=0\n";
        std::cout << "==============================================================\n";
        std::cout << "  ΔX = " << std::setprecision(3) << dx << " km\n";
        std::cout << "  ΔY = " << dy << " km\n";
        std::cout << "  ΔZ = " << dz << " km\n";
        std::cout << "  |Δr| = " << diff << " km\n\n";
        
        if (diff < 100) {
            std::cout << "✅ OTTIMO! Elementi AstDyS compatibili con Horizons @ T=0\n";
        } else if (diff < 10000) {
            std::cout << "⚠ Gli elementi sono diversi!\n";
            std::cout << "   AstDyS usa fit diverso da JPL#659\n";
        } else {
            std::cout << "❌ ERRORE CRITICO: elementi completamente diversi!\n";
            std::cout << "   Differenza: " << (int)(diff/1000) << ",000 km @ T=0\n";
        }
        
        std::cout << "\nNOTE:\n";
        std::cout << "• Se diff < 100 km: Gli elementi sono 'equivalenti'\n";
        std::cout << "• Se diff > 1000 km: Fit diversi (epoca, osservazioni, modello)\n";
        std::cout << "• JPL#659: fit 2021-May-24 su 9130 obs (1893-2021)\n";
        std::cout << "• AstDyS: fit continuo, aggiornato periodicamente\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n==============================================================\n";
    return 0;
}
