/**
 * @file test_astdys_vs_orbfit.cpp
 * @brief Confronto IOccultCalc vs OrbFit usando elementi AstDyS
 * 
 * Test:
 * 1. Elementi iniziali: AstDyS @ MJD 55918.807 (Dec 23, 2011)
 * 2. Propagazione forward a MJD 60300.0 (Jan 01, 2024) - ~12 anni
 * 3. Confronto con OrbFit (quando disponibile)
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <ioccultcalc/astdys_client.h>
#include <ioccultcalc/orbit_propagator.h>
#include <ioccultcalc/types.h>

using namespace ioccultcalc;

int main() {
    std::cout << "=================================================================\n";
    std::cout << "TEST: IOccultCalc con elementi AstDyS\n";
    std::cout << "=================================================================\n\n";
    
    std::cout << "Setup test:\n";
    std::cout << "  • Elementi: AstDyS @ MJD 55918.807 (Dec 23, 2011)\n";
    std::cout << "  • Propagazione: ~12 anni → MJD 60300.0 (Jan 01, 2024)\n";
    std::cout << "  • Modello: 16 asteroidi massivi + 9 pianeti + GR\n\n";
    
    try {
        // STEP 1: Download elementi da AstDyS
        std::cout << "--- STEP 1: Elementi iniziali AstDyS ---\n";
        AstDysClient astdys;
        auto elements_eq = astdys.getElements("433");
        
        double epoch0_mjd = elements_eq.epoch.toMJD();
        std::cout << "Epoca iniziale: MJD " << std::fixed << std::setprecision(6) 
                  << epoch0_mjd << "\n";
        std::cout << "Elementi equinoziali:\n";
        std::cout << "  a = " << elements_eq.a << " AU\n";
        std::cout << "  h = " << std::setprecision(9) << elements_eq.h << "\n";
        std::cout << "  k = " << elements_eq.k << "\n";
        std::cout << "  p = " << elements_eq.p << "\n";
        std::cout << "  q = " << elements_eq.q << "\n";
        std::cout << "  λ = " << elements_eq.lambda * 180.0 / M_PI << "°\n\n";
        
        // STEP 2: Converti a stato iniziale
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RK4;
        opts.stepSize = 0.1;
        opts.usePlanetaryPerturbations = true;
        
        OrbitPropagator prop(opts);
        auto state_0 = prop.elementsToState(elements_eq);
        
        std::cout << "Stato iniziale @ MJD " << epoch0_mjd << " (Eclittica J2000):\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << state_0.position.x << ", "
                  << state_0.position.y << ", "
                  << state_0.position.z << ") AU\n";
        std::cout << "  v = ("
                  << state_0.velocity.x << ", "
                  << state_0.velocity.y << ", "
                  << state_0.velocity.z << ") AU/d\n";
        double r0_mag = std::sqrt(state_0.position.x * state_0.position.x + 
                                   state_0.position.y * state_0.position.y +
                                   state_0.position.z * state_0.position.z);
        std::cout << "  |r| = " << std::setprecision(6) << r0_mag << " AU\n\n";
        
        // STEP 3: Propaga a MJD 60300.0
        std::cout << "--- STEP 2: Propagazione IOccultCalc ---\n";
        double target_mjd = 60300.0;  // 2024-Jan-01
        JulianDate target_epoch(target_mjd + 2400000.5);
        
        double dt_days = target_mjd - epoch0_mjd;
        std::cout << "Propagando " << std::setprecision(1) << dt_days << " giorni (~" 
                  << (dt_days/365.25) << " anni)...\n";
        std::cout << "(Con RK4 step=0.1d, ~" << (int)(dt_days/0.1) << " steps)\n";
        std::cout << "⏳ Tempo stimato: ~2-3 minuti...\n\n";
        
        auto state_final = prop.propagate(state_0, target_epoch);
        
        Vector3D r = state_final.position;
        Vector3D v = state_final.velocity;
        
        std::cout << "IOccultCalc @ MJD " << target_mjd << " (Eclittica J2000):\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << r.x << ", "
                  << r.y << ", "
                  << r.z << ") AU\n";
        std::cout << "  v = ("
                  << v.x << ", "
                  << v.y << ", "
                  << v.z << ") AU/d\n";
        double r_mag = std::sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
        std::cout << "  |r| = " << std::setprecision(6) << r_mag << " AU\n\n";
        
        // Conversione a frame Equatoriale per confronto futuro con OrbFit
        double eps = 23.439291111 * M_PI / 180.0;  // Obliquità J2000
        
        double rx_eq = r.x;
        double ry_eq = r.y * cos(eps) - r.z * sin(eps);
        double rz_eq = r.y * sin(eps) + r.z * cos(eps);
        
        double vx_eq = v.x;
        double vy_eq = v.y * cos(eps) - v.z * sin(eps);
        double vz_eq = v.y * sin(eps) + v.z * cos(eps);
        
        std::cout << "IOccultCalc @ MJD " << target_mjd << " (Equatoriale J2000):\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << rx_eq << ", "
                  << ry_eq << ", "
                  << rz_eq << ") AU\n";
        std::cout << "  v = ("
                  << vx_eq << ", "
                  << vy_eq << ", "
                  << vz_eq << ") AU/d\n\n";
        
        // STEP 4: Istruzioni per confronto con OrbFit
        std::cout << "=================================================================\n";
        std::cout << "RISULTATO IOccultCalc PRONTO\n";
        std::cout << "=================================================================\n\n";
        
        std::cout << "Per confrontare con OrbFit:\n\n";
        std::cout << "1. Verifica che OrbFit abbia generato 433_test.oel\n";
        std::cout << "2. Controlla i valori @ MJD " << target_mjd << "\n";
        std::cout << "3. Calcola differenza con i valori sopra (frame Equatoriale)\n\n";
        
        std::cout << "Comandi per estrarre risultato OrbFit:\n";
        std::cout << "  cd /Users/michelebigi/Astro/OrbFit/tests/orbfit/433-Eros\n";
        std::cout << "  grep 'CAR' 433_test.oel | head -2\n\n";
        
        std::cout << "Formato atteso OrbFit (Equatoriale J2000):\n";
        std::cout << "  CAR  rx ry rz  vx vy vz\n";
        std::cout << "  MJD  " << target_mjd << " TDT\n\n";
        
        std::cout << "Note:\n";
        std::cout << "• IOccultCalc usa 16 asteroidi massivi\n";
        std::cout << "• OrbFit usa 17 asteroidi massivi\n";
        std::cout << "• Differenze di 100-1000 km sono normali su 12 anni\n";
        std::cout << "• RK4 step fisso vs RA15 adattivo possono dare differenze\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n=================================================================\n";
    
    return 0;
}
