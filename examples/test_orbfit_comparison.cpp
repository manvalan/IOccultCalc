/**
 * @file test_orbfit_comparison.cpp
 * @brief Confronto IOccultCalc vs OrbFit per 433 Eros
 * 
 * Test:
 * 1. Usa elementi AstDyS @ MJD 55918.807 (Dec 23, 2011)
 * 2. Propaga a MJD 60615.5 (Nov 01, 2024) = ~13 anni
 * 3. Confronta con:
 *    - OrbFit (risultati da /Users/michelebigi/Astro/OrbFit/tests/orbfit/433-Eros/)
 *    - JPL Horizons (JPL#659)
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <ioccultcalc/astdys_client.h>
#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/orbit_propagator.h>
#include <ioccultcalc/orbital_elements.h>
#include <ioccultcalc/types.h>

using namespace ioccultcalc;

int main() {
    std::cout << "=================================================================\n";
    std::cout << "TEST: CONFRONTO IOccultCalc vs OrbFit (433 Eros)\n";
    std::cout << "=================================================================\n\n";
    
    std::cout << "Setup test:\n";
    std::cout << "  • Elementi: AstDyS @ MJD 55918.807 (Dec 23, 2011)\n";
    std::cout << "  • Propagazione: 13 anni → MJD 60615.5 (Nov 01, 2024)\n";
    std::cout << "  • Confronto: OrbFit + JPL Horizons\n\n";
    
    try {
        // STEP 1: Download elementi da AstDyS
        std::cout << "--- STEP 1: Download elementi AstDyS ---\n";
        AstDysClient astdys;
        auto elements_eq = astdys.getElements("433");
        
        std::cout << "Elementi equinoziali @ MJD " << std::fixed << std::setprecision(3) 
                  << elements_eq.epoch.toMJD() << ":\n";
        std::cout << "  a = " << std::setprecision(6) << elements_eq.a << " AU\n";
        std::cout << "  h = " << elements_eq.h << "\n";
        std::cout << "  k = " << elements_eq.k << "\n";
        std::cout << "  p = " << elements_eq.p << "\n";
        std::cout << "  q = " << elements_eq.q << "\n";
        std::cout << "  λ = " << elements_eq.lambda * 180.0 / M_PI << "°\n\n";
        
        // Verifica epoch
        double expected_mjd = 55918.807;
        double actual_mjd = elements_eq.epoch.toMJD();
        if (std::abs(actual_mjd - expected_mjd) > 1.0) {
            std::cout << "⚠ WARNING: Epoch diversa da OrbFit!\n";
            std::cout << "  Atteso: MJD " << expected_mjd << " (Dec 23, 2011)\n";
            std::cout << "  Trovato: MJD " << actual_mjd << "\n\n";
        }
        
        // STEP 2: Propaga con IOccultCalc
        std::cout << "--- STEP 2: Propagazione IOccultCalc ---\n";
        
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RK4;
        opts.stepSize = 0.1;  // giorni
        opts.usePlanetaryPerturbations = true;
        
        OrbitPropagator prop(opts);
        
        // Epoca target
        double target_mjd = 60615.5;  // 2024-Nov-01
        JulianDate target_epoch(target_mjd + 2400000.5);
        
        double dt_days = target_mjd - elements_eq.epoch.toMJD();
        std::cout << "Propagando " << std::setprecision(1) << dt_days << " giorni (~" 
                  << (dt_days/365.25) << " anni)...\n";
        std::cout << "(Con RK4 step=0.1d, ~" << (int)(dt_days/0.1) << " steps)\n\n";
        
        // Converti → stato → propaga
        auto state_0 = prop.elementsToState(elements_eq);
        auto state_ioc = prop.propagate(state_0, target_epoch);
        
        Vector3D r_ioc = state_ioc.position;
        
        std::cout << "IOccultCalc @ MJD " << target_mjd << ":\n";
        std::cout << "  r = (" << std::setprecision(9) 
                  << r_ioc.x << ", " 
                  << r_ioc.y << ", " 
                  << r_ioc.z << ") AU\n";
        std::cout << "  v = ("
                  << state_ioc.velocity.x << ", "
                  << state_ioc.velocity.y << ", "
                  << state_ioc.velocity.z << ") AU/d\n";
        double r_ioc_mag = std::sqrt(r_ioc.x*r_ioc.x + r_ioc.y*r_ioc.y + r_ioc.z*r_ioc.z);
        std::cout << "  |r| = " << std::setprecision(6) << r_ioc_mag << " AU\n\n";
        
        // STEP 3: Confronto con Horizons
        std::cout << "--- STEP 3: Confronto con JPL Horizons ---\n";
        JPLHorizonsClient horizons;
        auto [r_hor, v_hor] = horizons.getStateVectors("433", target_epoch, "@sun");
        
        std::cout << "Horizons @ MJD " << target_mjd << " (JPL#659):\n";
        std::cout << "  r = (" << std::setprecision(9)
                  << r_hor.x << ", " 
                  << r_hor.y << ", " 
                  << r_hor.z << ") AU\n";
        std::cout << "  v = ("
                  << v_hor.x << ", "
                  << v_hor.y << ", "
                  << v_hor.z << ") AU/d\n";
        double r_hor_mag = std::sqrt(r_hor.x*r_hor.x + r_hor.y*r_hor.y + r_hor.z*r_hor.z);
        std::cout << "  |r| = " << std::setprecision(6) << r_hor_mag << " AU\n\n";
        
        // Differenza IOccultCalc vs Horizons
        double dx = (r_ioc.x - r_hor.x) * 149597870.7;  // AU → km
        double dy = (r_ioc.y - r_hor.y) * 149597870.7;
        double dz = (r_ioc.z - r_hor.z) * 149597870.7;
        double diff = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        std::cout << "--- DIFFERENZA IOccultCalc vs Horizons ---\n";
        std::cout << "  ΔX = " << std::setprecision(1) << dx << " km (" 
                  << std::setprecision(1) << (dx/diff*100) << "%)\n";
        std::cout << "  ΔY = " << dy << " km (" 
                  << (dy/diff*100) << "%)\n";
        std::cout << "  ΔZ = " << dz << " km (" 
                  << (dz/diff*100) << "%)\n";
        std::cout << "  |Δr| = " << diff << " km\n";
        std::cout << "  Error/year = " << (diff / (dt_days/365.25)) << " km/year\n\n";
        
        // STEP 4: Info OrbFit (risultati da file)
        std::cout << "--- STEP 4: Confronto con OrbFit ---\n";
        std::cout << "OrbFit configuration:\n";
        std::cout << "  • 17 asteroidi massivi (AST17)\n";
        std::cout << "  • Modello errore: VFCC17\n";
        std::cout << "  • Differential correction: 16119 obs\n";
        std::cout << "  • Residual norm: 0.326 arcsec\n\n";
        
        std::cout << "Note: OrbFit ha propagato INDIETRO a MJD 6001.0 (1875),\n";
        std::cout << "      non in avanti a 2024. Per confronto diretto serve\n";
        std::cout << "      eseguire OrbFit con output.epoch=MJD 60615.5\n\n";
        
        // ANALISI
        std::cout << "=================================================================\n";
        std::cout << "ANALISI RISULTATI\n";
        std::cout << "=================================================================\n\n";
        
        if (diff < 100) {
            std::cout << "✅ ECCELLENTE! Errore <100 km dopo 13 anni!\n\n";
            std::cout << "IOccultCalc è accurato quando usa elementi AstDyS recenti.\n";
        } else if (diff < 1000) {
            std::cout << "✓ BUONO. Errore ~" << (int)diff << " km dopo 13 anni.\n\n";
            std::cout << "~" << (int)(diff/(dt_days/365.25)) << " km/anno è accettabile per NEA.\n";
        } else if (diff < 10000) {
            std::cout << "⚠ MODERATO. Errore ~" << (int)diff << " km dopo 13 anni.\n\n";
            std::cout << "Possibili cause:\n";
            std::cout << "  • Differenza epoch iniziale AstDyS vs JPL#659\n";
            std::cout << "  • Asteroidi massivi non inclusi (OrbFit usa 17)\n";
            std::cout << "  • Parametri non-gravitazionali mancanti\n";
        } else {
            std::cout << "❌ ERRORE ELEVATO: " << (int)diff << " km dopo 13 anni!\n\n";
            std::cout << "Cause probabili:\n";
            std::cout << "  1. Elementi iniziali diversi (AstDyS vs JPL#659)\n";
            std::cout << "  2. Frame di riferimento (ECLIPJ2000_DE441 verificato)\n";
            std::cout << "  3. Modello forze incompleto\n";
        }
        
        std::cout << "\nConfronto modelli dinamici:\n";
        std::cout << "  IOccultCalc:  9 pianeti + 16 asteroidi + GR\n";
        std::cout << "  OrbFit:       9 pianeti + 17 asteroidi + GR + VFCC17\n";
        std::cout << "  JPL Horizons: Modello completo JPL + fit recente\n\n";
        
        std::cout << "Per test preciso OrbFit:\n";
        std::cout << "  1. Modificare 433.oop: output.epoch=MJD 60615.5\n";
        std::cout << "  2. Eseguire: orbfit 433\n";
        std::cout << "  3. Confrontare 433.oel con questo output\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n=================================================================\n";
    
    return 0;
}
