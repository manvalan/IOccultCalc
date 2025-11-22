/**
 * @file test_orbfit_backward.cpp
 * @brief Confronto IOccultCalc vs OrbFit - Propagazione INDIETRO
 * 
 * Test:
 * 1. Usa elementi AstDyS @ MJD 55918.807 (Dec 23, 2011)
 * 2. Propaga INDIETRO a MJD 6001.0 (Apr 23, 1875) - 136 anni
 * 3. Confronta con risultato OrbFit 433.oel
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
    std::cout << "TEST: PROPAGAZIONE INDIETRO IOccultCalc vs OrbFit\n";
    std::cout << "=================================================================\n\n";
    
    std::cout << "Setup test:\n";
    std::cout << "  • Elementi: AstDyS @ MJD 55918.807 (Dec 23, 2011)\n";
    std::cout << "  • Propagazione: 136 anni INDIETRO → MJD 6001.0 (Apr 23, 1875)\n";
    std::cout << "  • Confronto: OrbFit 433.oel (17 asteroidi massivi)\n\n";
    
    try {
        // STEP 1: Download elementi da AstDyS
        std::cout << "--- STEP 1: Download elementi AstDyS ---\n";
        AstDysClient astdys;
        auto elements_eq = astdys.getElements("433");
        
        double epoch_mjd = elements_eq.epoch.toMJD();
        std::cout << "Elementi equinoziali @ MJD " << std::fixed << std::setprecision(3) 
                  << epoch_mjd << ":\n";
        std::cout << "  a = " << std::setprecision(6) << elements_eq.a << " AU\n";
        std::cout << "  h = " << elements_eq.h << "\n";
        std::cout << "  k = " << elements_eq.k << "\n";
        std::cout << "  p = " << elements_eq.p << "\n";
        std::cout << "  q = " << elements_eq.q << "\n";
        std::cout << "  λ = " << elements_eq.lambda * 180.0 / M_PI << "°\n\n";
        
        // STEP 2: Propaga INDIETRO con IOccultCalc
        std::cout << "--- STEP 2: Propagazione IOccultCalc INDIETRO ---\n";
        
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RK4;
        opts.stepSize = 0.1;  // giorni
        opts.usePlanetaryPerturbations = true;
        
        OrbitPropagator prop(opts);
        
        // Epoca target: MJD 6001.0 (1875-Apr-23)
        double target_mjd = 6001.0;
        JulianDate target_epoch(target_mjd + 2400000.5);
        
        double dt_days = target_mjd - epoch_mjd;
        std::cout << "Propagando " << std::setprecision(1) << dt_days << " giorni (~" 
                  << (dt_days/365.25) << " anni) INDIETRO...\n";
        std::cout << "(Con RK4 step=0.1d, ~" << (int)std::abs(dt_days/0.1) << " steps)\n";
        std::cout << "⏳ Questa propagazione richiederà ~10-15 minuti...\n\n";
        
        // Converti → stato → propaga
        auto state_0 = prop.elementsToState(elements_eq);
        auto state_ioc = prop.propagate(state_0, target_epoch);
        
        Vector3D r_ioc = state_ioc.position;
        Vector3D v_ioc = state_ioc.velocity;
        
        std::cout << "IOccultCalc @ MJD " << target_mjd << " (Eclittica J2000):\n";
        std::cout << "  r = (" << std::setprecision(12) 
                  << r_ioc.x << ", " 
                  << r_ioc.y << ", " 
                  << r_ioc.z << ") AU\n";
        std::cout << "  v = ("
                  << v_ioc.x << ", "
                  << v_ioc.y << ", "
                  << v_ioc.z << ") AU/d\n";
        double r_ioc_mag = std::sqrt(r_ioc.x*r_ioc.x + r_ioc.y*r_ioc.y + r_ioc.z*r_ioc.z);
        std::cout << "  |r| = " << std::setprecision(6) << r_ioc_mag << " AU\n\n";
        
        // STEP 3: Confronto con OrbFit
        std::cout << "--- STEP 3: Confronto con OrbFit ---\n";
        
        // Risultati OrbFit da 433.oel (EQUATORIALE J2000)
        // NOTA: OrbFit usa frame EQUATORIALE, IOccultCalc usa ECLITTICA!
        // Dati da /Users/michelebigi/Astro/OrbFit/tests/orbfit/433-Eros/433.oel
        double orbfit_rx = -0.109914;   // AU
        double orbfit_ry = -1.671713;   // AU
        double orbfit_rz = -0.201240;   // AU
        double orbfit_vx = 0.011705;    // AU/d
        double orbfit_vy = -0.002918;   // AU/d
        double orbfit_vz = 0.001512;    // AU/d
        
        std::cout << "OrbFit @ MJD " << target_mjd << " (Equatoriale J2000):\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << orbfit_rx << ", " 
                  << orbfit_ry << ", " 
                  << orbfit_rz << ") AU\n";
        std::cout << "  v = ("
                  << orbfit_vx << ", "
                  << orbfit_vy << ", "
                  << orbfit_vz << ") AU/d\n";
        double r_orbfit_mag = std::sqrt(orbfit_rx*orbfit_rx + orbfit_ry*orbfit_ry + orbfit_rz*orbfit_rz);
        std::cout << "  |r| = " << std::setprecision(6) << r_orbfit_mag << " AU\n\n";
        
        // ⚠️ CONVERSIONE FRAME NECESSARIA!
        // IOccultCalc: Eclittica J2000 (ECLIPJ2000_DE441)
        // OrbFit: Equatoriale J2000 (ICRF/J2000)
        // Conversione: obliquità J2000 = 23.439291111°
        double eps = 23.439291111 * M_PI / 180.0;
        
        // Eclittica → Equatoriale (rotazione intorno X)
        double r_ioc_eq_x = r_ioc.x;
        double r_ioc_eq_y = r_ioc.y * cos(eps) - r_ioc.z * sin(eps);
        double r_ioc_eq_z = r_ioc.y * sin(eps) + r_ioc.z * cos(eps);
        
        double v_ioc_eq_x = v_ioc.x;
        double v_ioc_eq_y = v_ioc.y * cos(eps) - v_ioc.z * sin(eps);
        double v_ioc_eq_z = v_ioc.y * sin(eps) + v_ioc.z * cos(eps);
        
        std::cout << "IOccultCalc @ MJD " << target_mjd << " (Equatoriale J2000 convertito):\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << r_ioc_eq_x << ", " 
                  << r_ioc_eq_y << ", " 
                  << r_ioc_eq_z << ") AU\n";
        std::cout << "  v = ("
                  << v_ioc_eq_x << ", "
                  << v_ioc_eq_y << ", "
                  << v_ioc_eq_z << ") AU/d\n\n";
        
        // Differenza IOccultCalc vs OrbFit (in frame equatoriale)
        double dx = (r_ioc_eq_x - orbfit_rx) * 149597870.7;  // AU → km
        double dy = (r_ioc_eq_y - orbfit_ry) * 149597870.7;
        double dz = (r_ioc_eq_z - orbfit_rz) * 149597870.7;
        double diff_pos = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        double dvx = (v_ioc_eq_x - orbfit_vx) * 149597870.7 / 86400.0;  // AU/d → m/s
        double dvy = (v_ioc_eq_y - orbfit_vy) * 149597870.7 / 86400.0;
        double dvz = (v_ioc_eq_z - orbfit_vz) * 149597870.7 / 86400.0;
        double diff_vel = std::sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
        
        std::cout << "--- DIFFERENZA IOccultCalc vs OrbFit ---\n";
        std::cout << "  ΔX = " << std::setprecision(1) << dx << " km (" 
                  << std::setprecision(1) << (dx/diff_pos*100) << "%)\n";
        std::cout << "  ΔY = " << dy << " km (" 
                  << (dy/diff_pos*100) << "%)\n";
        std::cout << "  ΔZ = " << dz << " km (" 
                  << (dz/diff_pos*100) << "%)\n";
        std::cout << "  |Δr| = " << diff_pos << " km\n";
        std::cout << "  |Δv| = " << diff_vel << " m/s\n\n";
        
        std::cout << "  Error/year = " << (diff_pos / std::abs(dt_days/365.25)) << " km/year\n\n";
        
        // ANALISI
        std::cout << "=================================================================\n";
        std::cout << "ANALISI RISULTATI\n";
        std::cout << "=================================================================\n\n";
        
        std::cout << "Confronto modelli:\n";
        std::cout << "  IOccultCalc:  9 pianeti + 16 asteroidi + GR (RK4, step=0.1d)\n";
        std::cout << "  OrbFit:       9 pianeti + 17 asteroidi + GR (RA15 multistep)\n\n";
        
        if (diff_pos < 1000) {
            std::cout << "✅ ECCELLENTE! Errore <1000 km dopo 136 anni!\n\n";
            std::cout << "IOccultCalc concorda con OrbFit a livello professionale.\n";
        } else if (diff_pos < 10000) {
            std::cout << "✓ BUONO. Errore ~" << (int)diff_pos << " km dopo 136 anni.\n\n";
            std::cout << "~" << (int)(diff_pos/(std::abs(dt_days)/365.25)) << " km/anno è ";
            std::cout << "accettabile per propagazione a lungo termine.\n\n";
            std::cout << "Differenze attese da:\n";
            std::cout << "  • Integratore: RK4 vs RA15 multistep\n";
            std::cout << "  • Asteroide mancante (#17 in OrbFit)\n";
            std::cout << "  • Step size: fisso 0.1d vs adattivo\n";
        } else if (diff_pos < 100000) {
            std::cout << "⚠ MODERATO. Errore ~" << (int)diff_pos << " km dopo 136 anni.\n\n";
            std::cout << "Possibili cause:\n";
            std::cout << "  • Integratore RK4 non ottimale per 136 anni\n";
            std::cout << "  • Step size troppo grande (0.1d)\n";
            std::cout << "  • 17° asteroide significativo per Eros\n";
        } else {
            std::cout << "❌ ERRORE ELEVATO: " << (int)diff_pos << " km dopo 136 anni!\n\n";
            std::cout << "Indica problemi con:\n";
            std::cout << "  1. Conversione frame eclittica↔equatoriale\n";
            std::cout << "  2. Modello forze incompleto\n";
            std::cout << "  3. Accumulo errore RK4 su 136 anni\n";
        }
        
        std::cout << "\nNote tecniche:\n";
        std::cout << "• Conversione frame applicata: Eclittica J2000 → Equatoriale J2000\n";
        std::cout << "• Obliquità J2000: 23.439291111°\n";
        std::cout << "• Propagazione su " << std::abs(dt_days) << " giorni (" 
                  << std::abs(dt_days)/365.25 << " anni)\n";
        std::cout << "• RK4 steps: ~" << (int)std::abs(dt_days/0.1) << " (ogni 0.1 giorni)\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n=================================================================\n";
    
    return 0;
}
