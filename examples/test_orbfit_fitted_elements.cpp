/**
 * @file test_orbfit_fitted_elements.cpp
 * @brief Confronto IOccultCalc vs OrbFit usando gli stessi elementi FITTED
 * 
 * Test:
 * 1. Elementi iniziali: CORRECTED da OrbFit dopo fit su 16119 osservazioni
 * 2. Epoca: MJD 55918.807 (Dec 23, 2011)
 * 3. Propagazione: 12 anni → MJD 60300.0 (Jan 01, 2024)
 * 4. Confronto diretto: stesso punto di partenza, stessa epoca target
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <ioccultcalc/orbit_propagator.h>
#include <ioccultcalc/types.h>

using namespace ioccultcalc;

int main() {
    std::cout << "=================================================================\n";
    std::cout << "TEST: IOccultCalc vs OrbFit con ELEMENTI FITTED\n";
    std::cout << "=================================================================\n\n";
    
    std::cout << "Setup test:\n";
    std::cout << "  • Elementi: OrbFit CORRECTED (fit su 16119 osservazioni)\n";
    std::cout << "  • Epoca iniziale: MJD 55918.807 (Dec 23, 2011)\n";
    std::cout << "  • Propagazione: 12 anni → MJD 60300.0 (Jan 01, 2024)\n";
    std::cout << "  • Confronto: Stesso punto partenza + stessa epoca target\n\n";
    
    try {
        // STEP 1: Elementi CORRECTED da OrbFit (dopo differential correction)
        std::cout << "--- STEP 1: Elementi iniziali FITTED da OrbFit ---\n";
        
        // Valori da 433.olg: "Corrected orbital elements"
        // Position vector @ MJD 55918.807 (Equatoriale J2000)
        double r0_eq_x = -2.23090281851654E-01;  // AU
        double r0_eq_y =  1.12387595809716E+00;  // AU
        double r0_eq_z =  8.60668959807618E-02;  // AU
        double v0_eq_x = -1.68320119058993E-02;  // AU/d
        double v0_eq_y = -4.36024765332142E-03;  // AU/d
        double v0_eq_z = -3.12870435158407E-03;  // AU/d
        
        double epoch0_mjd = 55918.807374655;
        JulianDate epoch0(epoch0_mjd + 2400000.5);
        
        std::cout << "Elementi OrbFit @ MJD " << std::fixed << std::setprecision(6) 
                  << epoch0_mjd << " (Equatoriale J2000):\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << r0_eq_x << ", "
                  << r0_eq_y << ", "
                  << r0_eq_z << ") AU\n";
        std::cout << "  v = ("
                  << v0_eq_x << ", "
                  << v0_eq_y << ", "
                  << v0_eq_z << ") AU/d\n";
        double r0_mag = std::sqrt(r0_eq_x*r0_eq_x + r0_eq_y*r0_eq_y + r0_eq_z*r0_eq_z);
        std::cout << "  |r| = " << std::setprecision(6) << r0_mag << " AU\n";
        std::cout << "  Residual norm: 0.326 arcsec (16119 osservazioni)\n\n";
        
        // Conversione Equatoriale → Eclittica (per IOccultCalc)
        double eps = 23.439291111 * M_PI / 180.0;  // Obliquità J2000
        
        double r0_ecl_x = r0_eq_x;
        double r0_ecl_y = r0_eq_y * cos(eps) + r0_eq_z * sin(eps);
        double r0_ecl_z = -r0_eq_y * sin(eps) + r0_eq_z * cos(eps);
        
        double v0_ecl_x = v0_eq_x;
        double v0_ecl_y = v0_eq_y * cos(eps) + v0_eq_z * sin(eps);
        double v0_ecl_z = -v0_eq_y * sin(eps) + v0_eq_z * cos(eps);
        
        std::cout << "Conversione a Eclittica J2000 (per IOccultCalc):\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << r0_ecl_x << ", "
                  << r0_ecl_y << ", "
                  << r0_ecl_z << ") AU\n";
        std::cout << "  v = ("
                  << v0_ecl_x << ", "
                  << v0_ecl_y << ", "
                  << v0_ecl_z << ") AU/d\n\n";
        
        // STEP 2: Propaga con IOccultCalc
        std::cout << "--- STEP 2: Propagazione IOccultCalc ---\n";
        
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RK4;
        opts.stepSize = 0.05;  // Ridotto da 0.1 a 0.05 per maggiore precisione
        opts.usePlanetaryPerturbations = true;
        
        OrbitPropagator prop(opts);
        
        // Crea stato iniziale
        OrbitState state_0;
        state_0.position = {r0_ecl_x, r0_ecl_y, r0_ecl_z};
        state_0.velocity = {v0_ecl_x, v0_ecl_y, v0_ecl_z};
        state_0.epoch = epoch0;
        
        // Epoca target
        double target_mjd = 60300.0;  // 2024-Jan-01
        JulianDate target_epoch(target_mjd + 2400000.5);
        
        double dt_days = target_mjd - epoch0_mjd;
        std::cout << "Propagando " << std::setprecision(1) << dt_days << " giorni (~" 
                  << (dt_days/365.25) << " anni)...\n";
        std::cout << "(Con RK4 step=0.05d, ~" << (int)(dt_days/0.05) << " steps)\n";
        std::cout << "⏳ Tempo stimato: ~4-6 minuti (step ridotto)...\n\n";
        
        auto state_final = prop.propagate(state_0, target_epoch);
        
        Vector3D r_ecl = state_final.position;
        Vector3D v_ecl = state_final.velocity;
        
        std::cout << "IOccultCalc @ MJD " << target_mjd << " (Eclittica J2000):\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << r_ecl.x << ", "
                  << r_ecl.y << ", "
                  << r_ecl.z << ") AU\n";
        std::cout << "  v = ("
                  << v_ecl.x << ", "
                  << v_ecl.y << ", "
                  << v_ecl.z << ") AU/d\n";
        double r_mag = std::sqrt(r_ecl.x*r_ecl.x + r_ecl.y*r_ecl.y + r_ecl.z*r_ecl.z);
        std::cout << "  |r| = " << std::setprecision(6) << r_mag << " AU\n\n";
        
        // Conversione finale a Equatoriale per confronto
        double r_eq_x = r_ecl.x;
        double r_eq_y = r_ecl.y * cos(eps) - r_ecl.z * sin(eps);
        double r_eq_z = r_ecl.y * sin(eps) + r_ecl.z * cos(eps);
        
        double v_eq_x = v_ecl.x;
        double v_eq_y = v_ecl.y * cos(eps) - v_ecl.z * sin(eps);
        double v_eq_z = v_ecl.y * sin(eps) + v_ecl.z * cos(eps);
        
        std::cout << "IOccultCalc @ MJD " << target_mjd << " (Equatoriale J2000):\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << r_eq_x << ", "
                  << r_eq_y << ", "
                  << r_eq_z << ") AU\n";
        std::cout << "  v = ("
                  << v_eq_x << ", "
                  << v_eq_y << ", "
                  << v_eq_z << ") AU/d\n\n";
        
        // STEP 3: Confronto con OrbFit
        std::cout << "--- STEP 3: Confronto con OrbFit ---\n";
        
        // Valori OrbFit @ MJD 60300.0 (da 433.oel)
        double orbfit_rx = 1.39933263576780E+00;
        double orbfit_ry = 4.05376107746071E-01;
        double orbfit_rz = 2.64807319940301E-01;
        double orbfit_vx = -6.88234140936674E-03;
        double orbfit_vy = 1.22202237086118E-02;
        double orbfit_vz = 2.28870868886193E-04;
        
        std::cout << "OrbFit @ MJD " << target_mjd << " (Equatoriale J2000):\n";
        std::cout << "  r = (" << std::setprecision(12)
                  << orbfit_rx << ", "
                  << orbfit_ry << ", "
                  << orbfit_rz << ") AU\n";
        std::cout << "  v = ("
                  << orbfit_vx << ", "
                  << orbfit_vy << ", "
                  << orbfit_vz << ") AU/d\n\n";
        
        // Differenze
        double dx = (r_eq_x - orbfit_rx) * 149597870.7;  // AU → km
        double dy = (r_eq_y - orbfit_ry) * 149597870.7;
        double dz = (r_eq_z - orbfit_rz) * 149597870.7;
        double diff_pos = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        double dvx = (v_eq_x - orbfit_vx) * 149597870.7 / 86400.0;  // AU/d → m/s
        double dvy = (v_eq_y - orbfit_vy) * 149597870.7 / 86400.0;
        double dvz = (v_eq_z - orbfit_vz) * 149597870.7 / 86400.0;
        double diff_vel = std::sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
        
        std::cout << "=================================================================\n";
        std::cout << "DIFFERENZA IOccultCalc - OrbFit\n";
        std::cout << "=================================================================\n\n";
        
        std::cout << "  ΔX = " << std::setprecision(1) << dx << " km (" 
                  << std::abs(dx/diff_pos*100) << "%)\n";
        std::cout << "  ΔY = " << dy << " km (" 
                  << std::abs(dy/diff_pos*100) << "%)\n";
        std::cout << "  ΔZ = " << dz << " km (" 
                  << std::abs(dz/diff_pos*100) << "%)\n";
        std::cout << "  |Δr| = " << diff_pos << " km\n";
        std::cout << "  |Δv| = " << std::setprecision(3) << diff_vel << " m/s\n\n";
        
        std::cout << "  Error/year = " << std::setprecision(1) 
                  << (diff_pos / (dt_days/365.25)) << " km/year\n\n";
        
        // ANALISI
        std::cout << "=================================================================\n";
        std::cout << "ANALISI RISULTATI\n";
        std::cout << "=================================================================\n\n";
        
        std::cout << "Condizioni test:\n";
        std::cout << "  • STESSO punto di partenza (elementi fitted OrbFit)\n";
        std::cout << "  • STESSA epoca target (MJD 60300.0)\n";
        std::cout << "  • IOccultCalc: 16 asteroidi + RK4 step=0.1d\n";
        std::cout << "  • OrbFit: 17 asteroidi + RA15 multistep\n\n";
        
        if (diff_pos < 1000) {
            std::cout << "✅ ECCELLENTE! Errore <1000 km dopo 12 anni!\n\n";
            std::cout << "IOccultCalc e OrbFit concordano a livello professionale.\n";
            std::cout << "Le differenze sono dovute solo a:\n";
            std::cout << "  • Integratore: RK4 fisso vs RA15 adattivo\n";
            std::cout << "  • Asteroide #17 (mancante in IOccultCalc)\n";
        } else if (diff_pos < 10000) {
            std::cout << "✓ BUONO. Errore ~" << (int)diff_pos << " km dopo 12 anni.\n\n";
            std::cout << "~" << (int)(diff_pos/(dt_days/365.25)) << " km/anno è accettabile.\n\n";
            std::cout << "Differenze attese da:\n";
            std::cout << "  • Integratore: RK4 step fisso vs RA15 multistep adattivo\n";
            std::cout << "  • Asteroide #17 in OrbFit (mancante in IOccultCalc)\n";
            std::cout << "  • Step size: 0.1d fisso può accumulare errore\n";
        } else if (diff_pos < 100000) {
            std::cout << "⚠ MODERATO. Errore ~" << (int)diff_pos << " km dopo 12 anni.\n\n";
            std::cout << "Possibili cause:\n";
            std::cout << "  • RK4 non ottimale per 12 anni (prova step=0.05d)\n";
            std::cout << "  • Asteroide #17 significativo per Eros\n";
            std::cout << "  • Conversione frame potrebbe avere imprecisioni\n";
        } else {
            std::cout << "❌ ERRORE ELEVATO: " << (int)diff_pos << " km dopo 12 anni!\n\n";
            std::cout << "Problemi probabili:\n";
            std::cout << "  1. Conversione frame eclittica↔equatoriale errata\n";
            std::cout << "  2. Modello forze significativamente diverso\n";
            std::cout << "  3. Bug nel propagatore IOccultCalc\n";
        }
        
        std::cout << "\nNote tecniche:\n";
        std::cout << "• Conversione frame: Equatoriale ↔ Eclittica J2000\n";
        std::cout << "• Obliquità J2000: 23.439291111°\n";
        std::cout << "• Propagazione: " << dt_days << " giorni (" 
                  << (dt_days/365.25) << " anni)\n";
        std::cout << "• Questo test elimina ambiguità su elementi iniziali!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n=================================================================\n";
    
    return 0;
}
