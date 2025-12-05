/**
 * @file test_astdyn_rkf78_ephemeris.cpp
 * @brief Tabella effemeridi 17030 con AstDyn RKF78 vs JPL Horizons
 */

#include "ioccultcalc/astdyn_interface.h"
#include "ioccultcalc/jpl_horizons_client.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

// Distanza angolare
double angDist(double ra1, double dec1, double ra2, double dec2) {
    double d_ra = (ra1-ra2) * M_PI/180.0;
    double d_dec = (dec1-dec2) * M_PI/180.0;
    double a = sin(d_dec/2)*sin(d_dec/2) + cos(dec1*M_PI/180)*cos(dec2*M_PI/180)*sin(d_ra/2)*sin(d_ra/2);
    return 2*atan2(sqrt(a),sqrt(1-a)) * 180/M_PI * 3600; // arcsec
}

int main() {
    std::cout << "================================================================\n";
    std::cout << "TABELLA EFFEMERIDI: AstDyn RKF78 vs JPL Horizons\n";
    std::cout << "Asteroide 17030 - Periodo: 24 Nov - 1 Dic 2025 (ogni 6h)\n";
    std::cout << "================================================================\n\n";
    
    // 1. Carica elementi da file .eq1
    std::cout << "Caricamento elementi orbitali da 17030_astdys.eq1...\n";
    AstDynPropagator prop;
    
    if (!prop.loadFromEQ1File("17030_astdys.eq1")) {
        std::cerr << "✗ Errore: file 17030_astdys.eq1 non trovato\n";
        std::cerr << "Scaricalo da: https://newton.spacedys.com/~astdys2/astdys/17030.eq1\n";
        return 1;
    }
    
    std::cout << "✓ Elementi orbitali caricati\n";
    std::cout << "  a = " << prop.getOrbitalElements().a << " AU\n";
    std::cout << "  e = " << prop.getOrbitalElements().e << "\n";
    std::cout << "  Epoca: MJD " << prop.getOrbitalElements().epoch_mjd << "\n\n";
    
    // 2. Setup JPL client
    JPLHorizonsClient jpl;
    jpl.setTimeout(30);
    
    // 3. Periodo: 24 Nov - 1 Dic, ogni 6h
    double start_mjd = 61003.0;  // 2025-11-24 00:00 UTC
    double end_mjd = 61010.0;    // 2025-12-01 00:00 UTC
    double step_hours = 6.0;
    double step_days = step_hours / 24.0;
    
    std::cout << "Propagazione periodo:\n";
    std::cout << "  Start: MJD " << std::fixed << std::setprecision(1) << start_mjd << "\n";
    std::cout << "  End:   MJD " << end_mjd << "\n";
    std::cout << "  Step:  " << step_hours << " ore\n\n";
    
    // Header tabella
    std::cout << "================================================================\n";
    std::cout << std::setw(5) << "Step" << " | "
              << std::setw(10) << "MJD" << " | "
              << std::setw(10) << "RA_AstDyn" << " | "
              << std::setw(10) << "Dec_AstDyn" << " | "
              << std::setw(10) << "RA_JPL" << " | "
              << std::setw(10) << "Dec_JPL" << " | "
              << std::setw(10) << "Error(\")" << "\n";
    std::cout << std::string(85, '-') << "\n";
    
    double max_err = 0, sum_err = 0;
    int n = 0;
    
    for (double mjd = start_mjd; mjd <= end_mjd; mjd += step_days) {
        try {
            // Propaga con AstDyn RKF78
            CartesianStateICRF state = prop.propagateToMJD(mjd);
            
            // Converti in RA/Dec
            double r = sqrt(state.position.x*state.position.x + 
                           state.position.y*state.position.y + 
                           state.position.z*state.position.z);
            double ra_astdyn = atan2(state.position.y, state.position.x) * 180/M_PI;
            if (ra_astdyn < 0) ra_astdyn += 360;
            double dec_astdyn = asin(state.position.z / r) * 180/M_PI;
            
            // Query JPL
            JulianDate epoch;
            epoch.jd = mjd + 2400000.5;
            auto [pos_jpl, vel_jpl] = jpl.getStateVectors("17030", epoch, "@geocenter");
            
            double r_jpl = sqrt(pos_jpl.x*pos_jpl.x + pos_jpl.y*pos_jpl.y + pos_jpl.z*pos_jpl.z);
            double ra_jpl = atan2(pos_jpl.y, pos_jpl.x) * 180/M_PI;
            if (ra_jpl < 0) ra_jpl += 360;
            double dec_jpl = asin(pos_jpl.z / r_jpl) * 180/M_PI;
            
            // Errore
            double err = angDist(ra_astdyn, dec_astdyn, ra_jpl, dec_jpl);
            
            // Stampa
            std::cout << std::setw(5) << n << " | "
                      << std::fixed << std::setprecision(4) << std::setw(10) << mjd << " | "
                      << std::setprecision(6) << std::setw(10) << ra_astdyn << " | "
                      << std::setw(10) << dec_astdyn << " | "
                      << std::setw(10) << ra_jpl << " | "
                      << std::setw(10) << dec_jpl << " | "
                      << std::setprecision(3) << std::setw(10) << err << "\n";
            
            if (err > max_err) max_err = err;
            sum_err += err;
            n++;
            
        } catch (const std::exception& e) {
            std::cout << std::setw(5) << n << " | " << std::setw(10) << mjd 
                      << " | ERROR: " << e.what() << "\n";
        }
    }
    
    std::cout << std::string(85, '-') << "\n\n";
    
    // Statistiche
    if (n > 0) {
        std::cout << "STATISTICHE:\n";
        std::cout << "  Punti: " << n << "\n";
        std::cout << "  Errore max:   " << std::fixed << std::setprecision(3) << max_err << " arcsec\n";
        std::cout << "  Errore medio: " << (sum_err/n) << " arcsec\n";
        
        if (max_err < 1.0) {
            std::cout << "\n✓✓✓ ECCELLENTE: JPL-grade precision\n";
        } else if (max_err < 5.0) {
            std::cout << "\n✓✓ OTTIMO\n";
        } else {
            std::cout << "\n⚠️ Errore elevato\n";
        }
    }
    
    std::cout << "\n================================================================\n";
    return 0;
}
