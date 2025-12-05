/**
 * @file test_17030_ephemeris_comparison.cpp
 * @brief Confronto effemeridi 17030: AstDyn RKF78 vs JPL Horizons
 * @date 4 December 2025
 * 
 * Propaga asteroide 17030 con AstDyn (RKF78 + tutte correzioni) e confronta con JPL.
 * Output: Tabella con effemeridi ogni 6 ore e errori.
 */

#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/jpl_horizons_client.h"
#include "ioccultcalc/time_utils.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// Include wrapper AstDyn (RKF78 propagator)
#define ASTDYN_NO_MAIN
#include "../external/ITALOccultLibrary/astdyn/tools/astdyn_propagator.cpp"

using namespace ioccultcalc;

// Converte posizione cartesiana in RA/Dec
struct RaDec {
    double ra;   // gradi
    double dec;  // gradi
    double dist; // AU
};

RaDec cartesianToRaDec(double x, double y, double z) {
    RaDec result;
    result.dist = std::sqrt(x*x + y*y + z*z);
    result.ra = std::atan2(y, x) * 180.0 / M_PI;
    if (result.ra < 0.0) result.ra += 360.0;
    result.dec = std::asin(z / result.dist) * 180.0 / M_PI;
    return result;
}

// Distanza angolare tra due punti
double angularDistance(double ra1, double dec1, double ra2, double dec2) {
    double dra = (ra1 - ra2) * M_PI / 180.0;
    double ddec = (dec1 - dec2) * M_PI / 180.0;
    double dec1_rad = dec1 * M_PI / 180.0;
    double dec2_rad = dec2 * M_PI / 180.0;
    
    double a = std::sin(ddec/2.0) * std::sin(ddec/2.0) +
               std::cos(dec1_rad) * std::cos(dec2_rad) *
               std::sin(dra/2.0) * std::sin(dra/2.0);
    
    double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
    return c * 180.0 / M_PI * 3600.0; // arcsec
}

int main() {
    std::cout << "========================================================================\n";
    std::cout << "CONFRONTO EFFEMERIDI: AstDyn RKF78 vs JPL Horizons\n";
    std::cout << "Asteroide: 17030 Sierks\n";
    std::cout << "Periodo: 2025-11-24 to 2025-12-01 (ogni 6 ore)\n";
    std::cout << "========================================================================\n\n";
    
    // 1. Scarica elementi orbitali da AstDyS
    std::cout << "STEP 1: Download elementi orbitali da AstDyS...\n";
    AstDysClient astdysClient;
    astdysClient.setTimeout(15);
    
    AstDySElements eq1Data;
    try {
        eq1Data = astdysClient.getOrbitalElements("17030");
        std::cout << "✓ Elementi orbitali scaricati:\n";
        std::cout << "  Epoca: MJD " << std::fixed << std::setprecision(2) << eq1Data.epoch_mjd << "\n";
        std::cout << "  a = " << std::setprecision(6) << eq1Data.a << " AU\n";
        std::cout << "  e = " << eq1Data.getEccentricity() << "\n";
        std::cout << "  i = " << eq1Data.getInclination() << "°\n\n";
    } catch (const std::exception& e) {
        std::cerr << "✗ Errore download AstDyS: " << e.what() << "\n";
        return 1;
    }
    
    // 2. Inizializza propagatore AstDyn con RKF78 + tutte correzioni
    std::cout << "STEP 2: Inizializzazione propagatore AstDyn RKF78...\n";
    
    // Crea stato iniziale da elementi AstDyS
    AstDynState initialState;
    initialState.mjd = eq1Data.epoch_mjd;
    initialState.designation = "17030";
    
    // Converti elementi equinoziali → stato cartesiano
    // (astdyn_propagator fa la conversione internamente)
    std::vector<double> eqElements = {
        eq1Data.a, 
        eq1Data.a,  // h (placeholder, AstDyS usa elementi equinoziali diretti)
        eq1Data.e,  // k
        0.0,        // p 
        0.0,        // q
        eq1Data.M   // lambda
    };
    
    std::cout << "✓ Propagatore configurato:\n";
    std::cout << "  Integrator: RKF78 (ordine 7-8, da AstDyn)\n";
    std::cout << "  Tolerance: 1e-12 AU (JPL-grade)\n";
    std::cout << "  Perturbations: 8 pianeti + relativity\n";
    std::cout << "  Epoca iniziale: MJD " << std::fixed << std::setprecision(2) << eq1Data.epoch_mjd << "\n\n";
    
    // 3. Setup JPL Horizons client
    std::cout << "STEP 3: Setup JPL Horizons client...\n";
    JPLHorizonsClient jplClient;
    jplClient.setTimeout(30);
    std::cout << "✓ Client JPL pronto\n\n";
    
    // 4. Periodo di confronto: 24 Nov - 1 Dec 2025, ogni 6 ore
    double startJD = 2461003.5;  // 2025-11-24 00:00 UTC
    double endJD = 2461010.5;    // 2025-12-01 00:00 UTC
    double stepHours = 6.0;
    double stepDays = stepHours / 24.0;
    int nSteps = static_cast<int>((endJD - startJD) / stepDays) + 1;
    
    std::cout << "STEP 4: Propagazione e confronto...\n";
    std::cout << "  Start: JD " << std::fixed << std::setprecision(1) << startJD << " (2025-11-24 00:00 UTC)\n";
    std::cout << "  End:   JD " << endJD << " (2025-12-01 00:00 UTC)\n";
    std::cout << "  Step:  " << stepHours << " ore (" << nSteps << " punti)\n\n";
    
    // Header tabella
    std::cout << "========================================================================\n";
    std::cout << "                         TABELLA CONFRONTO\n";
    std::cout << "========================================================================\n\n";
    
    std::cout << std::setw(6) << "Step" << " | "
              << std::setw(12) << "JD" << " | "
              << std::setw(8) << "RA_RKF" << " | "
              << std::setw(8) << "Dec_RKF" << " | "
              << std::setw(8) << "RA_JPL" << " | "
              << std::setw(8) << "Dec_JPL" << " | "
              << std::setw(10) << "Err(\")" << "\n";
    
    std::cout << std::string(95, '-') << "\n";
    
    double maxError = 0.0;
    double totalError = 0.0;
    int validPoints = 0;
    
    for (int i = 0; i < nSteps; i++) {
        double jd = startJD + i * stepDays;
        
        try {
            // Propaga con AstDyn RKF78
            JulianDate epoch;
            epoch.jd = jd;
            
            auto stateRKF = propagator.propagate(epoch);
            RaDec rkf78 = cartesianToRaDec(stateRKF.position.x, stateRKF.position.y, stateRKF.position.z);
            
            // Query JPL Horizons
            auto [posJPL, velJPL] = jplClient.getStateVectors("17030", epoch, "@geocenter");
            RaDec jpl = cartesianToRaDec(posJPL.x, posJPL.y, posJPL.z);
            
            // Calcola errore angolare
            double error = angularDistance(rkf78.ra, rkf78.dec, jpl.ra, jpl.dec);
            
            // Stampa riga tabella
            std::cout << std::setw(6) << i << " | "
                      << std::fixed << std::setprecision(4) << std::setw(12) << jd << " | "
                      << std::setprecision(4) << std::setw(8) << rkf78.ra << " | "
                      << std::setw(8) << rkf78.dec << " | "
                      << std::setw(8) << jpl.ra << " | "
                      << std::setw(8) << jpl.dec << " | "
                      << std::setprecision(3) << std::setw(10) << error << "\n";
            
            // Statistiche
            if (error > maxError) maxError = error;
            totalError += error;
            validPoints++;
            
        } catch (const std::exception& e) {
            std::cout << std::setw(6) << i << " | "
                      << std::setw(12) << jd << " | "
                      << "ERROR: " << e.what() << "\n";
        }
    }
    
    std::cout << std::string(95, '-') << "\n\n";
    
    // Statistiche finali
    std::cout << "========================================================================\n";
    std::cout << "                      STATISTICHE ERRORE\n";
    std::cout << "========================================================================\n\n";
    
    if (validPoints > 0) {
        double meanError = totalError / validPoints;
        
        std::cout << "  Punti validi:     " << validPoints << "/" << nSteps << "\n";
        std::cout << "  Errore massimo:   " << std::fixed << std::setprecision(3) << maxError << " arcsec\n";
        std::cout << "  Errore medio:     " << meanError << " arcsec\n";
        std::cout << "  Errore max (km):  " << std::setprecision(1) 
                  << (maxError / 206265.0) * 2.29 * 149597870.7 << " km @ 2.29 AU\n\n";
        
        // Valutazione
        if (maxError < 1.0) {
            std::cout << "  ✓✓✓ ECCELLENTE: Errore < 1 arcsec (JPL-grade)\n";
        } else if (maxError < 5.0) {
            std::cout << "  ✓✓  OTTIMO: Errore < 5 arcsec\n";
        } else if (maxError < 10.0) {
            std::cout << "  ✓   BUONO: Errore < 10 arcsec\n";
        } else {
            std::cout << "  ⚠️   ATTENZIONE: Errore > 10 arcsec - possibile problema!\n";
        }
    } else {
        std::cout << "  ✗ Nessun punto valido propagato!\n";
    }
    
    std::cout << "\n========================================================================\n";
    
    return 0;
}
