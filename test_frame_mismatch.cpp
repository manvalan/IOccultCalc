/**
 * @file test_frame_mismatch.cpp
 * @brief Test CRITICO: verifica differenza frame ECLIPJ2000_DE405 vs ICRF
 * 
 * Ipotesi: L'errore di 1943 km su Z è dovuto a mismatch tra:
 * - Horizons: ECLIPTIC + ICRF (moderno)
 * - IOccultCalc/SPICE: ECLIPJ2000_DE405 (vecchio DE405)
 */

#include <ioccultcalc/jpl_horizons_client.h>
#include <iostream>
#include <iomanip>

// SPICE include
extern "C" {
    #include "SpiceUsr.h"
}

using namespace ioccultcalc;

void printVector(const std::string& label, double x, double y, double z) {
    std::cout << label << ": ("
              << std::fixed << std::setprecision(8)
              << x << ", " << y << ", " << z << ")\n";
}

int main() {
    std::cout << "\n=================================================================\n";
    std::cout << "TEST CRITICO: Frame ECLIPJ2000_DE405 vs J2000 (ICRF)\n";
    std::cout << "=================================================================\n\n";
    
    // Carica kernel SPICE
    std::string home = getenv("HOME");
    std::string spkPath = home + "/.ioccultcalc/ephemerides/linux_m13000p17000.441";
    std::string framePath = home + "/.ioccultcalc/ephemerides/eclipj2000_de405.tf";
    
    furnsh_c(spkPath.c_str());
    furnsh_c(framePath.c_str());
    
    if (failed_c()) {
        char msg[1841];
        getmsg_c("SHORT", 1840, msg);
        std::cerr << "SPICE load error: " << msg << std::endl;
        return 1;
    }
    
    std::cout << "✓ Kernel SPICE caricati\n\n";
    
    // Test: posizione Terra @ J2000
    double et = 0.0;  // J2000.0 epoch
    double state_ecl[6], state_j2000[6];
    double lt;
    
    std::cout << "Test: Posizione Terra @ J2000.0\n\n";
    
    // Frame 1: ECLIPJ2000_DE405 (quello che usa IOccultCalc)
    spkezr_c("399", et, "ECLIPJ2000_DE405", "NONE", "10", state_ecl, &lt);
    std::cout << "Frame ECLIPJ2000_DE405 (DE405, vecchio):\n";
    printVector("  pos", state_ecl[0], state_ecl[1], state_ecl[2]);
    std::cout << "  |r| = " << sqrt(state_ecl[0]*state_ecl[0] + 
                                     state_ecl[1]*state_ecl[1] + 
                                     state_ecl[2]*state_ecl[2]) << " km\n\n";
    
    // Frame 2: J2000 (ICRF equatoriale)
    spkezr_c("399", et, "J2000", "NONE", "10", state_j2000, &lt);
    std::cout << "Frame J2000 (ICRF equatoriale):\n";
    printVector("  pos", state_j2000[0], state_j2000[1], state_j2000[2]);
    std::cout << "  |r| = " << sqrt(state_j2000[0]*state_j2000[0] + 
                                     state_j2000[1]*state_j2000[1] + 
                                     state_j2000[2]*state_j2000[2]) << " km\n\n";
    
    // Conversione manuale J2000 → Ecliptic usando obliquità standard
    double eps = 23.43928 * 3.141592653589793 / 180.0;  // Obliquità J2000
    double cos_eps = cos(eps);
    double sin_eps = sin(eps);
    
    double x_ecl = state_j2000[0];
    double y_ecl = state_j2000[1] * cos_eps + state_j2000[2] * sin_eps;
    double z_ecl = -state_j2000[1] * sin_eps + state_j2000[2] * cos_eps;
    
    std::cout << "Frame J2000 convertito manualmente a Eclittico:\n";
    printVector("  pos", x_ecl, y_ecl, z_ecl);
    std::cout << "  |r| = " << sqrt(x_ecl*x_ecl + y_ecl*y_ecl + z_ecl*z_ecl) << " km\n\n";
    
    // Differenza
    std::cout << "=================================================================\n";
    std::cout << "DIFFERENZE\n";
    std::cout << "=================================================================\n\n";
    
    double dx = state_ecl[0] - x_ecl;
    double dy = state_ecl[1] - y_ecl;
    double dz = state_ecl[2] - z_ecl;
    
    std::cout << "ECLIPJ2000_DE405 vs J2000+conversione:\n";
    printVector("  Δr", dx, dy, dz);
    std::cout << "  |Δr| = " << sqrt(dx*dx + dy*dy + dz*dz) << " km\n\n";
    
    std::cout << "Componenti:\n";
    std::cout << "  ΔX = " << dx << " km (" << fabs(dx/(dx*dx+dy*dy+dz*dz)*100) << "%)\n";
    std::cout << "  ΔY = " << dy << " km (" << fabs(dy/(dx*dx+dy*dy+dz*dz)*100) << "%)\n";
    std::cout << "  ΔZ = " << dz << " km (" << fabs(dz/(dx*dx+dy*dy+dz*dz)*100) << "%)\n\n";
    
    // Test con Eros @ 2024
    std::cout << "=================================================================\n";
    std::cout << "Test con 433 Eros @ 2024-Jan-01\n";
    std::cout << "=================================================================\n\n";
    
    // Scarica da Horizons
    JPLHorizonsClient horizons;
    JulianDate epoch;
    epoch.jd = 2460615.5;
    
    auto [pos_hor, vel_hor] = horizons.getStateVectors("433", epoch);
    std::cout << "Horizons (ECLIPTIC+ICRF):\n";
    std::cout << "  pos = (" << pos_hor.x << ", " << pos_hor.y << ", " << pos_hor.z << ") AU\n\n";
    
    // Leggi da SPICE in ECLIPJ2000_DE405
    et = (epoch.jd - 2451545.0) * 86400.0;
    double state_eros[6];
    spkezr_c("2000433", et, "ECLIPJ2000_DE405", "NONE", "10", state_eros, &lt);
    
    if (!failed_c()) {
        std::cout << "SPICE (ECLIPJ2000_DE405):\n";
        std::cout << "  pos = (" << state_eros[0]/149597870.7 << ", " 
                  << state_eros[1]/149597870.7 << ", " 
                  << state_eros[2]/149597870.7 << ") AU\n\n";
        
        // Differenza
        dx = state_eros[0]/149597870.7 - pos_hor.x;
        dy = state_eros[1]/149597870.7 - pos_hor.y;
        dz = state_eros[2]/149597870.7 - pos_hor.z;
        
        std::cout << "Differenza SPICE vs Horizons:\n";
        std::cout << "  Δr = (" << dx << ", " << dy << ", " << dz << ") AU\n";
        std::cout << "     = (" << dx*149597870.7 << ", " 
                  << dy*149597870.7 << ", " 
                  << dz*149597870.7 << ") km\n";
        std::cout << "  |Δr| = " << sqrt(dx*dx+dy*dy+dz*dz)*149597870.7 << " km\n\n";
        
        if (fabs(dz*149597870.7) > 1000) {
            std::cout << "⚠️  ATTENZIONE: Differenza >1000 km su Z!\n";
            std::cout << "    Questo spiega l'errore osservato di 1943 km su Z.\n";
        }
    } else {
        std::cout << "✗ Eros non disponibile in SPICE kernel\n";
    }
    
    return 0;
}
