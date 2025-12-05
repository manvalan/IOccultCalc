/**
 * @file test_17030_closeapproach.cpp
 * @brief Test close approach 17030 vs Gaia 3411546266140512128 - 28 Nov 2025
 */

#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/occultation_predictor.h"
#include "ioccultcalc/astdys_client.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

// Conversione arcsec in gradi
double arcsecToDegrees(double arcsec) {
    return arcsec / 3600.0;
}

// Distanza angolare tra due punti sulla sfera celeste
double angularDistance(double ra1, double dec1, double ra2, double dec2) {
    double d_ra = (ra1 - ra2) * M_PI / 180.0;
    double d_dec = (dec1 - dec2) * M_PI / 180.0;
    double a = std::sin(d_dec / 2.0) * std::sin(d_dec / 2.0) +
               std::cos(dec1 * M_PI / 180.0) * std::cos(dec2 * M_PI / 180.0) *
               std::sin(d_ra / 2.0) * std::sin(d_ra / 2.0);
    double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
    return c * 180.0 / M_PI * 3600.0; // ritorna arcsec
}

int main() {
    std::cout << "==================================================================\n";
    std::cout << "TEST CLOSE APPROACH: Asteroid 17030 vs Gaia 3411546266140512128\n";
    std::cout << "                     28 November 2025\n";
    std::cout << "==================================================================\n\n";
    
    // Stella target (da XML occultazione)
    double star_ra = 4.89440726;   // gradi
    double star_dec = 20.3316615;  // gradi
    double star_mag = 13.36;
    std::string gaia_id = "3411546266140512128";
    
    std::cout << "Stella target:\n";
    std::cout << "  Gaia DR3 ID: " << gaia_id << "\n";
    std::cout << "  RA:  " << std::fixed << std::setprecision(8) << star_ra << "° (J2000)\n";
    std::cout << "  Dec: " << star_dec << "° (J2000)\n";
    std::cout << "  Mag: " << std::setprecision(2) << star_mag << "\n\n";
    
    // Epoca evento
    double jd_event = 2460994.9306;  // 2025-11-28 10:20 UTC
    
    std::cout << "Epoca close approach:\n";
    std::cout << "  JD:  " << std::setprecision(4) << jd_event << "\n";
    std::cout << "  MJD: " << (jd_event - 2400000.5) << "\n";
    std::cout << "  Data: 2025-11-28 10:20:03 UTC\n\n";
    
    // Carica elementi orbitali da AstDyS
    std::cout << "Caricamento elementi orbitali 17030 da AstDyS...\n";
    AstDysClient client;
    client.setTimeout(10);
    
    try {
        auto eq1Data = client.getOrbitalElements("17030");
        
        std::cout << "✓ Elementi orbitali scaricati:\n";
        std::cout << "  Epoca: MJD " << std::setprecision(2) << eq1Data.epoch_mjd << "\n";
        std::cout << "  a = " << std::setprecision(6) << eq1Data.a << " AU\n";
        std::cout << "  e = " << eq1Data.getEccentricity() << "\n";
        std::cout << "  i = " << eq1Data.getInclination() << "°\n\n";
        
        // Converti in OrbitalElements
        OrbitalElements elements = eq1Data.toOrbitalElements();
        elements.designation = "17030";
        
        // Configura propagatore con massima precisione
        PropagationSettings settings;
        settings.integrator = Integrator::RKF78;
        settings.stepSize = 0.01;  // 14.4 minuti
        settings.tolerance = 1e-12;
        settings.usePlanetaryPerturbations = true;
        settings.numPerturbingBodies = 8;  // Tutti i pianeti
        settings.useRelativisticCorrections = true;
        
        OrbitPropagator propagator(elements, settings);
        
        std::cout << "Propagatore configurato:\n";
        std::cout << "  Integrator: RKF78 (ordine 7-8)\n";
        std::cout << "  Tolerance: 1e-12 AU\n";
        std::cout << "  Perturbations: 8 pianeti + relativity\n\n";
        
        // Propaga alla data dell'evento
        std::cout << "Propagazione all'epoca " << std::setprecision(4) << jd_event << "...\n";
        
        JulianDate epoch;
        epoch.jd = jd_event;
        
        auto state = propagator.propagate(epoch);
        
        // Converti stato in RA/Dec
        double r = std::sqrt(state.position.x * state.position.x +
                            state.position.y * state.position.y +
                            state.position.z * state.position.z);
        
        double ast_ra = std::atan2(state.position.y, state.position.x) * 180.0 / M_PI;
        if (ast_ra < 0.0) ast_ra += 360.0;
        
        double ast_dec = std::asin(state.position.z / r) * 180.0 / M_PI;
        
        std::cout << "\n✓ Posizione asteroide propagata:\n";
        std::cout << "  RA:  " << std::fixed << std::setprecision(8) << ast_ra << "° (J2000)\n";
        std::cout << "  Dec: " << ast_dec << "° (J2000)\n";
        std::cout << "  Distanza Terra: " << std::setprecision(6) << r << " AU\n\n";
        
        // Calcola distanza angolare
        double distance_arcsec = angularDistance(ast_ra, ast_dec, star_ra, star_dec);
        double distance_deg = distance_arcsec / 3600.0;
        
        std::cout << "==================================================================\n";
        std::cout << "RISULTATO CLOSE APPROACH:\n";
        std::cout << "==================================================================\n";
        std::cout << "  Distanza angolare: " << std::fixed << std::setprecision(2) 
                  << distance_arcsec << "\" (" << distance_deg << "°)\n";
        std::cout << "  ΔRA:  " << std::setprecision(8) << (ast_ra - star_ra) * 3600.0 << "\"\n";
        std::cout << "  ΔDec: " << (ast_dec - star_dec) * 3600.0 << "\"\n";
        
        // Diametro asteroide
        double diameter_km = 7.3;  // km (da XML)
        double distance_km = r * 149597870.7;  // AU → km
        double angular_diameter_arcsec = (diameter_km / distance_km) * 206265.0;
        
        std::cout << "\n  Diametro asteroide: " << diameter_km << " km\n";
        std::cout << "  Diametro angolare: " << std::setprecision(3) 
                  << angular_diameter_arcsec << "\"\n";
        
        if (distance_arcsec < angular_diameter_arcsec) {
            std::cout << "\n  ✓✓✓ OCCULTAZIONE CONFERMATA! ✓✓✓\n";
            std::cout << "  Stella dentro disco asteroide!\n";
        } else if (distance_arcsec < 2.0 * angular_diameter_arcsec) {
            std::cout << "\n  ⚠️  OCCULTAZIONE POSSIBILE (entro 2σ)\n";
        } else {
            std::cout << "\n  ✗ Nessuna occultazione (distanza > diametro)\n";
            std::cout << "  Rapporto: " << std::setprecision(1) 
                      << (distance_arcsec / angular_diameter_arcsec) << "× diametro\n";
        }
        
        std::cout << "==================================================================\n";
        
    } catch (const std::exception& e) {
        std::cerr << "✗ Errore: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
