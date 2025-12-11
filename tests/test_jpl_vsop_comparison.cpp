/**
 * @file test_jpl_vsop_comparison.cpp
 * @brief Confronto JPL DE441 vs VSOP87 per posizione Terra
 */

#include "ioccultcalc/spk_reader.h"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

using namespace ioccultcalc;

int main() {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  CONFRONTO JPL DE441 vs VSOP87 - POSIZIONE TERRA         ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    // Test date: 28 Nov 2025 12:00 UTC = MJD 61007.5 TDB
    double test_mjd = 61007.5;
    double jd_tdb = test_mjd + 2400000.5;
    
    std::cout << "Test date: MJD " << std::fixed << std::setprecision(5) << test_mjd 
              << " (JD " << jd_tdb << " TDB)\n\n";
    
    // 1. VSOP87 (AstDyn PlanetaryEphemeris)
    Eigen::Vector3d earth_vsop_helio_ecl = 
        astdyn::ephemeris::PlanetaryEphemeris::getPosition(
            astdyn::ephemeris::CelestialBody::EARTH, jd_tdb);
    
    Eigen::Vector3d sun_vsop_bary_ecl = 
        astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
    
    Eigen::Vector3d earth_vsop_bary_ecl = earth_vsop_helio_ecl - sun_vsop_bary_ecl;
    
    std::cout << "1. VSOP87 (AstDyn PlanetaryEphemeris):\n";
    std::cout << "   Terra baricentrica eclittica: [" << std::setprecision(10)
              << earth_vsop_bary_ecl[0] << ", " 
              << earth_vsop_bary_ecl[1] << ", " 
              << earth_vsop_bary_ecl[2] << "] AU\n\n";
    
    // 2. JPL DE441 (SPKReader)
    SPKReader spk;
    bool loaded = spk.ensureFileLoaded("DE441");
    
    if (!loaded) {
        std::cerr << "⚠ DE441 non disponibile, provo DE430...\n";
        loaded = spk.ensureFileLoaded("DE430");
    }
    
    if (loaded) {
        std::cout << "2. JPL " << spk.getVersion() << " (SPKReader):\n";
        
        // NAIF ID 3 = Earth
        // jpl_pleph: 12 = SSB (Solar System Barycenter)
        Vector3D earth_jpl_bary_ecl_vec;
        try {
            earth_jpl_bary_ecl_vec = spk.getPosition(3, jd_tdb, 12);  // SSB-centric (jpl_pleph ID 12)
        } catch (const std::exception& e) {
            std::cerr << "Errore lettura SPK con centro 12: " << e.what() << "\n";
            std::cerr << "Provando con centro 11 (Sun) e calcolo manuale...\n";
            // Fallback: usa centro 11 (Sun) e calcola baricentrico manualmente
            Vector3D earth_jpl_helio_ecl = spk.getPosition(3, jd_tdb, 11);
            // Per ottenere baricentrico, dobbiamo sottrarre la posizione del Sole rispetto a SSB
            Vector3D sun_bary_ecl = spk.getPosition(11, jd_tdb, 12);  // Sole rispetto a SSB
            earth_jpl_bary_ecl_vec.x = earth_jpl_helio_ecl.x - sun_bary_ecl.x;
            earth_jpl_bary_ecl_vec.y = earth_jpl_helio_ecl.y - sun_bary_ecl.y;
            earth_jpl_bary_ecl_vec.z = earth_jpl_helio_ecl.z - sun_bary_ecl.z;
        }
        
        Eigen::Vector3d earth_jpl_bary_ecl;
        earth_jpl_bary_ecl[0] = earth_jpl_bary_ecl_vec.x;
        earth_jpl_bary_ecl[1] = earth_jpl_bary_ecl_vec.y;
        earth_jpl_bary_ecl[2] = earth_jpl_bary_ecl_vec.z;
        
        std::cout << "   Terra baricentrica eclittica: [" << std::setprecision(10)
                  << earth_jpl_bary_ecl[0] << ", " 
                  << earth_jpl_bary_ecl[1] << ", " 
                  << earth_jpl_bary_ecl[2] << "] AU\n\n";
        
        // Confronto
        Eigen::Vector3d diff;
        diff[0] = earth_vsop_bary_ecl[0] - earth_jpl_bary_ecl[0];
        diff[1] = earth_vsop_bary_ecl[1] - earth_jpl_bary_ecl[1];
        diff[2] = earth_vsop_bary_ecl[2] - earth_jpl_bary_ecl[2];
        double diff_mag = diff.norm();
        
        std::cout << "3. CONFRONTO:\n";
        std::cout << "   Differenza: [" << diff[0] << ", " << diff[1] << ", " << diff[2] << "] AU\n";
        std::cout << "   Magnitudine: " << diff_mag << " AU = " 
                  << diff_mag * 1.495978707e8 << " km\n";
        
        if (diff_mag < 1e-6) {
            std::cout << "   ✓ COERENTI\n";
        } else if (diff_mag < 1e-4) {
            std::cout << "   ⚠ Piccola differenza (accettabile)\n";
            std::cout << "   Questo potrebbe spiegare parte della differenza con Horizons\n";
        } else {
            std::cout << "   ⚠ DIFFERENZA SIGNIFICATIVA\n";
            std::cout << "   Questo spiega la differenza con Horizons\n";
        }
        
        // Calcola impatto su RA/Dec
        std::cout << "\n4. IMPATTO SU RA/DEC:\n";
        std::cout << "   La differenza nella posizione Terra di " << diff_mag * 1.495978707e8 
                  << " km\n";
        std::cout << "   può causare una differenza angolare di ~" 
                  << std::atan(diff_mag / 2.3) * 180.0 / M_PI * 3600.0 
                  << " arcsec\n";
        std::cout << "   (assumendo distanza asteroide ~2.3 AU)\n";
        
    } else {
        std::cerr << "❌ JPL DE non disponibile\n";
        std::cerr << "   Installare DE441 o DE430 per confronto\n";
    }
    
    return 0;
}

