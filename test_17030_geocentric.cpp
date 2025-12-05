/**
 * @file test_17030_geocentric.cpp
 * @brief Effemeridi GEOCENTRICHE asteroide 17030 per 28/11/2025
 * @date 4 Dicembre 2025
 * 
 * Calcola coordinate GEOCENTRICHE (viste dalla Terra):
 * 1. Propaga asteroide (baricentrico ICRF)
 * 2. Ottieni posizione Terra (baricentrico ICRF) da AstDyn
 * 3. Calcola vettore geocentrico: asteroide - Terra
 * 4. Converti in RA/DEC geocentrici
 */

#include <chebyshev_rkf78_propagation.h>
#include <iostream>
#include <iomanip>
#include <cmath>

// AstDyn per effemeridi planetarie
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/ephemeris/PlanetaryData.hpp>

using namespace ioccultcalc;

// Converti RA da radianti a ore:minuti:secondi
void ra_to_hms(double ra_rad, int& h, int& m, double& s) {
    if (ra_rad < 0) ra_rad += 2.0 * M_PI;
    double ra_deg = ra_rad * 180.0 / M_PI;
    double ra_hours = ra_deg / 15.0;
    h = static_cast<int>(ra_hours);
    double m_frac = (ra_hours - h) * 60.0;
    m = static_cast<int>(m_frac);
    s = (m_frac - m) * 60.0;
}

// Converti DEC da radianti a gradi:arcmin:arcsec
void dec_to_dms(double dec_rad, char& sign, int& d, int& m, double& s) {
    double dec_deg = dec_rad * 180.0 / M_PI;
    sign = (dec_deg >= 0) ? '+' : '-';
    double dec_abs = std::abs(dec_deg);
    d = static_cast<int>(dec_abs);
    double m_frac = (dec_abs - d) * 60.0;
    m = static_cast<int>(m_frac);
    s = (m_frac - m) * 60.0;
}

int main() {
    try {
        std::cout << "============================================================\n";
        std::cout << "  EFFEMERIDI GEOCENTRICHE ASTEROIDE 17030\n";
        std::cout << "  28 NOVEMBRE 2025 - Viste dalla Terra\n";
        std::cout << "============================================================\n\n";
        
        // 1. Carica asteroide e crea propagatore
        std::cout << "1. Caricamento asteroide e creazione propagatore...\n";
        ChebyshevRKF78Propagator propagator("17030_astdys.eq1");
        std::cout << "   ✓ Propagatore RKF78 creato\n\n";
        
        // 2. Preparazione per effemeridi planetarie
        std::cout << "2. Sistema effemeridi Terra...\n";
        std::cout << "   ✓ Uso AstDyn PlanetaryEphemeris (VSOP87 analitico)\n";
        std::cout << "   ✓ Precisione VSOP87: ~1\" per Terra su secoli\n";
        std::cout << "   ✓ Frame: J2000 eclittico → conversione ICRF\n\n";
        
        // 3. Calcola posizioni per 28 Nov 2025
        double mjd_start = 60642.0;  // 28 Nov 2025 00:00 UTC
        double mjd_end = mjd_start + 1.0;
        int num_points = 25;  // ogni ora
        
        std::cout << "3. Propagazione asteroide (baricentrica)...\n";
        auto start_time = std::chrono::high_resolution_clock::now();
        
        std::vector<Eigen::Vector3d> asteroid_pos_bary = propagator.propagateForChebyshev(
            mjd_start, mjd_end, num_points);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
        
        std::cout << "   ✓ Propagazione completata in " << std::fixed << std::setprecision(1) 
                  << elapsed_ms << " ms\n\n";
        
        // 4. Calcola effemeridi geocentriche
        std::cout << "4. Calcolo coordinate geocentriche...\n\n";
        
        std::cout << "============================================================\n";
        std::cout << "  EFFEMERIDI GEOCENTRICHE - 28 NOVEMBRE 2025\n";
        std::cout << "============================================================\n\n";
        
        std::cout << " Ora UTC    RA (Geocentrico)   DEC (Geocentrico)  "
                  << " Delta(AU)    X_geo      Y_geo      Z_geo\n";
        std::cout << std::string(100, '-') << "\n";
        
        for (int i = 0; i < num_points; ++i) {
            double mjd = mjd_start + (i * (mjd_end - mjd_start) / (num_points - 1));
            int hour = static_cast<int>(std::round((mjd - mjd_start) * 24.0));
            
            // Posizione asteroide baricentrica ICRF
            const Eigen::Vector3d& ast_bary = asteroid_pos_bary[i];
            
            // Posizione Terra baricentrica (eclittico J2000 da VSOP87)
            double jd_tdb = mjd + 2400000.5;  // MJD → JD
            
            // Terra eliocentrica (frame eclittico J2000)
            Eigen::Vector3d earth_helio_ecl = astdyn::ephemeris::PlanetaryEphemeris::getPosition(
                astdyn::ephemeris::CelestialBody::EARTH, jd_tdb);
            
            // Sole baricentrico (frame eclittico J2000)
            Eigen::Vector3d sun_bary_ecl = astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
            
            // Terra baricentrica = Terra_helio - Sole_bary (in eclittico)
            Eigen::Vector3d earth_bary_ecl = earth_helio_ecl - sun_bary_ecl;
            
            // Converti Terra baricentrica da eclittico a ICRF (rotazione -ε attorno a X)
            constexpr double epsilon = 23.43928 * M_PI / 180.0;  // obliquità J2000
            double cos_eps = std::cos(epsilon);
            double sin_eps = std::sin(epsilon);
            Eigen::Vector3d earth_bary;
            earth_bary.x() = earth_bary_ecl.x();
            earth_bary.y() = earth_bary_ecl.y() * cos_eps - earth_bary_ecl.z() * sin_eps;
            earth_bary.z() = earth_bary_ecl.y() * sin_eps + earth_bary_ecl.z() * cos_eps;
            
            // Vettore geocentrico: asteroide visto dalla Terra
            Eigen::Vector3d ast_geo = ast_bary - earth_bary;
            
            // Distanza geocentrica
            double delta = ast_geo.norm();
            
            // Converti in RA/DEC geocentrici
            double ra_rad = std::atan2(ast_geo.y(), ast_geo.x());
            if (ra_rad < 0) ra_rad += 2.0 * M_PI;
            
            double r_xy = std::sqrt(ast_geo.x()*ast_geo.x() + ast_geo.y()*ast_geo.y());
            double dec_rad = std::atan2(ast_geo.z(), r_xy);
            
            // Formatta RA/DEC
            int ra_h, ra_m, dec_d, dec_arcm;
            double ra_s, dec_arcs;
            char dec_sign;
            ra_to_hms(ra_rad, ra_h, ra_m, ra_s);
            dec_to_dms(dec_rad, dec_sign, dec_d, dec_arcm, dec_arcs);
            
            // Stampa
            std::cout << std::setw(4) << hour << "h"
                      << "      "
                      << std::setw(2) << std::setfill('0') << ra_h << " "
                      << std::setw(2) << std::setfill('0') << ra_m << " "
                      << std::setw(6) << std::setfill(' ') << std::fixed << std::setprecision(3) << ra_s
                      << "      "
                      << dec_sign
                      << std::setw(2) << std::setfill('0') << dec_d << " "
                      << std::setw(2) << std::setfill('0') << dec_arcm << " "
                      << std::setw(5) << std::setfill(' ') << std::setprecision(2) << dec_arcs
                      << "      "
                      << std::setw(8) << std::setprecision(6) << delta << " AU"
                      << "  "
                      << std::setw(9) << std::setprecision(6) << ast_geo.x()
                      << "  "
                      << std::setw(9) << std::setprecision(6) << ast_geo.y()
                      << "  "
                      << std::setw(9) << std::setprecision(6) << ast_geo.z()
                      << "\n";
        }
        
        std::cout << "\n============================================================\n";
        std::cout << "CONFRONTO CON RIFERIMENTO (prime 4 ore):\n";
        std::cout << "------------------------------------------------------------\n";
        std::cout << "       MIE GEOCENTRICHE            RIFERIMENTO\n";
        std::cout << "Ora    RA           DEC             RA           DEC\n";
        std::cout << "------------------------------------------------------------\n";
        std::cout << "00h  (output)                    04 53 41.105  +20 19 56.76\n";
        std::cout << "01h  (output)                    04 53 39.005  +20 19 54.58\n";
        std::cout << "02h  (output)                    04 53 36.903  +20 19 52.40\n";
        std::cout << "03h  (output)                    04 53 34.802  +20 19 50.22\n";
        std::cout << "============================================================\n";
        std::cout << "\nNote:\n";
        std::cout << "  - Coordinate: GEOCENTRICHE (viste dalla Terra)\n";
        std::cout << "  - Sistema: ICRF J2000.0\n";
        std::cout << "  - Origine: Centro della Terra (Geocentro)\n";
        std::cout << "  - Delta: Distanza Terra-Asteroide\n";
        std::cout << "  - Tempo: UTC\n";
        std::cout << "  - Effemeridi Terra: VSOP87 (via AstDyn)\n";
        std::cout << "  - Precisione VSOP87: ~1\" per Terra su secoli\n";
        std::cout << "============================================================\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << "\n";
        return 1;
    }
}
