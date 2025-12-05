/**
 * @file test_17030_positions.cpp
 * @brief Stampa posizioni asteroide 17030 il 28/11/2025 ogni 2 ore
 * @date 4 Dicembre 2025
 */

#include "chebyshev_occultation_manager.h"
#include <iostream>
#include <iomanip>

using namespace ioccultcalc;

int main() {
    try {
        std::cout << "========================================\n";
        std::cout << "  Posizioni Asteroide 17030\n";
        std::cout << "  Data: 28 Novembre 2025\n";
        std::cout << "  Campionamento: ogni 2 ore\n";
        std::cout << "========================================\n\n";
        
        // MJD per 28 Nov 2025 00:00 UTC
        // 2025-11-28 = JD 2460642.5 = MJD 60642.0
        double mjd_target = 60642.0;  // 28 Nov 2025 00:00 UTC
        
        // Configurazione - propaga da inizio a fine della giornata
        OccultationSearchConfig config;
        config.start_mjd = mjd_target;         // 28 Nov 00:00
        config.end_mjd = mjd_target + 1.0;     // 29 Nov 00:00
        config.num_propagation_points = 200;
        config.num_chebyshev_coeffs = 12;      // Più accurato
        
        std::cout << "1. Caricamento asteroide 17030...\n";
        ChebyshevOccultationManager manager(config);
        
        if (!manager.loadAsteroidFromEQ1("17030_astdys.eq1")) {
            std::cerr << "ERRORE: Impossibile caricare 17030_astdys.eq1\n";
            return 1;
        }
        std::cout << "   ✓ Asteroide caricato\n\n";
        
        std::cout << "2. Propagazione e fitting Chebyshev...\n";
        if (!manager.propagateAndFit()) {
            std::cerr << "ERRORE: Propagazione fallita\n";
            return 1;
        }
        std::cout << "   ✓ Propagazione completata\n\n";
        
        std::cout << "========================================\n";
        std::cout << "  POSIZIONI 28 NOVEMBRE 2025\n";
        std::cout << "========================================\n\n";
        
        std::cout << std::setw(8) << "Ora UTC"
                  << std::setw(12) << "RA (°)"
                  << std::setw(12) << "DEC (°)"
                  << std::setw(18) << "RA (h m s)"
                  << std::setw(18) << "DEC (° ' \")"
                  << "\n";
        std::cout << std::string(68, '-') << "\n";
        
        // Stampa posizioni ogni 2 ore
        for (int hour = 0; hour <= 24; hour += 2) {
            double mjd = mjd_target + hour / 24.0;
            
            std::cout << "DEBUG: Richiedendo posizione a MJD " << mjd << std::endl;
            
            // Ottieni posizione dall'interpolazione Chebyshev
            Eigen::Vector3d pos_icrf;
            try {
                pos_icrf = manager.getPositionAtEpoch(mjd);
                std::cout << "DEBUG: Posizione ottenuta: " << pos_icrf.transpose() << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "ERRORE getPositionAtEpoch: " << e.what() << std::endl;
                return 1;
            }
            
            // Converti da coordinate cartesiane ICRF a RA/Dec
            double ra_rad = std::atan2(pos_icrf.y(), pos_icrf.x());
            if (ra_rad < 0) ra_rad += 2.0 * M_PI;
            
            double r_xy = std::sqrt(pos_icrf.x()*pos_icrf.x() + pos_icrf.y()*pos_icrf.y());
            double dec_rad = std::atan2(pos_icrf.z(), r_xy);
            
            // Converti in gradi
            double ra_deg = ra_rad * 180.0 / M_PI;
            double dec_deg = dec_rad * 180.0 / M_PI;
            
            // Converti RA in ore:minuti:secondi
            double ra_hours = ra_deg / 15.0;
            int ra_h = static_cast<int>(ra_hours);
            double ra_m_frac = (ra_hours - ra_h) * 60.0;
            int ra_m = static_cast<int>(ra_m_frac);
            double ra_s = (ra_m_frac - ra_m) * 60.0;
            
            // Converti DEC in gradi:minuti:secondi
            bool dec_neg = dec_deg < 0;
            double dec_abs = std::abs(dec_deg);
            int dec_d = static_cast<int>(dec_abs);
            double dec_m_frac = (dec_abs - dec_d) * 60.0;
            int dec_arcm = static_cast<int>(dec_m_frac);
            double dec_arcs = (dec_m_frac - dec_arcm) * 60.0;
            
            // Stampa
            std::cout << std::setw(6) << hour << "h"
                      << std::fixed << std::setprecision(6)
                      << std::setw(12) << ra_deg
                      << std::setw(12) << dec_deg
                      << "    "
                      << std::setw(2) << std::setfill('0') << ra_h << "h"
                      << std::setw(2) << std::setfill('0') << ra_m << "m"
                      << std::setw(5) << std::setfill(' ') << std::setprecision(2) << ra_s << "s"
                      << "    "
                      << (dec_neg ? "-" : "+")
                      << std::setw(2) << std::setfill('0') << dec_d << "°"
                      << std::setw(2) << std::setfill('0') << dec_arcm << "'"
                      << std::setw(5) << std::setfill(' ') << std::setprecision(2) << dec_arcs << "\""
                      << "\n";
        }
        
        std::cout << "\n========================================\n";
        std::cout << "Note:\n";
        std::cout << "  - Coordinate: ICRS (barycentric)\n";
        std::cout << "  - Epoca: J2000.0\n";
        std::cout << "  - Tempo: UTC\n";
        std::cout << "========================================\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << "\n";
        return 1;
    }
}
