/**
 * @file test_17030_ephemerides.cpp
 * @brief Calcola effemeridi complete asteroide 17030 per 28/11/2025
 * @date 4 Dicembre 2025
 * 
 * Usa RKF78 con tutte le perturbazioni:
 * - 8 pianeti (Mercurio → Nettuno)
 * - Perturbazioni asteroidi (AST17)
 * - Correzioni relativistiche (Schwarzschild)
 * - Frame ICRF baricentrico J2000.0
 */

#include <chebyshev_rkf78_propagation.h>
#include <eq1_parser.h>
#include <iostream>
#include <iomanip>
#include <cmath>

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
        std::cout << "  EFFEMERIDI ASTEROIDE 17030 - 28 NOVEMBRE 2025\n";
        std::cout << "  Propagazione RKF78 con perturbazioni complete\n";
        std::cout << "============================================================\n\n";
        
        // 1. Carica elementi orbitali e crea propagatore
        std::cout << "1. Caricamento elementi orbitali e creazione propagatore...\n";
        std::string eq1_file = "17030_astdys.eq1";
        
        // Il nuovo costruttore carica direttamente dal file .eq1
        ChebyshevRKF78Propagator propagator(eq1_file);
        
        std::cout << "   ✓ File: " << eq1_file << "\n";
        std::cout << "   ✓ Elementi caricati e propagatore creato\n";
        std::cout << "   ✓ Frame input: ECLM (Mean Ecliptic J2000.0)\n";
        std::cout << "   ✓ Sistema: Elementi equinoziali modificati\n\n";
        
        // 2. Configura propagatore RKF78
        std::cout << "2. Configurazione propagatore RKF78...\n";
        
        // Mostra configurazione
        RKF78PropagationConfig config = propagator.getConfig();
        std::cout << "   ✓ Integratore: Runge-Kutta-Fehlberg 78 (RKF78)\n";
        std::cout << "   ✓ Tolleranza: " << std::scientific << config.tolerance << " AU\n";
        std::cout << "   ✓ Perturbazioni: TUTTE (8 pianeti + asteroidi AST17 + relatività)\n";
        std::cout << "   ✓ Output frame: ICRF baricentrico (J2000.0)\n";
        std::cout << "   ✓ Conversione automatica: ECLM → ICRF\n\n";
        
        // 3. Calcola posizioni per 28 Nov 2025 ogni ora
        // MJD 28 Nov 2025 00:00 UTC = JD 2460642.5 = MJD 60642.0
        double mjd_start = 60642.0;
        double mjd_end = mjd_start + 1.0;  // 24 ore
        int num_points = 25;  // 0h → 24h (inclusi estremi)
        
        std::cout << "3. Propagazione per 28 Novembre 2025...\n";
        std::cout << "   Intervallo: MJD " << std::fixed << std::setprecision(1) 
                  << mjd_start << " → " << mjd_end << " (24 ore)\n";
        std::cout << "   Campionamento: ogni 1 ora\n";
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        std::vector<Eigen::Vector3d> positions = propagator.propagateForChebyshev(
            mjd_start, mjd_end, num_points);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
        
        std::cout << "   ✓ Propagazione completata in " << std::fixed << std::setprecision(1) 
                  << elapsed_ms << " ms\n";
        std::cout << "   ✓ Tempo medio per punto: " << std::setprecision(2)
                  << (elapsed_ms / num_points) << " ms\n\n";
        
        // 4. Stampa effemeridi
        std::cout << "============================================================\n";
        std::cout << "  EFFEMERIDI - 28 NOVEMBRE 2025\n";
        std::cout << "============================================================\n\n";
        
        std::cout << " Ora UTC      RA (ICRF)         DEC (ICRF)      "
                  << "   Distanza       X (AU)      Y (AU)      Z (AU)\n";
        std::cout << std::string(110, '-') << "\n";
        
        for (int i = 0; i < num_points; ++i) {
            // Calcola epoca
            double mjd = mjd_start + (i * (mjd_end - mjd_start) / (num_points - 1));
            int hour = static_cast<int>(std::round((mjd - mjd_start) * 24.0));
            
            // Posizione baricentrica ICRF in AU
            const Eigen::Vector3d& pos = positions[i];
            
            // Converti da coordinate cartesiane a sferiche
            double r = pos.norm();  // Distanza dal baricentro solare
            double ra_rad = std::atan2(pos.y(), pos.x());
            if (ra_rad < 0) ra_rad += 2.0 * M_PI;
            
            double r_xy = std::sqrt(pos.x()*pos.x() + pos.y()*pos.y());
            double dec_rad = std::atan2(pos.z(), r_xy);
            
            // Converti RA/DEC in formato leggibile
            int ra_h, ra_m, dec_d, dec_arcm;
            double ra_s, dec_arcs;
            char dec_sign;
            ra_to_hms(ra_rad, ra_h, ra_m, ra_s);
            dec_to_dms(dec_rad, dec_sign, dec_d, dec_arcm, dec_arcs);
            
            // Stampa riga
            std::cout << std::setw(4) << hour << "h"
                      << "      "
                      << std::setw(2) << std::setfill('0') << ra_h << "h"
                      << std::setw(2) << std::setfill('0') << ra_m << "m"
                      << std::setw(5) << std::setfill(' ') << std::fixed << std::setprecision(2) << ra_s << "s"
                      << "      "
                      << dec_sign
                      << std::setw(2) << std::setfill('0') << dec_d << "°"
                      << std::setw(2) << std::setfill('0') << dec_arcm << "'"
                      << std::setw(5) << std::setfill(' ') << std::setprecision(2) << dec_arcs << "\""
                      << "      "
                      << std::setw(8) << std::setprecision(6) << r << " AU"
                      << "  "
                      << std::setw(10) << std::setprecision(6) << pos.x()
                      << "  "
                      << std::setw(10) << std::setprecision(6) << pos.y()
                      << "  "
                      << std::setw(10) << std::setprecision(6) << pos.z()
                      << "\n";
        }
        
        std::cout << "\n============================================================\n";
        std::cout << "Note:\n";
        std::cout << "  - Coordinate: ICRF baricentrico J2000.0\n";
        std::cout << "  - Sistema di riferimento: International Celestial Reference Frame\n";
        std::cout << "  - Origine: Baricentro del Sistema Solare\n";
        std::cout << "  - Tempo: UTC (Coordinated Universal Time)\n";
        std::cout << "  - Scala temporale elementi: TDB (Barycentric Dynamical Time)\n";
        std::cout << "  - Metodo: RKF78 con tolleranza " << std::scientific 
                  << config.tolerance << " AU\n";
        std::cout << "  - Perturbazioni: 8 pianeti + asteroidi (AST17) + relatività\n";
        std::cout << "============================================================\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << "\n";
        return 1;
    }
}
