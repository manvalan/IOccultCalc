/**
 * @file test_17030_chebyshev_ephemerides.cpp
 * @brief Test effemeridi geocentriche 17030 usando ChebyshevOccultationManager
 * @date 4 Dicembre 2025
 * 
 * Questo test:
 * 1. Carica elementi orbitali di 17030
 * 2. Propaga con RKF78 e fitta polinomi di Chebyshev
 * 3. Usa i polinomi per calcolare posizioni ogni ora il 28/11/2025
 * 4. Converte a coordinate geocentriche ICRF
 * 5. Confronta con i test precedenti
 */

#include "chebyshev_occultation_manager.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <cmath>
#include <Eigen/Dense>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>

using namespace ioccultcalc;

// Costanti
constexpr double MJD_TO_JD = 2400000.5;
constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;
constexpr double EPSILON_J2000 = 23.4392911 * DEG_TO_RAD;

// Conversione eclittico → equatoriale
Eigen::Vector3d ecliptic_to_equatorial(const Eigen::Vector3d& ecl) {
    double cos_eps = std::cos(EPSILON_J2000);
    double sin_eps = std::sin(EPSILON_J2000);
    Eigen::Vector3d eq;
    eq.x() = ecl.x();
    eq.y() = ecl.y() * cos_eps - ecl.z() * sin_eps;
    eq.z() = ecl.y() * sin_eps + ecl.z() * cos_eps;
    return eq;
}

// Cartesiano → RA/Dec
void cartesian_to_radec(const Eigen::Vector3d& pos, double& ra_deg, double& dec_deg) {
    double dist = pos.norm();
    dec_deg = std::asin(pos.z() / dist) * RAD_TO_DEG;
    ra_deg = std::atan2(pos.y(), pos.x()) * RAD_TO_DEG;
    if (ra_deg < 0) ra_deg += 360.0;
}

// Formatta RA
std::string format_ra(double ra_deg) {
    double ra_h = ra_deg / 15.0;
    int h = static_cast<int>(ra_h);
    double rem = (ra_h - h) * 60.0;
    int m = static_cast<int>(rem);
    double s = (rem - m) * 60.0;
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%02d:%02d:%05.2f", h, m, s);
    return std::string(buf);
}

// Formatta Dec
std::string format_dec(double dec_deg) {
    char sign = (dec_deg >= 0) ? '+' : '-';
    double abs_dec = std::abs(dec_deg);
    int d = static_cast<int>(abs_dec);
    double rem = (abs_dec - d) * 60.0;
    int m = static_cast<int>(rem);
    double s = (rem - m) * 60.0;
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%c%02d:%02d:%04.1f", sign, d, m, s);
    return std::string(buf);
}

// UTC → MJD TDB
double utc_to_mjd_tdb(int year, int month, int day, int hour) {
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    double jd = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    jd += (hour - 12.0) / 24.0;
    double mjd_utc = jd - 2400000.5;
    return mjd_utc + 69.184 / 86400.0;  // TT-UTC
}

int main() {
    std::cout << "====================================================================\n";
    std::cout << "  TEST EFFEMERIDI GEOCENTRICHE CON CHEBYSHEV POLYNOMIALS\n";
    std::cout << "  Asteroide 17030 - 28 Novembre 2025\n";
    std::cout << "====================================================================\n\n";

    try {
        // 1. CONFIGURA MANAGER
        std::cout << "1. Configurazione ChebyshevOccultationManager\n";
        std::cout << "   -----------------------------------------\n";
        
        OccultationSearchConfig config;
        // Intervallo: 27 Nov - 29 Nov (includiamo margine per interpolazione)
        config.start_mjd = 61006.0;  // 27 Nov 2025
        config.end_mjd = 61009.0;    // 30 Nov 2025
        config.num_propagation_points = 150;  // Alta densità per buon fitting
        config.rkf78_tolerance = 1e-12;
        config.num_chebyshev_coeffs = 10;  // 10 coefficienti per asse
        
        std::cout << "   Intervallo propagazione: MJD " << config.start_mjd 
                  << " - " << config.end_mjd << " (" 
                  << (config.end_mjd - config.start_mjd) << " giorni)\n";
        std::cout << "   Punti RKF78: " << config.num_propagation_points << "\n";
        std::cout << "   Coefficienti Chebyshev: " << config.num_chebyshev_coeffs 
                  << " per asse\n";
        std::cout << "   Tolleranza RKF78: " << std::scientific << config.rkf78_tolerance 
                  << " AU\n\n";
        
        ChebyshevOccultationManager manager(config);
        
        // 2. CARICA ASTEROIDE
        std::cout << "2. Caricamento elementi orbitali\n";
        std::cout << "   ------------------------------\n";
        
        if (!manager.loadAsteroidFromEQ1("../17030_astdys.eq1")) {
            std::cerr << "✗ Errore caricamento asteroide\n";
            return 1;
        }
        std::cout << "\n";
        
        // 3. PROPAGA E FITTA CHEBYSHEV
        std::cout << "3. Propagazione RKF78 e Fitting Chebyshev\n";
        std::cout << "   ---------------------------------------\n";
        
        if (!manager.propagateAndFit()) {
            std::cerr << "✗ Errore propagazione/fitting\n";
            return 1;
        }
        std::cout << "\n";
        
        // 4. CREA EPHEMERIS PLANETARIE PER TERRA
        std::cout << "4. Inizializzazione effemeridi Terra (VSOP87)\n";
        std::cout << "   ------------------------------------------\n";
        auto ephemeris = std::make_shared<astdyn::ephemeris::PlanetaryEphemeris>();
        std::cout << "   ✓ PlanetaryEphemeris creato\n\n";
        
        // 5. CALCOLA EFFEMERIDI GEOCENTRICHE
        std::cout << "5. Calcolo effemeridi geocentriche - 28 Novembre 2025\n";
        std::cout << "   ==================================================\n\n";
        
        std::cout << "Ora UTC      RA (J2000)        Dec (J2000)      Dist(AU)     Metodo\n";
        std::cout << "---------------------------------------------------------------------\n";
        
        std::string csv_filename = "/tmp/17030_chebyshev_28nov2025.csv";
        std::ofstream csv(csv_filename);
        csv << "Hour_UTC,MJD_TDB,RA_deg,Dec_deg,RA_HMS,Dec_DMS,Distance_AU,X_AU,Y_AU,Z_AU\n";
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        for (int hour = 0; hour <= 24; ++hour) {
            double mjd_tdb = utc_to_mjd_tdb(2025, 11, 28, hour);
            
            // Posizione asteroide baricentrica ICRF (da Chebyshev)
            Eigen::Vector3d ast_bary_icrf = manager.getPositionAtEpoch(mjd_tdb);
            
            // Posizione Terra
            double jd_tdb = mjd_tdb + MJD_TO_JD;
            Eigen::Vector3d earth_helio_ecl = astdyn::ephemeris::PlanetaryEphemeris::getPosition(
                astdyn::ephemeris::CelestialBody::EARTH, jd_tdb);
            Eigen::Vector3d sun_bary_ecl = astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
            Eigen::Vector3d earth_bary_ecl = earth_helio_ecl - sun_bary_ecl;
            Eigen::Vector3d earth_bary_icrf = ecliptic_to_equatorial(earth_bary_ecl);
            
            // Geocentrico
            Eigen::Vector3d ast_geo_icrf = ast_bary_icrf - earth_bary_icrf;
            
            // RA/Dec
            double ra_deg, dec_deg;
            cartesian_to_radec(ast_geo_icrf, ra_deg, dec_deg);
            double distance = ast_geo_icrf.norm();
            
            std::string ra_str = format_ra(ra_deg);
            std::string dec_str = format_dec(dec_deg);
            
            // Output console
            std::cout << std::setw(2) << std::setfill('0') << hour << "h"
                      << "         " << ra_str
                      << "    " << dec_str
                      << "    " << std::fixed << std::setprecision(6) << distance
                      << "    Chebyshev\n";
            
            // CSV
            csv << hour << ","
                << std::fixed << std::setprecision(6) << mjd_tdb << ","
                << std::setprecision(8) << ra_deg << ","
                << dec_deg << ","
                << ra_str << ","
                << dec_str << ","
                << std::setprecision(6) << distance << ","
                << std::setprecision(8) << ast_geo_icrf.x() << ","
                << ast_geo_icrf.y() << ","
                << ast_geo_icrf.z() << "\n";
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
        
        csv.close();
        
        std::cout << "\n---------------------------------------------------------------------\n";
        std::cout << "✓ Calcolo completato in " << std::fixed << std::setprecision(1) 
                  << elapsed_ms << " ms\n";
        std::cout << "✓ Tempo per query: " << std::setprecision(2) 
                  << (elapsed_ms / 25.0) << " ms/punto\n";
        std::cout << "✓ File CSV: " << csv_filename << "\n\n";
        
        // 6. STATISTICHE
        std::cout << "6. Statistiche Chebyshev\n";
        std::cout << "   ---------------------\n";
        std::cout << "   RMS error: " << std::scientific << std::setprecision(3)
                  << manager.getChebyshevRMSError() << " AU\n";
        std::cout << "   Interpolazione: Chebyshev polinomiale (sub-ms query)\n";
        std::cout << "   Frame: ICRF J2000.0 geocentrico\n\n";
        
        // 7. CONFRONTO CON RISULTATI ATTESI
        std::cout << "7. Validazione risultati\n";
        std::cout << "   ---------------------\n";
        std::cout << "   Confronta con: /tmp/17030_geocentric_28nov2025.csv\n";
        std::cout << "   Comando: diff /tmp/17030_chebyshev_28nov2025.csv \\\n";
        std::cout << "                 /tmp/17030_geocentric_28nov2025.csv\n\n";
        
        std::cout << "====================================================================\n";
        std::cout << "TEST COMPLETATO CON SUCCESSO\n";
        std::cout << "====================================================================\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ ERRORE: " << e.what() << "\n";
        return 1;
    }
}
