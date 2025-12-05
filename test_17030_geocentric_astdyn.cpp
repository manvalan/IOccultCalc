/**
 * @file test_17030_geocentric_astdyn.cpp
 * @brief Coordinate equatoriali geocentriche per asteroide 17030
 * @date 4 Dicembre 2025
 * 
 * Data: 28/11/2025 dalle 00:00 alle 24:00 UTC (ogni ora)
 * Usa: AstDyn Propagator con RKF78, tutte le perturbazioni
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/io/parsers/OrbFitEQ1Parser.hpp>
#include <astdyn/core/Constants.hpp>

using namespace astdyn;
using namespace Eigen;

// Costanti
constexpr double MJD_TO_JD = 2400000.5;
constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Obliquità eclittica J2000
constexpr double EPSILON_J2000 = 23.4392911 * DEG_TO_RAD;

// Conversione UTC a MJD TDB (approssimata, ignora leap seconds ~70s)
double utc_to_mjd_tdb(int year, int month, int day, int hour, int min, double sec) {
    // Formula di Meeus per JD
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    double jd = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    jd += (hour - 12) / 24.0 + min / 1440.0 + sec / 86400.0;
    double mjd = jd - MJD_TO_JD;
    // TDB ≈ TT ≈ UTC + 69.184s (al 2025)
    return mjd + 69.184 / 86400.0;
}

// Rotazione ECLM J2000 → ICRF (equatoriale)
Vector3d ecliptic_to_equatorial(const Vector3d& ecl) {
    double cos_eps = std::cos(EPSILON_J2000);
    double sin_eps = std::sin(EPSILON_J2000);
    return Vector3d(
        ecl.x(),
        ecl.y() * cos_eps - ecl.z() * sin_eps,
        ecl.y() * sin_eps + ecl.z() * cos_eps
    );
}

// Conversione cartesiano equatoriale → RA/Dec
void cartesian_to_radec(const Vector3d& pos, double& ra_deg, double& dec_deg) {
    double dist_au = pos.norm();
    dec_deg = std::asin(pos.z() / dist_au) * RAD_TO_DEG;
    ra_deg = std::atan2(pos.y(), pos.x()) * RAD_TO_DEG;
    if (ra_deg < 0) ra_deg += 360.0;
}

// Formatta RA in HH:MM:SS.ss
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

// Formatta Dec in ±DD:MM:SS.s
std::string format_dec(double dec_deg) {
    char sign = dec_deg >= 0 ? '+' : '-';
    dec_deg = std::fabs(dec_deg);
    int d = static_cast<int>(dec_deg);
    double rem = (dec_deg - d) * 60.0;
    int m = static_cast<int>(rem);
    double s = (rem - m) * 60.0;
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%c%02d:%02d:%04.1f", sign, d, m, s);
    return std::string(buf);
}

int main() {
    std::cout << "=== Coordinate Equatoriali Geocentriche Asteroide 17030 ===\n";
    std::cout << "Data: 28 Novembre 2025, 00:00-24:00 UTC (ogni ora)\n";
    std::cout << "Integratore: RKF78, tolleranza 1e-12 AU\n";
    std::cout << "Perturbazioni: 8 pianeti + asteroidi + relatività\n";
    std::cout << "Frame output: ICRF J2000 geocentrico\n\n";

    // 1. Carica elementi orbitali da file .eq1
    std::string eq1_path = "../17030_astdys.eq1";
    // 2. Carica elementi orbitali
    io::parsers::OrbFitEQ1Parser parser;
    io::IOrbitParser::OrbitalElements elements;
    
    try {
        elements = parser.parse(eq1_path);
        std::cout << "✓ Elementi caricati: a=" << elements.semi_major_axis 
                  << " AU, e=" << elements.eccentricity << "\n";
    } catch (const std::exception& e) {
        std::cerr << "✗ Errore caricamento .eq1: " << e.what() << "\n";
        return 1;
    }

    // 3. Converti a KeplerianElements
    propagation::KeplerianElements kep_elements;
    kep_elements.semi_major_axis = elements.semi_major_axis;
    kep_elements.eccentricity = elements.eccentricity;
    kep_elements.inclination = elements.inclination;
    kep_elements.longitude_ascending_node = elements.longitude_asc_node;
    kep_elements.argument_perihelion = elements.argument_perihelion;
    kep_elements.mean_anomaly = elements.mean_anomaly;
    kep_elements.epoch_mjd_tdb = elements.epoch_mjd_tdb;
    kep_elements.gravitational_parameter = 1.32712440018e20 / 1.495978707e11 / 1.495978707e11 / 1.495978707e11 * 86400.0 * 86400.0;
    
    // 4. Crea integrator e ephemeris
    auto integrator = std::make_unique<propagation::RKF78Integrator>(0.1, 1e-12);
    auto ephemeris = std::make_shared<ephemeris::PlanetaryEphemeris>();
    
    // 5. Configura propagatore con tutte le perturbazioni
    propagation::PropagatorSettings settings;
    settings.include_planets = true;
    settings.include_asteroids = true;
    settings.include_relativity = true;
    settings.perturb_mercury = true;
    settings.perturb_venus = true;
    settings.perturb_earth = true;
    settings.perturb_mars = true;
    settings.perturb_jupiter = true;
    settings.perturb_saturn = true;
    settings.perturb_uranus = true;
    settings.perturb_neptune = true;
    
    propagation::Propagator propagator(std::move(integrator), ephemeris, settings);
    std::cout << "✓ Propagatore creato con tutte le perturbazioni\n\n";
    
    // 6. Output file CSV
    std::string csv_filename = "/tmp/17030_geocentric_28nov2025.csv";
    std::ofstream csv_file(csv_filename);
    if (!csv_file) {
        std::cerr << "✗ Errore creazione file CSV\n";
        return 1;
    }
    
    csv_file << "Hour_UTC,MJD_TDB,RA_deg,Dec_deg,RA_HMS,Dec_DMS,Distance_AU,X_AU,Y_AU,Z_AU\n";
    
    // 7. Header output console
    std::cout << "EFFEMERIDI GEOCENTRICHE - 28 NOVEMBRE 2025\n";
    std::cout << "==========================================\n\n";
    std::cout << "Ora UTC      RA (J2000)        Dec (J2000)      Dist(AU)\n";
    std::cout << "---------------------------------------------------------------\n";
    
    // 8. Loop per ogni ora
    auto start_time = std::chrono::high_resolution_clock::now();
    
    for (int hour = 0; hour <= 24; ++hour) {
        // Converti UTC → MJD TDB
        double mjd_tdb = utc_to_mjd_tdb(2025, 11, 28, hour, 0, 0.0);
        
        // Propaga asteroide (baricentrico eclittico J2000)
        auto kep_prop = propagator.propagate_keplerian(kep_elements, mjd_tdb);
        auto cart_ecl = propagation::keplerian_to_cartesian(kep_prop);
        
        // Converti asteroide a ICRF equatoriale
        Vector3d ast_bary_icrf = ecliptic_to_equatorial(cart_ecl.position);
        
        // Posizione Terra (eliocentrico eclittico J2000)
        double jd_tdb = mjd_tdb + MJD_TO_JD;
        Vector3d earth_helio_ecl = ephemeris::PlanetaryEphemeris::getPosition(
            ephemeris::CelestialBody::EARTH, jd_tdb);
        
        // Posizione Sole (baricentrico eclittico J2000)
        Vector3d sun_bary_ecl = ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
        
        // Terra baricentrica = Earth_helio - Sun_bary (eclittico)
        Vector3d earth_bary_ecl = earth_helio_ecl - sun_bary_ecl;
        
        // Converti Terra a ICRF equatoriale
        Vector3d earth_bary_icrf = ecliptic_to_equatorial(earth_bary_ecl);
        
        // Vettore geocentrico
        Vector3d ast_geo_icrf = ast_bary_icrf - earth_bary_icrf;
        
        // RA, Dec, distanza
        double ra_deg, dec_deg;
        cartesian_to_radec(ast_geo_icrf, ra_deg, dec_deg);
        double distance = ast_geo_icrf.norm();
        
        // Formatta output
        std::string ra_str = format_ra(ra_deg);
        std::string dec_str = format_dec(dec_deg);
        
        // Console output
        std::cout << std::setw(2) << std::setfill('0') << hour << "h"
                  << "         " << ra_str
                  << "    " << dec_str
                  << "    " << std::fixed << std::setprecision(6) << distance << "\n";
        
        // CSV output
        csv_file << hour << ","
                 << std::fixed << std::setprecision(6) << mjd_tdb << ","
                 << std::setprecision(8) << ra_deg << ","
                 << std::setprecision(8) << dec_deg << ","
                 << ra_str << ","
                 << dec_str << ","
                 << std::setprecision(6) << distance << ","
                 << std::setprecision(8) << ast_geo_icrf.x() << ","
                 << ast_geo_icrf.y() << ","
                 << ast_geo_icrf.z() << "\n";
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    csv_file.close();
    
    std::cout << "\n---------------------------------------------------------------\n";
    std::cout << "✓ Calcolo completato in " << std::fixed << std::setprecision(1) 
              << elapsed_ms << " ms\n";
    std::cout << "✓ File CSV salvato: " << csv_filename << "\n";
    std::cout << "\nNote:\n";
    std::cout << "  - Frame: ICRF J2000.0 geocentrico\n";
    std::cout << "  - Propagatore: AstDyn RKF78 (tolleranza 1e-12 AU)\n";
    std::cout << "  - Effemeridi Terra: VSOP87\n";
    std::cout << "  - Perturbazioni: 8 pianeti + relatività\n";
    
    return 0;
}
