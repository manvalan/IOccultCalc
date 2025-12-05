/**
 * @file test_17030_occultation_ucac4_552_011427.cpp
 * @brief Verifica occultazione asteroide 17030 con stella UCAC4 552-011427
 * 
 * OCCULTAZIONE PREVISTA:
 * - Data: 28 Novembre 2025
 * - Orario: 00:28:18 - 00:42:39 UT (finestra 14 minuti)
 * - Momento centrale: ~00:35:28 UT
 * - Stella: UCAC4 552-011427 (Gaia DR3 3411546266140512128)
 * - RA stella: 04h 53m 39.8661s = 73.416084°
 * - Dec stella: +20° 19' 53.981" = 20.331661°
 * - Mag V: 12.13, Mag R: 11.05
 * - Durata max: 0.5 sec
 * - Larghezza ombra: 7.9 km
 * 
 * OBIETTIVO:
 * Propagare l'asteroide durante la finestra dell'occultazione e verificare
 * la minima distanza dalla stella target.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <astdyn/io/parsers/OrbFitEQ1Parser.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/core/Constants.hpp>

using namespace std;
using namespace astdyn;
using Eigen::Vector3d;

// Costanti
constexpr double MJD_TO_JD = 2400000.5;
constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;
constexpr double EPSILON_J2000 = 23.4392911 * DEG_TO_RAD;

// Stella target
constexpr double STAR_RA_DEG = 73.416084;   // 04h 53m 39.8661s
constexpr double STAR_DEC_DEG = 20.331661;  // +20° 19' 53.981"

// Conversione UTC a MJD TDB
double utc_to_mjd_tdb(int year, int month, int day, int hour, int min, double sec) {
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    double jd = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    jd += (hour - 12) / 24.0 + min / 1440.0 + sec / 86400.0;
    double mjd_utc = jd - 2400000.5;
    return mjd_utc + 69.184 / 86400.0; // Approssimazione UTC->TDB
}

// Conversione eclittica J2000 -> equatoriale ICRF J2000
Vector3d ecliptic_to_equatorial(const Vector3d& ecl) {
    double cos_eps = cos(EPSILON_J2000);
    double sin_eps = sin(EPSILON_J2000);
    return Vector3d(
        ecl[0],
        cos_eps * ecl[1] - sin_eps * ecl[2],
        sin_eps * ecl[1] + cos_eps * ecl[2]
    );
}

// Conversione cartesiano -> RA/Dec
void cartesian_to_radec(const Vector3d& pos, double& ra, double& dec) {
    double r = pos.norm();
    dec = asin(pos[2] / r);
    ra = atan2(pos[1], pos[0]);
    if (ra < 0) ra += 2 * M_PI;
}

// Formatta RA in ore
string format_ra(double ra_deg) {
    double ra_hours = ra_deg / 15.0;
    int h = static_cast<int>(ra_hours);
    double rem = (ra_hours - h) * 60.0;
    int m = static_cast<int>(rem);
    double s = (rem - m) * 60.0;
    char buf[32];
    snprintf(buf, sizeof(buf), "%02dh%02dm%05.2fs", h, m, s);
    return string(buf);
}

// Formatta Dec in gradi
string format_dec(double dec_deg) {
    char sign = dec_deg >= 0 ? '+' : '-';
    dec_deg = fabs(dec_deg);
    int d = static_cast<int>(dec_deg);
    double rem = (dec_deg - d) * 60.0;
    int m = static_cast<int>(rem);
    double s = (rem - m) * 60.0;
    char buf[32];
    snprintf(buf, sizeof(buf), "%c%02d°%02d'%05.2f\"", sign, d, m, s);
    return string(buf);
}

// Distanza angolare haversine
double angular_distance(double ra1, double dec1, double ra2, double dec2) {
    double delta_ra = ra2 - ra1;
    double delta_dec = dec2 - dec1;
    
    double a = sin(delta_dec / 2) * sin(delta_dec / 2) +
               cos(dec1) * cos(dec2) * sin(delta_ra / 2) * sin(delta_ra / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    
    return c; // radianti
}

int main() {
    cout << "\n========================================\n";
    cout << "VERIFICA OCCULTAZIONE 17030 × UCAC4 552-011427\n";
    cout << "========================================\n\n";
    
    cout << "Stella Target: UCAC4 552-011427\n";
    cout << "  Gaia DR3 Source ID: 3411546266140512128\n";
    cout << "  RA:  " << format_ra(STAR_RA_DEG) << " (" << fixed << setprecision(6) 
         << STAR_RA_DEG << "°)\n";
    cout << "  Dec: " << format_dec(STAR_DEC_DEG) << " (" << STAR_DEC_DEG << "°)\n";
    cout << "  Mag V: 12.13, Mag R: 11.05\n\n";
    
    cout << "Occultazione Prevista:\n";
    cout << "  Data: 28 Novembre 2025\n";
    cout << "  Finestra: 00:28:18 - 00:42:39 UT\n";
    cout << "  Durata max: 0.5 sec\n";
    cout << "  Larghezza ombra: 7.9 km\n\n";
    
    try {
        // 1. Carica elementi orbitali
        cout << "1. Caricamento elementi asteroide 17030...\n";
        string eq1_path = "../17030_astdys.eq1";
        io::parsers::OrbFitEQ1Parser parser;
        io::IOrbitParser::OrbitalElements elements = parser.parse(eq1_path);
        
        cout << "   ✓ a=" << elements.semi_major_axis 
             << " AU, e=" << elements.eccentricity << "\n\n";
        
        // 2. Converti a Keplerian
        propagation::KeplerianElements kep_elements;
        kep_elements.semi_major_axis = elements.semi_major_axis;
        kep_elements.eccentricity = elements.eccentricity;
        kep_elements.inclination = elements.inclination;
        kep_elements.longitude_ascending_node = elements.longitude_asc_node;
        kep_elements.argument_perihelion = elements.argument_perihelion;
        kep_elements.mean_anomaly = elements.mean_anomaly;
        kep_elements.epoch_mjd_tdb = elements.epoch_mjd_tdb;
        kep_elements.gravitational_parameter = 1.32712440018e20 / 1.495978707e11 / 
                                               1.495978707e11 / 1.495978707e11 * 
                                               86400.0 * 86400.0;
        
        // 3. Crea propagatore
        auto integrator = make_unique<propagation::RKF78Integrator>(0.1, 1e-12);
        auto ephemeris = make_shared<ephemeris::PlanetaryEphemeris>();
        
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
        
        cout << "2. Propagazione durante finestra occultazione...\n";
        cout << "   Campionamento: ogni 60 secondi (14 minuti = 15 punti)\n\n";
        
        // 4. Propaga durante la finestra (00:28:00 - 00:43:00, ogni 60 sec)
        double min_distance_arcsec = 1e10;
        double best_time_mjd = 0;
        Vector3d best_ast_pos_geo;
        
        cout << "   UTC Time     Asteroid RA/Dec              Distance from Star\n";
        cout << "   ----------------------------------------------------------------\n";
        
        for (int i = 0; i <= 15; ++i) {
            int minutes = 28 + i;
            int hour = minutes / 60;
            int min = minutes % 60;
            
            double mjd_tdb = utc_to_mjd_tdb(2025, 11, 28, hour, min, 0.0);
            
            // Propaga asteroide (barycentric ecliptic)
            auto kep_prop = propagator.propagate_keplerian(kep_elements, mjd_tdb);
            auto cart_ecl = propagation::keplerian_to_cartesian(kep_prop);
            
            // Converti a ICRF equatoriale
            Vector3d ast_bary_icrf = ecliptic_to_equatorial(cart_ecl.position);
            
            // Ottieni posizione Terra
            double jd_tdb = mjd_tdb + MJD_TO_JD;
            Vector3d earth_helio_ecl = ephemeris::PlanetaryEphemeris::getPosition(
                ephemeris::CelestialBody::EARTH, jd_tdb);
            Vector3d sun_bary_ecl = ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
            Vector3d earth_bary_ecl = earth_helio_ecl - sun_bary_ecl;
            Vector3d earth_bary_icrf = ecliptic_to_equatorial(earth_bary_ecl);
            
            // Posizione geocentrica
            Vector3d ast_geo_icrf = ast_bary_icrf - earth_bary_icrf;
            
            // Converti in RA/Dec
            double ra_rad, dec_rad;
            cartesian_to_radec(ast_geo_icrf, ra_rad, dec_rad);
            double ra_deg = ra_rad * RAD_TO_DEG;
            double dec_deg = dec_rad * RAD_TO_DEG;
            
            // Calcola distanza dalla stella
            double star_ra_rad = STAR_RA_DEG * DEG_TO_RAD;
            double star_dec_rad = STAR_DEC_DEG * DEG_TO_RAD;
            double dist_rad = angular_distance(ra_rad, dec_rad, star_ra_rad, star_dec_rad);
            double dist_arcsec = dist_rad * RAD_TO_DEG * 3600.0;
            
            // Aggiorna minimo
            if (dist_arcsec < min_distance_arcsec) {
                min_distance_arcsec = dist_arcsec;
                best_time_mjd = mjd_tdb;
                best_ast_pos_geo = ast_geo_icrf;
            }
            
            // Output
            char time_buf[16];
            snprintf(time_buf, sizeof(time_buf), "   00:%02d:00", min);
            
            cout << time_buf << "     "
                 << format_ra(ra_deg) << " " << format_dec(dec_deg)
                 << "     " << setw(8) << fixed << setprecision(2) << dist_arcsec << " arcsec";
            
            if (dist_arcsec < 10.0) {
                cout << "  ← CLOSEST!";
            }
            cout << "\n";
        }
        
        cout << "\n========================================\n";
        cout << "RISULTATI\n";
        cout << "========================================\n\n";
        
        cout << "Distanza minima dalla stella: " << fixed << setprecision(3) 
             << min_distance_arcsec << " arcsec\n";
        cout << "Momento migliore: MJD " << setprecision(6) << best_time_mjd << "\n\n";
        
        // Analisi risultato
        double shadow_width_km = 7.9;
        double distance_earth_km = best_ast_pos_geo.norm() * 1.495978707e8; // AU to km
        double angular_size_asteroid_arcsec = (shadow_width_km / distance_earth_km) * RAD_TO_DEG * 3600.0;
        
        cout << "Analisi occultazione:\n";
        cout << "  - Distanza Terra-Asteroide: " << setprecision(2) 
             << (distance_earth_km / 1e6) << " milioni km\n";
        cout << "  - Dimensione angolare ombra: " << setprecision(3) 
             << angular_size_asteroid_arcsec << " arcsec\n";
        cout << "  - Distanza minima stella: " << min_distance_arcsec << " arcsec\n\n";
        
        if (min_distance_arcsec < angular_size_asteroid_arcsec) {
            cout << "✓ OCCULTAZIONE CONFERMATA!\n";
            cout << "  L'asteroide passa davanti alla stella.\n";
            cout << "  Margine: " << (angular_size_asteroid_arcsec - min_distance_arcsec) 
                 << " arcsec\n";
        } else if (min_distance_arcsec < angular_size_asteroid_arcsec * 2) {
            cout << "⚠ OCCULTAZIONE POSSIBILE (near miss)\n";
            cout << "  Distanza vicina alle dimensioni dell'ombra.\n";
            cout << "  Differenza: " << (min_distance_arcsec - angular_size_asteroid_arcsec) 
                 << " arcsec\n";
        } else {
            cout << "✗ OCCULTAZIONE NON CONFERMATA\n";
            cout << "  Distanza troppo grande rispetto alle dimensioni dell'ombra.\n";
            cout << "  Differenza: " << (min_distance_arcsec - angular_size_asteroid_arcsec) 
                 << " arcsec\n";
        }
        
        cout << "\n";
        
    } catch (const exception& e) {
        cerr << "❌ ERRORE: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
