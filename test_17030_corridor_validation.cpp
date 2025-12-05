/**
 * @file test_17030_corridor_validation.cpp
 * @brief Validazione completa dell'API Corridor per l'asteroide 17030
 * 
 * OBIETTIVO:
 * Verificare che TUTTE le stelle restituite dalla query Corridor siano effettivamente
 * alla distanza uguale o inferiore alla fascia (corridor width) specificata.
 * 
 * TEST:
 * 1. Carica elementi orbitali di 17030
 * 2. Propaga posizioni per creare il path
 * 3. Esegue query Corridor con larghezza specifica (es. 0.02 gradi = 72 arcsec)
 * 4. Per OGNI stella restituita:
 *    - Calcola distanza minima dal path usando haversine
 *    - Verifica che distanza <= corridor_width
 *    - Segnala eventuali violazioni
 * 5. Stampa statistiche complete
 * 
 * Data: 4 dicembre 2025
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <fstream>

// AstDyn headers
#include <astdyn/io/parsers/OrbFitEQ1Parser.hpp>
#include <astdyn/propagation/Propagator.hpp>
#include <astdyn/propagation/Integrator.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>
#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/core/Constants.hpp>

// IOC_GaiaLib headers
#include "ioc_gaialib/unified_gaia_catalog.h"
#include "ioc_gaialib/types.h"

using namespace std;
using namespace astdyn;
using namespace ioc::gaia;
using Eigen::Vector3d;

// === UTILITY FUNCTIONS ===

/**
 * Converte data UTC in MJD TDB
 */
double utc_to_mjd_tdb(int year, int month, int day, int hour, int minute = 0, int second = 0) {
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    double jd = jdn + (hour - 12) / 24.0 + minute / 1440.0 + second / 86400.0;
    double mjd_utc = jd - 2400000.5;
    double mjd_tdb = mjd_utc + 69.184 / 86400.0;
    return mjd_tdb;
}

/**
 * Converte coordinate eclittiche J2000 in ICRF equatoriali J2000
 * Rotazione intorno all'asse X di epsilon (obliquità dell'eclittica)
 */
void ecliptic_to_icrf(double x_ecl, double y_ecl, double z_ecl,
                       double& x_icrf, double& y_icrf, double& z_icrf) {
    const double epsilon = 23.4392911 * M_PI / 180.0; // obliquità J2000
    const double cos_eps = cos(epsilon);
    const double sin_eps = sin(epsilon);
    
    x_icrf = x_ecl;
    y_icrf = cos_eps * y_ecl - sin_eps * z_ecl;
    z_icrf = sin_eps * y_ecl + cos_eps * z_ecl;
}

/**
 * Converte coordinate cartesiane in RA/Dec
 */
void cartesian_to_radec(double x, double y, double z, double& ra, double& dec) {
    double r = sqrt(x*x + y*y + z*z);
    dec = asin(z / r);
    ra = atan2(y, x);
    if (ra < 0) ra += 2 * M_PI;
}

/**
 * Calcola la distanza angolare (great circle) tra due punti sulla sfera celeste
 * usando la formula haversine
 * 
 * @param ra1, dec1 Coordinate del primo punto [radianti]
 * @param ra2, dec2 Coordinate del secondo punto [radianti]
 * @return Distanza angolare [gradi]
 */
double haversine_distance(double ra1, double dec1, double ra2, double dec2) {
    double delta_ra = ra2 - ra1;
    double delta_dec = dec2 - dec1;
    
    double a = sin(delta_dec / 2) * sin(delta_dec / 2) +
               cos(dec1) * cos(dec2) * sin(delta_ra / 2) * sin(delta_ra / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    
    return c * 180.0 / M_PI; // Converti in gradi
}

/**
 * Calcola la distanza minima di un punto da un segmento sulla sfera celeste
 * 
 * @param ra_point, dec_point Coordinate del punto [radianti]
 * @param ra1, dec1 Coordinate inizio segmento [radianti]
 * @param ra2, dec2 Coordinate fine segmento [radianti]
 * @return Distanza minima [gradi]
 */
double distance_point_to_segment(double ra_point, double dec_point,
                                  double ra1, double dec1,
                                  double ra2, double dec2) {
    // Distanze dai due estremi
    double dist1 = haversine_distance(ra_point, dec_point, ra1, dec1);
    double dist2 = haversine_distance(ra_point, dec_point, ra2, dec2);
    
    // Lunghezza del segmento
    double seg_length = haversine_distance(ra1, dec1, ra2, dec2);
    
    // Se il segmento è molto corto, usa la distanza minima dai due estremi
    if (seg_length < 1e-10) {
        return min(dist1, dist2);
    }
    
    // Calcola la proiezione del punto sul segmento usando dot product
    // Converti in coordinate cartesiane
    double x1 = cos(dec1) * cos(ra1);
    double y1 = cos(dec1) * sin(ra1);
    double z1 = sin(dec1);
    
    double x2 = cos(dec2) * cos(ra2);
    double y2 = cos(dec2) * sin(ra2);
    double z2 = sin(dec2);
    
    double xp = cos(dec_point) * cos(ra_point);
    double yp = cos(dec_point) * sin(ra_point);
    double zp = sin(dec_point);
    
    // Vettore dal punto 1 al punto 2
    double vx = x2 - x1;
    double vy = y2 - y1;
    double vz = z2 - z1;
    
    // Vettore dal punto 1 al punto P
    double wx = xp - x1;
    double wy = yp - y1;
    double wz = zp - z1;
    
    // Parametro t della proiezione (0 <= t <= 1 se la proiezione è sul segmento)
    double dot_vw = vx*wx + vy*wy + vz*wz;
    double dot_vv = vx*vx + vy*vy + vz*vz;
    double t = dot_vw / dot_vv;
    
    // Se la proiezione è prima del punto 1, usa dist1
    if (t <= 0.0) {
        return dist1;
    }
    // Se la proiezione è dopo il punto 2, usa dist2
    if (t >= 1.0) {
        return dist2;
    }
    
    // La proiezione è sul segmento, calcola il punto proiettato
    double x_proj = x1 + t * vx;
    double y_proj = y1 + t * vy;
    double z_proj = z1 + t * vz;
    
    // Normalizza (per riportarlo sulla sfera unitaria)
    double norm = sqrt(x_proj*x_proj + y_proj*y_proj + z_proj*z_proj);
    x_proj /= norm;
    y_proj /= norm;
    z_proj /= norm;
    
    // Converti back in ra/dec
    double ra_proj = atan2(y_proj, x_proj);
    double dec_proj = asin(z_proj);
    
    // Calcola distanza dal punto al punto proiettato
    return haversine_distance(ra_point, dec_point, ra_proj, dec_proj);
}

/**
 * Calcola la distanza minima di un punto da un path composto da più segmenti
 * 
 * @param ra_point, dec_point Coordinate del punto [gradi]
 * @param path Vettore di punti del path (CelestialPoint in gradi)
 * @return Distanza minima [gradi]
 */
double distance_point_to_path(double ra_point, double dec_point,
                               const vector<CelestialPoint>& path) {
    double min_distance = 1e10;
    
    // Converti punto in radianti
    double ra_p_rad = ra_point * M_PI / 180.0;
    double dec_p_rad = dec_point * M_PI / 180.0;
    
    // Itera su tutti i segmenti del path
    for (size_t i = 0; i < path.size() - 1; ++i) {
        double ra1_rad = path[i].ra * M_PI / 180.0;
        double dec1_rad = path[i].dec * M_PI / 180.0;
        double ra2_rad = path[i+1].ra * M_PI / 180.0;
        double dec2_rad = path[i+1].dec * M_PI / 180.0;
        
        double dist = distance_point_to_segment(ra_p_rad, dec_p_rad,
                                                 ra1_rad, dec1_rad,
                                                 ra2_rad, dec2_rad);
        min_distance = min(min_distance, dist);
    }
    
    return min_distance;
}

// === MAIN TEST ===

int main(int argc, char* argv[]) {
    cout << "\n========================================\n";
    cout << "TEST VALIDAZIONE CORRIDOR API\n";
    cout << "Asteroide: (17030) 1999 CQ3\n";
    cout << "Data: 28 novembre 2025\n";
    cout << "========================================\n\n";
    
    try {
        // === 1. CARICA ELEMENTI ORBITALI ===
        cout << "1. Caricamento elementi orbitali...\n";
        
        string eq1_path = "../17030_astdys.eq1";
        io::parsers::OrbFitEQ1Parser parser;
        io::IOrbitParser::OrbitalElements orbital_elements;
        
        try {
            orbital_elements = parser.parse(eq1_path);
            cout << "   ✓ Elementi caricati:\n";
            cout << "     Semi-asse maggiore: " << orbital_elements.semi_major_axis << " AU\n";
            cout << "     Eccentricità: " << orbital_elements.eccentricity << "\n";
            cout << "     Epoca: MJD " << orbital_elements.epoch_mjd_tdb << "\n\n";
        } catch (const exception& e) {
            cerr << "   ✗ Errore caricamento: " << e.what() << "\n";
            return 1;
        }
        
        // === 2. CONVERTI IN KEPLERIAN ELEMENTS ===
        propagation::KeplerianElements kep_elements;
        kep_elements.semi_major_axis = orbital_elements.semi_major_axis;
        kep_elements.eccentricity = orbital_elements.eccentricity;
        kep_elements.inclination = orbital_elements.inclination;
        kep_elements.longitude_ascending_node = orbital_elements.longitude_asc_node;
        kep_elements.argument_perihelion = orbital_elements.argument_perihelion;
        kep_elements.mean_anomaly = orbital_elements.mean_anomaly;
        kep_elements.epoch_mjd_tdb = orbital_elements.epoch_mjd_tdb;
        kep_elements.gravitational_parameter = 1.32712440018e20 / 1.495978707e11 / 1.495978707e11 / 1.495978707e11 * 86400.0 * 86400.0;
        
        // === 3. PROPAGA POSIZIONI PER 28 NOV 2025 (24 ore) ===
        cout << "2. Propagazione posizioni per il 28 novembre 2025...\n";
        
        // Crea integrator e ephemeris
        auto integrator = make_unique<propagation::RKF78Integrator>(0.1, 1e-12);
        auto ephemeris = make_shared<ephemeris::PlanetaryEphemeris>();
        
        // Configurazione con tutte le perturbazioni
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
        
        // Crea propagatore
        propagation::Propagator propagator(std::move(integrator), ephemeris, settings);
        
        // Genera 25 epoche (0h - 24h, ogni ora)
        vector<double> mjds;
        vector<CelestialPoint> path;
        
        double start_mjd = utc_to_mjd_tdb(2025, 11, 28, 0);
        int num_points = 25;
        
        cout << "   Propagando " << num_points << " punti (ogni ora)...\n";
        auto t_start = chrono::high_resolution_clock::now();
        
        for (int i = 0; i < num_points; ++i) {
            double mjd = start_mjd + i / 24.0;
            mjds.push_back(mjd);
            
            // Propaga (barycentric ecliptic J2000)
            auto kep_prop = propagator.propagate_keplerian(kep_elements, mjd);
            auto cart_ecl = propagation::keplerian_to_cartesian(kep_prop);
            
            // Converti eclittica → ICRF
            double x_icrf, y_icrf, z_icrf;
            ecliptic_to_icrf(cart_ecl.position[0], cart_ecl.position[1], cart_ecl.position[2],
                            x_icrf, y_icrf, z_icrf);
            
            // Converti in RA/Dec
            double ra_rad, dec_rad;
            cartesian_to_radec(x_icrf, y_icrf, z_icrf, ra_rad, dec_rad);
            
            double ra_deg = ra_rad * 180.0 / M_PI;
            double dec_deg = dec_rad * 180.0 / M_PI;
            
            path.push_back(CelestialPoint(ra_deg, dec_deg));
        }
        
        auto t_end = chrono::high_resolution_clock::now();
        auto duration_ms = chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();
        
        cout << "   ✓ Propagazione completata in " << duration_ms << " ms\n";
        cout << "   ✓ Path creato con " << path.size() << " punti\n";
        cout << "     Prima posizione: RA " << fixed << setprecision(6) << path[0].ra 
             << "°, Dec " << path[0].dec << "°\n";
        cout << "     Ultima posizione: RA " << path[num_points-1].ra 
             << "°, Dec " << path[num_points-1].dec << "°\n\n";
        
        // === 4. INIZIALIZZA UNIFIED GAIA CATALOG ===
        cout << "3. Inizializzazione UnifiedGaiaCatalog...\n";
        
        string home = string(getenv("HOME"));
        string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
        
        string json_config = R"({
            "catalog_type": "multifile_v2",
            "multifile_directory": ")" + catalog_path + R"("
        })";
        
        if (!UnifiedGaiaCatalog::initialize(json_config)) {
            throw runtime_error("Impossibile inizializzare UnifiedGaiaCatalog. Verificare path: " + catalog_path);
        }
        
        auto& catalog = UnifiedGaiaCatalog::getInstance();
        auto info = catalog.getCatalogInfo();
        cout << "   ✓ Catalogo inizializzato: " << catalog_path << "\n";
        cout << "   ✓ Stelle totali: " << info.total_stars << "\n\n";
        
        // === 5. ESEGUI QUERY CORRIDOR ===
        // Per la fase 1 (screening candidate) usiamo un corridor più largo
        // per trovare tutte le stelle che potrebbero avere un closest approach interessante
        double corridor_width_deg = 0.05; // 180 arcsec (3 arcmin) - screening iniziale
        double max_magnitude = 18.0;
        
        cout << "4. Esecuzione query Corridor...\n";
        cout << "   Parametri:\n";
        cout << "     - Larghezza corridor: " << corridor_width_deg << " gradi ("
             << (corridor_width_deg * 3600) << " arcsec)\n";
        cout << "     - Magnitudine limite: " << max_magnitude << "\n";
        cout << "     - Punti del path: " << path.size() << "\n\n";
        
        CorridorQueryParams params;
        params.path = path;
        params.width = corridor_width_deg;
        params.max_magnitude = max_magnitude;
        params.min_parallax = -1.0; // nessun limite
        params.max_results = 0; // nessun limite
        
        t_start = chrono::high_resolution_clock::now();
        auto stars = catalog.queryCorridor(params);
        t_end = chrono::high_resolution_clock::now();
        duration_ms = chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();
        
        cout << "   ✓ Query completata in " << duration_ms << " ms\n";
        cout << "   ✓ Stelle trovate: " << stars.size() << "\n\n";
        
        // === 6. VALIDAZIONE: VERIFICA DISTANZE ===
        cout << "5. VALIDAZIONE DISTANZE DA PATH\n";
        cout << "   Verifico che TUTTE le " << stars.size() << " stelle siano entro " 
             << corridor_width_deg << " gradi...\n\n";
        
        size_t num_valid = 0;
        size_t num_invalid = 0;
        double max_distance_found = 0.0;
        double min_distance_found = 1e10;
        double sum_distances = 0.0;
        
        vector<pair<GaiaStar, double>> violations; // stelle fuori dal corridor con la loro distanza
        
        for (const auto& star : stars) {
            // Calcola distanza minima dal path
            double min_dist = distance_point_to_path(star.ra, star.dec, path);
            
            sum_distances += min_dist;
            max_distance_found = max(max_distance_found, min_dist);
            min_distance_found = min(min_distance_found, min_dist);
            
            // Verifica se è entro il corridor
            if (min_dist <= corridor_width_deg) {
                num_valid++;
            } else {
                num_invalid++;
                violations.push_back({star, min_dist});
            }
        }
        
        double avg_distance = sum_distances / stars.size();
        
        // === 7. STAMPA RISULTATI ===
        cout << "========================================\n";
        cout << "RISULTATI VALIDAZIONE\n";
        cout << "========================================\n\n";
        
        cout << "Stelle analizzate: " << stars.size() << "\n";
        cout << "Stelle VALIDE (dentro corridor): " << num_valid << "\n";
        cout << "Stelle INVALIDE (fuori corridor): " << num_invalid << "\n\n";
        
        cout << "Statistiche distanze:\n";
        cout << "  - Minima: " << fixed << setprecision(6) << min_distance_found 
             << " gradi (" << (min_distance_found * 3600) << " arcsec)\n";
        cout << "  - Massima: " << max_distance_found 
             << " gradi (" << (max_distance_found * 3600) << " arcsec)\n";
        cout << "  - Media: " << avg_distance 
             << " gradi (" << (avg_distance * 3600) << " arcsec)\n";
        cout << "  - Limite corridor: " << corridor_width_deg 
             << " gradi (" << (corridor_width_deg * 3600) << " arcsec)\n\n";
        
        if (num_invalid > 0) {
            cout << "⚠️  ATTENZIONE: " << num_invalid << " STELLE FUORI DAL CORRIDOR!\n\n";
            cout << "Prime 10 violazioni:\n";
            cout << "  ID             RA          Dec         Mag    Distanza[gradi]  Distanza[arcsec]  Eccesso[arcsec]\n";
            cout << "  " << string(110, '-') << "\n";
            
            for (size_t i = 0; i < min(violations.size(), size_t(10)); ++i) {
                const auto& star = violations[i].first;
                double dist = violations[i].second;
                double excess = (dist - corridor_width_deg) * 3600; // arcsec oltre il limite
                
                cout << "  " << setw(14) << star.source_id 
                     << "  " << setw(10) << fixed << setprecision(6) << star.ra
                     << "  " << setw(10) << star.dec
                     << "  " << setw(5) << setprecision(2) << star.phot_g_mean_mag
                     << "  " << setw(15) << setprecision(6) << dist
                     << "  " << setw(16) << setprecision(2) << (dist * 3600)
                     << "  " << setw(16) << setprecision(2) << excess << "\n";
            }
            cout << "\n";
        } else {
            cout << "✓ SUCCESSO: Tutte le stelle sono entro il corridor specificato!\n\n";
        }
        
        // === 8. SALVA RISULTATI IN FILE CSV ===
        string csv_path = "/tmp/corridor_validation_17030.csv";
        ofstream csv(csv_path);
        csv << "source_id,ra,dec,mag,min_distance_deg,min_distance_arcsec,within_corridor\n";
        
        for (const auto& star : stars) {
            double min_dist = distance_point_to_path(star.ra, star.dec, path);
            bool within = (min_dist <= corridor_width_deg);
            
            csv << star.source_id << ","
                << fixed << setprecision(8) << star.ra << ","
                << star.dec << ","
                << setprecision(3) << star.phot_g_mean_mag << ","
                << setprecision(8) << min_dist << ","
                << setprecision(3) << (min_dist * 3600) << ","
                << (within ? "YES" : "NO") << "\n";
        }
        csv.close();
        
        cout << "Risultati salvati in: " << csv_path << "\n\n";
        
        // === 9. CONCLUSIONE ===
        double success_rate = (num_valid * 100.0) / stars.size();
        
        cout << "========================================\n";
        cout << "TASSO DI SUCCESSO: " << fixed << setprecision(2) << success_rate << "%\n";
        cout << "========================================\n\n";
        
        // Shutdown catalog
        UnifiedGaiaCatalog::shutdown();
        
        return (num_invalid == 0) ? 0 : 1;
        
    } catch (const exception& e) {
        cerr << "❌ ERRORE: " << e.what() << endl;
        return 1;
    }
}
