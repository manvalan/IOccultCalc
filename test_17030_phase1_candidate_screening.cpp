/**
 * @file test_17030_phase1_candidate_screening.cpp
 * @brief FASE 1: Screening Stelle Candidate per Occultazioni
 * 
 * OBIETTIVO FASE 1:
 * Creare un PATH AD ALTA RISOLUZIONE dell'asteroide e identificare tutte
 * le stelle che hanno un closest approach interessante (potenziali occultazioni).
 * 
 * NON calcoliamo ancora l'occultazione precisa, solo:
 * 1. Path denso (es. ogni 10-30 secondi per 24 ore)
 * 2. Query Corridor con larghezza adeguata (es. 3-5 arcmin)
 * 3. Per ogni stella trovata: calcola closest approach
 * 4. Filtra stelle con closest approach < soglia (es. < 10 arcsec)
 * 
 * FASE 2 (successiva): Per le stelle filtrate, calcolo preciso dell'occultazione
 * 
 * Data: 4 Dicembre 2025
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

// Costanti
constexpr double MJD_TO_JD = 2400000.5;
constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;
constexpr double EPSILON_J2000 = 23.4392911 * DEG_TO_RAD;

// Struttura per punto del path
struct PathPoint {
    double mjd_tdb;
    double ra_deg;
    double dec_deg;
    Vector3d pos_geo_au;  // Posizione geocentrica [AU]
};

// Struttura per stella candidata
struct CandidateStar {
    uint64_t source_id;
    double ra_deg;
    double dec_deg;
    double mag;
    double closest_approach_arcsec;
    double closest_approach_mjd;
    int closest_segment_index;
};

// Conversione UTC a MJD TDB
double utc_to_mjd_tdb(int year, int month, int day, int hour, int min, double sec) {
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    double jd = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    jd += (hour - 12) / 24.0 + min / 1440.0 + sec / 86400.0;
    double mjd_utc = jd - 2400000.5;
    return mjd_utc + 69.184 / 86400.0;
}

// Conversione eclittica -> equatoriale ICRF
Vector3d ecliptic_to_equatorial(const Vector3d& ecl) {
    double cos_eps = cos(EPSILON_J2000);
    double sin_eps = sin(EPSILON_J2000);
    return Vector3d(ecl[0], cos_eps * ecl[1] - sin_eps * ecl[2], 
                    sin_eps * ecl[1] + cos_eps * ecl[2]);
}

// Conversione cartesiano -> RA/Dec
void cartesian_to_radec(const Vector3d& pos, double& ra, double& dec) {
    double r = pos.norm();
    dec = asin(pos[2] / r);
    ra = atan2(pos[1], pos[0]);
    if (ra < 0) ra += 2 * M_PI;
}

// Distanza angolare haversine
double haversine_distance(double ra1_rad, double dec1_rad, double ra2_rad, double dec2_rad) {
    double delta_ra = ra2_rad - ra1_rad;
    double delta_dec = dec2_rad - dec1_rad;
    double a = sin(delta_dec / 2) * sin(delta_dec / 2) +
               cos(dec1_rad) * cos(dec2_rad) * sin(delta_ra / 2) * sin(delta_ra / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return c * RAD_TO_DEG * 3600.0; // arcsec
}

// Calcola closest approach di una stella a un segmento del path
double closest_approach_to_segment(double star_ra_deg, double star_dec_deg,
                                    const PathPoint& p1, const PathPoint& p2,
                                    double& closest_mjd) {
    // Converti tutto in radianti
    double star_ra = star_ra_deg * DEG_TO_RAD;
    double star_dec = star_dec_deg * DEG_TO_RAD;
    double ra1 = p1.ra_deg * DEG_TO_RAD;
    double dec1 = p1.dec_deg * DEG_TO_RAD;
    double ra2 = p2.ra_deg * DEG_TO_RAD;
    double dec2 = p2.dec_deg * DEG_TO_RAD;
    
    // Distanze dai due estremi
    double dist1 = haversine_distance(star_ra, star_dec, ra1, dec1);
    double dist2 = haversine_distance(star_ra, star_dec, ra2, dec2);
    
    // Lunghezza segmento
    double seg_length = haversine_distance(ra1, dec1, ra2, dec2);
    
    if (seg_length < 0.01) { // Segmento molto corto
        closest_mjd = p1.mjd_tdb;
        return dist1;
    }
    
    // Converti in coordinate cartesiane 3D (sfera unitaria)
    double x1 = cos(dec1) * cos(ra1);
    double y1 = cos(dec1) * sin(ra1);
    double z1 = sin(dec1);
    
    double x2 = cos(dec2) * cos(ra2);
    double y2 = cos(dec2) * sin(ra2);
    double z2 = sin(dec2);
    
    double xs = cos(star_dec) * cos(star_ra);
    double ys = cos(star_dec) * sin(star_ra);
    double zs = sin(star_dec);
    
    // Vettore segmento e vettore da p1 a stella
    double vx = x2 - x1, vy = y2 - y1, vz = z2 - z1;
    double wx = xs - x1, wy = ys - y1, wz = zs - z1;
    
    // Parametro t della proiezione
    double dot_vw = vx*wx + vy*wy + vz*wz;
    double dot_vv = vx*vx + vy*vy + vz*vz;
    double t = dot_vw / dot_vv;
    
    if (t <= 0.0) {
        closest_mjd = p1.mjd_tdb;
        return dist1;
    }
    if (t >= 1.0) {
        closest_mjd = p2.mjd_tdb;
        return dist2;
    }
    
    // Interpola il tempo
    closest_mjd = p1.mjd_tdb + t * (p2.mjd_tdb - p1.mjd_tdb);
    
    // Punto proiettato
    double xp = x1 + t * vx;
    double yp = y1 + t * vy;
    double zp = z1 + t * vz;
    
    // Normalizza
    double norm = sqrt(xp*xp + yp*yp + zp*zp);
    xp /= norm; yp /= norm; zp /= norm;
    
    // Calcola distanza
    double dx = xs - xp, dy = ys - yp, dz = zs - zp;
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    return asin(min(1.0, dist / 2.0)) * 2.0 * RAD_TO_DEG * 3600.0; // arcsec
}

int main() {
    cout << "\n========================================\n";
    cout << "FASE 1: SCREENING STELLE CANDIDATE\n";
    cout << "Asteroide: (17030) 1999 CQ3\n";
    cout << "Data: 28 Novembre 2025\n";
    cout << "========================================\n\n";
    
    try {
        // === STEP 1: CARICA ELEMENTI ORBITALI ===
        cout << "STEP 1: Caricamento elementi orbitali...\n";
        string eq1_path = "../17030_astdys.eq1";
        io::parsers::OrbFitEQ1Parser parser;
        auto elements = parser.parse(eq1_path);
        
        cout << "  ✓ Semi-asse: " << elements.semi_major_axis << " AU\n";
        cout << "  ✓ Eccentricità: " << elements.eccentricity << "\n\n";
        
        // === STEP 2: CREA PROPAGATORE ===
        cout << "STEP 2: Inizializzazione propagatore RKF78...\n";
        
        propagation::KeplerianElements kep;
        kep.semi_major_axis = elements.semi_major_axis;
        kep.eccentricity = elements.eccentricity;
        kep.inclination = elements.inclination;
        kep.longitude_ascending_node = elements.longitude_asc_node;
        kep.argument_perihelion = elements.argument_perihelion;
        kep.mean_anomaly = elements.mean_anomaly;
        kep.epoch_mjd_tdb = elements.epoch_mjd_tdb;
        kep.gravitational_parameter = 1.32712440018e20 / pow(1.495978707e11, 3) * 
                                      pow(86400.0, 2);
        
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
        
        propagation::Propagator propagator(move(integrator), ephemeris, settings);
        cout << "  ✓ Propagatore configurato con tutte le perturbazioni\n\n";
        
        // === STEP 3: CREA PATH AD ALTA RISOLUZIONE ===
        cout << "STEP 3: Creazione path ad alta risoluzione...\n";
        cout << "  Periodo: 28 Nov 2025, 00:00 - 24:00 UT\n";
        
        // PARAMETRO CRITICO: Intervallo tra punti del path
        // Per screening iniziale: 30-60 secondi è sufficiente
        // Per calcolo occultazione precisa (Fase 2): 1-5 secondi
        int interval_seconds = 30;
        int total_seconds = 24 * 3600;
        int num_points = total_seconds / interval_seconds + 1;
        
        cout << "  Intervallo: " << interval_seconds << " secondi\n";
        cout << "  Punti totali: " << num_points << "\n";
        
        vector<PathPoint> path;
        path.reserve(num_points);
        
        auto t_start = chrono::high_resolution_clock::now();
        
        for (int i = 0; i < num_points; ++i) {
            double seconds = i * interval_seconds;
            int hour = static_cast<int>(seconds / 3600);
            int min = static_cast<int>((seconds - hour * 3600) / 60);
            double sec = seconds - hour * 3600 - min * 60;
            
            double mjd_tdb = utc_to_mjd_tdb(2025, 11, 28, hour, min, sec);
            
            // Propaga asteroide (barycentric ecliptic)
            auto kep_prop = propagator.propagate_keplerian(kep, mjd_tdb);
            auto cart_ecl = propagation::keplerian_to_cartesian(kep_prop);
            Vector3d ast_bary_icrf = ecliptic_to_equatorial(cart_ecl.position);
            
            // Posizione Terra
            double jd_tdb = mjd_tdb + MJD_TO_JD;
            Vector3d earth_helio_ecl = ephemeris::PlanetaryEphemeris::getPosition(
                ephemeris::CelestialBody::EARTH, jd_tdb);
            Vector3d sun_bary_ecl = ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd_tdb);
            Vector3d earth_bary_icrf = ecliptic_to_equatorial(earth_helio_ecl - sun_bary_ecl);
            
            // Geocentrico
            Vector3d ast_geo_icrf = ast_bary_icrf - earth_bary_icrf;
            
            // RA/Dec
            double ra_rad, dec_rad;
            cartesian_to_radec(ast_geo_icrf, ra_rad, dec_rad);
            
            PathPoint pt;
            pt.mjd_tdb = mjd_tdb;
            pt.ra_deg = ra_rad * RAD_TO_DEG;
            pt.dec_deg = dec_rad * RAD_TO_DEG;
            pt.pos_geo_au = ast_geo_icrf;
            
            path.push_back(pt);
        }
        
        auto t_end = chrono::high_resolution_clock::now();
        auto duration_ms = chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();
        
        cout << "  ✓ Path generato in " << duration_ms << " ms (" 
             << fixed << setprecision(2) << (duration_ms / 1000.0 / num_points * 1000) 
             << " ms/punto)\n";
        cout << "  ✓ Range RA:  " << setprecision(6) << path.front().ra_deg 
             << "° → " << path.back().ra_deg << "°\n";
        cout << "  ✓ Range Dec: " << path.front().dec_deg 
             << "° → " << path.back().dec_deg << "°\n\n";
        
        // === STEP 4: QUERY CORRIDOR PER STELLE CANDIDATE ===
        cout << "STEP 4: Query Corridor per stelle candidate...\n";
        
        string home = string(getenv("HOME"));
        string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
        
        string json_config = R"({
            "catalog_type": "multifile_v2",
            "multifile_directory": ")" + catalog_path + R"("
        })";
        
        if (!UnifiedGaiaCatalog::initialize(json_config)) {
            throw runtime_error("Impossibile inizializzare catalogo");
        }
        
        auto& catalog = UnifiedGaiaCatalog::getInstance();
        
        // PARAMETRO CRITICO: Larghezza corridor per screening
        // Troppo stretto: perdi candidate valide
        // Troppo largo: troppo tempo di calcolo nella Fase 2
        // Raccomandato: 2-5 arcmin (0.033 - 0.083 gradi)
        double corridor_width_deg = 0.05; // 3 arcmin = 180 arcsec
        double max_magnitude = 18.0;
        
        cout << "  Larghezza corridor: " << corridor_width_deg << "° (" 
             << (corridor_width_deg * 60) << " arcmin)\n";
        cout << "  Magnitudine limite: " << max_magnitude << "\n";
        
        // Costruisci parametri per queryCorridor
        CorridorQueryParams params;
        for (const auto& pt : path) {
            params.path.push_back(CelestialPoint(pt.ra_deg, pt.dec_deg));
        }
        params.width = corridor_width_deg;
        params.max_magnitude = max_magnitude;
        params.min_parallax = -1.0;
        params.max_results = 0;
        
        t_start = chrono::high_resolution_clock::now();
        auto stars = catalog.queryCorridor(params);
        t_end = chrono::high_resolution_clock::now();
        duration_ms = chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();
        
        cout << "  ✓ Query completata in " << duration_ms << " ms\n";
        cout << "  ✓ Stelle trovate: " << stars.size() << "\n\n";
        
        // === STEP 5: CALCOLA CLOSEST APPROACH PER OGNI STELLA ===
        cout << "STEP 5: Calcolo closest approach per ogni stella...\n";
        
        vector<CandidateStar> candidates;
        candidates.reserve(stars.size());
        
        t_start = chrono::high_resolution_clock::now();
        
        for (const auto& star : stars) {
            double min_distance = 1e10;
            double best_mjd = 0;
            int best_segment = -1;
            
            // Itera su tutti i segmenti del path
            for (size_t i = 0; i < path.size() - 1; ++i) {
                double closest_mjd;
                double dist = closest_approach_to_segment(star.ra, star.dec,
                                                          path[i], path[i+1],
                                                          closest_mjd);
                if (dist < min_distance) {
                    min_distance = dist;
                    best_mjd = closest_mjd;
                    best_segment = i;
                }
            }
            
            CandidateStar candidate;
            candidate.source_id = star.source_id;
            candidate.ra_deg = star.ra;
            candidate.dec_deg = star.dec;
            candidate.mag = star.phot_g_mean_mag;
            candidate.closest_approach_arcsec = min_distance;
            candidate.closest_approach_mjd = best_mjd;
            candidate.closest_segment_index = best_segment;
            
            candidates.push_back(candidate);
        }
        
        t_end = chrono::high_resolution_clock::now();
        duration_ms = chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();
        
        cout << "  ✓ Closest approach calcolato in " << duration_ms << " ms\n";
        cout << "  ✓ Tempo medio: " << fixed << setprecision(3) 
             << (duration_ms / double(stars.size())) << " ms/stella\n\n";
        
        // === STEP 6: FILTRA CANDIDATE INTERESSANTI ===
        cout << "STEP 6: Filtro stelle con closest approach interessante...\n";
        
        // PARAMETRO CRITICO: Soglia di closest approach
        // Per occultazioni: tipicamente < 2-5 arcsec
        // Per near-miss: fino a 10-20 arcsec
        double threshold_arcsec = 10.0;
        
        cout << "  Soglia closest approach: " << threshold_arcsec << " arcsec\n\n";
        
        // Ordina per closest approach
        sort(candidates.begin(), candidates.end(),
             [](const CandidateStar& a, const CandidateStar& b) {
                 return a.closest_approach_arcsec < b.closest_approach_arcsec;
             });
        
        vector<CandidateStar> filtered;
        for (const auto& c : candidates) {
            if (c.closest_approach_arcsec <= threshold_arcsec) {
                filtered.push_back(c);
            }
        }
        
        cout << "========================================\n";
        cout << "RISULTATI FASE 1\n";
        cout << "========================================\n\n";
        
        cout << "Stelle totali nel corridor: " << candidates.size() << "\n";
        cout << "Stelle candidate (CA < " << threshold_arcsec << " arcsec): " 
             << filtered.size() << "\n\n";
        
        if (filtered.size() > 0) {
            cout << "TOP 10 CANDIDATE per closest approach:\n";
            cout << "  #  Source ID            RA         Dec        Mag   CA[arcsec]  MJD\n";
            cout << "  -------------------------------------------------------------------------\n";
            
            for (size_t i = 0; i < min(filtered.size(), size_t(10)); ++i) {
                const auto& c = filtered[i];
                cout << "  " << setw(2) << (i+1)
                     << "  " << setw(19) << c.source_id
                     << "  " << setw(9) << fixed << setprecision(5) << c.ra_deg
                     << "  " << setw(9) << c.dec_deg
                     << "  " << setw(5) << setprecision(2) << c.mag
                     << "  " << setw(10) << setprecision(3) << c.closest_approach_arcsec
                     << "  " << setprecision(6) << c.closest_approach_mjd << "\n";
            }
            cout << "\n";
            
            cout << "PROSSIMI PASSI (FASE 2):\n";
            cout << "  1. Per ciascuna stella candidata, creare path ultra-denso\n";
            cout << "     nell'intervallo ±5 minuti attorno al closest approach\n";
            cout << "  2. Propagare con step di 1-5 secondi per calcolo preciso\n";
            cout << "  3. Calcolare geometria occultazione esatta\n";
            cout << "  4. Verificare se l'ombra passa sulla stella\n";
            cout << "  5. Calcolare chord, durata, e path dell'ombra\n\n";
        } else {
            cout << "Nessuna stella candidata trovata con CA < " << threshold_arcsec << " arcsec\n\n";
            cout << "Stelle più vicine trovate:\n";
            cout << "  #  Source ID            RA         Dec        Mag   CA[arcsec]\n";
            cout << "  -------------------------------------------------------------------\n";
            
            for (size_t i = 0; i < min(candidates.size(), size_t(5)); ++i) {
                const auto& c = candidates[i];
                cout << "  " << setw(2) << (i+1)
                     << "  " << setw(19) << c.source_id
                     << "  " << setw(9) << fixed << setprecision(5) << c.ra_deg
                     << "  " << setw(9) << c.dec_deg
                     << "  " << setw(5) << setprecision(2) << c.mag
                     << "  " << setw(10) << setprecision(3) << c.closest_approach_arcsec << "\n";
            }
            cout << "\n";
        }
        
        // === SALVA RISULTATI ===
        string csv_path = "/tmp/phase1_candidates_17030.csv";
        ofstream csv(csv_path);
        csv << "rank,source_id,ra,dec,mag,closest_approach_arcsec,closest_approach_mjd,segment_index\n";
        
        for (size_t i = 0; i < filtered.size(); ++i) {
            const auto& c = filtered[i];
            csv << (i+1) << "," << c.source_id << "," 
                << fixed << setprecision(8) << c.ra_deg << "," << c.dec_deg << ","
                << setprecision(3) << c.mag << ","
                << setprecision(6) << c.closest_approach_arcsec << ","
                << setprecision(8) << c.closest_approach_mjd << ","
                << c.closest_segment_index << "\n";
        }
        csv.close();
        
        cout << "Risultati salvati in: " << csv_path << "\n\n";
        
        UnifiedGaiaCatalog::shutdown();
        
    } catch (const exception& e) {
        cerr << "❌ ERRORE: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
