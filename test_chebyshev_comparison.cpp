/**
 * @file test_chebyshev_comparison.cpp
 * @brief Confronto 16 vs 32 vs 48 punti di controllo
 * @date 4 Dicembre 2025
 * 
 * Test comparativo per valutare:
 * - Affidabilità (numero candidati trovati)
 * - Performance (tempo di esecuzione)
 * - Trade-off ottimale
 */

#include "phase1_candidate_screening.h"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/io/parsers/OrbFitEQ1Parser.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Propagator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/Integrator.hpp"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/ephemeris/PlanetaryEphemeris.hpp"
#include "ioc_gaialib/unified_gaia_catalog.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>
#include <Eigen/Dense>

using namespace std;

// Costanti
constexpr double MJD_TO_JD = 2400000.5;
constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;
constexpr double RAD_TO_ARCSEC = RAD_TO_DEG * 3600.0;
constexpr double EPSILON_J2000 = 23.4392911 * DEG_TO_RAD;

struct Coord {
    double ra_rad;
    double dec_rad;
};

struct MinResult {
    double min_distance_arcsec;
    double min_time_mjd;
};

// Forward declarations
vector<double> fitChebyshev1D(const vector<double>& times, const vector<double>& values,
                             double t_min, double t_max, int degree);
double evaluateChebyshev(double t, double t_min, double t_max, const vector<double>& coeffs);
double evaluateChebyshevDerivative(double t, double t_min, double t_max, const vector<double>& coeffs);
Coord asteroidPosition(double mjd, double mjd_min, double mjd_max,
                       const vector<double>& ra_coeffs, const vector<double>& dec_coeffs);
double objectiveFunction(double mjd, const Coord& star_pos, double mjd_min, double mjd_max,
                        const vector<double>& ra_coeffs, const vector<double>& dec_coeffs);
double objectiveFunctionDerivative(double mjd, const Coord& star_pos, 
                                   double mjd_min, double mjd_max,
                                   const vector<double>& ra_coeffs, 
                                   const vector<double>& dec_coeffs);
MinResult findMinimumDistance(const Coord& star_pos, double mjd_min, double mjd_max,
                             const vector<double>& ra_coeffs, const vector<double>& dec_coeffs,
                             int coarse_steps);
void runTest(int NUM_CONTROL);

// ===================================================================
// IMPLEMENTAZIONE FUNZIONI
// ===================================================================

vector<double> fitChebyshev1D(const vector<double>& times, const vector<double>& values,
                             double t_min, double t_max, int degree) {
    int n = times.size();
    vector<double> coeffs(degree + 1, 0.0);
    
    vector<double> x_norm(n);
    for (int i = 0; i < n; ++i) {
        x_norm[i] = (2.0 * times[i] - (t_max + t_min)) / (t_max - t_min);
    }
    
    for (int k = 0; k <= degree; ++k) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            double T_k = cos(k * acos(x_norm[i]));
            sum += values[i] * T_k;
        }
        coeffs[k] = (k == 0) ? sum / n : 2.0 * sum / n;
    }
    
    return coeffs;
}

double evaluateChebyshev(double t, double t_min, double t_max, const vector<double>& coeffs) {
    double x = (2.0 * t - (t_max + t_min)) / (t_max - t_min);
    x = max(-1.0, min(1.0, x));
    
    int n = coeffs.size();
    if (n == 0) return 0.0;
    if (n == 1) return coeffs[0];
    
    double b_k_plus_2 = 0.0;
    double b_k_plus_1 = 0.0;
    
    for (int k = n - 1; k >= 1; --k) {
        double b_k = 2.0 * x * b_k_plus_1 - b_k_plus_2 + coeffs[k];
        b_k_plus_2 = b_k_plus_1;
        b_k_plus_1 = b_k;
    }
    
    return x * b_k_plus_1 - b_k_plus_2 + coeffs[0];
}

double evaluateChebyshevDerivative(double t, double t_min, double t_max, const vector<double>& coeffs) {
    int n = coeffs.size();
    if (n <= 1) return 0.0;
    
    vector<double> deriv_coeffs(n - 1);
    for (int k = n - 2; k >= 0; --k) {
        deriv_coeffs[k] = 0.0;
        for (int j = k + 1; j < n; j += 2) {
            deriv_coeffs[k] += 2.0 * j * coeffs[j];
        }
    }
    
    double dx_dt = 2.0 / (t_max - t_min);
    return evaluateChebyshev(t, t_min, t_max, deriv_coeffs) * dx_dt;
}

Coord asteroidPosition(double mjd, double mjd_min, double mjd_max,
                       const vector<double>& ra_coeffs, const vector<double>& dec_coeffs) {
    Coord pos;
    pos.ra_rad = evaluateChebyshev(mjd, mjd_min, mjd_max, ra_coeffs);
    pos.dec_rad = evaluateChebyshev(mjd, mjd_min, mjd_max, dec_coeffs);
    return pos;
}

double objectiveFunction(double mjd, const Coord& star_pos, double mjd_min, double mjd_max,
                        const vector<double>& ra_coeffs, const vector<double>& dec_coeffs) {
    Coord ast_pos = asteroidPosition(mjd, mjd_min, mjd_max, ra_coeffs, dec_coeffs);
    
    double sin_dec1 = sin(star_pos.dec_rad);
    double cos_dec1 = cos(star_pos.dec_rad);
    double sin_dec2 = sin(ast_pos.dec_rad);
    double cos_dec2 = cos(ast_pos.dec_rad);
    double cos_dra = cos(star_pos.ra_rad - ast_pos.ra_rad);
    
    double cos_theta = sin_dec1 * sin_dec2 + cos_dec1 * cos_dec2 * cos_dra;
    cos_theta = max(-1.0, min(1.0, cos_theta));
    
    return 1.0 - cos_theta;
}

double objectiveFunctionDerivative(double mjd, const Coord& star_pos, 
                                   double mjd_min, double mjd_max,
                                   const vector<double>& ra_coeffs, 
                                   const vector<double>& dec_coeffs) {
    Coord ast_pos = asteroidPosition(mjd, mjd_min, mjd_max, ra_coeffs, dec_coeffs);
    
    double dra_dt = evaluateChebyshevDerivative(mjd, mjd_min, mjd_max, ra_coeffs);
    double ddec_dt = evaluateChebyshevDerivative(mjd, mjd_min, mjd_max, dec_coeffs);
    
    double sin_dec1 = sin(star_pos.dec_rad);
    double cos_dec1 = cos(star_pos.dec_rad);
    double sin_dec2 = sin(ast_pos.dec_rad);
    double cos_dec2 = cos(ast_pos.dec_rad);
    double dra = star_pos.ra_rad - ast_pos.ra_rad;
    double cos_dra = cos(dra);
    double sin_dra = sin(dra);
    
    double dcos_theta_dra2 = -cos_dec1 * cos_dec2 * sin_dra;
    double dcos_theta_ddec2 = -sin_dec1 * cos_dec2 + cos_dec1 * sin_dec2 * cos_dra;
    
    double dcos_theta_dt = dcos_theta_dra2 * (-dra_dt) + dcos_theta_ddec2 * ddec_dt;
    
    return -dcos_theta_dt;
}

MinResult findMinimumDistance(const Coord& star_pos, double mjd_min, double mjd_max,
                             const vector<double>& ra_coeffs, const vector<double>& dec_coeffs,
                             int coarse_steps) {
    double step_size = (mjd_max - mjd_min) / coarse_steps;
    double min_objective = 2.0;
    double mjd_coarse = mjd_min;
    
    for (int i = 0; i <= coarse_steps; ++i) {
        double mjd = mjd_min + i * step_size;
        double current_objective = objectiveFunction(mjd, star_pos, mjd_min, mjd_max,
                                                     ra_coeffs, dec_coeffs);
        
        if (current_objective < min_objective) {
            min_objective = current_objective;
            mjd_coarse = mjd;
        }
    }
    
    double mjd_refined = mjd_coarse;
    constexpr int MAX_NEWTON_ITER = 10;
    constexpr double TOLERANCE = 1e-9;
    
    for (int iter = 0; iter < MAX_NEWTON_ITER; ++iter) {
        double f_deriv = objectiveFunctionDerivative(mjd_refined, star_pos, 
                                                      mjd_min, mjd_max,
                                                      ra_coeffs, dec_coeffs);
        
        if (abs(f_deriv) < TOLERANCE) {
            break;
        }
        
        double h = 1e-6;
        double f_deriv_plus = objectiveFunctionDerivative(mjd_refined + h, star_pos,
                                                          mjd_min, mjd_max,
                                                          ra_coeffs, dec_coeffs);
        double f_second_deriv = (f_deriv_plus - f_deriv) / h;
        
        if (abs(f_second_deriv) < 1e-12) {
            break;
        }
        
        double delta = -f_deriv / f_second_deriv;
        
        if (abs(delta) < TOLERANCE) {
            break;
        }
        
        mjd_refined += delta;
        mjd_refined = max(mjd_min, min(mjd_max, mjd_refined));
    }
    
    double final_objective = objectiveFunction(mjd_refined, star_pos, mjd_min, mjd_max,
                                               ra_coeffs, dec_coeffs);
    double cos_theta_min = 1.0 - final_objective;
    cos_theta_min = max(-1.0, min(1.0, cos_theta_min));
    
    MinResult result;
    result.min_distance_arcsec = acos(cos_theta_min) * RAD_TO_ARCSEC;
    result.min_time_mjd = mjd_refined;
    
    return result;
}

// ===================================================================
// FUNZIONE DI TEST PARAMETRIZZATA
// ===================================================================
void runTest(int NUM_CONTROL) {
    try {
        // Ottieni istanza catalogo
        auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
        
        // Carica elementi orbitali
        cout << "\n[1] Caricamento elementi orbitali 17030...\n";
        astdyn::io::parsers::OrbFitEQ1Parser parser;
        auto elements = parser.parse("17030_astdys.eq1");
        
        astdyn::propagation::KeplerianElements kep;
        kep.semi_major_axis = elements.semi_major_axis;
        kep.eccentricity = elements.eccentricity;
        kep.inclination = elements.inclination;
        kep.longitude_ascending_node = elements.longitude_asc_node;
        kep.argument_perihelion = elements.argument_perihelion;
        kep.mean_anomaly = elements.mean_anomaly;
        kep.epoch_mjd_tdb = elements.epoch_mjd_tdb;
        kep.gravitational_parameter = 1.32712440018e20 / pow(1.495978707e11, 3) * pow(86400.0, 2);
        cout << "✓ Elementi caricati\n\n";
        
        // Setup propagatore
        auto ephemeris = make_shared<astdyn::ephemeris::PlanetaryEphemeris>();
        auto integrator = make_unique<astdyn::propagation::RKF78Integrator>(0.5, 1e-9);
        astdyn::propagation::PropagatorSettings settings;
        settings.include_planets = false;
        auto propagator = make_unique<astdyn::propagation::Propagator>(
            move(integrator), ephemeris, settings);
        
        // STEP 1: Propaga N punti di controllo
        cout << "[2] Propagazione " << NUM_CONTROL << " punti di controllo...\n";
        double mjd_start = 61007.0;
        double mjd_end = 61008.0;
        
        vector<double> control_times;
        vector<double> control_ra, control_dec;
        
        auto t_prop_start = chrono::high_resolution_clock::now();
        
        for (int i = 0; i < NUM_CONTROL; ++i) {
            double fraction = static_cast<double>(i) / (NUM_CONTROL - 1);
            double mjd = mjd_start + fraction * (mjd_end - mjd_start);
            
            // Propaga posizione keplerian → cartesian (barycentric ecliptic)
            auto kep_prop = propagator->propagate_keplerian(kep, mjd);
            auto cart_ecl = astdyn::propagation::keplerian_to_cartesian(kep_prop);
            
            double x_ecl = cart_ecl.position[0];
            double y_ecl = cart_ecl.position[1];
            double z_ecl = cart_ecl.position[2];
            
            // Converti eclittica → equatoriale J2000
            double x_eq = x_ecl;
            double y_eq = y_ecl * cos(EPSILON_J2000) - z_ecl * sin(EPSILON_J2000);
            double z_eq = y_ecl * sin(EPSILON_J2000) + z_ecl * cos(EPSILON_J2000);
            
            // Ottieni posizione Terra (in realtà ci serve geocentrico, ma per 
            // Phase 1 usiamo approssimazione eliocentrica)
            double sun_dist = sqrt(x_eq*x_eq + y_eq*y_eq + z_eq*z_eq);
            double earth_x = -x_eq / sun_dist;
            double earth_y = -y_eq / sun_dist;
            double earth_z = -z_eq / sun_dist;
            
            double ra = atan2(earth_y, earth_x);
            double dec = atan2(earth_z, sqrt(earth_x*earth_x + earth_y*earth_y));
            
            if (ra < 0) ra += 2.0 * M_PI;
            
            control_times.push_back(mjd);
            control_ra.push_back(ra);
            control_dec.push_back(dec);
        }
        
        auto t_prop_end = chrono::high_resolution_clock::now();
        double prop_time = chrono::duration<double, milli>(t_prop_end - t_prop_start).count();
        cout << "✓ Propagazione completata in " << prop_time << " ms\n\n";
        
        // STEP 2: Fit Chebyshev
        cout << "[3] Fitting polinomi Chebyshev (grado " << (NUM_CONTROL - 1) << ")...\n";
        auto t_fit_start = chrono::high_resolution_clock::now();
        
        int degree = NUM_CONTROL - 1;
        auto ra_coeffs = fitChebyshev1D(control_times, control_ra, mjd_start, mjd_end, degree);
        auto dec_coeffs = fitChebyshev1D(control_times, control_dec, mjd_start, mjd_end, degree);
        
        auto t_fit_end = chrono::high_resolution_clock::now();
        double fit_time = chrono::duration<double, milli>(t_fit_end - t_fit_start).count();
        cout << "✓ Fitting completato in " << fit_time << " ms\n\n";
        
        // STEP 3: Query corridor
        cout << "[4] Query corridor (width=30 arcsec, mag<18)...\n";
        
        ioc::gaia::CorridorQueryParams params;
        for (int i = 0; i < NUM_CONTROL; ++i) {
            params.path.push_back(ioc::gaia::CelestialPoint(
                control_ra[i] * RAD_TO_DEG, 
                control_dec[i] * RAD_TO_DEG));
        }
        params.width = 0.0083;
        params.max_magnitude = 18.0;
        
        auto t_query_start = chrono::high_resolution_clock::now();
        auto stars = catalog.queryCorridor(params);
        auto t_query_end = chrono::high_resolution_clock::now();
        double query_time = chrono::duration<double, milli>(t_query_end - t_query_start).count();
        
        cout << "✓ Trovate " << stars.size() << " stelle in " << query_time << " ms\n\n";
        
        // STEP 4: Ottimizza closest approach
        cout << "[5] Ottimizzazione CA (1k grid + Newton-Raphson)...\n";
        auto t_opt_start = chrono::high_resolution_clock::now();
        
        vector<pair<MinResult, ioc::gaia::GaiaStar>> results;
        for (const auto& star : stars) {
            Coord star_coord;
            star_coord.ra_rad = star.ra * DEG_TO_RAD;
            star_coord.dec_rad = star.dec * DEG_TO_RAD;
            
            MinResult min_res = findMinimumDistance(star_coord, mjd_start, mjd_end,
                                                    ra_coeffs, dec_coeffs, 1000);
            results.push_back({min_res, star});
        }
        
        auto t_opt_end = chrono::high_resolution_clock::now();
        double opt_time = chrono::duration<double, milli>(t_opt_end - t_opt_start).count();
        
        cout << "✓ Ottimizzazione completata in " << opt_time << " ms\n";
        cout << "  Tempo medio per stella: " << (opt_time / stars.size()) << " ms\n\n";
        
        // Ordina per distanza
        sort(results.begin(), results.end(),
             [](const auto& a, const auto& b) {
                 return a.first.min_distance_arcsec < b.first.min_distance_arcsec;
             });
        
        // Mostra risultati
        cout << "══════════════════════════════════════════════════════════\n";
        cout << "CANDIDATI con CA < 15 arcsec:\n\n";
        
        int count = 0;
        bool found_target = false;
        for (const auto& [min_res, star] : results) {
            if (min_res.min_distance_arcsec < 15.0) {
                count++;
                cout << count << ". Source ID: " << star.source_id << "\n";
                cout << "   RA/Dec: " << fixed << setprecision(6) 
                     << star.ra << "° / " << star.dec << "°\n";
                cout << "   Mag: " << setprecision(2) << star.phot_g_mean_mag << "\n";
                cout << "   Closest Approach: " << setprecision(3) 
                     << min_res.min_distance_arcsec << " arcsec\n";
                cout << "   CA Epoch: MJD " << setprecision(6) 
                     << min_res.min_time_mjd << "\n";
                
                if (star.source_id == 3411546266140512128ULL) {
                    cout << "   ✓✓✓ TARGET STAR TROVATA! ✓✓✓\n";
                    found_target = true;
                }
                cout << "\n";
            }
        }
        
        // Riepilogo
        cout << "══════════════════════════════════════════════════════════\n";
        cout << "RIEPILOGO PERFORMANCE:\n\n";
        cout << "  Punti di controllo:         " << NUM_CONTROL << "\n";
        cout << "  Grado polinomio:            " << degree << "\n";
        cout << "  Propagazione:               " << setprecision(2) << prop_time << " ms\n";
        cout << "  Fitting Chebyshev:          " << fit_time << " ms\n";
        cout << "  Query corridor:             " << query_time << " ms\n";
        cout << "  Ottimizzazione CA:          " << opt_time << " ms\n";
        cout << "  ─────────────────────────────────────────────\n";
        double total_time = prop_time + fit_time + query_time + opt_time;
        cout << "  TOTALE:                     " << total_time << " ms\n\n";
        cout << "  Candidati trovati:          " << count << "\n";
        cout << "  Target star (3411546...):   " << (found_target ? "✓ SÌ" : "✗ NO") << "\n\n";
        
    } catch (const exception& e) {
        cerr << "\n❌ ERRORE: " << e.what() << "\n\n";
    }
}

// ===================================================================
// MAIN
// ===================================================================
int main() {
    cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    cout << "║  TEST COMPARATIVO: 16 vs 32 vs 48 PUNTI DI CONTROLLO    ║\n";
    cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    try {
        // Inizializza catalogo UNA VOLTA SOLA
        cout << "[Setup] Inizializzazione catalogo Gaia...\n";
        string home = string(getenv("HOME"));
        string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
        string json_config = R"({"catalog_type": "multifile_v2","multifile_directory": ")" + catalog_path + R"("})";
        
        if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) {
            throw runtime_error("Impossibile inizializzare catalogo");
        }
        cout << "✓ Catalogo inizializzato\n";
        
        // Test con 3 configurazioni
        vector<int> configs = {16, 32, 48};
        
        for (int num_points : configs) {
            cout << "\n\n" << string(70, '=') << "\n";
            cout << "TEST CON " << num_points << " PUNTI DI CONTROLLO\n";
            cout << string(70, '=') << "\n";
            
            runTest(num_points);
        }
        
        cout << "\n\n" << string(70, '=') << "\n";
        cout << "CONFRONTO COMPLETATO\n";
        cout << string(70, '=') << "\n\n";
        
        return 0;
        
    } catch (const exception& e) {
        cerr << "\n❌ ERRORE: " << e.what() << "\n\n";
        return 1;
    }
}
