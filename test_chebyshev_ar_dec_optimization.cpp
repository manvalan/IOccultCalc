/**
 * @file test_chebyshev_ar_dec_optimization.cpp
 * @brief Test approccio ottimizzato: Chebyshev su RA/Dec + Grid Search
 * @date 4 Dicembre 2025
 * 
 * STRATEGIA:
 * 1. Propaga 16 punti
 * 2. Calcola polinomi Chebyshev per RA e DEC (invece di X,Y,Z)
 * 3. Per ogni stella: ottimizza distanza angolare con grid search
 * 
 * VANTAGGI:
 * - Evita conversione cartesiana -> sferica ad ogni valutazione
 * - Distanza angolare calcolata direttamente su sfera
 * - Grid search trova minimum con precisione arbitraria
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

// Struttura per coordinate sferiche
struct Coord {
    double ra_rad;
    double dec_rad;
};

// Risultato ottimizzazione
struct MinResult {
    double min_distance_arcsec;
    double min_time_mjd;
};

// Forward declaration
void runTest(italoccult::UnifiedGaiaCatalog& catalog, int NUM_CONTROL);

// Conversione eclittica -> equatoriale
Eigen::Vector3d eclipticToEquatorial(const Eigen::Vector3d& ecl) {
    double cos_eps = std::cos(EPSILON_J2000);
    double sin_eps = std::sin(EPSILON_J2000);
    return Eigen::Vector3d(
        ecl[0],
        cos_eps * ecl[1] - sin_eps * ecl[2],
        sin_eps * ecl[1] + cos_eps * ecl[2]
    );
}

// Conversione cartesiano -> RA/Dec
void cartesianToRaDec(const Eigen::Vector3d& pos, double& ra_rad, double& dec_rad) {
    double r = pos.norm();
    dec_rad = std::asin(pos[2] / r);
    ra_rad = std::atan2(pos[1], pos[0]);
    if (ra_rad < 0) ra_rad += 2 * M_PI;
}

// Valutazione polinomio Chebyshev (algoritmo di Clenshaw)
double evaluateChebyshev(double x, const vector<double>& coeffs) {
    if (coeffs.empty()) return 0.0;
    
    double b_k = 0.0;
    double b_k_plus_1 = 0.0;
    
    for (int k = coeffs.size() - 1; k >= 1; --k) {
        double b_k_minus_1 = coeffs[k] + 2.0 * x * b_k - b_k_plus_1;
        b_k_plus_1 = b_k;
        b_k = b_k_minus_1;
    }
    
    return coeffs[0] + x * b_k - b_k_plus_1;
}

// Valutazione derivata del polinomio Chebyshev
// La derivata di T_n(x) è n * U_{n-1}(x) dove U è il polinomio di Chebyshev di seconda specie
double evaluateChebyshevDerivative(double x, const vector<double>& coeffs) {
    if (coeffs.size() <= 1) return 0.0;
    
    // Calcola i coefficienti della derivata
    // dT_0/dx = 0, dT_1/dx = 1, dT_n/dx = n * U_{n-1}(x)
    // Per Chebyshev: d/dx[sum(a_k T_k)] = sum(k * a_k * U_{k-1})
    // Semplificazione: usiamo la relazione di ricorrenza
    
    int n = coeffs.size();
    vector<double> deriv_coeffs(n - 1);
    
    // Per k >= 1: coefficiente della derivata è k * coeffs[k]
    // Ma dobbiamo esprimere in termini di T_k, non U_k
    // Usiamo: T'_n(x) = n * U_{n-1}(x) e U_{n-1}(x) può essere espresso come combinazione di T_k
    
    // Formula semplificata: per small degree usiamo differenze finite
    // Per ora, approccio diretto con formula ricorsiva
    
    deriv_coeffs[n-2] = 2.0 * (n - 1) * coeffs[n - 1];
    if (n > 2) {
        deriv_coeffs[n-3] = 2.0 * (n - 2) * coeffs[n - 2];
    }
    
    for (int k = n - 3; k >= 1; --k) {
        deriv_coeffs[k-1] = deriv_coeffs[k+1] + 2.0 * k * coeffs[k];
    }
    
    if (n > 1) {
        deriv_coeffs[0] = deriv_coeffs[1] / 2.0 + coeffs[1];
    }
    
    return evaluateChebyshev(x, deriv_coeffs);
}

// Posizione asteroide all'istante t
Coord asteroidPosition(double mjd, double mjd_min, double mjd_max,
                       const vector<double>& ra_coeffs,
                       const vector<double>& dec_coeffs) {
    // Normalizza tempo in [-1, 1]
    double x = (2.0 * mjd - (mjd_max + mjd_min)) / (mjd_max - mjd_min);
    
    Coord pos;
    pos.ra_rad = evaluateChebyshev(x, ra_coeffs);
    pos.dec_rad = evaluateChebyshev(x, dec_coeffs);
    
    return pos;
}

// Funzione obiettivo: 1 - cos(theta)
double objectiveFunction(double mjd, const Coord& star_pos,
                        double mjd_min, double mjd_max,
                        const vector<double>& ra_coeffs,
                        const vector<double>& dec_coeffs) {
    
    Coord ast_pos = asteroidPosition(mjd, mjd_min, mjd_max, ra_coeffs, dec_coeffs);
    
    double delta_ra = abs(star_pos.ra_rad - ast_pos.ra_rad);
    
    // Coseno sferico
    double cos_theta = sin(star_pos.dec_rad) * sin(ast_pos.dec_rad) +
                       cos(star_pos.dec_rad) * cos(ast_pos.dec_rad) * cos(delta_ra);
    
    // Clamp a [-1, 1]
    cos_theta = max(-1.0, min(1.0, cos_theta));
    
    return 1.0 - cos_theta;
}

// Derivata della funzione obiettivo rispetto al tempo
double objectiveFunctionDerivative(double mjd, const Coord& star_pos,
                                   double mjd_min, double mjd_max,
                                   const vector<double>& ra_coeffs,
                                   const vector<double>& dec_coeffs) {
    
    // Normalizza tempo
    double x = (2.0 * mjd - (mjd_max + mjd_min)) / (mjd_max - mjd_min);
    
    // Valuta posizione e derivate
    double ra_ast = evaluateChebyshev(x, ra_coeffs);
    double dec_ast = evaluateChebyshev(x, dec_coeffs);
    double dra_dx = evaluateChebyshevDerivative(x, ra_coeffs);
    double ddec_dx = evaluateChebyshevDerivative(x, dec_coeffs);
    
    // Converti da dx/dt
    double dx_dmjd = 2.0 / (mjd_max - mjd_min);
    double dra_dmjd = dra_dx * dx_dmjd;
    double ddec_dmjd = ddec_dx * dx_dmjd;
    
    // Delta RA con segno
    double delta_ra = star_pos.ra_rad - ra_ast;
    
    // Derivata di cos(theta) rispetto a mjd
    // cos(theta) = sin(dec_s)*sin(dec_a) + cos(dec_s)*cos(dec_a)*cos(delta_ra)
    // d/dt[cos(theta)] = sin(dec_s)*cos(dec_a)*ddec/dt + 
    //                    cos(dec_s)*[-sin(dec_a)*ddec/dt*cos(dra) - cos(dec_a)*sin(dra)*(-dra/dt)]
    
    double sin_dec_s = sin(star_pos.dec_rad);
    double cos_dec_s = cos(star_pos.dec_rad);
    double sin_dec_a = sin(dec_ast);
    double cos_dec_a = cos(dec_ast);
    double cos_dra = cos(delta_ra);
    double sin_dra = sin(delta_ra);
    
    double dcos_theta_dmjd = sin_dec_s * cos_dec_a * ddec_dmjd +
                             cos_dec_s * (-sin_dec_a * ddec_dmjd * cos_dra + 
                                          cos_dec_a * sin_dra * dra_dmjd);
    
    // Derivata di (1 - cos(theta))
    return -dcos_theta_dmjd;
}

// APPROCCIO IBRIDO: Grid Search + Newton-Raphson
MinResult findMinimumDistance(const Coord& star_pos,
                              double mjd_min, double mjd_max,
                              const vector<double>& ra_coeffs,
                              const vector<double>& dec_coeffs,
                              int coarse_steps = 1000) {
    
    // FASE 1: Grid search grossolana per identificare regione del minimo
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
    
    // FASE 2: Newton-Raphson refinement partendo dal minimo grossolano
    double mjd_refined = mjd_coarse;
    constexpr int MAX_NEWTON_ITER = 10;
    constexpr double TOLERANCE = 1e-9;  // Tolleranza in MJD (~0.1 millisecondi)
    
    for (int iter = 0; iter < MAX_NEWTON_ITER; ++iter) {
        double f_deriv = objectiveFunctionDerivative(mjd_refined, star_pos, 
                                                      mjd_min, mjd_max,
                                                      ra_coeffs, dec_coeffs);
        
        // Se derivata è quasi zero, siamo al minimo
        if (abs(f_deriv) < TOLERANCE) {
            break;
        }
        
        // Approssimazione derivata seconda con finite difference
        double h = 1e-6;  // 0.1 secondi
        double f_deriv_plus = objectiveFunctionDerivative(mjd_refined + h, star_pos,
                                                          mjd_min, mjd_max,
                                                          ra_coeffs, dec_coeffs);
        double f_second_deriv = (f_deriv_plus - f_deriv) / h;
        
        // Evita divisione per zero
        if (abs(f_second_deriv) < 1e-12) {
            break;
        }
        
        // Newton step: x_new = x_old - f'(x) / f''(x)
        double delta = f_deriv / f_second_deriv;
        mjd_refined -= delta;
        
        // Clamp entro i limiti
        mjd_refined = max(mjd_min, min(mjd_max, mjd_refined));
        
        // Converged?
        if (abs(delta) < TOLERANCE) {
            break;
        }
    }
    
    // Calcola distanza finale
    double final_objective = objectiveFunction(mjd_refined, star_pos, 
                                               mjd_min, mjd_max,
                                               ra_coeffs, dec_coeffs);
    double cos_theta_min = 1.0 - final_objective;
    cos_theta_min = max(-1.0, min(1.0, cos_theta_min));
    
    MinResult result;
    result.min_distance_arcsec = acos(cos_theta_min) * RAD_TO_ARCSEC;
    result.min_time_mjd = mjd_refined;
    
    return result;
}

// Fit polinomio Chebyshev 1D
vector<double> fitChebyshev1D(const vector<double>& times, const vector<double>& values,
                             double t_min, double t_max, int degree) {
    int n = times.size();
    vector<double> coeffs(degree + 1, 0.0);
    
    // Normalizza tempi in [-1, 1]
    vector<double> x_norm(n);
    for (int i = 0; i < n; ++i) {
        x_norm[i] = (2.0 * times[i] - (t_max + t_min)) / (t_max - t_min);
    }
    
    // Calcola coefficienti con proiezione discreta
    for (int k = 0; k <= degree; ++k) {
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            double T_k = cos(k * acos(x_norm[i]));  // T_k(x)
            sum += values[i] * T_k;
        }
        coeffs[k] = (k == 0) ? sum / n : 2.0 * sum / n;
    }
    
    return coeffs;
}

int main() {
    cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    cout << "║  TEST: Chebyshev RA/Dec + Grid Search Optimization      ║\n";
    cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    try {
        // Inizializza catalogo
        cout << "[1] Inizializzazione catalogo Gaia...\n";
        string home = string(getenv("HOME"));
        string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
        string json_config = R"({"catalog_type": "multifile_v2","multifile_directory": ")" + catalog_path + R"("})";
        
        if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) {
            throw runtime_error("Impossibile inizializzare catalogo");
        }
        cout << "✓ Catalogo inizializzato\n\n";
        
        // Carica elementi orbitali
        cout << "[2] Caricamento elementi orbitali 17030...\n";
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
        
        // STEP 1: Propaga 16 punti di controllo
        cout << "[3] Propagazione 16 punti di controllo...\n";
        constexpr int NUM_CONTROL = 16;
        double mjd_start = 61007.0;
        double mjd_end = 61008.0;
        
        vector<double> control_times;
        vector<double> control_ra, control_dec;
        
        auto t_start = chrono::high_resolution_clock::now();
        
        for (int i = 0; i < NUM_CONTROL; ++i) {
            double fraction = static_cast<double>(i) / (NUM_CONTROL - 1);
            double mjd = mjd_start + fraction * (mjd_end - mjd_start);
            
            auto kep_prop = propagator->propagate_keplerian(kep, mjd);
            auto cart_ecl = astdyn::propagation::keplerian_to_cartesian(kep_prop);
            Eigen::Vector3d ast_bary_icrf = eclipticToEquatorial(cart_ecl.position);
            
            double jd = mjd + MJD_TO_JD;
            Eigen::Vector3d earth_helio_ecl = 
                astdyn::ephemeris::PlanetaryEphemeris::getPosition(
                    astdyn::ephemeris::CelestialBody::EARTH, jd);
            Eigen::Vector3d sun_bary_ecl = 
                astdyn::ephemeris::PlanetaryEphemeris::getSunBarycentricPosition(jd);
            Eigen::Vector3d earth_bary_icrf = eclipticToEquatorial(earth_helio_ecl - sun_bary_ecl);
            Eigen::Vector3d ast_geo_icrf = ast_bary_icrf - earth_bary_icrf;
            
            double ra_rad, dec_rad;
            cartesianToRaDec(ast_geo_icrf, ra_rad, dec_rad);
            
            control_times.push_back(mjd);
            control_ra.push_back(ra_rad);
            control_dec.push_back(dec_rad);
        }
        
        auto t_prop = chrono::high_resolution_clock::now();
        double prop_time = chrono::duration<double, milli>(t_prop - t_start).count();
        cout << "✓ Propagazione completata in " << prop_time << " ms\n\n";
        
        // STEP 2: Fit polinomi Chebyshev per RA e DEC
        cout << "[4] Fitting polinomi Chebyshev (grado 15)...\n";
        auto t_fit_start = chrono::high_resolution_clock::now();
        
        vector<double> ra_coeffs = fitChebyshev1D(control_times, control_ra,
                                                   mjd_start, mjd_end, 15);
        vector<double> dec_coeffs = fitChebyshev1D(control_times, control_dec,
                                                    mjd_start, mjd_end, 15);
        
        auto t_fit_end = chrono::high_resolution_clock::now();
        double fit_time = chrono::duration<double, milli>(t_fit_end - t_fit_start).count();
        cout << "✓ Fitting completato in " << fit_time << " ms\n\n";
        
        // STEP 3: Query corridor con 16 punti
        cout << "[5] Query corridor (width=30 arcsec)...\n";
        auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
        
        ioc::gaia::CorridorQueryParams params;
        for (int i = 0; i < NUM_CONTROL; ++i) {
            params.path.push_back(ioc::gaia::CelestialPoint(
                control_ra[i] * RAD_TO_DEG, 
                control_dec[i] * RAD_TO_DEG));
        }
        params.width = 0.0083;  // 30 arcsec
        params.max_magnitude = 18.0;
        
        auto t_query_start = chrono::high_resolution_clock::now();
        auto stars = catalog.queryCorridor(params);
        auto t_query_end = chrono::high_resolution_clock::now();
        double query_time = chrono::duration<double, milli>(t_query_end - t_query_start).count();
        
        cout << "✓ Trovate " << stars.size() << " stelle in " << query_time << " ms\n\n";
        
        // STEP 4: Ottimizza closest approach per ogni stella
        cout << "[6] Ottimizzazione closest approach (Hybrid: 1k grid + Newton-Raphson)...\n";
        auto t_opt_start = chrono::high_resolution_clock::now();
        
        vector<pair<MinResult, ioc::gaia::GaiaStar>> results;
        for (const auto& star : stars) {
            Coord star_coord;
            star_coord.ra_rad = star.ra * DEG_TO_RAD;
            star_coord.dec_rad = star.dec * DEG_TO_RAD;
            
            // Usa 1000 steps per coarse grid (100x meno di prima)
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
        cout << "CANDIDATE con CA < 15 arcsec:\n\n";
        
        int count = 0;
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
                     << min_res.min_time_mjd << "\n\n";
            }
        }
        
        // Riepilogo performance
        cout << "══════════════════════════════════════════════════════════\n";
        cout << "RIEPILOGO PERFORMANCE:\n\n";
        cout << "  Propagazione (16 punti):    " << setprecision(2) << prop_time << " ms\n";
        cout << "  Fitting Chebyshev:          " << fit_time << " ms\n";
        cout << "  Query corridor:             " << query_time << " ms\n";
        cout << "  Ottimizzazione CA:          " << opt_time << " ms\n";
        cout << "  ─────────────────────────────────────────────\n";
        cout << "  TOTALE:                     " << (prop_time + fit_time + query_time + opt_time) 
             << " ms\n\n";
        
        cout << "✅ TEST COMPLETATO!\n\n";
        
    } catch (const exception& e) {
        cerr << "\n❌ ERRORE: " << e.what() << "\n\n";
        return 1;
    }
    
    return 0;
}
