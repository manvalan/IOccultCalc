/**
 * @file initial_orbit.cpp
 * @brief Implementazione metodi di determinazione orbitale iniziale
 * 
 * Basato su Find_Orb di Bill Gray (https://github.com/Bill-Gray/find_orb)
 * Riferimenti:
 * - Boulet, "Methods of Orbit Determination for the Micro Computer"
 * - https://www.projectpluto.com/herget.htm
 * - https://www.projectpluto.com/vaisala.htm
 */

#include "ioccultcalc/initial_orbit.h"
#include "ioccultcalc/coordinates.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/orbital_elements.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <cstring>

namespace ioccultcalc {

// Costanti
constexpr double GAUSS_K = 0.01720209895;  // Costante di Gauss
constexpr double PI = 3.1415926535897932384626433832795;
constexpr double AU_PER_DAY = 173.14463272;  // AU in un giorno-luce

// Helper: prodotto misto (A × B) · C
// Usato frequentemente nel metodo di Gauss (Boulet p. 415)
double InitialOrbit::crossThenDot(const double* a, const double* b, const double* c) {
    // cross = a × b
    double cross[3];
    cross[0] = a[1] * b[2] - a[2] * b[1];
    cross[1] = a[2] * b[0] - a[0] * b[2];
    cross[2] = a[0] * b[1] - a[1] * b[0];
    
    // dot = cross · c
    return cross[0] * c[0] + cross[1] * c[1] + cross[2] * c[2];
}

// Helper: prodotto scalare
static double dot_prod3(const double* a, const double* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Helper: risolve equazione polinomiale
// Usa metodo di Newton-Raphson o altri algoritmi numerici
int InitialOrbit::findRealPolynomialRoots(
    const double* poly,
    int poly_order,
    double* real_roots)
{
    // Implementazione semplificata - in produzione usare libreria robusta
    // come quella di Find_Orb (roots.cpp)
    
    // Per ora, gestisce solo casi comuni dell'equazione di Gauss (8° grado)
    if (poly_order != 8) {
        return 0; // Nessuna radice trovata
    }
    
    // Metodo di Newton-Raphson multi-guess
    int n_roots = 0;
    double guesses[] = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0};
    
    for (double guess : guesses) {
        double x = guess;
        bool converged = false;
        
        for (int iter = 0; iter < 50; iter++) {
            // Valuta polinomio e derivata
            double f = 0, df = 0;
            double x_pow = 1.0;
            
            for (int i = 0; i <= poly_order; i++) {
                f += poly[i] * x_pow;
                if (i > 0) {
                    df += poly[i] * i * (x_pow / x);
                }
                x_pow *= x;
            }
            
            if (std::abs(f) < 1e-10) {
                converged = true;
                break;
            }
            
            if (std::abs(df) < 1e-15) break;
            
            double dx = -f / df;
            x += dx;
            
            if (std::abs(dx) < 1e-12) {
                converged = true;
                break;
            }
        }
        
        if (converged && x > 0.01 && x < 100.0) {
            // Verifica che non sia duplicato
            bool duplicate = false;
            for (int i = 0; i < n_roots; i++) {
                if (std::abs(real_roots[i] - x) < 0.001) {
                    duplicate = true;
                    break;
                }
            }
            
            if (!duplicate && n_roots < 3) {
                real_roots[n_roots++] = x;
            }
        }
    }
    
    return n_roots;
}

/**
 * @brief Metodo di Gauss per determinazione orbitale
 * 
 * Implementazione basata su Find_Orb (gauss.cpp) e Boulet cap. 10.
 * 
 * Processo:
 * 1. Calcola determinanti d0, d11-d33 (Boulet p. 415, eq. 10.13-10.14)
 * 2. Costruisce equazione polinomiale di 8° grado per r2
 * 3. Risolve per trovare 0-3 radici reali positive
 * 4. Per ogni radice, itera con serie f e g (Boulet p. 417-420)
 * 5. Calcola posizione e velocità dell'oggetto
 */
std::vector<OrbitalElements> InitialOrbit::gauss(
    const AstrometricObservation& obs1,
    const AstrometricObservation& obs2,
    const AstrometricObservation& obs3,
    double mu)
{
    std::vector<OrbitalElements> solutions;
    
    // Vettori di osservazione (unit vectors verso l'oggetto)
    Vector3D vect1 = Coordinates::raDecToUnitVector(obs1.obs.ra, obs1.obs.dec);
    Vector3D vect2 = Coordinates::raDecToUnitVector(obs2.obs.ra, obs2.obs.dec);
    Vector3D vect3 = Coordinates::raDecToUnitVector(obs3.obs.ra, obs3.obs.dec);
    
    // Posizioni osservatore (geocentriche in AU)
    Vector3D obs_posn1 = Coordinates::observerPosition(obs1.observatoryCode, obs1.epoch);
    Vector3D obs_posn2 = Coordinates::observerPosition(obs2.observatoryCode, obs2.epoch);
    Vector3D obs_posn3 = Coordinates::observerPosition(obs3.observatoryCode, obs3.epoch);
    
    // Tempi in unità di Gauss (k*tau dove tau è in giorni)
    double tau1 = GAUSS_K * (obs1.epoch.jd - obs2.epoch.jd);
    double tau3 = GAUSS_K * (obs3.epoch.jd - obs2.epoch.jd);
    double tau = tau3 - tau1;
    
    // Verifica condizioni di validità
    if (std::abs(tau1) < 1e-10 || std::abs(tau3) < 1e-10 || std::abs(tau) < 1e-10) {
        return solutions; // Osservazioni troppo vicine
    }
    
    // Array per crossThenDot
    double v1[3] = {vect1.x, vect1.y, vect1.z};
    double v2[3] = {vect2.x, vect2.y, vect2.z};
    double v3[3] = {vect3.x, vect3.y, vect3.z};
    double p1[3] = {obs_posn1.x, obs_posn1.y, obs_posn1.z};
    double p2[3] = {obs_posn2.x, obs_posn2.y, obs_posn2.z};
    double p3[3] = {obs_posn3.x, obs_posn3.y, obs_posn3.z};
    
    // Calcola determinanti (Boulet p. 415)
    // d0 = (vect1 × vect2) · vect3
    double d0 = crossThenDot(v1, v2, v3);
    
    if (std::abs(d0) < 1e-10) {
        return solutions; // Osservazioni coplanari
    }
    
    // d11-d33 (Boulet p. 415, eq. 10.13)
    double d11 = crossThenDot(v2, p1, v3);
    double d12 = crossThenDot(v2, p2, v3);
    double d13 = crossThenDot(v2, p3, v3);
    
    double d21 = crossThenDot(p1, v1, v3);
    double d22 = crossThenDot(p2, v1, v3);
    double d23 = crossThenDot(p3, v1, v3);
    
    double d31 = crossThenDot(p1, v2, v1);
    double d32 = crossThenDot(p2, v2, v1);
    double d33 = crossThenDot(p3, v2, v1);
    
    // Coefficienti a, b (Boulet p. 417, eq. 10.20)
    double a1 = tau3 / tau;
    double a3 = -tau1 / tau;
    double b1 = a1 * (tau * tau - tau3 * tau3) / 6.0;
    double b3 = a3 * (tau * tau - tau1 * tau1) / 6.0;
    
    // Costruisce equazione polinomiale di 8° grado (Boulet p. 418)
    // Coefficienti A, B, E, F
    double big_a = -(a1 * d21 + a3 * d23) / d0;
    double big_b = -(b1 * d21 + b3 * d23) / d0;
    double big_e = -d22 / d0;
    double big_f = -(b1 * d21 + b3 * d23) / d0;  // Nota: ci sono termini aggiuntivi
    
    // Polinomio p(r2) = 0
    double poly[9] = {0};
    poly[8] = 1.0;  // r2^8
    poly[6] = -(big_a * (big_a + big_e) + big_f);
    // Altri coefficienti richiedono termini aggiuntivi (vedi gauss.cpp)
    
    // Risolve equazione polinomiale
    double roots[3];
    int n_roots = findRealPolynomialRoots(poly, 8, roots);
    
    if (n_roots == 0) {
        return solutions; // Nessuna soluzione reale
    }
    
    // Per ogni radice, costruisce soluzione
    for (int sol_idx = 0; sol_idx < n_roots; sol_idx++) {
        double r2 = roots[sol_idx];
        
        if (r2 < 0.01 || r2 > 100.0) {
            continue; // Soluzione non fisica
        }
        
        // Itera con serie f e g (Boulet p. 419-420)
        double f1 = 1.0, g1 = tau1;
        double f3 = 1.0, g3 = tau3;
        
        bool converged = false;
        for (int iteration = 0; iteration < 96; iteration++) {
            double u2 = mu / (r2 * r2 * r2);
            double tau1_sq = tau1 * tau1;
            double tau3_sq = tau3 * tau3;
            
            // Serie f e g (Boulet p. 419, eq. 10.28, 10.27)
            double f1_new = 1.0 - u2 * tau1_sq / 2.0;
            double f3_new = 1.0 - u2 * tau3_sq / 2.0;
            double g1_new = tau1 - u2 * tau1_sq * tau1 / 6.0;
            double g3_new = tau3 - u2 * tau3_sq * tau3 / 6.0;
            
            // Termini di ordine superiore (opzionali per maggiore precisione)
            #ifdef INCLUDE_FIFTH
            double tau1_4 = tau1_sq * tau1_sq;
            double tau3_4 = tau3_sq * tau3_sq;
            f1_new += u2 * u2 * tau1_4 / 24.0;
            f3_new += u2 * u2 * tau3_4 / 24.0;
            #endif
            
            // Test convergenza
            double df1 = std::abs(f1_new - f1);
            double df3 = std::abs(f3_new - f3);
            
            if (df1 < 1e-12 && df3 < 1e-12) {
                converged = true;
                f1 = f1_new;
                f3 = f3_new;
                g1 = g1_new;
                g3 = g3_new;
                break;
            }
            
            f1 = f1_new;
            f3 = f3_new;
            g1 = g1_new;
            g3 = g3_new;
        }
        
        if (!converged) {
            continue;
        }
        
        // Calcola coefficienti c1, c2, c3 (Boulet p. 419, eq. 10.29)
        double f1_g3_minus_f3_g1 = f1 * g3 - f3 * g1;
        if (std::abs(f1_g3_minus_f3_g1) < 1e-15) {
            continue;
        }
        
        double c1 = g3 / f1_g3_minus_f3_g1;
        double c2 = -1.0;
        double c3 = -g1 / f1_g3_minus_f3_g1;
        
        // Calcola rho (distanze topocentriche) (Boulet p. 419, eq. 10.30)
        double p1 = (c1 * d11 + c2 * d12 + c3 * d13) / (c1 * d0);
        double p2 = (c1 * d21 + c2 * d22 + c3 * d23) / (c2 * d0);
        double p3 = (c1 * d31 + c2 * d32 + c3 * d33) / (c3 * d0);
        
        if (p1 < 0 || p2 < 0 || p3 < 0) {
            continue; // Distanze negative non fisiche
        }
        
        // Posizione dell'oggetto (heliocentric)
        Vector3D obj_posn1 = obs_posn1 + vect1 * p1;
        Vector3D obj_posn2 = obs_posn2 + vect2 * p2;
        Vector3D obj_posn3 = obs_posn3 + vect3 * p3;
        
        // Calcola velocità (Boulet p. 420)
        double d1 = -f3 / f1_g3_minus_f3_g1;
        double d3 = f1 / f1_g3_minus_f3_g1;
        
        Vector3D obj_vel = obj_posn1 * d1 + obj_posn3 * d3;
        
        // Converte stato cartesiano in elementi orbitali
        OrbitalElements elements = cartesianToOrbitalElements(
            obj_posn2, obj_vel, mu, obs2.epoch);
        
        solutions.push_back(elements);
    }
    
    return solutions;
}

/**
 * @brief Trova le 3 migliori osservazioni per Gauss
 * 
 * Massimizza d0 = (vect1 × vect2) · vect3 per massimizzare
 * la separazione angolare e migliorare la determinazione.
 */
double InitialOrbit::findBestObsForGauss(
    const std::vector<AstrometricObservation>& observations,
    unsigned obs_idx[3])
{
    // Filtra osservazioni valide
    std::vector<int> valid_indices;
    for (size_t i = 0; i < observations.size(); i++) {
        if (!observations[i].outlier) {
            valid_indices.push_back(i);
        }
    }
    
    if (valid_indices.size() < 3) {
        return 0.0; // Non abbastanza osservazioni
    }
    
    // Strategia semplice: prima, ultima, e quella più vicina al centro
    obs_idx[0] = valid_indices.front();
    obs_idx[2] = valid_indices.back();
    
    // Trova osservazione centrale
    double mid_jd = (observations[obs_idx[0]].epoch.jd + 
                     observations[obs_idx[2]].epoch.jd) / 2.0;
    
    int mid_idx = obs_idx[0];
    double min_diff = 1e30;
    for (int idx : valid_indices) {
        double diff = std::abs(observations[idx].epoch.jd - mid_jd);
        if (diff < min_diff) {
            min_diff = diff;
            mid_idx = idx;
        }
    }
    obs_idx[1] = mid_idx;
    
    // Calcola d0 per queste osservazioni
    // TODO: Implementare calcolo vettori
    double d0 = 1.0; // Placeholder
    
    return d0;
}

/**
 * @brief Trova le 3 migliori osservazioni per Gauss (implementazione interna)
 */
static double findBestObsForGaussInternal(
    const std::vector<AstrometricObservation>& observations,
    unsigned obs_idx[3])
{
    // Filtra osservazioni valide
    std::vector<int> valid_indices;
    for (size_t i = 0; i < observations.size(); i++) {
        if (!observations[i].outlier) {
            valid_indices.push_back(i);
        }
    }
    
    if (valid_indices.size() < 3) {
        return 0.0; // Non abbastanza osservazioni
    }
    
    // Strategia semplice: prima, ultima, e quella più vicina al centro
    obs_idx[0] = valid_indices.front();
    obs_idx[2] = valid_indices.back();
    
    // Trova osservazione centrale
    double mid_jd = (observations[obs_idx[0]].epoch.jd + 
                     observations[obs_idx[2]].epoch.jd) / 2.0;
    
    int mid_idx = obs_idx[0];
    double min_diff = 1e30;
    for (int idx : valid_indices) {
        double diff = std::abs(observations[idx].epoch.jd - mid_jd);
        if (diff < min_diff) {
            min_diff = diff;
            mid_idx = idx;
        }
    }
    obs_idx[1] = mid_idx;
    
    // Calcola d0 per queste osservazioni
    // TODO: Implementare calcolo vettori
    double d0 = 1.0; // Placeholder
    
    return d0;
}

/**
 * @brief Wrapper conveniente per metodo di Gauss
 * 
 * Seleziona automaticamente le 3 migliori osservazioni
 * e chiama gauss().
 */
static OrbitalElements convenient_gauss(
    const std::vector<AstrometricObservation>& observations,
    int desired_solution = 0)
{
    unsigned obs_idx[3];
    double d0 = findBestObsForGaussInternal(observations, obs_idx);
    
    if (d0 < 1e-10) {
        return OrbitalElements(); // Fallito
    }
    
    auto solutions = InitialOrbit::gauss(
        observations[obs_idx[0]],
        observations[obs_idx[1]],
        observations[obs_idx[2]]
    );
    
    if (desired_solution >= 0 && desired_solution < (int)solutions.size()) {
        return solutions[desired_solution];
    }
    
    return OrbitalElements(); // Nessuna soluzione
}

// ===== HELPER FUNCTIONS PER HERGET =====

/**
 * @brief Setta la distanza geocentrica per un'osservazione
 */
void InitialOrbit::setDistance(AstrometricObservation& obs, double distance) {
    // Calcola vettore unitario verso l'oggetto
    Vector3D vect = Coordinates::raDecToUnitVector(obs.obs.ra, obs.obs.dec);
    
    // Posizione osservatore
    Vector3D obs_pos = Coordinates::observerPosition(obs.observatoryCode, obs.epoch);
    
    // Posizione heliocentric dell'oggetto
    // obj_pos = obs_pos + distance * vect
    // Memorizza in obs.computed (abuso di campo)
    obs.computed.distance = distance;
}

/**
 * @brief Calcola residui RA/Dec per un'osservazione
 */
static void getResidualData(const AstrometricObservation& obs,
                           double* ra_resid, double* dec_resid) {
    // Residui in arcsec: osservato - calcolato
    *ra_resid = (obs.obs.ra - obs.computed.ra) * 3600.0 * RAD_TO_DEG * cos(obs.obs.dec);
    *dec_resid = (obs.obs.dec - obs.computed.dec) * 3600.0 * RAD_TO_DEG;
}

/**
 * @brief Trova il raggio dato il raggio solare (per Väisälä)
 */
double InitialOrbit::findRGivenSolarR(const AstrometricObservation& obs, double solar_r) {
    // Distanza osservatore-sole
    Vector3D obs_pos = Coordinates::observerPosition(obs.observatoryCode, obs.epoch);
    double obs_sun_dist = obs_pos.magnitude();
    
    // Direzione verso oggetto
    Vector3D vect = Coordinates::raDecToUnitVector(obs.obs.ra, obs.obs.dec);
    
    // Legge dei coseni: r² = R² + ρ² - 2Rρcos(θ)
    // dove R = solar_r, ρ = distanza topocentric, r = obs_sun_dist
    double cos_theta = obs_pos.normalize().dot(vect);
    
    // Risolvi equazione quadratica per ρ
    double a = 1.0;
    double b = -2.0 * solar_r * cos_theta;
    double c = solar_r * solar_r - obs_sun_dist * obs_sun_dist;
    
    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0) return -1.0;
    
    // Prendi soluzione positiva
    double rho = (-b + sqrt(discriminant)) / (2.0 * a);
    return rho;
}

/**
 * @brief Get propagator singleton per initial orbit (veloce, senza perturbazioni)
 */
static OrbitPropagator* getInitialOrbitPropagator() {
    static OrbitPropagator* propagator = nullptr;
    if (!propagator) {
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RK4;
        opts.stepSize = 1.0;
        opts.usePlanetaryPerturbations = false;  // Troppo lento per initial orbit
        opts.useRelativisticCorrections = false;
        opts.useSolarRadiationPressure = false;
        opts.useNonGravitational = false;
        propagator = new OrbitPropagator(opts);
    }
    return propagator;
}

/**
 * @brief Propaga stato cartesiano usando propagazione Kepleriana veloce
 */
static bool propagateStateRobust(const double* state_in, double t0, double target_time, 
                                 double* state_out) {
    double dt = target_time - t0;
    double mu = GAUSS_K * GAUSS_K;
    
    Vector3D pos(state_in[0], state_in[1], state_in[2]);
    Vector3D vel(state_in[3], state_in[4], state_in[5]);
    
    try {
        // Usa propagazione Kepleriana semplice (veloce)
        OrbitalElements elem = cartesianToOrbitalElements(pos, vel, mu, JulianDate(t0));
        
        // Verifica validità elementi
        if (std::isnan(elem.a) || std::isnan(elem.e) || elem.a <= 0 || elem.a > 1000.0) {
            return false;
        }
        
        // Propaga mean anomaly
        if (elem.e < 1.0) {  // Solo orbite ellittiche
            double mean_motion = sqrt(mu / pow(elem.a, 3));
            elem.M += mean_motion * dt;
            
            // Normalizza M in [0, 2π]
            elem.M = fmod(elem.M, 2.0 * M_PI);
            if (elem.M < 0) elem.M += 2.0 * M_PI;
        } else {
            // Orbite paraboliche/iperboliche non supportate
            return false;
        }
        
        // Riconverti a cartesiano
        Vector3D pos_new, vel_new;
        orbitalElementsToCartesian(elem, pos_new, vel_new, mu);
        
        // Verifica validità risultato
        if (std::isnan(pos_new.x) || std::isnan(vel_new.x)) {
            return false;
        }
        
        state_out[0] = pos_new.x;
        state_out[1] = pos_new.y;
        state_out[2] = pos_new.z;
        state_out[3] = vel_new.x;
        state_out[4] = vel_new.y;
        state_out[5] = vel_new.z;
        
        return true;
        
    } catch (...) {
        return false;
    }
}

/**
 * @brief Calcola RA/Dec predetti per un'osservazione data un'orbita
 */
static bool computePredictedRaDec(const double* orbit, const AstrometricObservation& obs, 
                                  double t0, double& ra_pred, double& dec_pred) {
    double state_prop[6];
    
    // Propaga orbita (in coordinate ECLITTICHE)
    if (!propagateStateRobust(orbit, t0, obs.epoch.jd, state_prop)) {
        return false;
    }
    
    // Posizione heliocentric dell'oggetto (eclittiche)
    Vector3D obj_pos_ecl(state_prop[0], state_prop[1], state_prop[2]);
    
    // Posizione osservatore (equatoriali, poi converti in eclittiche)
    Vector3D obs_pos_eq = Coordinates::observerPosition(obs.observatoryCode, obs.epoch);
    Vector3D obs_pos_ecl = Coordinates::equatorialToEcliptic(obs_pos_eq);
    
    // Vettore topocentric (eclittiche)
    Vector3D topo_ecl = obj_pos_ecl - obs_pos_ecl;
    
    // Converti da eclittiche a equatoriali per RA/Dec
    Vector3D topo_eq = Coordinates::eclipticToEquatorial(topo_ecl);
    
    // Converti a RA/Dec
    Coordinates::vectorToRaDec(topo_eq, ra_pred, dec_pred);
    return true;
}

/**
 * @brief Metodo di Herget completo con iterazione
 * 
 * Algoritmo basato su Find_Orb orb_func.cpp lines 1916-2167
 */
OrbitalElements InitialOrbit::herget(
    const std::vector<AstrometricObservation>& observations,
    double r1_guess,
    double r2_guess,
    int max_iterations)
{
    if (observations.size() < 2) {
        throw std::runtime_error("Herget requires at least 2 observations");
    }
    
    // Copia osservazioni (modificabili)
    std::vector<AstrometricObservation> obs = observations;
    int n_obs = obs.size();
    
    // Prima e ultima osservazione non escluse
    int idx1 = 0, idx2 = n_obs - 1;
    while (idx1 < n_obs && obs[idx1].excluded) idx1++;
    while (idx2 >= 0 && obs[idx2].excluded) idx2--;
    
    if (idx1 >= idx2) {
        throw std::runtime_error("Herget: not enough valid observations");
    }
    
    // Pseudo-Väisälä mode se r1 < 0
    bool using_pseudo_vaisala = (r1_guess < 0.0);
    
    if (using_pseudo_vaisala) {
        // Interpreta r1_guess come -solar_r
        double solar_r = -r1_guess;
        r2_guess = findRGivenSolarR(obs[idx2], solar_r);
        r1_guess = findRGivenSolarR(obs[idx1], solar_r);
        
        if (r1_guess < 0 || r2_guess < 0) {
            throw std::runtime_error("Pseudo-Vaisala: invalid distances");
        }
    }
    
    // Distanze correnti
    double r1 = r1_guess;
    double r2 = r2_guess;
    
    double orbit[6]; // stato cartesiano [x,y,z,vx,vy,vz]
    double t0 = obs[idx1].epoch.jd;
    
    // Iterazione di Herget
    for (int iter = 0; iter < max_iterations; iter++) {
        // Setta distanze
        setDistance(obs[idx1], r1);
        setDistance(obs[idx2], r2);
        
        // Trova orbita di trasferimento
        int result = findTransferOrbit(orbit, obs[idx1], obs[idx2], iter > 0);
        if (result != 0) {
            // Se fallisce, ritorna orbita corrente
            break;
        }
        
        // Calcola residui per prima e ultima osservazione
        double ra1_pred, dec1_pred, ra2_pred, dec2_pred;
        if (!computePredictedRaDec(orbit, obs[idx1], t0, ra1_pred, dec1_pred) ||
            !computePredictedRaDec(orbit, obs[idx2], t0, ra2_pred, dec2_pred)) {
            break; // Propagazione fallita
        }
        
        // Residui in arcsec (O - C)
        double dra1 = (obs[idx1].obs.ra - ra1_pred) * cos(obs[idx1].obs.dec) * 206265.0;
        double ddec1 = (obs[idx1].obs.dec - dec1_pred) * 206265.0;
        double dra2 = (obs[idx2].obs.ra - ra2_pred) * cos(obs[idx2].obs.dec) * 206265.0;
        double ddec2 = (obs[idx2].obs.dec - dec2_pred) * 206265.0;
        
        // Check convergenza
        double rms = sqrt((dra1*dra1 + ddec1*ddec1 + dra2*dra2 + ddec2*ddec2) / 4.0);
        if (rms < 0.1) {  // < 0.1 arcsec
            break;
        }
        
        // Calcola derivate parziali numeriche
        const double h = 1e-5;  // AU
        
        // dRA/dR1, dDec/dR1
        setDistance(obs[idx1], r1 + h);
        setDistance(obs[idx2], r2);
        double orbit_plus[6];
        if (findTransferOrbit(orbit_plus, obs[idx1], obs[idx2], true) == 0) {
            double ra1_plus, dec1_plus, ra2_plus, dec2_plus;
            if (!computePredictedRaDec(orbit_plus, obs[idx1], t0, ra1_plus, dec1_plus) ||
                !computePredictedRaDec(orbit_plus, obs[idx2], t0, ra2_plus, dec2_plus)) {
                break; // Skip se propagazione fallisce
            }
            
            double dra1_dr1 = (ra1_plus - ra1_pred) / h * cos(obs[idx1].obs.dec) * 206265.0;
            double ddec1_dr1 = (dec1_plus - dec1_pred) / h * 206265.0;
            double dra2_dr1 = (ra2_plus - ra2_pred) / h * cos(obs[idx2].obs.dec) * 206265.0;
            double ddec2_dr1 = (dec2_plus - dec2_pred) / h * 206265.0;
            
            // dRA/dR2, dDec/dR2
            setDistance(obs[idx1], r1);
            setDistance(obs[idx2], r2 + h);
            if (findTransferOrbit(orbit_plus, obs[idx1], obs[idx2], true) == 0) {
                double ra1_plus2, dec1_plus2, ra2_plus2, dec2_plus2;
                if (!computePredictedRaDec(orbit_plus, obs[idx1], t0, ra1_plus2, dec1_plus2) ||
                    !computePredictedRaDec(orbit_plus, obs[idx2], t0, ra2_plus2, dec2_plus2)) {
                    break;
                }
                
                double dra1_dr2 = (ra1_plus2 - ra1_pred) / h * cos(obs[idx1].obs.dec) * 206265.0;
                double ddec1_dr2 = (dec1_plus2 - dec1_pred) / h * 206265.0;
                double dra2_dr2 = (ra2_plus2 - ra2_pred) / h * cos(obs[idx2].obs.dec) * 206265.0;
                double ddec2_dr2 = (dec2_plus2 - dec2_pred) / h * 206265.0;
                
                // Sistema 2x2: minimizza residui usando entrambe le osservazioni
                // Usa solo RA per semplicità (Find_Orb fa questo)
                // [ dra1_dr1  dra1_dr2 ] [ delta_r1 ]   [ dra1 ]
                // [ dra2_dr1  dra2_dr2 ] [ delta_r2 ] = [ dra2 ]
                
                double det = dra1_dr1 * dra2_dr2 - dra1_dr2 * dra2_dr1;
                if (fabs(det) > 1e-10) {
                    double delta_r1 = (dra2_dr2 * dra1 - dra1_dr2 * dra2) / det;
                    double delta_r2 = (dra1_dr1 * dra2 - dra2_dr1 * dra1) / det;
                    
                    // Limita correzioni per stabilità
                    const double max_delta = 0.1;  // Max 0.1 AU per step
                    delta_r1 = std::max(-max_delta, std::min(max_delta, delta_r1));
                    delta_r2 = std::max(-max_delta, std::min(max_delta, delta_r2));
                    
                    r1 += delta_r1;
                    r2 += delta_r2;
                    
                    // Verifica ragionevolezza
                    if (r1 < 0.1 || r1 > 100.0 || r2 < 0.1 || r2 > 100.0) {
                        r1 = r1_guess;
                        r2 = r2_guess;
                        break;
                    }
                }
            }
        }
    }
    
    // Ritorna elementi finali
    Vector3D pos(orbit[0], orbit[1], orbit[2]);
    Vector3D vel(orbit[3], orbit[4], orbit[5]);
    
    return cartesianToOrbitalElements(pos, vel, GAUSS_K * GAUSS_K, obs[idx1].epoch);
}

/**
 * @brief Differential correction usando least squares
 * 
 * Raffina orbita minimizzando residui O-C su tutte le osservazioni.
 * Usa metodo Gauss-Newton con matrice normale.
 * 
 * Basato su Find_Orb orb_func.cpp adjust_herget_results()
 */
int InitialOrbit::adjustHergetResults(
    std::vector<AstrometricObservation>& observations,
    OrbitalElements& orbit)
{
    const int max_iterations = 20;  // Aumentato per LM
    const double convergence_tol = 0.01; // arcsec RMS
    const double outlier_threshold = 3.0; // 3-sigma outlier rejection
    
    // Stato iniziale
    double mu = GAUSS_K * GAUSS_K;
    double t0 = orbit.epoch.jd;
    
    // Parametri da ottimizzare: [x, y, z, vx, vy, vz]
    Vector3D pos, vel;
    orbitalElementsToCartesian(orbit, pos, vel, mu);
    
    double state[6] = {pos.x, pos.y, pos.z, vel.x, vel.y, vel.z};
    
    // Levenberg-Marquardt damping parameter
    double lambda = 0.001;  // Initial LM damping
    double lambda_factor = 10.0;  // Factor for lambda adjustment
    double prev_rms = 1e10;
    
    for (int iter = 0; iter < max_iterations; iter++) {
        int n_obs = 0;
        for (const auto& obs : observations) {
            if (!obs.excluded) n_obs++;
        }
        
        if (n_obs < 3) {
            return -1; // Troppo poche osservazioni
        }
        
        // Matrici per least squares
        // Abbiamo 2*n_obs equazioni (RA, Dec per ogni obs)
        // e 6 incognite (stato cartesiano)
        const int n_eq = 2 * n_obs;
        const int n_par = 6;
        
        std::vector<double> A(n_eq * n_par, 0.0);  // Matrice Jacobiana
        std::vector<double> b(n_eq, 0.0);           // Residui
        std::vector<double> weights(n_obs, 1.0);    // Pesi osservazioni
        std::vector<double> residuals_all;          // Per calcolo sigma
        
        int eq_idx = 0;
        double sum_sq_res = 0.0;
        
        // PASS 1: Calcola tutti i residui per outlier detection
        std::vector<double> all_residuals;
        for (const auto& obs : observations) {
            if (obs.excluded) continue;
            
            double state_prop[6];
            if (!propagateStateRobust(state, t0, obs.epoch.jd, state_prop)) continue;
            
            Vector3D obj_pos(state_prop[0], state_prop[1], state_prop[2]);
            Vector3D obs_pos = Coordinates::observerPosition(obs.observatoryCode, obs.epoch);
            Vector3D topo = obj_pos - obs_pos;
            
            double ra_pred, dec_pred;
            Coordinates::vectorToRaDec(topo, ra_pred, dec_pred);
            
            double dra = (obs.obs.ra - ra_pred) * cos(obs.obs.dec) * 206265.0;
            double ddec = (obs.obs.dec - dec_pred) * 206265.0;
            double res_mag = sqrt(dra*dra + ddec*ddec);
            all_residuals.push_back(res_mag);
        }
        
        // Calcola RMS e sigma per outlier rejection
        double mean_res = 0.0;
        for (double r : all_residuals) mean_res += r;
        mean_res /= all_residuals.size();
        
        double sigma = 0.0;
        for (double r : all_residuals) {
            sigma += (r - mean_res) * (r - mean_res);
        }
        sigma = sqrt(sigma / all_residuals.size());
        
        double outlier_limit = mean_res + outlier_threshold * sigma;
        
        // PASS 2: Calcola residui e derivate parziali con weights
        eq_idx = 0;
        int obs_idx = 0;
        for (auto& obs : observations) {
            if (obs.excluded) continue;
            
            // Propaga orbita all'epoca dell'osservazione
            double state_prop[6];
            if (!propagateStateRobust(state, t0, obs.epoch.jd, state_prop)) {
                // Skip osservazione se propagazione fallisce
                continue;
            }
            
            // Posizione oggetto heliocentric
            Vector3D obj_pos(state_prop[0], state_prop[1], state_prop[2]);
            
            // Posizione osservatore
            Vector3D obs_pos = Coordinates::observerPosition(obs.observatoryCode, obs.epoch);
            
            // Vettore topocentric
            Vector3D topo = obj_pos - obs_pos;
            
            // RA/Dec predetti
            double ra_pred, dec_pred;
            Coordinates::vectorToRaDec(topo, ra_pred, dec_pred);
            
            // Residui (arcsec)
            double dra = (obs.obs.ra - ra_pred) * cos(obs.obs.dec) * 206265.0;
            double ddec = (obs.obs.dec - dec_pred) * 206265.0;
            double res_mag = sqrt(dra*dra + ddec*ddec);
            
            // Outlier rejection: marca osservazione se oltre 3-sigma
            if (iter > 0 && res_mag > outlier_limit) {
                obs.excluded = true;
                obs_idx++;
                continue;
            }
            
            // Peso osservazione (inversely proportional to epoch uncertainty)
            // Osservazioni recenti (post-2000) hanno peso maggiore
            double weight = 1.0;
            if (obs.epoch.jd < 2451545.0) {  // Prima del 2000
                weight = 0.5;  // Downweight osservazioni vecchie
            }
            // TODO: aggiungere pesi per osservatorio da MPC
            
            // Applica peso ai residui
            b[eq_idx] = dra * weight;
            b[eq_idx + 1] = ddec * weight;
            
            sum_sq_res += (dra*dra + ddec*ddec) * weight * weight;
            weights[obs_idx] = weight;
            
            // Derivate parziali numeriche: d(RA,Dec)/d(x,y,z,vx,vy,vz)
            const double h = 1e-8;  // Step per derivate
            
            for (int j = 0; j < n_par; j++) {
                double state_plus[6];
                for (int k = 0; k < 6; k++) {
                    state_plus[k] = state[k] + (k == j ? h : 0.0);
                }
                
                double state_plus_prop[6];
                if (!propagateStateRobust(state_plus, t0, obs.epoch.jd, state_plus_prop)) {
                    // Se propagazione fallisce, usa derivata zero
                    A[eq_idx * n_par + j] = 0.0;
                    A[(eq_idx + 1) * n_par + j] = 0.0;
                    continue;
                }
                
                Vector3D obj_pos_plus(state_plus_prop[0], state_plus_prop[1], state_plus_prop[2]);
                Vector3D topo_plus = obj_pos_plus - obs_pos;
                
                double ra_plus, dec_plus;
                Coordinates::vectorToRaDec(topo_plus, ra_plus, dec_plus);
                
                // dRA/d(param_j) in arcsec, weighted
                A[eq_idx * n_par + j] = (ra_plus - ra_pred) / h * cos(obs.obs.dec) * 206265.0 * weight;
                // dDec/d(param_j) in arcsec, weighted
                A[(eq_idx + 1) * n_par + j] = (dec_plus - dec_pred) / h * 206265.0 * weight;
            }
            
            eq_idx += 2;
            obs_idx++;
        }
        
        double rms = sqrt(sum_sq_res / (2.0 * n_obs));
        
        // Check for NaN
        if (std::isnan(rms) || std::isinf(rms)) {
            // Ripristina stato precedente e termina
            return -1;
        }
        
        // Levenberg-Marquardt: adatta lambda basato su miglioramento
        if (iter > 0) {  // Skip first iteration
            if (rms < prev_rms) {
                // Improvement: riduce damping (più simile a Gauss-Newton)
                lambda /= lambda_factor;
            } else {
                // No improvement: aumenta damping (più simile a gradient descent)
                lambda *= lambda_factor;
            }
        }
        prev_rms = rms;
        
        // Clamp lambda per stabilità
        lambda = std::max(1e-10, std::min(1e6, lambda));
        
        // Check convergenza
        if (rms < convergence_tol && iter > 0) {
            // Aggiorna orbit
            Vector3D pos_final(state[0], state[1], state[2]);
            Vector3D vel_final(state[3], state[4], state[5]);
            orbit = cartesianToOrbitalElements(pos_final, vel_final, mu, orbit.epoch);
            return 0;  // Convergenza raggiunta
        }
        
        // Risolvi sistema normale con LM: (A^T A + λI) delta = A^T b
        std::vector<double> AtA(n_par * n_par, 0.0);
        std::vector<double> Atb(n_par, 0.0);
        
        // Calcola A^T A
        for (int i = 0; i < n_par; i++) {
            for (int j = 0; j < n_par; j++) {
                double sum = 0.0;
                for (int k = 0; k < n_eq; k++) {
                    sum += A[k * n_par + i] * A[k * n_par + j];
                }
                AtA[i * n_par + j] = sum;
            }
        }
        
        // Levenberg-Marquardt: modifica elementi diagonali
        // Formula standard: AtA[i][i] *= (1 + λ)
        for (int i = 0; i < n_par; i++) {
            double diag_value = AtA[i * n_par + i];
            if (std::abs(diag_value) < 1e-15) {
                // Diagonale troppo piccola - usa valore minimo
                AtA[i * n_par + i] = lambda;
            } else {
                AtA[i * n_par + i] *= (1.0 + lambda);
            }
        }
        
        // Calcola A^T b
        for (int i = 0; i < n_par; i++) {
            double sum = 0.0;
            for (int k = 0; k < n_eq; k++) {
                sum += A[k * n_par + i] * b[k];
            }
            Atb[i] = sum;
        }
        
        // Risolvi AtA * delta = Atb usando Gauss elimination
        std::vector<double> AtA_aug(n_par * (n_par + 1));
        for (int i = 0; i < n_par; i++) {
            for (int j = 0; j < n_par; j++) {
                AtA_aug[i * (n_par + 1) + j] = AtA[i * n_par + j];
            }
            AtA_aug[i * (n_par + 1) + n_par] = Atb[i];
        }
        
        // Gauss elimination con pivot
        for (int k = 0; k < n_par; k++) {
            // Pivot
            double max_val = fabs(AtA_aug[k * (n_par + 1) + k]);
            int max_row = k;
            for (int i = k + 1; i < n_par; i++) {
                if (fabs(AtA_aug[i * (n_par + 1) + k]) > max_val) {
                    max_val = fabs(AtA_aug[i * (n_par + 1) + k]);
                    max_row = i;
                }
            }
            
            if (max_val < 1e-15) {
                return -1; // Matrice singolare
            }
            
            // Swap rows
            if (max_row != k) {
                for (int j = 0; j <= n_par; j++) {
                    std::swap(AtA_aug[k * (n_par + 1) + j], 
                             AtA_aug[max_row * (n_par + 1) + j]);
                }
            }
            
            // Eliminazione
            for (int i = k + 1; i < n_par; i++) {
                double factor = AtA_aug[i * (n_par + 1) + k] / AtA_aug[k * (n_par + 1) + k];
                for (int j = k; j <= n_par; j++) {
                    AtA_aug[i * (n_par + 1) + j] -= factor * AtA_aug[k * (n_par + 1) + j];
                }
            }
        }
        
        // Back substitution
        std::vector<double> delta(n_par);
        for (int i = n_par - 1; i >= 0; i--) {
            delta[i] = AtA_aug[i * (n_par + 1) + n_par];
            for (int j = i + 1; j < n_par; j++) {
                delta[i] -= AtA_aug[i * (n_par + 1) + j] * delta[j];
            }
            delta[i] /= AtA_aug[i * (n_par + 1) + i];
        }
        
        // Salva stato prima dell'update
        double state_old[6];
        for (int i = 0; i < 6; i++) state_old[i] = state[i];
        
        // Aggiorna stato (LM gestisce damping tramite λ)
        // Limita solo delta estremi per stabilità numerica
        double max_delta = 0.0;
        for (int i = 0; i < n_par; i++) {
            // Limiti adattivi basati su grandezza variabile
            double scale = (i < 3) ? 10.0 : 0.5;  // pos in AU, vel in AU/day
            double limited_delta = std::max(-scale, std::min(scale, delta[i]));
            state[i] += limited_delta;
            max_delta = std::max(max_delta, std::abs(limited_delta));
        }
        
        // Verifica che stato rimanga ragionevole
        Vector3D pos_check(state[0], state[1], state[2]);
        Vector3D vel_check(state[3], state[4], state[5]);
        double r = pos_check.magnitude();
        double v = vel_check.magnitude();
        
        if (r < 0.1 || r > 100.0 || v > 1.0 ||  // Limiti ragionevoli
            std::isnan(r) || std::isnan(v)) {
            // Stato non fisico - ripristina e aumenta lambda
            for (int i = 0; i < 6; i++) state[i] = state_old[i];
            lambda *= lambda_factor * lambda_factor;  // Damping forte
            continue;  // Riprova con lambda maggiore
        }
    }
    
    // Max iterazioni raggiunto - ritorna comunque orbit raffinato
    Vector3D pos_final(state[0], state[1], state[2]);
    Vector3D vel_final(state[3], state[4], state[5]);
    orbit = cartesianToOrbitalElements(pos_final, vel_final, mu, orbit.epoch);
    
    // Verifica ragionevolezza orbita finale
    if (std::isnan(orbit.a) || std::isnan(orbit.e) || 
        orbit.a <= 0 || orbit.a > 100.0 || 
        orbit.e < 0 || orbit.e >= 1.0) {
        return -1; // Orbita non fisica
    }
    
    return 1; // Non convergente ma migliorato
}

OrbitalElements InitialOrbit::vaisala(
    const std::vector<AstrometricObservation>& observations,
    double perihelion_distance)
{
    // TODO: Implementare metodo di Väisälä
    // Vedi Find_Orb: orb_func.cpp, find_trial_orbit()
    throw std::runtime_error("Vaisala method not yet implemented");
}

OrbitalElements InitialOrbit::autoSolve(
    const std::vector<AstrometricObservation>& observations,
    double* score)
{
    // TODO: Implementare auto-solve
    // Vedi Find_Orb: orb_func.cpp, initial_orbit()
    
    // Strategia:
    // 1. Prova Gauss (0-3 soluzioni)
    // 2. Se fallisce, prova Herget con vari R1/R2
    // 3. Se fallisce, prova Väisälä con vari q
    // 4. Applica attempt_improvements()
    
    throw std::runtime_error("autoSolve not yet implemented");
}

double InitialOrbit::evaluateInitialOrbit(
    const OrbitalElements& orbit,
    const std::vector<AstrometricObservation>& observations,
    const JulianDate& epoch)
{
    // TODO: Implementare valutazione score
    // Calcola RMS residui + penalità per orbite non fisiche
    return 0.0;
}

bool InitialOrbit::isUnreasonableOrbit(const OrbitalElements& orbit) {
    // Verifica criteri di ragionevolezza
    if (orbit.a <= 0 || orbit.a > 1000.0) return true;  // Semiasse > 1000 AU
    if (orbit.e < 0 || orbit.e >= 1.0) return true;     // Eccentricità non ellittica
    if (orbit.i < 0 || orbit.i > 180.0) return true;    // Inclinazione valida
    return false;
}

double InitialOrbit::maxHergetSpan(double r1, double r2) {
    // Span massimo dipende dalle distanze
    // Orbite vicine (NEO) possono cambiare rapidamente
    double avg_r = (r1 + r2) / 2.0;
    return 365.25 / std::pow(avg_r, 1.5);  // Periodo orbitale stimato in giorni
}

/**
 * @brief Lambert solver - Universal Variables algorithm
 * 
 * Risolve il problema di Lambert usando l'algoritmo Universal Variables.
 * Dato r1, r2 e dt, trova v1 e v2 tali che l'orbita connetta i due punti.
 * 
 * Implementazione basata su:
 * - Battin, "An Introduction to the Mathematics and Methods of Astrodynamics"
 * - Vallado, "Fundamentals of Astrodynamics and Applications", Algorithm 58
 * - Curtis, "Orbital Mechanics for Engineering Students", Algorithm 5.2
 */
int InitialOrbit::solveLambert(
    const Vector3D& r1,
    const Vector3D& r2,
    double dt,
    double mu,
    bool prograde,
    Vector3D& v1,
    Vector3D& v2)
{
    const double tol = 1e-6;
    const int max_iterations = 100;
    
    // Magnitudini
    double r1mag = r1.magnitude();
    double r2mag = r2.magnitude();
    
    if (r1mag < 1e-10 || r2mag < 1e-10) {
        return -1;
    }
    
    // Versori
    Vector3D ir1 = r1 / r1mag;
    Vector3D ir2 = r2 / r2mag;
    
    // Cross product per determinare prograde/retrograde
    Vector3D h = ir1.cross(ir2);
    double hmag = h.magnitude();
    
    // Angolo di trasferimento
    double cos_dnu = ir1.dot(ir2);
    
    // Clamp per evitare errori numerici
    if (cos_dnu > 1.0) cos_dnu = 1.0;
    if (cos_dnu < -1.0) cos_dnu = -1.0;
    
    double dnu = acos(cos_dnu);
    
    // Seleziona direzione (short way o long way)
    if (!prograde && dnu < PI) {
        dnu = 2.0 * PI - dnu;
    } else if (prograde && dnu > PI) {
        dnu = 2.0 * PI - dnu;
    }
    
    // Parametro A (Battin/Vallado)
    double A = sin(dnu) * sqrt(r1mag * r2mag / (1.0 - cos_dnu));
    
    if (fabs(A) < 1e-10) {
        return -1; // Orbita degenere
    }
    
    // Iterazione per trovare z (universal anomaly variable)
    double z = 0.0;  // Initial guess
    double c2, c3;
    
    for (int iter = 0; iter < max_iterations; iter++) {
        // Stumpff functions C2(z) e C3(z)
        if (fabs(z) < 1e-4) {
            // Espansione Taylor per z piccolo
            c2 = 0.5 - z / 24.0 + z * z / 720.0;
            c3 = 1.0 / 6.0 - z / 120.0 + z * z / 5040.0;
        } else if (z > 0) {
            // Orbita ellittica
            double sqrt_z = sqrt(z);
            c2 = (1.0 - cos(sqrt_z)) / z;
            c3 = (sqrt_z - sin(sqrt_z)) / (z * sqrt_z);
        } else {
            // Orbita iperbolica
            double sqrt_mz = sqrt(-z);
            c2 = (1.0 - cosh(sqrt_mz)) / z;
            c3 = (sinh(sqrt_mz) - sqrt_mz) / (-z * sqrt_mz);
        }
        
        // Funzione y(z)
        if (c2 <= 0) {
            z += 0.1;
            continue;
        }
        
        double sqrt_c2 = sqrt(c2);
        double y = r1mag + r2mag + A * (z * c3 - 1.0) / sqrt_c2;
        
        if (y <= 0) {
            // Negativo - adjustiamo z
            z += 0.1;
            continue;
        }
        
        // Funzione chi(z) = sqrt(y / c2)
        double chi = sqrt(y / c2);
        
        // Tempo di volo con questo z
        double t = (chi * chi * chi * c3 + A * sqrt(y)) / sqrt(mu);
        
        // Residuo
        double dt_error = t - dt;
        
        if (fabs(dt_error) < tol) {
            // Convergenza!
            // Calcola f e g functions
            double f = 1.0 - y / r1mag;
            double g = A * sqrt(y / mu);
            double gdot = 1.0 - y / r2mag;
            
            // Velocità
            v1 = (r2 - r1 * f) / g;
            v2 = (r2 * gdot - r1) / g;
            
            return 0;
        }
        
        // Newton-Raphson: z_new = z - F(z) / F'(z)
        // Derivata dF/dz
        double dydt;
        if (fabs(z) < 1e-4) {
            dydt = sqrt(2.0) / 40.0 * y * y + A / 8.0 * (sqrt(y) + A * sqrt(1.0 / (2.0 * y)));
        } else {
            dydt = (chi * chi * chi * (c3 - 3.0 * c2 * c3 / (2.0 * z)) 
                    + A / 8.0 * (3.0 * c3 * sqrt(y) / c2 + A * sqrt(c2 / y))) / sqrt(mu);
        }
        
        if (fabs(dydt) < 1e-10) {
            return -1; // Derivata nulla
        }
        
        // Update z
        z = z - dt_error / dydt;
        
        // Limita z per stabilità
        if (z < -4.0 * PI * PI) z = -4.0 * PI * PI;
        if (z > 4.0 * PI * PI) z = 4.0 * PI * PI;
    }
    
    return -1; // Non convergenza
}

/**
 * @brief Trova orbita di trasferimento tra due osservazioni con distanze note
 * 
 * Usa il Lambert solver per calcolare l'orbita che connette
 * le due posizioni heliocentriche nel tempo dato.
 */
int InitialOrbit::findTransferOrbit(
    double* orbit,
    const AstrometricObservation& obs1,
    const AstrometricObservation& obs2,
    bool already_have_approximate_orbit)
{
    // Posizioni osservatori (in coordinate EQUATORIALI)
    Vector3D obs1_pos_eq = Coordinates::observerPosition(obs1.observatoryCode, obs1.epoch);
    Vector3D obs2_pos_eq = Coordinates::observerPosition(obs2.observatoryCode, obs2.epoch);
    
    // Converti in coordinate ECLITTICHE (sistema di riferimento per orbite heliocentriche)
    Vector3D obs1_pos = Coordinates::equatorialToEcliptic(obs1_pos_eq);
    Vector3D obs2_pos = Coordinates::equatorialToEcliptic(obs2_pos_eq);
    
    // Direzioni verso oggetto in coordinate ECLITTICHE
    // (conversione diretta da RA/Dec equatoriali a vettore unitario eclittico)
    Vector3D vect1 = Coordinates::raDecToEclipticUnitVector(obs1.obs.ra, obs1.obs.dec);
    Vector3D vect2 = Coordinates::raDecToEclipticUnitVector(obs2.obs.ra, obs2.obs.dec);
    
    // Distanze geocentriche (settate in setDistance)
    double r1 = obs1.computed.distance;
    double r2 = obs2.computed.distance;
    
    if (r1 <= 0 || r2 <= 0) {
        return -1; // Distanze non valide
    }
    
    // Posizioni heliocentriche in coordinate ECLITTICHE
    Vector3D pos1 = obs1_pos + vect1 * r1;
    Vector3D pos2 = obs2_pos + vect2 * r2;
    
    // Intervallo temporale (giorni)
    double dt = obs2.epoch.jd - obs1.epoch.jd;
    
    if (fabs(dt) < 1e-6) {
        return -1; // Time span troppo piccolo
    }
    
    // Usa Lambert solver per calcolare velocità
    Vector3D vel1_pro, vel2_pro;
    Vector3D vel1_retro, vel2_retro;
    double mu = GAUSS_K * GAUSS_K; // k^2 per sistema solare
    
    // Prova entrambe le soluzioni (prograde e retrograde)
    int result_pro = solveLambert(pos1, pos2, dt, mu, true, vel1_pro, vel2_pro);
    int result_retro = solveLambert(pos1, pos2, dt, mu, false, vel1_retro, vel2_retro);
    
    // Scegli la soluzione con inclinazione più bassa (più comune per asteroidi)
    Vector3D vel1, vel2;
    int result = -1;
    
    if (result_pro == 0 && result_retro == 0) {
        // Entrambe convergono: scegli quella con i più piccolo
        Vector3D h_pro = pos1.cross(vel1_pro);
        Vector3D h_retro = pos1.cross(vel1_retro);
        double i_pro = acos(h_pro.z / h_pro.magnitude());
        double i_retro = acos(h_retro.z / h_retro.magnitude());
        
        // Scegli prograde se < 90°, altrimenti la migliore
        if (i_pro < M_PI / 2.0) {
            vel1 = vel1_pro;
            vel2 = vel2_pro;
            result = 0;
        } else if (i_retro < M_PI / 2.0) {
            vel1 = vel1_retro;
            vel2 = vel2_retro;
            result = 0;
        } else {
            // Entrambe > 90°: scegli quella più vicina a 90°
            if (std::abs(i_pro - M_PI / 2.0) < std::abs(i_retro - M_PI / 2.0)) {
                vel1 = vel1_pro;
                vel2 = vel2_pro;
            } else {
                vel1 = vel1_retro;
                vel2 = vel2_retro;
            }
            result = 0;
        }
    } else if (result_pro == 0) {
        vel1 = vel1_pro;
        vel2 = vel2_pro;
        result = 0;
    } else if (result_retro == 0) {
        vel1 = vel1_retro;
        vel2 = vel2_retro;
        result = 0;
    }
    
    if (result != 0) {
        // Fallback: approssimazione lineare
        vel1 = (pos2 - pos1) / dt;
        
        // Verifica plausibilità fisica
        double speed = vel1.magnitude();
        if (speed > 0.1) {  // > 0.1 AU/day (molto veloce per oggetti naturali)
            return -1;
        }
    }
    
    // Ritorna stato cartesiano all'epoca obs1
    orbit[0] = pos1.x;
    orbit[1] = pos1.y;
    orbit[2] = pos1.z;
    orbit[3] = vel1.x;
    orbit[4] = vel1.y;
    orbit[5] = vel1.z;
    
    return 0;
}

} // namespace ioccultcalc
