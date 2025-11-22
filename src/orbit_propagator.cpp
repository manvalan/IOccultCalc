/**
 * @file orbit_propagator.cpp
 * @brief Implementazione propagatore numerico orbite
 */

#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/jpl_ephemeris.h"
#include "ioccultcalc/spk_reader.h"
#include "ioccultcalc/spice_spk_reader.h"
#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/ra15_integrator.hpp"
#include <cmath>
#include <iostream>
#include <chrono>
#include <map>

namespace ioccultcalc {

// Costanti
static constexpr double AU_TO_KM = 149597870.7;
static constexpr double DAY_TO_SEC = 86400.0;
static constexpr double GM_SUN = 1.32712440041279419e11;  // km³/s²
static constexpr double GM_SUN_AU = GM_SUN / (AU_TO_KM * AU_TO_KM * AU_TO_KM) * (DAY_TO_SEC * DAY_TO_SEC);  // AU³/day²
static constexpr double C_LIGHT_AU_DAY = 299792.458 * DAY_TO_SEC / AU_TO_KM;  // AU/day

// GM pianeti in AU³/day² (valori da JPL DE441)
// https://ssd.jpl.nasa.gov/astro_par.html
struct PlanetaryGM {
    static constexpr double MERCURY = 2.2031868551e4 / (AU_TO_KM * AU_TO_KM * AU_TO_KM) * (DAY_TO_SEC * DAY_TO_SEC);
    static constexpr double VENUS   = 3.24858592000e5 / (AU_TO_KM * AU_TO_KM * AU_TO_KM) * (DAY_TO_SEC * DAY_TO_SEC);
    static constexpr double EARTH   = 3.98600435436e5 / (AU_TO_KM * AU_TO_KM * AU_TO_KM) * (DAY_TO_SEC * DAY_TO_SEC);
    static constexpr double MARS    = 4.2828375816e4 / (AU_TO_KM * AU_TO_KM * AU_TO_KM) * (DAY_TO_SEC * DAY_TO_SEC);
    static constexpr double JUPITER = 1.26712764100e8 / (AU_TO_KM * AU_TO_KM * AU_TO_KM) * (DAY_TO_SEC * DAY_TO_SEC);
    static constexpr double SATURN  = 3.7940585200e7 / (AU_TO_KM * AU_TO_KM * AU_TO_KM) * (DAY_TO_SEC * DAY_TO_SEC);
    static constexpr double URANUS  = 5.794556400e6 / (AU_TO_KM * AU_TO_KM * AU_TO_KM) * (DAY_TO_SEC * DAY_TO_SEC);
    static constexpr double NEPTUNE = 6.836535000e6 / (AU_TO_KM * AU_TO_KM * AU_TO_KM) * (DAY_TO_SEC * DAY_TO_SEC);
    static constexpr double MOON    = 4.902800076e3 / (AU_TO_KM * AU_TO_KM * AU_TO_KM) * (DAY_TO_SEC * DAY_TO_SEC);
};

// GM asteroidi principali in AU³/day² (da JPL Small-Body Database)
struct AsteroidGM {
    static constexpr double CERES   = 62.6284e-12;  // AU³/day²
    static constexpr double VESTA   = 17.8e-12;
    static constexpr double PALLAS  = 14.3e-12;
    static constexpr double HYGIEA  = 5.8e-12;
};

class OrbitPropagator::Impl {
public:
    JPLEphemerisReader jplReader;  // Fallback VSOP87
    SPKReader spkReader;            // JPL DE (alta precisione pianeti)
    SPICESPKReader asteroidReader;  // SPK asteroidi (CSPICE)
    bool useJPLDE;
    bool useAsteroidPerturbations;
    
    Impl() : useAsteroidPerturbations(false) {
        // Prova a caricare JPL DE441 (pianeti), fallback a DE430, poi VSOP87
        useJPLDE = spkReader.ensureFileLoaded("DE441");
        
        if (!useJPLDE) {
            std::cerr << "OrbitPropagator: DE441 not available, trying DE430..." << std::endl;
            useJPLDE = spkReader.ensureFileLoaded("DE430");
        }
        
        if (useJPLDE) {
            std::cerr << "OrbitPropagator: Using JPL " << spkReader.getVersion() 
                     << " (high precision planets)" << std::endl;
        } else {
            std::cerr << "OrbitPropagator: Using VSOP87 (fallback, ~100 km precision)" << std::endl;
            jplReader.loadFile("");  // Usa backend VSOP87
        }
        
        // Carica file asteroidi (codes_300ast contiene i primi 300 asteroidi inclusi AST17)
        useAsteroidPerturbations = asteroidReader.ensureFileLoaded("codes_300ast_20100725.bsp");
        if (useAsteroidPerturbations) {
            std::cerr << "OrbitPropagator: Asteroid perturbations enabled (AST17: 17 massive asteroids matching OrbFit)" << std::endl;
        } else {
            std::cerr << "OrbitPropagator: WARNING - Asteroid perturbations disabled (codes_300ast_20100725.bsp not found)" << std::endl;
        }
    }
    
    Vector3D getPlanetPosition(int naifId, double jd) {
        if (useJPLDE && spkReader.isLoaded()) {
            try {
                // Usa JPL DE per alta precisione (~1 km) - HELIOCENTRIC
                // DE441 fornisce posizioni eliocentriche (centro = 11 = Sun)
                return spkReader.getPosition(naifId, jd, 11);  // Heliocentric
            } catch (const std::exception& e) {
                std::cerr << "SPK read error, falling back to VSOP87: " << e.what() << std::endl;
                useJPLDE = false;  // Fallback permanente
            }
        }
        
        // Fallback: VSOP87 (~100 km precision)
        // Mappa NAIF ID a JPLBody per VSOP87
        static const std::map<int, JPLBody> naifToVSOP = {
            {1, JPLBody::MERCURY},
            {2, JPLBody::VENUS},
            {3, JPLBody::EARTH},
            {4, JPLBody::MARS},
            {5, JPLBody::JUPITER},
            {6, JPLBody::SATURN},
            {7, JPLBody::URANUS},
            {8, JPLBody::NEPTUNE}
        };
        
        auto it = naifToVSOP.find(naifId);
        if (it != naifToVSOP.end()) {
            Vector3D pos_km = jplReader.getPosition(it->second, jd);
            return Vector3D(pos_km.x / AU_TO_KM, pos_km.y / AU_TO_KM, pos_km.z / AU_TO_KM);
        }
        
        return Vector3D(0, 0, 0);
    }
};

OrbitPropagator::OrbitPropagator() 
    : pImpl(new Impl()), currentA1_(0), currentA2_(0), currentA3_(0) {
    // Opzioni di default
    options_.integrator = IntegratorType::RK4;
    options_.stepSize = 1.0;  // 1 giorno
    options_.usePlanetaryPerturbations = true;
    options_.useRelativisticCorrections = false;
}

OrbitPropagator::OrbitPropagator(const PropagatorOptions& options)
    : pImpl(new Impl()), options_(options), 
      currentA1_(0), currentA2_(0), currentA3_(0) {}

OrbitPropagator::~OrbitPropagator() {
    delete pImpl;
}

void OrbitPropagator::setOptions(const PropagatorOptions& options) {
    options_ = options;
}

OrbitState OrbitPropagator::elementsToState(const EquinoctialElements& elements) {
    // Converti elementi equinoziali in stato cartesiano (r, v)
    
    // Calcola anomalia media
    double omega_plus_Omega = atan2(elements.h, elements.k);
    double M = elements.lambda - omega_plus_Omega;
    
    // Eccentricità
    double e = sqrt(elements.h * elements.h + elements.k * elements.k);
    
    // Risolvi equazione di Keplero: E - e*sin(E) = M
    double E = M;
    for (int i = 0; i < 10; i++) {
        E = M + e * sin(E);
    }
    
    // Anomalia vera
    double nu = 2.0 * atan2(sqrt(1.0 + e) * sin(E / 2.0), 
                            sqrt(1.0 - e) * cos(E / 2.0));
    
    // Raggio
    double r = elements.a * (1.0 - e * cos(E));
    
    // Coordinate nel piano orbitale
    double x_orb = r * cos(nu);
    double y_orb = r * sin(nu);
    
    // Mean motion
    double n = sqrt(GM_SUN_AU / (elements.a * elements.a * elements.a));
    
    // Velocità nel piano orbitale
    double v_factor = n * elements.a / (1.0 - e * cos(E));
    double vx_orb = -v_factor * sin(E);
    double vy_orb = v_factor * sqrt(1.0 - e * e) * cos(E);
    
    // Trasforma al sistema inerziale usando elementi equinoziali
    double f = 1.0 + elements.p * elements.p + elements.q * elements.q;
    
    double m11 = (1.0 - elements.p * elements.p + elements.q * elements.q) / f;
    double m12 = 2.0 * elements.p * elements.q / f;
    double m21 = 2.0 * elements.p * elements.q / f;
    double m22 = (1.0 + elements.p * elements.p - elements.q * elements.q) / f;
    double m31 = -2.0 * elements.p / f;  // FIX: Changed sign for Z-axis
    double m32 = 2.0 * elements.q / f;   // FIX: Changed sign for Z-axis
    
    double cos_wp = cos(omega_plus_Omega);
    double sin_wp = sin(omega_plus_Omega);
    
    double x_ref = x_orb * cos_wp - y_orb * sin_wp;
    double y_ref = x_orb * sin_wp + y_orb * cos_wp;
    
    double vx_ref = vx_orb * cos_wp - vy_orb * sin_wp;
    double vy_ref = vx_orb * sin_wp + vy_orb * cos_wp;
    
    Vector3D position, velocity;
    position.x = m11 * x_ref + m12 * y_ref;
    position.y = m21 * x_ref + m22 * y_ref;
    position.z = m31 * x_ref + m32 * y_ref;
    
    velocity.x = m11 * vx_ref + m12 * vy_ref;
    velocity.y = m21 * vx_ref + m22 * vy_ref;
    velocity.z = m31 * vx_ref + m32 * vy_ref;
    
    return OrbitState(elements.epoch, position, velocity);
}

Vector3D OrbitPropagator::computeAcceleration(const JulianDate& jd,
                                             const Vector3D& pos,
                                             const Vector3D& vel) {
    // Accelerazione gravitazionale del Sole (frame heliocentric)
    double r = pos.magnitude();
    double r3 = r * r * r;
    
    Vector3D acc;
    acc.x = -GM_SUN_AU * pos.x / r3;
    acc.y = -GM_SUN_AU * pos.y / r3;
    acc.z = -GM_SUN_AU * pos.z / r3;
    
    // Perturbazioni planetarie
    if (options_.usePlanetaryPerturbations) {
        // Per ogni pianeta maggiore
        struct Planet {
            JPLBody body;
            double GM;
        };
        
        Planet planets[] = {
            {JPLBody::JUPITER, PlanetaryGM::JUPITER},
            {JPLBody::SATURN, PlanetaryGM::SATURN},
            {JPLBody::EARTH, PlanetaryGM::EARTH},
            {JPLBody::VENUS, PlanetaryGM::VENUS},
            {JPLBody::MARS, PlanetaryGM::MARS},
            {JPLBody::URANUS, PlanetaryGM::URANUS},
            {JPLBody::NEPTUNE, PlanetaryGM::NEPTUNE},
            {JPLBody::MERCURY, PlanetaryGM::MERCURY}  // Mercurio per completezza
        };
        
        for (const auto& planet : planets) {
            // Posizione pianeta (AU) - usa JPL DE se disponibile
            // Mappa JPLBody a NAIF ID
            int naifId = (planet.body == JPLBody::MERCURY) ? 1 :
                        (planet.body == JPLBody::VENUS) ? 2 :
                        (planet.body == JPLBody::EARTH) ? 3 :
                        (planet.body == JPLBody::MARS) ? 4 :
                        (planet.body == JPLBody::JUPITER) ? 5 :
                        (planet.body == JPLBody::SATURN) ? 6 :
                        (planet.body == JPLBody::URANUS) ? 7 :
                        (planet.body == JPLBody::NEPTUNE) ? 8 : 0;
            
            Vector3D planetPos = pImpl->getPlanetPosition(naifId, jd.jd);
            
            // Vettore asteroide->pianeta
            Vector3D d = planetPos - pos;
            double d_mag = d.magnitude();
            double d3 = d_mag * d_mag * d_mag;
            
            // Distanza pianeta dal Sole
            double rp = planetPos.magnitude();
            double rp3 = rp * rp * rp;
            
            // Perturbazione (termine diretto + termine indiretto)
            acc.x += planet.GM * (d.x / d3 - planetPos.x / rp3);
            acc.y += planet.GM * (d.y / d3 - planetPos.y / rp3);
            acc.z += planet.GM * (d.z / d3 - planetPos.z / rp3);
        }
        
        // TODO: Luna - DE441 non contiene corpo 301 (Moon) separatamente
        // Necessita file SPK planetario completo con satelliti
        
        // Perturbazioni asteroidali (usando CSPICE + file SPK sb441-n16)
        if (pImpl->useAsteroidPerturbations) {
            struct Asteroid {
                int naifId;
                double GM;
                const char* name;
            };
            
            // 17 asteroidi massivi AST17 (stesso set usato da OrbFit)
            // GM da Hilton 1997 (km³/s²), convertiti in AU³/day²
            // Matching OrbFit's AST17 configuration
            Asteroid asteroids[] = {
                {2000001, 62.6284e-12, "Ceres"},         // 1 Ceres (4.76E-10 Msun)
                {2000002, 14.3e-12, "Pallas"},           // 2 Pallas (1.08E-10 Msun)
                {2000003, 1.82e-12, "Juno"},             // 3 Juno (1.49E-11 Msun)
                {2000004, 17.8e-12, "Vesta"},            // 4 Vesta (1.34E-10 Msun)
                {2000006, 0.85e-12, "Hebe"},             // 6 Hebe (7.0E-12 Msun)
                {2000007, 0.79e-12, "Iris"},             // 7 Iris (6.48E-12 Msun)
                {2000010, 5.8e-12, "Hygiea"},            // 10 Hygiea (4.39E-11 Msun)
                {2000015, 0.38e-12, "Eunomia"},          // 15 Eunomia (2.9E-12 Msun)
                {2000016, 0.36e-12, "Psyche"},           // 16 Psyche (2.7E-12 Msun)
                {2000029, 0.13e-12, "Amphitrite"},       // 29 Amphitrite (1.0E-12 Msun)
                {2000052, 0.34e-12, "Europa"},           // 52 Europa (2.6E-12 Msun)
                {2000065, 0.15e-12, "Cybele"},           // 65 Cybele (1.1E-12 Msun)
                {2000087, 0.15e-12, "Sylvia"},           // 87 Sylvia (1.1E-12 Msun)
                {2000088, 0.16e-12, "Thisbe"},           // 88 Thisbe (1.2E-12 Msun)
                {2000511, 0.39e-12, "Davida"},           // 511 Davida (3.0E-12 Msun)
                {2000704, 0.34e-12, "Interamnia"},       // 704 Interamnia (2.6E-12 Msun)
                {2134340, 871.0e-12, "Pluto"}            // 134340 Pluto (6.58E-9 Msun)
            };
            
            for (const auto& ast : asteroids) {
                try {
                    // Posizione rispetto al Sole (heliocentric frame)
                    Vector3D astPos = pImpl->asteroidReader.getPosition(ast.naifId, jd.jd, 10);
                    
                    // Vettore asteroide_target -> asteroide_perturbatore
                    Vector3D d = astPos - pos;
                    double d_mag = d.magnitude();
                    
                    if (d_mag < 1e-10) continue;  // Skip se troppo vicino
                    
                    double d3 = d_mag * d_mag * d_mag;
                    
                    // Distanza perturbatore dal Sole
                    double rp = astPos.magnitude();
                    double rp3 = rp * rp * rp;
                    
                    // Perturbazione (termine diretto + termine indiretto)
                    acc.x += ast.GM * (d.x / d3 - astPos.x / rp3);
                    acc.y += ast.GM * (d.y / d3 - astPos.y / rp3);
                    acc.z += ast.GM * (d.z / d3 - astPos.z / rp3);
                    
                } catch (const std::exception& e) {
                    // Asteroide non disponibile o fuori range temporale, skip silently
                    static bool logged = false;
                    if (!logged) {
                        std::cerr << "Note: Some asteroids not available in SPK file" << std::endl;
                        logged = true;
                    }
                }
            }
        }
    }
    
    // Pressione radiazione solare (importante per NEAs e piccoli asteroidi)
    if (options_.useSolarRadiationPressure) {
        // Frame heliocentric: Sole è all'origine
        double r_sun = pos.magnitude();
        
        // Costante solare: 1361 W/m² a 1 AU
        constexpr double SOLAR_CONSTANT = 1361.0;  // W/m²
        constexpr double C_LIGHT_MS = 299792458.0;  // m/s
        constexpr double AU_TO_M = 149597870700.0;  // m
        
        // Pressione = L_sun / (4π c r²)
        double pressure = SOLAR_CONSTANT / (C_LIGHT_MS * r_sun * r_sun);  // N/m² normalizzato per AU²
        
        // Direzione radiale dal Sole (opposta a pos)
        Vector3D sunDir(-pos.x / r_sun, -pos.y / r_sun, -pos.z / r_sun);
        
        // Accelerazione SRP in AU/day²
        constexpr double M_TO_AU = 1.0 / AU_TO_M;
        constexpr double DAY_TO_S = 86400.0;
        double srp_factor = pressure * options_.areaToMassRatio * M_TO_AU * DAY_TO_S * DAY_TO_S;
        
        acc.x += srp_factor * sunDir.x;
        acc.y += srp_factor * sunDir.y;
        acc.z += srp_factor * sunDir.z;
    }
    
    // Correzioni relativistiche (post-Newtoniane)
    if (options_.useRelativisticCorrections) {
        // Frame heliocentric: usa direttamente pos
        double v2 = vel.x * vel.x + vel.y * vel.y + vel.z * vel.z;
        double rdotv = pos.x * vel.x + pos.y * vel.y + pos.z * vel.z;
        double c2 = C_LIGHT_AU_DAY * C_LIGHT_AU_DAY;
        
        // Termine relativistico (1PN)
        double factor = GM_SUN_AU / (r3 * c2);
        double term1 = 4.0 * GM_SUN_AU / r - v2;
        double term2 = 4.0 * rdotv;
        
        acc.x += factor * (term1 * pos.x + term2 * vel.x);
        acc.y += factor * (term1 * pos.y + term2 * vel.y);
        acc.z += factor * (term1 * pos.z + term2 * vel.z);
    }
    
    // Accelerazioni non gravitazionali (VFCC17 model: Yarkovsky, outgassing)
    if (options_.useNonGravitational && 
        (currentA1_ != 0 || currentA2_ != 0 || currentA3_ != 0)) {
        
        // Frame RSW (Radial-Transverse-Normal / RTN)
        // r̂ = radial (direzione Sole→asteroide normalizzata)
        Vector3D r_hat(pos.x / r, pos.y / r, pos.z / r);
        
        // Momento angolare h = r × v
        Vector3D h(pos.y * vel.z - pos.z * vel.y,
                   pos.z * vel.x - pos.x * vel.z,
                   pos.x * vel.y - pos.y * vel.x);
        
        // n̂ = normal (perpendicolare al piano orbitale)
        double h_mag = std::sqrt(h.x * h.x + h.y * h.y + h.z * h.z);
        if (h_mag > 1e-15) {  // Evita divisione per zero
            Vector3D n_hat(h.x / h_mag, h.y / h_mag, h.z / h_mag);
            
            // t̂ = transverse (n̂ × r̂, direzione moto nel piano orbitale)
            Vector3D t_hat(n_hat.y * r_hat.z - n_hat.z * r_hat.y,
                           n_hat.z * r_hat.x - n_hat.x * r_hat.z,
                           n_hat.x * r_hat.y - n_hat.y * r_hat.x);
            
            // Accelerazione non gravitazionale: a_ng = A1*r̂ + A2*t̂ + A3*n̂
            acc.x += currentA1_ * r_hat.x + currentA2_ * t_hat.x + currentA3_ * n_hat.x;
            acc.y += currentA1_ * r_hat.y + currentA2_ * t_hat.y + currentA3_ * n_hat.y;
            acc.z += currentA1_ * r_hat.z + currentA2_ * t_hat.z + currentA3_ * n_hat.z;
        }
    }
    
    return acc;
}

OrbitState OrbitPropagator::integrateRK4(const OrbitState& state0, double dt) {
    // Runge-Kutta 4° ordine per sistemi del 2° ordine
    // Notazione: y = posizione, v = velocità, a = accelerazione = y''
    
    const Vector3D& r0 = state0.position;
    const Vector3D& v0 = state0.velocity;
    const JulianDate& t0 = state0.epoch;
    
    // Stage 1: k1, l1
    Vector3D k1 = computeAcceleration(t0, r0, v0);  // k1 = a(t0, r0, v0)
    Vector3D l1 = v0;                                // l1 = v0
    
    // Stage 2: k2, l2 (punto medio usando k1, l1)
    JulianDate t2(t0.jd + dt / 2.0);
    Vector3D r2 = r0 + l1 * (dt / 2.0);              // r2 = r0 + v0*h/2
    Vector3D v2 = v0 + k1 * (dt / 2.0);              // v2 = v0 + a1*h/2
    Vector3D k2 = computeAcceleration(t2, r2, v2);   // k2 = a(t+h/2, r2, v2)
    Vector3D l2 = v2;                                // l2 = v2
    
    // Stage 3: k3, l3 (punto medio usando k2, l2)
    Vector3D r3 = r0 + l2 * (dt / 2.0);              // r3 = r0 + v2*h/2
    Vector3D v3 = v0 + k2 * (dt / 2.0);              // v3 = v0 + a2*h/2
    Vector3D k3 = computeAcceleration(t2, r3, v3);   // k3 = a(t+h/2, r3, v3)
    Vector3D l3 = v3;                                // l3 = v3
    
    // Stage 4: k4, l4 (punto finale usando k3, l3)
    JulianDate t4(t0.jd + dt);
    Vector3D r4 = r0 + l3 * dt;                      // r4 = r0 + v3*h
    Vector3D v4 = v0 + k3 * dt;                      // v4 = v0 + a3*h
    Vector3D k4 = computeAcceleration(t4, r4, v4);   // k4 = a(t+h, r4, v4)
    Vector3D l4 = v4;                                // l4 = v4
    
    // Combinazione finale RK4
    OrbitState result;
    result.epoch = t4;
    result.position = r0 + (l1 + l2 * 2.0 + l3 * 2.0 + l4) * (dt / 6.0);
    result.velocity = v0 + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
    
    return result;
}

OrbitState OrbitPropagator::integrateRA15(const OrbitState& state0, const JulianDate& targetEpoch) {
    // Configura opzioni RA15
    RA15Options ra15opts;
    ra15opts.llev = 10;                              // ss = 1e-10 (precisione massima)
    ra15opts.h_init = options_.stepSize;             // Step iniziale dall'utente
    ra15opts.eprk = options_.tolerance;              // Tolleranza convergenza
    ra15opts.lit1 = 20;                              // Max iterazioni primo step (OrbFit usa 20)
    ra15opts.lit2 = 20;                              // Max iterazioni step successivi (OrbFit usa 20)
    ra15opts.fixed_step = false;                     // Step adattivo
    ra15opts.max_steps = options_.maxSteps;
    ra15opts.verbose = false;                        // Quiet per produzione
    
    // Crea integratore RA15
    RA15Integrator integrator(ra15opts);
    
    // Funzione accelerazione per RA15 (lambda con capture di this)
    auto accel_func = [this](const JulianDate& t, const Vector3D& pos, const Vector3D& vel) -> Vector3D {
        return this->computeAcceleration(t, pos, vel);
    };
    
    // Integra (RA15 gestisce internamente loop e step adattivo)
    OrbitState result = integrator.integrate(state0, targetEpoch, accel_func);
    
    // Aggiorna statistiche
    auto stats = integrator.getStatistics();
    lastStats_.nSteps = stats.num_steps;
    lastStats_.nEvaluations = stats.num_evals;
    lastStats_.nRejections = stats.num_rejections;
    lastStats_.finalStepSize = stats.final_stepsize;
    lastStats_.avgIterations = stats.avg_iterations;
    
    return result;
}

OrbitState OrbitPropagator::propagate(const EquinoctialElements& initialElements,
                                     const JulianDate& targetEpoch) {
    // Salva parametri non gravitazionali dagli elementi
    currentA1_ = initialElements.A1;
    currentA2_ = initialElements.A2;
    currentA3_ = initialElements.A3;
    
    OrbitState initialState = elementsToState(initialElements);
    initialState.A1 = initialElements.A1;
    initialState.A2 = initialElements.A2;
    initialState.A3 = initialElements.A3;
    
    return propagate(initialState, targetEpoch);
}

OrbitState OrbitPropagator::propagate(const OrbitState& initialState,
                                     const JulianDate& targetEpoch) {
    auto start = std::chrono::high_resolution_clock::now();
    
    // Salva parametri non gravitazionali dallo stato
    currentA1_ = initialState.A1;
    currentA2_ = initialState.A2;
    currentA3_ = initialState.A3;
    
    double dt_total = targetEpoch.jd - initialState.epoch.jd;
    double dt_step = options_.stepSize;
    
    // Determina segno dello step
    if (dt_total < 0) dt_step = -dt_step;
    
    OrbitState currentState = initialState;
    int nSteps = 0;
    double maxError = 0;
    
    // Integrazione
    while (fabs(currentState.epoch.jd - targetEpoch.jd) > 1e-10) {
        double remaining = targetEpoch.jd - currentState.epoch.jd;
        
        // Step corrente (usa sempre step fisso per RK4, non variabile)
        double current_step = dt_step;
        
        // Ultimo step: aggiusta per arrivare esattamente al target
        if (fabs(remaining) < fabs(dt_step)) {
            current_step = remaining;
        }
        
        // Integra un passo
        if (options_.integrator == IntegratorType::RK4) {
            currentState = integrateRK4(currentState, current_step);
        } else if (options_.integrator == IntegratorType::RA15) {
            // RA15 gestisce internamente l'integrazione completa
            // Usa RA15 per l'intera propagazione in un colpo solo
            currentState = integrateRA15(initialState, targetEpoch);
            break;  // RA15 ha già fatto tutto
        } else {
            // Altri integratori TODO
            currentState = integrateRK4(currentState, current_step);
        }
        
        nSteps++;
        
        // Progress callback
        if (progressCallback_ && nSteps % 100 == 0) {
            double progress = fabs(currentState.epoch.jd - initialState.epoch.jd) / fabs(dt_total);
            progressCallback_(progress, currentState);
        }
        
        // Safety: max steps
        if (nSteps > options_.maxSteps) {
            std::cerr << "WARNING: Max steps reached in propagation\n";
            break;
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    // Salva statistiche
    lastStats_.nSteps = nSteps;
    lastStats_.finalStepSize = dt_step;
    lastStats_.maxError = maxError;
    lastStats_.computeTime = elapsed.count();
    
    return currentState;
}

std::vector<OrbitState> OrbitPropagator::propagateWithOutput(
    const OrbitState& initialState,
    const JulianDate& targetEpoch,
    double outputStep) {
    
    std::vector<OrbitState> results;
    results.push_back(initialState);
    
    double dt_total = targetEpoch.jd - initialState.epoch.jd;
    int nOutputs = static_cast<int>(fabs(dt_total / outputStep));
    
    for (int i = 1; i <= nOutputs; i++) {
        double jd = initialState.epoch.jd + i * outputStep;
        if ((dt_total > 0 && jd > targetEpoch.jd) ||
            (dt_total < 0 && jd < targetEpoch.jd)) {
            jd = targetEpoch.jd;
        }
        
        OrbitState state = propagate(initialState, JulianDate(jd));
        results.push_back(state);
        
        if (jd == targetEpoch.jd) break;
    }
    
    return results;
}

OrbitPropagator::PropagationStats OrbitPropagator::getLastStats() const {
    return lastStats_;
}

void OrbitPropagator::setProgressCallback(std::function<void(double, const OrbitState&)> callback) {
    progressCallback_ = callback;
}

} // namespace ioccultcalc
