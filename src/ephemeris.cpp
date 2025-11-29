#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/coordinates.h"
#include "ioccultcalc/spice_spk_reader.h"
#include "ioccultcalc/orbit_propagator.h"
#include <cmath>
#include <stdexcept>
#include <mutex>

namespace ioccultcalc {

// Thread-safe Earth position cache
static std::mutex s_cacheMutex;
static Ephemeris::EarthPositionFunc s_earthPositionCache = nullptr;

void Ephemeris::setEarthPositionCache(EarthPositionFunc cacheFunc) {
    std::lock_guard<std::mutex> lock(s_cacheMutex);
    s_earthPositionCache = cacheFunc;
}

void Ephemeris::clearEarthPositionCache() {
    std::lock_guard<std::mutex> lock(s_cacheMutex);
    s_earthPositionCache = nullptr;
}

bool Ephemeris::hasEarthPositionCache() {
    std::lock_guard<std::mutex> lock(s_cacheMutex);
    return s_earthPositionCache != nullptr;
}

// Costante gravitazionale gaussiana
constexpr double GAUSS_K = 0.01720209895; // AU^(3/2) / day

Ephemeris::Ephemeris() 
    : propagatorInitialized_(false) {}

Ephemeris::Ephemeris(const EquinoctialElements& elements) 
    : elements_(elements), propagatorInitialized_(false) {}

Ephemeris::~Ephemeris() = default;

void Ephemeris::setElements(const EquinoctialElements& elements) {
    elements_ = elements;
    // Reset propagator per nuovi elementi
    propagatorInitialized_ = false;
    propagator_.reset();
}

void Ephemeris::setOptions(const EphemerisOptions& options) {
    options_ = options;
    // Reset propagator se cambiano opzioni
    if (propagatorInitialized_) {
        propagatorInitialized_ = false;
        propagator_.reset();
    }
}

void Ephemeris::enablePerturbations(bool enable) {
    if (options_.usePerturbations != enable) {
        options_.usePerturbations = enable;
        propagatorInitialized_ = false;
        propagator_.reset();
    }
}

void Ephemeris::initializePropagator() {
    if (propagatorInitialized_) return;
    
    PropagatorOptions propOpts;
    propOpts.integrator = IntegratorType::RK4;
    propOpts.stepSize = options_.stepSize;
    propOpts.usePlanetaryPerturbations = true;  // Sempre attive quando usiamo propagator
    propOpts.useRelativisticCorrections = options_.useRelativisticCorrections;
    propOpts.tolerance = 1e-12;
    
    propagator_ = std::make_unique<OrbitPropagator>(propOpts);
    propagatorInitialized_ = true;
    
    // Reset cache ultimo stato propagato
    lastPropagatedEpoch_.jd = 0;
}

EphemerisData Ephemeris::compute(const JulianDate& jd) {
    EphemerisData data;
    data.jd = jd;
    
    // Propaga l'orbita (con o senza perturbazioni)
    if (options_.usePerturbations) {
        propagateOrbitWithPerturbations(jd, data.heliocentricPos, data.heliocentricVel);
    } else {
        propagateOrbit(jd, data.heliocentricPos, data.heliocentricVel);
    }
    
    // Posizione della Terra
    Vector3D earthPos = getEarthPosition(jd);
    Vector3D earthVel = getEarthVelocity(jd);
    
    // Posizione geocentrica dell'asteroide
    Vector3D geocentricVec = data.heliocentricPos - earthPos;
    data.distance = geocentricVec.magnitude();
    data.geocentricPos = Coordinates::cartesianToEquatorial(geocentricVec);
    
    // Posizione del Sole geocentrica
    Vector3D sunPos = getSunPosition(jd);
    
    // Elongazione solare
    Vector3D asteroidDir = geocentricVec.normalize();
    Vector3D sunDir = sunPos.normalize();
    data.elongation = acos(asteroidDir.dot(sunDir)) * RAD_TO_DEG;
    
    // Angolo di fase
    Vector3D toEarth = geocentricVec * -1.0;
    double cosPhase = data.heliocentricPos.normalize().dot(toEarth.normalize());
    data.phase = acos(cosPhase) * RAD_TO_DEG;
    
    // Magnitudine
    double r = data.heliocentricPos.magnitude(); // AU
    data.magnitude = calculateMagnitude(r, data.distance, data.phase);
    
    return data;
}

std::vector<EphemerisData> Ephemeris::computeRange(const JulianDate& startJD,
                                                   const JulianDate& endJD,
                                                   double stepDays) {
    std::vector<EphemerisData> results;
    
    for (double jd = startJD.jd; jd <= endJD.jd; jd += stepDays) {
        results.push_back(compute(JulianDate(jd)));
    }
    
    return results;
}

void Ephemeris::propagateOrbit(const JulianDate& targetJD,
                               Vector3D& helioPos, Vector3D& helioVel) {
    // Propaga usando elementi equinoziali
    
    // Tempo dall'epoca in giorni
    double dt = targetJD.jd - elements_.epoch.jd;
    
    // Mean motion
    double n = GAUSS_K / sqrt(elements_.a * elements_.a * elements_.a);
    
    // Mean longitude al tempo target
    double lambda_t = elements_.lambda + n * dt;
    while (lambda_t < 0) lambda_t += 2.0 * M_PI;
    while (lambda_t >= 2.0 * M_PI) lambda_t -= 2.0 * M_PI;
    
    // Calcola anomalia media
    double omega_plus_Omega = atan2(elements_.h, elements_.k);
    double M = lambda_t - omega_plus_Omega;
    
    // Eccentricità
    double e = sqrt(elements_.h * elements_.h + elements_.k * elements_.k);
    
    // Risolvi equazione di Keplero
    double E = solveKeplerEquation(M, e);
    
    // Anomalia vera
    double cosNu = (cos(E) - e) / (1.0 - e * cos(E));
    double sinNu = sqrt(1.0 - e * e) * sin(E) / (1.0 - e * cos(E));
    double nu = atan2(sinNu, cosNu);
    
    // Raggio vettore
    double r = elements_.a * (1.0 - e * cos(E));
    
    // Coordinate nel piano orbitale
    double x_orb = r * cos(nu);
    double y_orb = r * sin(nu);
    
    // Velocità nel piano orbitale
    double v_factor = GAUSS_K * sqrt(elements_.a) / r;
    double vx_orb = -v_factor * sin(E);
    double vy_orb = v_factor * sqrt(1.0 - e * e) * cos(E);
    
    // Converti elementi equinoziali in elementi classici
    // p = tan(i/2) * sin(Omega)
    // q = tan(i/2) * cos(Omega)
    double tan_half_i = std::sqrt(elements_.p * elements_.p + elements_.q * elements_.q);
    double i = 2.0 * std::atan(tan_half_i);  // Inclinazione
    double Omega = std::atan2(elements_.p, elements_.q);  // Longitudine nodo ascendente
    double omega = omega_plus_Omega - Omega;  // Argomento del periapside
    
    // Trasforma dal piano orbitale al frame eclittico usando rotazioni 3D
    // Sequenza: R_z(-Omega) * R_x(-i) * R_z(-omega)
    double cos_o = std::cos(omega);
    double sin_o = std::sin(omega);
    double cos_O = std::cos(Omega);
    double sin_O = std::sin(Omega);
    double cos_i = std::cos(i);
    double sin_i = std::sin(i);
    
    // Matrice di rotazione combinata (piano orbitale -> eclittico)
    double P11 = cos_O * cos_o - sin_O * sin_o * cos_i;
    double P12 = -cos_O * sin_o - sin_O * cos_o * cos_i;
    double P21 = sin_O * cos_o + cos_O * sin_o * cos_i;
    double P22 = -sin_O * sin_o + cos_O * cos_o * cos_i;
    double P31 = sin_i * sin_o;
    double P32 = sin_i * cos_o;
    
    // Posizione eclittica
    double x_ecl = P11 * x_orb + P12 * y_orb;
    double y_ecl = P21 * x_orb + P22 * y_orb;
    double z_ecl = P31 * x_orb + P32 * y_orb;
    
    // Velocità eclittica
    double vx_ecl = P11 * vx_orb + P12 * vy_orb;
    double vy_ecl = P21 * vx_orb + P22 * vy_orb;
    double vz_ecl = P31 * vx_orb + P32 * vy_orb;
    
    // Converti da eclittico a equatoriale (J2000)
    // Rotazione attorno asse X di ε = 23.4392811° (obliquità eclittica J2000)
    constexpr double OBLIQUITY_J2000 = 23.4392911 * M_PI / 180.0;  // rad
    double cos_eps = std::cos(OBLIQUITY_J2000);
    double sin_eps = std::sin(OBLIQUITY_J2000);
    
    helioPos.x = x_ecl;
    helioPos.y = y_ecl * cos_eps - z_ecl * sin_eps;
    helioPos.z = y_ecl * sin_eps + z_ecl * cos_eps;
    
    helioVel.x = vx_ecl;
    helioVel.y = vy_ecl * cos_eps - vz_ecl * sin_eps;
    helioVel.z = vy_ecl * sin_eps + vz_ecl * cos_eps;
}

void Ephemeris::propagateOrbitWithPerturbations(const JulianDate& targetJD,
                                                  Vector3D& helioPos, Vector3D& helioVel) {
    // Inizializza propagator se necessario
    if (!propagatorInitialized_) {
        initializePropagator();
    }
    
    OrbitState result;
    
    // OTTIMIZZAZIONE: Se abbiamo già propagato a un'epoca vicina, 
    // usa quello stato come punto di partenza invece di ripartire dall'epoca elementi
    if (lastPropagatedEpoch_.jd > 0) {
        double dt_from_last = std::abs(targetJD.jd - lastPropagatedEpoch_.jd);
        double dt_from_epoch = std::abs(targetJD.jd - elements_.epoch.jd);
        
        // Se l'ultimo stato è più vicino dell'epoca elementi, usa quello
        if (dt_from_last < dt_from_epoch && dt_from_last < 10.0) {  // max 10 giorni
            OrbitState lastState(lastPropagatedEpoch_, lastPropagatedPos_, lastPropagatedVel_);
            result = propagator_->propagate(lastState, targetJD);
        } else {
            result = propagator_->propagate(elements_, targetJD);
        }
    } else {
        result = propagator_->propagate(elements_, targetJD);
    }
    
    // Cache ultimo stato (in frame eclittico)
    lastPropagatedEpoch_ = targetJD;
    lastPropagatedPos_ = result.position;
    lastPropagatedVel_ = result.velocity;
    
    // NOTA: OrbitPropagator restituisce posizioni in frame eclittico
    // Dobbiamo convertire a equatoriale J2000 come fa propagateOrbit()
    
    // Converti da eclittico a equatoriale J2000
    constexpr double OBLIQUITY_J2000 = 23.4392911 * M_PI / 180.0;
    double cos_eps = std::cos(OBLIQUITY_J2000);
    double sin_eps = std::sin(OBLIQUITY_J2000);
    
    // result.position è in eclittico
    helioPos.x = result.position.x;
    helioPos.y = result.position.y * cos_eps - result.position.z * sin_eps;
    helioPos.z = result.position.y * sin_eps + result.position.z * cos_eps;
    
    helioVel.x = result.velocity.x;
    helioVel.y = result.velocity.y * cos_eps - result.velocity.z * sin_eps;
    helioVel.z = result.velocity.y * sin_eps + result.velocity.z * cos_eps;
}

double Ephemeris::solveKeplerEquation(double M, double e, double tolerance) {
    // Risolve E - e*sin(E) = M usando Newton-Raphson
    
    double E = M; // Prima approssimazione
    if (e > 0.8) {
        E = M_PI; // Migliore starting point per alte eccentricità
    }
    
    for (int i = 0; i < 100; i++) {
        double f = E - e * sin(E) - M;
        double fp = 1.0 - e * cos(E);
        double delta = f / fp;
        E -= delta;
        
        if (fabs(delta) < tolerance) {
            break;
        }
    }
    
    return E;
}

double Ephemeris::calculateMagnitude(double r, double delta, double phaseAngle) {
    // Formula HG per la magnitudine
    double phi1 = exp(-3.33 * pow(tan(phaseAngle * DEG_TO_RAD / 2.0), 0.63));
    double phi2 = exp(-1.87 * pow(tan(phaseAngle * DEG_TO_RAD / 2.0), 1.22));
    
    double H = elements_.H;
    double G = elements_.G;
    
    double mag = H + 5.0 * log10(r * delta) - 2.5 * log10((1.0 - G) * phi1 + G * phi2);
    
    return mag;
}

Vector3D Ephemeris::getSunPosition(const JulianDate& jd) {
    // Restituisce posizione Sole geocentrica (opposto della Terra eliocentrica)
    return getEarthPosition(jd) * -1.0;
}

Vector3D Ephemeris::getEarthPosition(const JulianDate& jd) {
    // THREAD-SAFE: Se è disponibile una cache pre-calcolata, usala
    // Questo evita problemi di thread-safety con SPICE
    {
        std::lock_guard<std::mutex> lock(s_cacheMutex);
        if (s_earthPositionCache != nullptr) {
            return s_earthPositionCache(jd.jd);
        }
    }
    
    // Prova a usare SPK (JPL DE441/DE440) per massima precisione
    static SPICESPKReader spkReader;
    static bool spkInitialized = false;
    
    if (!spkInitialized) {
        // Prova a caricare DE441, DE440s o DE440
        if (spkReader.ensureFileLoaded("de441.bsp") || 
            spkReader.ensureFileLoaded("de440s.bsp") ||
            spkReader.ensureFileLoaded("de440.bsp")) {
            spkInitialized = true;
        }
    }
    
    if (spkReader.isLoaded()) {
        try {
            // NAIF ID: 399 = Terra, 10 = Sole, 0 = SSB (barycenter)
            // NOTA: SPICE restituisce in ECLIPJ2000, dobbiamo convertire a EQUATORIALE J2000
            // per coerenza con propagateOrbit() che restituisce helioPos in equatoriale
            Vector3D pos_ecl;
            
            // Prova prima con Sole come centro (heliocentric)
            try {
                auto [pos, vel] = spkReader.getState(399, jd.jd, 10);
                pos_ecl = pos;
            } catch (...) {
                // Se fallisce, prova SSB e sottrai posizione Sole
                auto [earthPos, earthVel] = spkReader.getState(399, jd.jd, 0);
                auto [sunPos, sunVel] = spkReader.getState(10, jd.jd, 0);
                pos_ecl = earthPos - sunPos;
            }
            
            // Converti da eclittico a equatoriale J2000
            // Rotazione attorno asse X di ε = 23.4392911° (obliquità eclittica J2000)
            constexpr double OBLIQUITY_J2000 = 23.4392911 * M_PI / 180.0;
            double cos_eps = std::cos(OBLIQUITY_J2000);
            double sin_eps = std::sin(OBLIQUITY_J2000);
            
            Vector3D pos_eq;
            pos_eq.x = pos_ecl.x;
            pos_eq.y = pos_ecl.y * cos_eps - pos_ecl.z * sin_eps;
            pos_eq.z = pos_ecl.y * sin_eps + pos_ecl.z * cos_eps;
            
            return pos_eq;
        } catch (...) {
            // Fallback a formula analitica
        }
    }
    
    // FALLBACK: Formula analitica semplificata
    // NOTA: Questa ha precisione limitata (~0.4 AU di errore!)
    // Per occultazioni asteroidali serve JPL DE
    
    double T = (jd.jd - 2451545.0) / 36525.0;
    
    // Elementi orbitali medi Terra
    double L = 280.46646 + 36000.76983 * T + 0.0003032 * T * T;
    double M = 357.52911 + 35999.05029 * T - 0.0001537 * T * T;
    
    L = fmod(L, 360.0);
    if (L < 0) L += 360.0;
    M = fmod(M, 360.0);
    if (M < 0) M += 360.0;
    
    double M_rad = M * DEG_TO_RAD;
    double e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T * T;
    
    // Equazione del centro
    double C = (1.914602 - 0.004817 * T - 0.000014 * T * T) * sin(M_rad) +
               (0.019993 - 0.000101 * T) * sin(2.0 * M_rad) +
               0.000289 * sin(3.0 * M_rad);
    
    double sunLon = L + C;
    double R = (1.000001018 * (1.0 - e * e)) / (1.0 + e * cos(M_rad + C * DEG_TO_RAD));
    
    // ATTENZIONE: Questa formula è imprecisa!
    // Serve solo come fallback se SPK non disponibile
    double earthLon = sunLon + 180.0;
    if (earthLon >= 360.0) earthLon -= 360.0;
    double earthLon_rad = earthLon * DEG_TO_RAD;
    
    double x_ecl = R * cos(earthLon_rad);
    double y_ecl = R * sin(earthLon_rad);
    double z_ecl = 0.0;
    
    double eps = 23.4392911 * DEG_TO_RAD;
    
    Vector3D earthPos;
    earthPos.x = x_ecl;
    earthPos.y = y_ecl * cos(eps) - z_ecl * sin(eps);
    earthPos.z = y_ecl * sin(eps) + z_ecl * cos(eps);
    
    return earthPos;
}

Vector3D Ephemeris::getEarthVelocity(const JulianDate& jd) {
    // Prova prima SPK per velocità diretta
    static SPICESPKReader spkReader;
    static bool spkInitialized = false;
    
    if (!spkInitialized) {
        if (spkReader.ensureFileLoaded("de441.bsp") || 
            spkReader.ensureFileLoaded("de440s.bsp") ||
            spkReader.ensureFileLoaded("de440.bsp")) {
            spkInitialized = true;
        }
    }
    
    if (spkReader.isLoaded()) {
        try {
            auto [pos, vel_ecl] = spkReader.getState(399, jd.jd, 10);
            
            // Converti velocità da eclittico a equatoriale J2000
            constexpr double OBLIQUITY_J2000 = 23.4392911 * M_PI / 180.0;
            double cos_eps = std::cos(OBLIQUITY_J2000);
            double sin_eps = std::sin(OBLIQUITY_J2000);
            
            Vector3D vel_eq;
            vel_eq.x = vel_ecl.x;
            vel_eq.y = vel_ecl.y * cos_eps - vel_ecl.z * sin_eps;
            vel_eq.z = vel_ecl.y * sin_eps + vel_ecl.z * cos_eps;
            
            return vel_eq;
        } catch (...) {
            // Fallback a differenze finite
        }
    }
    
    // FALLBACK: Calcola velocità con differenze finite
    double dt = 0.1; // 0.1 giorni
    
    Vector3D pos1 = getEarthPosition(JulianDate(jd.jd - dt / 2.0));
    Vector3D pos2 = getEarthPosition(JulianDate(jd.jd + dt / 2.0));
    
    return (pos2 - pos1) * (1.0 / dt);
}

Vector3D Ephemeris::getEarthPositionWithCorrections(const JulianDate& jd, 
                                                     const Vector3D& observerPos) {
    // Ottiene posizione base da SPK con interpolazione
    Vector3D earthPos = getEarthPosition(jd);
    Vector3D earthVel = getEarthVelocity(jd);
    
    // Costanti fisiche
    constexpr double C_AU_PER_DAY = 173.1446326846693;  // Velocità luce in AU/day
    constexpr double GM_SUN = 0.000295912208286;        // GM☉ in AU³/day²
    constexpr double C_AU_PER_DAY_SQ = C_AU_PER_DAY * C_AU_PER_DAY;
    
    // === 1. ABERRAZIONE DELLA LUCE ===
    // Correzione per il tempo di propagazione della luce
    Vector3D toEarth = earthPos - observerPos;
    double distance = toEarth.magnitude();
    double lightTime = distance / C_AU_PER_DAY;
    
    // Iterazione per tempo luce (miglior accuratezza)
    Vector3D earthPosIterative = earthPos;
    for (int iter = 0; iter < 2; ++iter) {
        JulianDate jdRetarded(jd.jd - lightTime);
        earthPosIterative = getEarthPosition(jdRetarded);
        Vector3D toEarthIter = earthPosIterative - observerPos;
        lightTime = toEarthIter.magnitude() / C_AU_PER_DAY;
    }
    
    // Aberrazione stellare (velocità osservatore)
    Vector3D aberrationCorr = earthVel * (-lightTime);
    earthPos = earthPos + aberrationCorr;
    
    // === 2. CORREZIONI RELATIVISTICHE ===
    
    // 2a. Effetto Shapiro (ritardo gravitazionale)
    // Δt = (2*GM/c³) * ln[(r_e + r_o + d)/(r_e + r_o - d)]
    // dove r_e = distanza Terra-Sole, r_o = distanza Osservatore-Sole
    double r_earth = earthPos.magnitude();  // Distanza Terra-Sole
    double r_obs = observerPos.magnitude(); // Distanza Osservatore-Sole
    
    if (r_earth > 0.01 && distance > 0.001) {  // Evita divisioni per zero
        // Geometria del percorso luce
        double sum_distances = r_earth + r_obs;
        double arg1 = sum_distances + distance;
        double arg2 = std::abs(sum_distances - distance);
        
        if (arg1 > 0 && arg2 > 0) {
            // Ritardo Shapiro in giorni
            double shapiroDelay = (2.0 * GM_SUN / (C_AU_PER_DAY * C_AU_PER_DAY_SQ)) * 
                                 std::log(arg1 / arg2);
            
            // Correzione posizione: Δr = v * Δt
            Vector3D shapiroCorr = earthVel * shapiroDelay;
            earthPos = earthPos + shapiroCorr;
        }
    }
    
    // 2b. Deflexione gravitazionale (light bending)
    // Δθ = (4*GM/c²) / b, dove b = parametro impatto
    // Approssimazione: correzione proporzionale alla massa solare
    Vector3D sunToEarth = earthPos;  // Terra rispetto al Sole
    Vector3D sunToObs = observerPos; // Osservatore rispetto al Sole
    
    // Parametro impatto (distanza minima raggio luce dal Sole)
    Vector3D lightDir = (sunToObs - sunToEarth);
    double lightDist = lightDir.magnitude();
    
    if (lightDist > 0.001) {
        lightDir = lightDir * (1.0 / lightDist);
        
        // Proiezione della posizione Terra sulla direzione luce
        double projection = sunToEarth.x * lightDir.x + 
                          sunToEarth.y * lightDir.y + 
                          sunToEarth.z * lightDir.z;
        
        Vector3D closestPoint = lightDir * projection;
        Vector3D impactVector = sunToEarth - closestPoint;
        double impactParam = impactVector.magnitude();
        
        if (impactParam > 0.01) {  // Solo se non passa troppo vicino al Sole
            // Angolo di deflessione (radianti)
            double deflectionAngle = (4.0 * GM_SUN / C_AU_PER_DAY_SQ) / impactParam;
            
            // Correzione posizione (perpendicolare alla direzione luce)
            Vector3D deflectionDir = impactVector * (1.0 / impactParam);
            Vector3D deflectionCorr = deflectionDir * (deflectionAngle * distance);
            
            earthPos = earthPos + deflectionCorr;
        }
    }
    
    // === 3. INTERPOLAZIONE CHEBYSHEV (cache per query multiple) ===
    // Implementata direttamente in SPK reader per migliore efficienza
    
    return earthPos;
}

} // namespace ioccultcalc
