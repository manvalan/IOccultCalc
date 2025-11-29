/**
 * @file orbfit_force_model.cpp
 * @brief Implementazione modello di forze OrbFit-compatible
 * 
 * Traduzione fedele da force_model.f90 di OrbFit
 * Autore originale: A. Milani et al., Università di Pisa
 */

#include "ioccultcalc/orbfit_force_model.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>

namespace ioccultcalc {

using namespace OrbFitConstants;

// Obliquità eclittica J2000.0 (IAU 2006)
static constexpr double OBLIQUITY = 23.4392911 * M_PI / 180.0;
static const double COS_OBL = std::cos(OBLIQUITY);
static const double SIN_OBL = std::sin(OBLIQUITY);

// ============================================================================
// Costruttore / Distruttore
// ============================================================================

OrbFitForceModel::OrbFitForceModel(const OrbFitForceOptions& options)
    : options_(options), initialized_(false), cachedJD_(0) {
    
    // Calcola numero pianeti attivi (npla in OrbFit)
    numPlanets_ = 0;
    if (options_.includeMercury) numPlanets_++;
    if (options_.includeVenus) numPlanets_++;
    if (options_.includeEarth) numPlanets_++;
    if (options_.includeMars) numPlanets_++;
    if (options_.includeJupiter) numPlanets_++;
    if (options_.includeSaturn) numPlanets_++;
    if (options_.includeUranus) numPlanets_++;
    if (options_.includeNeptune) numPlanets_++;
    if (options_.includeMoon) numPlanets_++;
    if (options_.includePluto) numPlanets_++;
    
    // Numero totale corpi massicci
    int numAsteroids = options_.includeMassiveAsteroids ? 
                       std::min(options_.numMassiveAsteroids, NUM_MASSIVE_ASTEROIDS) : 0;
    numMass_ = numPlanets_ + numAsteroids;
    
    // Inizializza array masse GM
    GM_.resize(numMass_);
    
    // Masse planetarie (gm() in OrbFit) in AU³/day²
    int idx = 0;
    if (options_.includeMercury) GM_[idx++] = GM0 * MASS_MERCURY;
    if (options_.includeVenus)   GM_[idx++] = GM0 * MASS_VENUS;
    if (options_.includeEarth)   GM_[idx++] = GM0 * MASS_EARTH;
    if (options_.includeMars)    GM_[idx++] = GM0 * MASS_MARS;
    if (options_.includeJupiter) GM_[idx++] = GM0 * MASS_JUPITER;
    if (options_.includeSaturn)  GM_[idx++] = GM0 * MASS_SATURN;
    if (options_.includeUranus)  GM_[idx++] = GM0 * MASS_URANUS;
    if (options_.includeNeptune) GM_[idx++] = GM0 * MASS_NEPTUNE;
    if (options_.includeMoon)    GM_[idx++] = GM0 * MASS_MOON;
    if (options_.includePluto)   GM_[idx++] = GM0 * MASS_PLUTO;
    
    // Masse asteroidi
    for (int i = 0; i < numAsteroids; ++i) {
        // Escludi self-perturbation se richiesto
        if (MASSIVE_ASTEROIDS[i].number != options_.excludeAsteroidNumber) {
            GM_[idx++] = GM0 * MASSIVE_ASTEROIDS[i].mass;
        }
    }
    
    planetCache_.resize(numMass_);
}

OrbFitForceModel::~OrbFitForceModel() = default;

// ============================================================================
// Inizializzazione
// ============================================================================

bool OrbFitForceModel::initializeEphemeris(const std::string& ephemPath) {
    jplReader_ = std::make_unique<JPLEphemerisReader>();
    
    bool success = false;
    if (!ephemPath.empty()) {
        success = jplReader_->loadFile(ephemPath);
    } else {
        // Auto-download DE441
        success = jplReader_->downloadDE(JPLVersion::DE441);
    }
    
    if (success) {
        initialized_ = true;
        std::cerr << "OrbFitForceModel: Initialized with JPL DE" 
                  << static_cast<int>(jplReader_->getVersion()) << std::endl;
    }
    
    return success;
}

// ============================================================================
// Calcolo posizioni planetarie (PLANAST in OrbFit)
// ============================================================================

void OrbFitForceModel::computePlanetPositions(double jd) const {
    if (std::abs(jd - cachedJD_) < 1e-10) {
        return;  // Usa cache
    }
    
    if (!initialized_) {
        throw std::runtime_error("OrbFitForceModel: Ephemeris not initialized");
    }
    
    int idx = 0;
    
    // Ottieni posizioni da JPL (equatoriali) e converti in eclittiche
    auto getState = [&](JPLBody body, double gm) {
        PlanetaryState& state = planetCache_[idx];
        
        // JPL restituisce posizioni BARICENTRCHE, dobbiamo convertire in ELIOCENTRICHE
        Vector3D posEq = jplReader_->getPosition(body, jd);
        Vector3D velEq = jplReader_->getVelocity(body, jd);
        
        // Posizione del Sole rispetto al baricentro
        Vector3D sunPosEq = jplReader_->getPosition(JPLBody::SUN, jd);
        Vector3D sunVelEq = jplReader_->getVelocity(JPLBody::SUN, jd);
        
        // Converti in eliocentrico
        posEq = posEq - sunPosEq;
        velEq = velEq - sunVelEq;
        
        // Converti da km a AU e da km/s a AU/day
        constexpr double KM_TO_AU = 1.0 / 149597870.7;
        constexpr double KMS_TO_AUD = 86400.0 / 149597870.7;
        
        posEq = posEq * KM_TO_AU;
        velEq = velEq * KMS_TO_AUD;
        
        // Converti da equatoriale a eclittico (EQUM00 -> ECLM00 in OrbFit)
        state.position = equatorialToEcliptic(posEq);
        state.velocity = equatorialToEcliptic(velEq);
        state.distance = state.position.magnitude();
        state.GM = gm;
        
        idx++;
    };
    
    // Pianeti
    if (options_.includeMercury) getState(JPLBody::MERCURY, GM_[idx]);
    if (options_.includeVenus)   getState(JPLBody::VENUS, GM_[idx]);
    if (options_.includeEarth)   getState(JPLBody::EARTH, GM_[idx]);
    if (options_.includeMars)    getState(JPLBody::MARS, GM_[idx]);
    if (options_.includeJupiter) getState(JPLBody::JUPITER, GM_[idx]);
    if (options_.includeSaturn)  getState(JPLBody::SATURN, GM_[idx]);
    if (options_.includeUranus)  getState(JPLBody::URANUS, GM_[idx]);
    if (options_.includeNeptune) getState(JPLBody::NEPTUNE, GM_[idx]);
    if (options_.includeMoon)    getState(JPLBody::MOON, GM_[idx]);
    if (options_.includePluto)   getState(JPLBody::PLUTO, GM_[idx]);
    
    // TODO: Asteroidi massicci da SPK (sb441-n16.bsp)
    // Per ora usiamo posizioni approssimate da elementi orbitali
    
    cachedJD_ = jd;
}

// ============================================================================
// Accelerazione solare centrale
// ============================================================================

Vector3D OrbFitForceModel::solarAcceleration(const Vector3D& pos) const {
    // f = -GM0 * r / |r|³
    double r = pos.magnitude();
    if (r < 1e-15) return Vector3D(0, 0, 0);
    
    double r3 = r * r * r;
    return pos * (-GM0 / r3);
}

// ============================================================================
// Accelerazione planetaria (diretto + indiretto)
// ============================================================================

Vector3D OrbFitForceModel::planetaryAcceleration(const Vector3D& pos, double jd) const {
    computePlanetPositions(jd);
    
    Vector3D accel(0, 0, 0);
    
    for (int i = 0; i < numPlanets_; ++i) {
        const PlanetaryState& planet = planetCache_[i];
        double gm = planet.GM;
        
        // Vettore asteroide -> pianeta
        Vector3D delta = planet.position - pos;
        double d = delta.magnitude();
        if (d < 1e-15) continue;
        
        double d3 = d * d * d;
        double r3 = planet.distance * planet.distance * planet.distance;
        
        // Termine DIRETTO: attrazione pianeta -> asteroide
        // f += GM * (r_pla - r_ast) / |r_pla - r_ast|³
        accel = accel + delta * (gm / d3);
        
        // Termine INDIRETTO: correzione per riferimento eliocentrico (non inerziale)
        // f -= GM * r_pla / |r_pla|³
        accel = accel - planet.position * (gm / r3);
    }
    
    return accel;
}

// ============================================================================
// Accelerazione da asteroidi massicci
// ============================================================================

Vector3D OrbFitForceModel::asteroidAcceleration(const Vector3D& pos, double jd) const {
    if (!options_.includeMassiveAsteroids) {
        return Vector3D(0, 0, 0);
    }
    
    // TODO: Implementare lettura effemeridi asteroidi da SPK o .bep
    // Per ora ritorna zero (gli asteroidi sono già inclusi in alcuni DE)
    
    return Vector3D(0, 0, 0);
}

// ============================================================================
// GENREL: Relatività Schwarzschild (solo Sole)
// ============================================================================

Vector3D OrbFitForceModel::genrel(const Vector3D& pos, const Vector3D& vel) const {
    // Implementazione da genrel() in force_model.f90
    // Post-Newtonian ordine 1, solo termine solare
    //
    // a_rel = (GM/c²r³) * [(4GM/r - v²)*r + 4(r·v)*v]
    
    double r2 = pos.dot(pos);
    double r = std::sqrt(r2);
    if (r < 1e-15) return Vector3D(0, 0, 0);
    
    double v2 = vel.dot(vel);
    double rdotv = pos.dot(vel);
    
    double c2 = VLIGHT * VLIGHT;
    double GMoverR = GM0 / r;
    
    // Fattore comune: GM/(c² r³)
    double factor = GMoverR / (c2 * r2);
    
    // Bracket: 4*GM/r - v²
    double bracket = 4.0 * GMoverR - v2;
    
    // Accelerazione relativistica
    Vector3D aRel = pos * (factor * bracket) + vel * (4.0 * factor * rdotv);
    
    return aRel;
}

// ============================================================================
// EIHREL2: Relatività EIH completa (eq. 4-26 DESCANSO2)
// ============================================================================

Vector3D OrbFitForceModel::eihrel2(const Vector3D& pos, 
                                    const Vector3D& vel,
                                    double jd) const {
    // Implementazione fedele di eihrel2() in force_model.f90
    // Riferimento: DESCANSO Monograph No. 2, eq. (4-26)
    // http://descanso.jpl.nasa.gov/Monograph/series2/Descanso2_all.pdf
    
    computePlanetPositions(jd);
    
    Vector3D drgr(0, 0, 0);
    
    double c2 = VLIGHT * VLIGHT;
    double rast = pos.magnitude();
    if (rast < 1e-15) return Vector3D(0, 0, 0);
    
    // === Calcolo velocità baricentriche ===
    // vsun = -Σ GM[k] * v_pla[k] / GM_all
    double gmall = GM0;
    for (int k = 0; k < numPlanets_; ++k) {
        gmall += planetCache_[k].GM;
    }
    
    Vector3D vsun(0, 0, 0);
    for (int k = 0; k < numPlanets_; ++k) {
        vsun = vsun - planetCache_[k].velocity * (planetCache_[k].GM / gmall);
    }
    
    // vast = vsun + v (velocità baricentrica asteroide)
    Vector3D vast = vsun + vel;
    double vast2 = vast.dot(vast);
    double vsun2 = vsun.dot(vsun);
    
    // Velocità baricentriche pianeti
    std::vector<Vector3D> vpla(numPlanets_);
    for (int k = 0; k < numPlanets_; ++k) {
        vpla[k] = vsun + planetCache_[k].velocity;
    }
    
    // === CONTRIBUTI PLANETARI ===
    for (int j = 0; j < numPlanets_; ++j) {
        const PlanetaryState& planet = planetCache_[j];
        double gm_j = planet.GM;
        double r_j = planet.distance;
        
        // Distanza asteroide-pianeta
        Vector3D delta = pos - planet.position;
        double d_j = delta.magnitude();
        if (d_j < 1e-15) continue;
        
        // --- Prima riga: potenziali ---
        double term = -4.0 * GM0 / rast - GM0 / r_j;
        
        // Somma potenziali su asteroide
        for (int l = 0; l < numPlanets_; ++l) {
            Vector3D dl = pos - planetCache_[l].position;
            double d_l = dl.magnitude();
            if (d_l > 1e-15) term -= 4.0 * planetCache_[l].GM / d_l;
        }
        
        // Potenziali tra pianeti
        for (int k = 0; k < numPlanets_; ++k) {
            if (k == j) continue;
            Vector3D rjk = planetCache_[k].position - planet.position;
            double d_jk = rjk.magnitude();
            if (d_jk > 1e-15) term -= planetCache_[k].GM / d_jk;
        }
        
        // --- Seconda riga: velocità ---
        term += vast2 + 2.0 * vpla[j].dot(vpla[j]) - 4.0 * vast.dot(vpla[j]);
        
        // --- Terza riga: accelerazione Newtoniana sul pianeta j ---
        Vector3D g_j = planet.position * (-GM0 / (r_j * r_j * r_j));
        for (int k = 0; k < numPlanets_; ++k) {
            if (k == j) continue;
            Vector3D rjk = planetCache_[k].position - planet.position;
            double d_jk = rjk.magnitude();
            if (d_jk > 1e-15) {
                g_j = g_j + rjk * (planetCache_[k].GM / (d_jk * d_jk * d_jk));
            }
        }
        
        // Termine (n·v)² dove n = (r_ast - r_pla) / d
        double ndotv = delta.dot(vpla[j]) / d_j;
        term -= 1.5 * ndotv * ndotv;
        
        // Termine 0.5 * (r_pla - r_ast) · g_j
        term += 0.5 * (planet.position - pos).dot(g_j);
        
        // Primo contributo: term * GM_j * (r_pla - r_ast) / d³
        Vector3D contrib1 = (planet.position - pos) * (term * gm_j / (d_j * d_j * d_j));
        drgr = drgr + contrib1;
        
        // --- Quarta riga: termine dipendente dalla velocità ---
        Vector3D vastRel = vast - vpla[j];
        double bracket4 = delta.dot(vast * 4.0 - vpla[j] * 3.0);
        Vector3D contrib2 = vastRel * (gm_j * bracket4 / (d_j * d_j * d_j));
        drgr = drgr + contrib2;
        
        // --- Quinta riga: termine "third-body" ---
        Vector3D contrib3 = g_j * (3.5 * gm_j / d_j);
        drgr = drgr + contrib3;
    }
    
    // === CONTRIBUTO SOLARE (Schwarzschild-like) ===
    {
        // Prima riga
        double term = -4.0 * GM0 / rast;
        for (int l = 0; l < numPlanets_; ++l) {
            Vector3D dl = pos - planetCache_[l].position;
            double d_l = dl.magnitude();
            if (d_l > 1e-15) term -= 4.0 * planetCache_[l].GM / d_l;
        }
        for (int k = 0; k < numPlanets_; ++k) {
            term -= planetCache_[k].GM / planetCache_[k].distance;
        }
        
        // Seconda riga
        term += vast2 + 2.0 * vsun2 - 4.0 * vast.dot(vsun);
        
        // Terza riga: accelerazione Newtoniana sul Sole
        Vector3D g0(0, 0, 0);
        for (int k = 0; k < numPlanets_; ++k) {
            double r_k = planetCache_[k].distance;
            g0 = g0 + planetCache_[k].position * (planetCache_[k].GM / (r_k * r_k * r_k));
        }
        
        double ndotv = pos.dot(vsun) / rast;
        term -= 1.5 * ndotv * ndotv;
        term += 0.5 * (pos * (-1.0)).dot(g0);
        
        // Primo contributo solare
        Vector3D contrib1 = pos * (-term * GM0 / (rast * rast * rast));
        drgr = drgr + contrib1;
        
        // Quarta riga (solare)
        Vector3D vastRel = vast - vsun;
        double bracket4 = pos.dot(vast * 4.0 - vsun * 3.0);
        Vector3D contrib2 = vastRel * (GM0 * bracket4 / (rast * rast * rast));
        drgr = drgr + contrib2;
        
        // Quinta riga (solare)
        Vector3D contrib3 = g0 * (3.5 * GM0 / rast);
        drgr = drgr + contrib3;
    }
    
    // Dividi per c²
    return drgr * (1.0 / c2);
}

// ============================================================================
// J2SUN: Oblatezza del Sole
// ============================================================================

Vector3D OrbFitForceModel::j2sun(const Vector3D& pos) const {
    // Implementazione da j2sun() in force_model.f90
    // J2 del Sole = 2×10⁻⁷, raggio = 4.6527174×10⁻³ AU
    //
    // a_J2 = ratio * [r_x * (5z² - r²), r_y * (5z² - r²), r_z * (5z² - 3r²)]
    // dove ratio = 1.5 * GM0 * J2 * R²_sun / r⁵
    
    double r2 = pos.dot(pos);
    double r = std::sqrt(r2);
    if (r < 1e-15) return Vector3D(0, 0, 0);
    
    double r5 = r2 * r2 * r;
    double ratio = 1.5 * GM0 * SOLJ2 * (RADSUN * RADSUN / r2) / r5;
    
    double z2 = pos.z * pos.z;
    double brac1 = 5.0 * z2 - r2;    // Per x e y
    double brac2 = 5.0 * z2 - 3.0 * r2;  // Per z
    
    return Vector3D(ratio * pos.x * brac1,
                    ratio * pos.y * brac1,
                    ratio * pos.z * brac2);
}

// ============================================================================
// Accelerazione totale (FORCE in OrbFit)
// ============================================================================

Vector3D OrbFitForceModel::computeAcceleration(const Vector3D& pos,
                                                const Vector3D& vel,
                                                double jd) const {
    if (!initialized_) {
        const_cast<OrbFitForceModel*>(this)->initializeEphemeris();
    }
    
    Vector3D f(0, 0, 0);
    
    // 1. Attrazione solare centrale
    f = f + solarAcceleration(pos);
    
    // 2. Perturbazioni planetarie (diretto + indiretto)
    f = f + planetaryAcceleration(pos, jd);
    
    // 3. Perturbazioni asteroidi massicci
    f = f + asteroidAcceleration(pos, jd);
    
    // 4. Correzione relativistica
    if (options_.relativityLevel == 1) {
        // Solo Schwarzschild solare
        f = f + genrel(pos, vel);
    } else if (options_.relativityLevel >= 2) {
        // EIH completo
        f = f + eihrel2(pos, vel, jd);
        
        // J2 del Sole (solo con relatività completa)
        if (options_.includeJ2Sun) {
            f = f + j2sun(pos);
        }
    }
    
    return f;
}

Vector3D OrbFitForceModel::computeAccelerationDetailed(
    const Vector3D& pos,
    const Vector3D& vel,
    double jd,
    Vector3D& sunAccel,
    Vector3D& planetAccel,
    Vector3D& asteroidAccel,
    Vector3D& relativisticAccel) const {
    
    if (!initialized_) {
        const_cast<OrbFitForceModel*>(this)->initializeEphemeris();
    }
    
    sunAccel = solarAcceleration(pos);
    planetAccel = planetaryAcceleration(pos, jd);
    asteroidAccel = asteroidAcceleration(pos, jd);
    
    if (options_.relativityLevel == 1) {
        relativisticAccel = genrel(pos, vel);
    } else if (options_.relativityLevel >= 2) {
        relativisticAccel = eihrel2(pos, vel, jd);
        if (options_.includeJ2Sun) {
            relativisticAccel = relativisticAccel + j2sun(pos);
        }
    } else {
        relativisticAccel = Vector3D(0, 0, 0);
    }
    
    return sunAccel + planetAccel + asteroidAccel + relativisticAccel;
}

// ============================================================================
// Utility
// ============================================================================

void OrbFitForceModel::setOptions(const OrbFitForceOptions& options) {
    options_ = options;
    cachedJD_ = 0;  // Invalida cache
}

PlanetaryState OrbFitForceModel::getPlanetState(int planetIndex, double jd) const {
    computePlanetPositions(jd);
    
    if (planetIndex < 0 || planetIndex >= numPlanets_) {
        throw std::out_of_range("Planet index out of range");
    }
    
    return planetCache_[planetIndex];
}

Vector3D OrbFitForceModel::equatorialToEcliptic(const Vector3D& eq) {
    // Rotazione intorno all'asse X di -obliquità
    // x' = x
    // y' = y*cos(ε) + z*sin(ε)
    // z' = -y*sin(ε) + z*cos(ε)
    return Vector3D(eq.x,
                    eq.y * COS_OBL + eq.z * SIN_OBL,
                   -eq.y * SIN_OBL + eq.z * COS_OBL);
}

Vector3D OrbFitForceModel::eclipticToEquatorial(const Vector3D& ecl) {
    // Rotazione inversa
    return Vector3D(ecl.x,
                    ecl.y * COS_OBL - ecl.z * SIN_OBL,
                    ecl.y * SIN_OBL + ecl.z * COS_OBL);
}

} // namespace ioccultcalc
