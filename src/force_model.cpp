/**
 * @file force_model.cpp
 * @brief Implementazione modello completo di forze gravitazionali N-body
 */

#include "ioccultcalc/force_model.h"
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>

namespace ioccultcalc {

// Conversione AU/day
constexpr double AU_KM = 1.495978707e8;           // 1 AU in km
constexpr double DAY_SEC = 86400.0;               // 1 giorno in secondi
constexpr double AU_DAY_TO_KM_S = AU_KM / DAY_SEC; // conversione velocità

// ============================================================================
// BodyData
// ============================================================================

BodyData::BodyData() 
    : body(PerturbingBody::SUN), name("Unknown"), GM(0), GM_AU_day(0),
      radius(0), mass_ratio(0), isPlanet(false) {}

BodyData::BodyData(PerturbingBody b, const std::string& n, double gm, double rad)
    : body(b), name(n), GM(gm), radius(rad), isPlanet(true) {
    GM_AU_day = convertGM_to_AU_day(gm);
    mass_ratio = gm / PhysicalConstants::GM_SUN;
}

double BodyData::convertGM_to_AU_day(double GM_km_s) {
    // GM in km³/s² -> AU³/day²
    // AU³ = (AU_KM)³ km³
    // day² = (DAY_SEC)² s²
    return GM_km_s * (DAY_SEC * DAY_SEC) / (AU_KM * AU_KM * AU_KM);
}

// ============================================================================
// ForceModelConfig
// ============================================================================

ForceModelConfig::ForceModelConfig()
    : includeSun(true),
      includeMercury(true),
      includeVenus(true),
      includeEarth(true),
      includeMars(true),
      includeJupiter(true),
      includeSaturn(true),
      includeUranus(true),
      includeNeptune(true),
      includeMoon(true),
      includePluto(false),
      includeCeres(false),
      includePallas(false),
      includeVesta(false),
      includeRelativisticCorrection(false),
      minDistanceForPerturbation(0.0),
      useJ2Oblateness(false) {}

ForceModelConfig ForceModelConfig::fastConfig() {
    ForceModelConfig config;
    config.includeMercury = false;
    config.includeVenus = false;
    config.includeMars = false;
    config.includeUranus = false;
    config.includeNeptune = false;
    config.includeMoon = false;
    config.includePluto = false;
    config.includeCeres = false;
    config.includePallas = false;
    config.includeVesta = false;
    config.minDistanceForPerturbation = 1.0; // Ignora se > 1 AU distante
    return config; // Solo Giove, Saturno, Terra
}

ForceModelConfig ForceModelConfig::standardConfig() {
    ForceModelConfig config;
    config.includePluto = false;
    config.includeCeres = false;
    config.includePallas = false;
    config.includeVesta = false;
    return config; // Tutti i pianeti + Luna
}

ForceModelConfig ForceModelConfig::highPrecisionConfig() {
    ForceModelConfig config;
    config.includePluto = true;
    config.includeCeres = true;
    config.includePallas = true;
    config.includeVesta = true;
    config.includeRelativisticCorrection = true;
    return config; // Tutto + relatività
}

ForceModelConfig ForceModelConfig::fullConfig() {
    return highPrecisionConfig();
}

// ============================================================================
// ForceModel
// ============================================================================

ForceModel::ForceModel() : config_(ForceModelConfig::standardConfig()), jplInitialized_(false), spkInitialized_(false) {
    initializeBodyDatabase();
    // JPL viene inizializzato on-demand o tramite initializeJPL()
    // SPK per asteroidi viene caricato se necessario
}

ForceModel::ForceModel(const ForceModelConfig& config) : config_(config), jplInitialized_(false), spkInitialized_(false) {
    initializeBodyDatabase();
    // JPL viene inizializzato on-demand o tramite initializeJPL()
    // SPK per asteroidi viene caricato se necessario
}

void ForceModel::initializeBodyDatabase() {
    using namespace PhysicalConstants;
    
    bodyDatabase_[PerturbingBody::SUN] = 
        BodyData(PerturbingBody::SUN, "Sun", GM_SUN, 695700.0);
    
    bodyDatabase_[PerturbingBody::MERCURY] = 
        BodyData(PerturbingBody::MERCURY, "Mercury", GM_MERCURY, 2439.7);
    
    bodyDatabase_[PerturbingBody::VENUS] = 
        BodyData(PerturbingBody::VENUS, "Venus", GM_VENUS, 6051.8);
    
    bodyDatabase_[PerturbingBody::EARTH] = 
        BodyData(PerturbingBody::EARTH, "Earth", GM_EARTH, 6371.0);
    
    bodyDatabase_[PerturbingBody::MARS] = 
        BodyData(PerturbingBody::MARS, "Mars", GM_MARS, 3389.5);
    
    bodyDatabase_[PerturbingBody::JUPITER] = 
        BodyData(PerturbingBody::JUPITER, "Jupiter", GM_JUPITER, 69911.0);
    
    bodyDatabase_[PerturbingBody::SATURN] = 
        BodyData(PerturbingBody::SATURN, "Saturn", GM_SATURN, 58232.0);
    
    bodyDatabase_[PerturbingBody::URANUS] = 
        BodyData(PerturbingBody::URANUS, "Uranus", GM_URANUS, 25362.0);
    
    bodyDatabase_[PerturbingBody::NEPTUNE] = 
        BodyData(PerturbingBody::NEPTUNE, "Neptune", GM_NEPTUNE, 24622.0);
    
    bodyDatabase_[PerturbingBody::MOON] = 
        BodyData(PerturbingBody::MOON, "Moon", GM_MOON, 1737.4);
    
    bodyDatabase_[PerturbingBody::PLUTO] = 
        BodyData(PerturbingBody::PLUTO, "Pluto", GM_PLUTO, 1188.3);
    
    bodyDatabase_[PerturbingBody::CERES] = 
        BodyData(PerturbingBody::CERES, "(1) Ceres", GM_CERES, 476.2);
    bodyDatabase_[PerturbingBody::CERES].isPlanet = false;
    
    bodyDatabase_[PerturbingBody::PALLAS] = 
        BodyData(PerturbingBody::PALLAS, "(2) Pallas", GM_PALLAS, 272.5);
    bodyDatabase_[PerturbingBody::PALLAS].isPlanet = false;
    
    bodyDatabase_[PerturbingBody::VESTA] = 
        BodyData(PerturbingBody::VESTA, "(4) Vesta", GM_VESTA, 262.7);
    bodyDatabase_[PerturbingBody::VESTA].isPlanet = false;
}

Vector3D ForceModel::computeAcceleration(double jd, 
                                         const Vector3D& pos,
                                         const Vector3D& vel) const {
    Vector3D totalAccel(0, 0, 0);
    
    // Accelerazione centrale del Sole
    double r = pos.magnitude();
    if (r > 0) {
        double GM_sun_AU = bodyDatabase_.at(PerturbingBody::SUN).GM_AU_day;
        double r3 = r * r * r;
        totalAccel = pos * (-GM_sun_AU / r3);
    }
    
    // Perturbazioni planetarie
    auto activeBodies = getActiveBodies();
    for (const auto& body : activeBodies) {
        if (body == PerturbingBody::SUN) continue; // Già incluso
        
        Vector3D accel = computeSingleBodyAcceleration(body, jd, pos);
        totalAccel = totalAccel + accel;
    }
    
    // Correzione relativistica (opzionale)
    if (config_.includeRelativisticCorrection) {
        Vector3D relAccel = computeRelativisticCorrection(pos, vel);
        totalAccel = totalAccel + relAccel;
    }
    
    return totalAccel;
}

Vector3D ForceModel::computeAccelerationWithBreakdown(
    double jd,
    const Vector3D& pos,
    const Vector3D& vel,
    std::map<PerturbingBody, Vector3D>& contributions) const {
    
    contributions.clear();
    Vector3D totalAccel(0, 0, 0);
    
    // Sole (centrale)
    double r = pos.magnitude();
    if (r > 0) {
        double GM_sun_AU = bodyDatabase_.at(PerturbingBody::SUN).GM_AU_day;
        double r3 = r * r * r;
        Vector3D sunAccel = pos * (-GM_sun_AU / r3);
        contributions[PerturbingBody::SUN] = sunAccel;
        totalAccel = totalAccel + sunAccel;
    }
    
    // Perturbatori
    auto activeBodies = getActiveBodies();
    for (const auto& body : activeBodies) {
        if (body == PerturbingBody::SUN) continue;
        
        Vector3D accel = computeSingleBodyAcceleration(body, jd, pos);
        contributions[body] = accel;
        totalAccel = totalAccel + accel;
    }
    
    return totalAccel;
}

Vector3D ForceModel::computeSingleBodyAcceleration(PerturbingBody body,
                                                   double jd,
                                                   const Vector3D& asteroidPos) const {
    // Ottiene posizione del corpo perturbatore
    Vector3D bodyPos = getBodyPosition(body, jd);
    double GM = bodyDatabase_.at(body).GM_AU_day;
    
    // Verifica rilevanza (ottimizzazione)
    if (!isBodyRelevant(asteroidPos, bodyPos, GM)) {
        return Vector3D(0, 0, 0);
    }
    
    // Accelerazione punto-massa (metodo diretto)
    return pointMassAcceleration(asteroidPos, bodyPos, GM);
}

Vector3D ForceModel::pointMassAcceleration(const Vector3D& asteroidPos,
                                          const Vector3D& bodyPos,
                                          double GM) {
    // a = GM * [(r_body - r_ast) / |r_body - r_ast|³ - r_body / |r_body|³]
    // Primo termine: attrazione diretta
    // Secondo termine: correzione riferimento non inerziale (eliocentrico)
    
    Vector3D delta = bodyPos - asteroidPos;
    double delta_mag = delta.magnitude();
    
    if (delta_mag < 1e-10) return Vector3D(0, 0, 0); // Evita singolarità
    
    double delta3 = delta_mag * delta_mag * delta_mag;
    Vector3D directTerm = delta * (GM / delta3);
    
    // Termine indiretto (riferimento eliocentrico)
    double bodyDist = bodyPos.magnitude();
    if (bodyDist < 1e-10) return directTerm;
    
    double bodyDist3 = bodyDist * bodyDist * bodyDist;
    Vector3D indirectTerm = bodyPos * (GM / bodyDist3);
    
    return directTerm - indirectTerm;
}

// Converte PerturbingBody in JPLBody
static JPLBody toJPLBody(PerturbingBody body) {
    switch (body) {
        case PerturbingBody::MERCURY: return JPLBody::MERCURY;
        case PerturbingBody::VENUS: return JPLBody::VENUS;
        case PerturbingBody::EARTH: return JPLBody::EARTH;
        case PerturbingBody::MARS: return JPLBody::MARS;
        case PerturbingBody::JUPITER: return JPLBody::JUPITER;
        case PerturbingBody::SATURN: return JPLBody::SATURN;
        case PerturbingBody::URANUS: return JPLBody::URANUS;
        case PerturbingBody::NEPTUNE: return JPLBody::NEPTUNE;
        case PerturbingBody::MOON: return JPLBody::MOON;
        case PerturbingBody::PLUTO: return JPLBody::PLUTO;
        case PerturbingBody::SUN: return JPLBody::SUN;
        default: return JPLBody::SUN; // Fallback
    }
}

bool ForceModel::initializeJPL(JPLVersion version, const std::string& filepath) {
    if (jplInitialized_) {
        return true; // Già inizializzato
    }
    
    // TODO: JPL ephemeris not yet implemented - need to implement jpl_ephemeris.cpp first
    // jplReader_ = JPLEphemerisReader(version);
    (void)version; // Suppress unused warning
    
    bool success;
    if (!filepath.empty()) {
        // Carica da file specifico
        success = jplReader_.loadFile(filepath);
    } else {
        // Auto-download se necessario
        success = jplReader_.downloadDE(version);
    }
    
    if (success) {
        jplInitialized_ = true;
        
        // Aggiorna parametri GM dal file JPL
        JPLConstants constants = jplReader_.getConstants();
        bodyDatabase_[PerturbingBody::SUN].GM = constants.GM_Sun;
        bodyDatabase_[PerturbingBody::EARTH].GM = constants.GM_Earth;
        bodyDatabase_[PerturbingBody::MOON].GM = constants.GM_Moon;
        
        for (int i = 0; i < 9; ++i) {
            PerturbingBody body;
            switch (i) {
                case 0: body = PerturbingBody::MERCURY; break;
                case 1: body = PerturbingBody::VENUS; break;
                case 2: body = PerturbingBody::EARTH; break;
                case 3: body = PerturbingBody::MARS; break;
                case 4: body = PerturbingBody::JUPITER; break;
                case 5: body = PerturbingBody::SATURN; break;
                case 6: body = PerturbingBody::URANUS; break;
                case 7: body = PerturbingBody::NEPTUNE; break;
                case 8: body = PerturbingBody::PLUTO; break;
                default: continue;
            }
            bodyDatabase_[body].GM = constants.GM_planets[i];
            bodyDatabase_[body].GM_AU_day = BodyData::convertGM_to_AU_day(constants.GM_planets[i]);
        }
    }
    
    return success;
}

JPLVersion ForceModel::getJPLVersion() const {
    if (!jplInitialized_) {
        return JPLVersion::DE441; // Default
    }
    return jplReader_.getVersion();
}

Vector3D ForceModel::getBodyPosition(PerturbingBody body, double jd) const {
    // Verifica cache
    auto cacheKey = std::make_pair(body, jd);
    auto it = positionCache_.find(cacheKey);
    if (it != positionCache_.end()) {
        return it->second;
    }
    
    // Inizializza JPL se necessario (lazy initialization)
    if (!jplInitialized_) {
        const_cast<ForceModel*>(this)->initializeJPL();
    }
    
    Vector3D pos;
    
    if (body == PerturbingBody::SUN) {
        // Sole: otteniamo posizione baricentrica e convertiamo in eliocentrica
        // Per riferimento eliocentrico, Sole è all'origine
        pos = Vector3D(0, 0, 0);
    } else if (body == PerturbingBody::CERES || 
               body == PerturbingBody::PALLAS || 
               body == PerturbingBody::VESTA) {
        // Grandi asteroidi: usa SPK (SB441-N16.bsp)
        if (!spkInitialized_) {
            // Lazy initialization: carica SB441-N16.bsp
            std::string spkPath = "sb441-n16.bsp";  // Cercare in data directory
            const_cast<ForceModel*>(this)->spkInitialized_ = 
                const_cast<ForceModel*>(this)->spkReader_.ensureFileLoaded(spkPath);
        }
        
        if (spkInitialized_) {
            // NAIF IDs: Ceres=2000001, Pallas=2000002, Vesta=2000004
            int naifId = 0;
            if (body == PerturbingBody::CERES) naifId = 2000001;
            else if (body == PerturbingBody::PALLAS) naifId = 2000002;
            else if (body == PerturbingBody::VESTA) naifId = 2000004;
            
            try {
                // getPosition ritorna posizione eliocentrica (centerId=10)
                pos = spkReader_.getPosition(naifId, jd, 10);
            } catch (...) {
                // Se fallisce, restituisci origine (zero perturbation)
                pos = Vector3D(0, 0, 0);
            }
        } else {
            pos = Vector3D(0, 0, 0);
        }
    } else {
        // Pianeti e Luna: usa JPL DE
        JPLBody jplBody = toJPLBody(body);
        
        // Ottiene posizione barycentric da JPL (in km)
        Vector3D posKm = jplReader_.getPosition(jplBody, jd);
        
        // Converte in eliocentrico sottraendo posizione Sole
        Vector3D sunPosKm = jplReader_.getPosition(JPLBody::SUN, jd);
        Vector3D heliocentricKm = posKm - sunPosKm;
        
        // Converte da km ad AU
        constexpr double KM_TO_AU = 1.0 / AU_KM;
        pos = heliocentricKm * KM_TO_AU;
    }
    
    // Salva in cache
    positionCache_[cacheKey] = pos;
    
    return pos;
}

Vector3D ForceModel::getBodyVelocity(PerturbingBody body, double jd) const {
    // Verifica cache
    auto cacheKey = std::make_pair(body, jd);
    auto it = velocityCache_.find(cacheKey);
    if (it != velocityCache_.end()) {
        return it->second;
    }
    
    // Inizializza JPL se necessario
    if (!jplInitialized_) {
        const_cast<ForceModel*>(this)->initializeJPL();
    }
    
    Vector3D vel;
    
    if (body == PerturbingBody::SUN) {
        // Sole: velocità zero in riferimento eliocentrico
        vel = Vector3D(0, 0, 0);
    } else if (body == PerturbingBody::CERES || 
               body == PerturbingBody::PALLAS || 
               body == PerturbingBody::VESTA) {
        // Asteroidi: usa differenze finite (per ora)
        double dt = 0.001; // 0.001 giorni
        Vector3D pos_plus = getBodyPosition(body, jd + dt);
        Vector3D pos_minus = getBodyPosition(body, jd - dt);
        vel = (pos_plus - pos_minus) * (0.5 / dt);
    } else {
        // Pianeti: JPL fornisce velocità direttamente
        JPLBody jplBody = toJPLBody(body);
        
        // Ottiene velocità barycentric da JPL (in km/day)
        Vector3D velKm = jplReader_.getVelocity(jplBody, jd);
        
        // Converte in eliocentrico sottraendo velocità Sole
        Vector3D sunVelKm = jplReader_.getVelocity(JPLBody::SUN, jd);
        Vector3D heliocentricVelKm = velKm - sunVelKm;
        
        // Converte da km/day ad AU/day
        constexpr double KM_TO_AU = 1.0 / AU_KM;
        vel = heliocentricVelKm * KM_TO_AU;
    }
    
    // Salva in cache
    velocityCache_[cacheKey] = vel;
    
    return vel;
}

BodyData ForceModel::getBodyData(PerturbingBody body) const {
    auto it = bodyDatabase_.find(body);
    if (it != bodyDatabase_.end()) {
        return it->second;
    }
    return BodyData();
}

void ForceModel::setConfiguration(const ForceModelConfig& config) {
    config_ = config;
    clearCache(); // Invalida cache con nuova configurazione
}

double ForceModel::estimatePerturbationMagnitude(PerturbingBody body,
                                                 double jd,
                                                 const Vector3D& asteroidPos) const {
    Vector3D accel = computeSingleBodyAcceleration(body, jd, asteroidPos);
    return accel.magnitude();
}

std::vector<PerturbingBody> ForceModel::getActiveBodies() const {
    std::vector<PerturbingBody> bodies;
    
    if (config_.includeSun) bodies.push_back(PerturbingBody::SUN);
    if (config_.includeMercury) bodies.push_back(PerturbingBody::MERCURY);
    if (config_.includeVenus) bodies.push_back(PerturbingBody::VENUS);
    if (config_.includeEarth) bodies.push_back(PerturbingBody::EARTH);
    if (config_.includeMars) bodies.push_back(PerturbingBody::MARS);
    if (config_.includeJupiter) bodies.push_back(PerturbingBody::JUPITER);
    if (config_.includeSaturn) bodies.push_back(PerturbingBody::SATURN);
    if (config_.includeUranus) bodies.push_back(PerturbingBody::URANUS);
    if (config_.includeNeptune) bodies.push_back(PerturbingBody::NEPTUNE);
    if (config_.includeMoon) bodies.push_back(PerturbingBody::MOON);
    if (config_.includePluto) bodies.push_back(PerturbingBody::PLUTO);
    if (config_.includeCeres) bodies.push_back(PerturbingBody::CERES);
    if (config_.includePallas) bodies.push_back(PerturbingBody::PALLAS);
    if (config_.includeVesta) bodies.push_back(PerturbingBody::VESTA);
    
    return bodies;
}

Vector3D ForceModel::computeRelativisticCorrection(const Vector3D& pos,
                                                   const Vector3D& vel) const {
    // Correzione post-Newtoniana (Schwarzschild, ordine PN1)
    // a_rel = (GM/c²r³) * [4(GM/r) - v² ] * r + 4(r·v) * v
    
    using namespace PhysicalConstants;
    
    double r = pos.magnitude();
    if (r < 1e-10) return Vector3D(0, 0, 0);
    
    double v2 = vel.dot(vel);
    double rdotv = pos.dot(vel);
    
    // Converte unità: GM in AU³/day², c in AU/day
    double GM_sun_AU = bodyDatabase_.at(PerturbingBody::SUN).GM_AU_day;
    double c_AU_day = (SPEED_OF_LIGHT * DAY_SEC) / AU_KM; // ~173.14 AU/day
    double c2 = c_AU_day * c_AU_day;
    
    double r3 = r * r * r;
    double factor = GM_sun_AU / (c2 * r3);
    
    double coeff1 = 4.0 * GM_sun_AU / r - v2;
    Vector3D term1 = pos * (factor * coeff1);
    
    Vector3D term2 = vel * (4.0 * factor * rdotv);
    
    return term1 + term2;
}

bool ForceModel::isCacheInitialized() const {
    return !positionCache_.empty();
}

void ForceModel::precacheEphemerides(double jdStart, double jdEnd, double stepDays) {
    clearCache();
    
    auto bodies = getActiveBodies();
    
    for (double jd = jdStart; jd <= jdEnd; jd += stepDays) {
        for (const auto& body : bodies) {
            if (body == PerturbingBody::SUN) continue;
            
            // Forza calcolo e cache
            getBodyPosition(body, jd);
            getBodyVelocity(body, jd);
        }
    }
}

void ForceModel::clearCache() {
    positionCache_.clear();
    velocityCache_.clear();
}

bool ForceModel::isBodyRelevant(const Vector3D& asteroidPos,
                               const Vector3D& bodyPos,
                               double GM) const {
    if (config_.minDistanceForPerturbation <= 0) {
        return true; // Sempre includi
    }
    
    Vector3D delta = bodyPos - asteroidPos;
    double distance = delta.magnitude();
    
    // Stima accelerazione
    double accel_est = GM / (distance * distance);
    
    // Soglia: 1e-15 AU/day² (circa 0.01 km/anno²)
    constexpr double MIN_ACCEL = 1e-15;
    
    return accel_est > MIN_ACCEL;
}

// ============================================================================
// ForceModelAnalyzer
// ============================================================================

std::vector<ForceModelAnalyzer::PerturbationAnalysis> 
ForceModelAnalyzer::analyze(const ForceModel& model,
                            double jd,
                            const Vector3D& pos,
                            const Vector3D& vel) {
    std::map<PerturbingBody, Vector3D> contributions;
    model.computeAccelerationWithBreakdown(jd, pos, vel, contributions);
    
    // Calcola magnitudine totale perturbazioni (escluso Sole)
    double totalPerturbMag = 0;
    for (const auto& pair : contributions) {
        if (pair.first != PerturbingBody::SUN) {
            totalPerturbMag += pair.second.magnitude();
        }
    }
    
    std::vector<PerturbationAnalysis> results;
    
    for (const auto& pair : contributions) {
        PerturbationAnalysis analysis;
        analysis.body = pair.first;
        analysis.acceleration = pair.second;
        analysis.magnitude = pair.second.magnitude();
        
        if (totalPerturbMag > 0 && pair.first != PerturbingBody::SUN) {
            analysis.relativeContribution = 
                100.0 * analysis.magnitude / totalPerturbMag;
        } else {
            analysis.relativeContribution = 0;
        }
        
        Vector3D bodyPos = model.getBodyPosition(pair.first, jd);
        Vector3D delta = bodyPos - pos;
        analysis.distance = delta.magnitude();
        
        results.push_back(analysis);
    }
    
    // Ordina per magnitudine decrescente
    std::sort(results.begin(), results.end(),
              [](const PerturbationAnalysis& a, const PerturbationAnalysis& b) {
                  return a.magnitude > b.magnitude;
              });
    
    return results;
}

std::string ForceModelAnalyzer::generateReport(
    const std::vector<PerturbationAnalysis>& analysis) {
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    
    oss << "=== Force Model Perturbation Analysis ===\n\n";
    oss << std::setw(15) << "Body" 
        << std::setw(18) << "Accel (AU/day²)"
        << std::setw(12) << "Contrib %"
        << std::setw(15) << "Distance (AU)\n";
    oss << std::string(60, '-') << "\n";
    
    for (const auto& item : analysis) {
        std::string bodyName = "Unknown";
        switch (item.body) {
            case PerturbingBody::SUN: bodyName = "Sun"; break;
            case PerturbingBody::MERCURY: bodyName = "Mercury"; break;
            case PerturbingBody::VENUS: bodyName = "Venus"; break;
            case PerturbingBody::EARTH: bodyName = "Earth"; break;
            case PerturbingBody::MARS: bodyName = "Mars"; break;
            case PerturbingBody::JUPITER: bodyName = "Jupiter"; break;
            case PerturbingBody::SATURN: bodyName = "Saturn"; break;
            case PerturbingBody::URANUS: bodyName = "Uranus"; break;
            case PerturbingBody::NEPTUNE: bodyName = "Neptune"; break;
            case PerturbingBody::MOON: bodyName = "Moon"; break;
            case PerturbingBody::PLUTO: bodyName = "Pluto"; break;
            case PerturbingBody::CERES: bodyName = "Ceres"; break;
            case PerturbingBody::PALLAS: bodyName = "Pallas"; break;
            case PerturbingBody::VESTA: bodyName = "Vesta"; break;
        }
        
        oss << std::setw(15) << bodyName
            << std::setw(18) << std::scientific << item.magnitude
            << std::setw(12) << std::fixed << item.relativeContribution
            << std::setw(15) << item.distance << "\n";
    }
    
    return oss.str();
}

double ForceModelAnalyzer::estimateOmissionError(PerturbingBody body,
                                                double jd,
                                                const Vector3D& pos,
                                                double propagationDays) {
    // Stima errore = 0.5 * accel * t²
    // Formula semplificata per ordine di grandezza
    
    ForceModel model;
    double accelMag = model.estimatePerturbationMagnitude(body, jd, pos);
    
    // Converte in km
    double t_days = propagationDays;
    double error_AU = 0.5 * accelMag * t_days * t_days;
    double error_km = error_AU * AU_KM;
    
    return error_km;
}

} // namespace ioccultcalc
