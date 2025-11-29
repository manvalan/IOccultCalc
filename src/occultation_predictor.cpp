#include "ioccultcalc/occultation_predictor.h"
#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/coordinates.h"
#include "ioccultcalc/time_utils.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace ioccultcalc {

class OccultationPredictor::Impl {
public:
    EquinoctialElements asteroid;
    Ephemeris ephemeris;
    GaiaClient gaiaClient;
    AstDysClient astdysClient;
    
    double asteroidDiameter;    // km
    double orbitalUncertainty;  // km (1-sigma)
    
    // Chebyshev approximation support
    ChebyshevConfig chebyshevConfig;
    std::unique_ptr<ChebyshevCache> chebyshevCache;
    
    // Advanced options (opzionali, default false per retrocompatibilità)
    bool applyParallaxCorrection;   // Correzione parallasse stelle vicine
    double angularBufferArcmin;     // Buffer ricerca (arcmin)
    bool useAdaptiveTimestep;       // Timestep adattivo
    double coarseStepDays;          // Step grossolano per adaptive
    double fineStepMinutes;         // Step fine per adaptive
    
    Impl() : asteroidDiameter(0), orbitalUncertainty(100.0),
             applyParallaxCorrection(false), angularBufferArcmin(0.0),
             useAdaptiveTimestep(false), coarseStepDays(1.0), 
             fineStepMinutes(15.0) {
        // Default: Chebyshev disabilitato
        chebyshevConfig.enabled = false;
    }
};

OccultationPredictor::OccultationPredictor() : pImpl(new Impl()) {}

OccultationPredictor::~OccultationPredictor() = default;

void OccultationPredictor::setAsteroid(const EquinoctialElements& elements) {
    pImpl->asteroid = elements;
    pImpl->ephemeris.setElements(elements);
    
    if (elements.diameter > 0) {
        pImpl->asteroidDiameter = elements.diameter;
    }
}

void OccultationPredictor::loadAsteroidFromAstDyS(const std::string& designation) {
    EquinoctialElements elem = pImpl->astdysClient.getElements(designation);
    setAsteroid(elem);
}

void OccultationPredictor::setLocalAstDySDirectories(const std::string& eq1Dir, const std::string& rwoDir) {
    pImpl->astdysClient.setLocalEQ1Directory(eq1Dir);
    pImpl->astdysClient.setLocalRWODirectory(rwoDir);
}

void OccultationPredictor::setAsteroidDiameter(double diameter) {
    pImpl->asteroidDiameter = diameter;
}

void OccultationPredictor::setOrbitalUncertainty(double sigmaKm) {
    pImpl->orbitalUncertainty = sigmaKm;
}

void OccultationPredictor::setChebyshevConfig(const ChebyshevConfig& config) {
    pImpl->chebyshevConfig = config;
}

const ChebyshevConfig& OccultationPredictor::getChebyshevConfig() const {
    return pImpl->chebyshevConfig;
}

void OccultationPredictor::setParallaxCorrection(bool enabled) {
    pImpl->applyParallaxCorrection = enabled;
    if (enabled) {
        std::cout << "Parallax correction ENABLED (accurate for stars < 100 pc)" << std::endl;
    }
}

void OccultationPredictor::setAngularBuffer(double arcmin) {
    pImpl->angularBufferArcmin = arcmin;
    if (arcmin > 0) {
        std::cout << "Angular buffer set to " << arcmin << " arcmin" << std::endl;
    }
}

void OccultationPredictor::setAdaptiveTimestep(bool enabled, double coarseStepDays, 
                                               double fineStepMinutes) {
    pImpl->useAdaptiveTimestep = enabled;
    pImpl->coarseStepDays = coarseStepDays;
    pImpl->fineStepMinutes = fineStepMinutes;
    if (enabled) {
        std::cout << "Adaptive timestep ENABLED (coarse: " << coarseStepDays 
                  << " days, fine: " << fineStepMinutes << " min)" << std::endl;
    }
}

std::vector<OccultationEvent> OccultationPredictor::findOccultations(
    const JulianDate& startJD,
    const JulianDate& endJD,
    double maxMagnitude,
    double searchRadius,
    double minProbability) {
    
    std::vector<OccultationEvent> events;
    
    // Determina strategia computazionale
    bool useChebyshev = pImpl->chebyshevConfig.enabled;
    bool useAdaptive = pImpl->useAdaptiveTimestep;
    
    std::vector<EphemerisData> ephemerides;
    
    // ============================================================
    // FASE 1: Calcola effemeridi asteroide per tutto il periodo
    // ============================================================
    
    if (useChebyshev) {
        std::cout << "Using Chebyshev approximation for fast ephemeris computation" << std::endl;
        
        // Inizializza cache Chebyshev
        pImpl->chebyshevCache = std::make_unique<ChebyshevCache>(pImpl->chebyshevConfig);
        pImpl->chebyshevCache->initialize(pImpl->ephemeris, startJD, endJD);
        
        // Genera timesteps fini
        auto timesteps = pImpl->chebyshevCache->generateFineTimesteps();
        std::cout << "Generated " << timesteps.size() << " fine timesteps" << std::endl;
        
        // Valuta effemeridi usando Chebyshev
        for (const auto& jd : timesteps) {
            ephemerides.push_back(pImpl->chebyshevCache->evaluate(jd));
        }
        
    } else {
        // Metodo standard: integrazione RADAU timestep fisso
        double stepDays = useAdaptive ? pImpl->coarseStepDays : 1.0;
        ephemerides = pImpl->ephemeris.computeRange(startJD, endJD, stepDays);
    }
    
    std::cout << "Computed " << ephemerides.size() << " ephemeris points" << std::endl;
    
    // DEBUG: Stampa posizioni asteroide
    std::cout << "\nDEBUG: Asteroid path positions:\n";
    for (size_t i = 0; i < ephemerides.size(); i++) {
        const auto& eph = ephemerides[i];
        std::cout << "  [" << i << "] JD=" << std::fixed << std::setprecision(2) << eph.jd.jd 
                  << " RA=" << std::setprecision(4) << (eph.geocentricPos.ra * RAD_TO_DEG)
                  << "° Dec=" << (eph.geocentricPos.dec * RAD_TO_DEG) << "°\n";
    }
    std::cout << std::endl;
    
    // ============================================================
    // FASE 2: OTTIMIZZAZIONE CORRIDOR - Query stelle una sola volta
    // ============================================================
    
    std::vector<GaiaStar> candidateStars;
    
    // Costruisci percorso asteroide nel cielo
    std::vector<EquatorialCoordinates> asteroidPath;
    asteroidPath.reserve(ephemerides.size());
    
    for (const auto& eph : ephemerides) {
        EquatorialCoordinates pos;
        pos.ra = eph.geocentricPos.ra * RAD_TO_DEG;
        pos.dec = eph.geocentricPos.dec * RAD_TO_DEG;
        asteroidPath.push_back(pos);
    }
    
    // Semplifica il percorso se ha troppi punti (mantieni 1 punto ogni N)
    std::vector<EquatorialCoordinates> simplifiedPath;
    size_t pathStep = std::max(size_t(1), asteroidPath.size() / 50);  // Max ~50 punti
    for (size_t i = 0; i < asteroidPath.size(); i += pathStep) {
        simplifiedPath.push_back(asteroidPath[i]);
    }
    // Assicura che l'ultimo punto sia incluso
    if (simplifiedPath.back().ra != asteroidPath.back().ra ||
        simplifiedPath.back().dec != asteroidPath.back().dec) {
        simplifiedPath.push_back(asteroidPath.back());
    }
    
    if (simplifiedPath.size() >= 2) {
        // Usa query corridor per trovare TUTTE le stelle candidate in un colpo solo
        double corridorWidth = searchRadius * 2.0;  // Buffer per incertezze
        
        std::cout << "Corridor query: " << simplifiedPath.size() << " path points, width=" 
                  << corridorWidth << "°, mag<" << maxMagnitude << std::endl;
        
        candidateStars = pImpl->gaiaClient.queryAlongPath(simplifiedPath, corridorWidth, maxMagnitude);
        
        std::cout << "Found " << candidateStars.size() << " candidate stars along path" << std::endl;
    } else {
        // Fallback a singola query cone (periodo molto breve)
        double raDeg = ephemerides[0].geocentricPos.ra * RAD_TO_DEG;
        double decDeg = ephemerides[0].geocentricPos.dec * RAD_TO_DEG;
        candidateStars = pImpl->gaiaClient.queryCone(raDeg, decDeg, searchRadius * 3.0, maxMagnitude);
        std::cout << "Single cone query: " << candidateStars.size() << " candidate stars" << std::endl;
    }
    
    if (candidateStars.empty()) {
        std::cout << "No candidate stars found - returning empty" << std::endl;
        return events;
    }
    
    // ============================================================
    // FASE 3: Per ogni stella candidata, trova closest approach
    // ============================================================
    
    std::cout << "Analyzing " << candidateStars.size() << " candidate stars..." << std::endl;
    
    // DEBUG: Trova la stella target e stampa la separazione
    const std::string TARGET_STAR = "3411546266140512128";
    bool targetFound = false;
    
    for (const auto& star : candidateStars) {
        if (star.sourceId == TARGET_STAR) {
            std::cout << "\n*** DEBUG: TARGET STAR FOUND IN CORRIDOR! ***\n";
            std::cout << "    RA=" << star.pos.ra << "° Dec=" << star.pos.dec << "° G=" << star.phot_g_mean_mag << "\n";
            targetFound = true;
        }
        
        // Trova l'epoca di minima separazione scansionando le effemeridi
        double minSeparation = 1e9;
        size_t minIndex = 0;
        
        for (size_t i = 0; i < ephemerides.size(); i++) {
            const auto& eph = ephemerides[i];
            
            // Propaga stella all'epoca
            auto starPos = star.propagateToEpoch(eph.jd);
            
            // Calcola separazione angolare
            double separation = Coordinates::angularSeparation(eph.geocentricPos, starPos);
            
            if (separation < minSeparation) {
                minSeparation = separation;
                minIndex = i;
            }
        }
        
        double separationArcsec = minSeparation * RAD_TO_DEG * 3600.0;
        
        // DEBUG: Mostra stella target
        if (star.sourceId == TARGET_STAR) {
            std::cout << "    Min separation at epoch " << minIndex << ": " << separationArcsec << " arcsec\n";
            std::cout << "    Asteroid pos: RA=" << (ephemerides[minIndex].geocentricPos.ra * RAD_TO_DEG) 
                      << "° Dec=" << (ephemerides[minIndex].geocentricPos.dec * RAD_TO_DEG) << "°\n";
        }
        
        // Calcola threshold usando metodo LinOccult (geometria 3D)
        const auto& ephMin = ephemerides[minIndex];
        double asteroidDistanceKm = ephMin.distance * AU;
        
        const double R_EARTH_KM = 6371.0;
        const double MAX_SHADOW_RADIUS_KM = 3.0 * R_EARTH_KM;
        
        double uncertaintyKm = pImpl->orbitalUncertainty;
        double maxDistKm = uncertaintyKm + R_EARTH_KM;
        if (maxDistKm > MAX_SHADOW_RADIUS_KM) {
            maxDistKm = MAX_SHADOW_RADIUS_KM;
        }
        
        double thresholdArcsec = (maxDistKm / asteroidDistanceKm) * RAD_TO_DEG * 3600.0;
        if (thresholdArcsec < 10.0) thresholdArcsec = 10.0;
        
        // Aggiungi buffer angolare opzionale
        double effectiveThreshold = thresholdArcsec;
        if (pImpl->angularBufferArcmin > 0) {
            effectiveThreshold += pImpl->angularBufferArcmin * 60.0;
        }
        
        if (separationArcsec < effectiveThreshold) {
            // Refine closest approach con ricerca precisa
            JulianDate searchStart(ephemerides[std::max(0, (int)minIndex - 1)].jd.jd);
            JulianDate searchEnd(ephemerides[std::min(ephemerides.size()-1, minIndex + 1)].jd.jd);
            
            JulianDate caTime = findClosestApproach(star, searchStart, searchEnd);
            
            // Predici occultazione dettagliata
            OccultationEvent event = predictOccultation(star, caTime);
            
            if (event.probability >= minProbability) {
                events.push_back(event);
            }
        }
    }
    
    // Ordina per tempo
    std::sort(events.begin(), events.end(),
             [](const OccultationEvent& a, const OccultationEvent& b) {
                 return a.timeCA.jd < b.timeCA.jd;
             });
    
    std::cout << "Found " << events.size() << " occultation events" << std::endl;
    
    return events;
}

OccultationEvent OccultationPredictor::predictOccultation(
    const GaiaStar& star,
    const JulianDate& approximateTime) {
    
    OccultationEvent event;
    event.asteroid = pImpl->asteroid;
    event.star = star;
    
    // Trova closest approach preciso
    event.timeCA = findClosestApproach(star,
        JulianDate(approximateTime.jd - 0.1),
        JulianDate(approximateTime.jd + 0.1));
    
    // Calcola effemeridi al momento CA
    EphemerisData ephCA = pImpl->ephemeris.compute(event.timeCA);
    
    // Propaga posizione stella con correzione oblateness Terra
    auto starPosCA = star.propagateToEpoch(event.timeCA);
    
    // Applica correzione parallasse se abilitata (stelle vicine)
    if (pImpl->applyParallaxCorrection && star.parallax > 0.001) {
        // Parallasse > 1 mas → stella entro ~1000 pc
        // Correzione significativa per parallax > 10 mas (< 100 pc)
        
        // Posizione Terra barycenter (AU)
        Vector3D earthPos = Ephemeris::getEarthPosition(event.timeCA);
        
        // Distanza stella in AU
        double distance_pc = 1000.0 / star.parallax;  // mas → pc
        double distance_au = distance_pc * 206265.0;  // pc → AU
        
        // Correzione parallasse (effetto annuo)
        // ΔRA = -x_Earth / distance * sec(Dec)
        // ΔDec = -y_Earth / distance
        double dra = -earthPos.x / distance_au / cos(starPosCA.dec * DEG_TO_RAD);
        double ddec = -earthPos.y / distance_au;
        
        starPosCA.ra += dra * RAD_TO_DEG;
        starPosCA.dec += ddec * RAD_TO_DEG;
    }
    
    // Calcola geometria
    double separation, posAngle;
    calculateGeometry(ephCA, star, event.timeCA, separation, posAngle);
    
    event.closeApproachDistance = separation * RAD_TO_DEG * 3600.0; // arcsec
    event.positionAngle = posAngle * RAD_TO_DEG;
    
    // Dimensione angolare asteroide
    double asteroidAngularSize = 0;
    if (pImpl->asteroidDiameter > 0) {
        asteroidAngularSize = (pImpl->asteroidDiameter / (ephCA.distance * AU)) * RAD_TO_DEG * 3600.0;
    }
    
    // Calcola incertezza totale (asteroide + stella) - LinOccult method
    // Incertezza asteroide
    double asteroidUncertaintyArcsec = (pImpl->orbitalUncertainty / (ephCA.distance * AU)) * RAD_TO_DEG * 3600.0;
    
    // Incertezza stella dipendente da magnitudine (LinOccult CalculateTotalUncertainty)
    const double STAR_UNCERT_FAINT = 60.0 / 1000.0;   // 60 mas per Mv >= 9
    const double STAR_UNCERT_BRIGHT = 7.0 / 1000.0;   // 7 mas per Mv < 9
    const double PROPER_MOTION_UNCERT = 2.5 / 1000.0; // 2.5 mas/anno
    const double MJD_J2015 = 2457023.5;               // J2015.0 epoch
    
    double starMag = star.phot_g_mean_mag;
    double starUncertMas = (starMag < 9.0) ? STAR_UNCERT_BRIGHT : STAR_UNCERT_FAINT;
    
    // Aggiunge degradazione moto proprio dal 2015
    double yearsSince2015 = (event.timeCA.jd - MJD_J2015) / 365.25;
    starUncertMas += fabs(yearsSince2015) * PROPER_MOTION_UNCERT;
    
    double starUncertaintyArcsec = starUncertMas / 1000.0;  // mas → arcsec
    
    // Incertezza totale (somma quadratica)
    double totalUncertaintyArcsec = sqrt(asteroidUncertaintyArcsec * asteroidUncertaintyArcsec + 
                                         starUncertaintyArcsec * starUncertaintyArcsec);
    
    // Calcola probabilità con metodo LinOccult migliorato
    event.probability = calculateProbability(event.closeApproachDistance, 
                                            asteroidAngularSize,
                                            totalUncertaintyArcsec);
    
    // Durata massima
    if (pImpl->asteroidDiameter > 0) {
        // Velocità relativa asteroide-stella in km/s
        Vector3D relVel = ephCA.heliocentricVel - Ephemeris::getEarthVelocity(event.timeCA);
        double relVelKmS = relVel.magnitude() * AU / 86400.0;
        
        event.maxDuration = pImpl->asteroidDiameter / relVelKmS;
    }
    
    // Calcola shadow path
    event.shadowPath = calculateShadowPath(ephCA, star, event.timeCA, 120.0);
    
    // Incertezze (1-sigma)
    event.uncertaintyNorth = pImpl->orbitalUncertainty;
    event.uncertaintySouth = pImpl->orbitalUncertainty;
    
    // Genera ID evento
    event.eventId = pImpl->asteroid.designation + "_" + star.sourceId + "_" +
                   TimeUtils::jdToISO(event.timeCA);
    
    return event;
}

std::vector<ShadowPathPoint> OccultationPredictor::calculateShadowPath(
    const EphemerisData& asteroidPos,
    const GaiaStar& star,
    const JulianDate& centralTime,
    double timeSpanMinutes) {
    
    std::vector<ShadowPathPoint> path;
    
    double halfSpan = timeSpanMinutes / 2.0 / 1440.0; // Converti in giorni
    double stepMinutes = 1.0; // Step di 1 minuto
    double stepDays = stepMinutes / 1440.0;
    
    JulianDate startTime(centralTime.jd - halfSpan);
    JulianDate endTime(centralTime.jd + halfSpan);
    
    for (double jd = startTime.jd; jd <= endTime.jd; jd += stepDays) {
        JulianDate time(jd);
        
        // Effemeridi asteroide
        EphemerisData eph = pImpl->ephemeris.compute(time);
        
        // Posizione stella
        auto starPos = star.propagateToEpoch(time);
        
        // Vettore dalla Terra all'asteroide
        Vector3D asteroidVec = Coordinates::equatorialToCartesian(eph.geocentricPos);
        asteroidVec = asteroidVec * (eph.distance * AU); // Converti in km
        
        // Vettore dalla Terra alla stella (assumendo stella all'infinito)
        Vector3D starDir = Coordinates::equatorialToCartesian(starPos);
        starDir = starDir.normalize();
        
        // Piano fondamentale: perpendicolare alla linea di vista della stella
        // Trova il punto sulla Terra dove l'ombra interseca
        
        // Risolvi per il punto sulla superficie terrestre
        // La posizione del punto è dove la linea Terra-asteroide proiettata
        // lungo la direzione della stella interseca la superficie terrestre
        
        // Approssimazione: assumiamo che l'ombra viaggi lungo starDir
        // Punto ombra = posizione asteroide - t * starDir, dove t è tale che
        // il punto risultante è sulla superficie terrestre
        
        // Semplificazione: proiettiamo il vettore asteroide sul piano geocentrico
        double t = asteroidVec.dot(starDir);
        Vector3D shadowVec = asteroidVec - starDir * t;
        
        // Normalizza alla superficie terrestre
        double shadowDist = shadowVec.magnitude();
        if (shadowDist > EARTH_RADIUS * 0.1) { // Solo se non troppo vicino al centro
            Vector3D earthSurfacePoint = shadowVec.normalize() * EARTH_RADIUS;
            
            // Converti in coordinate geografiche
            GeographicCoordinates geoCoord = Coordinates::ecefToGeographic(earthSurfacePoint);
            
            ShadowPathPoint point;
            point.time = time;
            point.location = geoCoord;
            
            // Durata (già calcolata nell'evento)
            if (pImpl->asteroidDiameter > 0) {
                Vector3D relVel = eph.heliocentricVel - Ephemeris::getEarthVelocity(time);
                double relVelKmS = relVel.magnitude() * AU / 86400.0;
                point.duration = pImpl->asteroidDiameter / relVelKmS;
            }
            
            // Distanza dalla centerline (sempre 0 per questa semplificazione)
            point.centerlineDistance = 0;
            
            path.push_back(point);
        }
    }
    
    return path;
}

JulianDate OccultationPredictor::findClosestApproach(const GaiaStar& star,
                                                      const JulianDate& startJD,
                                                      const JulianDate& endJD) {
    // Ricerca binaria per trovare il momento di minima separazione
    
    auto calcSeparation = [&](const JulianDate& jd) -> double {
        EphemerisData eph = pImpl->ephemeris.compute(jd);
        auto starPos = star.propagateToEpoch(jd);
        return Coordinates::angularSeparation(eph.geocentricPos, starPos);
    };
    
    // Golden section search
    const double phi = (1.0 + sqrt(5.0)) / 2.0;
    const double resphi = 2.0 - phi;
    
    double a = startJD.jd;
    double b = endJD.jd;
    double tol = 1.0 / 86400.0; // 1 secondo
    
    double x1 = a + resphi * (b - a);
    double x2 = b - resphi * (b - a);
    double f1 = calcSeparation(JulianDate(x1));
    double f2 = calcSeparation(JulianDate(x2));
    
    while (fabs(b - a) > tol) {
        if (f1 < f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + resphi * (b - a);
            f1 = calcSeparation(JulianDate(x1));
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = b - resphi * (b - a);
            f2 = calcSeparation(JulianDate(x2));
        }
    }
    
    return JulianDate((a + b) / 2.0);
}

void OccultationPredictor::calculateGeometry(const EphemerisData& asteroidPos,
                                             const GaiaStar& star,
                                             const JulianDate& time,
                                             double& separation,
                                             double& posAngle) {
    auto starPos = star.propagateToEpoch(time);
    
    // Applica correzione oblateness Terra (LinOccult fac = 0.996647)
    // Ratio polar/equatorial Earth radius
    const double EARTH_FLATTENING = 0.003353;  // WGS84
    const double POLAR_EQUAT_RATIO = 1.0 - EARTH_FLATTENING;  // 0.996647
    
    // Converti coordinate stella in cartesiane e applica correzione
    Vector3D starDir = Coordinates::equatorialToCartesian(starPos);
    
    // Schiacciamento polare: scala componente z (DEC) con ratio
    // Questo compensa geometria ellissoidale Terra nel calcolo ombra
    starDir = Vector3D(starDir.x, starDir.y, starDir.z / POLAR_EQUAT_RATIO);
    starDir = starDir.normalize();
    
    // Riconverti in coordinate equatoriali corrette
    EquatorialCoordinates correctedStarPos = Coordinates::cartesianToEquatorial(starDir);
    
    separation = Coordinates::angularSeparation(asteroidPos.geocentricPos, correctedStarPos);
    posAngle = Coordinates::positionAngle(asteroidPos.geocentricPos, correctedStarPos);
}

double OccultationPredictor::calculateProbability(double separationArcsec,
                                                  double asteroidAngularSize,
                                                  double uncertaintyArcsec) {
    if (uncertaintyArcsec <= 0) {
        // Senza incertezza, probabilità binaria
        return (separationArcsec <= asteroidAngularSize / 2.0) ? 1.0 : 0.0;
    }
    
    // Metodo LinOccult: integra CDF gaussiana tra distance ± radius
    // CalculateProbability() linee 1691-1708
    
    double radius = asteroidAngularSize / 2.0;
    double sigma = uncertaintyArcsec;
    
    // Normalizza distanze
    double x1 = (separationArcsec + radius) / sigma;  // Bordo esterno
    double x2 = (separationArcsec - radius) / sigma;  // Bordo interno
    
    // CDF gaussiana standard: Φ(x) = 0.5 * (1 + erf(x/√2))
    auto gaussCDF = [](double x) -> double {
        return 0.5 * (1.0 + erf(x / sqrt(2.0)));
    };
    
    // Probabilità = |Φ(x1) - Φ(x2)| = probabilità nell'intervallo [x2, x1]
    double g1 = gaussCDF(x1);
    double g2 = gaussCDF(x2);
    double prob = fabs(g1 - g2);
    
    // Limita a [0, 1]
    if (prob < 0) prob = 0;
    if (prob > 1) prob = 1;
    
    return prob;
}

bool OccultationEvent::isVisibleFrom(const GeographicCoordinates& observer,
                                     double minElevationDeg) const {
    // Verifica se l'evento è visibile da una posizione
    
    // Calcola la posizione dell'asteroide nel sistema locale dell'osservatore
    // Questa è una implementazione semplificata
    
    // LST dell'osservatore
    double lst = TimeUtils::lst(timeCA, observer.longitude);
    
    // Hour angle
    double ha = lst - star.pos.ra;
    
    // Altitude
    double sinAlt = sin(observer.latitude * DEG_TO_RAD) * sin(star.pos.dec) +
                   cos(observer.latitude * DEG_TO_RAD) * cos(star.pos.dec) * cos(ha);
    double alt = asin(sinAlt) * RAD_TO_DEG;
    
    return alt >= minElevationDeg;
}

} // namespace ioccultcalc
