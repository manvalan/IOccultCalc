#include "ioccultcalc/occultation_predictor.h"
#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/coordinates.h"
#include "ioccultcalc/time_utils.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace ioccultcalc {

class OccultationPredictor::Impl {
public:
    EquinoctialElements asteroid;
    Ephemeris ephemeris;
    GaiaClient gaiaClient;
    AstDysClient astdysClient;
    
    double asteroidDiameter;    // km
    double orbitalUncertainty;  // km (1-sigma)
    
    Impl() : asteroidDiameter(0), orbitalUncertainty(100.0) {}
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

void OccultationPredictor::setAsteroidDiameter(double diameter) {
    pImpl->asteroidDiameter = diameter;
}

void OccultationPredictor::setOrbitalUncertainty(double sigmaKm) {
    pImpl->orbitalUncertainty = sigmaKm;
}

std::vector<OccultationEvent> OccultationPredictor::findOccultations(
    const JulianDate& startJD,
    const JulianDate& endJD,
    double maxMagnitude,
    double searchRadius,
    double minProbability) {
    
    std::vector<OccultationEvent> events;
    
    // Step temporale per la ricerca (1 giorno)
    double stepDays = 1.0;
    
    // Calcola effemeridi per l'intervallo
    auto ephemerides = pImpl->ephemeris.computeRange(startJD, endJD, stepDays);
    
    // Per ogni epoca, cerca stelle vicine
    for (const auto& eph : ephemerides) {
        // Converti posizione in gradi
        double raDeg = eph.geocentricPos.ra * RAD_TO_DEG;
        double decDeg = eph.geocentricPos.dec * RAD_TO_DEG;
        
        // Query stelle nella regione
        auto stars = pImpl->gaiaClient.queryCone(raDeg, decDeg, searchRadius, maxMagnitude);
        
        // Per ogni stella, verifica se c'è un'occultazione
        for (const auto& star : stars) {
            // Propaga stella all'epoca
            auto starPos = star.propagateToEpoch(eph.jd);
            
            // Calcola separazione angolare
            double separation = Coordinates::angularSeparation(eph.geocentricPos, starPos);
            double separationArcsec = separation * RAD_TO_DEG * 3600.0;
            
            // Stima dimensione angolare dell'asteroide
            double asteroidAngularSize = 0;
            if (pImpl->asteroidDiameter > 0 && eph.distance > 0) {
                // Dimensione angolare in arcsec
                asteroidAngularSize = (pImpl->asteroidDiameter / (eph.distance * AU)) * RAD_TO_DEG * 3600.0;
            }
            
            // Controllo preliminare: stella deve essere vicina
            double threshold = asteroidAngularSize + 3.0 * pImpl->orbitalUncertainty / (eph.distance * AU) * RAD_TO_DEG * 3600.0;
            if (threshold == 0) threshold = 10.0; // Default 10 arcsec
            
            if (separationArcsec < threshold) {
                // Trova momento di closest approach
                JulianDate caTime = findClosestApproach(star, 
                    JulianDate(eph.jd.jd - 0.5), 
                    JulianDate(eph.jd.jd + 0.5));
                
                // Predici occultazione dettagliata
                OccultationEvent event = predictOccultation(star, caTime);
                
                if (event.probability >= minProbability) {
                    events.push_back(event);
                }
            }
        }
    }
    
    // Ordina per tempo
    std::sort(events.begin(), events.end(),
             [](const OccultationEvent& a, const OccultationEvent& b) {
                 return a.timeCA.jd < b.timeCA.jd;
             });
    
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
    
    // Propaga posizione stella
    auto starPosCA = star.propagateToEpoch(event.timeCA);
    
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
    
    // Calcola probabilità
    double uncertaintyArcsec = (pImpl->orbitalUncertainty / (ephCA.distance * AU)) * RAD_TO_DEG * 3600.0;
    event.probability = calculateProbability(event.closeApproachDistance, 
                                            asteroidAngularSize,
                                            uncertaintyArcsec);
    
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
    
    separation = Coordinates::angularSeparation(asteroidPos.geocentricPos, starPos);
    posAngle = Coordinates::positionAngle(asteroidPos.geocentricPos, starPos);
    
    // DEBUG: Prime 3 stelle
    static int debugCount = 0;
    if (debugCount < 3) {
        double astRA = asteroidPos.geocentricPos.ra * 180.0 / M_PI;
        double astDec = asteroidPos.geocentricPos.dec * 180.0 / M_PI;
        double starRA = starPos.ra * 180.0 / M_PI;
        double starDec = starPos.dec * 180.0 / M_PI;
        std::cout << "[GEOM " << debugCount << "] JD=" << time.jd << " Ast=(RA=" << astRA 
                  << "° Dec=" << astDec << "°) Star=(RA=" << starRA 
                  << "° Dec=" << starDec << "°) Sep=" << (separation * 180.0 / M_PI * 3600.0) << "\"\n";
        debugCount++;
    }
}

double OccultationPredictor::calculateProbability(double separationArcsec,
                                                  double asteroidAngularSize,
                                                  double uncertaintyArcsec) {
    if (uncertaintyArcsec <= 0) {
        // Senza incertezza, probabilità binaria
        return (separationArcsec <= asteroidAngularSize / 2.0) ? 1.0 : 0.0;
    }
    
    // Distribuzione gaussiana
    // Probabilità che la vera posizione cada entro il raggio dell'asteroide
    
    double r_asteroid = asteroidAngularSize / 2.0;
    double sigma = uncertaintyArcsec;
    
    // Approssimazione: probabilità che la separazione reale sia < r_asteroid
    // Usando CDF della normale
    
    double z = (r_asteroid - separationArcsec) / sigma;
    
    // erf approximation
    double prob = 0.5 * (1.0 + erf(z / sqrt(2.0)));
    
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
