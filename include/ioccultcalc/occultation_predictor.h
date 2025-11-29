#ifndef IOCCULTCALC_OCCULTATION_PREDICTOR_H
#define IOCCULTCALC_OCCULTATION_PREDICTOR_H

#include "orbital_elements.h"
#include "ephemeris.h"
#include "gaia_client.h"
#include "chebyshev_approximation.h"
#include "types.h"
#include <string>
#include <vector>
#include <memory>

namespace ioccultcalc {

// Punto sulla traccia dell'occultazione
struct ShadowPathPoint {
    JulianDate time;
    GeographicCoordinates location;
    double duration;        // durata in secondi
    double centerlineDistance; // distanza dalla linea centrale (km)
    
    ShadowPathPoint() : duration(0), centerlineDistance(0) {}
};

// Informazioni complete su un'occultazione
struct OccultationEvent {
    std::string eventId;
    EquinoctialElements asteroid;
    GaiaStar star;
    
    JulianDate timeCA;              // Tempo di closest approach
    double closeApproachDistance;   // Distanza minima (arcsec)
    double positionAngle;            // Angolo di posizione (gradi)
    double probability;              // Probabilità di occultazione
    double maxDuration;              // Durata massima (secondi)
    
    std::vector<ShadowPathPoint> shadowPath;  // Traccia sulla Terra
    
    // Limiti di incertezza (1-sigma)
    double uncertaintyNorth;  // km
    double uncertaintySouth;  // km
    
    OccultationEvent() 
        : closeApproachDistance(0), positionAngle(0), 
          probability(0), maxDuration(0),
          uncertaintyNorth(0), uncertaintySouth(0) {}
    
    // Verifica se l'evento è visibile da una certa posizione
    bool isVisibleFrom(const GeographicCoordinates& observer, 
                      double minElevationDeg = 15.0) const;
};

class OccultationPredictor {
public:
    OccultationPredictor();
    ~OccultationPredictor();
    
    // Imposta l'asteroide da analizzare
    void setAsteroid(const EquinoctialElements& elements);
    
    // Carica l'asteroide da AstDyS
    void loadAsteroidFromAstDyS(const std::string& designation);
    
    // Configura directory locali per file AstDyS
    void setLocalAstDySDirectories(const std::string& eq1Dir, const std::string& rwoDir);
    
    // Cerca occultazioni in un intervallo temporale
    std::vector<OccultationEvent> findOccultations(
        const JulianDate& startJD,
        const JulianDate& endJD,
        double maxMagnitude = 12.0,     // Magnitudine limite stelle
        double searchRadius = 0.1,       // Raggio di ricerca in gradi
        double minProbability = 0.01     // Probabilità minima
    );
    
    // Calcola i dettagli di un'occultazione per una stella specifica
    OccultationEvent predictOccultation(const GaiaStar& star, 
                                       const JulianDate& approximateTime);
    
    // Calcola la traccia dell'ombra sulla Terra
    std::vector<ShadowPathPoint> calculateShadowPath(
        const EphemerisData& asteroidPos,
        const GaiaStar& star,
        const JulianDate& centralTime,
        double timeSpanMinutes = 120.0
    );
    
    // Imposta il diametro dell'asteroide (km) per calcoli più precisi
    void setAsteroidDiameter(double diameter);
    
    // Imposta l'errore orbitale (sigma in km) per calcolo incertezze
    void setOrbitalUncertainty(double sigmaKm);
    
    // Configura approssimazione Chebyshev per performance
    void setChebyshevConfig(const ChebyshevConfig& config);
    
    // Ottiene configurazione Chebyshev corrente
    const ChebyshevConfig& getChebyshevConfig() const;
    
    // ===== OPZIONI AVANZATE (OPZIONALI) =====
    
    // Abilita correzione parallasse per stelle vicine
    // Default: false (per retrocompatibilità)
    // Impatto: +0.5" accuratezza per stelle con parallax > 10 mas
    void setParallaxCorrection(bool enabled);
    
    // Imposta buffer angolare per ricerca stelle (arcmin)
    // Default: 0.0 (nessun buffer)
    // Raccomandato: 1.0 arcmin per evitare edge cases
    void setAngularBuffer(double arcmin);
    
    // Abilita timestep adattivo (scan grossolano + refine automatico)
    // Default: false (timestep fisso)
    // Impatto: 10× speedup per scan lunghi, cattura tutti eventi
    void setAdaptiveTimestep(bool enabled, double coarseStepDays = 1.0, 
                            double fineStepMinutes = 15.0);
    
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
    
    // Trova il momento di closest approach tra asteroide e stella
    JulianDate findClosestApproach(const GaiaStar& star,
                                   const JulianDate& startJD,
                                   const JulianDate& endJD);
    
    // Calcola la geometria dell'occultazione
    void calculateGeometry(const EphemerisData& asteroidPos,
                          const GaiaStar& star,
                          const JulianDate& time,
                          double& separation,
                          double& posAngle);
    
    // Calcola la probabilità di occultazione
    double calculateProbability(double separation, 
                               double asteroidAngularSize,
                               double uncertainty);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_OCCULTATION_PREDICTOR_H
