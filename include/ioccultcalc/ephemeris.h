#ifndef IOCCULTCALC_EPHEMERIS_H
#define IOCCULTCALC_EPHEMERIS_H

#include "orbital_elements.h"
#include "types.h"
#include <vector>
#include <functional>
#include <memory>

namespace ioccultcalc {

// Forward declarations
class OrbitPropagator;
struct PropagatorOptions;

// Dati effemeridali per un istante
struct EphemerisData {
    JulianDate jd;
    EquatorialCoordinates geocentricPos; // Posizione geocentrica
    Vector3D heliocentricPos;            // Posizione eliocentrica (AU)
    Vector3D heliocentricVel;            // Velocità eliocentrica (AU/day)
    double distance;                      // Distanza dalla Terra (AU)
    double elongation;                    // Elongazione solare (gradi)
    double phase;                         // Angolo di fase (gradi)
    double magnitude;                     // Magnitudine visuale
    
    EphemerisData() : distance(0), elongation(0), phase(0), magnitude(0) {}
};

/**
 * @brief Opzioni per il calcolo delle effemeridi
 */
struct EphemerisOptions {
    bool usePerturbations;           // Usa perturbazioni planetarie (più preciso, più lento)
    bool useRelativisticCorrections; // Correzioni relativistiche
    double stepSize;                  // Passo integrazione (giorni) - solo con perturbazioni
    
    EphemerisOptions() 
        : usePerturbations(false),
          useRelativisticCorrections(false),
          stepSize(0.1) {}
          
    // Preset configurazioni
    static EphemerisOptions fast() {
        return EphemerisOptions();  // Default: solo Keplero
    }
    
    static EphemerisOptions standard() {
        EphemerisOptions opts;
        opts.usePerturbations = true;
        return opts;
    }
    
    static EphemerisOptions highPrecision() {
        EphemerisOptions opts;
        opts.usePerturbations = true;
        opts.useRelativisticCorrections = true;
        opts.stepSize = 0.05;
        return opts;
    }
};

class Ephemeris {
public:
    Ephemeris();
    explicit Ephemeris(const EquinoctialElements& elements);
    ~Ephemeris();
    
    // Imposta gli elementi orbitali
    void setElements(const EquinoctialElements& elements);
    
    // Imposta opzioni (perturbazioni, etc.)
    void setOptions(const EphemerisOptions& options);
    const EphemerisOptions& getOptions() const { return options_; }
    
    // Abilita/disabilita perturbazioni (shortcut)
    void enablePerturbations(bool enable = true);
    bool hasPerturbations() const { return options_.usePerturbations; }
    
    // Calcola la posizione dell'asteroide per una data epoca
    EphemerisData compute(const JulianDate& jd);
    
    // Calcola effemeridi per un intervallo temporale
    std::vector<EphemerisData> computeRange(const JulianDate& startJD, 
                                           const JulianDate& endJD,
                                           double stepDays = 1.0);
    
    // Ottiene la posizione del Sole geocentrica
    static Vector3D getSunPosition(const JulianDate& jd);
    
    // Ottiene la posizione della Terra eliocentrica
    static Vector3D getEarthPosition(const JulianDate& jd);
    static Vector3D getEarthVelocity(const JulianDate& jd);
    
    // Ottiene posizione Terra con correzioni complete (aberrazione, relatività)
    static Vector3D getEarthPositionWithCorrections(const JulianDate& jd, 
                                                     const Vector3D& observerPos);
    
    // Thread-safe Earth position cache (for parallel processing)
    using EarthPositionFunc = std::function<Vector3D(double jd)>;
    static void setEarthPositionCache(EarthPositionFunc cacheFunc);
    static void clearEarthPositionCache();
    static bool hasEarthPositionCache();
    
private:
    EquinoctialElements elements_;
    EphemerisOptions options_;
    
    // Propagatore con perturbazioni (lazy init)
    std::unique_ptr<OrbitPropagator> propagator_;
    bool propagatorInitialized_;
    JulianDate lastPropagatedEpoch_;
    Vector3D lastPropagatedPos_;
    Vector3D lastPropagatedVel_;
    
    // Inizializza propagatore se necessario
    void initializePropagator();
    
    // Propaga l'orbita da epoca elementi a epoca target (Keplero puro)
    void propagateOrbit(const JulianDate& targetJD, 
                       Vector3D& helioPos, Vector3D& helioVel);
    
    // Propaga con perturbazioni (usa OrbitPropagator)
    void propagateOrbitWithPerturbations(const JulianDate& targetJD,
                                         Vector3D& helioPos, Vector3D& helioVel);
    
    // Risolve l'equazione di Keplero per anomalia eccentrica
    double solveKeplerEquation(double meanAnomaly, double eccentricity, 
                              double tolerance = 1e-12);
    
    // Calcola la magnitudine apparente
    double calculateMagnitude(double r, double delta, double phaseAngle);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_EPHEMERIS_H
