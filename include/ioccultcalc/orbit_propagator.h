/**
 * @file orbit_propagator.h
 * @brief Propagatore numerico di orbite con perturbazioni
 * 
 * Integrazione numerica delle equazioni del moto con:
 * - Runge-Kutta 4° ordine (RK4)
 * - Runge-Kutta-Fehlberg 7/8 (RKF78) per alta precisione
 * - Perturbazioni planetarie (n-body)
 * - Correzioni relativistiche (opzionali)
 * 
 * Ispirato a MERCURY (Chambers 1999) e Find_Orb (Project Pluto)
 */

#ifndef IOCCULTCALC_ORBIT_PROPAGATOR_H
#define IOCCULTCALC_ORBIT_PROPAGATOR_H

#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/types.h"
#include "ioccultcalc/force_model.h"
#include <vector>
#include <functional>

namespace ioccultcalc {

/**
 * @enum IntegratorType
 * @brief Tipo di integratore numerico
 */
enum class IntegratorType {
    RK4,        // Runge-Kutta 4° ordine (step fisso) - CONSIGLIATO: ΔE/E~10⁻¹⁴
    RKF78,      // Runge-Kutta-Fehlberg 7/8 (adattivo, alta precisione)
    GAUSS_RADAU, // Gauss-Radau (implicito, per problemi stiff)
    RA15        // ⚠ BUGGY: non conserva energia (ΔE/E~10⁻⁶, 323km errore/anno)
};

/**
 * @struct PropagatorOptions
 * @brief Opzioni per propagazione orbitale
 */
struct PropagatorOptions {
    IntegratorType integrator;      // Tipo di integratore
    double stepSize;                 // Step iniziale (giorni)
    double tolerance;                // Tolleranza errore per integratori adattivi
    bool usePlanetaryPerturbations; // Perturbazioni planetarie
    bool useRelativisticCorrections; // Correzioni relativistiche
    bool useSolarRadiationPressure;  // Pressione radiazione solare
    double areaToMassRatio;          // A/m ratio (m²/kg) per SRP
    bool useNonGravitational;        // Usa parametri non gravitazionali A1,A2,A3
    int maxSteps;                    // Max passi integrazione (safety)
    
    PropagatorOptions() 
        : integrator(IntegratorType::RK4),
          stepSize(0.1),  // 0.1 giorni = miglior trade-off velocità/precisione
          tolerance(1e-12),
          usePlanetaryPerturbations(true),
          useRelativisticCorrections(false),
          useSolarRadiationPressure(false),
          areaToMassRatio(0.01),  // Default per asteroide tipico
          useNonGravitational(false),
          maxSteps(100000) {}
};

/**
 * @struct OrbitState
 * @brief Stato orbitale completo (posizione + velocità)
 */
struct OrbitState {
    JulianDate epoch;
    Vector3D position;  // AU (eliocentrico)
    Vector3D velocity;  // AU/day (eliocentrico)
    
    // Parametri non gravitazionali (VFCC17 model)
    double A1;  // Radial (AU/day²)
    double A2;  // Transverse (AU/day²)
    double A3;  // Normal (AU/day²)
    
    OrbitState() : A1(0), A2(0), A3(0) {}
    OrbitState(const JulianDate& jd, const Vector3D& r, const Vector3D& v)
        : epoch(jd), position(r), velocity(v), A1(0), A2(0), A3(0) {}
    OrbitState(const JulianDate& jd, const Vector3D& r, const Vector3D& v,
               double a1, double a2, double a3)
        : epoch(jd), position(r), velocity(v), A1(a1), A2(a2), A3(a3) {}
};

/**
 * @class OrbitPropagator
 * @brief Propagatore numerico di orbite
 * 
 * Integra le equazioni del moto:
 * d²r/dt² = -μ/r³ * r + Σ perturbazioni
 * 
 * Perturbazioni incluse:
 * - Attrazione gravitazionale di 8 pianeti + Luna
 * - Correzioni relativistiche (PN formulation)
 * - Pressione radiazione solare (opzionale)
 */
class OrbitPropagator {
public:
    OrbitPropagator();
    explicit OrbitPropagator(const PropagatorOptions& options);
    ~OrbitPropagator();
    
    /**
     * @brief Propaga elementi da epoca iniziale a epoca target
     * 
     * @param initialElements Elementi all'epoca iniziale
     * @param targetEpoch Epoca target
     * @return Stato orbitale (r, v) all'epoca target
     */
    OrbitState propagate(const EquinoctialElements& initialElements,
                        const JulianDate& targetEpoch);
    
    /**
     * @brief Propaga stato orbitale
     * 
     * @param initialState Stato iniziale (r, v)
     * @param targetEpoch Epoca target
     * @return Stato finale
     */
    OrbitState propagate(const OrbitState& initialState,
                        const JulianDate& targetEpoch);
    
    /**
     * @brief Propaga con output intermedi
     * 
     * @param initialState Stato iniziale
     * @param targetEpoch Epoca finale
     * @param outputStep Step per output intermedi (giorni)
     * @return Vettore di stati agli istanti intermedi
     */
    std::vector<OrbitState> propagateWithOutput(const OrbitState& initialState,
                                                const JulianDate& targetEpoch,
                                                double outputStep = 1.0);
    
    /**
     * @brief Imposta opzioni propagatore
     */
    void setOptions(const PropagatorOptions& options);
    
    /**
     * @brief Ottiene statistiche ultima propagazione
     */
    struct PropagationStats {
        int nSteps = 0;              // Passi integrazione
        int nEvaluations = 0;        // Valutazioni forza (per RA15)
        int nRejections = 0;         // Step rifiutati (per RA15)
        double finalStepSize = 0.0;  // Step finale
        double maxError = 0.0;       // Errore max stimato
        double computeTime = 0.0;    // Tempo CPU (secondi)
        double avgIterations = 0.0;  // Media iterazioni/step (per RA15)
    };
    
    PropagationStats getLastStats() const;
    
    /**
     * @brief Callback per monitoraggio progresso
    /**
     * Chiamato ogni N steps (parametro del callback)
     * Utile per barre di progresso
     */
    void setProgressCallback(std::function<void(double, const OrbitState&)> callback);
    
    /**
     * @brief Conversione elementi orbitali → stato cartesiano
     * @param elements Elementi equinoziali
     * @return Stato (posizione, velocità)
     */
    OrbitState elementsToState(const EquinoctialElements& elements);
    
    /**
     * @brief Conversione stato cartesiano → elementi orbitali
     * @param state Stato (posizione, velocità)
     * @return Elementi equinoziali
     */
    EquinoctialElements stateToElements(const OrbitState& state);
    
private:
    class Impl;
    Impl* pImpl;
    
    PropagatorOptions options_;
    PropagationStats lastStats_;
    std::function<void(double, const OrbitState&)> progressCallback_;
    
    // Parametri non gravitazionali correnti (VFCC17 model)
    double currentA1_;
    double currentA2_;
    double currentA3_;
    
    // Metodi interni
    OrbitState integrateRK4(const OrbitState& state0, double dt);
    OrbitState integrateRKF78(const OrbitState& state0, double dt, double& errorEst);
    OrbitState integrateRA15(const OrbitState& state0, const JulianDate& targetEpoch);
    
public:  // ← TEMPORANEO per test validazione
    // Calcola accelerazione totale (gravitazione + perturbazioni)
    Vector3D computeAcceleration(const JulianDate& jd,
                                 const Vector3D& pos,
                                 const Vector3D& vel);
};

/**
 * @brief Test di accuratezza propagatore
 * 
 * Propaga orbita e confronta con effemeridi JPL DE
 * 
 * @param elements Elementi da testare
 * @param duration Durata propagazione (giorni)
 * @return Errore RMS in posizione (km)
 */
double testPropagatorAccuracy(const EquinoctialElements& elements,
                             double duration = 365.25);

} // namespace ioccultcalc

#endif // IOCCULTCALC_ORBIT_PROPAGATOR_H
