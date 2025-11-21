/**
 * @file numerical_integrator.h
 * @brief Integratori numerici per propagazione orbitale precisa
 * 
 * Implementa:
 * - RKF78: Runge-Kutta-Fehlberg 7(8) con step adattivo
 * - RK89: Runge-Kutta 8(9) di Verner per massima precisione
 * - DOPRI853: Dormand-Prince 8(5,3) con interpolazione densa
 * - Symplectic: Integratore simplettico per conservazione energia
 * 
 * Include equazioni variazionali per state transition matrix (STM)
 * per propagazione incertezze orbitali
 */

#ifndef IOCCULTCALC_NUMERICAL_INTEGRATOR_H
#define IOCCULTCALC_NUMERICAL_INTEGRATOR_H

#include "ioccultcalc/types.h"
#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/force_model.h"
#include <vector>
#include <functional>

namespace ioccultcalc {

/**
 * @struct OrbitalState
 * @brief Stato completo orbita (posizione + velocità)
 */
struct OrbitalState {
    JulianDate epoch;
    Vector3D position;  // AU
    Vector3D velocity;  // AU/giorno
    
    // State Transition Matrix (6x6) per propagazione incertezze
    std::vector<std::vector<double>> stm;
    
    OrbitalState();
    
    // Converte da elementi orbitali
    static OrbitalState fromElements(const EquinoctialElements& elements);
    
    // Converte a elementi orbitali
    EquinoctialElements toElements() const;
    
    // Energia orbitale
    double energy() const;
    
    // Momento angolare
    Vector3D angularMomentum() const;
};

/**
 * @typedef ForceFunction
 * @brief Funzione che calcola accelerazione
 * 
 * @param t Tempo (JD)
 * @param pos Posizione (AU)
 * @param vel Velocità (AU/giorno)
 * @return Accelerazione (AU/giorno²)
 */
using ForceFunction = std::function<Vector3D(double t, 
                                              const Vector3D& pos,
                                              const Vector3D& vel)>;

/**
 * @struct IntegratorOptions
 * @brief Opzioni per integratore numerico
 */
struct IntegratorOptions {
    double relTolerance;      // Tolleranza relativa (default: 1e-12)
    double absTolerance;      // Tolleranza assoluta (default: 1e-15)
    double minStep;           // Step minimo (giorni, default: 1e-6)
    double maxStep;           // Step massimo (giorni, default: 10.0)
    double initialStep;       // Step iniziale (giorni, default: 0.1)
    bool computeSTM;          // Calcola State Transition Matrix
    bool enforceEnergyConservation; // Correzione drift energia
    int maxSubsteps;          // Numero max valutazioni forza (default: 100000)
    
    IntegratorOptions();
};

/**
 * @struct IntegrationResult
 * @brief Risultato integrazione
 */
struct IntegrationResult {
    OrbitalState finalState;
    int numberOfSteps;
    int numberOfRejectedSteps;
    double energyError;       // Errore relativo energia
    bool success;
    std::string errorMessage;
    
    IntegrationResult();
};

/**
 * @class NumericalIntegrator
 * @brief Classe base per integratori numerici
 */
class NumericalIntegrator {
public:
    virtual ~NumericalIntegrator() = default;
    
    /**
     * @brief Propaga orbita da t0 a t1
     * @param initialState Stato iniziale
     * @param finalTime Tempo finale (JD)
     * @param forceFunc Funzione accelerazione
     * @param options Opzioni integrazione
     * @return Risultato integrazione
     */
    virtual IntegrationResult propagate(
        const OrbitalState& initialState,
        const JulianDate& finalTime,
        const ForceFunction& forceFunc,
        const IntegratorOptions& options = IntegratorOptions()) = 0;
    
    /**
     * @brief Propaga con output intermedi
     * @param initialState Stato iniziale
     * @param outputTimes Tempi di output
     * @param forceFunc Funzione accelerazione
     * @param options Opzioni
     * @return Vettore stati ai tempi richiesti
     */
    virtual std::vector<OrbitalState> propagateWithOutput(
        const OrbitalState& initialState,
        const std::vector<JulianDate>& outputTimes,
        const ForceFunction& forceFunc,
        const IntegratorOptions& options = IntegratorOptions()) = 0;
};

/**
 * @class RKF78Integrator
 * @brief Runge-Kutta-Fehlberg 7(8) con step adattivo
 * 
 * Metodo embedded di ordine 7 con stimatore errore di ordine 8
 * Molto preciso per orbite perturbate
 * 
 * Riferimento: Fehlberg (1968), NASA TR R-287
 */
class RKF78Integrator : public NumericalIntegrator {
public:
    IntegrationResult propagate(
        const OrbitalState& initialState,
        const JulianDate& finalTime,
        const ForceFunction& forceFunc,
        const IntegratorOptions& options = IntegratorOptions()) override;
    
    std::vector<OrbitalState> propagateWithOutput(
        const OrbitalState& initialState,
        const std::vector<JulianDate>& outputTimes,
        const ForceFunction& forceFunc,
        const IntegratorOptions& options = IntegratorOptions()) override;
    
private:
    // Coefficienti Butcher tableau RKF78
    static const double c[13];
    static const double a[13][12];
    static const double b7[13];  // Ordine 7
    static const double b8[13];  // Ordine 8
    
    // Singolo step
    bool step(double& t, Vector3D& pos, Vector3D& vel,
              double& h, const ForceFunction& forceFunc,
              const IntegratorOptions& options,
              double& error);
    
    // Interpolazione densa (per output intermedi)
    OrbitalState interpolate(double t_interp,
                            double t0, double t1,
                            const OrbitalState& state0,
                            const std::vector<Vector3D>& k_stages);
};

/**
 * @class DOPRI853Integrator
 * @brief Dormand-Prince 8(5,3) con interpolazione densa
 * 
 * Metodo molto efficiente con interpolazione cubica per output
 * Ordine 8 con stimatore locale ordine 5 e ordine 3
 * 
 * Riferimento: Hairer, Nørsett, Wanner (1993)
 */
class DOPRI853Integrator : public NumericalIntegrator {
public:
    IntegrationResult propagate(
        const OrbitalState& initialState,
        const JulianDate& finalTime,
        const ForceFunction& forceFunc,
        const IntegratorOptions& options = IntegratorOptions()) override;
    
    std::vector<OrbitalState> propagateWithOutput(
        const OrbitalState& initialState,
        const std::vector<JulianDate>& outputTimes,
        const ForceFunction& forceFunc,
        const IntegratorOptions& options = IntegratorOptions()) override;
};

/**
 * @class SymplecticIntegrator
 * @brief Integratore simplettico per conservazione hamiltoniana
 * 
 * Conserva esattamente energia + momento angolare a lungo termine
 * Ideale per propagazioni multi-orbitali (anni/secoli)
 * 
 * Implementa metodo Yoshida di ordine 6
 */
class SymplecticIntegrator : public NumericalIntegrator {
public:
    IntegrationResult propagate(
        const OrbitalState& initialState,
        const JulianDate& finalTime,
        const ForceFunction& forceFunc,
        const IntegratorOptions& options = IntegratorOptions()) override;
    
    std::vector<OrbitalState> propagateWithOutput(
        const OrbitalState& initialState,
        const std::vector<JulianDate>& outputTimes,
        const ForceFunction& forceFunc,
        const IntegratorOptions& options = IntegratorOptions()) override;
    
private:
    // Coefficienti Yoshida order 6
    static const double w[8];
};

/**
 * @class VariationalEquations
 * @brief Equazioni variazionali per State Transition Matrix
 * 
 * Propaga matrice 6x6 che mappa incertezze da t0 a t1:
 * δx(t) = Φ(t,t0) * δx(t0)
 * 
 * Usata per propagare ellissoidi di incertezza orbitale
 */
class VariationalEquations {
public:
    /**
     * @brief Calcola derivate STM
     * dΦ/dt = A(t) * Φ
     * dove A è la matrice Jacobiana del campo di forze
     * 
     * @param t Tempo
     * @param state Stato corrente
     * @param forceFunc Funzione forza
     * @return Derivata Φ̇ (linearizzata come vettore 36-elementi)
     */
    static std::vector<double> computeDerivatives(
        double t,
        const OrbitalState& state,
        const ForceFunction& forceFunc);
    
    /**
     * @brief Calcola matrice Jacobiana del campo di forze
     * A = [  0      I  ]
     *     [ ∂a/∂r  ∂a/∂v ]
     * 
     * @param t Tempo
     * @param pos Posizione
     * @param vel Velocità
     * @param forceFunc Funzione forza
     * @return Matrice 6x6
     */
    static std::vector<std::vector<double>> computeJacobian(
        double t,
        const Vector3D& pos,
        const Vector3D& vel,
        const ForceFunction& forceFunc);
    
    /**
     * @brief Propaga ellissoide di incertezza
     * 
     * Usa STM per mappare matrice covarianza da t0 a t1:
     * P(t1) = Φ(t1,t0) * P(t0) * Φ(t1,t0)^T
     * 
     * @param initialCovariance Matrice 6x6 covarianza iniziale
     * @param stm State Transition Matrix
     * @return Matrice covarianza finale
     */
    static std::vector<std::vector<double>> propagateCovariance(
        const std::vector<std::vector<double>>& initialCovariance,
        const std::vector<std::vector<double>>& stm);
    
    /**
     * @brief Estrae incertezza posizionale da covarianza
     * @param covariance Matrice 6x6
     * @param sigma Livello sigma (default: 1.0)
     * @return Semiassi ellissoide incertezza (AU)
     */
    static Vector3D extractPositionUncertainty(
        const std::vector<std::vector<double>>& covariance,
        double sigma = 1.0);
};

/**
 * @class AdaptiveStepControl
 * @brief Controllo step adattivo per integratori
 */
class AdaptiveStepControl {
public:
    /**
     * @brief Calcola nuovo step ottimale
     * 
     * Usa formula PI controller:
     * h_new = h * (tol / err)^(1/order) * safety_factor
     * 
     * @param currentStep Step corrente
     * @param error Errore stimato
     * @param tolerance Tolleranza richiesta
     * @param order Ordine metodo
     * @param options Opzioni integrazione
     * @return Nuovo step suggerito
     */
    static double computeNewStep(
        double currentStep,
        double error,
        double tolerance,
        int order,
        const IntegratorOptions& options);
    
    /**
     * @brief Verifica se step è accettabile
     * @param error Errore stimato
     * @param tolerance Tolleranza
     * @return true se step accettato
     */
    static bool acceptStep(double error, double tolerance);
    
private:
    static constexpr double SAFETY_FACTOR = 0.9;
    static constexpr double MIN_SCALE = 0.2;
    static constexpr double MAX_SCALE = 5.0;
};

/**
 * @namespace IntegratorUtils
 * @brief Utility per integrazione orbitale con modelli di forze
 */
namespace IntegratorUtils {
    /**
     * @brief Crea ForceFunction da ForceModel
     * 
     * Wrapper che converte ForceModel::computeAcceleration in ForceFunction
     * utilizzabile dagli integratori numerici
     * 
     * @param model Modello di forze (deve rimanere valido durante integrazione)
     * @return ForceFunction compatibile con integratori
     */
    ForceFunction createForceFunctionFromModel(const class ForceModel& model);
    
    /**
     * @brief Propaga orbita usando ForceModel completo
     * 
     * High-level function che combina integratore + force model
     * 
     * @param initialState Stato orbitale iniziale
     * @param finalTime Tempo finale
     * @param forceConfig Configurazione force model
     * @param integOptions Opzioni integratore
     * @return Risultato integrazione
     */
    IntegrationResult propagateWithForceModel(
        const OrbitalState& initialState,
        const JulianDate& finalTime,
        const class ForceModelConfig& forceConfig,
        const IntegratorOptions& integOptions = IntegratorOptions());
    
    /**
     * @brief Propaga con output intermedi usando ForceModel
     * 
     * @param initialState Stato iniziale
     * @param outputTimes Tempi di output
     * @param forceConfig Configurazione force model
     * @param integOptions Opzioni integratore
     * @return Stati ai tempi richiesti
     */
    std::vector<OrbitalState> propagateWithOutputAndForceModel(
        const OrbitalState& initialState,
        const std::vector<JulianDate>& outputTimes,
        const class ForceModelConfig& forceConfig,
        const IntegratorOptions& integOptions = IntegratorOptions());
}

} // namespace ioccultcalc

#endif // IOCCULTCALC_NUMERICAL_INTEGRATOR_H
