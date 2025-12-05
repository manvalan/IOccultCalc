#ifndef IOCCULTCALC_PROPAGATION_STRATEGY_H
#define IOCCULTCALC_PROPAGATION_STRATEGY_H

/**
 * @file propagation_strategy.h
 * @brief Strategia di propagazione a due fasi semplificata per IOccultCalc
 * 
 * Implementa una strategia ottimizzata usando solo AstDyn RKF78:
 * - FASE 1 (Screening): AstDyn RKF78 con tolleranza ridotta (1e-10) - veloce
 * - FASE 2 (Closest Approach): AstDyn RKF78 con tolleranza massima (1e-12) - preciso
 * 
 * @author Michele Bigi
 * @date 30 Novembre 2025
 */

#include "orbital_elements.h"
#include "astdyn_interface.h"
#include <string>
#include <memory>

namespace ioccultcalc {

/**
 * @brief Configurazione strategia propagazione semplificata
 */
struct PropagationConfig {
    // Fase 1: Screening candidati (AstDyn RKF78 con tolleranza ridotta)
    double screening_threshold_arcsec = 60.0;  // Soglia per promuovere a fase 2
    
    // Fase 2: Closest approach (AstDyn RKF78 con tolleranza massima)
    double rkf78_tolerance = 1e-12;            // Tolleranza integratore RKF78 FASE 2
    double convergence_seconds = 1.0;          // Convergenza golden section
    
    // Fitting orbitale - Controllo granulare
    enum class FittingMode {
        NEVER,           // Non usa mai orbit fitting (elementi originali)
        AUTO,           // Usa fitting se osservazioni disponibili
        ALWAYS_ATTEMPT  // Tenta sempre fitting (errore se osservazioni mancanti)
    };
    
    FittingMode orbit_fitting_mode = FittingMode::AUTO;  // Modalità orbit fitting
    double fitting_sigma_cutoff = 3.0;         // Soglia outlier σ per fitting  
    int max_fitting_iterations = 10;           // Max iterazioni differential correction
    int min_observations_for_fitting = 3;      // Minimo osservazioni per tentare fitting
    bool force_refit_each_prediction = false;  // Rifit elementi ad ogni closest approach
    
    // Retrocompatibilità (deprecated)
    bool enable_orbit_fitting = false;         // DEPRECATO: usa orbit_fitting_mode
    
    // Debug e performance
    bool verbose_timing = false;
    bool log_phase_transitions = false;
    bool log_fitting_details = false;          // Log dettagli orbit fitting
};

/**
 * @brief Strategia di propagazione a due fasi semplificata
 */
class TwoPhaseStrategy {
public:
    /**
     * @brief Coordinate equatoriali per propagazione
     */
    struct EquatorialCoords {
        double ra_deg;
        double dec_deg;
        double distance_au;
    };
    
    explicit TwoPhaseStrategy(const PropagationConfig& config = PropagationConfig());
    
    /**
     * @brief Imposta configurazione
     */
    void setConfig(const PropagationConfig& config);
    PropagationConfig getConfig() const { return config_; }
    
    /**
     * @brief Imposta elementi orbitali
     */
    void setElements(const AstDySElements& elements);
    void setElements(const OrbitalElements& elements);
    
    /**
     * @brief Carica osservazioni per orbit fitting (opzionale)
     */
    void setObservations(const std::vector<RWOObservation>& observations);
    
    /**
     * @brief Carica osservazioni da file RWO
     * @param rwo_file Path al file .rwo
     * @return true se caricamento riuscito
     */
    bool loadObservationsFromFile(const std::string& rwo_file);
    
    /**
     * @brief Carica osservazioni da AstDyS per un asteroide specifico
     * @param asteroid_number Numero asteroide
     * @return true se download e caricamento riusciti
     */
    bool loadObservationsFromAstDyS(int asteroid_number);
    
    bool hasObservations() const { return !observations_.empty(); }
    size_t getObservationCount() const { return observations_.size(); }
    
    /**
     * @brief Informazioni sulle osservazioni caricate
     */
    struct ObservationInfo {
        size_t total_observations = 0;
        double first_observation_mjd = 0.0;
        double last_observation_mjd = 0.0;
        double observation_arc_days = 0.0;
        std::string designation = "";
        int rejected_outliers = 0;
        double mean_rms_arcsec = 0.0;
        bool quality_acceptable = false;
    };
    
    /**
     * @brief Analizza qualità delle osservazioni caricate
     */
    ObservationInfo analyzeObservations() const;
    
    /**
     * @brief Valida e pulisce le osservazioni
     * @param max_rms_arcsec Soglia massima RMS per osservazione [arcsec]
     * @return Numero osservazioni rimosse come outliers
     */
    int cleanObservations(double max_rms_arcsec = 5.0);
    
    /**
     * @brief Controllo orbit fitting
     */
    void setFittingMode(PropagationConfig::FittingMode mode) { 
        config_.orbit_fitting_mode = mode; 
    }
    PropagationConfig::FittingMode getFittingMode() const { 
        return config_.orbit_fitting_mode; 
    }
    bool isUsingFittedElements() const { return use_fitted_elements_; }
    void resetFittedElements() { use_fitted_elements_ = false; }  // Forza re-fit prossima volta
    
    /**
     * @brief Trigger manuale orbit fitting
     * @return true se fitting successful
     */
    bool triggerOrbitFitting();
    
    /**
     * @brief FASE 1: Screening veloce con polinomi Chebyshev
     * Genera posizioni RKF78 precise, calcola polinomi Chebyshev, testa velocemente stelle
     * 
     * @param target_mjd Epoca target
     * @param ra_star RA stella [deg]  
     * @param dec_star Dec stella [deg]
     * @return Separazione angolare [arcsec]
     */
    double screenCandidate(double target_mjd, double ra_star, double dec_star);
    
    /**
     * @brief Ottieni posizione asteroide usando FASE 1 (Chebyshev)
     * @param target_mjd Epoca target
     * @return Coordinate RA/Dec [deg]
     */
    EquatorialCoords getChebyshevPosition(double target_mjd);
    
    /**
     * @brief Ottieni posizione asteroide usando RKF78 diretto
     * @param target_mjd Epoca target
     * @return Coordinate RA/Dec [deg]
     */
    EquatorialCoords getRKF78Position(double target_mjd);
    
    /**
     * @brief FASE 2: Closest approach preciso con AstDyn RKF78
     * Usa AstDyn RKF78 direttamente con tolleranza massima (1e-12) per massima precisione
     * 
     * @param target_mjd Epoca target approssimativa
     * @param ra_star RA stella [deg]
     * @param dec_star Dec stella [deg]  
     * @param search_window_hours Finestra ricerca [ore]
     * @return Tempo e separazione minima
     */
    struct CloseApproachResult {
        double closest_time_mjd;
        double minimum_separation_arcsec;
        double position_angle_deg;
        
        // Orbit fitting result (se abilitato)
        bool orbit_was_fitted = false;
        double orbit_improvement_arcsec = 0.0;
        int fitting_iterations = 0;
        double final_rms_arcsec = 0.0;
        
        // Performance
        int computation_steps;
        double computation_time_ms;
    };
    
    CloseApproachResult findClosestApproach(double target_mjd, 
                                          double ra_star, double dec_star,
                                          double search_window_hours = 2.0);
    
    /**
     * @brief Statistiche performance
     */
    struct PerformanceStats {
        int phase1_calls = 0;
        int phase2_calls = 0;
        double phase1_total_time_ms = 0.0;
        double phase2_total_time_ms = 0.0;
        int candidates_screened = 0;
        int candidates_promoted = 0;
        int orbits_fitted = 0;
        double total_orbit_improvement_arcsec = 0.0;
    };
    
    PerformanceStats getStats() const { return stats_; }
    void resetStats();
    
private:
    PropagationConfig config_;
    AstDySElements elements_;
    AstDySElements fitted_elements_;  // Elementi dopo fitting
    std::vector<RWOObservation> observations_;
    bool elements_set_ = false;
    bool use_fitted_elements_ = false;
    mutable PerformanceStats stats_;
    
    // Propagatori interni
    std::unique_ptr<class AstDynPropagator> astdyn_;
    
    // Polinomi di Chebyshev per screening veloce FASE 1
    struct ChebyshevPolynomials {
        std::vector<double> ra_coeffs;    // Coefficienti RA
        std::vector<double> dec_coeffs;   // Coefficienti Dec  
        double t_start_mjd;               // Inizio intervallo
        double t_end_mjd;                 // Fine intervallo
        int degree;                       // Grado polinomio
        bool valid = false;               // Polinomi calcolati
    };
    mutable ChebyshevPolynomials chebyshev_;
    
    /**
     * @brief FASE 1: Propagazione Kepleriana geocentrica veloce
     */
    EquatorialCoords propagateKeplerianGeocentric(double target_mjd);
    
    /**
     * @brief FASE 2: Propagazione AstDyn RKF78
     */
    EquatorialCoords propagateRKF78(double target_mjd);
    
    /**
     * @brief Orbit fitting con osservazioni
     */
    bool performOrbitFitting();
    
    /**
     * @brief Separazione angolare
     */
    double angularSeparation(double ra1_deg, double dec1_deg, 
                           double ra2_deg, double dec2_deg);
    
    /**
     * @brief FASE 1: Genera polinomi di Chebyshev da posizioni RKF78
     */
    void generateChebyshevPolynomials(double start_mjd, double end_mjd, int degree = 8);
    
    /**
     * @brief FASE 1: Valuta posizione usando polinomi di Chebyshev
     */
    EquatorialCoords evaluateChebyshev(double target_mjd) const;
    
    /**
     * @brief FASE 1: Verifica validità polinomi per epoca target
     */
    bool areChebyshevPolynomialsValid(double target_mjd) const;
    
    /**
     * @brief Conversione elementi per AstDyn
     */
    void updateAstDynElements();
    
    /**
     * @brief Metodi ausiliari per polinomi di Chebyshev
     */
    std::vector<double> calculateChebyshevCoeffs(const std::vector<double>& times, 
                                               const std::vector<double>& values, int degree);
    double evaluateChebyshevPolynomial(const std::vector<double>& coeffs, double x) const;
    double chebyshevBasisFunction(int n, double x) const;
    double interpolateLinear(const std::vector<double>& x, const std::vector<double>& y, 
                           double target_x) const;
};

/**
 * @brief Factory per configurazioni predefinite
 */
namespace propagation_presets {
    
    /**
     * @brief Configurazione veloce per survey
     * - Chebyshev fase 1, RKF78 fase 2
     * - NO orbit fitting (FittingMode::NEVER)
     */
    PropagationConfig createFastSurvey();
    
    /**
     * @brief Configurazione bilanciata
     * - Chebyshev fase 1, RKF78 fase 2
     * - Orbit fitting AUTO se disponibili osservazioni
     */
    PropagationConfig createBalanced();
    
    /**
     * @brief Configurazione per occultazioni finali
     * - Chebyshev fase 1, RKF78 fase 2 alta precisione
     * - Orbit fitting ALWAYS_ATTEMPT (richiede osservazioni)
     */
    PropagationConfig createPrecision();
    
    /**
     * @brief Configurazione per testing con fitting controllato
     * - Debug dettagliato orbit fitting
     * - Modalità AUTO con log completi
     */
    PropagationConfig createTestingWithFitting();
    
    /**
     * @brief Configurazione senza orbit fitting
     * - Usa solo elementi AstDyS originali
     * - FittingMode::NEVER con debug
     */
    PropagationConfig createNoFitting();
    
} // namespace propagation_presets

} // namespace ioccultcalc

#endif // IOCCULTCALC_PROPAGATION_STRATEGY_H
