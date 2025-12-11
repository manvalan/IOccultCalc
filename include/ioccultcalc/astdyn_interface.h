#ifndef IOCCULTCALC_ASTDYN_INTERFACE_H
#define IOCCULTCALC_ASTDYN_INTERFACE_H

/**
 * @file astdyn_interface.h
 * @brief Interfaccia opzionale a ITALOccultLibrary/AstDyn
 * 
 * Fornisce:
 * - Caricamento elementi orbitali da AstDyS (formato .eq1)
 * - Caricamento osservazioni RWO da AstDyS (formato .rwo)
 * - Propagazione con RKF78 ad alta precisione
 * - Fitting orbitale con osservazioni
 * - Calcolo residui O-C (Observed - Computed)
 * 
 * Integrazione opzionale: richiede ITALOccultLibrary/astdyn
 * 
 * @author Michele Bigi
 * @date 29 Novembre 2025
 * @version 1.0
 */

#include <string>
#include <vector>
#include <memory>
#include <optional>
#include "orbital_elements.h"
#include "observation.h"
#include "ioccultcalc/astdys_client.h"

namespace ioccultcalc {

/**
 * @struct AstDySElements
 * @brief Elementi orbitali nel formato AstDyS (equinoziali o kepleriani)
 */
struct AstDySElements {
    std::string name;              // Nome/designazione asteroide
    int number;                    // Numero asteroide (0 se non numerato)
    
    // Elementi kepleriani
    double a;                      // Semiasse maggiore [AU]
    double e;                      // Eccentricità
    double i;                      // Inclinazione [deg]
    double Omega;                  // Longitudine nodo ascendente [deg]
    double omega;                  // Argomento del perielio [deg]
    double M;                      // Anomalia media [deg]
    double epoch_mjd;              // Epoca [MJD TDT]
    
    // Parametri fisici
    double H;                      // Magnitudine assoluta
    double G;                      // Slope parameter
    
    // Matrice covarianza (opzionale)
    bool has_covariance;
    std::vector<double> covariance; // Triangolo superiore 6x6 (21 elementi)
    
    // Conversione a OrbitalElements di IOccultCalc
    OrbitalElements toOrbitalElements() const;
    
    // Crea da file .eq1
    static AstDySElements fromFile(const std::string& filename);
    
    // Download da AstDyS
    static AstDySElements download(int asteroid_number);
    static AstDySElements download(const std::string& designation);
};

/**
 * @struct RWOObservation
 * @brief Osservazione astrometrica dal formato RWO di AstDyS
 */
struct RWOObservation {
    std::string designation;       // Designazione oggetto
    double mjd_utc;                // Epoca [MJD UTC]
    double ra_deg;                 // RA osservata [deg]
    double dec_deg;                // Dec osservata [deg]
    double ra_sigma_arcsec;        // Errore RA [arcsec]
    double dec_sigma_arcsec;       // Errore Dec [arcsec]
    std::string obs_code;          // Codice osservatorio MPC
    double magnitude;              // Magnitudine osservata
    std::string band;              // Banda fotometrica
    
    // Residui O-C (calcolati dal fitter)
    double ra_residual_arcsec;     // O-C in RA*cos(dec) [arcsec]
    double dec_residual_arcsec;    // O-C in Dec [arcsec]
    double chi_squared;            // Chi² normalizzato
    bool is_outlier;               // Flag outlier (> 3σ)
    
    // Conversione a Observation di IOccultCalc
    AstrometricObservation toObservation() const;
    
    // Carica da file .rwo
    static std::vector<RWOObservation> fromFile(const std::string& filename);
    
    // Download da AstDyS
    static std::vector<RWOObservation> download(int asteroid_number);
    static std::vector<RWOObservation> download(const std::string& designation);
};

/**
 * @struct AstDynOrbitFitResult
 * @brief Risultato del fitting orbitale con statistiche (AstDyn)
 */
struct AstDynOrbitFitResult {
    AstDySElements fitted_elements;       // Elementi migliorati
    
    // Statistiche globali
    int n_observations;                    // Numero osservazioni totali
    int n_used;                            // Osservazioni usate (non outlier)
    int n_outliers;                        // Outlier rilevati (> 3σ)
    
    double rms_ra_arcsec;                  // RMS residui RA [arcsec]
    double rms_dec_arcsec;                 // RMS residui Dec [arcsec]
    double rms_total_arcsec;               // RMS totale [arcsec]
    
    double chi2;                           // Chi²
    double chi2_reduced;                   // Chi² ridotto
    
    double mean_ra_residual;               // Media residui RA [arcsec]
    double mean_dec_residual;              // Media residui Dec [arcsec]
    double max_residual_arcsec;            // Residuo massimo [arcsec]
    
    double time_span_days;                 // Arco temporale [giorni]
    double first_obs_mjd;                  // Prima osservazione [MJD]
    double last_obs_mjd;                   // Ultima osservazione [MJD]
    
    // Osservazioni con residui calcolati
    std::vector<RWOObservation> observations;
    
    // Propagazione
    int propagation_steps;                 // Passi integratore
    double propagation_time_sec;           // Tempo calcolo [sec]
    
    // Metodo usato
    std::string method;                    // "AstDyn-RKF78", "IOccultCalc-RK4", etc.
    
    // Qualità del fit
    bool is_good_fit() const {
        return (chi2_reduced < 2.0 && rms_total_arcsec < 1.0);
    }
    
    // Report testuale
    std::string toReport() const;
};

/**
 * @class AstDynPropagator
 * @brief Interfaccia al propagatore RKF78 di ITALOccultLibrary
 * 
 * Usa il propagatore ad alta precisione con:
 * - RKF78 (ordine 7, stima errore ordine 8)
 * - Perturbazioni 8 pianeti (Simon et al. 1994)
 * - Perturbazioni AST17 (16 asteroidi maggiori)
 * - Correzioni relativistiche (Schwarzschild)
 * - Step adattivo con controllo automatico
 */
class AstDynPropagator {
public:
    /**
     * @brief Costruttore con tolleranza
     * @param tolerance Tolleranza errore per step adattivo (default: 1e-12)
     */
    explicit AstDynPropagator(double tolerance = 1e-12);
    
    ~AstDynPropagator();
    
    // Configurazione
    void setTolerance(double tol);
    void setStepLimits(double h_min_days, double h_max_days);
    void usePlanetPerturbations(bool enable);
    void useAsteroidPerturbations(bool enable); // AST17
    void useRelativisticCorrections(bool enable);
    
    /**
     * @brief Propaga elementi orbitali da epoca iniziale a finale
     * @param elements Elementi all'epoca iniziale
     * @param target_mjd Epoca target [MJD]
     * @return Elementi propagati all'epoca target
     */
    AstDySElements propagate(const AstDySElements& elements, 
                             double target_mjd);
    
    /**
     * @brief Calcola effemeridi (RA, Dec) per data epoca
     * @param elements Elementi orbitali
     * @param mjd_utc Epoca [MJD UTC]
     * @return Coppia (RA [deg], Dec [deg])
     */
    std::pair<double, double> getRADec(const AstDySElements& elements,
                                        double mjd_utc);
    
    /**
     * @brief Calcola residui O-C per lista di osservazioni
     * @param elements Elementi orbitali
     * @param observations Osservazioni da testare
     * @return Osservazioni con residui calcolati
     */
    std::vector<RWOObservation> computeResiduals(
        const AstDySElements& elements,
        const std::vector<RWOObservation>& observations);
    
    // Statistiche ultima propagazione
    int getLastStepsAccepted() const;
    int getLastStepsRejected() const;
    double getLastMinStep() const;
    double getLastMaxStep() const;
    
private:
    class Impl;
    std::unique_ptr<Impl> pimpl_;
};

/**
 * @class AstDynOrbitFitter
 * @brief Fit orbitale con osservazioni usando differential correction
 * 
 * Implementa:
 * - Propagazione RKF78 ad alta precisione
 * - Calcolo residui O-C
 * - Outlier detection (3σ)
 * - Differential correction (iterativa)
 * - Matrice covarianza elementi
 */
class AstDynOrbitFitter {
public:
    /**
     * @brief Costruttore
     * @param tolerance Tolleranza propagatore (default: 1e-12)
     */
    explicit AstDynOrbitFitter(double tolerance = 1e-12);
    
    ~AstDynOrbitFitter();
    
    // Configurazione
    void setOutlierThreshold(double sigma);          // Default: 3.0
    void setMaxIterations(int max_iter);             // Default: 20
    void setConvergenceTolerance(double tol_au);     // Default: 1e-6 AU
    void setVerbose(bool verbose);
    
    /**
     * @brief Fit orbita con osservazioni
     * @param initial_elements Elementi iniziali (epoch = epoca prima osservazione)
     * @param observations Osservazioni RWO da fittare
     * @return Risultato con elementi migliorati e statistiche
     */
    AstDynOrbitFitResult fit(const AstDySElements& initial_elements,
                       const std::vector<RWOObservation>& observations);
    
    /**
     * @brief Solo calcolo residui (no fitting)
     * @param elements Elementi orbitali
     * @param observations Osservazioni da testare
     * @return Statistiche residui
     */
    AstDynOrbitFitResult computeResidualsOnly(
        const AstDySElements& elements,
        const std::vector<RWOObservation>& observations);
    
private:
    class Impl;
    std::unique_ptr<Impl> pimpl_;
};

/**
 * @brief Client per download dati da AstDyS
 * 
 * AstDySClient è un alias di AstDysClient definito in astdys_client.h
 * per mantenere compatibilità con il codice esistente.
 */
// RIMUOVI la dichiarazione class AstDySClient { ... }; (linee 264-300)

// Alias per compatibilità
using AstDySClient = AstDysClient;

/**
 * @brief Utility per conversione tra formati
 */
namespace astdyn_utils {
    
    // Converti elementi IOccultCalc ↔ AstDyS
    AstDySElements toAstDySElements(const OrbitalElements& elem);
    OrbitalElements fromAstDySElements(const AstDySElements& elem);
    
    // Converti osservazioni IOccultCalc ↔ RWO
    RWOObservation toRWOObservation(const AstrometricObservation& obs);
    AstrometricObservation fromRWOObservation(const RWOObservation& rwo);
    
    // Formato → stringa
    std::string formatResidual(double arcsec);
    std::string formatRMS(double arcsec);
    std::string formatChi2(double chi2, int ndf);
    
    // Parser file locali
    AstDySElements parseEQ1File(const std::string& filename);
    std::vector<RWOObservation> parseRWOFile(const std::string& filename);
    
} // namespace astdyn_utils

} // namespace ioccultcalc

#endif // IOCCULTCALC_ASTDYN_INTERFACE_H
