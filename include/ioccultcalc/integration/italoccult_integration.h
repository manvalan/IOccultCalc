/**
 * @file italoccult_integration.h
 * @brief Integrazione ITALOccultLibrary in IOccultCalc
 * @author IOccultCalc Team
 * @date 1 Dicembre 2025
 * 
 * Questo modulo fornisce l'interfaccia per usare ITALOccultLibrary
 * (wrapper AstDyn con frame conversion) all'interno di IOccultCalc.
 * 
 * Sostituisce/completa orbfit_wrapper.h con propagazione ad alta precisione.
 */

#ifndef IOCCULTCALC_ITALOCCULT_INTEGRATION_H
#define IOCCULTCALC_ITALOCCULT_INTEGRATION_H

#include <italoccultlib/eq1_parser.h>
#include <italoccultlib/orbital_conversions.h>
#include <italoccultlib/astdyn_wrapper.h>
#include <string>
#include <vector>
#include <memory>

namespace ioccultcalc {

/**
 * @struct AsteroidState
 * @brief Stato completo di un asteroide per IOccultCalc
 */
struct AsteroidState {
    std::string name;               ///< Nome/numero asteroide
    double epoch_mjd_tdb;           ///< Epoca [MJD TDB]
    
    // Posizione e velocità ICRF (J2000 equatoriale)
    double x_au, y_au, z_au;        ///< Posizione [AU]
    double vx_au_day, vy_au_day, vz_au_day;  ///< Velocità [AU/day]
    
    // Parametri orbitali (opzionali, per reference)
    double a_au;                    ///< Semiasse maggiore [AU]
    double e;                       ///< Eccentricità
    double i_deg;                   ///< Inclinazione [gradi]
    
    AsteroidState() 
        : epoch_mjd_tdb(0.0), 
          x_au(0.0), y_au(0.0), z_au(0.0),
          vx_au_day(0.0), vy_au_day(0.0), vz_au_day(0.0),
          a_au(0.0), e(0.0), i_deg(0.0) {}
};

/**
 * @class ITALOccultIntegration
 * @brief Classe di integrazione per ITALOccultLibrary in IOccultCalc
 * 
 * Fornisce metodi high-level per:
 * - Caricare asteroidi da file .eq1
 * - Propagare orbite con perturbazioni complete
 * - Convertire stati in formato IOccultCalc
 * 
 * Esempio d'uso:
 * @code
 *   ITALOccultIntegration integrator;
 *   
 *   // Carica asteroide da AstDyS
 *   integrator.loadAsteroidFromEQ1("17030.eq1");
 *   
 *   // Propaga a epoca target
 *   auto state = integrator.propagateToEpoch(61007.0);
 *   
 *   // Usa lo stato per calcoli occultazione
 *   std::cout << "Posizione: " << state.x_au << ", " 
 *             << state.y_au << ", " << state.z_au << " AU\n";
 * @endcode
 */
class ITALOccultIntegration {
public:
    /**
     * @brief Costruttore con configurazione propagazione
     * @param use_high_accuracy Se true, usa configurazione JPL-compliant
     */
    explicit ITALOccultIntegration(bool use_high_accuracy = true);
    
    /**
     * @brief Distruttore
     */
    ~ITALOccultIntegration();
    
    /**
     * @brief Carica asteroide da file .eq1 (formato AstDyS/OrbFit)
     * @param filepath Percorso file .eq1
     * @return true se caricamento riuscito
     */
    bool loadAsteroidFromEQ1(const std::string& filepath);
    
    /**
     * @brief Imposta elementi orbitali manualmente
     * @param name Nome asteroide
     * @param a Semiasse maggiore [AU]
     * @param e Eccentricità
     * @param i Inclinazione [gradi]
     * @param Omega Longitudine nodo ascendente [gradi]
     * @param omega Argomento perielio [gradi]
     * @param M Anomalia media [gradi]
     * @param epoch_mjd_tdb Epoca [MJD TDB]
     */
    void setOrbitalElements(const std::string& name,
                           double a, double e, double i,
                           double Omega, double omega, double M,
                           double epoch_mjd_tdb);
    
    /**
     * @brief Propaga asteroide a epoca target
     * @param target_mjd_tdb Epoca target [MJD TDB]
     * @return Stato dell'asteroide in formato IOccultCalc
     * @throws std::runtime_error Se asteroide non caricato
     */
    AsteroidState propagateToEpoch(double target_mjd_tdb);
    
    /**
     * @brief Propaga a multiple epoche (efficiente per sequenze)
     * @param target_mjd_tdbs Vector di epoche target [MJD TDB]
     * @return Vector di stati corrispondenti
     */
    std::vector<AsteroidState> propagateToEpochs(const std::vector<double>& target_mjd_tdbs);
    
    /**
     * @brief Ottiene epoca corrente elementi
     * @return Epoca [MJD TDB]
     */
    double getCurrentEpoch() const;
    
    /**
     * @brief Ottiene nome asteroide corrente
     * @return Nome
     */
    std::string getAsteroidName() const;
    
    /**
     * @brief Verifica se asteroide caricato
     * @return true se elementi validi disponibili
     */
    bool isInitialized() const;
    
    /**
     * @brief Ottiene informazioni ultima propagazione
     * @return Stringa con statistiche (tempo, step, ecc.)
     */
    std::string getLastPropagationInfo() const;
    
    /**
     * @brief Configura livello accuratezza
     * @param high_accuracy Se true, usa massima precisione (lento)
     */
    void setHighAccuracy(bool high_accuracy);

private:
    std::unique_ptr<AstDynWrapper> wrapper_;
    bool high_accuracy_;
    std::string asteroid_name_;
    double current_epoch_;
    bool initialized_;
    
    /// Converte CartesianStateICRF in AsteroidState
    AsteroidState convertToAsteroidState(const CartesianStateICRF& icrf_state,
                                        const std::string& name) const;
};

/**
 * @brief Funzione helper: carica e propaga in un solo passo
 * @param eq1_file File .eq1 da caricare
 * @param target_mjd_tdb Epoca target [MJD TDB]
 * @param high_accuracy Usa massima precisione (default: true)
 * @return Stato dell'asteroide
 * @throws std::runtime_error In caso di errori
 */
AsteroidState quickPropagateFromEQ1(const std::string& eq1_file,
                                    double target_mjd_tdb,
                                    bool high_accuracy = true);

} // namespace ioccultcalc

#endif // IOCCULTCALC_ITALOCCULT_INTEGRATION_H
