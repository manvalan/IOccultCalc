/**
 * @file chebyshev_occultation_manager.h
 * @brief Manager per ricerca occultazioni con propagazione Chebyshev RKF78
 * @author IOccultCalc Development Team
 * @date 4 Dicembre 2025
 * 
 * Questa classe gestisce l'intero workflow per la ricerca di occultazioni:
 * 1. Carica elementi orbitali da file .eq1 all'epoca t0
 * 2. Propaga con RKF78 ad alta precisione all'epoca del test
 *    (usa ChebyshevRKF78Propagator di AstDyn - contiene già propagatore integrato!)
 * 3. Fitta polinomi di Chebyshev per query veloci
 * 4. Crea corridor path e cerca stelle candidate
 * 5. Calcola closest approach per ogni stella
 * 6. Filtra e ordina le occultazioni candidate
 * 
 * Integra:
 * - AstDyn ChebyshevRKF78Propagator (propagazione RKF78 + fitting integrati)
 * - ChebyshevApproximation (valutazione polinomi)
 * - Corridor API (ricerca stelle)
 * - Closest approach detector
 * 
 * IMPORTANTE: ChebyshevRKF78Propagator contiene GIÀ:
 * - RKF78 integrator (tolerance 1e-12 AU)
 * - Tutte le perturbazioni (8 pianeti + asteroidi + relatività)
 * - Conversione frame ECLM → ICRF automatica
 * - Output: coordinate barycentriche ICRF J2000.0 in AU
 */

#ifndef CHEBYSHEV_OCCULTATION_MANAGER_H
#define CHEBYSHEV_OCCULTATION_MANAGER_H

#include <string>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <eq1_parser.h>

// Forward declarations per tipi non usati nella parte pubblica
namespace ioccultcalc {
    struct CartesianStateICRF;
}

namespace ioccultcalc {

/**
 * @brief Configurazione per la ricerca di occultazioni
 */
struct OccultationSearchConfig {
    // Parametri temporali
    double start_mjd;           ///< Inizio ricerca (MJD TDB)
    double end_mjd;             ///< Fine ricerca (MJD TDB)
    
    // Parametri propagazione
    size_t num_propagation_points; ///< Numero punti RKF78 (default: 100)
    double rkf78_tolerance;     ///< Tolleranza integrazione (default: 1e-12 AU)
    
    // Parametri Chebyshev
    size_t num_chebyshev_coeffs; ///< Coefficienti per asse (default: 8)
    
    // Parametri corridor
    double corridor_width_km;   ///< Larghezza corridor (default: 100 km)
    double search_step_days;    ///< Passo campionamento path (default: 0.1 giorni)
    
    // Filtri stelle
    double max_magnitude;       ///< Magnitudine massima (default: 16.0)
    double min_altitude_deg;    ///< Altitudine minima osservatore (default: 15°)
    
    // Parametri closest approach
    double threshold_arcsec;    ///< Soglia distanza occultazione (default: 2.0")
    
    // Observer location
    double observer_lat_deg;    ///< Latitudine osservatore (default: 40° Nord)
    double observer_lon_deg;    ///< Longitudine osservatore (default: 0° Greenwich)
    double observer_alt_m;      ///< Altitudine osservatore (default: 0 m slm)
    
    /**
     * @brief Costruttore con valori di default
     */
    OccultationSearchConfig();
    
    /**
     * @brief Preset per ricerca veloce (survey)
     */
    static OccultationSearchConfig fastSurvey();
    
    /**
     * @brief Preset per ricerca ad alta precisione
     */
    static OccultationSearchConfig highPrecision();
};

/**
 * @brief Informazioni su una stella candidata
 */
struct StarCandidate {
    std::string catalog_id;     ///< ID catalogo (es: "Gaia DR3 123456789")
    double ra_deg;              ///< RA (gradi, ICRS)
    double dec_deg;             ///< Dec (gradi, ICRS)
    double magnitude;           ///< Magnitudine apparente
    double pmra_mas_yr;         ///< Proper motion RA (mas/year)
    double pmdec_mas_yr;        ///< Proper motion Dec (mas/year)
    
    // Risultati closest approach
    double closest_distance_arcsec; ///< Distanza minima (arcsec)
    double closest_epoch_mjd;       ///< Epoca closest approach (MJD TDB)
    double position_angle_deg;      ///< Angolo posizione (gradi da Nord)
    double relative_velocity_km_s;  ///< Velocità relativa (km/s)
    
    bool is_occultation;        ///< true se distanza < threshold
};

/**
 * @brief Risultati della ricerca di occultazioni
 */
struct OccultationSearchResults {
    std::string asteroid_name;
    double asteroid_diameter_km;
    EquinoctialElements elements;
    
    std::vector<StarCandidate> candidates;  ///< Tutte le stelle nel corridor
    std::vector<StarCandidate> occultations; ///< Solo occultazioni (< threshold)
    
    // Statistiche
    size_t total_stars_searched;
    size_t stars_in_corridor;
    double search_duration_seconds;
    double propagation_time_ms;
    double fitting_time_ms;
    double corridor_query_time_ms;
    double closest_approach_time_ms;
    
    // Accuratezza Chebyshev fitting
    double chebyshev_rms_error_au;
    double chebyshev_max_error_au;
};

/**
 * @brief Manager principale per ricerca occultazioni con Chebyshev
 */
class ChebyshevOccultationManager {
public:
    /**
     * @brief Costruttore
     * @param config Configurazione ricerca
     */
    explicit ChebyshevOccultationManager(const OccultationSearchConfig& config);
    
    /**
     * @brief Distruttore
     */
    ~ChebyshevOccultationManager();
    
    // Disable copy (contiene unique_ptr)
    ChebyshevOccultationManager(const ChebyshevOccultationManager&) = delete;
    ChebyshevOccultationManager& operator=(const ChebyshevOccultationManager&) = delete;
    
    /**
     * @brief Carica asteroide da file .eq1 e inizializza propagatore
     * @param eq1_file Path al file .eq1 (formato AstDyS/OrbFit)
     * @return true se caricamento riuscito
     * 
     * Internamente crea ChebyshevRKF78Propagator che:
     * - Carica elementi orbitali dal file .eq1
     * - Inizializza RKF78 integrator (tolerance 1e-12 AU)
     * - Configura tutte le perturbazioni (8 pianeti + asteroidi + relatività)
     * - Prepara conversione frame ECLM → ICRF
     */
    bool loadAsteroidFromEQ1(const std::string& eq1_file);
    
    /**
     * @brief Propaga orbita con RKF78 e fitta Chebyshev
     * @return true se propagazione e fitting riusciti
     * 
     * Steps:
     * 1. Usa ChebyshevRKF78Propagator::propagateForChebyshev() per ottenere
     *    posizioni ad alta precisione (già con tutte le correzioni applicate!)
     * 2. Fitta polinomi di Chebyshev ai punti propagati con ChebyshevApproximation
     * 3. Calcola errore RMS del fitting
     * 
     * NOTA: Tutti i dati ritornati da propagateForChebyshev() sono:
     * - Frame: ICRF J2000.0 (barycentric)
     * - Unità: AU (posizioni), AU/day (velocità)
     * - Accuratezza: 0.7 km vs JPL Horizons
     */
    bool propagateAndFit();
    
    /**
     * @brief Cerca stelle nel corridor path
     * @return Numero di stelle trovate
     * 
     * Steps:
     * 1. Campiona path dell'asteroide (ogni search_step_days)
     * 2. Per ogni punto, query Gaia nel raggio corridor_width_km
     * 3. Applica filtri (magnitudine, altitudine)
     * 4. Deduplicazione stelle
     */
    size_t searchStarsInCorridor();
    
    /**
     * @brief Calcola closest approach per tutte le stelle
     * @return Numero di occultazioni trovate (distanza < threshold)
     * 
     * Steps:
     * 1. Per ogni stella, trova epoca di closest approach
     * 2. Usa Chebyshev per query veloce (sub-microsecondo)
     * 3. Calcola distanza minima, position angle, velocità
     * 4. Classifica come occultazione se distanza < threshold
     */
    size_t computeClosestApproaches();
    
    /**
     * @brief Esegue ricerca completa (propagate + search + closest approach)
     * @return Risultati completi della ricerca
     */
    OccultationSearchResults performFullSearch();
    
    /**
     * @brief Ottieni posizione asteroide all'epoca specificata
     * @param epoch_mjd Epoca (MJD TDB)
     * @return Posizione baricentrica ICRF (AU)
     * 
     * Usa il fitting Chebyshev per query veloce (deve aver chiamato propagateAndFit())
     */
    Eigen::Vector3d getPositionAtEpoch(double epoch_mjd) const;
    
    /**
     * @brief Ottieni posizione + velocità all'epoca specificata
     * @param epoch_mjd Epoca (MJD TDB)
     * @return State baricentrico ICRF (posizione AU, velocità AU/day)
     */
    CartesianStateICRF getStateAtEpoch(double epoch_mjd) const;
    
    /**
     * @brief Ottieni lista stelle candidate
     */
    const std::vector<StarCandidate>& getCandidates() const { return candidates_; }
    
    /**
     * @brief Ottieni lista occultazioni (solo distanza < threshold)
     */
    std::vector<StarCandidate> getOccultations() const;
    
    /**
     * @brief Ottieni configurazione attuale
     */
    const OccultationSearchConfig& getConfig() const { return config_; }
    
    /**
     * @brief Ottieni nome asteroide caricato
     */
    std::string getAsteroidName() const;
    
    /**
     * @brief Ottieni elementi orbitali caricati
     */
    EquinoctialElements getOrbitalElements() const;
    
    /**
     * @brief Ottieni accuratezza Chebyshev (RMS error)
     * @return RMS error in AU (deve aver chiamato propagateAndFit())
     */
    double getChebyshevRMSError() const;
    
    /**
     * @brief Esporta risultati in formato JSON
     * @param output_file Path file output
     * @return true se esportazione riuscita
     */
    bool exportResultsJSON(const std::string& output_file) const;
    
    /**
     * @brief Esporta risultati in formato OOP (IOTA prediction)
     * @param output_file Path file output
     * @return true se esportazione riuscita
     */
    bool exportResultsOOP(const std::string& output_file) const;
    
    /**
     * @brief Stampa riepilogo risultati su console
     */
    void printSummary() const;

private:
    // Implementation details (PIMPL pattern)
    class Impl;
    std::unique_ptr<Impl> pimpl_;
    
    // Configuration
    OccultationSearchConfig config_;
    
    // Results
    std::vector<StarCandidate> candidates_;
    
    // Internal state
    bool asteroid_loaded_;
    bool propagation_done_;
    bool fitting_done_;
    bool corridor_searched_;
};

} // namespace ioccultcalc

#endif // CHEBYSHEV_OCCULTATION_MANAGER_H
