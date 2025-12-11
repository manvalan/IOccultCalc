/**
 * @file phase2_occultation_geometry.h
 * @brief Phase 2: Calcolo geometria precisa occultazione
 * @date 4 Dicembre 2025
 * 
 * OBIETTIVO:
 * ===========
 * Data la lista di candidati da Phase 1, calcola con precisione massima:
 * - Percorso dell'ombra sulla superficie terrestre
 * - Istante esatto del closest approach
 * - Durata massima dell'occultazione
 * - Chord length e position angle
 * - Ellisse di incertezza
 * 
 * STRATEGIA:
 * ===========
 * Per ogni candidato:
 * 1. Crea path ultra-denso attorno al CA (±5 min, step 1-5 sec)
 * 2. Usa propagazione precisa (RKF78 con perturbazioni planetarie)
 * 3. Corregge parallasse e aberrazione
 * 4. Proietta ombra su ellissoide terrestre (WGS84)
 * 5. Calcola parametri osservativi per diverse località
 * 
 * DIFFERENZE DA PHASE 1:
 * ======================
 * - Intervallo temporale: ±5 minuti (vs 24 ore)
 * - Risoluzione: 1-5 secondi (vs 1.5 ore)
 * - Perturbazioni: TUTTE attive (vs disabilitate)
 * - Correzioni: parallasse, aberrazione, proper motion
 * - Output: geometria completa per report
 */

#ifndef PHASE2_OCCULTATION_GEOMETRY_H
#define PHASE2_OCCULTATION_GEOMETRY_H

#include "phase1_candidate_screening.h"
#include <vector>
#include <string>
#include <memory>

// Forward declarations
namespace astdyn {
    namespace propagation {
        struct KeplerianElements;
        class Propagator;
    }
}

// Include completo per usare KeplerianElements nei membri struct
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"

namespace ioc {
namespace gaia {
    class GaiaStar;
}
}

// ═══════════════════════════════════════════════════════════════
// STRUTTURE DATI
// ═══════════════════════════════════════════════════════════════

/**
 * @brief Punto sul percorso dell'ombra sulla Terra
 */
struct ShadowPathPoint {
    double time_mjd_utc;        ///< Istante (MJD UTC)
    double latitude_deg;         ///< Latitudine (WGS84)
    double longitude_deg;        ///< Longitudine (WGS84)
    double speed_km_s;          ///< Velocità dell'ombra (km/s)
    double position_angle_deg;   ///< Position angle (gradi Est da Nord)
};

/**
 * @brief Geometria dell'occultazione per un osservatore
 */
struct ObserverGeometry {
    std::string site_name;           ///< Nome località
    double latitude_deg;              ///< Latitudine osservatore
    double longitude_deg;             ///< Longitudine osservatore
    double elevation_m;               ///< Elevazione sul livello del mare
    
    double time_ca_mjd_utc;          ///< Istante closest approach (UTC)
    double miss_distance_mas;         ///< Miss distance (milliarcsec)
    double max_duration_sec;          ///< Durata massima occultazione
    double sun_altitude_deg;          ///< Altitudine Sole (per visibilità)
    double moon_separation_deg;       ///< Separazione angolare Luna
    double target_altitude_deg;       ///< Altitudine target
    double target_azimuth_deg;        ///< Azimut target
    
    bool is_in_shadow_path;          ///< True se osservatore nel path
    double distance_from_centerline_km; ///< Distanza dalla linea centrale
};

/**
 * @brief Ellisse di incertezza 1-sigma
 */
struct UncertaintyEllipse {
    double semi_major_axis_mas;   ///< Semi-asse maggiore (mas)
    double semi_minor_axis_mas;   ///< Semi-asse minore (mas)
    double position_angle_deg;    ///< Angolo di posizione (deg)
    double time_uncertainty_sec;  ///< Incertezza temporale (sec)
};

/**
 * @brief Evento di occultazione completo
 */
struct OccultationEvent {
    // Identificazione
    uint64_t star_source_id;         ///< Gaia source_id
    std::string asteroid_name;        ///< Nome asteroide
    int asteroid_number;              ///< Numero asteroide
    
    // Dati stella
    double star_ra_deg;               ///< RA stella J2000
    double star_dec_deg;              ///< Dec stella J2000
    double star_magnitude;            ///< Magnitudine G
    double star_pm_ra_mas_yr;        ///< Proper motion RA
    double star_pm_dec_mas_yr;       ///< Proper motion Dec
    
    // Geometria globale
    double time_ca_mjd_utc;          ///< Istante CA geocentrico (UTC)
    double closest_approach_mas;      ///< Minima distanza geocentrica (mas)
    double max_duration_sec;          ///< Durata massima teorica
    double chord_length_km;           ///< Lunghezza chord
    double shadow_width_km;           ///< Larghezza ombra (diametro asteroide proiettato)
    double position_angle_deg;        ///< PA del moto relativo
    
    // Path dell'ombra sulla Terra
    std::vector<ShadowPathPoint> shadow_path;  ///< Traccia sulla superficie
    double path_length_km;            ///< Lunghezza totale percorso
    double path_duration_sec;         ///< Durata totale attraversamento Terra
    
    // Geometria per osservatori specifici
    std::vector<ObserverGeometry> observer_predictions;
    
    // Incertezza
    UncertaintyEllipse uncertainty;
    
    // Dati propagazione
    double asteroid_distance_au;      ///< Distanza asteroide dalla Terra
    double star_distance_au;          ///< Distanza stella (se nota)
    double sun_target_elongation_deg; ///< Elongazione dal Sole
    
    // Verifica occultazione (LinOccult method)
    bool is_occultation;              ///< True se c'è occultazione (shadow_threshold_factor = 3.0)
    double occultation_probability;   ///< Probabilità occultazione (CDF method)
    double shadow_threshold_arcsec;   ///< Soglia shadow (angular_radius + 3*uncertainty)
    double angular_radius_arcsec;      ///< Raggio angolare asteroide
    double total_uncertainty_arcsec;  ///< Incertezza totale (orbitale + stellare)
    
    // Quality flags
    bool high_confidence;             ///< True se geometria affidabile
    double snr;                       ///< Signal-to-noise ratio stima
    std::string notes;                ///< Note aggiuntive
};

// ═══════════════════════════════════════════════════════════════
// CONFIGURAZIONE PHASE 2
// ═══════════════════════════════════════════════════════════════

/**
 * @brief Configurazione calcolo Phase 2
 */
struct Phase2Config {
    // ORBITAL FITTING (nuova funzionalità)
    bool refine_orbit_from_observations = false; ///< Scarica RWO e rifà fit orbitale (DISATTIVATO per ora)
    std::string mpc_code = "";                    ///< Codice MPC asteroide (es: "17030")
    int observation_arc_days = 365;               ///< Arco osservativo per fit (giorni)
    bool use_all_available_observations = false;  ///< Usa tutte osservazioni disponibili
    int max_observations_for_fit = 20;            ///< Numero massimo osservazioni da usare (le più recenti)
    
    // Parametri fitting
    bool fit_planetary_perturbations = true;      ///< Include perturbazioni nel fit
    bool fit_relativistic_effects = true;         ///< Include relatività nel fit
    bool fit_asteroid_perturbations = false;      ///< Perturbazioni da Ceres, Vesta, etc
    double fit_tolerance = 1e-12;                 ///< Tolleranza integrazione fit
    int max_fit_iterations = 50;                  ///< Iterazioni massime fit
    double convergence_threshold = 1e-9;          ///< Soglia convergenza fit
    
    // Parametri temporali
    double time_window_minutes = 5.0;    ///< Finestra attorno CA (±minuti)
    double time_step_seconds = 1.0;      ///< Risoluzione temporale (sec)
    
    // Parametri propagazione finale
    bool use_planetary_perturbations = true;  ///< Abilita perturbazioni
    bool use_relativistic_effects = true;     ///< Effetti relativistici
    double integrator_tolerance = 1e-12;      ///< Tolleranza integrazione
    
    // Correzioni astrometriche
    bool apply_parallax = true;          ///< Correggi parallasse
    bool apply_aberration = true;        ///< Correggi aberrazione
    bool apply_proper_motion = true;     ///< Applica proper motion stella
    double proper_motion_epoch = 2016.0; ///< Epoca proper motion Gaia
    
    // Output
    bool compute_uncertainty = true;     ///< Calcola ellisse incertezza
    bool compute_shadow_path = true;     ///< Calcola path su Terra
    int shadow_path_points = 100;        ///< Numero punti path
    
    // Osservatori
    std::vector<ObserverGeometry> observer_sites;  ///< Località da calcolare
};

/**
 * @brief Risultati fitting orbitale
 */
struct OrbitalFitResults {
    bool fit_successful = false;
    int num_observations_used = 0;
    double rms_residuals_arcsec = 0.0;
    double max_residual_arcsec = 0.0;
    double chi_squared = 0.0;
    int iterations_performed = 0;
    astdyn::propagation::KeplerianElements refined_elements;
    Eigen::MatrixXd covariance_matrix;  ///< Matrice covarianza 6x6
    std::string fit_notes;
};

/**
 * @brief Risultati Phase 2
 */
struct Phase2Results {
    std::vector<OccultationEvent> events;  ///< Eventi calcolati
    int successful_calculations = 0;
    int failed_calculations = 0;
    double total_computation_time_ms = 0.0;
    std::string error_messages;
    
    // Risultati fitting orbitale (se abilitato)
    bool orbit_refined = false;
    OrbitalFitResults orbital_fit;
};

// ═══════════════════════════════════════════════════════════════
// CLASSE PRINCIPALE
// ═══════════════════════════════════════════════════════════════

/**
 * @brief Calcola geometria precisa delle occultazioni
 * 
 * Prende i candidati da Phase1 e calcola:
 * - Geometria dettagliata dell'evento
 * - Path dell'ombra sulla Terra
 * - Predizioni per osservatori specifici
 */
class Phase2OccultationGeometry {
public:
    Phase2OccultationGeometry();
    ~Phase2OccultationGeometry();
    
    // Non copiabile
    Phase2OccultationGeometry(const Phase2OccultationGeometry&) = delete;
    Phase2OccultationGeometry& operator=(const Phase2OccultationGeometry&) = delete;
    
    /**
     * @brief Carica elementi orbitali da file EQ1
     * @param eq1_path Percorso file
     * @return true se successo
     */
    bool loadAsteroidFromEQ1(const std::string& eq1_path);
    
    /**
     * @brief Carica elementi orbitali dal database JSON locale
     * @param asteroid_number Numero dell'asteroide (es: 17030)
     * @param json_path Path al file JSON (default: ~/.ioccultcalc/data/all_numbered_asteroids.json)
     * @return true se caricamento riuscito
     */
    bool loadAsteroidFromJSON(int asteroid_number, const std::string& json_path = "");
    
    /**
     * @brief Imposta elementi orbitali direttamente
     */
    void setOrbitalElements(const astdyn::propagation::KeplerianElements& elements);
    
    /**
     * @brief Ottiene elementi orbitali correnti
     */
    const astdyn::propagation::KeplerianElements& getOrbitalElements() const;
    
    /**
     * @brief Imposta diametro asteroide (km) per calcoli LinOccult
     * @param diameter_km Diametro in km (0 = usa placeholder)
     */
    void setAsteroidDiameter(double diameter_km);
    
    /**
     * @brief Imposta incertezza orbitale (km) per calcoli LinOccult
     * @param uncertainty_km Incertezza orbitale in km (1-sigma)
     */
    void setOrbitalUncertainty(double uncertainty_km);
    
    /**
     * @brief Calcola geometria per lista candidati
     * 
     * @param candidates Lista candidati da Phase 1
     * @param config Configurazione calcolo
     * @return Risultati con eventi calcolati
     */
    Phase2Results calculateGeometry(
        const std::vector<ioccultcalc::CandidateStar>& candidates,
        const Phase2Config& config);
    
    /**
     * @brief Calcola geometria per singolo candidato
     */
    OccultationEvent calculateSingleEvent(
        const ioccultcalc::CandidateStar& candidate,
        const Phase2Config& config);
    
    /**
     * @brief Aggiunge sito osservatorio alle predizioni
     */
    void addObserverSite(const std::string& name, 
                        double lat_deg, double lon_deg, double elev_m);
    
    /**
     * @brief Pulisce lista osservatori
     */
    void clearObserverSites();
    
private:
    class Impl;
    std::unique_ptr<Impl> pimpl_;
};

#endif // PHASE2_OCCULTATION_GEOMETRY_H
