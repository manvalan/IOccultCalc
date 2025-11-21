/**
 * @file force_model.h
 * @brief Modello completo di forze gravitazionali per propagazione orbitale
 * 
 * Include perturbazioni di:
 * - Sole (corpo centrale)
 * - 8 Pianeti (Mercurio - Nettuno)
 * - Luna
 * - Plutone
 * - Maggiori asteroidi: (1) Ceres, (2) Pallas, (4) Vesta
 * - Correzioni relativistiche (opzionali)
 * 
 * Utilizza VSOP87D per posizioni planetarie precise
 */

#ifndef IOCCULTCALC_FORCE_MODEL_H
#define IOCCULTCALC_FORCE_MODEL_H

#include "ioccultcalc/types.h"
#include "ioccultcalc/vsop87.h"
#include <vector>
#include <string>
#include <memory>

namespace ioccultcalc {

/**
 * @enum PerturbingBody
 * @brief Corpi celesti che possono perturbare l'orbita
 */
enum class PerturbingBody {
    SUN,
    MERCURY,
    VENUS,
    EARTH,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    MOON,
    PLUTO,
    CERES,      // (1) Ceres - pianeta nano
    PALLAS,     // (2) Pallas - grande asteroide
    VESTA       // (4) Vesta - grande asteroide
};

/**
 * @struct BodyData
 * @brief Dati fisici di un corpo celeste
 */
struct BodyData {
    PerturbingBody body;
    std::string name;
    double GM;              // Parametro gravitazionale (km³/s²)
    double GM_AU_day;       // Parametro gravitazionale (AU³/day²)
    double radius;          // Raggio medio (km)
    double mass_ratio;      // Massa relativa al Sole
    bool isPlanet;
    
    BodyData();
    BodyData(PerturbingBody b, const std::string& n, double gm, double rad);
    
    // Converte GM da km³/s² a AU³/day²
    static double convertGM_to_AU_day(double GM_km_s);
};

/**
 * @struct ForceModelConfig
 * @brief Configurazione modello di forze
 */
struct ForceModelConfig {
    // Quali corpi includere
    bool includeSun;              // Sempre true (corpo centrale)
    bool includeMercury;
    bool includeVenus;
    bool includeEarth;
    bool includeMars;
    bool includeJupiter;
    bool includeSaturn;
    bool includeUranus;
    bool includeNeptune;
    bool includeMoon;             // Luna (importante per asteroidi NEA)
    bool includePluto;            // Rilevante per TNO
    bool includeCeres;            // (1) Ceres
    bool includePallas;           // (2) Pallas
    bool includeVesta;            // (4) Vesta
    
    // Correzioni relativistiche
    bool includeRelativisticCorrection;  // Schwarzschild (PN1)
    
    // Ottimizzazione
    double minDistanceForPerturbation;   // AU - ignora corpi troppo lontani
    bool useJ2Oblateness;                // J2 del Sole (trascurabile)
    
    ForceModelConfig();
    
    // Configurazioni predefinite
    static ForceModelConfig fastConfig();      // Solo giganti gassosi
    static ForceModelConfig standardConfig();  // Tutti i pianeti
    static ForceModelConfig highPrecisionConfig(); // Pianeti + Luna + grandi asteroidi
    static ForceModelConfig fullConfig();      // Tutto incluso
};

/**
 * @class ForceModel
 * @brief Modello completo di forze gravitazionali N-body
 * 
 * Calcola accelerazione totale su un corpo minore dovuta a:
 * - Attrazione gravitazionale del Sole (centrale)
 * - Perturbazioni planetarie (tutti i pianeti)
 * - Perturbazioni di grandi asteroidi
 * - Correzioni relativistiche post-Newtoniane
 * 
 * Usa coordinate eliocentriche baricentriche (SSB)
 */
class ForceModel {
public:
    ForceModel();
    explicit ForceModel(const ForceModelConfig& config);
    
    /**
     * @brief Calcola accelerazione gravitazionale totale
     * 
     * @param jd Epoca (Julian Date TDB)
     * @param pos Posizione eliocentrica dell'asteroide (AU)
     * @param vel Velocità eliocentrica dell'asteroide (AU/day)
     * @return Accelerazione totale (AU/day²)
     */
    Vector3D computeAcceleration(double jd, 
                                 const Vector3D& pos,
                                 const Vector3D& vel) const;
    
    /**
     * @brief Calcola accelerazione con breakdown per componente
     * 
     * Utile per diagnostica e analisi contributi
     * 
     * @param jd Epoca
     * @param pos Posizione asteroide
     * @param vel Velocità asteroide
     * @param[out] contributions Mappa body->accelerazione
     * @return Accelerazione totale
     */
    Vector3D computeAccelerationWithBreakdown(
        double jd,
        const Vector3D& pos,
        const Vector3D& vel,
        std::map<PerturbingBody, Vector3D>& contributions) const;
    
    /**
     * @brief Calcola accelerazione di un singolo corpo perturbatore
     * 
     * @param body Corpo perturbatore
     * @param jd Epoca
     * @param asteroidPos Posizione asteroide (AU, eliocentrica)
     * @return Accelerazione dovuta a body (AU/day²)
     */
    Vector3D computeSingleBodyAcceleration(PerturbingBody body,
                                           double jd,
                                           const Vector3D& asteroidPos) const;
    
    /**
     * @brief Ottiene posizione di un corpo perturbatore
     * 
     * @param body Corpo
     * @param jd Epoca
     * @return Posizione eliocentrica (AU)
     */
    Vector3D getBodyPosition(PerturbingBody body, double jd) const;
    
    /**
     * @brief Ottiene velocità di un corpo perturbatore
     * 
     * @param body Corpo
     * @param jd Epoca
     * @return Velocità eliocentrica (AU/day)
     */
    Vector3D getBodyVelocity(PerturbingBody body, double jd) const;
    
    /**
     * @brief Ottiene dati fisici di un corpo
     * 
     * @param body Corpo
     * @return Dati fisici (GM, raggio, etc.)
     */
    BodyData getBodyData(PerturbingBody body) const;
    
    /**
     * @brief Imposta configurazione
     */
    void setConfiguration(const ForceModelConfig& config);
    
    /**
     * @brief Ottiene configurazione corrente
     */
    const ForceModelConfig& getConfiguration() const { return config_; }
    
    /**
     * @brief Stima magnitudine perturbazione di un corpo
     * 
     * Utile per decidere quali corpi includere
     * 
     * @param body Corpo
     * @param jd Epoca
     * @param asteroidPos Posizione asteroide
     * @return Magnitudine accelerazione (AU/day²)
     */
    double estimatePerturbationMagnitude(PerturbingBody body,
                                         double jd,
                                         const Vector3D& asteroidPos) const;
    
    /**
     * @brief Lista corpi attivi nella configurazione corrente
     * 
     * @return Vector di corpi che vengono inclusi nel calcolo
     */
    std::vector<PerturbingBody> getActiveBodies() const;
    
    /**
     * @brief Calcola correzione relativistica post-Newtoniana (PN1)
     * 
     * Termine Schwarzschild dovuto al campo gravitazionale solare
     * 
     * @param pos Posizione (AU)
     * @param vel Velocità (AU/day)
     * @return Accelerazione relativistica (AU/day²)
     */
    Vector3D computeRelativisticCorrection(const Vector3D& pos,
                                          const Vector3D& vel) const;
    
    /**
     * @brief Verifica se cache VSOP87 è inizializzata
     */
    bool isCacheInitialized() const;
    
    /**
     * @brief Pre-carica posizioni planetarie per intervallo
     * 
     * Ottimizzazione: pre-calcola VSOP87 per evitare ricalcoli
     * 
     * @param jdStart Inizio intervallo
     * @param jdEnd Fine intervallo
     * @param stepDays Step di cache (default: 1 giorno)
     */
    void precacheEphemerides(double jdStart, double jdEnd, double stepDays = 1.0);
    
    /**
     * @brief Pulisce cache
     */
    void clearCache();
    
private:
    ForceModelConfig config_;
    
    // VSOP87 engine per posizioni planetarie
    mutable std::unique_ptr<VSOP87Calculator> vsop87_;
    mutable std::unique_ptr<ELP2000Calculator> elp2000_;
    
    // Cache per posizioni planetarie (ottimizzazione)
    mutable std::map<std::pair<PerturbingBody, double>, Vector3D> positionCache_;
    mutable std::map<std::pair<PerturbingBody, double>, Vector3D> velocityCache_;
    
    // Database parametri fisici
    void initializeBodyDatabase();
    std::map<PerturbingBody, BodyData> bodyDatabase_;
    
    // Metodi di supporto
    Vector3D computePlanetaryAcceleration(PerturbingBody body,
                                         double jd,
                                         const Vector3D& asteroidPos) const;
    
    Vector3D computeAsteroidAcceleration(PerturbingBody body,
                                        double jd,
                                        const Vector3D& asteroidPos) const;
    
    // Funzione di accelerazione punto-massa
    static Vector3D pointMassAcceleration(const Vector3D& asteroidPos,
                                         const Vector3D& bodyPos,
                                         double GM);
    
    // Verifica se corpo è sufficientemente vicino per perturbazione significativa
    bool isBodyRelevant(const Vector3D& asteroidPos,
                       const Vector3D& bodyPos,
                       double GM) const;
};

/**
 * @class ForceModelAnalyzer
 * @brief Analizza contributi delle varie perturbazioni
 * 
 * Utility per diagnostica e ottimizzazione modello di forze
 */
class ForceModelAnalyzer {
public:
    struct PerturbationAnalysis {
        PerturbingBody body;
        Vector3D acceleration;
        double magnitude;          // Magnitudine (AU/day²)
        double relativeContribution; // % rispetto a perturbazione totale
        double distance;           // Distanza dal corpo (AU)
    };
    
    /**
     * @brief Analizza tutte le perturbazioni per un dato stato
     * 
     * @param model Modello di forze
     * @param jd Epoca
     * @param pos Posizione asteroide
     * @param vel Velocità asteroide
     * @return Vector di analisi per ogni corpo
     */
    static std::vector<PerturbationAnalysis> analyze(
        const ForceModel& model,
        double jd,
        const Vector3D& pos,
        const Vector3D& vel);
    
    /**
     * @brief Stampa report testuale analisi
     */
    static std::string generateReport(const std::vector<PerturbationAnalysis>& analysis);
    
    /**
     * @brief Stima errore accumulato se un corpo viene omesso
     * 
     * @param body Corpo da omettere
     * @param jd Epoca iniziale
     * @param pos Posizione iniziale
     * @param propagationDays Giorni di propagazione
     * @return Errore stimato in posizione (km)
     */
    static double estimateOmissionError(PerturbingBody body,
                                       double jd,
                                       const Vector3D& pos,
                                       double propagationDays);
};

// Costanti fisiche
namespace PhysicalConstants {
    // Parametri gravitazionali (km³/s²) - JPL DE441
    constexpr double GM_SUN      = 1.32712440041279419e11;
    constexpr double GM_MERCURY  = 2.2031868551e4;
    constexpr double GM_VENUS    = 3.2485859882e5;
    constexpr double GM_EARTH    = 3.9860043543e5;
    constexpr double GM_MOON     = 4.9028000661e3;
    constexpr double GM_MARS     = 4.282837362e4;
    constexpr double GM_JUPITER  = 1.26712764100000e8;
    constexpr double GM_SATURN   = 3.79405852000e7;
    constexpr double GM_URANUS   = 5.794556400e6;
    constexpr double GM_NEPTUNE  = 6.836527100580e6;
    constexpr double GM_PLUTO    = 8.696138177608748e2;
    
    // Grandi asteroidi (km³/s²)
    constexpr double GM_CERES    = 62.6284;      // (1) Ceres
    constexpr double GM_PALLAS   = 14.3;         // (2) Pallas  
    constexpr double GM_VESTA    = 17.8;         // (4) Vesta
    
    // Sistema Terra-Luna
    constexpr double GM_EARTH_MOON = GM_EARTH + GM_MOON;
    
    // Velocità della luce (km/s)
    constexpr double SPEED_OF_LIGHT = 299792.458;
}

} // namespace ioccultcalc

#endif // IOCCULTCALC_FORCE_MODEL_H
