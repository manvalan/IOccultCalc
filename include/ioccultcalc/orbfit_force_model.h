/**
 * @file orbfit_force_model.h
 * @brief Modello di forze OrbFit-compatible per propagazione orbitale ad alta precisione
 * 
 * Implementazione fedele del modello di forze di OrbFit (Milani et al.)
 * Riferimento: http://adams.dm.unipi.it/orbfit/
 * 
 * Componenti del modello:
 * 1. Attrazione solare centrale (Kepler)
 * 2. Perturbazioni planetarie (diretto + indiretto)
 * 3. Perturbazioni asteroidi massicci (Ceres, Pallas, Vesta, etc.)
 * 4. Relatività generale EIH (Einstein-Infeld-Hoffmann) - eq. (4-26) DESCANSO2
 * 5. J2 del Sole (oblatezza)
 * 
 * Tutte le coordinate sono ECLITTICHE J2000.0 (come OrbFit)
 * Unità: AU, giorni, masse solari
 */

#ifndef IOCCULTCALC_ORBFIT_FORCE_MODEL_H
#define IOCCULTCALC_ORBFIT_FORCE_MODEL_H

#include "types.h"
#include "jpl_ephemeris.h"
#include <vector>
#include <array>
#include <string>
#include <memory>

namespace ioccultcalc {

/**
 * @brief Costanti fisiche fondamentali (IAU 2015 / DE441)
 * 
 * Valori identici a quelli usati da OrbFit per massima compatibilità
 */
namespace OrbFitConstants {
    // Velocità della luce in AU/day
    constexpr double VLIGHT = 173.14463267424034;  // c in AU/day (299792.458 km/s)
    
    // Masse in unità solari (GM/GM_Sun) - da DE441
    constexpr double GM0 = 2.959122082855911e-4;   // GM_Sun in AU³/day²
    
    // Rapporti massa/Sole (da jpllis2 in OrbFit)
    constexpr double MASS_MERCURY = 1.6601367952719304e-7;
    constexpr double MASS_VENUS   = 2.4478383396645447e-6;
    constexpr double MASS_EARTH   = 3.0034896149157645e-6;  // Terra sola
    constexpr double MASS_MOON    = 3.6943433937890996e-8;
    constexpr double MASS_EMB     = 3.0404234083651406e-6;  // Terra+Luna
    constexpr double MASS_MARS    = 3.2271514450538743e-7;
    constexpr double MASS_JUPITER = 9.5479193599427834e-4;
    constexpr double MASS_SATURN  = 2.8588598066820652e-4;
    constexpr double MASS_URANUS  = 4.3662440433515473e-5;
    constexpr double MASS_NEPTUNE = 5.1513890204661145e-5;
    constexpr double MASS_PLUTO   = 7.3505989612085371e-9;
    
    // Earth/Moon mass ratio
    constexpr double EMRAT = 81.30056907419062;
    
    // Raggio del Sole in AU (per J2)
    constexpr double RADSUN = 4.6527174e-3;  // 695700 km / AU
    
    // J2 del Sole
    constexpr double SOLJ2 = 2.0e-7;
    
    // Asteroidi massicci - masse in unità solari (da CPV.abe)
    // Lista OrbFit standard: 17 asteroidi (AST17)
    struct MassiveAsteroid {
        int number;
        const char* name;
        double mass;  // GM in unità solari
    };
    
    constexpr MassiveAsteroid MASSIVE_ASTEROIDS[] = {
        {1,   "Ceres",        4.72e-10},
        {2,   "Pallas",       1.03e-10},
        {4,   "Vesta",        1.35e-10},
        {10,  "Hygiea",       4.45e-11},
        {704, "Interamnia",   1.90e-11},
        {511, "Davida",       1.80e-11},
        {65,  "Cybele",       1.78e-11},
        {52,  "Europa",       1.59e-11},
        {87,  "Sylvia",       1.48e-11},
        {31,  "Euphrosyne",   1.39e-11},
        {15,  "Eunomia",      1.31e-11},
        {3,   "Juno",         1.14e-11},
        {16,  "Psyche",       1.13e-11},
        {48,  "Doris",        1.12e-11},
        {29,  "Amphitrite",   1.00e-11},
        {324, "Bamberga",     9.5e-12},
        {7,   "Iris",         8.6e-12}
    };
    constexpr int NUM_MASSIVE_ASTEROIDS = 17;
}

/**
 * @brief Opzioni per il modello di forze OrbFit
 */
struct OrbFitForceOptions {
    // Pianeti da includere
    bool includeMercury = true;
    bool includeVenus = true;
    bool includeEarth = true;
    bool includeMars = true;
    bool includeJupiter = true;
    bool includeSaturn = true;
    bool includeUranus = true;
    bool includeNeptune = true;
    bool includeMoon = true;
    bool includePluto = false;
    
    // Asteroidi massicci
    bool includeMassiveAsteroids = true;
    int numMassiveAsteroids = 17;  // Quanti dei 17 standard
    
    // Correzioni relativistiche (icrel in OrbFit)
    // 0 = nessuna, 1 = solo Sole (genrel), 2 = EIH completo (eihrel2)
    int relativityLevel = 2;
    
    // J2 del Sole (solo se icrel >= 2)
    bool includeJ2Sun = true;
    
    // Self-perturbation: escludere asteroide target dalla lista perturbatori
    int excludeAsteroidNumber = 0;  // 0 = nessuna esclusione
    
    // Costruttori helper
    static OrbFitForceOptions standardConfig() {
        return OrbFitForceOptions();  // Default = OrbFit standard
    }
    
    static OrbFitForceOptions fastConfig() {
        OrbFitForceOptions opt;
        opt.includeMercury = false;
        opt.includeMars = false;
        opt.includeUranus = false;
        opt.includeNeptune = false;
        opt.includeMoon = false;
        opt.includePluto = false;
        opt.includeMassiveAsteroids = false;
        opt.relativityLevel = 1;  // Solo Schwarzschild solare
        opt.includeJ2Sun = false;
        return opt;
    }
    
    static OrbFitForceOptions highPrecisionConfig() {
        OrbFitForceOptions opt;
        opt.includePluto = true;
        opt.relativityLevel = 2;
        opt.includeJ2Sun = true;
        return opt;
    }
};

/**
 * @brief Stato planetario per calcolo forze
 */
struct PlanetaryState {
    Vector3D position;   // Posizione eliocentrica eclittica [AU]
    Vector3D velocity;   // Velocità eliocentrica eclittica [AU/day]
    double distance;     // Distanza dal Sole [AU]
    double GM;           // Parametro gravitazionale [AU³/day²]
};

/**
 * @brief Modello di forze OrbFit-compatible
 * 
 * Implementazione fedele del modello force_model.f90 di OrbFit
 * per propagazione orbitale ad alta precisione.
 */
class OrbFitForceModel {
public:
    /**
     * @brief Costruttore
     * @param options Configurazione del modello
     */
    explicit OrbFitForceModel(const OrbFitForceOptions& options = OrbFitForceOptions());
    
    /**
     * @brief Distruttore
     */
    ~OrbFitForceModel();
    
    /**
     * @brief Inizializza effemeridi JPL
     * @param ephemPath Percorso file DE (auto-download se vuoto)
     * @return true se inizializzazione riuscita
     */
    bool initializeEphemeris(const std::string& ephemPath = "");
    
    /**
     * @brief Calcola accelerazione totale su asteroide
     * 
     * Implementa la funzione FORCE() di OrbFit
     * 
     * @param pos Posizione eliocentrica eclittica dell'asteroide [AU]
     * @param vel Velocità eliocentrica eclittica dell'asteroide [AU/day]
     * @param jd Data giuliana (MJD in OrbFit, ma convertiamo internamente)
     * @return Accelerazione [AU/day²]
     */
    Vector3D computeAcceleration(const Vector3D& pos, 
                                 const Vector3D& vel, 
                                 double jd) const;
    
    /**
     * @brief Calcola accelerazione con breakdown per componente
     * @param[out] sunAccel Accelerazione solare
     * @param[out] planetAccel Somma perturbazioni planetarie
     * @param[out] asteroidAccel Somma perturbazioni asteroidi
     * @param[out] relativisticAccel Correzione relativistica
     * @return Accelerazione totale
     */
    Vector3D computeAccelerationDetailed(const Vector3D& pos,
                                         const Vector3D& vel,
                                         double jd,
                                         Vector3D& sunAccel,
                                         Vector3D& planetAccel,
                                         Vector3D& asteroidAccel,
                                         Vector3D& relativisticAccel) const;
    
    /**
     * @brief Ottiene posizione/velocità pianeta
     * @param planetIndex Indice pianeta (1=Mercurio, ..., 9=Plutone)
     * @param jd Data giuliana
     * @return Stato planetario (eliocentrico eclittico)
     */
    PlanetaryState getPlanetState(int planetIndex, double jd) const;
    
    /**
     * @brief Imposta opzioni
     */
    void setOptions(const OrbFitForceOptions& options);
    
    /**
     * @brief Ottiene opzioni correnti
     */
    const OrbFitForceOptions& getOptions() const { return options_; }
    
    /**
     * @brief Verifica se effemeridi sono inizializzate
     */
    bool isInitialized() const { return initialized_; }

private:
    OrbFitForceOptions options_;
    bool initialized_;
    
    // Effemeridi JPL
    std::unique_ptr<JPLEphemerisReader> jplReader_;
    
    // Cache posizioni planetarie
    mutable std::vector<PlanetaryState> planetCache_;
    mutable double cachedJD_;
    
    // Numero di corpi massicci attivi
    int numPlanets_;  // npla in OrbFit
    int numMass_;     // nmass = npla + iatrue
    
    // Masse GM in AU³/day² (gm() array in OrbFit)
    std::vector<double> GM_;
    
    // === Metodi privati - implementazione OrbFit ===
    
    /**
     * @brief Calcola posizioni planetarie (PLANAST in OrbFit)
     */
    void computePlanetPositions(double jd) const;
    
    /**
     * @brief Accelerazione solare centrale (Kepler)
     */
    Vector3D solarAcceleration(const Vector3D& pos) const;
    
    /**
     * @brief Accelerazione planetaria diretta + indiretta
     * 
     * f += GM[i] * (r_pla - r_ast) / |r_pla - r_ast|³  [diretto]
     * f -= GM[i] * r_pla / |r_pla|³                    [indiretto]
     */
    Vector3D planetaryAcceleration(const Vector3D& pos, double jd) const;
    
    /**
     * @brief Accelerazione da asteroidi massicci
     */
    Vector3D asteroidAcceleration(const Vector3D& pos, double jd) const;
    
    /**
     * @brief Correzione relativistica Schwarzschild (GENREL in OrbFit)
     * Solo termine solare, ordine PN1
     */
    Vector3D genrel(const Vector3D& pos, const Vector3D& vel) const;
    
    /**
     * @brief Correzione relativistica EIH completa (EIHREL2 in OrbFit)
     * Eq. (4-26) da DESCANSO Monograph No. 2
     * Include contributi di tutti i pianeti
     */
    Vector3D eihrel2(const Vector3D& pos, 
                     const Vector3D& vel,
                     double jd) const;
    
    /**
     * @brief J2 del Sole (J2SUN in OrbFit)
     */
    Vector3D j2sun(const Vector3D& pos) const;
    
    /**
     * @brief Converte da equatoriale a eclittico
     */
    static Vector3D equatorialToEcliptic(const Vector3D& eq);
    
    /**
     * @brief Converte da eclittico a equatoriale
     */
    static Vector3D eclipticToEquatorial(const Vector3D& ecl);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_ORBFIT_FORCE_MODEL_H
