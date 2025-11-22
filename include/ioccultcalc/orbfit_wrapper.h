/**
 * @file orbfit_wrapper.h
 * @brief Wrapper C++ per integrare le librerie Fortran di OrbFit
 * 
 * Questo modulo fornisce un'interfaccia C++ moderna per utilizzare
 * gli algoritmi di orbit determination e propagazione di OrbFit.
 * 
 * OrbFit è sviluppato da Andrea Milani et al. (Università di Pisa)
 * Ref: https://adams.dm.unipi.it/orbfit/
 */

#ifndef IOCCULTCALC_ORBFIT_WRAPPER_H
#define IOCCULTCALC_ORBFIT_WRAPPER_H

#include "ioccultcalc/observation.h"
#include "ioccultcalc/orbital_elements.h"
#include <vector>
#include <string>

namespace ioccultcalc {

/**
 * @brief Wrapper per OrbFit - metodi di orbit determination
 * 
 * Questa classe fornisce accesso ai metodi di preliminary orbit determination
 * implementati in OrbFit (Gauss, Vaisala, etc.) tramite interfaccia C++.
 */
class OrbFitWrapper {
public:
    /**
     * @brief Risultato di orbit determination con OrbFit
     */
    struct Solution {
        OrbitalElements elements;
        double topocentric_distance; // Distanza topocentric (AU)
        double rms_residual;         // RMS residui (arcsec)
        int root_index;              // Indice radice polinomio (per Gauss)
        bool is_valid;               // Soluzione fisica valida
    };

    /**
     * @brief Metodo di Gauss di OrbFit
     * 
     * Implementa il metodo di Gauss con solver polinomiale di 8° grado
     * ottimizzato di OrbFit. Più robusto della nostra implementazione.
     * 
     * @param obs1 Prima osservazione
     * @param obs2 Seconda osservazione  
     * @param obs3 Terza osservazione
     * @param debug Stampa informazioni di debug
     * @return Vector di soluzioni (può ritornare multiple radici)
     */
    static std::vector<Solution> gauss(
        const AstrometricObservation& obs1,
        const AstrometricObservation& obs2,
        const AstrometricObservation& obs3,
        bool debug = false
    );

    /**
     * @brief Metodo di Vaisala di OrbFit
     * 
     * Algoritmo alternativo a Gauss, spesso più stabile per NEO
     * e oggetti con geometria difficile.
     * 
     * @param obs1 Prima osservazione
     * @param obs2 Seconda osservazione
     * @param heliocentric_distance Distanza eliocentrica iniziale (AU)
     * @param debug Stampa informazioni di debug
     * @return Soluzione singola
     */
    static Solution vaisala(
        const AstrometricObservation& obs1,
        const AstrometricObservation& obs2,
        double heliocentric_distance,
        bool debug = false
    );

    /**
     * @brief Propagazione orbitale - DISABLED
     * 
     * NOTE: propag() requires orbit_elem type which is too complex to wrap.
     * Use native RA15Integrator or TwoBodyPropagator for orbit propagation.
     * OrbFit wrapper only provides preliminary orbit determination.
     */
    // static bool propagate(...); // DISABLED

    /**
     * @brief Differential corrector - NOT YET IMPLEMENTED
     * 
     * Requires complex file I/O and configuration setup.
     * Use native differential_correction() for now.
     */
    // static bool fitObservations(...); // NOT IMPLEMENTED

    /**
     * @brief Inizializza OrbFit (carica dati, ephemerides, etc.)
     * 
     * Deve essere chiamato una volta prima di usare altri metodi.
     * Carica:
     * - JPL ephemerides
     * - Observatory codes (OBSCODE.dat)
     * - Earth orientation parameters
     * - Asteroid masses
     * 
     * @param orbfit_lib_path Path alla directory lib/ di OrbFit
     * @return true se successo
     */
    static bool initialize(const std::string& orbfit_lib_path);

    /**
     * @brief Cleanup risorse OrbFit
     */
    static void cleanup();

private:
    static bool s_initialized;
    static std::string s_lib_path;
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_ORBFIT_WRAPPER_H
