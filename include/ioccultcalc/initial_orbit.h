#ifndef IOCCULTCALC_INITIAL_ORBIT_H
#define IOCCULTCALC_INITIAL_ORBIT_H

#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/observation.h"
#include <vector>

namespace ioccultcalc {

/**
 * @brief Metodi per determinazione orbitale iniziale
 * 
 * Implementazione dei metodi classici per determinare un'orbita preliminare
 * da osservazioni astrometriche, basata su Find_Orb di Bill Gray.
 */
class InitialOrbit {
public:
    /**
     * @brief Metodo di Gauss per determinazione orbitale da 3 osservazioni
     * 
     * Il metodo di Gauss risolve il problema della determinazione orbitale
     * usando esattamente 3 osservazioni. Può restituire da 0 a 3 soluzioni.
     * 
     * Riferimenti:
     * - Boulet, "Methods of Orbit Determination", cap. 10
     * - Find_Orb: gauss.cpp
     * 
     * @param obs1 Prima osservazione
     * @param obs2 Osservazione centrale (tipicamente al centro dell'arco)
     * @param obs3 Ultima osservazione
     * @param mu Parametro gravitazionale (default: k^2 per sistema solare)
     * @return Vector di soluzioni (0-3 elementi orbitali validi)
     * 
     * Note:
     * - Funziona meglio per archi di qualche settimana a qualche mese
     * - Le 3 osservazioni devono essere non-coplanari (d0 != 0)
     * - Per oggetti eliocentrici solamente (non satelliti)
     */
    static std::vector<OrbitalElements> gauss(
        const AstrometricObservation& obs1,
        const AstrometricObservation& obs2,
        const AstrometricObservation& obs3,
        double mu = 0.01720209895 * 0.01720209895); // GAUSS_K^2
    
    /**
     * @brief Metodo di Herget per determinazione orbitale iterativa
     * 
     * Il metodo di Herget assume distanze R1 e R2 per le prime/ultime
     * osservazioni e le raffina iterativamente per minimizzare i residui.
     * 
     * Riferimenti:
     * - https://www.projectpluto.com/herget.htm
     * - Find_Orb: orb_func.cpp, herget_method()
     * 
     * @param observations Set completo di osservazioni
     * @param r1_guess Guess iniziale per la distanza della prima osservazione (AU)
     * @param r2_guess Guess iniziale per la distanza dell'ultima osservazione (AU)
     * @param max_iterations Numero massimo di iterazioni (default: 5)
     * @return Elementi orbitali se convergenza, altrimenti elementi nulli
     * 
     * Note:
     * - Più robusto di Gauss per archi lunghi
     * - Converge anche con guess iniziali non ottimali
     * - Tipicamente seguito da adjust_herget_results() per raffinamento
     */
    static OrbitalElements herget(
        const std::vector<AstrometricObservation>& observations,
        double r1_guess,
        double r2_guess,
        int max_iterations = 5);
    
    /**
     * @brief Raffina i risultati del metodo di Herget
     * 
     * Dopo il metodo di Herget, questa funzione applica un fit least-squares
     * con distanze fissate per migliorare la soluzione.
     * 
     * @param observations Osservazioni (modificate con distanze calcolate)
     * @param orbit Orbita da raffinare (in input/output)
     * @return 0 se successo, codice errore altrimenti
     */
    static int adjustHergetResults(
        std::vector<AstrometricObservation>& observations,
        OrbitalElements& orbit);
    
    /**
     * @brief Metodo di Väisälä per archi molto brevi
     * 
     * Assume che l'oggetto sia al perielio e richiede un guess della
     * distanza perielica.
     * 
     * Riferimenti:
     * - https://www.projectpluto.com/vaisala.htm
     * - Find_Orb: orb_func.cpp, find_trial_orbit()
     * 
     * @param observations Set di osservazioni (arco breve)
     * @param perihelion_distance Distanza perielica (AU)
     * @return Elementi orbitali se convergenza
     * 
     * Note:
     * - Usato solo per archi di pochi giorni
     * - Richiede guess ragionevole di q (tipicamente 1-3 AU)
     */
    static OrbitalElements vaisala(
        const std::vector<AstrometricObservation>& observations,
        double perihelion_distance);
    
    /**
     * @brief Trova automaticamente l'orbita migliore provando vari metodi
     * 
     * Strategia Auto-Solve di Find_Orb:
     * 1. Prova Gauss (tutte le soluzioni)
     * 2. Se fallisce, prova Herget con vari guess di R1/R2
     * 3. Se fallisce, prova Väisälä con vari guess di q
     * 4. Applica attempt_improvements() alla soluzione migliore
     * 
     * @param observations Set completo di osservazioni
     * @param score Output: score della soluzione (RMS residui)
     * @return Elementi orbitali migliori trovati
     * 
     * Note:
     * - Metodo più robusto, converge nella maggior parte dei casi
     * - Può richiedere diversi secondi per archi difficili
     */
    static OrbitalElements autoSolve(
        const std::vector<AstrometricObservation>& observations,
        double* score = nullptr);
    
    /**
     * @brief Valuta la qualità di un'orbita iniziale
     * 
     * Calcola uno score basato su RMS residui e altre metriche.
     * Usato da auto_solve() per confrontare soluzioni.
     * 
     * @param orbit Orbita da valutare
     * @param observations Osservazioni
     * @param epoch Epoca dell'orbita
     * @return Score (più basso = migliore)
     */
    static double evaluateInitialOrbit(
        const OrbitalElements& orbit,
        const std::vector<AstrometricObservation>& observations,
        const JulianDate& epoch);
    
    /**
     * @brief Lambert solver - Risolve il problema di Lambert (transfer orbit)
     * 
     * Calcola i vettori velocità iniziale e finale per un'orbita che connette
     * due posizioni r1 e r2 in un tempo dato dt.
     * 
     * Implementazione basata su Universal Variables algorithm (Battin, 1999).
     * 
     * Riferimenti:
     * - Battin, "An Introduction to the Mathematics and Methods of Astrodynamics"
     * - Vallado, "Fundamentals of Astrodynamics and Applications", cap. 7.6
     * - Curtis, "Orbital Mechanics for Engineering Students", cap. 5.3
     * 
     * @param r1 Vettore posizione iniziale (AU)
     * @param r2 Vettore posizione finale (AU)
     * @param dt Tempo di volo (giorni)
     * @param mu Parametro gravitazionale (default: k^2)
     * @param prograde True per orbita prograde, false per retrograde
     * @param v1 Output: vettore velocità iniziale (AU/giorno)
     * @param v2 Output: vettore velocità finale (AU/giorno)
     * @return 0 se successo, -1 se fallimento
     * 
     * Note:
     * - Può avere soluzioni multiple (orbite corte/lunghe)
     * - Il parametro prograde seleziona quale soluzione
     * - Converge in 5-15 iterazioni per la maggior parte dei casi
     * - Tollera archi da poche ore a diversi anni
     */
    static int solveLambert(
        const Vector3D& r1,
        const Vector3D& r2,
        double dt,
        double mu,
        bool prograde,
        Vector3D& v1,
        Vector3D& v2);

private:
    // Helper functions (internal)
    
    /**
     * @brief Trova l'orbita di trasferimento tra due osservazioni
     * 
     * Dato un'orbita iniziale, la integra e minimizza la differenza
     * con la posizione osservata finale usando derivate parziali.
     */
    static int findTransferOrbit(
        double* orbit,
        const AstrometricObservation& obs1,
        const AstrometricObservation& obs2,
        bool already_have_approximate_orbit);
    
    /**
     * @brief Setta le distanze per un'osservazione
     */
    static void setDistance(
        AstrometricObservation& obs,
        double distance);
    
    /**
     * @brief Trova il raggio data la distanza solare
     */
    static double findRGivenSolarR(
        const AstrometricObservation& obs,
        double solar_r);
    
    /**
     * @brief Risolve equazione polinomiale (per Gauss)
     */
    static int findRealPolynomialRoots(
        const double* poly,
        int poly_order,
        double* real_roots);
    
    /**
     * @brief Calcola prodotto misto (A × B) · C
     */
    static double crossThenDot(
        const double* a,
        const double* b,
        const double* c);
    
    /**
     * @brief Seleziona le 3 migliori osservazioni per Gauss
     * 
     * Massimizza d0 = (obs1 × obs2) · obs3
     */
    static double findBestObsForGauss(
        const std::vector<AstrometricObservation>& observations,
        unsigned obs_idx[3]);
    
    /**
     * @brief Tenta miglioramenti iterativi su un'orbita
     * 
     * Alterna passi di Herget e full-step fino a convergenza.
     */
    static double attemptImprovements(
        OrbitalElements& orbit,
        std::vector<AstrometricObservation>& observations);
    
    /**
     * @brief Verifica se un'orbita è ragionevole
     */
    static bool isUnreasonableOrbit(const OrbitalElements& orbit);
    
    /**
     * @brief Calcola lo span massimo per metodo di Herget
     */
    static double maxHergetSpan(double r1, double r2);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_INITIAL_ORBIT_H
