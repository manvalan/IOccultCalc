#ifndef IOCCULTCALC_CHEBYSHEV_APPROXIMATION_H
#define IOCCULTCALC_CHEBYSHEV_APPROXIMATION_H

#include "types.h"
#include "ephemeris.h"
#include <vector>
#include <memory>

namespace ioccultcalc {

/**
 * Configurazione approssimazione Chebyshev
 */
struct ChebyshevConfig {
    bool enabled = false;              // Abilita Chebyshev approximation
    int order = 11;                    // Ordine polinomio (7-15, LinOccult usa 11)
    double intervalDays = 1.0;         // Intervallo per segmento (0.5-2.0 giorni)
    double fineTimestepMinutes = 15.0; // Timestep fine per scanning
    bool precisionCheck = false;       // Valida errore vs integrazione diretta
    double maxErrorKm = 10.0;          // Errore massimo accettabile (km)
};

/**
 * Singolo segmento di approssimazione Chebyshev
 * Rappresenta un polinomio per intervallo temporale [t0, t1]
 */
class ChebyshevSegment {
public:
    ChebyshevSegment(const JulianDate& t0, const JulianDate& t1, int order);
    
    /**
     * Calcola coefficienti Chebyshev per questo segmento
     * Usa punti campionati da ephemeris.compute()
     */
    void computeCoefficients(Ephemeris& ephemeris);
    
    /**
     * Valuta posizione/velocità al tempo t usando interpolazione Chebyshev
     * Restituisce posizione geocentrica in coordinate equatoriali
     */
    EphemerisData evaluate(const JulianDate& t) const;
    
    /**
     * Verifica se tempo t è contenuto in questo segmento
     */
    bool contains(const JulianDate& t) const;
    
    /**
     * Errore stimato dell'approssimazione (km)
     */
    double estimatedError() const { return m_estimatedError; }
    
    /**
     * Intervallo temporale
     */
    double startJD() const { return m_t0.jd; }
    double endJD() const { return m_t1.jd; }
    
private:
    JulianDate m_t0, m_t1;    // Intervallo temporale
    int m_order;              // Ordine polinomio
    double m_estimatedError;  // Errore stimato (km)
    
    // Coefficienti Chebyshev per X, Y, Z (posizione)
    std::vector<double> m_coeffsX;
    std::vector<double> m_coeffsY;
    std::vector<double> m_coeffsZ;
    
    // Coefficienti per VX, VY, VZ (velocità) - opzionale
    std::vector<double> m_coeffsVX;
    std::vector<double> m_coeffsVY;
    std::vector<double> m_coeffsVZ;
    
    /**
     * Normalizza tempo t ∈ [t0, t1] → x ∈ [-1, 1]
     */
    double normalizeTime(const JulianDate& t) const;
    
    /**
     * Valuta polinomio Chebyshev: Σ cᵢ × Tᵢ(x)
     */
    double evaluatePolynomial(const std::vector<double>& coeffs, double x) const;
    
    /**
     * Calcola polinomi Chebyshev T₀, T₁, ..., Tₙ in x
     * Usa ricorsione: Tₙ₊₁(x) = 2x×Tₙ(x) - Tₙ₋₁(x)
     */
    void evaluateChebyshevBasis(double x, std::vector<double>& T) const;
};

/**
 * Cache di segmenti Chebyshev per intervallo temporale esteso
 * Gestisce pre-computazione e lookup efficiente
 */
class ChebyshevCache {
public:
    explicit ChebyshevCache(const ChebyshevConfig& config);
    
    /**
     * Inizializza cache per intervallo temporale
     * Pre-computa tutti i segmenti necessari
     */
    void initialize(Ephemeris& ephemeris, 
                   const JulianDate& startJD,
                   const JulianDate& endJD);
    
    /**
     * Valuta posizione al tempo t usando segmento appropriato
     * Restituisce nullptr se t fuori range
     */
    EphemerisData evaluate(const JulianDate& t) const;
    
    /**
     * Genera timestep fine per scanning
     * Restituisce vettore di JD con risoluzione fineTimestepMinutes
     */
    std::vector<JulianDate> generateFineTimesteps() const;
    
    /**
     * Statistiche cache
     */
    int segmentCount() const { return m_segments.size(); }
    double maxError() const { return m_maxError; }
    double avgError() const { return m_avgError; }
    bool isInitialized() const { return m_initialized; }
    
    /**
     * Configurazione
     */
    const ChebyshevConfig& config() const { return m_config; }
    
private:
    ChebyshevConfig m_config;
    std::vector<std::shared_ptr<ChebyshevSegment>> m_segments;
    JulianDate m_startJD, m_endJD;
    bool m_initialized;
    double m_maxError, m_avgError;
    
    /**
     * Trova segmento contenente tempo t
     * Usa binary search per O(log n)
     */
    std::shared_ptr<ChebyshevSegment> findSegment(const JulianDate& t) const;
};

/**
 * Helper per fitting coefficienti Chebyshev
 * Usa trasformata discreta di Chebyshev
 */
class ChebyshevFitter {
public:
    /**
     * Calcola coefficienti Chebyshev dato set di punti campionati
     * 
     * @param order Ordine polinomio
     * @param t Vettore di tempi normalizzati in [-1, 1]
     * @param y Vettore di valori campionati
     * @return Coefficienti c₀, c₁, ..., cₙ
     */
    static std::vector<double> fit(int order, 
                                   const std::vector<double>& t,
                                   const std::vector<double>& y);
    
    /**
     * Calcola punti di campionamento Chebyshev ottimali in [-1, 1]
     * Usa nodi: xᵢ = cos(π(2i+1)/(2n+2)) per i=0..n
     */
    static std::vector<double> chebyshevNodes(int order);
    
    /**
     * Valuta errore approssimazione confrontando con punti test
     */
    static double evaluateError(const std::vector<double>& coeffs,
                                const std::vector<double>& t_test,
                                const std::vector<double>& y_test);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_CHEBYSHEV_APPROXIMATION_H
