#ifndef IOCCULTCALC_CHEBYSHEV_DETECTOR_H
#define IOCCULTCALC_CHEBYSHEV_DETECTOR_H

#include "types.h"
#include "ephemeris.h"
#include "chebyshev_approximation.h"
#include <vector>
#include <memory>
#include <functional>

namespace ioccultcalc {

/**
 * Polinomio Chebyshev 1D con supporto per derivata e ricerca radici
 */
class ChebyshevPolynomial {
public:
    ChebyshevPolynomial() : m_t0(0), m_t1(1) {}
    ChebyshevPolynomial(double t0, double t1, const std::vector<double>& coeffs);
    
    /**
     * Valuta polinomio al tempo t (in JD)
     */
    double evaluate(double t) const;
    
    /**
     * Calcola derivata (restituisce nuovo polinomio)
     * Usa: T'ₙ(x) = n × U_{n-1}(x) dove Uₙ è Chebyshev 2a specie
     * Ma più semplice: d/dt p(t) con cambio scala
     */
    ChebyshevPolynomial derivative() const;
    
    /**
     * Trova radici del polinomio in [t0, t1]
     * Usa metodo companion matrix + eigenvalues
     */
    std::vector<double> findRoots() const;
    
    /**
     * Trova minimi locali in [t0, t1]
     * Calcola radici della derivata e verifica segno derivata seconda
     */
    std::vector<double> findMinima() const;
    
    /**
     * Sottrae costante dal polinomio
     */
    ChebyshevPolynomial subtract(double value) const;
    
    /**
     * Moltiplica per altro polinomio
     */
    ChebyshevPolynomial multiply(const ChebyshevPolynomial& other) const;
    
    /**
     * Somma con altro polinomio
     */
    ChebyshevPolynomial add(const ChebyshevPolynomial& other) const;
    
    /**
     * Intervallo temporale
     */
    double t0() const { return m_t0; }
    double t1() const { return m_t1; }
    
    /**
     * Coefficienti
     */
    const std::vector<double>& coeffs() const { return m_coeffs; }
    int order() const { return m_coeffs.size() - 1; }
    
private:
    double m_t0, m_t1;           // Intervallo temporale [JD]
    std::vector<double> m_coeffs; // Coefficienti Chebyshev
    
    /**
     * Normalizza tempo t ∈ [t0, t1] → x ∈ [-1, 1]
     */
    double normalizeTime(double t) const;
    
    /**
     * De-normalizza x ∈ [-1, 1] → t ∈ [t0, t1]
     */
    double denormalizeTime(double x) const;
    
    /**
     * Valuta polinomio in x normalizzato
     */
    double evaluateNormalized(double x) const;
};

/**
 * Segmento Chebyshev per RA/Dec (coordinate angolari)
 */
class ChebyshevRADecSegment {
public:
    ChebyshevRADecSegment(double t0, double t1, int order);
    
    /**
     * Calcola coefficienti da ephemeris
     */
    void computeCoefficients(Ephemeris& ephemeris);
    
    /**
     * Valuta RA, Dec al tempo t
     */
    void evaluate(double t, double& ra, double& dec) const;
    
    /**
     * Restituisce polinomi RA e Dec
     */
    const ChebyshevPolynomial& raPolynomial() const { return m_raCheb; }
    const ChebyshevPolynomial& decPolynomial() const { return m_decCheb; }
    
    /**
     * Intervallo
     */
    double t0() const { return m_t0; }
    double t1() const { return m_t1; }
    bool contains(double t) const { return t >= m_t0 && t <= m_t1; }
    
private:
    double m_t0, m_t1;
    int m_order;
    ChebyshevPolynomial m_raCheb;
    ChebyshevPolynomial m_decCheb;
};

/**
 * Candidato occultazione trovato da detector
 */
struct OccultationCandidate {
    double jd;              // Tempo minima distanza angolare
    double minDistArcsec;   // Distanza minima in arcsec
    double starRA;          // RA stella (rad)
    double starDec;         // Dec stella (rad)
    int starIndex;          // Indice stella nel catalogo
};

/**
 * Detector occultazioni basato su intersezione Chebyshev
 * 
 * Algoritmo:
 * 1. Interpola traiettoria asteroide (RA, Dec) con Chebyshev
 * 2. Per ogni stella S, calcola polinomio D²(t) = (RA(t)-RA_s)² + (Dec(t)-Dec_s)²
 * 3. Trova minimi di D²(t) (radici di dD²/dt)
 * 4. Se min(D²) < threshold² → candidato occultazione
 */
class ChebyshevOccultationDetector {
public:
    struct Config {
        int order = 11;                    // Ordine polinomio Chebyshev
        double segmentDays = 1.0;          // Durata segmento (giorni)
        double thresholdArcsec = 300.0;    // Threshold ricerca (arcsec) - ~5 arcmin
        double refinementArcsec = 60.0;    // Threshold per refine (arcsec)
        bool verbose = false;
    };
    
    explicit ChebyshevOccultationDetector(const Config& config);
    
    /**
     * Inizializza detector per intervallo temporale
     * Pre-calcola segmenti Chebyshev per l'asteroide
     */
    void initialize(Ephemeris& ephemeris, double startJD, double endJD);
    
    /**
     * Cerca occultazioni con lista di stelle
     * Restituisce candidati ordinati per tempo
     */
    std::vector<OccultationCandidate> findCandidates(
        const std::vector<std::pair<double, double>>& stars  // (RA, Dec) in radianti
    ) const;
    
    /**
     * Cerca occultazione con singola stella
     * Restituisce tutti i passaggi entro threshold
     */
    std::vector<OccultationCandidate> findCandidatesForStar(
        double starRA, double starDec, int starIndex
    ) const;
    
    /**
     * Valuta distanza angolare asteroide-stella al tempo t
     */
    double angularDistance(double t, double starRA, double starDec) const;
    
    /**
     * Statistiche
     */
    int segmentCount() const { return m_segments.size(); }
    bool isInitialized() const { return m_initialized; }
    
private:
    Config m_config;
    std::vector<std::shared_ptr<ChebyshevRADecSegment>> m_segments;
    double m_startJD, m_endJD;
    bool m_initialized;
    
    /**
     * Trova segmento contenente tempo t
     */
    std::shared_ptr<ChebyshevRADecSegment> findSegment(double t) const;
    
    /**
     * Calcola polinomio distanza² da stella per un segmento
     * D²(t) = (RA(t) - RA_s)² × cos²(Dec_avg) + (Dec(t) - Dec_s)²
     */
    ChebyshevPolynomial computeDistanceSquared(
        const ChebyshevRADecSegment& segment,
        double starRA, double starDec
    ) const;
    
    /**
     * Refine tempo minimo usando Newton-Raphson
     */
    double refineMinimum(double t_approx, double starRA, double starDec) const;
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_CHEBYSHEV_DETECTOR_H
