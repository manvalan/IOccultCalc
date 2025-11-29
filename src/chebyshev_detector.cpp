#include "ioccultcalc/chebyshev_detector.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace ioccultcalc {

// ============================================================================
// ChebyshevPolynomial Implementation
// ============================================================================

ChebyshevPolynomial::ChebyshevPolynomial(double t0, double t1, 
                                         const std::vector<double>& coeffs)
    : m_t0(t0), m_t1(t1), m_coeffs(coeffs) {
}

double ChebyshevPolynomial::normalizeTime(double t) const {
    return 2.0 * (t - m_t0) / (m_t1 - m_t0) - 1.0;
}

double ChebyshevPolynomial::denormalizeTime(double x) const {
    return m_t0 + (x + 1.0) * (m_t1 - m_t0) / 2.0;
}

double ChebyshevPolynomial::evaluateNormalized(double x) const {
    // Algoritmo di Clenshaw per stabilità numerica
    int n = m_coeffs.size() - 1;
    if (n < 0) return 0.0;
    
    double b_k2 = 0.0;
    double b_k1 = 0.0;
    
    for (int k = n; k >= 1; k--) {
        double b_k = 2.0 * x * b_k1 - b_k2 + m_coeffs[k];
        b_k2 = b_k1;
        b_k1 = b_k;
    }
    
    return x * b_k1 - b_k2 + m_coeffs[0];
}

double ChebyshevPolynomial::evaluate(double t) const {
    double x = normalizeTime(t);
    return evaluateNormalized(x);
}

ChebyshevPolynomial ChebyshevPolynomial::derivative() const {
    // La derivata di un polinomio Chebyshev p(x) = Σ cₖ Tₖ(x)
    // è p'(x) = Σ c'ₖ Tₖ(x) dove:
    // c'ₙ₋₁ = 2n cₙ
    // c'ₖ = c'ₖ₊₂ + 2(k+1) cₖ₊₁  per k = n-2, ..., 0
    
    int n = m_coeffs.size() - 1;
    if (n <= 0) {
        return ChebyshevPolynomial(m_t0, m_t1, {0.0});
    }
    
    std::vector<double> dcoeffs(n, 0.0);
    
    if (n >= 1) {
        dcoeffs[n-1] = 2.0 * n * m_coeffs[n];
    }
    
    for (int k = n - 2; k >= 0; k--) {
        dcoeffs[k] = (k + 2 < n ? dcoeffs[k + 2] : 0.0) + 2.0 * (k + 1) * m_coeffs[k + 1];
    }
    
    // Correggere per il cambio di scala: d/dt = d/dx × dx/dt = d/dx × 2/(t1-t0)
    double scale = 2.0 / (m_t1 - m_t0);
    for (double& c : dcoeffs) {
        c *= scale;
    }
    
    return ChebyshevPolynomial(m_t0, m_t1, dcoeffs);
}

std::vector<double> ChebyshevPolynomial::findRoots() const {
    // Metodo numerico: sampling + bisection
    // Per polinomi Chebyshev, campiona densamente e usa bisection
    // dove c'è cambio di segno
    
    int n = m_coeffs.size();
    if (n == 0) return {};
    if (n == 1) return {};  // Costante
    
    // Numero di campioni basato sull'ordine del polinomio
    int numSamples = std::max(50, 10 * n);
    double step = (m_t1 - m_t0) / numSamples;
    
    std::vector<double> roots;
    
    double prevT = m_t0;
    double prevVal = evaluate(m_t0);
    
    for (int i = 1; i <= numSamples; i++) {
        double t = m_t0 + i * step;
        double val = evaluate(t);
        
        // Cambio di segno → c'è una radice
        if (prevVal * val < 0) {
            // Bisection per trovare la radice
            double lo = prevT, hi = t;
            double loVal = prevVal, hiVal = val;
            
            for (int iter = 0; iter < 50; iter++) {
                double mid = (lo + hi) / 2.0;
                double midVal = evaluate(mid);
                
                if (std::abs(midVal) < 1e-14 || (hi - lo) < 1e-12) {
                    roots.push_back(mid);
                    break;
                }
                
                if (loVal * midVal < 0) {
                    hi = mid;
                    hiVal = midVal;
                } else {
                    lo = mid;
                    loVal = midVal;
                }
            }
        }
        
        prevT = t;
        prevVal = val;
    }
    
    // Rimuovi duplicati
    std::sort(roots.begin(), roots.end());
    auto last = std::unique(roots.begin(), roots.end(), 
        [](double a, double b) { return std::abs(a - b) < 1e-10; });
    roots.erase(last, roots.end());
    
    return roots;
}

std::vector<double> ChebyshevPolynomial::findMinima() const {
    // I minimi sono dove la derivata è zero e la derivata seconda è positiva
    ChebyshevPolynomial dp = derivative();
    std::vector<double> criticalPoints = dp.findRoots();
    
    ChebyshevPolynomial d2p = dp.derivative();
    
    std::vector<double> minima;
    for (double t : criticalPoints) {
        if (d2p.evaluate(t) > 0) {
            minima.push_back(t);
        }
    }
    
    // Aggiungi anche i bordi come potenziali minimi
    double v0 = evaluate(m_t0);
    double v1 = evaluate(m_t1);
    
    // Controlla se i bordi sono minimi locali
    // (sarebbero minimi se la funzione sta crescendo/decrescendo appropriatamente)
    
    return minima;
}

ChebyshevPolynomial ChebyshevPolynomial::subtract(double value) const {
    std::vector<double> newCoeffs = m_coeffs;
    if (!newCoeffs.empty()) {
        newCoeffs[0] -= value;  // T₀(x) = 1, quindi sottraiamo da c₀
    }
    return ChebyshevPolynomial(m_t0, m_t1, newCoeffs);
}

ChebyshevPolynomial ChebyshevPolynomial::add(const ChebyshevPolynomial& other) const {
    if (std::abs(m_t0 - other.m_t0) > 1e-10 || std::abs(m_t1 - other.m_t1) > 1e-10) {
        throw std::runtime_error("Cannot add polynomials with different intervals");
    }
    
    size_t maxSize = std::max(m_coeffs.size(), other.m_coeffs.size());
    std::vector<double> result(maxSize, 0.0);
    
    for (size_t i = 0; i < m_coeffs.size(); i++) {
        result[i] += m_coeffs[i];
    }
    for (size_t i = 0; i < other.m_coeffs.size(); i++) {
        result[i] += other.m_coeffs[i];
    }
    
    return ChebyshevPolynomial(m_t0, m_t1, result);
}

ChebyshevPolynomial ChebyshevPolynomial::multiply(const ChebyshevPolynomial& other) const {
    if (std::abs(m_t0 - other.m_t0) > 1e-10 || std::abs(m_t1 - other.m_t1) > 1e-10) {
        throw std::runtime_error("Cannot multiply polynomials with different intervals");
    }
    
    // Moltiplicazione di polinomi Chebyshev:
    // Tₘ(x) × Tₙ(x) = ½[Tₘ₊ₙ(x) + T|ₘ₋ₙ|(x)]
    
    int n1 = m_coeffs.size();
    int n2 = other.m_coeffs.size();
    int n_result = n1 + n2 - 1;
    
    std::vector<double> result(n_result, 0.0);
    
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            double prod = m_coeffs[i] * other.m_coeffs[j];
            int sum_idx = i + j;
            int diff_idx = std::abs(i - j);
            
            // Tᵢ × Tⱼ = ½(Tᵢ₊ⱼ + T|ᵢ₋ⱼ|)
            if (sum_idx < n_result) {
                result[sum_idx] += 0.5 * prod;
            }
            result[diff_idx] += 0.5 * prod;
        }
    }
    
    return ChebyshevPolynomial(m_t0, m_t1, result);
}

// ============================================================================
// ChebyshevRADecSegment Implementation
// ============================================================================

ChebyshevRADecSegment::ChebyshevRADecSegment(double t0, double t1, int order)
    : m_t0(t0), m_t1(t1), m_order(order) {
}

void ChebyshevRADecSegment::computeCoefficients(Ephemeris& ephemeris) {
    // Genera nodi Chebyshev
    std::vector<double> nodes;
    for (int i = 0; i <= m_order; i++) {
        double theta = M_PI * (2.0 * i + 1.0) / (2.0 * (m_order + 1));
        nodes.push_back(cos(theta));
    }
    
    // Campiona RA e Dec ai nodi
    std::vector<double> raValues, decValues;
    
    double prevRA = 0.0;
    bool first = true;
    
    for (double node : nodes) {
        double jd = m_t0 + (node + 1.0) * (m_t1 - m_t0) / 2.0;
        EphemerisData eph = ephemeris.compute(JulianDate(jd));
        
        double ra = eph.geocentricPos.ra;
        double dec = eph.geocentricPos.dec;
        
        // Gestisci wrap-around RA (0 → 2π)
        if (!first) {
            while (ra - prevRA > M_PI) ra -= 2 * M_PI;
            while (prevRA - ra > M_PI) ra += 2 * M_PI;
        }
        prevRA = ra;
        first = false;
        
        raValues.push_back(ra);
        decValues.push_back(dec);
    }
    
    // Fit coefficienti Chebyshev
    auto fitCoeffs = [this, &nodes](const std::vector<double>& values) {
        std::vector<double> coeffs(m_order + 1, 0.0);
        int n = nodes.size();
        
        for (int j = 0; j <= m_order; j++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                double T_j = (j == 0) ? 1.0 : (j == 1) ? nodes[i] : 0.0;
                if (j >= 2) {
                    double T_prev2 = 1.0, T_prev1 = nodes[i];
                    for (int k = 2; k <= j; k++) {
                        T_j = 2.0 * nodes[i] * T_prev1 - T_prev2;
                        T_prev2 = T_prev1;
                        T_prev1 = T_j;
                    }
                }
                sum += values[i] * T_j;
            }
            coeffs[j] = (j == 0) ? (sum / n) : (2.0 * sum / n);
        }
        return coeffs;
    };
    
    m_raCheb = ChebyshevPolynomial(m_t0, m_t1, fitCoeffs(raValues));
    m_decCheb = ChebyshevPolynomial(m_t0, m_t1, fitCoeffs(decValues));
}

void ChebyshevRADecSegment::evaluate(double t, double& ra, double& dec) const {
    ra = m_raCheb.evaluate(t);
    dec = m_decCheb.evaluate(t);
}

// ============================================================================
// ChebyshevOccultationDetector Implementation
// ============================================================================

ChebyshevOccultationDetector::ChebyshevOccultationDetector(const Config& config)
    : m_config(config), m_startJD(0), m_endJD(0), m_initialized(false) {
}

void ChebyshevOccultationDetector::initialize(Ephemeris& ephemeris, 
                                              double startJD, double endJD) {
    m_startJD = startJD;
    m_endJD = endJD;
    m_segments.clear();
    
    if (m_config.verbose) {
        std::cout << "ChebyshevOccultationDetector: Initializing" << std::endl;
        std::cout << "  Interval: JD " << startJD << " - " << endJD 
                  << " (" << (endJD - startJD) << " days)" << std::endl;
        std::cout << "  Order: " << m_config.order << std::endl;
        std::cout << "  Segment size: " << m_config.segmentDays << " days" << std::endl;
    }
    
    // Crea segmenti RA/Dec
    for (double jd = startJD; jd < endJD; jd += m_config.segmentDays) {
        double jdEnd = std::min(jd + m_config.segmentDays, endJD);
        
        auto segment = std::make_shared<ChebyshevRADecSegment>(
            jd, jdEnd, m_config.order);
        segment->computeCoefficients(ephemeris);
        m_segments.push_back(segment);
    }
    
    m_initialized = true;
    
    if (m_config.verbose) {
        std::cout << "  Segments created: " << m_segments.size() << std::endl;
    }
}

std::shared_ptr<ChebyshevRADecSegment> 
ChebyshevOccultationDetector::findSegment(double t) const {
    // Binary search
    int left = 0;
    int right = m_segments.size() - 1;
    
    while (left <= right) {
        int mid = (left + right) / 2;
        if (m_segments[mid]->contains(t)) {
            return m_segments[mid];
        }
        if (t < m_segments[mid]->t0()) {
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }
    return nullptr;
}

double ChebyshevOccultationDetector::angularDistance(double t, 
                                                      double starRA, 
                                                      double starDec) const {
    auto segment = findSegment(t);
    if (!segment) return 1e10;
    
    double ra, dec;
    segment->evaluate(t, ra, dec);
    
    // Distanza angolare con formula del coseno sferico
    double cosDist = sin(dec) * sin(starDec) + 
                     cos(dec) * cos(starDec) * cos(ra - starRA);
    cosDist = std::max(-1.0, std::min(1.0, cosDist));
    
    return acos(cosDist);  // radianti
}

ChebyshevPolynomial ChebyshevOccultationDetector::computeDistanceSquared(
    const ChebyshevRADecSegment& segment,
    double starRA, double starDec) const {
    
    // D²(t) ≈ (RA(t) - RA_s)² × cos²(Dec_avg) + (Dec(t) - Dec_s)²
    // Approssimazione valida per piccoli angoli
    
    const auto& raCheb = segment.raPolynomial();
    const auto& decCheb = segment.decPolynomial();
    
    // Stima Dec media per fattore cos²
    double t_mid = (segment.t0() + segment.t1()) / 2.0;
    double dec_mid = decCheb.evaluate(t_mid);
    double cosDec2 = cos(dec_mid) * cos(dec_mid);
    
    // ΔRA = RA(t) - RA_s  → sottrai costante
    ChebyshevPolynomial deltaRA = raCheb.subtract(starRA);
    
    // ΔDec = Dec(t) - Dec_s
    ChebyshevPolynomial deltaDec = decCheb.subtract(starDec);
    
    // D² = ΔRA² × cos²(Dec) + ΔDec²
    ChebyshevPolynomial deltaRA2 = deltaRA.multiply(deltaRA);
    ChebyshevPolynomial deltaDec2 = deltaDec.multiply(deltaDec);
    
    // Scala ΔRA² per cos²(Dec)
    std::vector<double> scaledCoeffs = deltaRA2.coeffs();
    for (double& c : scaledCoeffs) {
        c *= cosDec2;
    }
    ChebyshevPolynomial scaledDeltaRA2(segment.t0(), segment.t1(), scaledCoeffs);
    
    return scaledDeltaRA2.add(deltaDec2);
}

double ChebyshevOccultationDetector::refineMinimum(double t_approx, 
                                                    double starRA, 
                                                    double starDec) const {
    // Newton-Raphson per trovare minimo preciso
    // Minimizza D(t) trovando zero di dD/dt
    
    const double h = 1e-6;  // 0.1 secondi
    const int maxIter = 10;
    const double tol = 1e-8;  // ~1 ms
    
    double t = t_approx;
    
    for (int iter = 0; iter < maxIter; iter++) {
        double D = angularDistance(t, starRA, starDec);
        double D_plus = angularDistance(t + h, starRA, starDec);
        double D_minus = angularDistance(t - h, starRA, starDec);
        
        double dD = (D_plus - D_minus) / (2 * h);
        double d2D = (D_plus - 2*D + D_minus) / (h * h);
        
        if (std::abs(d2D) < 1e-20) break;
        
        double dt = -dD / d2D;
        
        // Limita step
        dt = std::max(-0.01, std::min(0.01, dt));  // max 15 minuti
        
        t += dt;
        
        if (std::abs(dt) < tol) break;
    }
    
    // Assicurati che t sia nel range
    t = std::max(m_startJD, std::min(m_endJD, t));
    
    return t;
}

std::vector<OccultationCandidate> 
ChebyshevOccultationDetector::findCandidatesForStar(double starRA, 
                                                     double starDec, 
                                                     int starIndex) const {
    std::vector<OccultationCandidate> candidates;
    
    double thresholdRad = m_config.thresholdArcsec / 206265.0;
    double thresholdRad2 = thresholdRad * thresholdRad;
    
    // Per ogni segmento
    for (const auto& segment : m_segments) {
        // Calcola polinomio D²(t)
        ChebyshevPolynomial distSq = computeDistanceSquared(*segment, starRA, starDec);
        
        // Trova minimi
        std::vector<double> minima = distSq.findMinima();
        
        // Aggiungi anche i bordi come potenziali candidati
        minima.push_back(segment->t0());
        minima.push_back(segment->t1());
        
        for (double t : minima) {
            if (t < segment->t0() || t > segment->t1()) continue;
            
            // Refina il minimo
            double t_refined = refineMinimum(t, starRA, starDec);
            
            // Calcola distanza precisa
            double dist = angularDistance(t_refined, starRA, starDec);
            double distArcsec = dist * 206265.0;
            
            if (dist < thresholdRad) {
                OccultationCandidate cand;
                cand.jd = t_refined;
                cand.minDistArcsec = distArcsec;
                cand.starRA = starRA;
                cand.starDec = starDec;
                cand.starIndex = starIndex;
                candidates.push_back(cand);
            }
        }
    }
    
    // Rimuovi duplicati (eventi troppo vicini nel tempo)
    std::sort(candidates.begin(), candidates.end(), 
              [](const auto& a, const auto& b) { return a.jd < b.jd; });
    
    std::vector<OccultationCandidate> unique;
    for (const auto& c : candidates) {
        if (unique.empty() || (c.jd - unique.back().jd) > 0.001) {  // > 1.4 minuti
            unique.push_back(c);
        } else if (c.minDistArcsec < unique.back().minDistArcsec) {
            unique.back() = c;  // Mantieni il più vicino
        }
    }
    
    return unique;
}

std::vector<OccultationCandidate> 
ChebyshevOccultationDetector::findCandidates(
    const std::vector<std::pair<double, double>>& stars) const {
    
    if (!m_initialized) {
        throw std::runtime_error("Detector not initialized");
    }
    
    std::vector<OccultationCandidate> allCandidates;
    
    for (size_t i = 0; i < stars.size(); i++) {
        auto starCandidates = findCandidatesForStar(
            stars[i].first, stars[i].second, static_cast<int>(i));
        
        allCandidates.insert(allCandidates.end(), 
                            starCandidates.begin(), starCandidates.end());
    }
    
    // Ordina per tempo
    std::sort(allCandidates.begin(), allCandidates.end(),
              [](const auto& a, const auto& b) { return a.jd < b.jd; });
    
    return allCandidates;
}

} // namespace ioccultcalc
