#include "ioccultcalc/chebyshev_approximation.h"
#include "ioccultcalc/coordinates.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace ioccultcalc {

// ============================================================================
// ChebyshevSegment Implementation
// ============================================================================

ChebyshevSegment::ChebyshevSegment(const JulianDate& t0, const JulianDate& t1, int order)
    : m_t0(t0), m_t1(t1), m_order(order), m_estimatedError(0.0) {
    
    if (order < 5 || order > 25) {
        throw std::runtime_error("Chebyshev order must be 5-25");
    }
    
    // Alloca spazio per coefficienti
    m_coeffsX.resize(order + 1);
    m_coeffsY.resize(order + 1);
    m_coeffsZ.resize(order + 1);
    m_coeffsVX.resize(order + 1);
    m_coeffsVY.resize(order + 1);
    m_coeffsVZ.resize(order + 1);
}

void ChebyshevSegment::computeCoefficients(Ephemeris& ephemeris) {
    // Genera nodi di campionamento Chebyshev ottimali
    std::vector<double> nodes = ChebyshevFitter::chebyshevNodes(m_order);
    
    // Campiona posizioni e velocità ai nodi
    std::vector<double> samplesX, samplesY, samplesZ;
    std::vector<double> samplesVX, samplesVY, samplesVZ;
    
    for (double node : nodes) {
        // Converti nodo [-1, 1] → tempo JD
        double jd = m_t0.jd + (node + 1.0) * (m_t1.jd - m_t0.jd) / 2.0;
        
        // Calcola ephemeris con integrazione esatta
        EphemerisData eph = ephemeris.compute(JulianDate(jd));
        
        // Estrai posizione geocentrica (assumendo già in km)
        Vector3D pos = Coordinates::equatorialToCartesian(eph.geocentricPos);
        pos = pos * (eph.distance * AU);  // Converti in km
        
        samplesX.push_back(pos.x);
        samplesY.push_back(pos.y);
        samplesZ.push_back(pos.z);
        
        // Velocità (km/s)
        Vector3D vel = eph.heliocentricVel * (AU / 86400.0);
        samplesVX.push_back(vel.x);
        samplesVY.push_back(vel.y);
        samplesVZ.push_back(vel.z);
    }
    
    // Fit coefficienti Chebyshev per ogni componente
    m_coeffsX = ChebyshevFitter::fit(m_order, nodes, samplesX);
    m_coeffsY = ChebyshevFitter::fit(m_order, nodes, samplesY);
    m_coeffsZ = ChebyshevFitter::fit(m_order, nodes, samplesZ);
    m_coeffsVX = ChebyshevFitter::fit(m_order, nodes, samplesVX);
    m_coeffsVY = ChebyshevFitter::fit(m_order, nodes, samplesVY);
    m_coeffsVZ = ChebyshevFitter::fit(m_order, nodes, samplesVZ);
    
    // Stima errore usando punti test intermedi
    std::vector<double> testNodes;
    std::vector<double> testX, testY, testZ;
    
    for (int i = 0; i < m_order; i++) {
        double nodeTest = -1.0 + (2.0 * i + 1.5) / (m_order + 1);
        testNodes.push_back(nodeTest);
        
        double jd = m_t0.jd + (nodeTest + 1.0) * (m_t1.jd - m_t0.jd) / 2.0;
        EphemerisData eph = ephemeris.compute(JulianDate(jd));
        Vector3D pos = Coordinates::equatorialToCartesian(eph.geocentricPos);
        pos = pos * (eph.distance * AU);
        
        testX.push_back(pos.x);
        testY.push_back(pos.y);
        testZ.push_back(pos.z);
    }
    
    double errX = ChebyshevFitter::evaluateError(m_coeffsX, testNodes, testX);
    double errY = ChebyshevFitter::evaluateError(m_coeffsY, testNodes, testY);
    double errZ = ChebyshevFitter::evaluateError(m_coeffsZ, testNodes, testZ);
    
    m_estimatedError = sqrt(errX*errX + errY*errY + errZ*errZ);
}

EphemerisData ChebyshevSegment::evaluate(const JulianDate& t) const {
    if (!contains(t)) {
        throw std::runtime_error("Time outside segment range");
    }
    
    // Normalizza tempo
    double x = normalizeTime(t);
    
    // Valuta polinomi per posizione
    double posX = evaluatePolynomial(m_coeffsX, x);
    double posY = evaluatePolynomial(m_coeffsY, x);
    double posZ = evaluatePolynomial(m_coeffsZ, x);
    
    // Valuta polinomi per velocità
    double velX = evaluatePolynomial(m_coeffsVX, x);
    double velY = evaluatePolynomial(m_coeffsVY, x);
    double velZ = evaluatePolynomial(m_coeffsVZ, x);
    
    // Converti back in coordinate equatoriali
    Vector3D posKm(posX, posY, posZ);
    double distance = posKm.magnitude();
    Vector3D posDir = posKm.normalize();
    
    EphemerisData eph;
    eph.jd = t;
    eph.geocentricPos = Coordinates::cartesianToEquatorial(posDir);
    eph.distance = distance / AU;  // km → AU
    
    Vector3D velKmS(velX, velY, velZ);
    eph.heliocentricVel = velKmS * (86400.0 / AU);  // km/s → AU/day
    
    return eph;
}

bool ChebyshevSegment::contains(const JulianDate& t) const {
    return t.jd >= m_t0.jd && t.jd <= m_t1.jd;
}

double ChebyshevSegment::normalizeTime(const JulianDate& t) const {
    // Mappa [t0, t1] → [-1, 1]
    return 2.0 * (t.jd - m_t0.jd) / (m_t1.jd - m_t0.jd) - 1.0;
}

double ChebyshevSegment::evaluatePolynomial(const std::vector<double>& coeffs, double x) const {
    // Valuta Σ cᵢ × Tᵢ(x) usando ricorsione Clenshaw
    // Più stabile numericamente della formula diretta
    
    int n = coeffs.size() - 1;
    if (n < 0) return 0.0;
    
    double b_k2 = 0.0;
    double b_k1 = 0.0;
    
    for (int k = n; k >= 1; k--) {
        double b_k = 2.0 * x * b_k1 - b_k2 + coeffs[k];
        b_k2 = b_k1;
        b_k1 = b_k;
    }
    
    return x * b_k1 - b_k2 + coeffs[0];
}

void ChebyshevSegment::evaluateChebyshevBasis(double x, std::vector<double>& T) const {
    T.resize(m_order + 1);
    
    if (m_order >= 0) T[0] = 1.0;
    if (m_order >= 1) T[1] = x;
    
    // Ricorsione: Tₙ₊₁(x) = 2x×Tₙ(x) - Tₙ₋₁(x)
    for (int n = 2; n <= m_order; n++) {
        T[n] = 2.0 * x * T[n-1] - T[n-2];
    }
}

// ============================================================================
// ChebyshevCache Implementation
// ============================================================================

ChebyshevCache::ChebyshevCache(const ChebyshevConfig& config)
    : m_config(config), m_initialized(false), m_maxError(0.0), m_avgError(0.0) {
}

void ChebyshevCache::initialize(Ephemeris& ephemeris, 
                               const JulianDate& startJD,
                               const JulianDate& endJD) {
    
    m_startJD = startJD;
    m_endJD = endJD;
    m_segments.clear();
    
    std::cout << "Initializing Chebyshev cache:" << std::endl;
    std::cout << "  Interval: JD " << startJD.jd << " - " << endJD.jd 
              << " (" << (endJD.jd - startJD.jd) << " days)" << std::endl;
    std::cout << "  Order: " << m_config.order << std::endl;
    std::cout << "  Segment size: " << m_config.intervalDays << " days" << std::endl;
    
    // Crea segmenti
    int segmentCount = 0;
    double totalError = 0.0;
    
    for (double jd = startJD.jd; jd < endJD.jd; jd += m_config.intervalDays) {
        double jdEnd = std::min(jd + m_config.intervalDays, endJD.jd);
        
        auto segment = std::make_shared<ChebyshevSegment>(
            JulianDate(jd), JulianDate(jdEnd), m_config.order);
        
        segment->computeCoefficients(ephemeris);
        
        double error = segment->estimatedError();
        totalError += error;
        m_maxError = std::max(m_maxError, error);
        
        m_segments.push_back(segment);
        segmentCount++;
        
        if (segmentCount % 10 == 0) {
            std::cout << "  Computed " << segmentCount << " segments..." << std::endl;
        }
    }
    
    m_avgError = totalError / segmentCount;
    m_initialized = true;
    
    std::cout << "Chebyshev cache initialized:" << std::endl;
    std::cout << "  Segments: " << segmentCount << std::endl;
    std::cout << "  Max error: " << m_maxError << " km" << std::endl;
    std::cout << "  Avg error: " << m_avgError << " km" << std::endl;
    
    if (m_maxError > m_config.maxErrorKm) {
        std::cerr << "WARNING: Max Chebyshev error " << m_maxError 
                  << " km exceeds threshold " << m_config.maxErrorKm << " km" << std::endl;
    }
}

EphemerisData ChebyshevCache::evaluate(const JulianDate& t) const {
    if (!m_initialized) {
        throw std::runtime_error("Chebyshev cache not initialized");
    }
    
    auto segment = findSegment(t);
    if (!segment) {
        throw std::runtime_error("Time outside cache range");
    }
    
    return segment->evaluate(t);
}

std::vector<JulianDate> ChebyshevCache::generateFineTimesteps() const {
    std::vector<JulianDate> timesteps;
    
    double stepDays = m_config.fineTimestepMinutes / 1440.0;
    
    for (double jd = m_startJD.jd; jd <= m_endJD.jd; jd += stepDays) {
        timesteps.push_back(JulianDate(jd));
    }
    
    return timesteps;
}

std::shared_ptr<ChebyshevSegment> ChebyshevCache::findSegment(const JulianDate& t) const {
    // Binary search per O(log n)
    int left = 0;
    int right = m_segments.size() - 1;
    
    while (left <= right) {
        int mid = (left + right) / 2;
        
        if (m_segments[mid]->contains(t)) {
            return m_segments[mid];
        }
        
        if (t.jd < m_segments[mid]->startJD()) {
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }
    
    return nullptr;
}

// ============================================================================
// ChebyshevFitter Implementation
// ============================================================================

std::vector<double> ChebyshevFitter::fit(int order, 
                                        const std::vector<double>& t,
                                        const std::vector<double>& y) {
    
    if (t.size() != y.size() || t.size() < order + 1) {
        throw std::runtime_error("Insufficient data points for Chebyshev fit");
    }
    
    std::vector<double> coeffs(order + 1, 0.0);
    int n = t.size();
    
    // Discrete Chebyshev Transform (DCT)
    // cⱼ = (2/n) × Σᵢ yᵢ × Tⱼ(xᵢ) con correzione per j=0
    
    for (int j = 0; j <= order; j++) {
        double sum = 0.0;
        
        for (int i = 0; i < n; i++) {
            double T_j = 0.0;
            
            // Calcola Tⱼ(xᵢ) usando ricorsione
            if (j == 0) {
                T_j = 1.0;
            } else if (j == 1) {
                T_j = t[i];
            } else {
                double T_prev2 = 1.0;
                double T_prev1 = t[i];
                for (int k = 2; k <= j; k++) {
                    T_j = 2.0 * t[i] * T_prev1 - T_prev2;
                    T_prev2 = T_prev1;
                    T_prev1 = T_j;
                }
            }
            
            sum += y[i] * T_j;
        }
        
        coeffs[j] = (j == 0) ? (sum / n) : (2.0 * sum / n);
    }
    
    return coeffs;
}

std::vector<double> ChebyshevFitter::chebyshevNodes(int order) {
    std::vector<double> nodes;
    
    // Nodi di Chebyshev: xᵢ = cos(π(2i+1)/(2(order+1)))
    // Ottimali per minimizzare oscillazione di Runge
    
    for (int i = 0; i <= order; i++) {
        double theta = M_PI * (2.0 * i + 1.0) / (2.0 * (order + 1));
        double node = cos(theta);
        nodes.push_back(node);
    }
    
    return nodes;
}

double ChebyshevFitter::evaluateError(const std::vector<double>& coeffs,
                                     const std::vector<double>& t_test,
                                     const std::vector<double>& y_test) {
    
    if (t_test.size() != y_test.size()) {
        throw std::runtime_error("Test data size mismatch");
    }
    
    double sumSquaredError = 0.0;
    
    for (size_t i = 0; i < t_test.size(); i++) {
        // Valuta polinomio ai punti test
        double y_approx = 0.0;
        
        for (size_t j = 0; j < coeffs.size(); j++) {
            double T_j = 0.0;
            
            if (j == 0) {
                T_j = 1.0;
            } else if (j == 1) {
                T_j = t_test[i];
            } else {
                double T_prev2 = 1.0;
                double T_prev1 = t_test[i];
                for (size_t k = 2; k <= j; k++) {
                    T_j = 2.0 * t_test[i] * T_prev1 - T_prev2;
                    T_prev2 = T_prev1;
                    T_prev1 = T_j;
                }
            }
            
            y_approx += coeffs[j] * T_j;
        }
        
        double error = y_approx - y_test[i];
        sumSquaredError += error * error;
    }
    
    return sqrt(sumSquaredError / t_test.size());  // RMS error
}

} // namespace ioccultcalc
