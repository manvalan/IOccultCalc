#include "ioccultcalc/chebyshev_approximation.h"
#include "ioccultcalc/ephemeris.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace ioccultcalc {

ChebyshevApproximation::ChebyshevApproximation(const ChebyshevConfig& cfg)
    : config_(cfg) {}

bool ChebyshevApproximation::generate(Ephemeris& ephemeris, double startJD, double endJD) {
    segments_.clear();
    if (endJD <= startJD) return false;
    
    double currentJD = startJD;
    while (currentJD < endJD) {
        double segmentEnd = std::min(currentJD + config_.segmentDays, endJD);
        ChebyshevSegment segment;
        if (!generateSegment(ephemeris, currentJD, segmentEnd, segment)) return false;
        segments_.push_back(segment);
        currentJD = segmentEnd;
    }
    
    // DEBUG: Verifica ordine segmenti
    std::cout << "[SEGMENTS] Generated " << segments_.size() << " segments:\n";
    for (size_t i = 0; i < std::min(segments_.size(), size_t(3)); i++) {
        std::cout << "  Segment " << i << ": [" << std::fixed << std::setprecision(2) 
                  << segments_[i].startJD << ", " << segments_[i].endJD << "] mid=" << segments_[i].midJD << "\n";
    }
    
    return !segments_.empty();
}

bool ChebyshevApproximation::generateSegment(Ephemeris& ephemeris, double startJD, double endJD, ChebyshevSegment& segment) {
    segment.startJD = startJD;
    segment.endJD = endJD;
    segment.midJD = (startJD + endJD) / 2.0;
    segment.halfSpan = (endJD - startJD) / 2.0;
    
    int n = config_.order + 1;
    std::vector<double> xVals(n), yVals(n), zVals(n);
    
    std::cout << "[CHEB FIT] Segment JD " << startJD << " to " << endJD 
              << " (mid=" << segment.midJD << ", order=" << config_.order << ")\n";
    
    for (int i = 0; i < n; i++) {
        double theta = M_PI * (i + 0.5) / n;
        double x = std::cos(theta);  // FIXED: rimosso il segno negativo per allineare con DCT
        double jd = segment.midJD + x * segment.halfSpan;
        
        JulianDate epoch;
        epoch.jd = jd;
        EphemerisData eph = ephemeris.compute(epoch);
        
        double px, py, pz;
        if (config_.geocentric) {
            // NOTA: eph.geocentricPos.ra/dec sono GIÀ in radianti (vedi types.h)
            double ra_rad = eph.geocentricPos.ra;
            double dec_rad = eph.geocentricPos.dec;
            double r = eph.distance;
            px = r * std::cos(dec_rad) * std::cos(ra_rad);
            py = r * std::cos(dec_rad) * std::sin(ra_rad);
            pz = r * std::sin(dec_rad);
            
            // DEBUG: Mostra primi e ultimi 3 punti del fit
            if (i < 3 || i >= n-3) {
                double ra_deg = ra_rad * 180.0 / M_PI;
                double dec_deg = dec_rad * 180.0 / M_PI;
                std::cout << "  [PROPAGATOR " << i << "] JD=" << std::fixed << std::setprecision(2) << jd 
                          << " (x=" << std::setprecision(3) << x << ", t=" << ((jd - segment.midJD) / segment.halfSpan) << ")"
                          << " RA=" << std::setprecision(4) << ra_deg << "° Dec=" << dec_deg 
                          << "° xyz=(" << std::setprecision(6) << px << ", " << py << ", " << pz << ") AU\n";
            }
        } else {
            throw std::runtime_error("Solo geocentrico supportato");
        }
        xVals[i] = px;
        yVals[i] = py;
        zVals[i] = pz;
    }
    
    segment.coeffX = computeChebyshevCoefficients(xVals);
    segment.coeffY = computeChebyshevCoefficients(yVals);
    segment.coeffZ = computeChebyshevCoefficients(zVals);
    
    // VERIFICA: Test al centro e agli estremi
    static int testCount = 0;
    if (testCount < 1) {
        std::cout << "  [SAMPLES] xVals[0]=" << std::setprecision(6) << xVals[0] 
                  << " xVals[n-1]=" << xVals[n-1] << "\n";
        std::cout << "  [COEFFS] coeffX[0]=" << segment.coeffX[0] 
                  << " coeffX[1]=" << segment.coeffX[1] << "\n";
        // Test t=-1 (inizio), t=0 (centro), t=+1 (fine)
        for (double t : {-1.0, 0.0, 1.0}) {
            double x_eval = evaluateChebyshevPolynomial(segment.coeffX, t);
            double y_eval = evaluateChebyshevPolynomial(segment.coeffY, t);
            double z_eval = evaluateChebyshevPolynomial(segment.coeffZ, t);
            double dist = std::sqrt(x_eval*x_eval + y_eval*y_eval + z_eval*z_eval);
            double ra_eval = std::atan2(y_eval, x_eval) * 180.0 / M_PI;
            if (ra_eval < 0) ra_eval += 360.0;
            double dec_eval = std::asin(z_eval / dist) * 180.0 / M_PI;
            double jd_test = segment.midJD + t * segment.halfSpan;
            
            std::cout << "  [VERIFY t=" << std::setw(4) << t << "] JD=" << std::fixed << std::setprecision(2) << jd_test
                      << " RA=" << std::setprecision(4) << ra_eval << "° Dec=" << dec_eval << "°\n";
        }
        testCount++;
    }
    
    return true;
}

std::vector<double> ChebyshevApproximation::computeChebyshevCoefficients(const std::vector<double>& y) const {
    int n = y.size();
    std::vector<double> coeffs(config_.order + 1, 0.0);
    
    for (int j = 0; j <= config_.order; j++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            double theta = M_PI * (i + 0.5) / n;
            sum += y[i] * std::cos(j * theta);
        }
        coeffs[j] = (j == 0) ? (sum / n) : (2.0 * sum / n);
    }
    return coeffs;
}

bool ChebyshevApproximation::evaluate(double jd, double& ra, double& dec, double& distance) const {
    const ChebyshevSegment* segment = nullptr;
    static int segSelectCount = 0;
    for (const auto& seg : segments_) {
        if (jd >= seg.startJD && jd <= seg.endJD) {
            segment = &seg;
            if (segSelectCount < 3) {
                std::cout << "    [SEG SELECT " << segSelectCount << "] JD=" << std::fixed << std::setprecision(2) << jd 
                          << " → segment [" << seg.startJD << ", " << seg.endJD << "] mid=" << seg.midJD << "\n";
                segSelectCount++;
            }
            break;
        }
    }
    if (!segment) return false;
    
    double t = normalizeTime(jd, *segment);
    double x = evaluateChebyshevPolynomial(segment->coeffX, t);
    double y = evaluateChebyshevPolynomial(segment->coeffY, t);
    double z = evaluateChebyshevPolynomial(segment->coeffZ, t);
    
    distance = std::sqrt(x*x + y*y + z*z);
    if (distance == 0.0) return false;
    
    ra = std::atan2(y, x) * 180.0 / M_PI;
    if (ra < 0.0) ra += 360.0;
    dec = std::asin(z / distance) * 180.0 / M_PI;
    
    // DEBUG: Mostra valutazioni chiave (inizio, centro, fine del primo segmento)
    static int evalCount = 0;
    static bool testedMid = false, testedEnd = false;
    bool shouldPrint = (evalCount < 3) || 
                       (!testedMid && fabs(jd - segment->midJD) < 0.01) ||
                       (!testedEnd && fabs(jd - segment->endJD) < 0.01);
    
    if (shouldPrint) {
        std::cout << "  [CHEBYSHEV " << evalCount << "] JD=" << std::fixed << std::setprecision(2) << jd 
                  << " RA=" << std::setprecision(4) << ra << "° Dec=" << dec 
                  << "° xyz=(" << std::setprecision(6) << x << ", " << y << ", " << z << ") AU";
        if (fabs(jd - segment->midJD) < 0.01) {
            std::cout << " ← MID";
            testedMid = true;
        }
        if (fabs(jd - segment->endJD) < 0.01) {
            std::cout << " ← END";
            testedEnd = true;
        }
        std::cout << "\n";
        evalCount++;
    }
    
    return true;
}

bool ChebyshevApproximation::evaluateStateVector(double jd, double pos[3], double vel[3]) const {
    const ChebyshevSegment* segment = nullptr;
    for (const auto& seg : segments_) {
        if (jd >= seg.startJD && jd <= seg.endJD) {
            segment = &seg;
            break;
        }
    }
    if (!segment) return false;
    
    double t = normalizeTime(jd, *segment);
    pos[0] = evaluateChebyshevPolynomial(segment->coeffX, t);
    pos[1] = evaluateChebyshevPolynomial(segment->coeffY, t);
    pos[2] = evaluateChebyshevPolynomial(segment->coeffZ, t);
    vel[0] = evaluateChebyshevDerivative(segment->coeffX, t) / segment->halfSpan;
    vel[1] = evaluateChebyshevDerivative(segment->coeffY, t) / segment->halfSpan;
    vel[2] = evaluateChebyshevDerivative(segment->coeffZ, t) / segment->halfSpan;
    return true;
}

double ChebyshevApproximation::evaluateChebyshevPolynomial(const std::vector<double>& coeffs, double x) const {
    if (coeffs.empty()) return 0.0;
    if (coeffs.size() == 1) return coeffs[0];
    
    double b_k2 = 0.0, b_k1 = 0.0;
    for (int k = coeffs.size() - 1; k >= 0; k--) {
        double b_k = coeffs[k] + 2.0 * x * b_k1 - b_k2;
        b_k2 = b_k1;
        b_k1 = b_k;
    }
    return b_k1 - x * b_k2;
}

double ChebyshevApproximation::evaluateChebyshevDerivative(const std::vector<double>& coeffs, double x) const {
    if (coeffs.size() <= 1) return 0.0;
    std::vector<double> derivCoeffs(coeffs.size() - 1);
    for (size_t n = 1; n < coeffs.size(); n++) {
        derivCoeffs[n - 1] = n * coeffs[n];
    }
    return evaluateChebyshevPolynomial(derivCoeffs, x);
}

double ChebyshevApproximation::estimateMaxError(Ephemeris& ephemeris) const {
    double maxError = 0.0;
    for (const auto& seg : segments_) {
        for (int i = 0; i < 20; i++) {
            double fraction = (i + 0.3) / 20.0;
            double jd = seg.startJD + fraction * (seg.endJD - seg.startJD);
            
            JulianDate epoch;
            epoch.jd = jd;
            EphemerisData exact = ephemeris.compute(epoch);
            
            double ra_approx, dec_approx, dist_approx;
            if (!evaluate(jd, ra_approx, dec_approx, dist_approx)) continue;
            
            double dRA = (exact.geocentricPos.ra - ra_approx) * std::cos(exact.geocentricPos.dec * M_PI / 180.0);
            double dDec = exact.geocentricPos.dec - dec_approx;
            double error = std::sqrt(dRA*dRA + dDec*dDec) * 3600.0;
            maxError = std::max(maxError, error);
        }
    }
    return maxError;
}

double ChebyshevApproximation::normalizeTime(double jd, const ChebyshevSegment& segment) const {
    return (jd - segment.midJD) / segment.halfSpan;
}

} // namespace ioccultcalc
