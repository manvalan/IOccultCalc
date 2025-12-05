#include "ioccultcalc/chebyshev_detector.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace ioccultcalc {

ChebyshevOccultationDetector::ChebyshevOccultationDetector(const Config& cfg)
    : config_(cfg), useAstDys_(false) {}

bool ChebyshevOccultationDetector::initialize(Ephemeris& ephemeris, double startJD, double endJD) {
    ChebyshevConfig chebConfig;
    chebConfig.order = config_.order;
    chebConfig.segmentDays = config_.segmentDays;
    chebConfig.geocentric = true;
    
    approximation_ = std::make_unique<ChebyshevApproximation>(chebConfig);
    bool success = approximation_->generate(ephemeris, startJD, endJD);
    
    if (success) {
        startJD_ = startJD;
        endJD_ = endJD;
        useAstDys_ = false;  // Usa approximation_ generata internamente
    }
    return success;
}

bool ChebyshevOccultationDetector::initializeFromAstDyn(const std::string& designation, 
                                                        double startJD, 
                                                        double endJD) {
    if (config_.verbose) {
        std::cout << "[Chebyshev AstDyn] Initializing from AstDyn for " << designation 
                  << " [JD " << startJD << " - " << endJD << "]\n";
    }
    
    // TODO: Implementare con nuove API AstDyn ChebyshevRKF78Propagator
    // Workflow:
    // 1. Carica file .eq1 per designation (es: "17030.eq1")
    // 2. Crea ChebyshevRKF78Propagator propagator(eq1_file)
    // 3. auto positions = propagator.propagateForChebyshev(startJD, endJD, num_points)
    // 4. Fit Chebyshev: approx.fit(positions, startJD, endJD)
    // 5. Salva coefficienti in astdysSegments_
    
    if (config_.verbose) {
        std::cerr << "[Chebyshev AstDyn] ✗ Not yet implemented - need .eq1 file integration\n";
    }
    return false;
}

bool ChebyshevOccultationDetector::initializeFromAstDynMulti(const std::string& designation,
                                                             double startJD,
                                                             double endJD,
                                                             double segmentDays) {
    if (config_.verbose) {
        std::cout << "[Chebyshev AstDyn Multi] Initializing from AstDyn for " << designation 
                  << " [JD " << startJD << " - " << endJD << "] with " 
                  << segmentDays << "-day segments\n";
    }
    
    // TODO: Implementare con nuove API AstDyn ChebyshevRKF78Propagator
    // Workflow multi-segmento:
    // 1. Divide intervallo [startJD, endJD] in segmenti di segmentDays giorni
    // 2. Per ogni segmento:
    //    - Propaga con RKF78: propagator.propagateForChebyshev(seg_start, seg_end, num_points)
    //    - Fit Chebyshev locale
    //    - Salva in astdysSegments_[i]
    // 3. Set useAstDys_ = true
    
    if (config_.verbose) {
        std::cerr << "[Chebyshev AstDyn Multi] ✗ Not yet implemented - need .eq1 file integration\n";
    }
    return false;
}

const ChebyshevEphemeris* ChebyshevOccultationDetector::findAstDysSegment(double jd) const {
    for (const auto& segment : astdysSegments_) {
        if (segment.isValid(jd)) {
            return &segment;
        }
    }
    return nullptr;
}

std::vector<OccultationCandidate> ChebyshevOccultationDetector::findCandidates(
    const std::vector<std::pair<double, double>>& stars) const {
    
    // Verifica inizializzazione (classica O AstDyn)
    if (!useAstDys_ && !approximation_) return {};
    if (useAstDys_ && astdysSegments_.empty()) return {};
    
    if (config_.verbose) {
        const char* method = useAstDys_ ? "AstDyn" : "Classic";
        std::cout << "[Chebyshev " << method << "] Searching " << stars.size() 
                  << " stars with threshold=" << config_.thresholdArcsec << " arcsec\n";
    }
    
    std::vector<OccultationCandidate> candidates;
    candidates.reserve(stars.size() / 100);
    
    for (size_t starIdx = 0; starIdx < stars.size(); starIdx++) {
        double ra = stars[starIdx].first;
        double dec = stars[starIdx].second;
        
        OccultationCandidate candidate = findCandidate(starIdx, ra, dec);
        
        if (config_.verbose && starIdx < 5) {
            std::cout << "[Chebyshev] Star " << starIdx << ": RA=" << ra << " Dec=" << dec 
                      << " minDist=" << candidate.minDistArcsec << " arcsec JD=" << candidate.jd << "\n";
        }
        
        // DEBUG: Se minDist è enorme (> 1 milione arcsec = 277 gradi), qualcosa è rotto
        if (starIdx < 3 && candidate.minDistArcsec > 100000.0) {
            std::cerr << "[Chebyshev WARNING] Star " << starIdx << " has HUGE distance: " 
                      << candidate.minDistArcsec << " arcsec - approximation may be broken!\n";
        }
        
        if (candidate.minDistArcsec < config_.thresholdArcsec) {
            if (config_.verbose) {
                std::cout << "[Chebyshev] CANDIDATE FOUND! Star " << starIdx 
                          << " minDist=" << candidate.minDistArcsec << " arcsec at JD=" << candidate.jd << "\n";
            }
            candidates.push_back(candidate);
        }
    }
    
    if (config_.verbose) {
        std::cout << "[Chebyshev] Found " << candidates.size() << " candidates\n";
    }
    
    return candidates;
}

OccultationCandidate ChebyshevOccultationDetector::findCandidate(int starIndex, double ra, double dec) const {
    OccultationCandidate result;
    result.starIndex = starIndex;
    result.jd = 0.0;
    result.minDistArcsec = 1e10;
    result.uncertaintyArcsec = 0.0;
    
    // DEBUG: Stella target specifica (Gaia 3411546266140512128)
    bool isTargetStar = (std::abs(ra - 4.89440726) < 0.001 && std::abs(dec - 20.3316615) < 0.001);
    
    if (isTargetStar && config_.verbose) {
        std::cout << "\n[CHEBYSHEV DEBUG] ★★★ STELLA TARGET TROVATA ★★★\n";
        std::cout << "  Star Index: " << starIndex << "\n";
        std::cout << "  RA:  " << std::fixed << std::setprecision(8) << ra << "° (input)\n";
        std::cout << "  Dec: " << dec << "° (input)\n";
        std::cout << "  Intervallo JD: " << std::setprecision(4) << startJD_ << " - " << endJD_ << "\n";
        std::cout << "  Mode: " << (useAstDys_ ? "AstDyn" : "Classic") << "\n";
    }
    
    // Verifica che sia inizializzato (classico O AstDyn)
    if (!useAstDys_ && !approximation_) {
        if (config_.verbose) {
            std::cout << "[Chebyshev] ERROR: not initialized (no approximation)\n";
        }
        return result;
    }
    
    if (useAstDys_ && astdysSegments_.empty()) {
        if (config_.verbose) {
            std::cout << "[Chebyshev] ERROR: not initialized (no AstDyn segments)\n";
        }
        return result;
    }
    
    double jdStart = startJD_;
    double jdEnd = endJD_;
    int nSamples = 100;
    
    double minDist = 1e10;
    double minJD = (jdStart + jdEnd) / 2.0;
    
    // DEBUG: Campionamento per stella target
    int debugSample = 0;
    
    for (int i = 0; i < nSamples; i++) {
        double jd = jdStart + (jdEnd - jdStart) * i / (nSamples - 1);
        
        double ast_ra, ast_dec, ast_dist;
        
        // USA ASTDYN se disponibile, altrimenti approximation_ classica
        if (useAstDys_) {
            // Trova segmento valido per questa JD
            const ChebyshevEphemeris* segment = findAstDysSegment(jd);
            if (!segment) continue;  // JD fuori da tutti i segmenti
            
            // Valuta RA/Dec da Chebyshev AstDyn
            auto radec = segment->getRADecDist(jd);
            ast_ra = radec[0];
            ast_dec = radec[1];
            ast_dist = radec[2];
        } else {
            // Valuta con approximation_ classica
            if (!approximation_->evaluate(jd, ast_ra, ast_dec, ast_dist)) continue;
        }
        
        double dist = angularDistance(ast_ra, ast_dec, ra, dec);
        
        // DEBUG: Stampa primi 10 campioni per stella target
        if (isTargetStar && config_.verbose && debugSample < 10) {
            std::cout << "  [Sample " << debugSample << "] JD=" << std::setprecision(4) << jd 
                      << " Ast(RA=" << std::setprecision(6) << ast_ra << "°, Dec=" << ast_dec 
                      << "°) Dist=" << std::setprecision(2) << dist << "\"\n";
            debugSample++;
        }
        
        // DEBUG: Cerca intorno a JD evento atteso (2461007.9306)
        if (isTargetStar && config_.verbose && std::abs(jd - 2461007.9306) < 0.1) {
            std::cout << "  [NEAR EVENT] JD=" << std::setprecision(4) << jd 
                      << " Ast(RA=" << std::setprecision(6) << ast_ra << "°, Dec=" << ast_dec 
                      << "°) Dist=" << std::setprecision(2) << dist << "\"\n";
        }
        
        if (dist < minDist) {
            minDist = dist;
            minJD = jd;
        }
    }
    
    if (isTargetStar && config_.verbose) {
        std::cout << "  [RESULT] Min Distance: " << std::setprecision(2) << minDist 
                  << "\" @ JD " << std::setprecision(4) << minJD << "\n";
        std::cout << "  Threshold: " << config_.thresholdArcsec << "\"\n";
        std::cout << "  Passed: " << (minDist < config_.thresholdArcsec ? "YES ✓" : "NO ✗") << "\n\n";
    }
    
    if (minDist < config_.refinementArcsec) {
        jdStart = minJD - 0.1;
        jdEnd = minJD + 0.1;
        nSamples = 50;
        
        for (int i = 0; i < nSamples; i++) {
            double jd = jdStart + (jdEnd - jdStart) * i / (nSamples - 1);
            
            double ast_ra, ast_dec, ast_dist;
            if (!approximation_->evaluate(jd, ast_ra, ast_dec, ast_dist)) continue;
            
            double dist = angularDistance(ast_ra, ast_dec, ra, dec);
            
            if (dist < minDist) {
                minDist = dist;
                minJD = jd;
            }
        }
    }
    
    result.jd = minJD;
    result.minDistArcsec = minDist;
    result.uncertaintyArcsec = 10.0;
    return result;
}

double ChebyshevOccultationDetector::angularDistance(double ra1, double dec1, double ra2, double dec2) const {
    double ra1_rad = ra1 * M_PI / 180.0;
    double dec1_rad = dec1 * M_PI / 180.0;
    double ra2_rad = ra2 * M_PI / 180.0;
    double dec2_rad = dec2 * M_PI / 180.0;
    
    double dra = ra2_rad - ra1_rad;
    double ddec = dec2_rad - dec1_rad;
    
    double a = sin(ddec/2.0) * sin(ddec/2.0) +
               cos(dec1_rad) * cos(dec2_rad) *
               sin(dra/2.0) * sin(dra/2.0);
    
    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
    return c * 180.0 / M_PI * 3600.0;
}

} // namespace ioccultcalc
