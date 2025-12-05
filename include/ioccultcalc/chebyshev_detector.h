#ifndef IOCCULTCALC_CHEBYSHEV_DETECTOR_H
#define IOCCULTCALC_CHEBYSHEV_DETECTOR_H

#include "chebyshev_approximation.h"
#include "astdys_client.h"
#include <vector>
#include <memory>

namespace ioccultcalc {

class Ephemeris;

struct OccultationCandidate {
    int starIndex;
    double jd;
    double minDistArcsec;
    double uncertaintyArcsec;
};

class ChebyshevOccultationDetector {
public:
    struct Config {
        int order = 11;
        double segmentDays = 1.0;
        double thresholdArcsec = 300.0;
        double refinementArcsec = 60.0;
        bool verbose = false;
    };
    
    explicit ChebyshevOccultationDetector(const Config& cfg);
    
    // Inizializzazione CLASSICA: genera Chebyshev da Ephemeris interno
    bool initialize(Ephemeris& ephemeris, double startJD, double endJD);
    
    // Inizializzazione ASTDYN: usa coefficienti Chebyshev da libreria AstDyn
    // Scarica dati da AstDyS (sito web) tramite AstDysClient
    // VANTAGGI: precisione sub-arcsec, velocità 1000x, coerenza con AstDyS
    bool initializeFromAstDyn(const std::string& designation, double startJD, double endJD);
    
    // Inizializzazione ASTDYN AVANZATA: multi-segmento per intervalli lunghi
    bool initializeFromAstDynMulti(const std::string& designation, 
                                   double startJD, 
                                   double endJD,
                                   double segmentDays = 7.0);
    
    std::vector<OccultationCandidate> findCandidates(const std::vector<std::pair<double, double>>& stars) const;
    
    OccultationCandidate findCandidate(int starIndex, double ra, double dec) const;

private:
    Config config_;
    std::unique_ptr<ChebyshevApproximation> approximation_;
    std::vector<ChebyshevEphemeris> astdysSegments_;  // Segmenti da AstDys
    bool useAstDys_;  // Flag: true = usa astdysSegments_, false = usa approximation_
    double startJD_;
    double endJD_;
    
    double angularDistance(double ra1, double dec1, double ra2, double dec2) const;
    
    // Helper: trova segmento AstDys valido per data JD
    const ChebyshevEphemeris* findAstDysSegment(double jd) const;
};

} // namespace ioccultcalc

#endif
