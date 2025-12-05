#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/orbit_fitter.h"
#include "ioccultcalc/observation.h"
#include "ioccultcalc/time_utils.h"

using namespace ioccultcalc;

// Parse osservazioni da file .rwo
std::vector<AstrometricObservation> parseRWO(const std::string& filename) {
    std::vector<AstrometricObservation> observations;
    std::ifstream file(filename);
    
    if (!file.good()) {
        throw std::runtime_error("File " + filename + " non trovato!");
    }
    
    std::string line;
    bool inHeader = true;
    
    while (std::getline(file, line)) {
        // Skip header
        if (line.find("END_OF_HEADER") != std::string::npos) {
            inHeader = false;
            continue;
        }
        if (inHeader || line.empty() || line[0] == '!') {
            continue;
        }
        
        // Parse observation line
        // Format: designation obs_type yyyy mm dd.ddddd ra_h ra_m ra_s dec_d dec_m dec_s mag obs_code
        if (line.length() < 80) continue;
        
        try {
            std::istringstream iss(line);
            std::string designation, obsType, flag1, flag2;
            int year, month;
            double day;
            double ra_h, ra_m, ra_s, dec_d, dec_m, dec_s;
            
            // Parse basic fields
            iss >> designation >> obsType >> flag1 >> flag2 >> year >> month >> day;
            
            // Skip accuracy field
            std::string accuracy;
            iss >> accuracy;
            
            // Parse RA (HH MM SS.sss)
            iss >> ra_h >> ra_m >> ra_s;
            
            // Skip RA accuracy fields
            for (int i = 0; i < 4; i++) iss >> accuracy;
            
            // Parse sign and Dec
            std::string dec_sign;
            iss >> dec_sign;
            dec_d = std::stod(dec_sign.substr(1)); // Skip +/- sign
            iss >> dec_m >> dec_s;
            
            // Convert to decimal degrees
            double ra_deg = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0;
            double dec_deg = dec_d + dec_m/60.0 + dec_s/3600.0;
            if (dec_sign[0] == '-') dec_deg = -dec_deg;
            
            // Convert to JD
            JulianDate epoch;
            int hour = (int)((day - (int)day) * 24.0);
            int minute = (int)(((day - (int)day) * 24.0 - hour) * 60.0);
            double second = ((((day - (int)day) * 24.0 - hour) * 60.0) - minute) * 60.0;
            
            epoch = TimeUtils::calendarToJD(year, month, (int)day, hour, minute, second);
            
            // Create observation
            AstrometricObservation obs;
            obs.epoch = epoch;
            obs.obs.ra = ra_deg * DEG_TO_RAD;
            obs.obs.dec = dec_deg * DEG_TO_RAD;
            obs.raError = 0.5;  // 0.5 arcsec
            obs.decError = 0.5;
            obs.observatoryCode = ""; // Parsed later if needed
            
            observations.push_back(obs);
            
        } catch (const std::exception& e) {
            // Skip malformed lines
            continue;
        }
    }
    
    file.close();
    return observations;
}

int main() {
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "ORBIT FITTING: (17030) Sierks con osservazioni RWO\n";
    std::cout << "═══════════════════════════════════════════════════════════\n\n";
    
    // ================================================================
    // STEP 1: Carica elementi iniziali da AstDyS
    // ================================================================
    std::cout << "STEP 1: Caricamento elementi iniziali AstDyS\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    
    EquinoctialElements astElem;
    std::ifstream eq1File("17030_astdys.eq1");
    
    if (!eq1File.good()) {
        std::cerr << "❌ File 17030_astdys.eq1 non trovato!\n";
        return 1;
    }
    
    std::string line;
    bool found = false;
    while (std::getline(eq1File, line)) {
        if (line.find(" EQU ") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            double a, h, k, p, q, lambda;
            if (iss >> tag >> a >> h >> k >> p >> q >> lambda) {
                astElem.a = a;
                astElem.h = h;
                astElem.k = k;
                astElem.p = p;
                astElem.q = q;
                astElem.lambda = lambda;
                found = true;
            }
        } else if (found && line.find(" MJD ") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            double mjd;
            if (iss >> tag >> mjd) {
                astElem.epoch.jd = mjd + 2400000.5;
                break;
            }
        }
    }
    eq1File.close();
    
    if (!found) {
        std::cerr << "❌ Elementi non trovati nel file!\n";
        return 1;
    }
    
    auto kepInit = astElem.toKeplerian();
    std::cout << "Elementi iniziali (AstDyS):\n";
    std::cout << "  Epoca = JD " << std::fixed << std::setprecision(1) << kepInit.epoch.jd << "\n";
    std::cout << "  a = " << std::setprecision(8) << kepInit.a << " AU\n";
    std::cout << "  e = " << std::setprecision(6) << kepInit.e << "\n";
    std::cout << "  i = " << std::setprecision(4) << (kepInit.i * RAD_TO_DEG) << "°\n";
    std::cout << "  Ω = " << (kepInit.Omega * RAD_TO_DEG) << "°\n";
    std::cout << "  ω = " << (kepInit.omega * RAD_TO_DEG) << "°\n";
    std::cout << "  M = " << (kepInit.M * RAD_TO_DEG) << "°\n\n";
    
    // ================================================================
    // STEP 2: Carica osservazioni da file .rwo
    // ================================================================
    std::cout << "STEP 2: Caricamento osservazioni da 17030_astdys.rwo\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    
    std::vector<AstrometricObservation> observations;
    
    try {
        observations = parseRWO("17030_astdys.rwo");
        std::cout << "✓ Caricate " << observations.size() << " osservazioni\n";
        
        if (!observations.empty()) {
            std::cout << "  Prima osservazione: JD " << std::fixed << std::setprecision(2) 
                      << observations.front().epoch.jd << "\n";
            std::cout << "  Ultima osservazione: JD " << observations.back().epoch.jd << "\n";
            
            double span = observations.back().epoch.jd - observations.front().epoch.jd;
            std::cout << "  Arco temporale: " << std::setprecision(1) << span << " giorni"
                      << " (" << std::setprecision(2) << (span/365.25) << " anni)\n\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "❌ Errore caricamento osservazioni: " << e.what() << "\n";
        std::cerr << "Continuo senza fitting...\n\n";
    }
    
    // ================================================================
    // STEP 3: Orbit Fitting (se abbiamo osservazioni)
    // ================================================================
    OrbitalElements fittedElements = kepInit;
    bool fittingSuccess = false;
    
    if (observations.size() >= 3) {
        std::cout << "STEP 3: Orbit Fitting con " << observations.size() << " osservazioni\n";
        std::cout << "─────────────────────────────────────────────────────────\n";
        
        try {
            OrbitFitter fitter;
            
            // Configura fitter
            OrbitFitOptions opts;
            opts.maxIterations = 50;
            opts.convergenceThreshold = 0.01;
            opts.rejectOutliers = true;
            opts.outlierSigma = 3.0;
            opts.includeJupiter = true;
            opts.includeSaturn = true;
            
            std::cout << "Configurazione fitter:\n";
            std::cout << "  Max iterazioni: " << opts.maxIterations << "\n";
            std::cout << "  Convergenza: " << opts.convergenceThreshold << " arcsec\n";
            std::cout << "  Outlier rejection: " << (opts.rejectOutliers ? "SÌ" : "NO") << "\n";
            std::cout << "  Perturbatori: Jupiter, Saturn\n\n";
            
            std::cout << "Fitting in corso...\n";
            
            // Crea ObservationSet
            ObservationSet obsSet;
            obsSet.objectDesignation = "17030";
            obsSet.observations = observations;
            obsSet.firstObservation = observations.front().epoch;
            obsSet.lastObservation = observations.back().epoch;
            obsSet.numberOfObservations = observations.size();
            
            // Esegui fitting
            auto result = fitter.fit(kepInit, obsSet, opts);
            
            if (result.converged) {
                fittedElements = result.fittedElements;
                fittingSuccess = true;
                
                std::cout << "✓ Fitting CONVERGITO!\n";
                std::cout << "  Iterazioni: " << result.iterations << "\n";
                std::cout << "  RMS residui: " << std::fixed << std::setprecision(3) 
                          << result.rmsResidual << " arcsec\n";
                std::cout << "  Chi-squared: " << std::setprecision(2) << result.chi2 << "\n\n";
                
                std::cout << "Elementi fitted:\n";
                std::cout << "  Epoca = JD " << std::fixed << std::setprecision(1) << fittedElements.epoch.jd << "\n";
                std::cout << "  a = " << std::setprecision(8) << fittedElements.a << " AU\n";
                std::cout << "  e = " << std::setprecision(6) << fittedElements.e << "\n";
                std::cout << "  i = " << std::setprecision(4) << (fittedElements.i * RAD_TO_DEG) << "°\n";
                std::cout << "  Ω = " << (fittedElements.Omega * RAD_TO_DEG) << "°\n";
                std::cout << "  ω = " << (fittedElements.omega * RAD_TO_DEG) << "°\n";
                std::cout << "  M = " << (fittedElements.M * RAD_TO_DEG) << "°\n\n";
                
                std::cout << "Differenze (Fitted - Iniziale):\n";
                std::cout << "  Δa = " << std::scientific << std::setprecision(3) 
                          << (fittedElements.a - kepInit.a) << " AU\n";
                std::cout << "  Δe = " << (fittedElements.e - kepInit.e) << "\n";
                std::cout << "  Δi = " << std::fixed << std::setprecision(6) 
                          << ((fittedElements.i - kepInit.i) * RAD_TO_DEG) << "°\n";
                std::cout << "  ΔΩ = " << ((fittedElements.Omega - kepInit.Omega) * RAD_TO_DEG) << "°\n";
                std::cout << "  Δω = " << ((fittedElements.omega - kepInit.omega) * RAD_TO_DEG) << "°\n\n";
                
            } else {
                std::cout << "⚠️  Fitting NON convergito\n";
                std::cout << "  Iterazioni: " << result.iterations << "\n";
                std::cout << "  RMS residui: " << std::fixed << std::setprecision(3) 
                          << result.rmsResidual << " arcsec\n\n";
                std::cout << "Uso elementi iniziali...\n\n";
            }
            
        } catch (const std::exception& e) {
            std::cout << "❌ Errore durante fitting: " << e.what() << "\n";
            std::cout << "Uso elementi iniziali...\n\n";
        }
    } else {
        std::cout << "STEP 3: Orbit Fitting SALTATO\n";
        std::cout << "─────────────────────────────────────────────────────────\n";
        std::cout << "⚠️  Servono almeno 3 osservazioni, ne abbiamo " << observations.size() << "\n\n";
    }
    
    // ================================================================
    // STEP 4: Propaga a epoca target (28 Nov 2025)
    // ================================================================
    std::cout << "STEP 4: Propagazione a epoca target (28 Nov 2025)\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    
    JulianDate targetEpoch;
    targetEpoch.jd = 2461007.5;  // 28 Nov 2025
    
    std::cout << "Epoca target: JD " << std::fixed << std::setprecision(1) << targetEpoch.jd << "\n";
    std::cout << "             (2025-11-28 00:00 UTC)\n\n";
    
    // Propaga con elementi iniziali
    std::cout << "a) Propagazione con elementi INIZIALI (AstDyS):\n";
    EquinoctialElements eqInit = kepInit.toEquinoctial();
    Ephemeris eph1(eqInit);
    eph1.enableNumericalPropagation(true, IntegratorType::RKF78);
    
    PropagatorOptions opts1;
    opts1.integrator = IntegratorType::RKF78;
    opts1.stepSize = 0.01;
    opts1.tolerance = 1e-12;
    opts1.usePlanetaryPerturbations = true;
    eph1.setPropagatorOptions(opts1);
    
    EphemerisData data1 = eph1.compute(targetEpoch);
    
    std::cout << "   RA  = " << std::fixed << std::setprecision(6) 
              << (data1.geocentricPos.ra * RAD_TO_DEG) << "°\n";
    std::cout << "   Dec = " << (data1.geocentricPos.dec * RAD_TO_DEG) << "°\n";
    std::cout << "   Dist = " << std::setprecision(8) << data1.distance << " AU\n\n";
    
    // Propaga con elementi fitted (se disponibili)
    if (fittingSuccess) {
        std::cout << "b) Propagazione con elementi FITTED:\n";
        EquinoctialElements eqFit = fittedElements.toEquinoctial();
        Ephemeris eph2(eqFit);
        eph2.enableNumericalPropagation(true, IntegratorType::RKF78);
        eph2.setPropagatorOptions(opts1);
        
        EphemerisData data2 = eph2.compute(targetEpoch);
        
        std::cout << "   RA  = " << std::fixed << std::setprecision(6) 
                  << (data2.geocentricPos.ra * RAD_TO_DEG) << "°\n";
        std::cout << "   Dec = " << (data2.geocentricPos.dec * RAD_TO_DEG) << "°\n";
        std::cout << "   Dist = " << std::setprecision(8) << data2.distance << " AU\n\n";
        
        std::cout << "Differenza (Fitted - Iniziale):\n";
        double deltaRA = (data2.geocentricPos.ra - data1.geocentricPos.ra) * RAD_TO_DEG * 3600.0;
        double deltaDec = (data2.geocentricPos.dec - data1.geocentricPos.dec) * RAD_TO_DEG * 3600.0;
        std::cout << "   ΔRA  = " << std::fixed << std::setprecision(2) << deltaRA << " arcsec\n";
        std::cout << "   ΔDec = " << deltaDec << " arcsec\n";
        std::cout << "   Separazione = " << std::setprecision(1) 
                  << sqrt(deltaRA*deltaRA + deltaDec*deltaDec) << " arcsec\n\n";
    }
    
    // ================================================================
    // STEP 5: Confronto con stella target
    // ================================================================
    std::cout << "STEP 5: Confronto con stella target\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    std::cout << "Stella: Gaia DR3 3411546266140512128\n";
    std::cout << "  RA  = 73.416106°\n";
    std::cout << "  Dec = 20.331662°\n\n";
    
    double star_ra = 73.416106;
    double star_dec = 20.331662;
    
    double sep1 = sqrt(pow(data1.geocentricPos.ra * RAD_TO_DEG - star_ra, 2) +
                       pow(data1.geocentricPos.dec * RAD_TO_DEG - star_dec, 2)) * 3600.0;
    
    std::cout << "Separazione asteroide-stella:\n";
    std::cout << "  Elementi iniziali: " << std::fixed << std::setprecision(1) << sep1 << " arcsec\n";
    
    if (fittingSuccess) {
        EquinoctialElements eqFit = fittedElements.toEquinoctial();
        Ephemeris eph2(eqFit);
        eph2.enableNumericalPropagation(true, IntegratorType::RKF78);
        eph2.setPropagatorOptions(opts1);
        EphemerisData data2 = eph2.compute(targetEpoch);
        
        double sep2 = sqrt(pow(data2.geocentricPos.ra * RAD_TO_DEG - star_ra, 2) +
                          pow(data2.geocentricPos.dec * RAD_TO_DEG - star_dec, 2)) * 3600.0;
        
        std::cout << "  Elementi fitted:   " << std::fixed << std::setprecision(1) << sep2 << " arcsec\n";
        std::cout << "\n  Miglioramento: " << std::setprecision(1) << (sep1 - sep2) << " arcsec\n";
    }
    
    std::cout << "\n═══════════════════════════════════════════════════════════\n";
    std::cout << "TEST COMPLETATO\n";
    std::cout << "═══════════════════════════════════════════════════════════\n";
    
    return 0;
}
