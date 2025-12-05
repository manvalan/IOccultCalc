/**
 * @file test_star_screening.cpp
 * @brief Test unitario per la classe StarScreening
 */

#include <iostream>
#include <ioccultcalc/star_screening.h>
#include <ioccultcalc/astdys_client.h>
#include <ioccultcalc/time_utils.h>
#include <ioccultcalc/chebyshev_astrometry.h>

using namespace ioccultcalc;

int main() {
    std::cout << "🧪 Test StarScreening Class\n";
    std::cout << "============================\n\n";
    
    try {
        // Test 1: Configurazione base
        std::cout << "Test 1: Basic configuration...\n";
        
        ScreeningConfig config;
        config.stepSize = 1.0;  // 1 giorno
        config.candidateThreshold = 30.0;  // 30 arcsec
        config.useDetailedLogging = true;
        
        StarScreening screening(config);
        
        // Verifica stato iniziale
        std::cout << "Initial status:\n" << screening.getStatusString() << "\n\n";
        
        // Test 2: Setup elementi orbitali (asteroid 17030)
        std::cout << "Test 2: Setting up asteroid elements...\n";
        
        // Elementi semplificati per test (17030 in approssimazione)
        EquinoctialElements elements;
        elements.epoch = JulianDate(2460611.5);  // 2025-01-01
        elements.a = 3.175;     // AU
        elements.lambda = 2.0;  // Longitudine media
        elements.k = 0.02;      // e*cos(omega)
        elements.h = 0.03;      // e*sin(omega)
        elements.p = 0.025;     // tan(i/2)*sin(Omega)
        elements.q = 0.015;     // tan(i/2)*cos(Omega)
        elements.name = "17030";
        
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RK4;
        opts.stepSize = 0.1;
        opts.usePlanetaryPerturbations = true;
        
        screening.setAsteroidElements(elements, opts);
        std::cout << "✓ Asteroid elements set\n";
        
        // Test 3: Setup stella
        std::cout << "\nTest 3: Setting up star coordinates...\n";
        double starRA = 15.5 * M_PI / 180.0;   // 15.5 gradi
        double starDec = 45.2 * M_PI / 180.0;  // 45.2 gradi
        
        screening.setStarCoordinates(starRA, starDec);
        std::cout << "✓ Star coordinates set: RA=" << (starRA * 180.0 / M_PI) 
                 << "°, Dec=" << (starDec * 180.0 / M_PI) << "°\n";
        
        // Test 4: Setup proprietà asteroide
        std::cout << "\nTest 4: Setting up asteroid properties...\n";
        screening.setAsteroidProperties(25.0, 13.5);  // 25 km, H=13.5
        std::cout << "✓ Asteroid properties set: D=25 km, H=13.5\n";
        
        // Test 5: Verifica configurazione completa
        std::cout << "\nTest 5: Configuration verification...\n";
        std::cout << "Final status:\n" << screening.getStatusString() << "\n\n";
        
        if (!screening.isConfigured()) {
            std::cerr << "❌ Configuration incomplete!\n";
            return 1;
        }
        
        // Test 6: Screening singola epoca
        std::cout << "Test 6: Single epoch screening...\n";
        JulianDate testEpoch(2460611.5 + 10.0);  // +10 giorni
        
        ScreeningResult result = screening.screenSingleEpoch(testEpoch);
        
        std::cout << "Single epoch result:\n";
        std::cout << "  Epoch: " << TimeUtils::jdToISO(result.epoch) << "\n";
        std::cout << "  Angular separation: " << result.angularSeparation << " arcsec\n";
        std::cout << "  Distance: " << result.asteroidDistance << " AU\n";
        std::cout << "  Is candidate: " << (result.isCandidate ? "Yes" : "No") << "\n";
        std::cout << "  Is occultation: " << (result.isOccultation ? "Yes" : "No") << "\n";
        std::cout << "  Probability: " << (result.probability * 100.0) << "%\n";
        std::cout << "  Used Chebyshev: " << (result.usedChebyshev ? "Yes" : "No") << "\n\n";
        
        // Test 7: Screening intervallo breve
        std::cout << "Test 7: Time interval screening (short)...\n";
        
        JulianDate startTime(2460611.5);
        JulianDate endTime(2460616.5);  // 5 giorni
        
        auto candidates = screening.screenTimeInterval(startTime, endTime);
        
        std::cout << "Interval screening results:\n";
        std::cout << "  Candidates found: " << candidates.size() << "\n";
        
        const auto& stats = screening.getLastStatistics();
        std::cout << "  Total calculations: " << stats.totalCalculations << "\n";
        std::cout << "  Computation time: " << stats.computationTime << " sec\n";
        std::cout << "  Average error: " << stats.averageError << " arcsec\n";
        
        if (!candidates.empty()) {
            std::cout << "\n  Best candidate:\n";
            const auto& best = candidates[0];
            std::cout << "    Date: " << TimeUtils::jdToISO(best.epoch) << "\n";
            std::cout << "    Separation: " << best.angularSeparation << " arcsec\n";
            std::cout << "    Probability: " << (best.probability * 100.0) << "%\n";
        }
        
        std::cout << "\n";
        
        // Test 8: Test con modello Chebyshev (opzionale)
        std::cout << "Test 8: Chebyshev model integration (optional)...\n";
        
        // Crea un modello Chebyshev di test
        ChebyshevAstrometry chebyshev(3, 3);  // Ordine 3x3
        
        // Coefficienti di test (campo simulato)
        ChebyshevCoefficients coeffs;
        coeffs.alpha_base = starRA;   // Centro campo = posizione stella
        coeffs.delta_base = starDec;
        coeffs.max_order_u = 3;
        coeffs.max_order_v = 3;
        
        // Coefficienti di test (identità per test base)
        coeffs.a_coefficients.resize(4, std::vector<double>(4, 0.0));
        coeffs.b_coefficients.resize(4, std::vector<double>(4, 0.0));
        coeffs.a_coefficients[0][0] = 0.0;  // Termine costante
        coeffs.b_coefficients[0][0] = 0.0;
        coeffs.a_coefficients[1][0] = 0.01 * M_PI / 180.0;  // Termine lineare
        coeffs.b_coefficients[0][1] = 0.01 * M_PI / 180.0;
        
        chebyshev.setCoefficients(coeffs);
        
        if (chebyshev.isValid()) {
            screening.setChebyshevModel(&chebyshev);
            std::cout << "✓ Chebyshev model set and valid\n";
            
            // Ripeti screening con Chebyshev
            ScreeningResult chebyResult = screening.screenSingleEpoch(testEpoch);
            
            std::cout << "  Chebyshev result:\n";
            std::cout << "    Angular separation: " << chebyResult.angularSeparation << " arcsec\n";
            std::cout << "    Astrometric error: " << chebyResult.astrometricError << " arcsec\n";
            std::cout << "    Used Chebyshev: " << (chebyResult.usedChebyshev ? "Yes" : "No") << "\n";
            
        } else {
            std::cout << "⚠️  Chebyshev model not valid (expected for test)\n";
        }
        
        std::cout << "\n✅ All tests completed successfully!\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Test failed: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}