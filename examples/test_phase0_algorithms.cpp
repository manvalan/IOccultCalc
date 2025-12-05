/**
 * @brief Test del sistema Phase 0 con diversi threshold
 */
#include <iostream>
#include <ioccultcalc/star_screening.h>
#include <ioccultcalc/astdys_client.h>
#include <ioccultcalc/time_utils.h>

using namespace ioccultcalc;

int main() {
    std::cout << "🔬 Phase 0 Algorithm Test Suite\n\n";
    
    try {
        // Setup basic configuration
        ScreeningConfig config;
        config.useFastPreScreening = true;
        config.fastThreshold = 120.0;        // 2 arcmin threshold
        config.fastGridPoints = 20;
        config.useApproximateDistance = true;
        config.useDetailedLogging = true;
        
        // Get asteroid elements
        AstDysClient client;
        auto elements = client.getElements("433"); // Eros
        
        // Setup screening system
        StarScreening screening(config);
        
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RK4;
        opts.stepSize = 0.1;
        
        screening.setAsteroidElements(elements, opts);
        
        // Test times
        JulianDate startTime = TimeUtils::isoToJD("2025-01-01");
        JulianDate endTime = TimeUtils::isoToJD("2025-01-31");
        
        // Test 1: Very distant star (should be rejected by Phase 0)
        std::cout << "📊 Test 1: Distant star (RA=45°, Dec=23°)\n";
        screening.setStarCoordinates(45.5 * M_PI / 180.0, 23.8 * M_PI / 180.0);
        
        auto fastResult1 = screening.performFastPreScreening(startTime, endTime);
        std::cout << "   Result: " << (fastResult1.isCandidate ? "CANDIDATE" : "REJECTED") << "\n";
        std::cout << "   Min distance: " << fastResult1.dMinApprox << " arcsec\n";
        std::cout << "   Scan time: " << fastResult1.scanTime * 1000 << " ms\n\n";
        
        // Test 2: Complete 3-phase analysis
        std::cout << "📊 Test 2: Complete 3-Phase Analysis\n";
        auto completeResult = screening.performCompleteAnalysis(startTime, endTime, 150000.0, 15.0);
        std::cout << "   Occultation possible: " << (completeResult.occultationPossible ? "YES" : "NO") << "\n";
        std::cout << "   Min distance: " << completeResult.dMinArcsec << " arcsec\n\n";
        
        // Test 3: Performance comparison
        std::cout << "📊 Test 3: Performance Comparison\n";
        
        auto start = std::chrono::high_resolution_clock::now();
        screening.performFastPreScreening(startTime, endTime);
        auto end = std::chrono::high_resolution_clock::now();
        double phase0Time = std::chrono::duration<double>(end - start).count();
        
        start = std::chrono::high_resolution_clock::now();
        screening.screenTimeInterval(startTime, endTime);
        end = std::chrono::high_resolution_clock::now();
        double fullTime = std::chrono::duration<double>(end - start).count();
        
        std::cout << "   Phase 0 time: " << phase0Time * 1000 << " ms\n";
        std::cout << "   Full screening time: " << fullTime * 1000 << " ms\n";
        std::cout << "   Speedup: " << (fullTime / phase0Time) << "x\n\n";
        
        std::cout << "✅ Phase 0 tests completed successfully!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Test failed: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}