// Test corridor optimization for occultation prediction
// Compare performance with and without corridor query

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstdlib>
#include "ioccultcalc/occultation_predictor.h"
#include "ioccultcalc/time_utils.h"
#include "ioccultcalc/astdys_client.h"

using namespace ioccultcalc;

int main() {
    std::cout << "=== Test Corridor Optimization for Occultation Prediction ===\n\n";
    
    // Test con (17030) Sierks - evento noto 28 Nov 2025 ~00:35 UT
    std::cout << "Test case: (17030) Sierks occultation\n";
    std::cout << "Expected: 28 Nov 2025 ~00:35 UT\n";
    std::cout << "Target star: Gaia DR3 3411546266140512128 (G=12.1)\n\n";
    
    // Carica elementi orbitali da AstDyS
    AstDysClient astdys;
    EquinoctialElements sierks;
    
    try {
        std::cout << "Loading orbital elements for (17030) Sierks...\n";
        sierks = astdys.getElements("17030");
        std::cout << "  a = " << sierks.a << " AU\n";
        std::cout << "  h = " << sierks.h << "\n";
        std::cout << "  k = " << sierks.k << "\n";
        std::cout << "  Epoch = JD " << sierks.epoch.jd << "\n\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR loading elements: " << e.what() << "\n";
        
        // Fallback: usa elementi approssimati
        std::cout << "Using approximate elements...\n";
        sierks.a = 2.77;
        sierks.h = 0.05;
        sierks.k = 0.12;
        sierks.p = 0.05;
        sierks.q = 0.08;
        sierks.lambda = 2.3;
        sierks.epoch = JulianDate(2460600.5);  // ~ Nov 2024
        sierks.designation = "17030";
    }
    
    // Periodo di ricerca: 25-30 Nov 2025
    JulianDate startJD = TimeUtils::calendarToJD(2025, 11, 25, 0.0);
    JulianDate endJD = TimeUtils::calendarToJD(2025, 11, 30, 0.0);
    
    std::cout << "Search period: 25-30 Nov 2025\n";
    std::cout << "  Start JD = " << std::fixed << std::setprecision(2) << startJD.jd << "\n";
    std::cout << "  End JD   = " << endJD.jd << "\n\n";
    
    // Configura predictor
    OccultationPredictor predictor;
    predictor.setAsteroid(sierks);
    predictor.setAsteroidDiameter(5.0);  // ~5 km stimato
    predictor.setOrbitalUncertainty(50.0);  // 50 km 1-sigma
    
    // Parametri di ricerca
    double maxMagnitude = 14.0;    // Stelle fino a mag 14
    double searchRadius = 0.5;     // 0.5° radius
    double minProbability = 0.001; // 0.1% minimo
    
    std::cout << "Search parameters:\n";
    std::cout << "  Max magnitude: " << maxMagnitude << "\n";
    std::cout << "  Search radius: " << searchRadius << "°\n";
    std::cout << "  Min probability: " << (minProbability * 100) << "%\n\n";
    
    // Esegui ricerca
    std::cout << "=== Running Occultation Search ===\n\n";
    
    auto t_start = std::chrono::high_resolution_clock::now();
    
    auto events = predictor.findOccultations(startJD, endJD, maxMagnitude, 
                                             searchRadius, minProbability);
    
    auto t_end = std::chrono::high_resolution_clock::now();
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
    
    std::cout << "\n=== Results ===\n\n";
    std::cout << "Search completed in " << duration_ms << " ms\n";
    std::cout << "Found " << events.size() << " occultation events\n\n";
    
    // Mostra eventi trovati
    if (!events.empty()) {
        std::cout << "Events:\n";
        for (size_t i = 0; i < events.size(); i++) {
            const auto& ev = events[i];
            
            int year, month, day, hour, minute;
            double second;
            TimeUtils::jdToCalendar(ev.timeCA, year, month, day, hour, minute, second);
            
            std::cout << "\n[" << (i+1) << "] Event at " << year << "-" 
                      << std::setw(2) << std::setfill('0') << month << "-"
                      << std::setw(2) << std::setfill('0') << day << " "
                      << std::setw(2) << std::setfill('0') << hour << ":"
                      << std::setw(2) << std::setfill('0') << minute << ":"
                      << std::fixed << std::setprecision(1) << second << " UT\n";
            std::cout << std::setfill(' ');
            std::cout << "    Star: " << ev.star.sourceId << "\n";
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "    Star RA/Dec: " << ev.star.pos.ra << "° / " << ev.star.pos.dec << "°\n";
            std::cout << std::setprecision(2);
            std::cout << "    Star mag: " << ev.star.phot_g_mean_mag << "\n";
            std::cout << std::setprecision(3);
            std::cout << "    Close approach: " << ev.closeApproachDistance << " arcsec\n";
            std::cout << "    Probability: " << (ev.probability * 100) << "%\n";
            std::cout << "    Max duration: " << ev.maxDuration << " sec\n";
            
            // Check se è l'evento target
            if (ev.star.sourceId == "3411546266140512128") {
                std::cout << "\n    *** THIS IS THE TARGET SIERKS EVENT! ***\n";
            }
        }
    } else {
        std::cout << "No events found.\n";
        std::cout << "\nDEBUG: Check that:\n";
        std::cout << "  1. Gaia catalog is properly initialized\n";
        std::cout << "  2. Ephemeris computation is correct\n";
        std::cout << "  3. Search parameters are reasonable\n";
    }
    
    std::cout << "\n=== Test Complete ===\n";
    return 0;
}
