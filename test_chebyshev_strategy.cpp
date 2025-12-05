/**
 * @file test_chebyshev_strategy.cpp
 * @brief Test implementazione FASE 1 con polinomi di Chebyshev
 * @author Michele Bigi
 * @date 30 Novembre 2025
 */

#include "ioccultcalc/propagation_strategy.h"
#include "ioccultcalc/astdyn_interface.h"
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace ioccultcalc;

int main() {
    std::cout << "=== TEST STRATEGIA CHEBYSHEV FASE 1 ===" << std::endl;
    
    try {
        // Elementi orbitali 17030 Sierks
        AstDySElements elements;
        elements.a = 2.794282840;      // Semi-asse maggiore [AU]
        elements.e = 0.120508410;      // Eccentricità
        elements.i = 8.740394000;      // Inclinazione [gradi]
        elements.Omega = 165.542610;   // Nodo ascendente [gradi]
        elements.omega = 343.121070;   // Argomento perielio [gradi]
        elements.M = 128.656770;       // Anomalia media [gradi] @ epoca
        elements.epoch_mjd = 60000.0;  // 2023-02-25.0 TDB
        
        std::cout << "Asteroide: 17030 Sierks" << std::endl;
        std::cout << "Epoca elementi: MJD " << std::fixed << std::setprecision(1) << elements.epoch_mjd << std::endl;
        
        // Inizializza strategia
        PropagationConfig config;
        config.rkf78_tolerance = 1e-12;
        config.verbose_timing = true;
        
        TwoPhaseStrategy strategy(config);
        strategy.setElements(elements);
        
        // Epoca test: 28 Nov 2025 00:35 UTC (MJD 60006.024306)
        double target_mjd = 60006.024306;
        std::cout << "\nEpoca target: MJD " << std::fixed << std::setprecision(6) << target_mjd << std::endl;
        
        // Stella Gaia DR3 3411546266140512128
        double ra_star = 40.8652440;   // RA [deg]
        double dec_star = 11.8825320;  // Dec [deg]
        std::cout << "Stella: RA=" << std::fixed << std::setprecision(6) << ra_star 
                  << "°, Dec=" << dec_star << "°" << std::endl;
        
        // Test FASE 1: Screening con Chebyshev
        std::cout << "\n=== FASE 1: Screening Chebyshev ===" << std::endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        double separation_fase1 = strategy.screenCandidate(target_mjd, ra_star, dec_star);
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        std::cout << "Separazione FASE 1: " << std::fixed << std::setprecision(2) 
                  << separation_fase1 << "\" (tempo: " << duration.count() << " μs)" << std::endl;
        
        // Test FASE 2: Closest approach con RKF78
        std::cout << "\n=== FASE 2: Closest Approach RKF78 ===" << std::endl;
        
        start = std::chrono::high_resolution_clock::now();
        auto result = strategy.findClosestApproach(target_mjd, ra_star, dec_star, 2.0);
        end = std::chrono::high_resolution_clock::now();
        
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        std::cout << "Tempo closest approach: MJD " << std::fixed << std::setprecision(6) 
                  << result.closest_time_mjd << std::endl;
        std::cout << "Separazione minima: " << std::fixed << std::setprecision(2) 
                  << result.minimum_separation_arcsec << "\" (tempo: " << duration.count() << " ms)" << std::endl;
        std::cout << "Position angle: " << std::fixed << std::setprecision(1) 
                  << result.position_angle_deg << "°" << std::endl;
        
        // Statistiche performance
        auto stats = strategy.getStats();
        std::cout << "\n=== STATISTICHE PERFORMANCE ===" << std::endl;
        std::cout << "Chiamate FASE 1: " << stats.phase1_calls << std::endl;
        std::cout << "Chiamate FASE 2: " << stats.phase2_calls << std::endl;
        std::cout << "Tempo totale FASE 1: " << std::fixed << std::setprecision(3) 
                  << stats.phase1_total_time_ms << " ms" << std::endl;
        std::cout << "Tempo totale FASE 2: " << std::fixed << std::setprecision(3) 
                  << stats.phase2_total_time_ms << " ms" << std::endl;
        
        // Verifica accuratezza: confronto con valutazione diretta FASE 2
        std::cout << "\n=== VERIFICA ACCURATEZZA ===" << std::endl;
        
        // Usa FASE 2 direttamente come riferimento (è già RKF78 accurato)  
        double reference_separation = result.minimum_separation_arcsec;
        
        std::cout << "FASE 2 (riferimento): " << std::fixed << std::setprecision(2) 
                  << reference_separation << "\"" << std::endl;
        std::cout << "Differenza FASE 1-FASE 2: " << std::fixed << std::setprecision(3) 
                  << std::abs(separation_fase1 - reference_separation) << "\"" << std::endl;
        
        if (std::abs(separation_fase1 - reference_separation) < 10.0) {  // 10" tolleranza per Chebyshev
            std::cout << "✓ FASE 1 Chebyshev accettabile (<10\" errore)" << std::endl;
        } else {
            std::cout << "✗ FASE 1 Chebyshev errore troppo grande" << std::endl;
        }
        
        std::cout << "\n=== TEST COMPLETATO ===" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << std::endl;
        return 1;
    }
}