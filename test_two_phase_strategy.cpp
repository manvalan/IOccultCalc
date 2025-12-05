/**
 * @file test_two_phase_strategy.cpp
 * @brief Test della strategia di propagazione a due fasi semplificata
 * @author Michele Bigi
 * @date 30 Novembre 2025
 */

#include "include/ioccultcalc/propagation_strategy.h"
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace ioccultcalc;

int main() {
    std::cout << "================================================================\n";
    std::cout << " TEST STRATEGIA DUE FASI SEMPLIFICATA\n";
    std::cout << " Asteroide 17030 Sierks - 28 Nov 2025 00:35 UTC\n";
    std::cout << "================================================================\n\n";
    
    // Elementi asteroide 17030 da AstDyS (corretti dal test precedente)
    AstDySElements elem;
    elem.name = "17030";
    elem.number = 17030;
    elem.epoch_mjd = 61000.0; // 2025-11-21
    elem.a = 3.17547;
    elem.e = 0.0454207;
    elem.i = 2.9046;
    elem.Omega = 104.185;
    elem.omega = 102.150;
    elem.M = 225.123;  // Questa potrebbe essere il problema
    elem.H = 13.33;
    elem.G = 0.15;
    elem.has_covariance = false;
    
    double target_mjd = 61007.024306; // 28 Nov 2025 00:35 UTC
    double star_ra = 73.416167;       // Stella target
    double star_dec = 20.332083;
    
    std::cout << "📂 ELEMENTI ASTEROIDE 17030:\n";
    std::cout << "   Epoca: MJD " << elem.epoch_mjd << " (2025-11-21)\n";
    std::cout << "   a = " << elem.a << " AU\n";
    std::cout << "   e = " << elem.e << "\n";
    std::cout << "   i = " << elem.i << "°\n\n";
    
    std::cout << "🎯 TARGET:\n";
    std::cout << "   Data: MJD " << target_mjd << " (2025-11-28 00:35 UTC)\n";
    std::cout << "   Stella: RA=" << star_ra << "°, Dec=" << star_dec << "°\n\n";
    
    // Test diverse configurazioni
    std::cout << "================================================================\n";
    std::cout << " TEST CONFIGURAZIONI\n";
    std::cout << "================================================================\n\n";
    
    auto test_config = [&](const std::string& name, const PropagationConfig& config) {
        std::cout << "🔧 " << name << ":\n";
        
        TwoPhaseStrategy strategy(config);
        strategy.setElements(elem);
        
        auto start_total = std::chrono::high_resolution_clock::now();
        
        // FASE 1: Screening
        std::cout << "   FASE 1 (Screening Kepleriano)...\n";
        double separation_phase1 = strategy.screenCandidate(target_mjd, star_ra, star_dec);
        
        std::cout << "   → Separazione: " << std::fixed << std::setprecision(2) 
                  << separation_phase1 << "\"\n";
        
        // Check se passa la soglia
        bool promoted = separation_phase1 < config.screening_threshold_arcsec;
        std::cout << "   → Soglia: " << config.screening_threshold_arcsec 
                  << "\" → " << (promoted ? "✅ PROMOSSO" : "❌ SCARTATO") << "\n\n";
        
        if (promoted) {
            // FASE 2: Closest approach preciso
            std::cout << "   FASE 2 (Closest Approach RKF78)...\n";
            auto result = strategy.findClosestApproach(target_mjd, star_ra, star_dec, 2.0);
            
            std::cout << "   → Tempo CA: MJD " << std::fixed << std::setprecision(6) 
                      << result.closest_time_mjd << "\n";
            std::cout << "   → Separazione min: " << std::setprecision(2) 
                      << result.minimum_separation_arcsec << "\"\n";
            std::cout << "   → Position angle: " << std::setprecision(1) 
                      << result.position_angle_deg << "°\n";
            std::cout << "   → Steps: " << result.computation_steps << "\n";
            std::cout << "   → Tempo calcolo: " << result.computation_time_ms << " ms\n";
            
            if (result.orbit_was_fitted) {
                std::cout << "   → 🔄 Orbit fitting: " 
                          << result.fitting_iterations << " iter, RMS=" 
                          << result.final_rms_arcsec << "\"\n";
            }
        }
        
        auto end_total = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_total - start_total);
        
        // Statistiche
        auto stats = strategy.getStats();
        std::cout << "\n   📊 STATISTICHE:\n";
        std::cout << "      Fase 1: " << stats.phase1_calls << " calls, " 
                  << stats.phase1_total_time_ms << " ms\n";
        std::cout << "      Fase 2: " << stats.phase2_calls << " calls, " 
                  << stats.phase2_total_time_ms << " ms\n";
        std::cout << "      Totale: " << total_time.count() << " ms\n";
        std::cout << "      Candidati: " << stats.candidates_screened 
                  << " screened, " << stats.candidates_promoted << " promoted\n\n";
    };
    
    // Test configurazione veloce
    test_config("CONFIGURAZIONE VELOCE (Survey)", propagation_presets::createFastSurvey());
    
    // Test configurazione bilanciata
    test_config("CONFIGURAZIONE BILANCIATA", propagation_presets::createBalanced());
    
    // Test configurazione precisione
    test_config("CONFIGURAZIONE PRECISIONE", propagation_presets::createPrecision());
    
    // Test performance mass screening
    std::cout << "================================================================\n";
    std::cout << " TEST PERFORMANCE SCREENING MASSA\n";
    std::cout << "================================================================\n\n";
    
    PropagationConfig fast_config = propagation_presets::createFastSurvey();
    TwoPhaseStrategy fast_strategy(fast_config);
    fast_strategy.setElements(elem);
    
    const int num_stars = 1000;
    std::cout << "🚀 Screening " << num_stars << " stelle candidat...\n";
    
    auto start_mass = std::chrono::high_resolution_clock::now();
    
    int promoted = 0;
    for (int i = 0; i < num_stars; ++i) {
        // Stelle casuali attorno alla posizione target
        double ra_offset = (i % 100 - 50) * 0.01;  // ±0.5°
        double dec_offset = ((i / 100) % 100 - 50) * 0.01;
        
        double test_ra = star_ra + ra_offset;
        double test_dec = star_dec + dec_offset;
        
        double sep = fast_strategy.screenCandidate(target_mjd, test_ra, test_dec);
        if (sep < fast_config.screening_threshold_arcsec) {
            promoted++;
        }
    }
    
    auto end_mass = std::chrono::high_resolution_clock::now();
    auto mass_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_mass - start_mass);
    
    auto final_stats = fast_strategy.getStats();
    
    std::cout << "\n✅ RISULTATI SCREENING MASSA:\n";
    std::cout << "   Stelle processate: " << num_stars << "\n";
    std::cout << "   Candidate promosse: " << promoted << " (" 
              << (100.0 * promoted / num_stars) << "%)\n";
    std::cout << "   Tempo totale: " << mass_time.count() << " ms\n";
    std::cout << "   Velocità: " << std::fixed << std::setprecision(1) 
              << (num_stars * 1000.0 / mass_time.count()) << " stelle/sec\n";
    std::cout << "   Tempo medio per stella: " << std::setprecision(3)
              << (final_stats.phase1_total_time_ms / final_stats.phase1_calls) 
              << " ms\n\n";
    
    std::cout << "💡 ANALISI:\n";
    std::cout << "   • FASE 1 (Kepleriano): Ultra-veloce per screening\n";
    std::cout << "   • FASE 2 (RKF78): Precisione per closest approach\n";
    std::cout << "   • Orbit fitting: Opzionale per migliorare accuracy\n";
    std::cout << "   • Scalabilità: >1000 stelle/sec per survey\n\n";
    
    std::cout << "✅ Test strategia a due fasi completato!\n";
    return 0;
}
