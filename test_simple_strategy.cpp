#include "include/ioccultcalc/propagation_strategy.h"
#include <iostream>
#include <iomanip>

using namespace ioccultcalc;

int main() {
    std::cout << "Test strategia semplificato\n\n";
    
    // Elementi corretti dal file AstDyS .eq1
    AstDySElements elem;
    elem.name = "17030";
    elem.number = 17030;
    elem.epoch_mjd = 61000.0;
    elem.a = 3.1754732060579491;
    elem.e = 0.0454207; // Calcolato da h e k
    elem.i = 2.9046;
    elem.Omega = 104.185;
    elem.omega = 102.150;
    elem.M = 99.035;  // Dal file MPC che funzionava
    elem.H = 13.29;
    elem.G = 0.15;
    elem.has_covariance = false;
    
    PropagationConfig config = propagation_presets::createBalanced();
    config.screening_threshold_arcsec = 100.0; // Soglia più permissiva
    
    TwoPhaseStrategy strategy(config);
    strategy.setElements(elem);
    
    double target_mjd = 61007.024306;
    double star_ra = 73.416167;
    double star_dec = 20.332083;
    
    std::cout << "Target: MJD " << target_mjd << "\n";
    std::cout << "Stella: RA=" << star_ra << "°, Dec=" << star_dec << "°\n\n";
    
    // Test FASE 1
    double sep = strategy.screenCandidate(target_mjd, star_ra, star_dec);
    std::cout << "FASE 1 separazione: " << std::fixed << std::setprecision(1) << sep << "\"\n";
    
    if (sep < config.screening_threshold_arcsec) {
        std::cout << "✅ Candidato promosso alla FASE 2\n\n";
        
        auto result = strategy.findClosestApproach(target_mjd, star_ra, star_dec);
        std::cout << "FASE 2 risultato:\n";
        std::cout << "  Tempo CA: MJD " << std::setprecision(6) << result.closest_time_mjd << "\n";
        std::cout << "  Separazione min: " << std::setprecision(1) << result.minimum_separation_arcsec << "\"\n";
        std::cout << "  Position angle: " << result.position_angle_deg << "°\n";
        std::cout << "  Tempo calcolo: " << result.computation_time_ms << " ms\n";
    } else {
        std::cout << "❌ Candidato scartato (soglia=" << config.screening_threshold_arcsec << "\")\n";
    }
    
    auto stats = strategy.getStats();
    std::cout << "\nStatistiche: " << stats.phase1_calls << " fase1, " << stats.phase2_calls << " fase2\n";
    
    return 0;
}
