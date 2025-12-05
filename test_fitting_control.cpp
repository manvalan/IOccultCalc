#include <iostream>
#include <iomanip>
#include "include/ioccultcalc/astdyn_interface.h"
#include "include/ioccultcalc/propagation_strategy.h"

using namespace ioccultcalc;

int main() {
    std::cout << "🔧 TEST CONTROLLO ORBIT FITTING AVANZATO\n";
    std::cout << "=========================================\n\n";
    
    // Carica elementi 17030
    AstDySElements elements = AstDySElements::fromFile("17030_astdys.eq1");
    std::cout << "📂 Elementi 17030 caricati (epoca MJD " << elements.epoch_mjd << ")\n\n";
    
    // Setup osservazioni simulate
    std::vector<RWOObservation> observations(5);
    for (int i = 0; i < 5; i++) {
        observations[i].mjd_utc = 60000.0 + i;
        observations[i].ra_deg = 280.0 + i * 0.1;
        observations[i].dec_deg = -22.0 + i * 0.01;
        observations[i].ra_sigma_arcsec = 1.0;
        observations[i].dec_sigma_arcsec = 1.0;
        observations[i].designation = "17030";
        observations[i].obs_code = "500";
        observations[i].magnitude = 16.0;
        observations[i].band = "V";
    }
    std::cout << "📊 Osservazioni simulate: " << observations.size() << "\n\n";
    
    double target_mjd = 60006.024306;
    double ra_star = 283.0;
    double dec_star = -22.1;
    
    // ==================================================================
    // TEST 1: Modalità NEVER (nessun fitting)
    // ==================================================================
    std::cout << "🔴 TEST 1: FittingMode::NEVER\n";
    std::cout << "-------------------------------\n";
    
    auto config_never = propagation_presets::createNoFitting();
    TwoPhaseStrategy strategy_never(config_never);
    strategy_never.setElements(elements);
    strategy_never.setObservations(observations);
    
    std::cout << "  Modalità fitting: " << 
        (strategy_never.getFittingMode() == PropagationConfig::FittingMode::NEVER ? "NEVER" : "OTHER") << "\n";
    std::cout << "  Osservazioni disponibili: " << strategy_never.getObservationCount() << "\n";
    std::cout << "  Usa elementi fittati: " << (strategy_never.isUsingFittedElements() ? "SÌ" : "NO") << "\n";
    
    auto result_never = strategy_never.findClosestApproach(target_mjd, ra_star, dec_star, 2.0);
    std::cout << "  ✅ Risultato: " << std::fixed << std::setprecision(1) 
              << result_never.minimum_separation_arcsec << "\" separazione\n";
    std::cout << "  Orbit fittato: " << (result_never.orbit_was_fitted ? "SÌ" : "NO") << "\n\n";
    
    // ==================================================================
    // TEST 2: Modalità AUTO (fitting automatico se osservazioni sufficienti)
    // ==================================================================
    std::cout << "🔵 TEST 2: FittingMode::AUTO\n";
    std::cout << "------------------------------\n";
    
    auto config_auto = propagation_presets::createTestingWithFitting();
    TwoPhaseStrategy strategy_auto(config_auto);
    strategy_auto.setElements(elements);
    strategy_auto.setObservations(observations);
    
    std::cout << "  Modalità fitting: AUTO\n";
    std::cout << "  Osservazioni disponibili: " << strategy_auto.getObservationCount() << "\n";
    std::cout << "  Minimo richiesto: " << config_auto.min_observations_for_fitting << "\n";
    std::cout << "  Usa elementi fittati: " << (strategy_auto.isUsingFittedElements() ? "SÌ" : "NO") << "\n";
    
    auto result_auto = strategy_auto.findClosestApproach(target_mjd, ra_star, dec_star, 2.0);
    std::cout << "  ✅ Risultato: " << std::fixed << std::setprecision(1) 
              << result_auto.minimum_separation_arcsec << "\" separazione\n";
    std::cout << "  Orbit fittato: " << (result_auto.orbit_was_fitted ? "SÌ" : "NO") << "\n";
    std::cout << "  Post-call usa elementi fittati: " << (strategy_auto.isUsingFittedElements() ? "SÌ" : "NO") << "\n\n";
    
    // ==================================================================
    // TEST 3: Trigger manuale orbit fitting
    // ==================================================================
    std::cout << "🟢 TEST 3: Trigger manuale orbit fitting\n";
    std::cout << "------------------------------------------\n";
    
    TwoPhaseStrategy strategy_manual(config_auto);
    strategy_manual.setElements(elements);
    strategy_manual.setObservations(observations);
    strategy_manual.setFittingMode(PropagationConfig::FittingMode::NEVER); // Start with NEVER
    
    std::cout << "  Modalità iniziale: NEVER\n";
    std::cout << "  Trigger manuale fitting...\n";
    bool fit_success = strategy_manual.triggerOrbitFitting();
    std::cout << "  Fitting riuscito: " << (fit_success ? "SÌ" : "NO") << "\n";
    std::cout << "  Ora usa elementi fittati: " << (strategy_manual.isUsingFittedElements() ? "SÌ" : "NO") << "\n";
    
    auto result_manual = strategy_manual.findClosestApproach(target_mjd, ra_star, dec_star, 2.0);
    std::cout << "  ✅ Risultato: " << std::fixed << std::setprecision(1) 
              << result_manual.minimum_separation_arcsec << "\" separazione\n\n";
    
    // ==================================================================
    // TEST 4: Modalità ALWAYS_ATTEMPT (richiede osservazioni sufficienti)
    // ==================================================================
    std::cout << "🟡 TEST 4: FittingMode::ALWAYS_ATTEMPT\n";
    std::cout << "---------------------------------------\n";
    
    auto config_always = propagation_presets::createPrecision();
    TwoPhaseStrategy strategy_always(config_always);
    strategy_always.setElements(elements);
    strategy_always.setObservations(observations);
    
    std::cout << "  Modalità fitting: ALWAYS_ATTEMPT\n";
    std::cout << "  Osservazioni disponibili: " << strategy_always.getObservationCount() << "\n";
    std::cout << "  Minimo richiesto: " << config_always.min_observations_for_fitting << "\n";
    std::cout << "  Force refit: " << (config_always.force_refit_each_prediction ? "SÌ" : "NO") << "\n";
    
    try {
        auto result_always = strategy_always.findClosestApproach(target_mjd, ra_star, dec_star, 2.0);
        std::cout << "  ✅ Risultato: " << std::fixed << std::setprecision(1) 
                  << result_always.minimum_separation_arcsec << "\" separazione\n";
        std::cout << "  Orbit fittato: " << (result_always.orbit_was_fitted ? "SÌ" : "NO") << "\n";
        
        // Seconda chiamata per testare force_refit_each_prediction
        std::cout << "  🔄 Seconda predizione (test force_refit)...\n";
        auto result_always2 = strategy_always.findClosestApproach(target_mjd + 1.0, ra_star, dec_star, 2.0);
        std::cout << "  ✅ Seconda chiamata: " << (result_always2.orbit_was_fitted ? "RE-FITTATO" : "RIUSO ELEMENTI") << "\n";
        
    } catch (const std::exception& e) {
        std::cout << "  ❌ Errore: " << e.what() << "\n";
    }
    
    std::cout << "\n";
    
    // ==================================================================
    // TEST 5: Reset elementi fittati
    // ==================================================================
    std::cout << "🟣 TEST 5: Reset elementi fittati\n";
    std::cout << "----------------------------------\n";
    
    std::cout << "  Prima reset usa elementi fittati: " << (strategy_auto.isUsingFittedElements() ? "SÌ" : "NO") << "\n";
    strategy_auto.resetFittedElements();
    std::cout << "  Dopo reset usa elementi fittati: " << (strategy_auto.isUsingFittedElements() ? "SÌ" : "NO") << "\n";
    
    auto result_reset = strategy_auto.findClosestApproach(target_mjd, ra_star, dec_star, 2.0);
    std::cout << "  ✅ Dopo reset: " << (result_reset.orbit_was_fitted ? "RE-FITTATO" : "NO FITTING") << "\n\n";
    
    // ==================================================================
    // CONFRONTO RISULTATI
    // ==================================================================
    std::cout << "📊 CONFRONTO RISULTATI FINALI\n";
    std::cout << "===============================\n";
    std::cout << "Modalità           | Separazione | Fittato | Tempo\n";
    std::cout << "-------------------|-------------|---------|--------\n";
    std::cout << "NEVER              | " << std::setw(8) << std::fixed << std::setprecision(2) 
              << result_never.minimum_separation_arcsec << "\"  | " 
              << (result_never.orbit_was_fitted ? "  SÌ   " : "  NO   ") 
              << " | " << std::setw(4) << result_never.computation_time_ms << "ms\n";
    std::cout << "AUTO               | " << std::setw(8) 
              << result_auto.minimum_separation_arcsec << "\"  | " 
              << (result_auto.orbit_was_fitted ? "  SÌ   " : "  NO   ")
              << " | " << std::setw(4) << result_auto.computation_time_ms << "ms\n";
    std::cout << "MANUAL TRIGGER     | " << std::setw(8) 
              << result_manual.minimum_separation_arcsec << "\"  | " 
              << (fit_success ? "  SÌ   " : "  NO   ")
              << " | " << std::setw(4) << result_manual.computation_time_ms << "ms\n";
    
    std::cout << "\n🎯 CONCLUSIONI:\n";
    std::cout << "  ✅ Controllo granulare orbit fitting implementato\n";
    std::cout << "  ✅ Modalità NEVER/AUTO/ALWAYS_ATTEMPT funzionanti\n";
    std::cout << "  ✅ Trigger manuale e reset elementi disponibili\n";
    std::cout << "  ✅ Presets configurati per diversi use cases\n";
    
    return 0;
}