/**
 * Esempio Avanzato: Propagatore con Osservazioni Opzionali
 * ========================================================
 * 
 * Dimostra come usare il propagatore con controllo granulare delle osservazioni
 * per l'orbit fitting, inclusi caricamento da file e validazione qualità.
 */

#include <iostream>
#include <iomanip>
#include "ioccultcalc/propagation_strategy.h"

using namespace std;
using namespace ioccultcalc;

void testObservationManagement() {
    cout << "🧪 TEST GESTIONE OSSERVAZIONI" << endl;
    cout << "=============================" << endl;
    
    // Configurazione con orbit fitting abilitato
    PropagationConfig config = propagation_presets::createPrecision();
    config.log_fitting_details = true;
    config.verbose_timing = true;
    
    TwoPhaseStrategy propagator(config);
    
    cout << "\n1️⃣ CARICAMENTO OSSERVAZIONI" << endl;
    cout << "   Modalità fitting: " << 
        (config.orbit_fitting_mode == PropagationConfig::FittingMode::ALWAYS_ATTEMPT ? "SEMPRE" :
         config.orbit_fitting_mode == PropagationConfig::FittingMode::AUTO ? "AUTO" : "MAI") << endl;
    
    // Esempio 1: Caricamento da file
    cout << "\n📁 Test caricamento da file RWO..." << endl;
    bool loaded_from_file = propagator.loadObservationsFromFile("examples/433.rwo");
    
    if (!loaded_from_file) {
        cout << "⚠️  File RWO non trovato, provo download da AstDyS..." << endl;
        
        // Esempio 2: Download da AstDyS
        bool loaded_from_astdys = propagator.loadObservationsFromAstDyS(433);  // Eros
        
        if (!loaded_from_astdys) {
            cout << "❌ Download fallito, creo osservazioni mock..." << endl;
            
            // Esempio 3: Osservazioni create manualmente
            vector<RWOObservation> mock_observations;
            
            for (int i = 0; i < 10; ++i) {
                RWOObservation obs;
                obs.designation = "433";
                obs.mjd_utc = 60000.0 + i * 10.0;  // Ogni 10 giorni
                obs.ra_deg = 150.0 + i * 0.1;
                obs.dec_deg = -20.0 + i * 0.05;
                obs.rms_ra_arcsec = 0.5;
                obs.rms_dec_arcsec = 0.4;
                mock_observations.push_back(obs);
            }
            
            propagator.setObservations(mock_observations);
            cout << "✅ Create " << mock_observations.size() << " osservazioni mock" << endl;
        }
    }
    
    cout << "\n2️⃣ ANALISI QUALITÀ OSSERVAZIONI" << endl;
    
    if (propagator.hasObservations()) {
        auto info = propagator.analyzeObservations();
        
        cout << "📊 Informazioni osservazioni:" << endl;
        cout << "   Totale: " << info.total_observations << " osservazioni" << endl;
        cout << "   Arco: " << fixed << setprecision(1) << info.observation_arc_days << " giorni" << endl;
        cout << "   Periodo: MJD " << setprecision(1) << info.first_observation_mjd 
             << " → " << info.last_observation_mjd << endl;
        cout << "   RMS medio: " << setprecision(2) << info.mean_rms_arcsec << " arcsec" << endl;
        cout << "   Qualità: " << (info.quality_acceptable ? "✅ BUONA" : "⚠️  SCARSA") << endl;
        
        // Pulizia outliers
        cout << "\n🧹 Pulizia outliers..." << endl;
        int removed = propagator.cleanObservations(3.0);  // Soglia 3 arcsec
        cout << "   Rimosse: " << removed << " osservazioni con RMS > 3\"" << endl;
        cout << "   Rimaste: " << propagator.getObservationCount() << " osservazioni" << endl;
        
    } else {
        cout << "❌ Nessuna osservazione disponibile" << endl;
    }
    
    cout << "\n3️⃣ TEST MODALITÀ FITTING" << endl;
    
    // Elementi di esempio per Eros
    AstDySElements elements;
    elements.name = "433 Eros";
    elements.number = 433;
    elements.epoch_mjd = 60000.0;
    elements.a = 1.458;     // AU
    elements.e = 0.223;
    elements.i = 10.83;     // deg
    elements.omega = 178.8; // deg
    elements.Omega = 304.4; // deg
    elements.M = 320.1;     // deg
    
    propagator.setElements(elements);
    
    // Test diverse modalità di fitting
    vector<PropagationConfig::FittingMode> test_modes = {
        PropagationConfig::FittingMode::NEVER,
        PropagationConfig::FittingMode::AUTO,
        PropagationConfig::FittingMode::ALWAYS_ATTEMPT
    };
    
    for (auto mode : test_modes) {
        cout << "\n🔧 Test modalità: ";
        switch (mode) {
            case PropagationConfig::FittingMode::NEVER:
                cout << "NEVER (no fitting)" << endl;
                break;
            case PropagationConfig::FittingMode::AUTO:
                cout << "AUTO (fitting se possibile)" << endl;
                break;
            case PropagationConfig::FittingMode::ALWAYS_ATTEMPT:
                cout << "ALWAYS_ATTEMPT (fitting obbligatorio)" << endl;
                break;
        }
        
        propagator.setFittingMode(mode);
        propagator.resetFittedElements();  // Forza re-fit
        
        // Test closest approach con fitting
        double target_mjd = 60010.0;
        double ra_star = 150.5;  // deg
        double dec_star = -20.2; // deg
        
        try {
            auto result = propagator.findClosestApproach(target_mjd, ra_star, dec_star, 4.0);
            
            cout << "   📍 Closest approach: " << fixed << setprecision(3) 
                 << result.minimum_separation_arcsec << " arcsec" << endl;
            cout << "   🕐 Tempo: MJD " << setprecision(6) << result.closest_time_mjd << endl;
            
            if (result.orbit_was_fitted) {
                cout << "   ✅ Orbit fitting: " << result.fitting_iterations << " iterazioni" << endl;
                cout << "   📈 Miglioramento: " << setprecision(2) 
                     << result.orbit_improvement_arcsec << " arcsec" << endl;
                cout << "   📊 RMS finale: " << result.final_rms_arcsec << " arcsec" << endl;
            } else {
                cout << "   ➖ Orbit fitting: NON ESEGUITO" << endl;
            }
            
        } catch (const exception& e) {
            cout << "   ❌ Errore: " << e.what() << endl;
        }
    }
    
    cout << "\n4️⃣ CONTROLLO MANUALE FITTING" << endl;
    
    // Reset e trigger manuale
    propagator.setFittingMode(PropagationConfig::FittingMode::NEVER);
    propagator.resetFittedElements();
    
    cout << "📍 Stato attuale: " << (propagator.isUsingFittedElements() ? 
                                   "Usando elementi fitted" : "Usando elementi originali") << endl;
    
    cout << "\n🔧 Trigger manuale orbit fitting..." << endl;
    bool fitting_success = propagator.triggerOrbitFitting();
    
    if (fitting_success) {
        cout << "✅ Fitting completato con successo" << endl;
        cout << "📍 Nuovo stato: " << (propagator.isUsingFittedElements() ? 
                                     "Usando elementi fitted" : "Usando elementi originali") << endl;
    } else {
        cout << "❌ Fitting fallito o non eseguibile" << endl;
    }
    
    cout << "\n5️⃣ STATISTICHE PERFORMANCE" << endl;
    
    auto stats = propagator.getStats();
    
    cout << "📊 Chiamate FASE 1: " << stats.phase1_calls << endl;
    cout << "📊 Chiamate FASE 2: " << stats.phase2_calls << endl;
    cout << "📊 Candidati testati: " << stats.candidates_screened << endl;
    cout << "📊 Candidati promossi: " << stats.candidates_promoted << endl;
    cout << "📊 Orbit fitting eseguiti: " << stats.orbits_fitted << endl;
    
    if (stats.orbits_fitted > 0) {
        cout << "📊 Miglioramento medio: " << fixed << setprecision(2) 
             << (stats.total_orbit_improvement_arcsec / stats.orbits_fitted) 
             << " arcsec" << endl;
    }
    
    if (stats.phase1_calls > 0) {
        cout << "⚡ Tempo medio FASE 1: " << setprecision(3) 
             << (stats.phase1_total_time_ms / stats.phase1_calls) << " ms" << endl;
    }
    
    if (stats.phase2_calls > 0) {
        cout << "⚡ Tempo medio FASE 2: " << setprecision(1) 
             << (stats.phase2_total_time_ms / stats.phase2_calls) << " ms" << endl;
    }
    
    cout << "\n🎉 TEST COMPLETATO!" << endl;
}

int main() {
    cout << "🚀 PROPAGATORE CON OSSERVAZIONI OPZIONALI" << endl;
    cout << "=========================================" << endl;
    cout << endl;
    
    try {
        testObservationManagement();
        
        cout << "\n✅ Tutti i test completati con successo!" << endl;
        cout << "\n📖 RIASSUNTO FUNZIONALITÀ:" << endl;
        cout << "   • Caricamento osservazioni da file .rwo" << endl;
        cout << "   • Download automatico da AstDyS" << endl;
        cout << "   • Analisi qualità osservazioni" << endl;
        cout << "   • Pulizia automatica outliers" << endl;
        cout << "   • Controllo granulare fitting (NEVER/AUTO/ALWAYS)" << endl;
        cout << "   • Trigger manuale orbit fitting" << endl;
        cout << "   • Statistiche performance dettagliate" << endl;
        
    } catch (const exception& e) {
        cout << "❌ Errore durante test: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}