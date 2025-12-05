/**
 * Esempio Semplice: Propagatore con Osservazioni
 * ==============================================
 * 
 * Mostra l'uso base del propagatore con osservazioni per orbit fitting.
 */

#include <iostream>
#include <vector>
#include "ioccultcalc/propagation_strategy.h"

using namespace std;
using namespace ioccultcalc;

int main() {
    cout << "📡 ESEMPIO: Propagatore con Osservazioni" << endl;
    cout << "=======================================" << endl;
    
    // 1. Configurazione con fitting AUTO
    PropagationConfig config = propagation_presets::createBalanced();
    config.orbit_fitting_mode = PropagationConfig::FittingMode::AUTO;
    config.log_fitting_details = true;
    
    TwoPhaseStrategy propagator(config);
    
    // 2. Elementi orbitali di esempio (Eros)
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
    
    cout << "✅ Elementi orbitali caricati: " << elements.name << endl;
    
    // 3. Prova caricamento osservazioni (opzionale)
    cout << "\n📁 Caricamento osservazioni..." << endl;
    
    bool has_observations = false;
    
    // Opzione A: Da file locale
    if (propagator.loadObservationsFromFile("examples/433.rwo")) {
        cout << "✅ Caricate da file locale" << endl;
        has_observations = true;
    }
    // Opzione B: Download da AstDyS  
    else if (propagator.loadObservationsFromAstDyS(433)) {
        cout << "✅ Scaricate da AstDyS" << endl;
        has_observations = true;
    }
    // Opzione C: Nessuna osservazione (usa elementi originali)
    else {
        cout << "ℹ️  Nessuna osservazione - uso elementi originali" << endl;
    }
    
    // 4. Informazioni osservazioni (se disponibili)
    if (has_observations) {
        auto info = propagator.analyzeObservations();
        cout << "📊 " << info.total_observations << " osservazioni, "
             << "arco " << (int)info.observation_arc_days << " giorni" << endl;
        cout << "🎯 Qualità: " << (info.quality_acceptable ? "BUONA" : "SCARSA") 
             << " (RMS: " << info.mean_rms_arcsec << "\")" << endl;
        
        // Pulizia outliers
        int removed = propagator.cleanObservations(3.0);
        if (removed > 0) {
            cout << "🧹 Rimossi " << removed << " outliers" << endl;
        }
    }
    
    // 5. Test propagazione
    cout << "\n🎯 Test propagazione..." << endl;
    
    double target_mjd = 60010.0;   // Data target
    double ra_star = 150.5;        // RA stella [deg]
    double dec_star = -20.2;       // Dec stella [deg]
    
    cout << "📅 Epoca: MJD " << target_mjd << endl;
    cout << "⭐ Stella: RA " << ra_star << "°, Dec " << dec_star << "°" << endl;
    
    // FASE 1: Screening veloce
    cout << "\n🔍 FASE 1 - Screening veloce..." << endl;
    double separation = propagator.screenCandidate(target_mjd, ra_star, dec_star);
    cout << "📏 Separazione: " << separation << " arcsec" << endl;
    
    // FASE 2: Closest approach preciso (se candidato promettente)
    if (separation < 3600.0) {  // Se < 1 grado
        cout << "\n🎯 FASE 2 - Closest approach preciso..." << endl;
        
        try {
            auto result = propagator.findClosestApproach(target_mjd, ra_star, dec_star, 4.0);
            
            cout << "📍 Approccio più vicino:" << endl;
            cout << "   🕐 Tempo: MJD " << result.closest_time_mjd << endl;
            cout << "   📏 Distanza: " << result.minimum_separation_arcsec << " arcsec" << endl;
            cout << "   📐 Angolo: " << result.position_angle_deg << "°" << endl;
            
            // Informazioni orbit fitting
            if (result.orbit_was_fitted) {
                cout << "\n✨ Orbit fitting eseguito:" << endl;
                cout << "   🔄 Iterazioni: " << result.fitting_iterations << endl;
                cout << "   📈 Miglioramento: " << result.orbit_improvement_arcsec << " arcsec" << endl;
                cout << "   📊 RMS finale: " << result.final_rms_arcsec << " arcsec" << endl;
            } else {
                cout << "\nℹ️  Orbit fitting non eseguito (modalità: ";
                switch (propagator.getFittingMode()) {
                    case PropagationConfig::FittingMode::NEVER:
                        cout << "NEVER)";
                        break;
                    case PropagationConfig::FittingMode::AUTO:
                        cout << "AUTO, osservazioni insufficienti)";
                        break;
                    case PropagationConfig::FittingMode::ALWAYS_ATTEMPT:
                        cout << "SEMPRE, ma fallito)";
                        break;
                }
                cout << endl;
            }
            
            cout << "\n⚡ Performance: " << result.computation_time_ms << " ms" << endl;
            
        } catch (const exception& e) {
            cout << "❌ Errore FASE 2: " << e.what() << endl;
        }
        
    } else {
        cout << "❌ Candidato troppo distante per FASE 2" << endl;
    }
    
    // 6. Statistiche finali
    cout << "\n📊 Statistiche:" << endl;
    auto stats = propagator.getStats();
    cout << "   Screening: " << stats.candidates_screened << " candidati" << endl;
    cout << "   Promossi: " << stats.candidates_promoted << " alla FASE 2" << endl;
    cout << "   Orbit fitting: " << stats.orbits_fitted << " eseguiti" << endl;
    
    cout << "\n🎉 Esempio completato!" << endl;
    
    return 0;
}