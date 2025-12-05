/**
 * @file test_phase1_phase2_integration.cpp
 * @brief Test integrazione Phase 1 → Phase 2
 * 
 * Verifica il flusso completo:
 * 1. Phase1: carica elementi da JSON, trova stelle candidate
 * 2. Phase2: usa stessi elementi e analizza solo le candidate
 * 
 * @date 5 Dicembre 2025
 */

#include "phase1_candidate_screening.h"
#include "phase2_occultation_geometry.h"
#include "ioc_gaialib/unified_gaia_catalog.h"

#include <iostream>
#include <iomanip>
#include <chrono>

using namespace ioccultcalc;
using namespace ioc::gaia;

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║     TEST INTEGRAZIONE PHASE 1 → PHASE 2                     ║\n";
    std::cout << "║     Asteroide: 17030 Sierks                                  ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // INIZIALIZZAZIONE CATALOGO GAIA
    // ═══════════════════════════════════════════════════════════════
    
    std::cout << "=== INIZIALIZZAZIONE CATALOGO GAIA ===\n";
    
    const char* home = std::getenv("HOME");
    if (!home) {
        std::cerr << "Errore: HOME non definito\n";
        return 1;
    }
    
    std::string json_config = R"({
        "catalog_type": "multifile_v2",
        "multifile_directory": ")" + std::string(home) + "/.catalog/gaia_mag18_v2_multifile" + R"("
    })";
    
    try {
        UnifiedGaiaCatalog::initialize(json_config);
        std::cout << "✓ Catalogo Gaia inizializzato\n\n";
    } catch (const std::exception& e) {
        std::cerr << "✗ Errore inizializzazione catalogo: " << e.what() << "\n";
        return 1;
    }
    
    // ═══════════════════════════════════════════════════════════════
    // PHASE 1: SCREENING CANDIDATE
    // ═══════════════════════════════════════════════════════════════
    
    std::cout << "=== PHASE 1: SCREENING CANDIDATE ===\n";
    
    Phase1CandidateScreening phase1;
    
    // Carica elementi orbitali dal database JSON (NON da EQ1)
    if (!phase1.loadAsteroidFromJSON(17030)) {
        std::cerr << "✗ Errore caricamento elementi orbitali\n";
        return 1;
    }
    
    // Configura Phase 1 per 25 Nov - 2 Dic 2025 (dal preset)
    Phase1Config config1;
    
    // Intervallo temporale dal preset: JD 2460638.5 - 2460645.5
    // MJD = JD - 2400000.5
    config1.start_mjd_tdb = 2460638.5 - 2400000.5;  // 25 Nov 2025 00:00 UTC
    config1.end_mjd_tdb = 2460645.5 - 2400000.5;    // 2 Dic 2025 00:00 UTC (7 giorni)
    
    // Parametri screening dal preset
    config1.corridor_width_deg = 0.0083;  // 30 arcsec
    config1.max_magnitude = 16.0;         // Dal preset
    config1.closest_approach_threshold_arcsec = 15.0;  // Solo CA < 15"
    
    auto start = std::chrono::high_resolution_clock::now();
    
    try {
        Phase1Results results1 = phase1.screenCandidates(config1);
        
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
        
        std::cout << "\n--- Risultati Phase 1 ---\n";
        std::cout << "Path points:        " << results1.num_path_points << "\n";
        std::cout << "Stelle nel corridor: " << results1.num_stars_in_corridor << "\n";
        std::cout << "Candidate filtrate: " << results1.num_candidates_filtered << "\n";
        std::cout << "Tempo totale:       " << std::fixed << std::setprecision(1) << elapsed_ms << " ms\n";
        
        // Mostra candidate trovate
        std::cout << "\n--- Candidate per Phase 2 ---\n";
        for (size_t i = 0; i < results1.candidates.size(); i++) {
            const auto& c = results1.candidates[i];
            std::cout << i+1 << ") Gaia DR3 " << c.source_id << "\n";
            std::cout << "   RA=" << std::fixed << std::setprecision(6) << c.ra_deg << "° ";
            std::cout << "Dec=" << c.dec_deg << "°\n";
            std::cout << "   Gmag=" << std::setprecision(2) << c.phot_g_mean_mag;
            std::cout << "  CA=" << std::setprecision(3) << c.closest_approach_arcsec << "\"\n";
            std::cout << "   MJD CA=" << std::setprecision(5) << c.closest_approach_mjd << "\n";
            std::cout << "\n";
        }
        
        if (results1.candidates.empty()) {
            std::cout << "Nessuna candidata trovata. Fine test.\n";
            return 0;
        }
        
        // ═══════════════════════════════════════════════════════════════
        // PHASE 2: GEOMETRIA PRECISA (solo candidate da Phase 1)
        // ═══════════════════════════════════════════════════════════════
        
        std::cout << "\n=== PHASE 2: GEOMETRIA PRECISA ===\n";
        std::cout << "Analisi delle " << results1.candidates.size() << " candidate trovate...\n\n";
        
        Phase2OccultationGeometry phase2;
        
        // USA GLI STESSI ELEMENTI DELLA PHASE 1!
        phase2.setOrbitalElements(phase1.getOrbitalElements());
        std::cout << "✓ Elementi orbitali passati da Phase 1 a Phase 2\n";
        
        // Configura Phase 2
        Phase2Config config2;
        config2.time_window_minutes = 5.0;    // ±5 min attorno CA
        config2.time_step_seconds = 1.0;       // 1 sec risoluzione
        config2.use_planetary_perturbations = true;
        config2.use_relativistic_effects = true;
        config2.refine_orbit_from_observations = false;  // NO fitting (disabilitato)
        
        // Aggiungi osservatori italiani
        phase2.addObserverSite("Roma", 41.9028, 12.4964, 21);
        phase2.addObserverSite("Milano", 45.4642, 9.1900, 120);
        phase2.addObserverSite("Bologna", 44.4949, 11.3426, 54);
        
        // Calcola geometria per le candidate
        start = std::chrono::high_resolution_clock::now();
        
        Phase2Results results2 = phase2.calculateGeometry(results1.candidates, config2);
        
        end = std::chrono::high_resolution_clock::now();
        elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
        
        std::cout << "\n--- Risultati Phase 2 ---\n";
        std::cout << "Eventi calcolati:   " << results2.successful_calculations << "\n";
        std::cout << "Errori:             " << results2.failed_calculations << "\n";
        std::cout << "Tempo totale:       " << std::fixed << std::setprecision(1) << elapsed_ms << " ms\n";
        
        // Mostra dettagli eventi
        for (size_t i = 0; i < results2.events.size(); i++) {
            const auto& evt = results2.events[i];
            std::cout << "\n═══ EVENTO " << i+1 << " ═══\n";
            std::cout << "Stella:      Gaia DR3 " << evt.star_source_id << "\n";
            std::cout << "Magnitudine: " << std::setprecision(2) << evt.star_magnitude << " G\n";
            std::cout << "Closest approach:\n";
            std::cout << "  MJD UTC:   " << std::setprecision(6) << evt.time_ca_mjd_utc << "\n";
            std::cout << "  Distanza:  " << std::setprecision(1) << evt.closest_approach_mas << " mas\n";
            std::cout << "  Durata max: " << std::setprecision(2) << evt.max_duration_sec << " sec\n";
            std::cout << "  PA:        " << std::setprecision(1) << evt.position_angle_deg << "°\n";
            
            if (!evt.observer_predictions.empty()) {
                std::cout << "\nPredizioni osservatori:\n";
                for (const auto& obs : evt.observer_predictions) {
                    std::cout << "  " << obs.site_name << ": ";
                    if (obs.is_in_shadow_path) {
                        std::cout << "NEL PATH! ";
                    }
                    std::cout << "Alt=" << std::setprecision(0) << obs.target_altitude_deg << "°";
                    std::cout << " Az=" << obs.target_azimuth_deg << "°\n";
                }
            }
        }
        
        std::cout << "\n✓ Test completato con successo!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "✗ Errore: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
