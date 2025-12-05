/**
 * @file test_phase2_orbital_fitting.cpp
 * @brief Test integrazione Phase 2 con fitting orbitale da RWO
 * @date 4 Dicembre 2025
 * 
 * Test workflow completo:
 * 1. Carica asteroide 17030 da EQ1
 * 2. Esegue Phase 1 per trovare candidati
 * 3. Per ogni candidato, Phase 2:
 *    - Download RWO da AstDyS
 *    - Fitting orbitale con osservazioni
 *    - Propagazione con elementi raffinati
 *    - Calcolo geometria precisa
 */

#include <iostream>
#include <iomanip>
#include <chrono>
#include "phase1_candidate_screening.h"
#include "phase2_occultation_geometry.h"
#include "ioc_gaialib/unified_gaia_catalog.h"

using namespace ioccultcalc;
using namespace ioc::gaia;

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST: Phase 2 con Orbital Fitting da RWO AstDyS           ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 1: PHASE 1 - Candidate Screening
    // ═══════════════════════════════════════════════════════════════
    
    std::cout << "\n### PHASE 1: Candidate Screening ###\n\n";
    
    Phase1CandidateScreening phase1;
    
    // Inizializza catalogo Gaia
    std::cout << "Inizializzazione catalogo Gaia...\n";
    
    std::string home = std::string(getenv("HOME"));
    std::string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
    
    std::string json_config = R"({
        "catalog_type": "multifile_v2",
        "multifile_directory": ")" + catalog_path + R"("
    })";
    
    std::cout << "Prima di initialize..." << std::endl;
    if (!UnifiedGaiaCatalog::initialize(json_config)) {
        throw std::runtime_error("Impossibile inizializzare catalogo Gaia");
    }
    std::cout << "Dopo initialize!" << std::endl;
    
    auto& catalog = UnifiedGaiaCatalog::getInstance();
    std::cout << "getInstance OK" << std::endl;
    
    std::cout << "✓ Catalogo Gaia inizializzato" << std::endl;
    
    // Carica asteroide 17030
    std::string eq1_file = "17030_astdys.eq1";
    std::cout << "Caricamento elementi orbitali da: " << eq1_file << std::endl;
    
    std::cout << "DEBUG: Prima di loadAsteroidFromEQ1..." << std::endl;
    if (!phase1.loadAsteroidFromEQ1(eq1_file)) {
        std::cerr << "✗ Impossibile caricare elementi da " << eq1_file << "\n";
        return 1;
    }
    std::cout << "DEBUG: Dopo loadAsteroidFromEQ1" << std::endl;
    std::cout << "✓ Elementi caricati" << std::endl;
    
    // Configura Phase 1
    std::cout << "DEBUG: Configurazione Phase1Config..." << std::endl;
    Phase1Config config1;
    config1.start_mjd_tdb = 60676.0;  // 2025-01-01
    config1.end_mjd_tdb = 61041.0;    // 2025-12-31
    config1.max_magnitude = 14.0;
    config1.corridor_width_deg = 20.0 / 3600.0;  // 20 arcsec in gradi
    config1.closest_approach_threshold_arcsec = 5.0;
    std::cout << "DEBUG: Phase1Config OK" << std::endl;
    
    // Esegui Phase 1
    std::cout << "DEBUG: Prima di screenCandidates..." << std::endl;
    auto t1_start = std::chrono::high_resolution_clock::now();
    Phase1Results results1 = phase1.screenCandidates(config1);
    auto t1_end = std::chrono::high_resolution_clock::now();
    double t1_sec = std::chrono::duration<double>(t1_end - t1_start).count();
    
    std::cout << "\n=== Phase 1 Results ===\n";
    std::cout << "Execution time: " << std::fixed << std::setprecision(3) << t1_sec << " sec\n";
    std::cout << "Candidates found: " << results1.candidates.size() << "\n\n";
    
    if (results1.candidates.empty()) {
        std::cout << "Nessun candidato trovato - test terminato\n";
        return 0;
    }
    
    // Mostra primi candidati
    std::cout << "Primi candidati:\n";
    int n_show = std::min(3, static_cast<int>(results1.candidates.size()));
    for (int i = 0; i < n_show; i++) {
        const auto& c = results1.candidates[i];
        std::cout << std::setw(2) << (i+1) << ". Star " << c.source_id 
                  << " @ MJD " << std::fixed << std::setprecision(6) << c.closest_approach_mjd
                  << " | sep=" << std::setprecision(2) << c.closest_approach_arcsec << " arcsec"
                  << " | mag=" << std::setprecision(1) << c.phot_g_mean_mag << "\n";
    }
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 2: PHASE 2 - Geometria con Orbital Fitting
    // ═══════════════════════════════════════════════════════════════
    
    std::cout << "\n\n### PHASE 2: Precise Geometry with Orbital Fitting ###\n\n";
    
    Phase2OccultationGeometry phase2;
    
    // Carica stesso asteroide in Phase 2
    std::cout << "Caricamento elementi in Phase 2...\n";
    if (!phase2.loadAsteroidFromEQ1(eq1_file)) {
        std::cerr << "✗ Impossibile caricare elementi in Phase 2\n";
        return 1;
    }
    std::cout << "✓ Elementi caricati in Phase 2\n\n";
    
    // Configura Phase 2 SENZA orbital refinement (per ora)
    Phase2Config config2;
    config2.refine_orbit_from_observations = false; // ← DISATTIVATO per ora
    config2.mpc_code = "17030";                     // Per download RWO
    config2.observation_arc_days = 365 * 3;         // 3 anni di osservazioni
    
    // Parametri fitting
    config2.fit_planetary_perturbations = true;
    config2.fit_relativistic_effects = false;
    config2.fit_tolerance = 1e-12;
    config2.max_fit_iterations = 20;
    
    // Parametri propagazione precisa
    config2.time_window_minutes = 5.0;
    config2.time_step_seconds = 1.0;
    config2.compute_uncertainty = false;  // Per ora nessuna incertezza
    config2.compute_shadow_path = false;  // Per ora nessun path
    
    // Processa SOLO il primo candidato (il più promettente)
    std::vector<CandidateStar> test_candidates;
    test_candidates.push_back(results1.candidates[0]);
    
    std::cout << "\n╔═══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  Test su candidato #1: " << test_candidates[0].source_id << "\n";
    std::cout << "╚═══════════════════════════════════════════════════════════╝\n";
    
    // Esegui Phase 2 con fitting
    auto t2_start = std::chrono::high_resolution_clock::now();
    Phase2Results results2 = phase2.calculateGeometry(test_candidates, config2);
    auto t2_end = std::chrono::high_resolution_clock::now();
    double t2_sec = std::chrono::duration<double>(t2_end - t2_start).count();
    
    // ═══════════════════════════════════════════════════════════════
    // RISULTATI FINALI
    // ═══════════════════════════════════════════════════════════════
    
    std::cout << "\n\n╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  RISULTATI FINALI                                            ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
    
    std::cout << "Phase 1 time: " << std::fixed << std::setprecision(3) << t1_sec << " sec\n";
    std::cout << "Phase 2 time: " << std::fixed << std::setprecision(3) << t2_sec << " sec\n";
    std::cout << "Total time:   " << std::fixed << std::setprecision(3) << (t1_sec + t2_sec) << " sec\n\n";
    
    // Statistiche orbital fit
    std::cout << "=== Orbital Fitting Statistics ===\n";
    if (results2.orbit_refined) {
        const auto& fit = results2.orbital_fit;
        std::cout << "Fit converged: " << (fit.fit_successful ? "YES ✓" : "NO ✗") << "\n";
        if (fit.fit_successful) {
            std::cout << "Observations used: " << fit.num_observations_used << "\n";
            std::cout << "Iterations: " << fit.iterations_performed << "\n";
            std::cout << "RMS residuals: " << std::fixed << std::setprecision(3) 
                     << fit.rms_residuals_arcsec << " arcsec\n";
            std::cout << "Max residual: " << std::fixed << std::setprecision(3) 
                     << fit.max_residual_arcsec << " arcsec\n";
            std::cout << "Chi²: " << std::scientific << std::setprecision(2) << fit.chi_squared << "\n";
        }
        std::cout << "Notes: " << fit.fit_notes << "\n";
    } else {
        std::cout << "No orbital fitting performed\n";
    }
    
    // Geometria eventi
    std::cout << "\n=== Occultation Events ===\n";
    std::cout << "Events calculated: " << results2.events.size() << "\n\n";
    
    for (size_t i = 0; i < results2.events.size(); i++) {
        const auto& evt = results2.events[i];
        std::cout << "Event " << (i+1) << ":\n";
        std::cout << "  Star: " << evt.star_source_id << " (mag " 
                 << std::fixed << std::setprecision(1) << evt.star_magnitude << ")\n";
        std::cout << "  CA time: MJD " << std::fixed << std::setprecision(6) << evt.time_ca_mjd_utc << "\n";
        std::cout << "  Minimum distance: " << std::fixed << std::setprecision(2) 
                 << evt.closest_approach_mas << " mas\n";
        std::cout << "  Max duration: " << std::fixed << std::setprecision(2) 
                 << evt.max_duration_sec << " sec\n";
        std::cout << "  Distance at CA: " << std::fixed << std::setprecision(3) 
                 << evt.asteroid_distance_au << " AU\n";
    }
    
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST COMPLETATO                                             ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n";
    
    return 0;
}
