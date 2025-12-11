/**
 * @file test_astdyn_integration.cpp
 * @brief Test completo integrazione ITALOccultLibrary/AstDyn
 * 
 * Verifica:
 * 1. AstDynPropagationHelper funziona correttamente
 * 2. OccultationPredictor usa AstDyn
 * 3. OccultationSearchAstDyn wrapper unificato
 * 4. Phase1 e Phase2 con AstDyn
 * 5. Conversioni tra formati
 * 
 * Compilazione:
 *   cd build && make test_astdyn_integration
 * 
 * Esecuzione:
 *   ./tests/test_astdyn_integration
 */

#include "ioccultcalc/astdyn_propagation_helper.h"
#include "ioccultcalc/astdyn_interface.h"
#include "ioccultcalc/occultation_predictor.h"
#include "ioccultcalc/occultation_search_astdyn.h"
#include "phase1_candidate_screening.h"
#include "phase2_occultation_geometry.h"
#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/time_utils.h"
#include "ioc_gaialib/unified_gaia_catalog.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <chrono>
 
 using namespace ioccultcalc;
 
 // Costanti
 constexpr double EPSILON = 1e-6;  // Tolleranza per test
 
 // Colori per output
 #define RESET   "\033[0m"
 #define GREEN   "\033[32m"
 #define RED     "\033[31m"
 #define YELLOW  "\033[33m"
 #define BLUE    "\033[34m"
 
 // Helper per stampa risultati
 void printTestHeader(const std::string& name) {
     std::cout << "\n" << BLUE << "═══════════════════════════════════════════════════════\n";
     std::cout << "TEST: " << name << "\n";
     std::cout << "═══════════════════════════════════════════════════════" << RESET << "\n\n";
 }
 
 void printTestResult(const std::string& name, bool passed) {
     if (passed) {
         std::cout << GREEN << "✓ PASS: " << name << RESET << "\n";
     } else {
         std::cout << RED << "✗ FAIL: " << name << RESET << "\n";
     }
 }
 
 // ============================================================================
 // TEST 1: AstDynPropagationHelper
 // ============================================================================
 
 bool testAstDynPropagationHelper() {
     printTestHeader("AstDynPropagationHelper");
     
     bool allPassed = true;
     
     try {
         // Test 1.1: Singleton funziona
         auto& helper1 = AstDynPropagationHelper::getInstance();
         auto& helper2 = AstDynPropagationHelper::getInstance();
         
         if (&helper1 != &helper2) {
             printTestResult("Singleton instance", false);
             allPassed = false;
         } else {
             printTestResult("Singleton instance", true);
         }
         
         // Test 1.2: Download elementi da AstDyS
         std::cout << "  Downloading elements for asteroid 433 (Eros)...\n";
         AstDySElements elements;
         try {
             elements = AstDySClient::downloadElements(433);
             printTestResult("Download elements from AstDyS", true);
             std::cout << "    Asteroide: " << elements.name << "\n";
             std::cout << "    a = " << std::fixed << std::setprecision(8) << elements.a << " AU\n";
             std::cout << "    e = " << elements.e << "\n";
             std::cout << "    epoch = " << elements.epoch_mjd << " MJD TDB\n";
         } catch (const std::exception& e) {
             printTestResult("Download elements from AstDyS", false);
             std::cerr << "    Errore: " << e.what() << "\n";
             allPassed = false;
             return allPassed;  // Non possiamo continuare senza elementi
         }
         
         // Test 1.3: Propagazione base
         std::cout << "  Testing propagation...\n";
         double target_mjd = elements.epoch_mjd + 10.0;  // +10 giorni
         AstDySElements propagated = helper1.propagate(elements, target_mjd);
         
         // Verifica che elementi siano validi
         if (propagated.a > 0 && propagated.e >= 0 && propagated.e < 1.0) {
             printTestResult("Basic propagation", true);
             std::cout << "    Propagated a = " << propagated.a << " AU\n";
         } else {
             printTestResult("Basic propagation", false);
             allPassed = false;
         }
         
         // Test 1.4: Calcolo RA/Dec
         std::cout << "  Testing RA/Dec calculation...\n";
         auto radec = helper1.getRADec(elements, elements.epoch_mjd);
         
         if (radec.first >= 0 && radec.first < 360 && 
             radec.second >= -90 && radec.second <= 90) {
             printTestResult("RA/Dec calculation", true);
             std::cout << "    RA = " << std::fixed << std::setprecision(6) 
                      << radec.first << "°\n";
             std::cout << "    Dec = " << radec.second << "°\n";
         } else {
             printTestResult("RA/Dec calculation", false);
             allPassed = false;
         }
         
         // Test 1.5: Propagazione range
         std::cout << "  Testing range propagation...\n";
         auto range = helper1.propagateRange(elements, elements.epoch_mjd, 
                                            elements.epoch_mjd + 5.0, 1.0);
         
         if (range.size() == 6) {  // 0, 1, 2, 3, 4, 5 giorni
             printTestResult("Range propagation", true);
             std::cout << "    Generated " << range.size() << " points\n";
         } else {
             printTestResult("Range propagation", false);
             std::cout << "    Expected 6 points, got " << range.size() << "\n";
             allPassed = false;
         }
         
         // Test 1.6: Configurazione tolleranza
         std::cout << "  Testing tolerance configuration...\n";
         helper1.setTolerance(1e-10);
         printTestResult("Tolerance configuration", true);
         
     } catch (const std::exception& e) {
         std::cerr << RED << "  Exception: " << e.what() << RESET << "\n";
         allPassed = false;
     }
     
     return allPassed;
 }
 
 // ============================================================================
 // TEST 2: OccultationPredictor con AstDyn
 // ============================================================================
 
 bool testOccultationPredictorAstDyn() {
     printTestHeader("OccultationPredictor with AstDyn");
     
     bool allPassed = true;
     
     try {
         OccultationPredictor predictor;
         
         // Test 2.1: Caricamento asteroide da AstDyS
         std::cout << "  Loading asteroid from AstDyS...\n";
         try {
             predictor.loadAsteroidFromAstDyS("433");  // Eros
             printTestResult("Load asteroid from AstDyS", true);
         } catch (const std::exception& e) {
             printTestResult("Load asteroid from AstDyS", false);
             std::cerr << "    Errore: " << e.what() << "\n";
             allPassed = false;
             return allPassed;
         }
         
         // Test 2.2: Configurazione parametri
         predictor.setAsteroidDiameter(16.84);  // km (Eros)
         predictor.setOrbitalUncertainty(50.0);  // km
         printTestResult("Set asteroid parameters", true);
         
         // Test 2.3: Ricerca occultazioni (intervallo breve per test)
         std::cout << "  Searching for occultations (short interval)...\n";
         JulianDate startJD = TimeUtils::isoToJD("2026-01-01");
         JulianDate endJD = TimeUtils::isoToJD("2026-01-02");  // Solo 1 giorno
         
         auto events = predictor.findOccultations(startJD, endJD, 18.0, 0.1, 0.001);
         
         printTestResult("Find occultations", true);
         std::cout << "    Found " << events.size() << " potential events\n";
         
         // Test 2.4: Predizione singola occultazione (se trovata)
         if (!events.empty()) {
             std::cout << "  Testing single event prediction...\n";
             auto event = predictor.predictOccultation(events[0].star, events[0].timeCA);
             
             if (event.probability >= 0 && event.probability <= 1.0) {
                 printTestResult("Single event prediction", true);
                 std::cout << "    Probability: " << std::fixed << std::setprecision(4)
                          << event.probability << "\n";
                 std::cout << "    Separation: " << event.closeApproachDistance << " arcsec\n";
             } else {
                 printTestResult("Single event prediction", false);
                 allPassed = false;
             }
         } else {
             std::cout << YELLOW << "  ⚠ No events found (this is OK for short interval)" << RESET << "\n";
         }
         
     } catch (const std::exception& e) {
         std::cerr << RED << "  Exception: " << e.what() << RESET << "\n";
         allPassed = false;
     }
     
     return allPassed;
 }
 
 // ============================================================================
 // TEST 3: OccultationSearchAstDyn
 // ============================================================================
 
 bool testOccultationSearchAstDyn() {
     printTestHeader("OccultationSearchAstDyn");
     
     bool allPassed = true;
     
     try {
         OccultationSearchAstDyn searcher;
         
         // Test 3.1: Caricamento asteroide
         std::cout << "  Loading asteroid...\n";
         if (searcher.loadAsteroid(433)) {  // Eros
             printTestResult("Load asteroid", true);
             std::cout << "    Name: " << searcher.getAsteroidName() << "\n";
             std::cout << "    Number: " << searcher.getAsteroidNumber() << "\n";
         } else {
             printTestResult("Load asteroid", false);
             allPassed = false;
             return allPassed;
         }
         
         // Test 3.2: Configurazione
         searcher.setAsteroidDiameter(16.84);
         searcher.setOrbitalUncertainty(50.0);
         printTestResult("Set configuration", true);
         
         // Test 3.3: Ricerca semplificata
         std::cout << "  Running simplified search...\n";
         double start_mjd = 61007.0;  // ~2026-01-01
         double end_mjd = 61008.0;    // ~2026-01-02
         
         auto results = searcher.search(start_mjd, end_mjd, 18.0);
         
         printTestResult("Simplified search", true);
         std::cout << "    Candidates found: " << results.num_candidates_found << "\n";
         std::cout << "    Stars in corridor: " << results.num_stars_in_corridor << "\n";
         std::cout << "    Events calculated: " << results.num_events_calculated << "\n";
         std::cout << "    Total time: " << std::fixed << std::setprecision(2)
                  << results.total_time_ms << " ms\n";
         
         // Test 3.4: Ricerca completa con configurazione
         std::cout << "  Running full search with config...\n";
         OccultationSearchConfig config;
         config.start_mjd_tdb = start_mjd;
         config.end_mjd_tdb = end_mjd;
         config.max_magnitude = 18.0;
         config.phase1_config = Phase1Config::conservative();
         config.phase2_config.time_window_minutes = 5.0;
         config.phase2_config.time_step_seconds = 1.0;
         
         auto results2 = searcher.search(config);
         
         printTestResult("Full search with config", true);
         std::cout << "    Events: " << results2.events.size() << "\n";
         
     } catch (const std::exception& e) {
         std::cerr << RED << "  Exception: " << e.what() << RESET << "\n";
         allPassed = false;
     }
     
     return allPassed;
 }
 
 // ============================================================================
 // TEST 4: Phase1 e Phase2 con AstDyn
 // ============================================================================
 
 bool testPhase1Phase2AstDyn() {
     printTestHeader("Phase1 and Phase2 with AstDyn");
     
     bool allPassed = true;
     
     try {
         // Test 4.1: Phase1
         std::cout << "  Testing Phase1CandidateScreening...\n";
         Phase1CandidateScreening phase1;
         
         if (phase1.loadAsteroidFromJSON(433)) {  // Eros
             printTestResult("Phase1 load asteroid", true);
         } else {
             printTestResult("Phase1 load asteroid", false);
             allPassed = false;
             return allPassed;
         }
         
         // Configura Phase1
         // phase1.setPropagatorTolerance(1e-9); // Non esiste un overload che accetta int
         // phase1.setUsePropagationHelper(false); // Non esiste un overload che accetta int
         
         Phase1Config phase1_config;
         phase1_config.start_mjd_tdb = 61007.0;
         phase1_config.end_mjd_tdb = 61008.0;
         phase1_config.max_magnitude = 18.0;
         phase1_config.corridor_width_deg = 0.0083;  // 30 arcsec
         
         auto phase1_results = phase1.screenCandidates(phase1_config);
         
         printTestResult("Phase1 screening", true);
         std::cout << "    Path points: " << phase1_results.num_path_points << "\n";
         std::cout << "    Stars in corridor: " << phase1_results.num_stars_in_corridor << "\n";
         std::cout << "    Candidates: " << phase1_results.num_candidates_filtered << "\n";
         std::cout << "    Phase1 time: " << std::fixed << std::setprecision(2)
                  << phase1_results.propagation_time_ms + 
                     phase1_results.corridor_query_time_ms +
                     phase1_results.closest_approach_calc_time_ms << " ms\n";
         
         // Test 4.2: Phase2 (solo se ci sono candidati)
         if (phase1_results.num_candidates_filtered > 0) {
             std::cout << "  Testing Phase2OccultationGeometry...\n";
             Phase2OccultationGeometry phase2;
             
             if (phase2.loadAsteroidFromJSON(433)) {
                 printTestResult("Phase2 load asteroid", true);
             } else {
                 printTestResult("Phase2 load asteroid", false);
                 allPassed = false;
                 return allPassed;
             }
             
             // Configura Phase2
             // phase2.setPropagatorTolerance(1e-12); // Non esiste un overload che accetta int
             // phase2.setPerturbations(true, false, true); // Non esiste un overload che accetta int
             
             Phase2Config phase2_config;
             phase2_config.time_window_minutes = 5.0;
             phase2_config.time_step_seconds = 1.0;
             phase2_config.use_planetary_perturbations = true;
             
             // Usa solo primi 3 candidati per test veloce
             std::vector<CandidateStar> test_candidates;
             for (size_t i = 0; i < std::min(3UL, phase1_results.candidates.size()); i++) {
                 test_candidates.push_back(phase1_results.candidates[i]);
             }
             
             auto phase2_results = phase2.calculateGeometry(test_candidates, phase2_config);
             
             printTestResult("Phase2 geometry calculation", true);
             std::cout << "    Events calculated: " << phase2_results.successful_calculations << "\n";
             std::cout << "    Events failed: " << phase2_results.failed_calculations << "\n";
             std::cout << "    Computation time: " << std::fixed << std::setprecision(2)
                      << phase2_results.total_computation_time_ms << " ms\n";
         } else {
             std::cout << YELLOW << "  ⚠ No candidates found, skipping Phase2 test" << RESET << "\n";
         }
         
     } catch (const std::exception& e) {
         std::cerr << RED << "  Exception: " << e.what() << RESET << "\n";
         allPassed = false;
     }
     
     return allPassed;
 }
 
 // ============================================================================
 // TEST 5: Conversioni tra formati
 // ============================================================================
 
 bool testFormatConversions() {
     printTestHeader("Format Conversions");
     
     bool allPassed = true;
     
     try {
         // Test 5.1: OrbitalElements -> AstDySElements
         std::cout << "  Testing OrbitalElements -> AstDySElements...\n";
         OrbitalElements orb;
         orb.designation = "Test";
         orb.a = 2.5;
         orb.e = 0.2;
         orb.i = 10.0 * DEG_TO_RAD;
         orb.Omega = 20.0 * DEG_TO_RAD;
         orb.omega = 30.0 * DEG_TO_RAD;
         orb.M = 40.0 * DEG_TO_RAD;
         orb.epoch.jd = 2450000.0;
         
         AstDySElements astdys = AstDynPropagationHelper::convertFromOrbital(orb);
         
         if (std::abs(astdys.a - 2.5) < EPSILON &&
             std::abs(astdys.e - 0.2) < EPSILON &&
             std::abs(astdys.i - 10.0) < EPSILON) {
             printTestResult("OrbitalElements -> AstDySElements", true);
         } else {
             printTestResult("OrbitalElements -> AstDySElements", false);
             allPassed = false;
         }
         
         // Test 5.2: EquinoctialElements -> AstDySElements
         std::cout << "  Testing EquinoctialElements -> AstDySElements...\n";
         EquinoctialElements eq;
         eq.designation = "Test";
         eq.a = 2.5;
         eq.h = 0.1;
         eq.k = 0.15;
         eq.p = 0.05;
         eq.q = 0.08;
         eq.lambda = 100.0 * DEG_TO_RAD;
         eq.epoch.jd = 2450000.0;
         
         AstDySElements astdys2 = AstDynPropagationHelper::convertFromEquinoctial(eq);
         
         if (std::abs(astdys2.a - 2.5) < EPSILON) {
             printTestResult("EquinoctialElements -> AstDySElements", true);
         } else {
             printTestResult("EquinoctialElements -> AstDySElements", false);
             allPassed = false;
         }
         
         // Test 5.3: Utility functions
         std::cout << "  Testing utility functions...\n";
         AstDySElements test_elem = astdys;
         test_elem.epoch_mjd = 60000.0;
         
         auto quick_prop = astdyn_utils::quickPropagate(test_elem, 60010.0);
         auto quick_radec = astdyn_utils::quickRADec(test_elem, 60000.0);
         
         if (quick_prop.a > 0 && quick_radec.first >= 0) {
             printTestResult("Utility functions", true);
         } else {
             printTestResult("Utility functions", false);
             allPassed = false;
         }
         
     } catch (const std::exception& e) {
         std::cerr << RED << "  Exception: " << e.what() << RESET << "\n";
         allPassed = false;
     }
     
     return allPassed;
 }
 
 // ============================================================================
 // MAIN
 // ============================================================================
 
int main() {
    std::cout << GREEN << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST INTEGRAZIONE ITALOccultLibrary/AstDyn              ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << RESET << "\n";
    
    // Inizializza UnifiedGaiaCatalog (richiesto per Phase1 e OccultationSearchAstDyn)
    std::cout << "Initializing UnifiedGaiaCatalog...\n";
    std::string gaia_config = R"({
        "catalog_type": "online_esa",
        "timeout_seconds": 60,
        "log_level": "warning"
    })";
    
    if (!ioc::gaia::UnifiedGaiaCatalog::initialize(gaia_config)) {
        std::cerr << YELLOW << "⚠ Warning: Failed to initialize UnifiedGaiaCatalog (online mode)\n";
        std::cerr << "  Some tests may fail. Continuing anyway...\n" << RESET << "\n";
    } else {
        std::cout << GREEN << "✓ UnifiedGaiaCatalog initialized (online mode)\n" << RESET << "\n";
    }
    
    int testsPassed = 0;
    int testsTotal = 5;
    
    auto start = std::chrono::high_resolution_clock::now();
     
     // Esegui tutti i test
     if (testAstDynPropagationHelper()) testsPassed++;
     if (testOccultationPredictorAstDyn()) testsPassed++;
     if (testOccultationSearchAstDyn()) testsPassed++;
     if (testPhase1Phase2AstDyn()) testsPassed++;
     if (testFormatConversions()) testsPassed++;
     
     auto end = std::chrono::high_resolution_clock::now();
     double elapsed = std::chrono::duration<double>(end - start).count();
     
     // Riepilogo
     std::cout << "\n";
     std::cout << BLUE << "═══════════════════════════════════════════════════════\n";
     std::cout << "RIEPILOGO TEST\n";
     std::cout << "═══════════════════════════════════════════════════════" << RESET << "\n\n";
     
     std::cout << "Test passati: " << GREEN << testsPassed << "/" << testsTotal << RESET << "\n";
     std::cout << "Tempo totale: " << std::fixed << std::setprecision(2) << elapsed << " secondi\n\n";
     
    // Cleanup
    ioc::gaia::UnifiedGaiaCatalog::shutdown();
    
    if (testsPassed == testsTotal) {
        std::cout << GREEN << "✓ TUTTI I TEST PASSATI!" << RESET << "\n\n";
        return 0;
    } else {
        std::cout << RED << "✗ ALCUNI TEST FALLITI" << RESET << "\n\n";
        return 1;
    }
}