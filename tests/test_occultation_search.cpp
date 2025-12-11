/**
 * @file test_occultation_search.cpp
 * @brief Test completo per OccultationSearchAstDyn
 * 
 * Verifica:
 * 1. Caricamento asteroidi da AstDyS
 * 2. Ricerca semplificata
 * 3. Ricerca completa con configurazione
 * 4. Validazione risultati
 * 
 * Compilazione:
 *   cd build && make test_occultation_search
 * 
 * Esecuzione:
 *   ./tests/test_occultation_search
 */

#include "ioccultcalc/occultation_search_astdyn.h"
#include "ioc_gaialib/unified_gaia_catalog.h"
#include "ioccultcalc/time_utils.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>

using namespace ioccultcalc;

// Colori per output
#define RESET   "\033[0m"
#define GREEN   "\033[32m"
#define RED     "\033[31m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define CYAN    "\033[36m"

void printHeader(const std::string& title) {
    std::cout << "\n" << BLUE << "═══════════════════════════════════════════════════════\n";
    std::cout << title << "\n";
    std::cout << "═══════════════════════════════════════════════════════" << RESET << "\n\n";
}

void printTest(const std::string& name, bool passed) {
    if (passed) {
        std::cout << GREEN << "✓ " << name << RESET << "\n";
    } else {
        std::cout << RED << "✗ " << name << RESET << "\n";
    }
}

void printInfo(const std::string& label, const std::string& value) {
    std::cout << "  " << CYAN << label << ": " << RESET << value << "\n";
}

void printValue(const std::string& label, double value, int precision = 2) {
    std::cout << "  " << CYAN << label << ": " << RESET 
              << std::fixed << std::setprecision(precision) << value << "\n";
}

// ============================================================================
// TEST 1: Caricamento asteroide
// ============================================================================

bool testLoadAsteroid() {
    printHeader("TEST 1: Caricamento Asteroide");
    
    bool allPassed = true;
    
    try {
        OccultationSearchAstDyn searcher;
        
        // Test 1.1: Carica Eros (433)
        std::cout << "Caricamento asteroide 433 (Eros)...\n";
        if (searcher.loadAsteroid(433)) {
            printTest("Load asteroid by number", true);
            printInfo("Name", searcher.getAsteroidName());
            printValue("Number", searcher.getAsteroidNumber());
            
            // Verifica elementi orbitali
            auto elements = searcher.getAsteroidElements();
            printValue("Semi-major axis (AU)", elements.a, 6);
            printValue("Eccentricity", elements.e, 6);
            printValue("Inclination (deg)", elements.i, 4);
            printValue("Epoch (MJD TDB)", elements.epoch_mjd, 2);
        } else {
            printTest("Load asteroid by number", false);
            allPassed = false;
        }
        
        // Test 1.2: Carica per designazione
        std::cout << "\nCaricamento asteroide per designazione '433'...\n";
        OccultationSearchAstDyn searcher2;
        if (searcher2.loadAsteroid("433")) {
            printTest("Load asteroid by designation", true);
        } else {
            printTest("Load asteroid by designation", false);
            allPassed = false;
        }
        
        // Test 1.3: Verifica hasAsteroid
        std::cout << "\nVerifica stato caricamento...\n";
        if (searcher.hasAsteroid()) {
            printTest("hasAsteroid() returns true", true);
        } else {
            printTest("hasAsteroid() returns true", false);
            allPassed = false;
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Exception: " << e.what() << RESET << "\n";
        allPassed = false;
    }
    
    return allPassed;
}

// ============================================================================
// TEST 2: Ricerca semplificata
// ============================================================================

bool testSimplifiedSearch() {
    printHeader("TEST 2: Ricerca Semplificata");
    
    bool allPassed = true;
    
    try {
        OccultationSearchAstDyn searcher;
        
        // Carica asteroide
        if (!searcher.loadAsteroid(433)) {
            std::cerr << RED << "Failed to load asteroid" << RESET << "\n";
            return false;
        }
        
        // Configura asteroide
        searcher.setAsteroidDiameter(16.84);  // Eros diameter in km
        searcher.setOrbitalUncertainty(50.0);  // 50 km uncertainty
        
        // Intervallo temporale: ~1 giorno a gennaio 2026
        double start_mjd = 61007.0;  // 2026-01-01
        double end_mjd = 61008.0;    // 2026-01-02
        double max_magnitude = 18.0;
        
        std::cout << "Esecuzione ricerca semplificata...\n";
        std::cout << "  Periodo: MJD " << start_mjd << " - " << end_mjd << "\n";
        std::cout << "  Magnitudine limite: " << max_magnitude << "\n";
        
        auto t_start = std::chrono::high_resolution_clock::now();
        auto results = searcher.search(start_mjd, end_mjd, max_magnitude);
        auto t_end = std::chrono::high_resolution_clock::now();
        
        double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        
        printTest("Simplified search completed", true);
        printValue("Total time (ms)", elapsed_ms);
        printValue("Candidates found", results.num_candidates_found);
        printValue("Stars in corridor", results.num_stars_in_corridor);
        printValue("Events calculated", results.num_events_calculated);
        printValue("Phase1 time (ms)", results.phase1_time_ms);
        printValue("Phase2 time (ms)", results.phase2_time_ms);
        
        // Verifica risultati
        if (results.num_candidates_found >= 0) {
            printTest("Valid candidate count", true);
        } else {
            printTest("Valid candidate count", false);
            allPassed = false;
        }
        
        if (results.events.size() == static_cast<size_t>(results.num_events_calculated)) {
            printTest("Event count consistency", true);
        } else {
            printTest("Event count consistency", false);
            allPassed = false;
        }
        
        // Mostra eventi trovati
        if (results.events.size() > 0) {
            std::cout << "\n" << GREEN << "Eventi trovati:" << RESET << "\n";
            for (size_t i = 0; i < std::min(results.events.size(), size_t(5)); ++i) {
                const auto& event = results.events[i];
                std::cout << "  Evento " << (i+1) << ":\n";
                printValue("  Time CA (JD)", event.timeCA.jd, 5);
                printInfo("  Star ID", event.eventId);
                printValue("  Separation (arcsec)", event.closeApproachDistance, 2);
                printValue("  Probability", event.probability, 4);
                printValue("  Max duration (sec)", event.maxDuration, 2);
            }
            if (results.events.size() > 5) {
                std::cout << "  ... e altri " << (results.events.size() - 5) << " eventi\n";
            }
        } else {
            std::cout << YELLOW << "  Nessun evento trovato (normale per intervallo breve)" << RESET << "\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Exception: " << e.what() << RESET << "\n";
        allPassed = false;
    }
    
    return allPassed;
}

// ============================================================================
// TEST 3: Ricerca completa con configurazione
// ============================================================================

bool testFullSearch() {
    printHeader("TEST 3: Ricerca Completa con Configurazione");
    
    bool allPassed = true;
    
    try {
        OccultationSearchAstDyn searcher;
        
        // Carica asteroide
        if (!searcher.loadAsteroid(433)) {
            std::cerr << RED << "Failed to load asteroid" << RESET << "\n";
            return false;
        }
        
        // Configura asteroide
        searcher.setAsteroidDiameter(16.84);
        searcher.setOrbitalUncertainty(50.0);
        
        // Crea configurazione completa
        OccultationSearchConfig config;
        config.start_mjd_tdb = 61007.0;
        config.end_mjd_tdb = 61010.0;  // 3 giorni
        config.max_magnitude = 18.0;
        config.asteroid_diameter_km = 16.84;
        config.orbital_uncertainty_km = 50.0;
        config.min_probability = 0.01;
        
        // Configurazione Phase1
        config.phase1_config = Phase1Config::conservative();
        config.phase1_config.start_mjd_tdb = config.start_mjd_tdb;
        config.phase1_config.end_mjd_tdb = config.end_mjd_tdb;
        config.phase1_config.max_magnitude = config.max_magnitude;
        config.phase1_config.corridor_width_deg = 0.0083;  // 30 arcsec
        
        // Configurazione Phase2
        config.phase2_config.time_window_minutes = 10.0;
        config.phase2_config.time_step_seconds = 1.0;
        
        std::cout << "Esecuzione ricerca completa...\n";
        std::cout << "  Periodo: MJD " << config.start_mjd_tdb << " - " << config.end_mjd_tdb << "\n";
        std::cout << "  Magnitudine limite: " << config.max_magnitude << "\n";
        std::cout << "  Diametro asteroide: " << config.asteroid_diameter_km << " km\n";
        std::cout << "  Incertezza orbitale: " << config.orbital_uncertainty_km << " km\n";
        
        auto t_start = std::chrono::high_resolution_clock::now();
        auto results = searcher.search(config);
        auto t_end = std::chrono::high_resolution_clock::now();
        
        double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        
        printTest("Full search completed", true);
        printValue("Total time (ms)", elapsed_ms);
        printValue("Candidates found", results.num_candidates_found);
        printValue("Stars in corridor", results.num_stars_in_corridor);
        printValue("Events calculated", results.num_events_calculated);
        printValue("Events failed", results.num_events_failed);
        printValue("Phase1 time (ms)", results.phase1_time_ms);
        printValue("Phase2 time (ms)", results.phase2_time_ms);
        printValue("Total time (ms)", results.total_time_ms);
        
        // Verifica performance
        if (results.phase1_time_ms > 0 && results.phase2_time_ms > 0) {
            printTest("Performance timing valid", true);
        } else {
            printTest("Performance timing valid", false);
            allPassed = false;
        }
        
        // Mostra statistiche
        if (results.num_candidates_found > 0) {
            double avg_time_per_candidate = results.phase2_time_ms / results.num_candidates_found;
            printValue("Avg time per candidate (ms)", avg_time_per_candidate, 2);
        }
        
        // Mostra eventi
        if (results.events.size() > 0) {
            std::cout << "\n" << GREEN << "Eventi trovati: " << results.events.size() << RESET << "\n";
            
            // Conta eventi con durata > 1 secondo
            int long_duration_count = 0;
            for (const auto& event : results.events) {
                if (event.maxDuration > 1.0) {
                    long_duration_count++;
                }
            }
            printValue("Long duration events (>1 sec)", long_duration_count);
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Exception: " << e.what() << RESET << "\n";
        allPassed = false;
    }
    
    return allPassed;
}

// ============================================================================
// TEST 4: Test con asteroide diverso
// ============================================================================

bool testDifferentAsteroid() {
    printHeader("TEST 4: Test con Asteroide Diverso");
    
    bool allPassed = true;
    
    try {
        OccultationSearchAstDyn searcher;
        
        // Prova con un asteroide più grande (Ceres)
        std::cout << "Caricamento asteroide 1 (Ceres)...\n";
        if (searcher.loadAsteroid(1)) {
            printTest("Load Ceres", true);
            
            searcher.setAsteroidDiameter(939.4);  // Ceres diameter in km
            searcher.setOrbitalUncertainty(10.0);  // 10 km uncertainty (più preciso)
            
            // Intervallo breve
            double start_mjd = 61007.0;
            double end_mjd = 61008.0;
            
            std::cout << "Esecuzione ricerca per Ceres...\n";
            auto results = searcher.search(start_mjd, end_mjd, 18.0);
            
            printTest("Search for Ceres", true);
            printValue("Candidates found", results.num_candidates_found);
            printValue("Events calculated", results.num_events_calculated);
            
        } else {
            printTest("Load Ceres", false);
            allPassed = false;
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Exception: " << e.what() << RESET << "\n";
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
    std::cout << "║  TEST OCCULTATION SEARCH                                 ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << RESET << "\n";
    
    // Inizializza UnifiedGaiaCatalog
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
    int testsTotal = 4;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // Esegui tutti i test
    if (testLoadAsteroid()) testsPassed++;
    if (testSimplifiedSearch()) testsPassed++;
    if (testFullSearch()) testsPassed++;
    if (testDifferentAsteroid()) testsPassed++;
    
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    
    // Cleanup
    ioc::gaia::UnifiedGaiaCatalog::shutdown();
    
    // Riepilogo
    std::cout << "\n";
    std::cout << BLUE << "═══════════════════════════════════════════════════════\n";
    std::cout << "RIEPILOGO TEST\n";
    std::cout << "═══════════════════════════════════════════════════════" << RESET << "\n\n";
    
    std::cout << "Test passati: " << GREEN << testsPassed << "/" << testsTotal << RESET << "\n";
    std::cout << "Tempo totale: " << std::fixed << std::setprecision(2) << elapsed << " secondi\n\n";
    
    if (testsPassed == testsTotal) {
        std::cout << GREEN << "✓ TUTTI I TEST PASSATI!" << RESET << "\n\n";
        return 0;
    } else {
        std::cout << RED << "✗ ALCUNI TEST FALLITI" << RESET << "\n\n";
        return 1;
    }
}

