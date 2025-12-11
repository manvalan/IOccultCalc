/**
 * @file test_occultation_search_17030.cpp
 * @brief Test specifico per asteroide 17030 (26-30 novembre 2025)
 * 
 * Verifica occultazione per stella GAIA:3411546266140512128
 * 
 * Compilazione:
 *   cd build && make test_occultation_search_17030
 * 
 * Esecuzione:
 *   ./tests/test_occultation_search_17030
 */

#include "ioccultcalc/occultation_search_astdyn.h"
#include "ioccultcalc/time_utils.h"
#include "ioccultcalc/astdyn_propagation_helper.h"
#include "ioccultcalc/types.h"
#include "ioc_gaialib/unified_gaia_catalog.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <fstream>
#include <cstdlib>
#include <filesystem>

using namespace ioccultcalc;

// Colori per output
#define RESET   "\033[0m"
#define GREEN   "\033[32m"
#define RED     "\033[31m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define CYAN    "\033[36m"
#define MAGENTA "\033[35m"

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

int main() {
    std::cout << GREEN << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST OCCULTATION SEARCH - ASTEROIDE 17030               ║\n";
    std::cout << "║  26-30 Novembre 2025                                     ║\n";
    std::cout << "║  Target: GAIA:3411546266140512128                        ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << RESET << "\n";
    
    // Inizializza UnifiedGaiaCatalog con catalogo locale
    std::cout << "Initializing UnifiedGaiaCatalog with local catalog...\n";
    
    std::string home = std::string(getenv("HOME"));
    std::string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
    
    std::string json_config = R"({
        "catalog_type": "multifile_v2",
        "multifile_directory": ")" + catalog_path + R"(",
        "max_cached_chunks": 50,
        "log_level": "warning"
    })";
    
    std::cout << "  Using catalog: " << catalog_path << "\n";
    
    if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) {
        std::cerr << RED << "✗ Failed to initialize UnifiedGaiaCatalog\n" << RESET;
        std::cerr << "  Check if catalog exists at: " << catalog_path << "\n";
        return 1;
    }
    
    // Ottieni l'istanza solo dopo l'inizializzazione riuscita
    auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
    auto info = catalog.getCatalogInfo();
    std::cout << GREEN << "✓ UnifiedGaiaCatalog initialized\n" << RESET;
    std::cout << "  Catalog: " << info.catalog_name << "\n";
    std::cout << "  Total stars: " << info.total_stars << "\n";
    std::cout << "  Magnitude limit: " << info.magnitude_limit << "\n";
    std::cout << "  Online: " << (info.is_online ? "Yes" : "No") << "\n\n";
    
    bool allPassed = true;
    
    try {
        // Converti date in MJD TDB
        // 26 novembre 2025 00:00 UTC
        JulianDate startJD = TimeUtils::calendarToJD(2025, 11, 26, 0, 0, 0);
        // 30 novembre 2025 23:59 UTC
        JulianDate endJD = TimeUtils::calendarToJD(2025, 11, 30, 23, 59, 0);
        
        double start_mjd = startJD.jd - 2400000.5;
        double end_mjd = endJD.jd - 2400000.5;
        
        printHeader("Configurazione Test");
        printInfo("Asteroide", "17030");
        printValue("Start date", start_mjd, 5);
        printValue("End date", end_mjd, 5);
        printInfo("Target star", "GAIA:3411546266140512128");
        
        // Crea searcher
        OccultationSearchAstDyn searcher;
        
        // Carica asteroide 17030
        printHeader("Caricamento Asteroide 17030");
        if (!searcher.loadAsteroid(17030)) {
            std::cerr << RED << "✗ Failed to load asteroid 17030" << RESET << "\n";
            ioc::gaia::UnifiedGaiaCatalog::shutdown();
            return 1;
        }
        printTest("Load asteroid 17030", true);
        
        auto elements = searcher.getAsteroidElements();
        printInfo("Name", elements.name);
        printValue("Number", elements.number);
        printValue("Semi-major axis (AU)", elements.a, 6);
        printValue("Eccentricity", elements.e, 6);
        printValue("Inclination (deg)", elements.i, 4);
        printValue("Epoch (MJD TDB)", elements.epoch_mjd, 2);
        
        // Configura asteroide (diametro tipico per un asteroide di questa dimensione)
        // 17030 è un asteroide di medie dimensioni, assumiamo ~10-20 km
        searcher.setAsteroidDiameter(15.0);  // km
        searcher.setOrbitalUncertainty(50.0);  // 50 km uncertainty
        
        // Esegui ricerca
        printHeader("Esecuzione Ricerca");
        std::cout << "Cercando occultazioni dal " << start_mjd << " al " << end_mjd << " MJD TDB...\n";
        std::cout << "  (26-30 novembre 2025, ~5 giorni)\n\n";
        
        auto t_start = std::chrono::high_resolution_clock::now();
        
        OccultationSearchConfig config;
        config.start_mjd_tdb = start_mjd;
        config.end_mjd_tdb = end_mjd;
        config.max_magnitude = 18.0;
        config.asteroid_diameter_km = 15.0;
        config.orbital_uncertainty_km = 50.0;
        config.min_probability = 0.01;
        
        // Configurazione Phase1
        config.phase1_config = Phase1Config::conservative();
        config.phase1_config.start_mjd_tdb = start_mjd;
        config.phase1_config.end_mjd_tdb = end_mjd;
        config.phase1_config.max_magnitude = 18.0;
        config.phase1_config.corridor_width_deg = 0.0083;  // 30 arcsec
        
        // Configurazione Phase2
        config.phase2_config.time_window_minutes = 10.0;
        config.phase2_config.time_step_seconds = 1.0;
        
        std::cout << "Avvio ricerca (Phase1 + Phase2)...\n";
        std::cout << "  Phase1: screening candidati (query Gaia corridor)\n";
        std::cout << "  Phase2: calcolo geometria precisa\n";
        std::cout << "  Questo può richiedere tempo se il catalogo è online...\n\n";
        std::cout.flush();
        
        auto results = searcher.search(config);
        
        auto t_end = std::chrono::high_resolution_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        
        printHeader("Risultati Ricerca");
        printValue("Total time (ms)", elapsed_ms);
        printValue("Candidates found", results.num_candidates_found);
        printValue("Stars in corridor", results.num_stars_in_corridor);
        printValue("Events calculated", results.num_events_calculated);
        printValue("Events failed", results.num_events_failed);
        printValue("Phase1 time (ms)", results.phase1_time_ms);
        printValue("Phase2 time (ms)", results.phase2_time_ms);
        printValue("Total time (ms)", results.total_time_ms);
        
        // Calcola distanza esplicita tra asteroide e stella target
        printHeader("Calcolo Distanza Asteroide-Stella");
        uint64_t target_star_id = 3411546266140512128ULL;
        std::string target_star_id_str = std::to_string(target_star_id);
        
        // Ottieni stella dal catalogo
        auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
        auto star_opt = catalog.queryBySourceId(target_star_id);
        
        if (star_opt.has_value()) {
            const auto& star = star_opt.value();
            std::cout << "Stella target trovata nel catalogo:\n";
            printInfo("  Source ID", std::to_string(star.source_id));
            printValue("  RA (deg)", star.ra, 6);
            printValue("  Dec (deg)", star.dec, 6);
            printValue("  G magnitude", star.phot_g_mean_mag, 2);
            
            // Propaga asteroide a metà periodo
            double mid_mjd = (start_mjd + end_mjd) / 2.0;
            auto& helper = AstDynPropagationHelper::getInstance();
            auto asteroid_ra_dec = helper.getRADec(elements, mid_mjd);
            
            std::cout << "\nPosizione asteroide a MJD " << std::fixed << std::setprecision(5) << mid_mjd << ":\n";
            printValue("  RA (deg)", asteroid_ra_dec.first, 6);
            printValue("  Dec (deg)", asteroid_ra_dec.second, 6);
            
            // Calcola distanza angolare (haversine)
            double ra1_rad = star.ra * DEG_TO_RAD;
            double dec1_rad = star.dec * DEG_TO_RAD;
            double ra2_rad = asteroid_ra_dec.first * DEG_TO_RAD;
            double dec2_rad = asteroid_ra_dec.second * DEG_TO_RAD;
            
            double delta_ra = ra2_rad - ra1_rad;
            double delta_dec = dec2_rad - dec1_rad;
            
            double a = std::sin(delta_dec / 2) * std::sin(delta_dec / 2) +
                       std::cos(dec1_rad) * std::cos(dec2_rad) * 
                       std::sin(delta_ra / 2) * std::sin(delta_ra / 2);
            double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
            double separation_arcsec = c * RAD_TO_DEG * 3600.0;
            
            std::cout << "\nDistanza angolare:\n";
            printValue("  Separation (arcsec)", separation_arcsec, 3);
            printValue("  Separation (arcmin)", separation_arcsec / 60.0, 3);
            printValue("  Separation (deg)", separation_arcsec / 3600.0, 6);
            
            // Verifica se è nel corridoio
            double corridor_width_arcsec = config.phase1_config.corridor_width_deg * 3600.0;
            bool in_corridor = separation_arcsec <= corridor_width_arcsec;
            std::cout << "\n  Corridor width: " << corridor_width_arcsec << " arcsec\n";
            std::cout << "  " << (in_corridor ? std::string(GREEN) + "✓" : std::string(RED) + "✗") << RESET 
                     << " Stella " << (in_corridor ? "dentro" : "fuori") << " corridoio\n";
        } else {
            std::cout << YELLOW << "⚠ Stella target non trovata nel catalogo" << RESET << "\n";
        }
        
        std::cout << "\n";
        
        // Cerca l'evento per la stella target
        printHeader("Verifica Evento Target");
        bool found_target = false;
        
        std::cout << "Cercando evento per stella GAIA:" << target_star_id << "...\n";
        std::cout << "  Total events found: " << results.events.size() << "\n\n";
        
        // Cerca prima se la stella è tra gli eventi
        for (size_t i = 0; i < results.events.size(); ++i) {
            const auto& event = results.events[i];
            
            // Verifica se è la stella target
            bool is_target = (event.star.sourceId == target_star_id_str || 
                            event.eventId == target_star_id_str ||
                            event.eventId.find(target_star_id_str) != std::string::npos);
            
            if (is_target) {
                std::cout << GREEN << ">>> TARGET STAR FOUND! <<<" << RESET << "\n";
                std::cout << MAGENTA << "Evento " << (i+1) << ":" << RESET << "\n";
                printInfo("  Event ID", event.eventId);
                printInfo("  Star ID", event.star.sourceId);
                printValue("  Time CA (JD)", event.timeCA.jd, 5);
                
                // Converti JD in data leggibile
                JulianDate jd = event.timeCA;
                int year, month, day, hour, minute;
                double second;
                TimeUtils::jdToCalendar(jd, year, month, day, hour, minute, second);
                std::cout << "  " << CYAN << "Time CA (UTC): " << RESET 
                         << year << "-" << std::setfill('0') << std::setw(2) << month 
                         << "-" << std::setw(2) << day << " "
                         << std::setw(2) << hour << ":" << std::setw(2) << minute 
                         << ":" << std::setw(2) << static_cast<int>(second) << "\n";
                
                printValue("  Separation (arcsec)", event.closeApproachDistance, 3);
                printValue("  Probability", event.probability, 4);
                printValue("  Max duration (sec)", event.maxDuration, 2);
                printValue("  Position angle (deg)", event.positionAngle, 2);
                found_target = true;
                std::cout << "\n";
                break;
            }
        }
        
        // Mostra i primi 10 eventi per debug
        if (!found_target && results.events.size() > 0) {
            std::cout << YELLOW << "Target star not found. Showing first 10 events for debug:" << RESET << "\n\n";
            for (size_t i = 0; i < std::min(results.events.size(), size_t(10)); ++i) {
                const auto& event = results.events[i];
                std::cout << MAGENTA << "Evento " << (i+1) << ":" << RESET << "\n";
                printInfo("  Event ID", event.eventId);
                printInfo("  Star ID", event.star.sourceId);
                printValue("  Separation (arcsec)", event.closeApproachDistance, 3);
                std::cout << "\n";
            }
        }
        
        if (found_target) {
            printTest("Target star event found", true);
        } else {
            printTest("Target star event found", false);
            std::cout << YELLOW << "  ⚠ Evento per stella target non trovato" << RESET << "\n";
            std::cout << "  Totale eventi trovati: " << results.events.size() << "\n";
            allPassed = false;
        }
        
        // Statistiche generali
        if (results.events.size() > 0) {
            printHeader("Statistiche Eventi");
            int high_prob_count = 0;
            int long_duration_count = 0;
            double min_separation = 1e9;
            double max_probability = 0.0;
            
            for (const auto& event : results.events) {
                if (event.probability > 0.5) high_prob_count++;
                if (event.maxDuration > 1.0) long_duration_count++;
                if (event.closeApproachDistance < min_separation) {
                    min_separation = event.closeApproachDistance;
                }
                if (event.probability > max_probability) {
                    max_probability = event.probability;
                }
            }
            
            printValue("Total events", results.events.size());
            printValue("High probability (>0.5)", high_prob_count);
            printValue("Long duration (>1 sec)", long_duration_count);
            printValue("Min separation (arcsec)", min_separation, 3);
            printValue("Max probability", max_probability, 4);
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Exception: " << e.what() << RESET << "\n";
        allPassed = false;
    }
    
    // Cleanup
    ioc::gaia::UnifiedGaiaCatalog::shutdown();
    
    // Riepilogo finale
    printHeader("Riepilogo");
    if (allPassed) {
        std::cout << GREEN << "✓ TEST COMPLETATO CON SUCCESSO!" << RESET << "\n\n";
        return 0;
    } else {
        std::cout << YELLOW << "⚠ TEST COMPLETATO CON AVVISI" << RESET << "\n\n";
        return 1;
    }
}

