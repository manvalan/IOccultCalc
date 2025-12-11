/**
 * @file test_occultation_17030_28nov.cpp
 * @brief Test ricerca occultazione asteroide 17030 - 28 novembre 2025
 */

#include "ioccultcalc/occultation_search_astdyn.h"
#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/time_utils.h"
#include "ioccultcalc/astdyn_propagation_helper.h"
#include "ioc_gaialib/unified_gaia_catalog.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <cstdlib>

using namespace ioccultcalc;
using namespace ioc::gaia;

// Colori per output
const char* GREEN = "\033[32m";
const char* YELLOW = "\033[33m";
const char* RED = "\033[31m";
const char* BLUE = "\033[34m";
const char* CYAN = "\033[36m";
const char* RESET = "\033[0m";

void printHeader(const std::string& text) {
    std::cout << "\n" << CYAN << "═══════════════════════════════════════════════════════\n";
    std::cout << text << "\n";
    std::cout << "═══════════════════════════════════════════════════════" << RESET << "\n";
}

void printInfo(const std::string& label, const std::string& value) {
    std::cout << "  " << CYAN << label << ": " << RESET << value << "\n";
}

void printValue(const std::string& label, double value, int precision = 6) {
    std::cout << "  " << CYAN << label << ": " << RESET 
              << std::fixed << std::setprecision(precision) << value << "\n";
}

int main() {
    std::cout << GREEN << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST OCCULTATION SEARCH - ASTEROIDE 17030               ║\n";
    std::cout << "║  28 Novembre 2025 (finestra ristretta)                  ║\n";
    std::cout << "║  Target: GAIA:3411546266140512128                        ║\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n" << RESET;
    
    // Inizializza Gaia Catalog
    std::cout << "\nInitializing UnifiedGaiaCatalog with local catalog...\n";
    
    std::string home = std::getenv("HOME") ? std::getenv("HOME") : "";
    std::string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
    
    // Risolvi symlink se presente
    if (std::filesystem::exists(catalog_path) && std::filesystem::is_symlink(catalog_path)) {
        catalog_path = std::filesystem::canonical(catalog_path);
    }
    
    std::string json_config = R"({
        "catalog_type": "multifile_v2",
        "multifile_directory": ")" + catalog_path + R"(",
        "max_cached_chunks": 50,
        "log_level": "warning"
    })";
    
    if (!UnifiedGaiaCatalog::initialize(json_config)) {
        std::cerr << RED << "❌ Failed to initialize UnifiedGaiaCatalog" << RESET << "\n";
        return 1;
    }
    
    auto& catalog = UnifiedGaiaCatalog::getInstance();
    auto info = catalog.getCatalogInfo();
    std::cout << GREEN << "✓ UnifiedGaiaCatalog initialized" << RESET << "\n";
    std::cout << "  Catalog: " << info.catalog_name << "\n";
    std::cout << "  Total stars: " << info.total_stars << "\n";
    std::cout << "  Magnitude limit: " << info.magnitude_limit << "\n";
    std::cout << "  Online: " << (info.is_online ? "Yes" : "No") << "\n\n";
    
    // Configurazione test
    printHeader("Configurazione Test");
    
    int asteroid_number = 17030;
    std::string target_star = "GAIA:3411546266140512128";
    
    // 28 novembre 2025: MJD 61007.0 (00:00 UTC) a 61007.99931 (23:59 UTC)
    double start_mjd = 61007.0;  // 28 Nov 2025 00:00 UTC
    double end_mjd = 61008.0;    // 29 Nov 2025 00:00 UTC
    
    printInfo("Asteroide", std::to_string(asteroid_number));
    printValue("Start date (MJD)", start_mjd, 5);
    printValue("End date (MJD)", end_mjd, 5);
    printInfo("Target star", target_star);
    
    // Carica asteroide
    printHeader("Caricamento Asteroide 17030");
    
    OccultationSearchAstDyn search;
    search.setGaiaCatalog(&catalog);
    
    if (!search.loadAsteroid(asteroid_number)) {
        std::cerr << RED << "❌ Failed to load asteroid " << asteroid_number << RESET << "\n";
        return 1;
    }
    
    std::cout << GREEN << "✓ Load asteroid " << asteroid_number << RESET << "\n";
    auto elements = search.getAsteroidElements();
    printValue("Semi-major axis (AU)", elements.a, 6);
    printValue("Eccentricity", elements.e, 6);
    printValue("Inclination (deg)", elements.i * 180.0 / M_PI, 4);
    printValue("Epoch (MJD TDB)", elements.epoch_mjd, 2);
    
    // Configurazione ricerca
    OccultationSearchConfig config;
    config.start_mjd_tdb = start_mjd;
    config.end_mjd_tdb = end_mjd;
    config.max_magnitude = 18.0;
    
    // Phase1: configurazione conservativa per catturare tutti i candidati
    config.phase1_config = Phase1Config::conservative();
    config.phase1_config.start_mjd_tdb = start_mjd;
    config.phase1_config.end_mjd_tdb = end_mjd;
    config.phase1_config.max_magnitude = 18.0;
    
    // Phase2: configurazione precisa
    config.phase2_config.time_window_minutes = 10.0;  // ±10 minuti
    config.min_probability = 0.01;  // 1% minimo
    
    // Esegui ricerca
    printHeader("Esecuzione Ricerca");
    
    std::cout << "Cercando occultazioni dal " << std::fixed << std::setprecision(2) 
              << start_mjd << " al " << end_mjd << " MJD TDB...\n";
    std::cout << "  (28 novembre 2025, ~24 ore)\n\n";
    
    std::cout << "Avvio ricerca (Phase1 + Phase2)...\n";
    std::cout << "  Phase1: screening candidati (query Gaia corridor)\n";
    std::cout << "  Phase2: calcolo geometria precisa\n\n";
    
    auto results = search.search(config);
    
    // Risultati
    printHeader("Risultati Ricerca");
    
    std::cout << "Tempo totale: " << std::fixed << std::setprecision(2) 
              << results.total_time_ms << " ms\n";
    std::cout << "  Phase1: " << results.phase1_time_ms << " ms\n";
    std::cout << "  Phase2: " << results.phase2_time_ms << " ms\n\n";
    
    std::cout << "Stelle nel corridoio: " << results.num_stars_in_corridor << "\n";
    std::cout << "Candidati Phase1: " << results.num_candidates_found << "\n";
    std::cout << "Eventi calcolati: " << results.num_events_calculated << "\n";
    std::cout << "Eventi falliti: " << results.num_events_failed << "\n\n";
    
    // Cerca evento per stella target
    printHeader("Eventi Trovati");
    
    if (results.events.empty()) {
        std::cout << YELLOW << "⚠ Nessun evento trovato" << RESET << "\n";
    } else {
        std::cout << "Trovati " << results.events.size() << " eventi:\n\n";
        
        for (size_t i = 0; i < results.events.size(); ++i) {
            const auto& event = results.events[i];
            std::cout << "Evento " << (i + 1) << ":\n";
            printValue("  Time CA (MJD UTC)", event.timeCA.jd - 2400000.5, 8);
            printValue("  Closest Approach (arcsec)", event.closeApproachDistance, 3);
            printValue("  Position Angle (deg)", event.positionAngle, 2);
            printValue("  Max Duration (sec)", event.maxDuration, 2);
            std::cout << "\n";
        }
    }
    
    // Verifica stella target
    printHeader("Verifica Stella Target");
    
    uint64_t target_source_id = std::stoull(target_star.substr(5));  // Rimuovi "GAIA:"
    bool found_target = false;
    
    for (const auto& event : results.events) {
        // Verifica se l'evento corrisponde alla stella target
        // Nota: potrebbe essere necessario verificare event.star_source_id o event.eventId
        if (event.eventId == std::to_string(target_source_id)) {
            found_target = true;
            std::cout << GREEN << "✓ Evento trovato per stella target!" << RESET << "\n\n";
            printValue("Time CA (MJD UTC)", event.timeCA.jd - 2400000.5, 8);
            printValue("Closest Approach (arcsec)", event.closeApproachDistance, 3);
            printValue("Position Angle (deg)", event.positionAngle, 2);
            printValue("Max Duration (sec)", event.maxDuration, 2);
            
            // Converti MJD UTC a data/ora
            ioccultcalc::JulianDate jd_utc(event.timeCA.jd);
            int year, month, day, hour, minute;
            double second;
            ioccultcalc::TimeUtils::jdToCalendar(jd_utc, year, month, day, hour, minute, second);
            
            std::cout << "\nData/ora evento:\n";
            std::cout << "  " << std::setfill('0') << std::setw(2) << day << "/"
                      << std::setw(2) << month << "/" << year << " "
                      << std::setw(2) << hour << ":" << std::setw(2) << minute << ":"
                      << std::setw(2) << (int)second << " UTC\n";
            
            if (event.closeApproachDistance < 0.1) {
                std::cout << GREEN << "\n✓ OCCULTAZIONE CONFERMATA!" << RESET << "\n";
            } else {
                std::cout << YELLOW << "\n⚠ Avvicinamento ravvicinato ma non occultazione" << RESET << "\n";
            }
            break;
        }
    }
    
    if (!found_target) {
        std::cout << YELLOW << "⚠ Nessun evento trovato per stella target " << target_star << RESET << "\n";
        
        // Verifica se la stella è nel catalogo
        auto star_opt = catalog.queryBySourceId(target_source_id);
        if (star_opt.has_value()) {
            const auto& star = star_opt.value();
            std::cout << "\nStella trovata nel catalogo:\n";
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
            
            // Calcola distanza angolare
            double ra1_rad = star.ra * M_PI / 180.0;
            double dec1_rad = star.dec * M_PI / 180.0;
            double ra2_rad = asteroid_ra_dec.first * M_PI / 180.0;
            double dec2_rad = asteroid_ra_dec.second * M_PI / 180.0;
            
            double delta_ra = ra2_rad - ra1_rad;
            double delta_dec = dec2_rad - dec1_rad;
            
            double a = std::sin(delta_dec / 2) * std::sin(delta_dec / 2) +
                       std::cos(dec1_rad) * std::cos(dec2_rad) * 
                       std::sin(delta_ra / 2) * std::sin(delta_ra / 2);
            double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
            double separation_arcsec = c * 180.0 / M_PI * 3600.0;
            
            std::cout << "\nDistanza angolare:\n";
            printValue("  Separation (arcsec)", separation_arcsec, 3);
            printValue("  Separation (arcmin)", separation_arcsec / 60.0, 3);
        } else {
            std::cout << RED << "\n❌ Stella non trovata nel catalogo" << RESET << "\n";
        }
    }
    
    std::cout << "\n";
    return 0;
}

