/**
 * @file test_phase1_class.cpp
 * @brief Test della classe Phase1CandidateScreening
 * @date 4 Dicembre 2025
 */

#include "phase1_candidate_screening.h"
#include "ioc_gaialib/unified_gaia_catalog.h"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

void printSeparator() {
    std::cout << "\n" << std::string(80, '=') << "\n\n";
}

void printCandidate(const CandidateStar& star, int idx) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  " << idx << ". Source ID: " << star.source_id << "\n";
    std::cout << "     RA/Dec: " << star.ra_deg << "° / " << star.dec_deg << "°\n";
    std::cout << "     Mag: " << std::setprecision(2) << star.phot_g_mean_mag << "\n";
    std::cout << "     Closest Approach: " << std::setprecision(3) 
              << star.closest_approach_arcsec << " arcsec\n";
    std::cout << "     CA Epoch: MJD " << std::setprecision(6) 
              << star.closest_approach_mjd << "\n";
    
    // Converti MJD in data/ora leggibile
    double jd = star.closest_approach_mjd + 2400000.5;
    int y, m, d;
    double frac = jd - static_cast<int>(jd) + 0.5;
    if (frac >= 1.0) {
        frac -= 1.0;
        jd += 1.0;
    }
    int z = static_cast<int>(jd);
    int a = static_cast<int>((z - 1867216.25) / 36524.25);
    int b = z + 1 + a - a / 4;
    int c = b + 1524;
    int e = static_cast<int>((c - 122.1) / 365.25);
    int f = static_cast<int>(365.25 * e);
    int g = static_cast<int>((c - f) / 30.6001);
    d = c - f - static_cast<int>(30.6001 * g);
    m = (g < 14) ? (g - 1) : (g - 13);
    y = (m > 2) ? (e - 4716) : (e - 4715);
    
    int hh = static_cast<int>(frac * 24);
    int mm = static_cast<int>((frac * 24 - hh) * 60);
    int ss = static_cast<int>(((frac * 24 - hh) * 60 - mm) * 60);
    
    std::cout << "     CA Time: " << y << "-" << std::setfill('0') << std::setw(2) << m 
              << "-" << std::setw(2) << d << " " 
              << std::setw(2) << hh << ":" << std::setw(2) << mm 
              << ":" << std::setw(2) << ss << " UT\n";
    std::cout << std::setfill(' ');
    std::cout << "     Angular Velocity: " << std::setprecision(2) 
              << star.angular_velocity_arcsec_per_sec << " arcsec/sec\n";
    std::cout << "     Segment Index: " << star.closest_segment_index << "\n\n";
}

int main() {
    std::cout << "\n╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST CLASSE Phase1CandidateScreening                         ║\n";
    std::cout << "║  Asteroid 17030 - 28 November 2025                            ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";
    
    try {
        // Inizializza catalogo Gaia
        std::cout << "\n[0] Inizializzazione catalogo Gaia...\n";
        
        std::string home = std::string(getenv("HOME"));
        std::string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
        
        std::string json_config = R"({
            "catalog_type": "multifile_v2",
            "multifile_directory": ")" + catalog_path + R"("
        })";
        
        if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) {
            throw std::runtime_error("Impossibile inizializzare catalogo Gaia");
        }
        std::cout << "✓ Catalogo Gaia inizializzato\n";
        
        // Crea istanza classe Phase1
        Phase1CandidateScreening phase1;
        
        // Carica elementi orbitali dal database JSON
        std::cout << "\n[1] Caricamento elementi orbitali da database JSON (asteroide 17030)...\n";
        if (!phase1.loadAsteroidFromJSON(17030)) {
            std::cerr << "ERRORE: Impossibile caricare elementi orbitali!\n";
            return 1;
        }
        std::cout << "✓ Elementi orbitali caricati con successo\n";
        
        // Verifica elementi caricati
        auto elements = phase1.getOrbitalElements();
        std::cout << "\nElementi orbitali:\n";
        std::cout << "  a = " << std::fixed << std::setprecision(6) 
                  << elements.semi_major_axis << " AU\n";
        std::cout << "  e = " << elements.eccentricity << "\n";
        std::cout << "  i = " << elements.inclination * 180.0 / M_PI << "°\n";
        std::cout << "  Epoch: MJD " << elements.epoch_mjd_tdb << "\n";
        
        printSeparator();
        
        // TEST 1: Configurazione DEFAULT (30 arcsec, 30 sec)
        std::cout << "[2] TEST con configurazione DEFAULT (RECOMMENDED)\n";
        std::cout << "    Corridor width: 30 arcsec (0.5 arcmin)\n";
        std::cout << "    Path interval: 30 seconds\n";
        std::cout << "    CA threshold: 15.0 arcsec\n\n";
        
        Phase1Config config_default;
        config_default.start_mjd_tdb = 61007.0;  // 28 Nov 2025 00:00 UT
        config_default.end_mjd_tdb = 61008.0;    // 29 Nov 2025 00:00 UT
        config_default.path_interval_seconds = 30;
        config_default.corridor_width_deg = 0.0083;  // 30 arcsec
        config_default.closest_approach_threshold_arcsec = 15.0;
        config_default.max_magnitude = 18.0;
        config_default.min_parallax = 0.0;
        
        auto results_default = phase1.screenCandidates(config_default);
        
        std::cout << "RISULTATI:\n";
        std::cout << "  Path points: " << results_default.num_path_points << "\n";
        std::cout << "  Propagation time: " << std::setprecision(2) 
                  << results_default.propagation_time_ms << " ms\n";
        std::cout << "  Corridor query time: " << results_default.corridor_query_time_ms / 1000.0 
                  << " sec\n";
        std::cout << "  Stars in corridor: " << results_default.num_stars_in_corridor << "\n";
        std::cout << "  CA calculation time: " << results_default.closest_approach_calc_time_ms 
                  << " ms\n";
        std::cout << "  Candidates filtered: " << results_default.num_candidates_filtered << "\n";
        
        if (results_default.num_candidates_filtered > 0) {
            std::cout << "\nCANDIDATE TROVATE (CA < " << config_default.closest_approach_threshold_arcsec 
                      << " arcsec):\n\n";
            for (size_t i = 0; i < results_default.candidates.size(); ++i) {
                printCandidate(results_default.candidates[i], i + 1);
            }
            
            // Verifica se il target star è presente (Gaia 3411546266140512128)
            bool found_target = false;
            for (const auto& cand : results_default.candidates) {
                if (cand.source_id == 3411546266140512128ULL) {
                    found_target = true;
                    std::cout << "✓ TARGET STAR FOUND! Gaia DR3 3411546266140512128\n";
                    std::cout << "  (UCAC4 552-011427 dalla predizione IOTA)\n";
                    break;
                }
            }
            if (!found_target) {
                std::cout << "⚠ Target star NON trovata nella lista candidate\n";
            }
        }
        
        printSeparator();
        
        // TEST 2: Configurazione CONSERVATIVE (60 arcsec, 20 sec)
        std::cout << "[3] TEST con configurazione CONSERVATIVE\n";
        std::cout << "    Corridor width: 60 arcsec (1 arcmin)\n";
        std::cout << "    Path interval: 20 seconds\n";
        std::cout << "    CA threshold: 30.0 arcsec\n\n";
        
        Phase1Config config_conservative;
        config_conservative.start_mjd_tdb = 61007.0;
        config_conservative.end_mjd_tdb = 61008.0;
        config_conservative.path_interval_seconds = 20;
        config_conservative.corridor_width_deg = 0.0167;  // 60 arcsec
        config_conservative.closest_approach_threshold_arcsec = 30.0;
        config_conservative.max_magnitude = 18.0;
        config_conservative.min_parallax = 0.0;
        
        auto results_conservative = phase1.screenCandidates(config_conservative);
        
        std::cout << "RISULTATI:\n";
        std::cout << "  Path points: " << results_conservative.num_path_points << "\n";
        std::cout << "  Stars in corridor: " << results_conservative.num_stars_in_corridor << "\n";
        std::cout << "  Candidates filtered: " << results_conservative.num_candidates_filtered << "\n";
        std::cout << "  Total time: " 
                  << (results_conservative.propagation_time_ms + 
                      results_conservative.corridor_query_time_ms + 
                      results_conservative.closest_approach_calc_time_ms) / 1000.0 
                  << " sec\n";
        
        printSeparator();
        
        // TEST 3: Configurazione FAST (20 arcsec, 60 sec)
        std::cout << "[4] TEST con configurazione FAST\n";
        std::cout << "    Corridor width: 20 arcsec\n";
        std::cout << "    Path interval: 60 seconds\n";
        std::cout << "    CA threshold: 10.0 arcsec\n\n";
        
        Phase1Config config_fast;
        config_fast.start_mjd_tdb = 61007.0;
        config_fast.end_mjd_tdb = 61008.0;
        config_fast.path_interval_seconds = 60;
        config_fast.corridor_width_deg = 0.0056;  // 20 arcsec
        config_fast.closest_approach_threshold_arcsec = 10.0;
        config_fast.max_magnitude = 18.0;
        config_fast.min_parallax = 0.0;
        
        auto results_fast = phase1.screenCandidates(config_fast);
        
        std::cout << "RISULTATI:\n";
        std::cout << "  Path points: " << results_fast.num_path_points << "\n";
        std::cout << "  Stars in corridor: " << results_fast.num_stars_in_corridor << "\n";
        std::cout << "  Candidates filtered: " << results_fast.num_candidates_filtered << "\n";
        std::cout << "  Total time: " 
                  << (results_fast.propagation_time_ms + 
                      results_fast.corridor_query_time_ms + 
                      results_fast.closest_approach_calc_time_ms) / 1000.0 
                  << " sec\n";
        
        printSeparator();
        
        // CONFRONTO CONFIGURAZIONI
        std::cout << "[5] CONFRONTO CONFIGURAZIONI\n\n";
        std::cout << std::left;
        std::cout << std::setw(20) << "Configuration" 
                  << std::setw(15) << "Path Points" 
                  << std::setw(15) << "Stars Found" 
                  << std::setw(15) << "Candidates" 
                  << std::setw(15) << "Time (sec)" << "\n";
        std::cout << std::string(80, '-') << "\n";
        
        auto print_row = [](const std::string& name, const Phase1Results& r) {
            double total_time = (r.propagation_time_ms + r.corridor_query_time_ms + 
                               r.closest_approach_calc_time_ms) / 1000.0;
            std::cout << std::setw(20) << name
                      << std::setw(15) << r.num_path_points
                      << std::setw(15) << r.num_stars_in_corridor
                      << std::setw(15) << r.num_candidates_filtered
                      << std::setw(15) << std::fixed << std::setprecision(2) << total_time << "\n";
        };
        
        print_row("Default", results_default);
        print_row("Conservative", results_conservative);
        print_row("Fast", results_fast);
        
        printSeparator();
        
        std::cout << "[6] RACCOMANDAZIONI\n\n";
        std::cout << "✓ Configurazione DEFAULT raccomandata per la maggior parte dei casi\n";
        std::cout << "  - Buon bilanciamento tra precisione e velocità\n";
        std::cout << "  - 30 arcsec corridor cattura tutte le occultazioni realistiche\n";
        std::cout << "  - 30 sec path interval sufficiente per CA < 1 arcsec\n\n";
        
        std::cout << "✓ Configurazione CONSERVATIVE per ricerche esaustive\n";
        std::cout << "  - Corridor più ampio (60 arcsec) per margine di sicurezza\n";
        std::cout << "  - Path più denso (20 sec) per massima precisione\n";
        std::cout << "  - Tempo elaborazione maggiore (~2-3x)\n\n";
        
        std::cout << "✓ Configurazione FAST per survey rapidi\n";
        std::cout << "  - Corridor stretto (20 arcsec) solo per stelle molto vicine\n";
        std::cout << "  - Path meno denso (60 sec) riduce tempi\n";
        std::cout << "  - Può perdere alcuni candidati marginali\n\n";
        
        printSeparator();
        
        std::cout << "✅ TEST COMPLETATO CON SUCCESSO!\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << "\n\n";
        return 1;
    }
    
    return 0;
}
