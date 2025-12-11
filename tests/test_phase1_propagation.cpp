/**
 * @file test_phase1_propagation.cpp
 * @brief Test propagazione Phase1 per asteroide 17030
 * 
 * Verifica che Phase1 propaga correttamente gli elementi
 */

#include "phase1_candidate_screening.h"
#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/astdyn_propagation_helper.h"
#include "ioccultcalc/time_utils.h"
#include "ioc_gaialib/unified_gaia_catalog.h"
#include "../external/ITALOccultLibrary/astdyn/include/astdyn/propagation/OrbitalElements.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
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

void printValue(const std::string& label, double value, int precision = 6) {
    std::cout << "  " << CYAN << label << ": " << RESET 
              << std::fixed << std::setprecision(precision) << value << "\n";
}

void printVector(const std::string& label, const Eigen::Vector3d& v, int precision = 6) {
    std::cout << "  " << CYAN << label << ": " << RESET 
              << "[" << std::fixed << std::setprecision(precision) 
              << v[0] << ", " << v[1] << ", " << v[2] << "] AU\n";
}

int main() {
    std::cout << GREEN << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST PROPAGAZIONE PHASE1 - ASTEROIDE 17030         ║\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n" << RESET;
    
    // Inizializza UnifiedGaiaCatalog (necessario per Phase1)
    std::string home = std::getenv("HOME") ? std::getenv("HOME") : "";
    std::string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
    std::string json_config = R"({ "catalog_type": "multifile_v2", "multifile_directory": ")" + catalog_path + R"(", "max_cached_chunks": 50, "log_level": "warning" })";
    if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) {
        std::cerr << RED << "✗ Errore inizializzazione UnifiedGaiaCatalog" << RESET << "\n";
        return 1;
    }
    
    try {
        // Date di test: 26-30 novembre 2025
        JulianDate startJD = TimeUtils::calendarToJD(2025, 11, 26, 0, 0, 0);
        JulianDate endJD = TimeUtils::calendarToJD(2025, 11, 30, 23, 59, 0);
        double start_mjd = startJD.jd - 2400000.5;
        double end_mjd = endJD.jd - 2400000.5;
        double mid_mjd = (start_mjd + end_mjd) / 2.0;
        
        printHeader("Configurazione Test");
        printValue("Start MJD", start_mjd, 5);
        printValue("End MJD", end_mjd, 5);
        printValue("Mid MJD", mid_mjd, 5);
        std::cout << "  " << CYAN << "Mid date: " << RESET << "2025-11-28 ~12:00\n\n";
        
        // Carica elementi
        printHeader("Caricamento Elementi");
        AstDySElements elements = AstDySClient::downloadElements(17030);
        if (elements.a <= 0) {
            std::cerr << RED << "✗ Errore caricamento elementi" << RESET << "\n";
            return 1;
        }
        std::cout << GREEN << "✓ Elementi caricati\n" << RESET;
        printValue("Epoch (MJD TDB)", elements.epoch_mjd, 5);
        printValue("a (AU)", elements.a, 8);
        printValue("e", elements.e, 8);
        
        // Crea Phase1CandidateScreening
        printHeader("Phase1CandidateScreening");
        Phase1CandidateScreening phase1;
        
        // Carica elementi in Phase1 usando AstDyS (come fa OccultationSearchAstDyn)
        // Converti AstDySElements in KeplerianElements
        astdyn::propagation::KeplerianElements kep;
        kep.semi_major_axis = elements.a;
        kep.eccentricity = elements.e;
        kep.inclination = elements.i * M_PI / 180.0;
        kep.longitude_ascending_node = elements.Omega * M_PI / 180.0;
        kep.argument_perihelion = elements.omega * M_PI / 180.0;
        kep.mean_anomaly = elements.M * M_PI / 180.0;
        kep.epoch_mjd_tdb = elements.epoch_mjd;
        kep.gravitational_parameter = 
            1.32712440018e20 / std::pow(1.495978707e11, 3) * std::pow(86400.0, 2);
        
        phase1.setOrbitalElements(kep);
        std::cout << GREEN << "✓ Asteroide caricato in Phase1\n" << RESET;
        
        // Configura Phase1
        Phase1Config config;
        config.start_mjd_tdb = start_mjd;
        config.end_mjd_tdb = end_mjd;
        config.path_interval_seconds = 30;
        config.corridor_width_deg = 0.0083;  // 30 arcsec
        config.max_magnitude = 18.0;
        config.closest_approach_threshold_arcsec = 15.0;
        
        printHeader("Esecuzione Phase1 ScreenCandidates");
        auto results = phase1.screenCandidates(config);
        
        std::cout << "  " << CYAN << "Risultati Phase1:\n" << RESET;
        printValue("  Numero punti path", results.num_path_points);
        printValue("  Stelle nel corridoio", results.num_stars_in_corridor);
        printValue("  Candidati filtrati", results.num_candidates_filtered);
        printValue("  Tempo propagazione (ms)", results.propagation_time_ms, 2);
        
        const auto& path = results.path;
        
        // Mostra alcuni punti del path
        std::cout << "\n  " << CYAN << "Primi 3 punti del path:\n" << RESET;
        for (size_t i = 0; i < std::min(path.size(), size_t(3)); ++i) {
            std::cout << "    Punto " << (i+1) << ":\n";
            printValue("      MJD TDB", path[i].mjd_tdb, 8);
            printValue("      RA (deg)", path[i].ra_deg, 8);
            printValue("      Dec (deg)", path[i].dec_deg, 8);
            printValue("      Distance Earth (AU)", path[i].distance_earth_au, 8);
        }
        
        // Verifica punto centrale
        if (path.size() > 0) {
            size_t mid_idx = path.size() / 2;
            std::cout << "\n  " << CYAN << "Punto centrale (indice " << mid_idx << "):\n" << RESET;
            printValue("    MJD TDB", path[mid_idx].mjd_tdb, 8);
            printValue("    RA (deg)", path[mid_idx].ra_deg, 8);
            printValue("    Dec (deg)", path[mid_idx].dec_deg, 8);
            
            // Confronta con propagazione diretta
            printHeader("Confronto con Propagazione Diretta");
            auto& helper = AstDynPropagationHelper::getInstance();
            auto radec_direct = helper.getRADec(elements, path[mid_idx].mjd_tdb);
            
            std::cout << "  " << CYAN << "Phase1 path:\n" << RESET;
            printValue("    RA (deg)", path[mid_idx].ra_deg, 8);
            printValue("    Dec (deg)", path[mid_idx].dec_deg, 8);
            
            std::cout << "  " << CYAN << "Propagazione diretta:\n" << RESET;
            printValue("    RA (deg)", radec_direct.first, 8);
            printValue("    Dec (deg)", radec_direct.second, 8);
            
            double diff_ra = std::abs(path[mid_idx].ra_deg - radec_direct.first);
            double diff_dec = std::abs(path[mid_idx].dec_deg - radec_direct.second);
            
            // Normalizza differenza RA (può essere > 180°)
            if (diff_ra > 180.0) diff_ra = 360.0 - diff_ra;
            
            std::cout << "\n  " << CYAN << "Differenze:\n" << RESET;
            printValue("    RA (arcsec)", diff_ra * 3600.0, 3);
            printValue("    Dec (arcsec)", diff_dec * 3600.0, 3);
            
            if (diff_ra * 3600.0 > 1.0 || diff_dec * 3600.0 > 1.0) {
                std::cout << "\n" << RED << "✗ DISCREPANZA TROVATA!\n" << RESET;
                std::cout << "  La propagazione Phase1 differisce dalla propagazione diretta\n";
            } else {
                std::cout << "\n" << GREEN << "✓ Propagazione Phase1 corretta\n" << RESET;
            }
        }
        
        // Verifica tutti i punti del path
        printHeader("Verifica Tutti i Punti del Path");
        auto& helper = AstDynPropagationHelper::getInstance();
        int discrepancies = 0;
        for (size_t i = 0; i < path.size(); ++i) {
            auto radec_direct = helper.getRADec(elements, path[i].mjd_tdb);
            double diff_ra = std::abs(path[i].ra_deg - radec_direct.first);
            double diff_dec = std::abs(path[i].dec_deg - radec_direct.second);
            
            // Normalizza differenza RA
            if (diff_ra > 180.0) diff_ra = 360.0 - diff_ra;
            
            if (diff_ra * 3600.0 > 1.0 || diff_dec * 3600.0 > 1.0) {
                discrepancies++;
                if (discrepancies <= 3) {  // Mostra solo i primi 3
                    std::cout << "  " << RED << "✗ Punto " << i << ":\n" << RESET;
                    printValue("    Phase1 RA", path[i].ra_deg, 8);
                    printValue("    Direct RA", radec_direct.first, 8);
                    printValue("    Diff RA (arcsec)", diff_ra * 3600.0, 3);
                    printValue("    Phase1 Dec", path[i].dec_deg, 8);
                    printValue("    Direct Dec", radec_direct.second, 8);
                    printValue("    Diff Dec (arcsec)", diff_dec * 3600.0, 3);
                }
            }
        }
        
        if (discrepancies > 0) {
            std::cout << "\n  " << RED << "✗ Trovate " << discrepancies 
                      << " discrepanze su " << path.size() << " punti\n" << RESET;
        } else {
            std::cout << "\n  " << GREEN << "✓ Tutti i punti del path sono corretti\n" << RESET;
        }
        
        // Analizza closest approach delle stelle
        printHeader("Analisi Closest Approach");
        if (results.all_stars.size() > 0) {
            std::vector<double> ca_values;
            for (const auto& star : results.all_stars) {
                ca_values.push_back(star.closest_approach_arcsec);
            }
            std::sort(ca_values.begin(), ca_values.end());
            
            std::cout << "  " << CYAN << "Statistiche closest approach (" 
                      << ca_values.size() << " stelle):\n" << RESET;
            printValue("    Min (arcsec)", ca_values[0], 3);
            printValue("    Max (arcsec)", ca_values.back(), 3);
            printValue("    Mediana (arcsec)", ca_values[ca_values.size()/2], 3);
            printValue("    Soglia filtro (arcsec)", config.closest_approach_threshold_arcsec, 1);
            
            int below_threshold = 0;
            for (double ca : ca_values) {
                if (ca <= config.closest_approach_threshold_arcsec) below_threshold++;
            }
            std::cout << "  " << CYAN << "Stelle sotto soglia: " << RESET 
                      << below_threshold << " / " << ca_values.size() << "\n";
            
            // Mostra prime 10 stelle con CA più piccolo
            std::cout << "\n  " << CYAN << "Prime 10 stelle (CA più piccolo):\n" << RESET;
            for (size_t i = 0; i < std::min(size_t(10), ca_values.size()); ++i) {
                std::cout << "    " << (i+1) << ": " << std::fixed << std::setprecision(3) 
                          << ca_values[i] << " arcsec\n";
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Exception: " << e.what() << RESET << "\n";
        ioc::gaia::UnifiedGaiaCatalog::shutdown();
        return 1;
    }
    
    ioc::gaia::UnifiedGaiaCatalog::shutdown();
    std::cout << "\n" << GREEN << "✓ Test completato\n" << RESET;
    return 0;
}

