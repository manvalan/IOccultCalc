/**
 * Test 17030 Sierks - Pipeline completa Phase1
 * Usa l'API pubblica screenCandidates con catalogo Gaia
 * CARICANDO ELEMENTI DA FILE EQ1 (formato corretto)
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include "phase1_candidate_screening.h"
#include "ioc_gaialib/unified_gaia_catalog.h"

// Per accedere alle costanti AstDyn
#include "astdyn/propagation/OrbitalElements.hpp"
#include "astdyn/core/Constants.hpp"

using namespace ioccultcalc;
using namespace ioc::gaia;

constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;

std::string formatRA(double ra_deg) {
    double hours = ra_deg / 15.0;
    int h = static_cast<int>(hours);
    double mins = (hours - h) * 60.0;
    int m = static_cast<int>(mins);
    double secs = (mins - m) * 60.0;
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << h << "h "
        << std::setw(2) << m << "m "
        << std::fixed << std::setprecision(2) << std::setw(5) << secs << "s";
    return oss.str();
}

std::string formatDec(double dec_deg) {
    char sign = dec_deg >= 0 ? '+' : '-';
    dec_deg = std::abs(dec_deg);
    int d = static_cast<int>(dec_deg);
    double arcmins = (dec_deg - d) * 60.0;
    int am = static_cast<int>(arcmins);
    double arcsecs = (arcmins - am) * 60.0;
    std::ostringstream oss;
    oss << sign << std::setfill('0') << std::setw(2) << d << "d "
        << std::setw(2) << am << "' "
        << std::fixed << std::setprecision(2) << std::setw(5) << arcsecs << "\"";
    return oss.str();
}

int main() {
    std::cout << "=========================================================\n";
    std::cout << "  TEST 17030 SIERKS - PIPELINE COMPLETA PHASE1\n";
    std::cout << "  (Caricamento da file EQ1 originale)\n";
    std::cout << "=========================================================\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // INIZIALIZZAZIONE CATALOGO GAIA
    // ═══════════════════════════════════════════════════════════════
    
    std::cout << "=== INIZIALIZZAZIONE CATALOGO GAIA ===\n\n";
    
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
        if (!UnifiedGaiaCatalog::initialize(json_config)) {
            std::cerr << "Errore inizializzazione catalogo Gaia\n";
            return 1;
        }
        std::cout << "Catalogo Gaia inizializzato con successo!\n\n";
    } catch (const std::exception& e) {
        std::cerr << "Errore: " << e.what() << "\n";
        return 1;
    }
    
    // ═══════════════════════════════════════════════════════════════
    // CARICAMENTO ELEMENTI ORBITALI DA FILE EQ1
    // ═══════════════════════════════════════════════════════════════
    
    // Intervallo: 28 Nov 2025 (MJD 60607) - un giorno
    double start_mjd = 60607.0;
    double end_mjd = 60608.0;
    
    std::cout << "Intervallo ricerca: MJD " << std::fixed << std::setprecision(1) 
              << start_mjd << " - " << end_mjd << "\n";
    std::cout << "  (28-29 Nov 2025)\n\n";
    
    try {
        // Crea oggetto Phase1CandidateScreening
        Phase1CandidateScreening screener;
        
        // CARICA DA FILE EQ1 - questo converte correttamente da equinoziali a kepleriani!
        std::cout << "Caricamento elementi orbitali da 17030_astdys.eq1...\n";
        if (!screener.loadAsteroidFromEQ1("17030_astdys.eq1")) {
            std::cerr << "Errore caricamento elementi orbitali da file EQ1\n";
            return 1;
        }
        
        std::cout << "Elementi orbitali caricati: " << (screener.hasOrbitalElements() ? "SI" : "NO") << "\n";
        
        // Mostra gli elementi convertiti
        const auto& kep = screener.getOrbitalElements();
        std::cout << "\nElementi orbitali (convertiti da equinoziali):\n";
        std::cout << "  Epoca:     MJD " << std::fixed << std::setprecision(1) << kep.epoch_mjd_tdb << "\n";
        std::cout << "  a:         " << std::setprecision(7) << kep.semi_major_axis << " AU\n";
        std::cout << "  e:         " << kep.eccentricity << "\n";
        std::cout << "  i:         " << (kep.inclination * RAD_TO_DEG) << " deg\n";
        std::cout << "  Omega:     " << (kep.longitude_ascending_node * RAD_TO_DEG) << " deg\n";
        std::cout << "  omega:     " << (kep.argument_perihelion * RAD_TO_DEG) << " deg\n";
        std::cout << "  M:         " << (kep.mean_anomaly * RAD_TO_DEG) << " deg\n\n";
        
        // Configura Phase1 - ricerca ampia per test
        Phase1Config config = Phase1Config::survey();  // Configurazione ampia
        config.start_mjd_tdb = start_mjd;
        config.end_mjd_tdb = end_mjd;
        config.max_magnitude = 16.0;  // Solo stelle brillanti per test
        config.corridor_width_deg = 0.1;  // 6 arcmin, ampio
        config.closest_approach_threshold_arcsec = 300.0;  // 5 arcmin threshold
        
        std::cout << "Configurazione Phase1:\n";
        std::cout << "  Intervallo path:     " << config.path_interval_seconds << " sec\n";
        std::cout << "  Larghezza corridor:  " << config.corridor_width_deg << " deg = " 
                  << (config.corridor_width_deg * 3600) << " arcsec\n";
        std::cout << "  Magnitudine limite:  " << config.max_magnitude << "\n";
        std::cout << "  Soglia CA:           " << config.closest_approach_threshold_arcsec << " arcsec\n\n";
        
        std::cout << "Esecuzione screenCandidates...\n";
        
        // Esegui screening
        auto results = screener.screenCandidates(config);
        
        std::cout << "\n=== RISULTATI PHASE1 ===\n\n";
        
        // Statistiche
        std::cout << "Statistiche:\n";
        std::cout << "  Punti nel path:        " << results.num_path_points << "\n";
        std::cout << "  Stelle nel corridor:   " << results.num_stars_in_corridor << "\n";
        std::cout << "  Candidate filtrate:    " << results.num_candidates_filtered << "\n\n";
        
        std::cout << "Tempi di esecuzione:\n";
        std::cout << "  Propagazione:          " << std::fixed << std::setprecision(2) 
                  << results.propagation_time_ms << " ms\n";
        std::cout << "  Query corridor:        " << results.corridor_query_time_ms << " ms\n";
        std::cout << "  Calcolo CA:            " << results.closest_approach_calc_time_ms << " ms\n\n";
        
        // Path dell'asteroide
        std::cout << "=== PATH ASTEROIDE ===\n\n";
        std::cout << std::setw(4) << "Idx" << "  "
                  << std::setw(14) << "MJD" << "  "
                  << std::setw(12) << "RA (deg)" << "  "
                  << std::setw(12) << "Dec (deg)" << "  "
                  << std::setw(12) << "Dist (AU)" << "\n";
        std::cout << std::string(60, '-') << "\n";
        
        // JPL reference per 28 Nov 2025 00:00 UTC
        double jpl_ra = 70.54869;
        double jpl_dec = 20.66593;
        
        for (size_t i = 0; i < results.path.size(); ++i) {
            const auto& pt = results.path[i];
            std::cout << std::setw(4) << i << "  "
                      << std::fixed << std::setprecision(6) << std::setw(14) << pt.mjd_tdb << "  "
                      << std::setprecision(5) << std::setw(12) << pt.ra_deg << "  "
                      << std::setw(12) << pt.dec_deg << "  "
                      << std::setprecision(8) << std::setw(12) << pt.distance_earth_au << "\n";
        }
        
        // Confronto con JPL (primo punto)
        if (!results.path.empty()) {
            const auto& pt = results.path[0];
            double delta_ra = pt.ra_deg - jpl_ra;
            double delta_dec = pt.dec_deg - jpl_dec;
            double sep = std::sqrt(std::pow(delta_ra * std::cos(pt.dec_deg * DEG_TO_RAD), 2) + 
                                   std::pow(delta_dec, 2));
            double sep_arcsec = sep * 3600.0;
            
            std::cout << "\n=== CONFRONTO CON JPL HORIZONS ===\n\n";
            std::cout << "IOccultCalc (primo punto):\n";
            std::cout << "  RA:  " << std::fixed << std::setprecision(5) << pt.ra_deg 
                      << " deg = " << formatRA(pt.ra_deg) << "\n";
            std::cout << "  Dec: " << pt.dec_deg << " deg = " << formatDec(pt.dec_deg) << "\n";
            std::cout << "\nJPL Horizons (28 Nov 2025 00:00 UTC):\n";
            std::cout << "  RA:  " << jpl_ra << " deg = " << formatRA(jpl_ra) << "\n";
            std::cout << "  Dec: " << jpl_dec << " deg = " << formatDec(jpl_dec) << "\n";
            std::cout << "\nDifferenze:\n";
            std::cout << "  Delta RA:     " << std::setprecision(5) << delta_ra << " deg = " 
                      << std::setprecision(1) << (delta_ra * 3600) << " arcsec\n";
            std::cout << "  Delta Dec:    " << std::setprecision(5) << delta_dec << " deg = " 
                      << std::setprecision(1) << (delta_dec * 3600) << " arcsec\n";
            std::cout << "  Separazione:  " << std::setprecision(5) << sep << " deg = "
                      << std::setprecision(1) << sep_arcsec << " arcsec\n";
            
            // Valutazione
            std::cout << "\n=== VALUTAZIONE ACCURATEZZA ===\n\n";
            if (sep_arcsec < 10.0) {
                std::cout << "ECCELLENTE: Errore < 10 arcsec\n";
            } else if (sep_arcsec < 60.0) {
                std::cout << "BUONO: Errore < 1 arcmin\n";
            } else if (sep_arcsec < 300.0) {
                std::cout << "ACCETTABILE: Errore < 5 arcmin\n";
            } else if (sep_arcsec < 3600.0) {
                std::cout << "ATTENZIONE: Errore < 1 deg\n";
            } else {
                std::cout << "ERRORE SIGNIFICATIVO: Errore > 1 deg\n";
            }
        }
        
        // Stelle candidate
        if (!results.candidates.empty()) {
            std::cout << "\n=== STELLE CANDIDATE ===\n\n";
            std::cout << std::setw(20) << "Gaia Source ID" << "  "
                      << std::setw(10) << "RA" << "  "
                      << std::setw(10) << "Dec" << "  "
                      << std::setw(6) << "Gmag" << "  "
                      << std::setw(12) << "CA (arcsec)" << "  "
                      << std::setw(14) << "CA MJD" << "\n";
            std::cout << std::string(82, '-') << "\n";
            
            for (const auto& star : results.candidates) {
                std::cout << std::setw(20) << star.source_id << "  "
                          << std::fixed << std::setprecision(5) << std::setw(10) << star.ra_deg << "  "
                          << std::setw(10) << star.dec_deg << "  "
                          << std::setprecision(2) << std::setw(6) << star.phot_g_mean_mag << "  "
                          << std::setprecision(2) << std::setw(12) << star.closest_approach_arcsec << "  "
                          << std::setprecision(6) << std::setw(14) << star.closest_approach_mjd << "\n";
            }
        } else {
            std::cout << "\n=== NESSUNA STELLA CANDIDATA TROVATA ===\n";
            std::cout << "(Questo può essere normale se non ci sono occultazioni nel periodo)\n";
        }
        
        // Tutte le stelle nel corridor (per debug)
        if (!results.all_stars.empty() && results.all_stars.size() <= 50) {
            std::cout << "\n=== TUTTE LE STELLE NEL CORRIDOR (prime 50) ===\n\n";
            std::cout << std::setw(20) << "Gaia Source ID" << "  "
                      << std::setw(10) << "RA" << "  "
                      << std::setw(10) << "Dec" << "  "
                      << std::setw(6) << "Gmag" << "  "
                      << std::setw(12) << "CA (arcsec)" << "\n";
            std::cout << std::string(70, '-') << "\n";
            
            int count = 0;
            for (const auto& star : results.all_stars) {
                if (count++ >= 50) break;
                std::cout << std::setw(20) << star.source_id << "  "
                          << std::fixed << std::setprecision(5) << std::setw(10) << star.ra_deg << "  "
                          << std::setw(10) << star.dec_deg << "  "
                          << std::setprecision(2) << std::setw(6) << star.phot_g_mean_mag << "  "
                          << std::setprecision(2) << std::setw(12) << star.closest_approach_arcsec << "\n";
            }
        } else if (results.all_stars.size() > 50) {
            std::cout << "\n(Troppe stelle nel corridor per elencarle tutte: " 
                      << results.all_stars.size() << ")\n";
        }
        
    } catch (const std::exception& ex) {
        std::cerr << "\nERRORE: " << ex.what() << "\n";
        return 1;
    }
    
    std::cout << "\n=========================================================\n";
    return 0;
}
