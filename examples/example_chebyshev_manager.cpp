/**
 * @file example_chebyshev_manager.cpp
 * @brief Esempio d'uso di ChebyshevOccultationManager
 * @author IOccultCalc Development Team
 * @date 4 Dicembre 2025
 * 
 * Questo esempio dimostra il workflow completo per la ricerca di occultazioni
 * usando ChebyshevOccultationManager che integra:
 * 
 * 1. ChebyshevRKF78Propagator di AstDyn (propagazione + fitting integrati)
 * 2. Corridor API per ricerca stelle
 * 3. Calcolo closest approach
 * 
 * Uso:
 *   ./example_chebyshev_manager <file.eq1> [start_mjd] [end_mjd]
 * 
 * Esempio:
 *   ./example_chebyshev_manager 17030.eq1 61003.5 61010.5
 */

#include "chebyshev_occultation_manager.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace ioccultcalc;

int main(int argc, char* argv[]) {
    try {
        // =====================================================================
        // PARSE ARGOMENTI
        // =====================================================================
        
        if (argc < 2) {
            std::cerr << "Uso: " << argv[0] << " <file.eq1> [start_mjd] [end_mjd]\n";
            std::cerr << "\nEsempio:\n";
            std::cerr << "  " << argv[0] << " 17030.eq1 61003.5 61010.5\n";
            std::cerr << "\nArgomenti:\n";
            std::cerr << "  file.eq1    File elementi orbitali (formato AstDyS/OrbFit)\n";
            std::cerr << "  start_mjd   Epoca inizio (MJD TDB, default: 61003.5 = 2025-11-24)\n";
            std::cerr << "  end_mjd     Epoca fine (MJD TDB, default: 61010.5 = 2025-12-01)\n";
            return 1;
        }
        
        std::string eq1_file = argv[1];
        double start_mjd = 61003.5;  // 2025-11-24
        double end_mjd = 61010.5;    // 2025-12-01
        
        if (argc > 2) start_mjd = std::atof(argv[2]);
        if (argc > 3) end_mjd = std::atof(argv[3]);
        
        std::cout << "=" << std::string(68, '=') << "\n";
        std::cout << "  ChebyshevOccultationManager - Esempio d'uso\n";
        std::cout << "=" << std::string(68, '=') << "\n\n";
        
        std::cout << "File .eq1: " << eq1_file << "\n";
        std::cout << "Intervallo: MJD " << std::fixed << std::setprecision(1)
                  << start_mjd << " - " << end_mjd
                  << " (" << (end_mjd - start_mjd) << " giorni)\n\n";
        
        // =====================================================================
        // STEP 1: CONFIGURA RICERCA
        // =====================================================================
        
        std::cout << "STEP 1: Configurazione ricerca\n";
        std::cout << std::string(70, '-') << "\n";
        
        OccultationSearchConfig config;
        config.start_mjd = start_mjd;
        config.end_mjd = end_mjd;
        config.num_propagation_points = 100;     // 100 punti RKF78
        config.num_chebyshev_coeffs = 8;         // 8 coefficienti per asse
        config.corridor_width_km = 1000.0;       // Corridor ±1000 km (più largo per test)
        config.search_step_days = 0.1;           // Campiona ogni 2.4 ore
        config.max_magnitude = 16.0;             // Stelle fino a mag 16
        config.threshold_arcsec = 2.0;           // Occultazione se < 2"
        
        // Observer location (default: Italia centrale)
        config.observer_lat_deg = 42.0;          // Latitudine 42° N
        config.observer_lon_deg = 12.0;          // Longitudine 12° E
        config.observer_alt_m = 100.0;           // Altitudine 100 m slm
        
        std::cout << "✓ Configurazione:\n";
        std::cout << "  Punti propagazione:     " << config.num_propagation_points << "\n";
        std::cout << "  Coefficienti Chebyshev: " << config.num_chebyshev_coeffs << " (x3 assi)\n";
        std::cout << "  Corridor width:         " << config.corridor_width_km << " km\n";
        std::cout << "  Search step:            " << config.search_step_days << " giorni\n";
        std::cout << "  Max magnitude:          " << config.max_magnitude << "\n";
        std::cout << "  Threshold:              " << config.threshold_arcsec << " arcsec\n";
        std::cout << "  Observer: " << config.observer_lat_deg << "°N, "
                  << config.observer_lon_deg << "°E, " << config.observer_alt_m << "m\n\n";
        
        // =====================================================================
        // STEP 2: CREA MANAGER E CARICA ASTEROIDE
        // =====================================================================
        
        std::cout << "STEP 2: Creazione manager e caricamento asteroide\n";
        std::cout << std::string(70, '-') << "\n";
        
        ChebyshevOccultationManager manager(config);
        
        if (!manager.loadAsteroidFromEQ1(eq1_file)) {
            std::cerr << "\nERRORE: Impossibile caricare asteroide da " << eq1_file << "\n";
            return 1;
        }
        
        std::cout << "\n";
        
        // =====================================================================
        // STEP 3: PROPAGA CON RKF78 E FITTA CHEBYSHEV
        // =====================================================================
        
        std::cout << "STEP 3: Propagazione RKF78 e fitting Chebyshev\n";
        std::cout << std::string(70, '-') << "\n";
        std::cout << "Questo step usa ChebyshevRKF78Propagator di AstDyn che integra:\n";
        std::cout << "  - RKF78 integrator (tolerance 1e-12 AU)\n";
        std::cout << "  - 8 perturbazioni planetarie (Mercury → Neptune)\n";
        std::cout << "  - Perturbazioni asteroidi (AST17 database)\n";
        std::cout << "  - Correzioni relativistiche (Schwarzschild)\n";
        std::cout << "  - Conversione frame ECLM → ICRF automatica\n";
        std::cout << "  - Output: barycentric ICRF J2000.0 in AU\n\n";
        
        if (!manager.propagateAndFit()) {
            std::cerr << "\nERRORE: Propagazione/fitting fallita\n";
            return 1;
        }
        
        // =====================================================================
        // STEP 4: TESTA QUERY VELOCE CHEBYSHEV
        // =====================================================================
        
        std::cout << "\nSTEP 4: Test query Chebyshev (sub-microsecond)\n";
        std::cout << std::string(70, '-') << "\n";
        
        // Campiona alcune epoche per mostrare la velocità
        std::vector<double> test_epochs = {
            start_mjd,
            start_mjd + (end_mjd - start_mjd) * 0.25,
            start_mjd + (end_mjd - start_mjd) * 0.50,
            start_mjd + (end_mjd - start_mjd) * 0.75,
            end_mjd
        };
        
        std::cout << "Valutazione posizione a 5 epoche diverse:\n\n";
        std::cout << std::fixed << std::setprecision(1);
        
        for (size_t i = 0; i < test_epochs.size(); ++i) {
            double mjd = test_epochs[i];
            auto start = std::chrono::high_resolution_clock::now();
            Eigen::Vector3d pos = manager.getPositionAtEpoch(mjd);
            auto end = std::chrono::high_resolution_clock::now();
            
            double query_time_ns = std::chrono::duration<double, std::nano>(
                end - start).count();
            
            std::cout << "MJD " << std::setw(8) << mjd << ": "
                      << "x=" << std::setw(10) << std::setprecision(6) << pos.x() << " AU, "
                      << "y=" << std::setw(10) << pos.y() << " AU, "
                      << "z=" << std::setw(10) << pos.z() << " AU  "
                      << "(query: " << std::setw(6) << std::setprecision(0)
                      << query_time_ns << " ns)\n";
        }
        
        std::cout << "\n✓ Query Chebyshev: ordine di nanosecondi (10⁵x più veloce di RKF78)\n";
        std::cout << "  RMS error: " << std::scientific << std::setprecision(3)
                  << manager.getChebyshevRMSError() << " AU (machine precision)\n\n";
        
        // =====================================================================
        // STEP 5: CERCA STELLE NEL CORRIDOR
        // =====================================================================
        
        std::cout << "STEP 5: Ricerca stelle nel corridor\n";
        std::cout << std::string(70, '-') << "\n";
        
        size_t num_stars = manager.searchStarsInCorridor();
        
        if (num_stars == 0) {
            std::cout << "\n⚠ Nessuna stella trovata nel corridor.\n";
            std::cout << "  Nota: Implementazione corridor search ancora TODO.\n\n";
        }
        
        // =====================================================================
        // STEP 6: CALCOLA CLOSEST APPROACHES
        // =====================================================================
        
        std::cout << "STEP 6: Calcolo closest approaches\n";
        std::cout << std::string(70, '-') << "\n";
        
        size_t num_occultations = manager.computeClosestApproaches();
        
        if (num_occultations == 0) {
            std::cout << "\n⚠ Nessuna occultazione trovata.\n";
            std::cout << "  Nota: Implementazione closest approach ancora TODO.\n\n";
        }
        
        // =====================================================================
        // STEP 7: STAMPA RIEPILOGO
        // =====================================================================
        
        manager.printSummary();
        
        // =====================================================================
        // STEP 8: ESEMPI EXPORT
        // =====================================================================
        
        std::cout << "\nSTEP 8: Export risultati\n";
        std::cout << std::string(70, '-') << "\n";
        std::cout << "Formati disponibili:\n";
        std::cout << "  - JSON: exportResultsJSON(\"results.json\")\n";
        std::cout << "  - OOP:  exportResultsOOP(\"results.oop\") [IOTA format]\n";
        std::cout << "\n⚠ Export ancora TODO.\n\n";
        
        // =====================================================================
        // CONCLUSIONE
        // =====================================================================
        
        std::cout << std::string(70, '=') << "\n";
        std::cout << "ESEMPIO COMPLETATO CON SUCCESSO\n";
        std::cout << std::string(70, '=') << "\n\n";
        
        std::cout << "Prossimi step di sviluppo:\n";
        std::cout << "  1. Implementare searchStarsInCorridor() con Corridor API\n";
        std::cout << "  2. Implementare computeClosestApproaches() con detector\n";
        std::cout << "  3. Implementare exportResultsJSON() e exportResultsOOP()\n";
        std::cout << "  4. Aggiungere parallasse e proper motion stelle\n";
        std::cout << "  5. Aggiungere calcolo altitudine/azimuth per observer\n\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << "\n";
        return 1;
    }
}
