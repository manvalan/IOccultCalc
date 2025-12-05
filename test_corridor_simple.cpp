/**
 * @file test_corridor_simple.cpp
 * @brief Test semplice query corridor UnifiedGaiaCatalog
 * @date 4 Dicembre 2025
 * 
 * Test: Stelle tra RA 4h30m DEC +22°45' e RA 4h35m DEC +22°35'
 * Magnitudine: <= 12.0
 */

#include "ioc_gaialib/unified_gaia_catalog.h"
#include "ioc_gaialib/types.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>

int main() {
    try {
        std::cout << "========================================\n";
        std::cout << "  Test Query Corridor - UnifiedGaiaCatalog\n";
        std::cout << "========================================\n\n";
        
        // 1. Inizializza catalogo
        std::cout << "1. Inizializzazione catalogo...\n";
        std::string home = std::getenv("HOME");
        std::string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
        
        std::string json_config = R"({
            "catalog_type": "multifile_v2",
            "multifile_directory": ")" + catalog_path + R"(",
            "cache_size_mb": 512,
            "prefetch_enabled": true,
            "log_level": "info"
        })";
        
        std::cout << "   Path: " << catalog_path << "\n";
        
        if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) {
            std::cerr << "ERRORE: Impossibile inizializzare catalogo\n";
            return 1;
        }
        
        std::cout << "   ✓ Catalogo inizializzato\n\n";
        
        // 2. Ottieni istanza
        auto& catalog = ioc::gaia::UnifiedGaiaCatalog::getInstance();
        
        // 3. Definisci corridor path
        // RA 4h30m = 4.5 * 15 = 67.5 gradi
        // RA 4h35m = 4.583333 * 15 = 68.75 gradi
        // DEC +22°45' = 22.75 gradi
        // DEC +22°35' = 22.583333 gradi
        
        std::cout << "2. Creazione corridor path...\n";
        std::vector<ioc::gaia::CelestialPoint> path;
        
        ioc::gaia::CelestialPoint p1, p2;
        p1.ra = 67.5;      // 4h30m
        p1.dec = 22.75;    // +22°45'
        
        p2.ra = 68.75;     // 4h35m
        p2.dec = 22.583333; // +22°35'
        
        path.push_back(p1);
        path.push_back(p2);
        
        std::cout << "   P1: RA=" << std::fixed << std::setprecision(2) 
                  << p1.ra << "° DEC=" << p1.dec << "°\n";
        std::cout << "   P2: RA=" << p2.ra << "° DEC=" << p2.dec << "°\n";
        
        // 4. Configura parametri query
        ioc::gaia::CorridorQueryParams params;
        params.path = path;
        params.width = 0.5;  // ±0.5 gradi = ±30 arcmin (molto largo)
        params.max_magnitude = 12.0;
        params.min_parallax = -1.0;  // no limit
        params.max_results = 0;      // no limit
        
        std::cout << "   Width: ±" << params.width << "° (±" 
                  << (params.width * 60.0) << " arcmin)\n";
        std::cout << "   Max magnitude: " << params.max_magnitude << "\n\n";
        
        // 5. Esegui query
        std::cout << "3. Query corridor...\n";
        auto start = std::chrono::high_resolution_clock::now();
        
        std::vector<ioc::gaia::GaiaStar> stars = catalog.queryCorridor(params);
        
        auto end = std::chrono::high_resolution_clock::now();
        double duration_ms = std::chrono::duration<double, std::milli>(end - start).count();
        
        std::cout << "   ✓ Query completata in " << std::fixed << std::setprecision(1) 
                  << duration_ms << " ms\n";
        std::cout << "   Stelle trovate: " << stars.size() << "\n\n";
        
        // 6. Stampa risultati
        if (stars.empty()) {
            std::cout << "⚠ Nessuna stella trovata nel corridor.\n";
        } else {
            std::cout << "========================================\n";
            std::cout << "  STELLE TROVATE\n";
            std::cout << "========================================\n\n";
            
            std::cout << std::setw(5) << "N" 
                      << std::setw(10) << "RA (°)"
                      << std::setw(10) << "DEC (°)"
                      << std::setw(8) << "Mag"
                      << "  Gaia DR3 Source ID\n";
            std::cout << std::string(60, '-') << "\n";
            
            for (size_t i = 0; i < stars.size(); ++i) {
                const auto& star = stars[i];
                
                // Converti RA in ore
                double ra_hours = star.ra / 15.0;
                int ra_h = static_cast<int>(ra_hours);
                double ra_m = (ra_hours - ra_h) * 60.0;
                
                // Converti DEC in gradi/minuti
                int dec_d = static_cast<int>(star.dec);
                double dec_m = std::abs((star.dec - dec_d) * 60.0);
                
                std::cout << std::setw(5) << (i+1)
                          << std::setw(10) << std::fixed << std::setprecision(4) << star.ra
                          << std::setw(10) << std::setprecision(4) << star.dec
                          << std::setw(8) << std::setprecision(2) << star.phot_g_mean_mag
                          << "  " << star.source_id
                          << " (RA: " << ra_h << "h" << std::setprecision(1) << ra_m << "m"
                          << " DEC: " << (star.dec >= 0 ? "+" : "") << dec_d << "°" 
                          << std::setprecision(1) << dec_m << "')\n";
            }
            
            std::cout << "\n========================================\n";
            std::cout << "Totale: " << stars.size() << " stelle\n";
            std::cout << "========================================\n";
        }
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << "\n";
        return 1;
    }
}
