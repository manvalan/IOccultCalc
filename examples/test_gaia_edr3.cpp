/**
 * @file test_gaia_edr3.cpp
 * @brief Test query Gaia EDR3 vs DR3
 */

#include "ioccultcalc/gaia_client.h"
#include <iostream>
#include <iomanip>

using namespace ioccultcalc;

int main() {
    std::cout << "=== TEST GAIA EDR3 vs DR3 ===\n\n";
    
    // Coordinate test: Regione Hyades
    double ra = 67.0;    // Circa 4h 28m
    double dec = 15.0;   // +15°
    double radius = 0.5; // 0.5 gradi
    double mag_limit = 12.0;
    
    std::cout << "Regione test:\n";
    std::cout << "  RA:  " << ra << "°\n";
    std::cout << "  Dec: " << dec << "°\n";
    std::cout << "  Raggio: " << radius << "°\n";
    std::cout << "  Mag limite: " << mag_limit << "\n\n";
    
    try {
        // Test 1: Query con EDR3
        std::cout << "--- Query EDR3 ---\n";
        GaiaClient clientEDR3;
        clientEDR3.setGaiaVersion(GaiaVersion::EDR3);
        clientEDR3.setMaxRows(10);  // Limita a 10 stelle per test
        
        auto starsEDR3 = clientEDR3.queryCone(ra, dec, radius, mag_limit);
        std::cout << "Trovate " << starsEDR3.size() << " stelle con EDR3\n\n";
        
        if (!starsEDR3.empty()) {
            std::cout << "Prime 3 stelle EDR3:\n";
            for (size_t i = 0; i < std::min(size_t(3), starsEDR3.size()); ++i) {
                const auto& s = starsEDR3[i];
                std::cout << "  [" << i+1 << "] ID: " << s.sourceId << "\n";
                std::cout << "      RA:  " << std::fixed << std::setprecision(6) 
                          << s.pos.ra * 180.0 / M_PI << "°\n";
                std::cout << "      Dec: " << s.pos.dec * 180.0 / M_PI << "°\n";
                std::cout << "      Mag: " << std::setprecision(2) << s.phot_g_mean_mag << "\n";
                std::cout << "      PM RA:  " << std::setprecision(3) << s.pmra << " mas/yr\n";
                std::cout << "      PM Dec: " << s.pmdec << " mas/yr\n\n";
            }
        }
        
        // Test 2: Query con DR3
        std::cout << "--- Query DR3 ---\n";
        GaiaClient clientDR3;
        clientDR3.setGaiaVersion(GaiaVersion::DR3);
        clientDR3.setMaxRows(10);
        
        auto starsDR3 = clientDR3.queryCone(ra, dec, radius, mag_limit);
        std::cout << "Trovate " << starsDR3.size() << " stelle con DR3\n\n";
        
        if (!starsDR3.empty()) {
            std::cout << "Prime 3 stelle DR3:\n";
            for (size_t i = 0; i < std::min(size_t(3), starsDR3.size()); ++i) {
                const auto& s = starsDR3[i];
                std::cout << "  [" << i+1 << "] ID: " << s.sourceId << "\n";
                std::cout << "      RA:  " << std::fixed << std::setprecision(6) 
                          << s.pos.ra * 180.0 / M_PI << "°\n";
                std::cout << "      Dec: " << s.pos.dec * 180.0 / M_PI << "°\n";
                std::cout << "      Mag: " << std::setprecision(2) << s.phot_g_mean_mag << "\n";
                std::cout << "      PM RA:  " << std::setprecision(3) << s.pmra << " mas/yr\n";
                std::cout << "      PM Dec: " << s.pmdec << " mas/yr\n\n";
            }
        }
        
        // Confronto
        std::cout << "--- Confronto ---\n";
        std::cout << "Stelle EDR3: " << starsEDR3.size() << "\n";
        std::cout << "Stelle DR3:  " << starsDR3.size() << "\n";
        
        if (starsEDR3.size() == starsDR3.size()) {
            std::cout << "✓ Stesso numero di stelle trovate\n";
        } else {
            std::cout << "⚠️  Numero diverso di stelle (differenze nel catalogo)\n";
        }
        
        std::cout << "\n✓ Test completato con successo!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
