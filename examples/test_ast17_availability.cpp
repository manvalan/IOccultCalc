// Test per verificare la disponibilità di tutti i 17 asteroidi AST17
// @ MJD 55918.807 (Dec 23, 2011) - epoca del test OrbFit

#include <iostream>
#include <iomanip>
#include "ioccultcalc/spice_spk_reader.h"

using namespace ioccultcalc;

int main() {
    try {
        SPICESPKReader spk;
        
        // Carica file asteroidi (codes_300ast contiene i primi 300 asteroidi)
        bool loaded = spk.ensureFileLoaded("codes_300ast_20100725.bsp");
        if (!loaded) {
            std::cerr << "❌ ERRORE: Impossibile caricare codes_300ast_20100725.bsp\n";
            std::cerr << "   Verificare che sia presente in ~/.ioccultcalc/ephemerides/\n";
            return 1;
        }
        std::cout << "✓ codes_300ast_20100725.bsp caricato con CSPICE\n\n";
        
        // MJD 55380.0 (July 1, 2010) - dentro la copertura del file asteroidi
        double mjd = 55380.0;
        double jd = mjd + 2400000.5;
        std::cout << "NOTA: Usando MJD " << mjd << " (July 1, 2010)\n";
        std::cout << "      perché codes_300ast copre fino a luglio 2010\n\n";
        
        // Lista completa AST17 (matching OrbFit)
        struct AsteroidInfo {
            int naifId;
            int number;
            std::string name;
            double GM_AU3day2;  // AU³/day²
        };
        
        AsteroidInfo asteroids[] = {
            {2000001, 1, "Ceres", 62.6284e-12},
            {2000002, 2, "Pallas", 14.3e-12},
            {2000003, 3, "Juno", 1.82e-12},
            {2000004, 4, "Vesta", 17.8e-12},
            {2000006, 6, "Hebe", 0.85e-12},
            {2000007, 7, "Iris", 0.79e-12},
            {2000010, 10, "Hygiea", 5.8e-12},
            {2000015, 15, "Eunomia", 0.38e-12},
            {2000016, 16, "Psyche", 0.36e-12},
            {2000029, 29, "Amphitrite", 0.13e-12},
            {2000052, 52, "Europa", 0.34e-12},
            {2000065, 65, "Cybele", 0.15e-12},
            {2000087, 87, "Sylvia", 0.15e-12},
            {2000088, 88, "Thisbe", 0.16e-12},
            {2000511, 511, "Davida", 0.39e-12},
            {2000704, 704, "Interamnia", 0.34e-12},
            {2134340, 134340, "Pluto", 871.0e-12}
        };
        
        std::cout << "=== TEST DISPONIBILITÀ AST17 @ MJD " << mjd << " ===\n\n";
        std::cout << std::setw(8) << "NAIF ID" << " "
                  << std::setw(6) << "Num" << " "
                  << std::setw(12) << "Nome" << " "
                  << std::setw(15) << "GM (AU³/d²)" << " "
                  << std::setw(10) << "Status" << "\n";
        std::cout << std::string(65, '-') << "\n";
        
        int available = 0;
        int missing = 0;
        
        for (const auto& ast : asteroids) {
            try {
                // Prova a leggere posizione @ epoca test (barycentric frame)
                Vector3D pos = spk.getPosition(ast.naifId, jd, 0);  // SSB center
                double r = pos.magnitude();
                
                std::cout << std::setw(8) << ast.naifId << " "
                          << std::setw(6) << ast.number << " "
                          << std::setw(12) << ast.name << " "
                          << std::scientific << std::setprecision(2)
                          << std::setw(15) << ast.GM_AU3day2 << " "
                          << std::fixed << std::setprecision(3)
                          << "✓ OK (r=" << r << " AU)\n";
                available++;
                
            } catch (const std::exception& e) {
                std::cout << std::setw(8) << ast.naifId << " "
                          << std::setw(6) << ast.number << " "
                          << std::setw(12) << ast.name << " "
                          << std::scientific << std::setprecision(2)
                          << std::setw(15) << ast.GM_AU3day2 << " "
                          << "✗ MISSING\n";
                missing++;
            }
        }
        
        std::cout << "\n=== RIEPILOGO ===\n";
        std::cout << "Asteroidi disponibili: " << available << "/17\n";
        std::cout << "Asteroidi mancanti: " << missing << "/17\n\n";
        
        if (available == 17) {
            std::cout << "✅ PERFETTO! Tutti i 17 asteroidi AST17 sono disponibili!\n";
            std::cout << "   IOccultCalc ora ha la stessa configurazione di OrbFit.\n";
        } else if (available >= 16) {
            std::cout << "⚠️  Quasi tutti disponibili. Possibile mancanza solo per tempi lontani.\n";
        } else {
            std::cout << "❌ ATTENZIONE: Alcuni asteroidi mancano in sb441-n16.bsp\n";
            std::cout << "   Potrebbe essere necessario un file SPK più completo.\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
