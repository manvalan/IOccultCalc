/**
 * Test caricamento file .eq1 e .rwo locali
 * Verifica che AstDysClient carichi correttamente elementi e osservazioni da file
 */

#include "ioccultcalc/astdys_client.h"
#include <iostream>
#include <iomanip>

using namespace ioccultcalc;

int main() {
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "   TEST: Caricamento File Locali AstDyS (.eq1 + .rwo)\n";
    std::cout << "═══════════════════════════════════════════════════════════\n\n";
    
    try {
        AstDysClient client;
        
        // Configura directory locali
        std::string dataDir = "test_astdys_download";
        std::cout << "📂 Configurazione directory:\n";
        std::cout << "   .eq1: " << dataDir << "\n";
        std::cout << "   .rwo: " << dataDir << "\n\n";
        
        client.setLocalEQ1Directory(dataDir);
        client.setLocalRWODirectory(dataDir);
        
        // Test asteroidi disponibili
        std::vector<std::string> asteroids = {"433", "10", "203"};
        
        for (const auto& ast : asteroids) {
            std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
            std::cout << "Asteroide (" << ast << ")\n";
            std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n";
            
            // 1. Carica elementi orbitali (.eq1)
            std::cout << "1️⃣  Caricamento elementi orbitali (.eq1)...\n";
            EquinoctialElements elem = client.getElements(ast);
            
            std::cout << "   ✓ Elementi caricati:\n";
            std::cout << "     a     = " << std::fixed << std::setprecision(6) << elem.a << " AU\n";
            std::cout << "     h     = " << elem.h << "\n";
            std::cout << "     k     = " << elem.k << "\n";
            std::cout << "     p     = " << elem.p << "\n";
            std::cout << "     q     = " << elem.q << "\n";
            std::cout << "     λ     = " << (elem.lambda * 180.0 / M_PI) << "°\n";
            std::cout << "     Epoca = MJD " << (elem.epoch.jd - 2400000.5) << "\n";
            std::cout << "     H     = " << elem.H << " mag\n";
            std::cout << "     G     = " << elem.G << "\n\n";
            
            // 2. Carica osservazioni (.rwo)
            std::cout << "2️⃣  Caricamento osservazioni (.rwo)...\n";
            auto observations = client.getObservations(ast);
            
            std::cout << "   ✓ Osservazioni caricate: " << observations.size() << " righe\n";
            if (!observations.empty()) {
                std::cout << "     Prima osservazione:\n";
                std::cout << "     " << observations[0].substr(0, 80) << "...\n";
                std::cout << "     Ultima osservazione:\n";
                std::cout << "     " << observations.back().substr(0, 80) << "...\n";
            }
            std::cout << "\n";
        }
        
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        std::cout << "✅ TEST COMPLETATO CON SUCCESSO\n";
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << "\n\n";
        return 1;
    }
}
