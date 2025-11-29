#include <iostream>
#include <iomanip>
#include "ioccultcalc/spice_spk_reader.h"
#include "ioccultcalc/types.h"

using namespace ioccultcalc;

int main() {
    std::cout << "=== TEST SPK CERES ===\n\n";
    
    SPICESPKReader reader;
    
    // Prova a caricare file asteroidi (sb441-n16 copre 2025!)
    if (!reader.ensureFileLoaded("sb441-n16.bsp")) {
        std::cerr << "Errore caricamento SPK asteroidi sb441\n";
        // Prova anche codes_300ast
        if (!reader.ensureFileLoaded("codes_300ast_20100725.bsp")) {
            std::cerr << "Errore caricamento SPK asteroidi\n";
            return 1;
        }
    }
    
    // NAIF ID Ceres = 2000001
    double jd = 2460646.5;  // 29 Nov 2025
    
    std::cout << "Data: JD " << std::fixed << std::setprecision(1) << jd << "\n";
    
    try {
        Vector3D pos = reader.getPosition(2000001, jd, 10);  // 10 = Sun
        std::cout << "\nPosizione Ceres da SPK (eliocentrica):\n";
        std::cout << "  X = " << std::setprecision(6) << pos.x << " AU\n";
        std::cout << "  Y = " << pos.y << " AU\n";
        std::cout << "  Z = " << pos.z << " AU\n";
        
        // Riferimento JPL eclittico
        double x_jpl = 1.8630;
        double y_jpl = -2.1630;
        double z_jpl = -0.4890;
        
        std::cout << "\nRiferimento JPL (eclittico):\n";
        std::cout << "  X = " << x_jpl << " AU\n";
        std::cout << "  Y = " << y_jpl << " AU\n";
        std::cout << "  Z = " << z_jpl << " AU\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Errore: " << e.what() << std::endl;
    }
    
    return 0;
}
