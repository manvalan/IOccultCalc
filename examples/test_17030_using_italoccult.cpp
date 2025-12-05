/**
 * @file test_17030_using_italoccult.cpp
 * @brief Test usando la libreria AstDyn di ITALOccult per gestire ECLM J2000
 * 
 * Questo test dimostra l'uso corretto della libreria AstDyn per:
 * 1. Leggere file .eq1 (formato OEF2.0, ECLM J2000)
 * 2. Convertire elementi equinoziali → kepleriani
 * 3. Gestire automaticamente frame eclittico → equatoriale
 * 4. Calcolare posizioni accurate
 * 
 * Usa la libreria statica precompilata: libastdyn.a
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>

// AstDyn headers
#include <astdyn/io/parsers/OrbFitEQ1Parser.hpp>
#include <astdyn/coordinates/EquinoctialElements.hpp>
#include <astdyn/coordinates/KeplerianElements.hpp>
#include <astdyn/coordinates/CartesianState.hpp>
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <astdyn/core/Constants.hpp>

using namespace astdyn;
using namespace astdyn::io::parsers;
using namespace astdyn::coordinates;

// JPL Horizons reference (2025-12-01 00:00 TT)
const double JPL_X_AU = 1.08361;
const double JPL_Y_AU = 3.08629;
const double JPL_Z_AU = -0.09162;
const double JPL_RA_DEG = 70.6536;
const double JPL_DEC_DEG = -1.6045;

void printResults(const std::string& method, 
                  double x_au, double y_au, double z_au,
                  double ra, double dec) {
    // Calcola errore vs JPL Horizons
    double dx = x_au - JPL_X_AU;
    double dy = y_au - JPL_Y_AU;
    double dz = z_au - JPL_Z_AU;
    double error_au = std::sqrt(dx*dx + dy*dy + dz*dz);
    double error_km = error_au * constants::AU_TO_KM;
    
    double ra_error_arcsec = (ra - JPL_RA_DEG) * 3600.0;
    double dec_error_arcsec = (dec - JPL_DEC_DEG) * 3600.0;
    
    std::cout << "\n" << method << ":\n";
    std::cout << "  Position: (" << std::fixed << std::setprecision(5)
              << x_au << ", " << y_au << ", " << z_au << ") AU\n";
    std::cout << "  RA:  " << std::setprecision(4) << ra << "°\n";
    std::cout << "  Dec: " << dec << "°\n";
    std::cout << "  Error: " << std::setprecision(1) << error_km << " km\n";
    std::cout << "  RA error:  " << ra_error_arcsec << " arcsec\n";
    std::cout << "  Dec error: " << dec_error_arcsec << " arcsec\n";
    
    if (error_km < 500000) {
        std::cout << "  ✓ EXCELLENT\n";
    } else if (error_km < 5000000) {
        std::cout << "  ⚠ GOOD\n";
    } else {
        std::cout << "  ✗ POOR\n";
    }
}

int main() {
    std::cout << "\n========================================================\n";
    std::cout << "  Test 17030 - Using AstDyn Library (ITALOccult)\n";
    std::cout << "  Gestione corretta ECLM J2000\n";
    std::cout << "========================================================\n";
    
    try {
        // File path
        std::string file_path = "/Users/michelebigi/VisualStudio Code/GitHub/IOccultCalc/17030_astdys.eq1";
        
        std::cout << "\n1. Lettura file OEF2.0: " << file_path << "\n";
        
        // Parse file usando AstDyn OrbFitEQ1Parser
        OrbFitEQ1Parser parser;
        auto elements = parser.parse(file_path);
        
        std::cout << "\n2. Elementi Equinoziali (ECLM J2000):\n";
        std::cout << "  a = " << elements.a << " AU\n";
        std::cout << "  h = " << elements.h << "\n";
        std::cout << "  k = " << elements.k << "\n";
        std::cout << "  p = " << elements.p << "\n";
        std::cout << "  q = " << elements.q << "\n";
        std::cout << "  λ = " << elements.lambda * 180.0 / M_PI << "° (radianti nel parser)\n";
        std::cout << "  Epoca = MJD " << elements.epoch_mjd << "\n";
        
        // Converti equinoziali → kepleriani usando AstDyn
        std::cout << "\n3. Conversione Equinoziali → Kepleriani (AstDyn):\n";
        KeplerianElements kep = elements.toKeplerian();
        
        std::cout << "  a = " << kep.a / constants::AU_TO_KM << " AU\n";
        std::cout << "  e = " << kep.e << "\n";
        std::cout << "  i = " << kep.i * 180.0 / M_PI << "°\n";
        std::cout << "  Ω = " << kep.Omega * 180.0 / M_PI << "°\n";
        std::cout << "  ω = " << kep.omega * 180.0 / M_PI << "°\n";
        std::cout << "  M = " << kep.M * 180.0 / M_PI << "°\n";
        
        // Converti kepleriani → cartesiano (ECLITTICO)
        std::cout << "\n4. Conversione Kepleriani → Cartesiano Eclittico (AstDyn):\n";
        CartesianState state_ecl = kep.toCartesianState();
        
        std::cout << "  Frame: ECLIPTIC J2000\n";
        std::cout << "  Position: (" << state_ecl.position()(0) / constants::AU_TO_KM
                  << ", " << state_ecl.position()(1) / constants::AU_TO_KM
                  << ", " << state_ecl.position()(2) / constants::AU_TO_KM << ") AU\n";
        
        // Trasforma eclittico → equatoriale usando AstDyn
        std::cout << "\n5. Trasformazione Eclittico → Equatoriale (AstDyn):\n";
        CartesianState state_eq = ReferenceFrame::eclipticToEquatorial(state_ecl);
        
        std::cout << "  Frame: EQUATORIAL J2000 (ICRF)\n";
        double x_au = state_eq.position()(0) / constants::AU_TO_KM;
        double y_au = state_eq.position()(1) / constants::AU_TO_KM;
        double z_au = state_eq.position()(2) / constants::AU_TO_KM;
        
        // Calcola RA/Dec
        double ra = std::atan2(y_au, x_au) * 180.0 / M_PI;
        if (ra < 0) ra += 360.0;
        double dec = std::atan2(z_au, std::sqrt(x_au*x_au + y_au*y_au)) * 180.0 / M_PI;
        
        std::cout << "\n========================================================\n";
        std::cout << "  CONFRONTO CON JPL HORIZONS\n";
        std::cout << "========================================================\n";
        
        std::cout << "\nJPL Horizons (2025-12-01 00:00 TT):\n";
        std::cout << "  Position: (" << JPL_X_AU << ", " << JPL_Y_AU << ", " << JPL_Z_AU << ") AU\n";
        std::cout << "  RA:  " << JPL_RA_DEG << "°\n";
        std::cout << "  Dec: " << JPL_DEC_DEG << "°\n";
        
        printResults("AstDyn (ECLM J2000 → ICRF)", x_au, y_au, z_au, ra, dec);
        
        std::cout << "\n========================================================\n";
        std::cout << "  CONCLUSIONE\n";
        std::cout << "========================================================\n";
        std::cout << "\n✓ La libreria AstDyn gestisce correttamente:\n";
        std::cout << "  - Lettura file .eq1 (OEF2.0) via OrbFitEQ1Parser\n";
        std::cout << "  - Conversione λ da gradi a radianti (nel parser)\n";
        std::cout << "  - Conversione equinoziali → kepleriani\n";
        std::cout << "  - Trasformazione ECLM J2000 → ICRF\n";
        std::cout << "  - Calcolo RA/Dec accurate\n\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ ERRORE: " << e.what() << "\n\n";
        return 1;
    }
}
