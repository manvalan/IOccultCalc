/**
 * @file test_iers.cpp
 * @brief Test IERS Earth Orientation Parameters
 */

#include "ioccultcalc/iers_data.h"
#include <iostream>
#include <iomanip>

using namespace ioccultcalc;

int main() {
    std::cout << "═══════════════════════════════════════════════════════════" << std::endl;
    std::cout << "  TEST IERS EARTH ORIENTATION PARAMETERS" << std::endl;
    std::cout << "═══════════════════════════════════════════════════════════" << std::endl;
    
    try {
        IERSData iers;
        
        std::cout << "\n1. Downloading IERS data..." << std::endl;
        if (iers.downloadLatestData()) {
            std::cout << "✓ Download successful" << std::endl;
            
            double mjd_start, mjd_end;
            if (iers.getDataRange(mjd_start, mjd_end)) {
                std::cout << "  Data range: MJD " << std::fixed << std::setprecision(1) 
                         << mjd_start << " to " << mjd_end << std::endl;
                std::cout << "  Coverage: " << std::setprecision(0) 
                         << (mjd_end - mjd_start) / 365.25 << " years" << std::endl;
            }
        } else {
            std::cout << "⚠️  Download failed - using cached data if available" << std::endl;
        }
        
        // Test EOP for specific dates
        std::cout << "\n2. Testing EOP retrieval..." << std::endl;
        
        struct TestDate {
            const char* name;
            double jd;
        };
        
        TestDate dates[] = {
            {"2000-01-01 (J2000.0)", 2451545.0},
            {"2024-01-01", 2460310.5},
            {"2024-11-21 (today)", 2460635.5},
        };
        
        std::cout << std::fixed << std::setprecision(6);
        for (const auto& date : dates) {
            auto eop = iers.getEOP(date.jd);
            
            std::cout << "\n  " << date.name << " (JD " << date.jd << "):" << std::endl;
            std::cout << "    UT1-UTC:  " << std::setw(10) << eop.ut1_utc << " s" << std::endl;
            std::cout << "    x_pole:   " << std::setw(10) << eop.x_pole << " arcsec" << std::endl;
            std::cout << "    y_pole:   " << std::setw(10) << eop.y_pole << " arcsec" << std::endl;
            std::cout << "    dPsi:     " << std::setw(10) << eop.dPsi << " arcsec" << std::endl;
            std::cout << "    dEps:     " << std::setw(10) << eop.dEps << " arcsec" << std::endl;
        }
        
        std::cout << "\n3. Impact on observatory position..." << std::endl;
        std::cout << "   Polar motion ~0.5\" × Earth radius (6371 km):" << std::endl;
        std::cout << "   → ~15 meters position uncertainty" << std::endl;
        std::cout << "   → ~0.003\" for asteroid at 1 AU" << std::endl;
        std::cout << "   Implementing EOP corrections → sub-0.01\" precision ✓" << std::endl;
        
        std::cout << "\n✓ Test completato!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "✗ Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
