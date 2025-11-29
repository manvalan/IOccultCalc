/**
 * @file test_sierks_quick.cpp
 * @brief Test veloce per verificare posizione Sierks e Ceres
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/astdys_client.h"

using namespace ioccultcalc;

void testAsteroid(const std::string& name, const std::string& designation,
                   double jd_target, double ref_ra = -999, double ref_dec = -999) {
    std::cout << "=== " << name << " (" << designation << ") ===\n";
    
    AstDysClient client;
    EquinoctialElements elements;
    
    try {
        elements = client.getElements(designation);
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << "\n";
        return;
    }
    
    std::cout << "Elementi orbitali:\n";
    std::cout << "  a      = " << std::fixed << std::setprecision(8) << elements.a << " AU\n";
    std::cout << "  lambda = " << elements.lambda << " rad (" 
              << (elements.lambda * 180.0/M_PI) << " deg)\n";
    std::cout << "  epoch  = JD " << std::setprecision(1) << elements.epoch.jd << "\n";
    std::cout << "  Delta epoca: " << (jd_target - elements.epoch.jd) << " giorni\n\n";
    
    Ephemeris eph(elements);
    eph.setOptions(EphemerisOptions::fast());
    
    EphemerisData data = eph.compute(JulianDate(jd_target));
    
    std::cout << "Posizione geocentrica (Keplero puro):\n";
    std::cout << "  RA  = " << std::setprecision(4) << data.geocentricPos.ra << " deg\n";
    std::cout << "  Dec = " << std::showpos << data.geocentricPos.dec << std::noshowpos << " deg\n";
    
    if (ref_ra > -900) {
        double dra = (data.geocentricPos.ra - ref_ra) * 3600.0;
        double ddec = (data.geocentricPos.dec - ref_dec) * 3600.0;
        std::cout << "\nRiferimento: RA=" << ref_ra << "°, Dec=" << std::showpos << ref_dec << std::noshowpos << "°\n";
        std::cout << "Errore RA:  " << std::showpos << dra << std::noshowpos << " arcsec\n";
        std::cout << "Errore Dec: " << std::showpos << ddec << std::noshowpos << " arcsec\n";
    }
    
    std::cout << "\n";
}

int main() {
    double jd_oggi = 2460646.5;  // 29 Nov 2025
    
    std::cout << "Test posizioni asteroidi - JD " << std::fixed << std::setprecision(1) << jd_oggi << "\n";
    std::cout << "(29 Nov 2025)\n\n";
    
    // Ceres - JPL: RA≈7°, Dec≈-9.5° il 29 Nov 2025
    testAsteroid("Ceres", "1", jd_oggi, 7.0, -9.5);
    
    // Sierks - Stella IOTA: RA=73.416°, Dec=+20.332°
    // Occultazione prevista 28 Nov 2025 00:35 UT
    double jd_sierks = 2460645.524;  // 28 Nov 2025 00:35 UT
    testAsteroid("Sierks", "17030", jd_sierks, 73.416, 20.332);
    
    return 0;
}
