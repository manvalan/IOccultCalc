/**
 * @file test_keplero_debug.cpp
 * @brief Test con elementi MPC (Kepleriani osculanti)
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/orbital_elements.h"
#include <nlohmann/json.hpp>

using namespace ioccultcalc;
using json = nlohmann::json;

// Legge elementi Kepleriani da MPC JSON per un asteroide specifico
OrbitalElements readMPCElements(const std::string& number) {
    std::string path = std::string(getenv("HOME")) + "/.ioccultcalc/data/all_numbered_asteroids.json";
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Cannot open MPC file: " + path);
    }
    
    json data = json::parse(file);
    
    std::string target = "(" + number + ")";
    
    for (const auto& item : data) {
        if (item.contains("Number") && item["Number"].get<std::string>() == target) {
            OrbitalElements elem;
            
            elem.a = item["a"].get<double>();
            elem.e = item["e"].get<double>();
            elem.i = item["i"].get<double>() * M_PI / 180.0;  // gradi -> rad
            elem.Omega = item["Node"].get<double>() * M_PI / 180.0;
            elem.omega = item["Peri"].get<double>() * M_PI / 180.0;
            elem.M = item["M"].get<double>() * M_PI / 180.0;
            elem.epoch = JulianDate(item["Epoch"].get<double>());
            
            if (item.contains("H")) elem.H = item["H"].get<double>();
            if (item.contains("G")) elem.G = item["G"].get<double>();
            if (item.contains("Name")) elem.name = item["Name"].get<std::string>();
            elem.designation = number;
            
            return elem;
        }
    }
    
    throw std::runtime_error("Asteroid not found: " + number);
}

int main() {
    std::cout << "=== TEST CON ELEMENTI MPC (Kepleriani osculanti) ===\n\n";
    
    // Leggi Ceres da MPC
    OrbitalElements kepCeres = readMPCElements("1");
    
    std::cout << "Elementi MPC Ceres:\n";
    std::cout << "  a     = " << std::fixed << std::setprecision(7) << kepCeres.a << " AU\n";
    std::cout << "  e     = " << kepCeres.e << "\n";
    std::cout << "  i     = " << (kepCeres.i * 180/M_PI) << "°\n";
    std::cout << "  Ω     = " << (kepCeres.Omega * 180/M_PI) << "°\n";
    std::cout << "  ω     = " << (kepCeres.omega * 180/M_PI) << "°\n";
    std::cout << "  M     = " << (kepCeres.M * 180/M_PI) << "°\n";
    std::cout << "  Epoca = JD " << std::setprecision(1) << kepCeres.epoch.jd << "\n\n";
    
    // Converti in equinoziali
    EquinoctialElements eqCeres = kepCeres.toEquinoctial();
    
    std::cout << "Convertito in equinoziali:\n";
    std::cout << "  a = " << std::setprecision(7) << eqCeres.a << "\n";
    std::cout << "  h = " << eqCeres.h << " (e·sin(ω̃))\n";
    std::cout << "  k = " << eqCeres.k << " (e·cos(ω̃))\n";
    std::cout << "  p = " << eqCeres.p << "\n";
    std::cout << "  q = " << eqCeres.q << "\n";
    std::cout << "  λ = " << eqCeres.lambda << " rad\n\n";
    
    // Target: 29 Nov 2025
    JulianDate jd_target(2460646.5);
    double dt = jd_target.jd - kepCeres.epoch.jd;
    std::cout << "Target: JD " << std::setprecision(1) << jd_target.jd << " (Δt = " << dt << " giorni)\n\n";
    
    // Calcola con elementi MPC convertiti
    Ephemeris eph(eqCeres);
    eph.setOptions(EphemerisOptions::fast());
    auto data = eph.compute(jd_target);
    
    std::cout << "Posizione geocentrica (MPC → Equinoziali → Keplero):\n";
    std::cout << "  RA  = " << std::setprecision(4) << data.geocentricPos.ra << "°\n";
    std::cout << "  Dec = " << std::showpos << data.geocentricPos.dec << std::noshowpos << "°\n\n";
    
    // Riferimento JPL
    std::cout << "=== RIFERIMENTO JPL ===\n";
    std::cout << "  Ceres 29 Nov 2025: RA ≈ 6.8°, Dec ≈ -9.6°\n\n";
    
    double err_ra = (data.geocentricPos.ra - 6.8) * 3600.0;
    double err_dec = (data.geocentricPos.dec - (-9.6)) * 3600.0;
    std::cout << "  Errore: ΔRA = " << std::showpos << err_ra << "\", ΔDec = " << err_dec << std::noshowpos << "\"\n";
    
    // ===== DEBUG: Posizione Terra e Asteroide =====
    std::cout << "\n=== DEBUG POSIZIONI ===\n";
    Vector3D earthPos = Ephemeris::getEarthPosition(jd_target);
    std::cout << "Terra eliocentrica (equatoriale J2000):\n";
    std::cout << "  X = " << std::setprecision(6) << earthPos.x << " AU\n";
    std::cout << "  Y = " << earthPos.y << " AU\n";
    std::cout << "  Z = " << earthPos.z << " AU\n";
    
    std::cout << "Asteroide eliocentrico (da ephemeris):\n";
    std::cout << "  X = " << data.heliocentricPos.x << " AU\n";
    std::cout << "  Y = " << data.heliocentricPos.y << " AU\n";
    std::cout << "  Z = " << data.heliocentricPos.z << " AU\n";
    
    // JPL Horizons Ceres 29 Nov 2025:
    // Heliocentric (ecliptic): X≈1.86, Y≈-2.16, Z≈-0.49 AU
    // Ma noi usiamo equatoriale, quindi Z è diverso
    std::cout << "\nRiferimento JPL Ceres helio (ECLITTICO): X≈1.86, Y≈-2.16, Z≈-0.49 AU\n";
    
    // Converti nostra posizione da equatoriale a eclittico per confronto
    constexpr double EPS = 23.4392911 * M_PI / 180.0;
    double x_ecl = data.heliocentricPos.x;
    double y_ecl = data.heliocentricPos.y * cos(EPS) + data.heliocentricPos.z * sin(EPS);
    double z_ecl = -data.heliocentricPos.y * sin(EPS) + data.heliocentricPos.z * cos(EPS);
    std::cout << "Nostra posizione (convertita in ECLITTICO): X=" << x_ecl 
              << ", Y=" << y_ecl << ", Z=" << z_ecl << " AU\n";
    
    // ===== TEST PROPAGAZIONE IN AVANTI =====
    std::cout << "\n=== TEST PROPAGAZIONE IN AVANTI (dall'epoca) ===\n";
    JulianDate jd_futuro(2461100.5);  // 100 giorni dopo l'epoca
    double dt_fut = jd_futuro.jd - kepCeres.epoch.jd;
    std::cout << "Target: JD " << std::setprecision(1) << jd_futuro.jd << " (Δt = +" << dt_fut << " giorni)\n";
    
    auto data_fut = eph.compute(jd_futuro);
    std::cout << "  RA  = " << std::setprecision(4) << data_fut.geocentricPos.ra << "°\n";
    std::cout << "  Dec = " << std::showpos << data_fut.geocentricPos.dec << std::noshowpos << "°\n";
    std::cout << "  (Non abbiamo riferimento JPL per questa data futura)\n";
    
    // ===== SIERKS =====
    std::cout << "\n=== TEST SIERKS ===\n";
    OrbitalElements kepSierks = readMPCElements("17030");
    
    std::cout << "Elementi MPC Sierks:\n";
    std::cout << "  a     = " << std::setprecision(7) << kepSierks.a << " AU\n";
    std::cout << "  e     = " << kepSierks.e << "\n";
    std::cout << "  i     = " << (kepSierks.i * 180/M_PI) << "°\n";
    std::cout << "  Ω     = " << (kepSierks.Omega * 180/M_PI) << "°\n";
    std::cout << "  ω     = " << (kepSierks.omega * 180/M_PI) << "°\n";
    std::cout << "  M     = " << (kepSierks.M * 180/M_PI) << "°\n";
    std::cout << "  Epoca = JD " << std::setprecision(1) << kepSierks.epoch.jd << "\n\n";
    
    // Converti e calcola
    EquinoctialElements eqSierks = kepSierks.toEquinoctial();
    
    JulianDate jd_sierks(2460645.524);  // 28 Nov 2025 00:35 UT
    
    Ephemeris ephS(eqSierks);
    ephS.setOptions(EphemerisOptions::fast());
    auto dataS = ephS.compute(jd_sierks);
    
    std::cout << "Posizione geocentrica Sierks (28 Nov 2025 00:35 UT):\n";
    std::cout << "  RA  = " << std::setprecision(4) << dataS.geocentricPos.ra << "°\n";
    std::cout << "  Dec = " << std::showpos << dataS.geocentricPos.dec << std::noshowpos << "°\n\n";
    
    // Stella IOTA
    double star_ra = 73.416106;
    double star_dec = 20.331662;
    
    double sep_ra = (dataS.geocentricPos.ra - star_ra) * 3600.0;
    double sep_dec = (dataS.geocentricPos.dec - star_dec) * 3600.0;
    double sep = std::sqrt(sep_ra*sep_ra + sep_dec*sep_dec);
    
    std::cout << "Stella IOTA: RA=" << star_ra << "°, Dec=" << std::showpos << star_dec << std::noshowpos << "°\n";
    std::cout << "Separazione: " << std::setprecision(1) << sep << " arcsec\n";
    std::cout << "  ΔRA = " << std::showpos << sep_ra << std::noshowpos << " arcsec\n";
    std::cout << "  ΔDec = " << std::showpos << sep_dec << std::noshowpos << " arcsec\n";
    
    return 0;
}
