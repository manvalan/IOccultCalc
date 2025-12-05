/**
 * @file test_17030_simple_conversion.cpp
 * @brief Test semplice: verificare se gli elementi AstDyS, convertiti in 
 *        coordinate cartesiane, corrispondono a JPL Horizons
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ioccultcalc/orbital_elements.h>
#include <ioccultcalc/time_utils.h>
#include <ioccultcalc/jpl_horizons_client.h>

using namespace ioccultcalc;

bool loadAstDysElements(const std::string& filename, EquinoctialElements& elements) {
    std::ifstream file(filename);
    if (!file.good()) {
        return false;
    }
    
    std::string line;
    bool foundEqu = false, foundMjd = false;
    
    while (std::getline(file, line)) {
        if (line.find(" EQU ") != std::string::npos && !foundEqu) {
            std::istringstream iss(line);
            std::string tag;
            double a, h, k, p, q, lambda;
            if (iss >> tag >> a >> h >> k >> p >> q >> lambda) {
                elements.a = a;
                elements.h = h;
                elements.k = k;
                elements.p = p;
                elements.q = q;
                elements.lambda = lambda;
                foundEqu = true;
            }
        }
        if (line.find(" MJD ") != std::string::npos && foundEqu && !foundMjd) {
            std::istringstream iss(line);
            std::string tag;
            double mjd;
            if (iss >> tag >> mjd) {
                elements.epoch.jd = mjd + 2400000.5;
                foundMjd = true;
                break;
            }
        }
    }
    file.close();
    
    return foundEqu && foundMjd;
}

int main() {
    std::cout << "===== TEST CONVERSIONE ELEMENTI AstDyS =====\n\n";
    
    try {
        EquinoctialElements equElements;
        if (!loadAstDysElements("17030_astdys.eq1", equElements)) {
            std::cerr << "Errore caricamento elementi\n";
            return 1;
        }
        
        std::cout << "Elementi Equinoziali (AstDyS):\n";
        std::cout << "  a      = " << std::fixed << std::setprecision(8) << equElements.a << " AU\n";
        std::cout << "  h      = " << equElements.h << "\n";
        std::cout << "  k      = " << equElements.k << "\n";
        std::cout << "  p      = " << equElements.p << "\n";
        std::cout << "  q      = " << equElements.q << "\n";
        std::cout << "  lambda = " << equElements.lambda << " rad = " 
                  << (equElements.lambda * 180.0 / M_PI) << "°\n";
        std::cout << "  Epoca  = JD " << equElements.epoch.jd << "\n\n";
        
        // Converti in Keplerian
        auto kepElements = equElements.toKeplerian();
        std::cout << "Elementi Keplerian:\n";
        std::cout << "  a = " << std::setprecision(8) << kepElements.a << " AU\n";
        std::cout << "  e = " << kepElements.e << "\n";
        std::cout << "  i = " << std::setprecision(4) << (kepElements.i * 180.0 / M_PI) << "°\n";
        std::cout << "  Ω = " << (kepElements.Omega * 180.0 / M_PI) << "°\n";
        std::cout << "  ω = " << (kepElements.omega * 180.0 / M_PI) << "°\n";
        std::cout << "  M = " << (kepElements.M * 180.0 / M_PI) << "°\n\n";
        
        // Converti in coordinate cartesiane usando il propagatore
        // Usiamo la formula di Keplero pura
        double E = kepElements.M;  // Prima approssimazione
        for (int i = 0; i < 10; i++) {
            E = kepElements.M + kepElements.e * sin(E);  // Equazione di Keplero
        }
        double nu = 2.0 * atan2(sqrt(1.0 + kepElements.e) * sin(E/2.0),
                               sqrt(1.0 - kepElements.e) * cos(E/2.0));
        
        double r = kepElements.a * (1.0 - kepElements.e * cos(E));
        
        // Coordinate in piano orbitale
        double x_orb = r * cos(nu);
        double y_orb = r * sin(nu);
        double z_orb = 0.0;
        
        // Rotazioni per portare al frame equatoriale J2000
        // Rotazione 3 (per Ω): attorno all'asse z
        double xr1 = x_orb * cos(kepElements.Omega) - y_orb * sin(kepElements.Omega);
        double yr1 = x_orb * sin(kepElements.Omega) + y_orb * cos(kepElements.Omega);
        double zr1 = z_orb;
        
        // Rotazione 1 (per i): attorno all'asse x
        double xr2 = xr1;
        double yr2 = yr1 * cos(kepElements.i) - zr1 * sin(kepElements.i);
        double zr2 = yr1 * sin(kepElements.i) + zr1 * cos(kepElements.i);
        
        // Rotazione 3 (per ω): attorno all'asse z
        double x_eq = xr2 * cos(kepElements.omega) - yr2 * sin(kepElements.omega);
        double y_eq = xr2 * sin(kepElements.omega) + yr2 * cos(kepElements.omega);
        double z_eq = zr2;
        
        std::cout << "Posizione calcolata da elementi Keplerian (conversione pura):\n";
        std::cout << "  X = " << std::setprecision(10) << x_eq << " AU\n";
        std::cout << "  Y = " << y_eq << " AU\n";
        std::cout << "  Z = " << z_eq << " AU\n";
        std::cout << "  |r| = " << sqrt(x_eq*x_eq + y_eq*y_eq + z_eq*z_eq) << " AU\n\n";
        
        // Scarica da JPL
        JPLHorizonsClient horizons;
        horizons.setTimeout(60);
        
        auto [posJPL, velJPL] = horizons.getStateVectors("17030", equElements.epoch, "@sun");
        
        std::cout << "Posizione JPL Horizons (heliocentric):\n";
        std::cout << "  X = " << std::setprecision(10) << posJPL.x << " AU\n";
        std::cout << "  Y = " << posJPL.y << " AU\n";
        std::cout << "  Z = " << posJPL.z << " AU\n";
        std::cout << "  |r| = " << posJPL.magnitude() << " AU\n\n";
        
        // Differenza
        std::cout << "CONFRONTO:\n";
        std::cout << "  ΔX = " << std::scientific << std::setprecision(6) 
                  << (x_eq - posJPL.x) << " AU = "
                  << std::fixed << (x_eq - posJPL.x) * 149597870.7 << " km\n";
        std::cout << "  ΔY = " << std::scientific << (y_eq - posJPL.y) 
                  << " AU = " << std::fixed << (y_eq - posJPL.y) * 149597870.7 << " km\n";
        std::cout << "  ΔZ = " << std::scientific << (z_eq - posJPL.z) 
                  << " AU = " << std::fixed << (z_eq - posJPL.z) * 149597870.7 << " km\n\n";
        
        double delta_r = sqrt((x_eq - posJPL.x)*(x_eq - posJPL.x) +
                            (y_eq - posJPL.y)*(y_eq - posJPL.y) +
                            (z_eq - posJPL.z)*(z_eq - posJPL.z));
        std::cout << "  |Δpos| = " << std::scientific << delta_r 
                  << " AU = " << std::fixed << delta_r * 149597870.7 << " km\n";
        std::cout << "  Percentuale: " << (delta_r / posJPL.magnitude() * 100.0) << "%\n\n";
        
        // Coordinate sferiche
        std::cout << "COORDINATE SFERICHE:\n";
        double ra_calc = atan2(y_eq, x_eq) * 180.0 / M_PI;
        double dec_calc = asin(z_eq / sqrt(x_eq*x_eq + y_eq*y_eq + z_eq*z_eq)) * 180.0 / M_PI;
        
        double ra_jpl = atan2(posJPL.y, posJPL.x) * 180.0 / M_PI;
        double dec_jpl = asin(posJPL.z / posJPL.magnitude()) * 180.0 / M_PI;
        
        std::cout << "  Calc: RA = " << ra_calc << "°,  Dec = " << dec_calc << "°\n";
        std::cout << "  JPL:  RA = " << ra_jpl << "°,  Dec = " << dec_jpl << "°\n";
        std::cout << "  ΔRA = " << (ra_calc - ra_jpl) * 3600.0 << " arcsec\n";
        std::cout << "  ΔDec = " << (dec_calc - dec_jpl) * 3600.0 << " arcsec\n";
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
