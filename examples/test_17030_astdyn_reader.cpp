/**
 * @file test_17030_astdyn_reader.cpp
 * @brief Test usando il reader AstDyn dalla libreria ITALOccultLibrary
 * 
 * Questo test carica gli elementi con AstDyn (che sa gestire correttamente
 * i formato J2000 Equinoctial Mean di AstDyS) e confronta con JPL Horizons
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/time_utils.h>

using namespace ioccultcalc;

// Struttura semplice per elementi AstDyS
struct AstDySElemsSimple {
    double a, e, i, Omega, omega, M;
    double epoch_mjd;
    std::string name;
};

int main() {
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "   TEST: Asteroide 17030 - AstDyn Conversion vs JPL Horizons\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";
    
    try {
        // ═══════════════════════════════════════════════════════════════
        // STEP 1: Carica elementi AstDyS (equinoziali J2000)
        // ═══════════════════════════════════════════════════════════════
        std::cout << "STEP 1: Caricamento elementi AstDyS (J2000 Equinoctial Mean)\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        // Leggi file
        std::ifstream file("17030_astdys.eq1");
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open 17030_astdys.eq1");
        }
        
        AstDySElemsSimple elem;
        double h_eq, k_eq, p_eq, q_eq, lambda_eq;
        std::string line;
        bool found_name = false, found_equ = false, found_mjd = false;
        
        while (std::getline(file, line)) {
            if (line.find(" EQU ") != std::string::npos && !found_equ) {
                std::istringstream iss(line);
                std::string tag;
                iss >> tag >> elem.a >> h_eq >> k_eq >> p_eq >> q_eq >> lambda_eq;
                found_equ = true;
            }
            if (line.find(" MJD ") != std::string::npos && found_equ && !found_mjd) {
                std::istringstream iss(line);
                std::string tag;
                iss >> tag >> elem.epoch_mjd;
                found_mjd = true;
                break;
            }
            if (!found_name && line[0] != '!' && line[0] != ' ' && 
                line.find("format") == std::string::npos &&
                line.find("rectype") == std::string::npos) {
                elem.name = line;
                found_name = true;
            }
        }
        file.close();
        
        if (!found_equ || !found_mjd) {
            throw std::runtime_error("Could not parse AstDyS file properly");
        }
        
        std::cout << "\n✓ Elementi equinoziali caricati:\n";
        std::cout << "  a:      " << std::fixed << std::setprecision(8) << elem.a << " AU\n";
        std::cout << "  h:      " << h_eq << " (e*sin(ω+Ω))\n";
        std::cout << "  k:      " << k_eq << " (e*cos(ω+Ω))\n";
        std::cout << "  p:      " << p_eq << " (tan(i/2)*sin(Ω))\n";
        std::cout << "  q:      " << q_eq << " (tan(i/2)*cos(Ω))\n";
        std::cout << "  λ:      " << lambda_eq << " rad\n";
        std::cout << "  Epoca:  " << std::setprecision(1) << elem.epoch_mjd << " MJD\n\n";
        
        // ═══════════════════════════════════════════════════════════════
        // STEP 2: Converti equinoziali → Kepleriani (MANUALMENTE come AstDyn)
        // ═══════════════════════════════════════════════════════════════
        std::cout << "STEP 2: Conversione equinoziali → Kepleriani (metodo AstDyn)\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        // Formula di conversione corretta per J2000 Equinoctial Mean:
        elem.e = sqrt(h_eq * h_eq + k_eq * k_eq);
        elem.i = 2.0 * atan(sqrt(p_eq * p_eq + q_eq * q_eq));
        elem.Omega = atan2(p_eq, q_eq);
        if (elem.Omega < 0) elem.Omega += 2.0 * M_PI;
        
        double omega_plus_Omega = atan2(h_eq, k_eq);
        if (omega_plus_Omega < 0) omega_plus_Omega += 2.0 * M_PI;
        
        elem.omega = omega_plus_Omega - elem.Omega;
        if (elem.omega < 0) elem.omega += 2.0 * M_PI;
        
        elem.M = lambda_eq - omega_plus_Omega;
        while (elem.M < 0) elem.M += 2.0 * M_PI;
        while (elem.M >= 2.0 * M_PI) elem.M -= 2.0 * M_PI;
        
        std::cout << "\n✓ Elementi Kepleriani convertiti:\n";
        std::cout << "  a:  " << std::setprecision(8) << elem.a << " AU\n";
        std::cout << "  e:  " << elem.e << "\n";
        std::cout << "  i:  " << std::setprecision(4) << (elem.i * 180.0 / M_PI) << "°\n";
        std::cout << "  Ω:  " << (elem.Omega * 180.0 / M_PI) << "°\n";
        std::cout << "  ω:  " << (elem.omega * 180.0 / M_PI) << "°\n";
        std::cout << "  M:  " << (elem.M * 180.0 / M_PI) << "°\n\n";
        
        // ═══════════════════════════════════════════════════════════════
        // STEP 3: Converti elementi Keplerian → cartesiane (2-body Kepler)
        // ═══════════════════════════════════════════════════════════════
        std::cout << "STEP 3: Propagazione elementi → coordinate cartesiane\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        // Risolvi equazione di Keplero: E - e*sin(E) = M
        double E = elem.M + elem.e * sin(elem.M);
        for (int iter = 0; iter < 10; iter++) {
            E = elem.M + elem.e * sin(E);
        }
        
        // Anomalia vera
        double nu = 2.0 * atan2(sqrt(1.0 + elem.e) * sin(E / 2.0),
                               sqrt(1.0 - elem.e) * cos(E / 2.0));
        
        // Distanza
        double r = elem.a * (1.0 - elem.e * cos(E));
        
        // Coordinate nel piano orbitale
        double x_orb = r * cos(nu);
        double y_orb = r * sin(nu);
        
        // Rotazioni per portare al frame eclittico J2000
        // R_z(Ω) * R_x(i) * R_z(ω)
        double x_ecl = (cos(elem.Omega) * cos(elem.omega) - 
                       sin(elem.Omega) * sin(elem.omega) * cos(elem.i)) * x_orb +
                       (-cos(elem.Omega) * sin(elem.omega) - 
                        sin(elem.Omega) * cos(elem.omega) * cos(elem.i)) * y_orb;
        
        double y_ecl = (sin(elem.Omega) * cos(elem.omega) + 
                       cos(elem.Omega) * sin(elem.omega) * cos(elem.i)) * x_orb +
                       (-sin(elem.Omega) * sin(elem.omega) + 
                        cos(elem.Omega) * cos(elem.omega) * cos(elem.i)) * y_orb;
        
        double z_ecl = sin(elem.omega) * sin(elem.i) * x_orb + 
                      cos(elem.omega) * sin(elem.i) * y_orb;
        
        double r_calc_ecl = sqrt(x_ecl * x_ecl + y_ecl * y_ecl + z_ecl * z_ecl);
        
        std::cout << "\n✓ Coordinate eclittiche (equinoziali J2000):\n";
        std::cout << "  X:  " << std::setprecision(10) << x_ecl << " AU\n";
        std::cout << "  Y:  " << y_ecl << " AU\n";
        std::cout << "  Z:  " << z_ecl << " AU\n";
        std::cout << "  |r|: " << r_calc_ecl << " AU\n\n";
        
        // ═══════════════════════════════════════════════════════════════
        // STEP 3.5: Rotazione ECLITTICO → EQUATORIALE (ICRF J2000)
        // ═══════════════════════════════════════════════════════════════
        std::cout << "STEP 3.5: Rotazione eclittico → equatoriale ICRF\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        // Obliquità eclittica J2000
        constexpr double eps = 23.4392911 * M_PI / 180.0;  // 0.40909302 rad
        double cos_eps = cos(eps);   // 0.91748
        double sin_eps = sin(eps);   // 0.39776
        
        // Rotazione attorno asse X:
        // [x_eq]   [1      0        0    ] [x_ecl]
        // [y_eq] = [0   cos_ε  -sin_ε  ] [y_ecl]
        // [z_eq]   [0   sin_ε   cos_ε  ] [z_ecl]
        
        double x_eq = x_ecl;
        double y_eq = cos_eps * y_ecl - sin_eps * z_ecl;
        double z_eq = sin_eps * y_ecl + cos_eps * z_ecl;
        
        double r_calc = sqrt(x_eq * x_eq + y_eq * y_eq + z_eq * z_eq);
        
        std::cout << "\n✓ Coordinate equatoriali ICRF (J2000):\n";
        std::cout << "  X:  " << std::setprecision(10) << x_eq << " AU\n";
        std::cout << "  Y:  " << y_eq << " AU\n";
        std::cout << "  Z:  " << z_eq << " AU\n";
        std::cout << "  |r|: " << r_calc << " AU\n\n";
        
        // ═══════════════════════════════════════════════════════════════
        // STEP 4: Scarica da JPL Horizons
        // ═══════════════════════════════════════════════════════════════
        std::cout << "STEP 4: Scaricamento da JPL Horizons\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        JPLHorizonsClient horizons;
        horizons.setTimeout(60);
        
        JulianDate epoch_jd;
        epoch_jd.jd = elem.epoch_mjd + 2400000.5;
        
        auto [pos_jpl, vel_jpl] = horizons.getStateVectors("17030", epoch_jd, "@sun");
        
        std::cout << "\n✓ Coordinate JPL Horizons (eliocentriche):\n";
        std::cout << "  X:  " << std::setprecision(10) << pos_jpl.x << " AU\n";
        std::cout << "  Y:  " << pos_jpl.y << " AU\n";
        std::cout << "  Z:  " << pos_jpl.z << " AU\n";
        std::cout << "  |r|: " << pos_jpl.magnitude() << " AU\n\n";
        
        // ═══════════════════════════════════════════════════════════════
        // STEP 5: Confronto
        // ═══════════════════════════════════════════════════════════════
        std::cout << "STEP 5: Confronto AstDyn (manuale) vs JPL\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        double delta_x = x_eq - pos_jpl.x;
        double delta_y = y_eq - pos_jpl.y;
        double delta_z = z_eq - pos_jpl.z;
        double delta_r = sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
        
        std::cout << "\n✓ Differenze cartesiane (ICRF):\n";
        std::cout << "  ΔX:    " << std::scientific << std::setprecision(6) 
                  << delta_x << " AU = " << std::fixed 
                  << (delta_x * 149597870.7) << " km\n";
        std::cout << "  ΔY:    " << std::scientific << delta_y 
                  << " AU = " << std::fixed << (delta_y * 149597870.7) << " km\n";
        std::cout << "  ΔZ:    " << std::scientific << delta_z 
                  << " AU = " << std::fixed << (delta_z * 149597870.7) << " km\n";
        std::cout << "  |Δr|:  " << std::scientific << delta_r 
                  << " AU = " << std::fixed << (delta_r * 149597870.7) << " km\n\n";
        
        // Coordinate sferiche (usando coordinate equatoriali)
        double ra_calc = atan2(y_eq, x_eq) * 180.0 / M_PI;
        double dec_calc = asin(z_eq / r_calc) * 180.0 / M_PI;
        
        double ra_jpl = atan2(pos_jpl.y, pos_jpl.x) * 180.0 / M_PI;
        double dec_jpl = asin(pos_jpl.z / pos_jpl.magnitude()) * 180.0 / M_PI;
        
        double delta_ra_arcsec = (ra_calc - ra_jpl) * 3600.0;
        double delta_dec_arcsec = (dec_calc - dec_jpl) * 3600.0;
        
        std::cout << "✓ Coordinate sferiche:\n";
        std::cout << "  AstDyn: RA=" << std::setprecision(6) << ra_calc 
                  << "°,  Dec=" << dec_calc << "°\n";
        std::cout << "  JPL:    RA=" << ra_jpl << "°,  Dec=" << dec_jpl << "°\n";
        std::cout << "  ΔRA:    " << delta_ra_arcsec << " arcsec\n";
        std::cout << "  ΔDec:   " << delta_dec_arcsec << " arcsec\n\n";
        
        // ═══════════════════════════════════════════════════════════════
        // STEP 6: Valutazione
        // ═══════════════════════════════════════════════════════════════
        std::cout << "STEP 6: Valutazione accuratezza\n";
        std::cout << "───────────────────────────────────────────────────────────────\n\n";
        
        double delta_ang = sqrt(delta_ra_arcsec * delta_ra_arcsec + 
                               delta_dec_arcsec * delta_dec_arcsec);
        
        if (delta_ang < 0.1) {
            std::cout << "✓✓✓ ECCEZIONALE: < 0.1 arcsec\n";
        } else if (delta_ang < 1.0) {
            std::cout << "✓✓  ECCELLENTE: < 1 arcsec\n";
        } else if (delta_ang < 10.0) {
            std::cout << "✓   OTTIMO: < 10 arcsec\n";
        } else if (delta_ang < 60.0) {
            std::cout << "⚠   BUONO: < 1 arcmin\n";
        } else {
            std::cout << "✗   PROBLEMATICO: > 1 arcmin\n";
        }
        
        std::cout << "\n";
        std::cout << "═══════════════════════════════════════════════════════════════\n";
        std::cout << "✓ TEST COMPLETATO\n";
        std::cout << "═══════════════════════════════════════════════════════════════\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ ERRORE: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
