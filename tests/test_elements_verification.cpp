/**
 * @file test_elements_verification.cpp
 * @brief Verifica elementi equinoctiali eclittici medi scaricati da AstDyS
 */

#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/astdyn_propagation_helper.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

int main() {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  VERIFICA ELEMENTI EQUINOCTIALI ECLITTICI MEDI          ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    try {
        // Scarica elementi direttamente
        AstDysClient client;
        EquinoctialElements eq = client.getElements("17030");
        
        std::cout << "Elementi Equinoctiali scaricati da AstDyS:\n";
        std::cout << "  Designation: " << eq.designation << "\n";
        std::cout << "  Name: " << eq.name << "\n";
        std::cout << "  Epoch (JD): " << std::fixed << std::setprecision(5) << eq.epoch.jd << "\n";
        std::cout << "  Epoch (MJD): " << (eq.epoch.jd - 2400000.5) << "\n";
        std::cout << "  a (AU): " << std::setprecision(10) << eq.a << "\n";
        std::cout << "  h: " << eq.h << "\n";
        std::cout << "  k: " << eq.k << "\n";
        std::cout << "  p: " << eq.p << "\n";
        std::cout << "  q: " << eq.q << "\n";
        std::cout << "  lambda (deg): " << eq.lambda * 180.0 / M_PI << "\n";
        std::cout << "  H: " << eq.H << "\n";
        std::cout << "  G: " << eq.G << "\n\n";
        
        // Converti in AstDySElements
        AstDySElements astdys = AstDynPropagationHelper::convertFromEquinoctial(eq);
        
        std::cout << "Elementi convertiti in AstDySElements (Kepleriani):\n";
        std::cout << "  Epoch (MJD TDB): " << std::setprecision(5) << astdys.epoch_mjd << "\n";
        std::cout << "  a (AU): " << std::setprecision(10) << astdys.a << "\n";
        std::cout << "  e: " << astdys.e << "\n";
        std::cout << "  i (deg): " << astdys.i << "\n";
        std::cout << "  Omega (deg): " << astdys.Omega << "\n";
        std::cout << "  omega (deg): " << astdys.omega << "\n";
        std::cout << "  M (deg): " << astdys.M << "\n\n";
        
        // Verifica formule di conversione
        std::cout << "Verifica conversioni:\n";
        double e_calc = sqrt(eq.h * eq.h + eq.k * eq.k);
        std::cout << "  e calcolato: " << e_calc << " (dovrebbe essere " << astdys.e << ")\n";
        
        double i_calc = 2.0 * atan(sqrt(eq.p * eq.p + eq.q * eq.q));
        std::cout << "  i calcolato: " << i_calc * 180.0 / M_PI << " deg (dovrebbe essere " << astdys.i << ")\n";
        
        double Omega_calc = atan2(eq.p, eq.q);
        if (Omega_calc < 0) Omega_calc += 2.0 * M_PI;
        std::cout << "  Omega calcolato: " << Omega_calc * 180.0 / M_PI << " deg (dovrebbe essere " << astdys.Omega << ")\n";
        
        double omega_plus_Omega = atan2(eq.h, eq.k);
        double omega_calc = omega_plus_Omega - Omega_calc;
        if (omega_calc < 0) omega_calc += 2.0 * M_PI;
        std::cout << "  omega calcolato: " << omega_calc * 180.0 / M_PI << " deg (dovrebbe essere " << astdys.omega << ")\n";
        
        double lambda_deg = eq.lambda * 180.0 / M_PI;
        double M_calc = lambda_deg - omega_plus_Omega * 180.0 / M_PI;
        while (M_calc < 0) M_calc += 360.0;
        while (M_calc >= 360.0) M_calc -= 360.0;
        std::cout << "  M calcolato: " << M_calc << " deg (dovrebbe essere " << astdys.M << ")\n";
        
        // Verifica che gli elementi siano in eclittico medio
        std::cout << "\nNOTA: Gli elementi equinoctiali di AstDyS sono in:\n";
        std::cout << "  - Sistema di riferimento: Eclittico medio J2000\n";
        std::cout << "  - Epoca: " << astdys.epoch_mjd << " MJD TDB\n";
        std::cout << "  - AstDyn si aspetta elementi kepleriani in eclittico J2000\n";
        std::cout << "  - La conversione dovrebbe preservare il sistema di riferimento\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

