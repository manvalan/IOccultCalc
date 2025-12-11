/**
 * @file test_17030_position_verification.cpp
 * @brief Verifica posizione asteroide 17030 il 28 novembre 2025
 */

#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/astdyn_propagation_helper.h"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "=" << std::string(70, '=') << "\n";
    std::cout << "VERIFICA POSIZIONE ASTEROIDE 17030\n";
    std::cout << "=" << std::string(70, '=') << "\n\n";
    
    // Data target: 28 novembre 2025
    // JD 2460311.5 = MJD 60311.0
    double target_mjd = 60311.0;
    
    std::cout << "Data target: 28 novembre 2025\n";
    std::cout << "MJD: " << std::fixed << std::setprecision(6) << target_mjd << "\n\n";
    
    // Carica elementi orbitali
    std::cout << "Caricamento elementi orbitali da AstDyS...\n";
    ioccultcalc::AstDysClient client;
    auto elements = client.getElements("17030");
    
    std::cout << "Elementi caricati:\n";
    std::cout << "  Epoca: MJD " << elements.epoch.toMJD() << "\n";
    std::cout << "  a = " << elements.a << " AU\n";
    std::cout << "  e = " << elements.e << "\n";
    std::cout << "  i = " << elements.i * 180.0 / M_PI << " deg\n";
    std::cout << "  Omega = " << elements.Omega * 180.0 / M_PI << " deg\n";
    std::cout << "  omega = " << elements.omega * 180.0 / M_PI << " deg\n";
    std::cout << "  lambda = " << elements.lambda * 180.0 / M_PI << " deg\n\n";
    
    // Propaga e calcola RA/Dec
    std::cout << "Propagazione a MJD " << target_mjd << "...\n";
    auto& helper = ioccultcalc::AstDynPropagationHelper::getInstance();
    
    // Converti elementi equinoziali in AstDySElements
    auto astdys = ioccultcalc::AstDynPropagationHelper::convertFromEquinoctial(elements);
    
    // Propaga
    auto astdys_prop = helper.propagate(astdys, target_mjd);
    
    // Calcola RA/Dec
    auto [ra, dec] = helper.getRADec(astdys_prop, target_mjd);
    
    std::cout << "\nCoordinate calcolate:\n";
    std::cout << "  RA: " << std::fixed << std::setprecision(6) << ra << "°\n";
    std::cout << "  Dec: " << dec << "°\n\n";
    
    // Coordinate stella attesa
    double ra_star = 73.416109;  // 4h 53m 39.866s
    double dec_star = 20.331661;  // +20° 19' 53.981"
    
    std::cout << "Coordinate stella UCAC4 552-011427:\n";
    std::cout << "  RA: " << ra_star << "°\n";
    std::cout << "  Dec: " << dec_star << "°\n\n";
    
    // Calcola separazione angolare
    double ra1_rad = ra * M_PI / 180.0;
    double dec1_rad = dec * M_PI / 180.0;
    double ra2_rad = ra_star * M_PI / 180.0;
    double dec2_rad = dec_star * M_PI / 180.0;
    
    double cos_sep = std::sin(dec1_rad) * std::sin(dec2_rad) + 
                     std::cos(dec1_rad) * std::cos(dec2_rad) * 
                     std::cos(ra1_rad - ra2_rad);
    cos_sep = std::max(-1.0, std::min(1.0, cos_sep));
    double sep_rad = std::acos(cos_sep);
    double sep_deg = sep_rad * 180.0 / M_PI;
    double sep_arcsec = sep_deg * 3600.0;
    
    std::cout << "Separazione angolare:\n";
    std::cout << "  " << sep_deg << "° (" << sep_arcsec << " arcsec)\n\n";
    
    if (sep_arcsec < 30) {
        std::cout << "✓ Separazione CORRETTA! La stella è vicina al percorso.\n";
    } else {
        std::cout << "❌ Separazione TROPPO GRANDE! C'è un errore nelle coordinate.\n";
        std::cout << "   Verificare:\n";
        std::cout << "   1. Propagazione dall'epoca originale\n";
        std::cout << "   2. Conversione baricentrico -> geocentrico\n";
        std::cout << "   3. Sistema di riferimento (ICRF J2000)\n";
    }
    
    return 0;
}

