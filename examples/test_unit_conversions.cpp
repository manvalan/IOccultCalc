/**
 * @file test_unit_conversions.cpp
 * @brief Test SPECIFICO conversioni gradi <-> radianti
 * 
 * Questo test verifica TUTTE le conversioni di unità nel pipeline
 * per trovare errori di conversione gradi/radianti
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <ioccultcalc/orbital_elements.h>

using namespace ioccultcalc;

int main() {
    std::cout << "===================================================================\n";
    std::cout << "  TEST: Conversioni Gradi <-> Radianti nel Pipeline\n";
    std::cout << "===================================================================\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 1: Leggi file ed estrai elementi equinoziali
    // ═══════════════════════════════════════════════════════════════
    std::cout << "STEP 1: Parsing file 17030_astdys.eq1\n";
    std::cout << "-------------------------------------------------------------------\n";
    
    std::ifstream file("17030_astdys.eq1");
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open 17030_astdys.eq1\n";
        return 1;
    }
    
    double a_eq, h_eq, k_eq, p_eq, q_eq, lambda_eq, mjd_epoch;
    std::string line;
    bool found_equ = false, found_mjd = false;
    
    while (std::getline(file, line)) {
        if (line.find(" EQU ") != std::string::npos && !found_equ) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag >> a_eq >> h_eq >> k_eq >> p_eq >> q_eq >> lambda_eq;
            found_equ = true;
        }
        if (line.find(" MJD ") != std::string::npos && found_equ && !found_mjd) {
            std::istringstream iss(line);
            std::string tag;
            iss >> tag >> mjd_epoch;
            found_mjd = true;
            break;
        }
    }
    file.close();
    
    if (!found_equ || !found_mjd) {
        std::cerr << "ERROR: Could not parse 17030_astdys.eq1\n";
        return 1;
    }
    
    std::cout << "File content (raw values):\n";
    std::cout << "  a        = " << std::scientific << std::setprecision(15) << a_eq << " AU\n";
    std::cout << "  h        = " << h_eq << "\n";
    std::cout << "  k        = " << k_eq << "\n";
    std::cout << "  p        = " << p_eq << "\n";
    std::cout << "  q        = " << q_eq << "\n";
    std::cout << "  λ (RAW)  = " << lambda_eq << " <-- UNITÀ???\n";
    std::cout << "  MJD      = " << mjd_epoch << " (TDT)\n";
    
    std::cout << "\n⚠ DOMANDA CRITICA: λ = " << lambda_eq 
              << " è in RADIANTI o GRADI?\n";
    
    // Converte λ da radianti a gradi per vedere quale ha senso
    double lambda_deg_from_rad = lambda_eq * 180.0 / M_PI;
    double lambda_rad_from_deg = lambda_eq * M_PI / 180.0;
    
    std::cout << "\nInterpretazioni possibili:\n";
    std::cout << "  Se λ è in RADIANTI:\n";
    std::cout << "    λ = " << lambda_eq << " rad = " << lambda_deg_from_rad 
              << "° (normalizzato: " << fmod(lambda_deg_from_rad, 360.0) << "°)\n";
    std::cout << "    Questo ha senso? " << (fmod(lambda_deg_from_rad, 360.0) > 360 ? "NO - troppo grande!" : "Forse") << "\n";
    
    std::cout << "\n  Se λ è in GRADI:\n";
    std::cout << "    λ = " << lambda_eq << "° = " << lambda_rad_from_deg 
              << " rad\n";
    std::cout << "    Questo ha senso? " << (lambda_eq < 360 ? "SÌ!" : "NO") << "\n";
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 2: Crea EquinoctialElements e converti
    // ═══════════════════════════════════════════════════════════════
    std::cout << "\n" << std::string(63, '=') << "\n";
    std::cout << "STEP 2: Creazione EquinoctialElements\n";
    std::cout << "-------------------------------------------------------------------\n";
    
    // TENTATIVO 1: λ in radianti (come nel file)
    std::cout << "\nTENTATIVO 1: Assumi λ in RADIANTI\n";
    
    EquinoctialElements eq1;
    eq1.a = a_eq;
    eq1.h = h_eq;
    eq1.k = k_eq;
    eq1.p = p_eq;
    eq1.q = q_eq;
    eq1.lambda = lambda_eq;  // Come nel file (radianti?)
    eq1.epoch.jd = mjd_epoch + 2400000.5;
    eq1.designation = "17030";
    
    std::cout << "  λ = " << std::scientific << eq1.lambda << " rad\n";
    std::cout << "  λ = " << std::fixed << std::defaultfloat 
              << (eq1.lambda * 180.0 / M_PI) << "° (normalizzato: " 
              << fmod(eq1.lambda * 180.0 / M_PI, 360.0) << "°)\n";
    
    OrbitalElements kep1 = eq1.toKeplerian();
    
    std::cout << "\nElementi Keplerian risultanti:\n";
    std::cout << "  a = " << std::setprecision(8) << kep1.a << " AU\n";
    std::cout << "  e = " << kep1.e << "\n";
    std::cout << "  i = " << kep1.i << " rad = " << (kep1.i * 180.0 / M_PI) << "°\n";
    std::cout << "  Ω = " << kep1.Omega << " rad = " << (kep1.Omega * 180.0 / M_PI) << "°\n";
    std::cout << "  ω = " << kep1.omega << " rad = " << (kep1.omega * 180.0 / M_PI) << "°\n";
    std::cout << "  M = " << kep1.M << " rad = " << (kep1.M * 180.0 / M_PI) << "°\n";
    
    // TENTATIVO 2: λ in gradi (convertito a radianti)
    std::cout << "\n" << std::string(63, '-') << "\n";
    std::cout << "TENTATIVO 2: Assumi λ in GRADI (converti a radianti)\n";
    
    EquinoctialElements eq2;
    eq2.a = a_eq;
    eq2.h = h_eq;
    eq2.k = k_eq;
    eq2.p = p_eq;
    eq2.q = q_eq;
    eq2.lambda = lambda_eq * M_PI / 180.0;  // Converti da gradi a radianti
    eq2.epoch.jd = mjd_epoch + 2400000.5;
    eq2.designation = "17030";
    
    std::cout << "  λ (input) = " << lambda_eq << "° (convertito in radianti)\n";
    std::cout << "  λ = " << std::scientific << eq2.lambda << " rad\n";
    std::cout << "  λ = " << std::fixed << std::defaultfloat 
              << (eq2.lambda * 180.0 / M_PI) << "° (normalizzato: " 
              << fmod(eq2.lambda * 180.0 / M_PI, 360.0) << "°)\n";
    
    OrbitalElements kep2 = eq2.toKeplerian();
    
    std::cout << "\nElementi Keplerian risultanti:\n";
    std::cout << "  a = " << std::setprecision(8) << kep2.a << " AU\n";
    std::cout << "  e = " << kep2.e << "\n";
    std::cout << "  i = " << kep2.i << " rad = " << (kep2.i * 180.0 / M_PI) << "°\n";
    std::cout << "  Ω = " << kep2.Omega << " rad = " << (kep2.Omega * 180.0 / M_PI) << "°\n";
    std::cout << "  ω = " << kep2.omega << " rad = " << (kep2.omega * 180.0 / M_PI) << "°\n";
    std::cout << "  M = " << kep2.M << " rad = " << (kep2.M * 180.0 / M_PI) << "°\n";
    
    // ═══════════════════════════════════════════════════════════════
    // STEP 3: Confronto con valori attesi dal documento
    // ═══════════════════════════════════════════════════════════════
    std::cout << "\n" << std::string(63, '=') << "\n";
    std::cout << "STEP 3: Confronto con documento ITALOccultLibrary\n";
    std::cout << "-------------------------------------------------------------------\n";
    
    std::cout << "\nValori dal documento (per i DATI DI TEST):\n";
    std::cout << "  e = 0.045407\n";
    std::cout << "  i = 2.9046° = 0.050707 rad\n";
    std::cout << "  Ω = 94.058° = 1.6417 rad\n";
    std::cout << "  ω = 110.28° = 1.9245 rad\n";
    std::cout << "  M = 25.45° = 0.4441 rad\n";
    
    std::cout << "\nValori dal nostro file (TENTATIVO 1 - λ in radianti):\n";
    std::cout << "  e = " << std::setprecision(6) << kep1.e << "\n";
    std::cout << "  i = " << (kep1.i * 180.0 / M_PI) << "° = " << kep1.i << " rad\n";
    std::cout << "  Ω = " << (kep1.Omega * 180.0 / M_PI) << "° = " << kep1.Omega << " rad\n";
    std::cout << "  ω = " << (kep1.omega * 180.0 / M_PI) << "° = " << kep1.omega << " rad\n";
    std::cout << "  M = " << (kep1.M * 180.0 / M_PI) << "° = " << kep1.M << " rad\n";
    
    std::cout << "\nValori dal nostro file (TENTATIVO 2 - λ in gradi):\n";
    std::cout << "  e = " << kep2.e << "\n";
    std::cout << "  i = " << (kep2.i * 180.0 / M_PI) << "° = " << kep2.i << " rad\n";
    std::cout << "  Ω = " << (kep2.Omega * 180.0 / M_PI) << "° = " << kep2.Omega << " rad\n";
    std::cout << "  ω = " << (kep2.omega * 180.0 / M_PI) << "° = " << kep2.omega << " rad\n";
    std::cout << "  M = " << (kep2.M * 180.0 / M_PI) << "° = " << kep2.M << " rad\n";
    
    std::cout << "\n" << std::string(63, '=') << "\n";
    std::cout << "CONCLUSIONE:\n";
    std::cout << "-------------------------------------------------------------------\n";
    
    double error1_M = fabs((kep1.M * 180.0 / M_PI) - 25.45);
    double error2_M = fabs((kep2.M * 180.0 / M_PI) - 101.99);  // Aspetto M diverso per file reale
    
    std::cout << "\nSe λ è in RADIANTI (TENTATIVO 1):\n";
    std::cout << "  Errore in M: " << error1_M << "° (vs atteso ~25.45°)\n";
    
    std::cout << "\nSe λ è in GRADI (TENTATIVO 2):\n";
    std::cout << "  M = " << (kep2.M * 180.0 / M_PI) << "° (file reale - confronta con Horizons)\n";
    
    return 0;
}
