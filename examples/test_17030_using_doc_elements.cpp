/**
 * @file test_17030_using_doc_elements.cpp
 * @brief Test usando gli elementi dal documento ITALOccultLibrary per verificare se il nostro algoritmo è corretto
 */

#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "   TEST: Usando elementi dal documento ITALOccultLibrary\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // FASE 1: Elementi equinoziali DAL DOCUMENTO
    // ═══════════════════════════════════════════════════════════════
    std::cout << "FASE 1: Elementi equinoziali (DAL DOCUMENTO)\n";
    std::cout << "───────────────────────────────────────────────────────────────\n";
    
    // Dal documento ITALOccultLibrary (Tabella FASE 0 INPUT):
    double a = 3.17547321;
    double h = -0.01896287;
    double k = -0.04127282;
    double p = 0.02458228;
    double q = -0.00620313;
    double lambda_deg = 229.790880;  // gradi dal documento!
    double lambda = lambda_deg * M_PI / 180.0;  // converti a radianti
    
    std::cout << "Dal documento:\n";
    std::cout << "  a = " << std::fixed << std::setprecision(10) << a << " AU\n";
    std::cout << "  h = " << h << "\n";
    std::cout << "  k = " << k << "\n";
    std::cout << "  p = " << p << "\n";
    std::cout << "  q = " << q << "\n";
    std::cout << "  λ = " << lambda_deg << "° = " << lambda << " rad\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // FASE 2: Conversione equinoziali → Kepleriani
    // ═══════════════════════════════════════════════════════════════
    std::cout << "FASE 2: Conversione equinoziali → Kepleriani\n";
    std::cout << "───────────────────────────────────────────────────────────────\n";
    
    double e = sqrt(h*h + k*k);
    double i = 2.0 * atan(sqrt(p*p + q*q));
    double Omega = atan2(p, q);
    if (Omega < 0) Omega += 2.0 * M_PI;
    double varpi = atan2(h, k);
    if (varpi < 0) varpi += 2.0 * M_PI;
    double omega = varpi - Omega;
    if (omega < 0) omega += 2.0 * M_PI;
    double M = lambda - varpi;
    while (M < 0) M += 2.0 * M_PI;
    while (M >= 2.0 * M_PI) M -= 2.0 * M_PI;
    
    std::cout << "✓ Calcolati Kepleriani:\n";
    std::cout << "  a = " << a << " AU\n";
    std::cout << "  e = " << std::setprecision(8) << e << "\n";
    std::cout << "  i = " << (i * 180.0 / M_PI) << "°\n";
    std::cout << "  Ω = " << (Omega * 180.0 / M_PI) << "°\n";
    std::cout << "  ω = " << (omega * 180.0 / M_PI) << "°\n";
    std::cout << "  M = " << (M * 180.0 / M_PI) << "°\n\n";
    
    std::cout << "✓ Attesi dal documento:\n";
    std::cout << "  a = 3.17547321 AU\n";
    std::cout << "  e = 0.045407\n";
    std::cout << "  i = 2.9046°\n";
    std::cout << "  Ω = 94.058°  ← ATTENZIONE: Il nostro Ω calcolato da p,q è DIVERSO!\n";
    std::cout << "  ω = 110.28°\n";
    std::cout << "  M = 25.45°\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // FASE 3: Risoluzione equazione di Keplero
    // ═══════════════════════════════════════════════════════════════
    std::cout << "FASE 3: Equazione di Keplero\n";
    std::cout << "───────────────────────────────────────────────────────────────\n";
    
    // Nota: Usiamo M calcolato da λ - ϖ
    double E = M + e * sin(M);
    for (int iter = 0; iter < 10; iter++) {
        E = M + e * sin(E);
    }
    
    double nu = 2.0 * atan2(sqrt(1.0 + e) * sin(E/2.0), sqrt(1.0 - e) * cos(E/2.0));
    double r = a * (1.0 - e * cos(E));
    
    std::cout << "✓ Anomalia eccentrica: E = " << E << " rad = " << (E * 180.0 / M_PI) << "°\n";
    std::cout << "  Attesa: E = 0.4461 rad = 25.56°\n\n";
    std::cout << "✓ Anomalia vera: ν = " << nu << " rad = " << (nu * 180.0 / M_PI) << "°\n";
    std::cout << "  Attesa: ν = 0.4593 rad = 26.32°\n\n";
    std::cout << "✓ Distanza: r = " << std::setprecision(6) << r << " AU\n";
    std::cout << "  Attesa: r = 3.0454 AU\n\n";
    
    // ═════════════════════════════════════════════════════════════════
    // FASE 4: Posizione nel piano orbitale
    // ═════════════════════════════════════════════════════════════════
    std::cout << "FASE 4: Coordinate nel piano orbitale\n";
    std::cout << "───────────────────────────────────────────────────────────────\n";
    
    double x_orb = r * cos(nu);
    double y_orb = r * sin(nu);
    double z_orb = 0.0;
    
    std::cout << "✓ (x_orb, y_orb, z_orb) = (" << x_orb << ", " << y_orb << ", " << z_orb << ") AU\n";
    std::cout << "  Attesa: (2.7393, 1.3337, 0) AU\n\n";
    
    // ═════════════════════════════════════════════════════════════════
    // FASE 5: Rotazioni Euleriane
    // ═════════════════════════════════════════════════════════════════
    std::cout << "FASE 5: Rotazioni Euleriane\n";
    std::cout << "───────────────────────────────────────────────────────────────\n";
    
    double cos_Omega = cos(Omega);
    double sin_Omega = sin(Omega);
    double cos_omega = cos(omega);
    double sin_omega = sin(omega);
    double cos_i = cos(i);
    double sin_i = sin(i);
    
    double x_ecl = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * x_orb +
                   (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * y_orb;
    
    double y_ecl = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * x_orb +
                   (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * y_orb;
    
    double z_ecl = sin_omega * sin_i * x_orb + cos_omega * sin_i * y_orb;
    
    std::cout << "✓ (x_ecl, y_ecl, z_ecl) = (" << x_ecl << ", " << y_ecl << ", " << z_ecl << ") AU\n";
    std::cout << "  Attesa: (-1.9899, -2.3073, 0.1093) AU\n\n";
    
    // ═════════════════════════════════════════════════════════════════
    // FASE 6: Rotazione eclittica → ICRF
    // ═════════════════════════════════════════════════════════════════
    std::cout << "FASE 6: Rotazione eclittica → ICRF\n";
    std::cout << "───────────────────────────────────────────────────────────────\n";
    
    double epsilon_rad = 0.40909302;  // obliquità in radianti
    double cos_epsilon = cos(epsilon_rad);
    double sin_epsilon = sin(epsilon_rad);
    
    double x_eq = x_ecl;
    double y_eq = y_ecl * cos_epsilon - z_ecl * sin_epsilon;
    double z_eq = y_ecl * sin_epsilon + z_ecl * cos_epsilon;
    
    std::cout << "✓ (x_eq, y_eq, z_eq) = (" << x_eq << ", " << y_eq << ", " << z_eq << ") AU\n";
    std::cout << "  Attesa: (-1.9899, -2.1599, -0.8174) AU\n\n";
    
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    double dx = x_eq - (-1.9899);
    double dy = y_eq - (-2.1599);
    double dz = z_eq - (-0.8174);
    
    if (fabs(dx) < 0.01 && fabs(dy) < 0.01 && fabs(dz) < 0.01) {
        std::cout << "✅ RISULTATI CORRETTI! (Errori < 0.01 AU)\n";
    } else {
        std::cout << "⚠️  Errori: Δx=" << dx << ", Δy=" << dy << ", Δz=" << dz << "\n";
    }
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    
    return 0;
}
