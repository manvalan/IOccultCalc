/**
 * @file test_17030_phase_by_phase_debug.cpp
 * @brief Debug passo dopo passo - Confronta i nostri calcoli con il documento
 */

#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "   DEBUG: Confronto valori calcolati vs documento ITALOccultLibrary\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // FASE 1: Elementi equinoziali caricati
    // ═══════════════════════════════════════════════════════════════
    std::cout << "FASE 1: Elementi equinoziali\n";
    std::cout << "───────────────────────────────────────────────────────────────\n";
    
    double a = 3.17547321;
    double h = -0.01896287;
    double k = -0.04127282;
    double p = 0.02458228;
    double q = -0.00620313;
    double lambda = 74.46741573;  // rad
    
    std::cout << "Caricati:\n";
    std::cout << "  a = " << std::fixed << std::setprecision(10) << a << " AU\n";
    std::cout << "  h = " << h << "\n";
    std::cout << "  k = " << k << "\n";
    std::cout << "  p = " << p << "\n";
    std::cout << "  q = " << q << "\n";
    std::cout << "  λ = " << lambda << " rad\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // FASE 2: Conversione equinoziali → Kepleriani
    // ═══════════════════════════════════════════════════════════════
    std::cout << "FASE 2: Conversione equinoziali → Kepleriani\n";
    std::cout << "───────────────────────────────────────────────────────────────\n";
    
    // Eccentricità
    double e_calc = sqrt(h*h + k*k);
    std::cout << "✓ Eccentricità:\n";
    std::cout << "  Calcolata: e = √(" << h << "² + " << k << "²) = " 
              << std::setprecision(8) << e_calc << "\n";
    std::cout << "  Attesa:    e = 0.045407\n";
    std::cout << "  Match: " << (fabs(e_calc - 0.045407) < 0.000001 ? "✓" : "✗") << "\n\n";
    
    // Inclinazione
    double tan_half_i = sqrt(p*p + q*q);
    double i_calc = 2.0 * atan(tan_half_i);
    std::cout << "✓ Inclinazione:\n";
    std::cout << "  Calcolata: i = 2·arctan(√(" << p << "² + " << q << "²))\n";
    std::cout << "           = " << i_calc << " rad = " << (i_calc * 180.0 / M_PI) 
              << "°\n";
    std::cout << "  Attesa:    i = 2.9046°\n";
    std::cout << "  Match: " << (fabs(i_calc * 180.0 / M_PI - 2.9046) < 0.01 ? "✓" : "✗") << "\n\n";
    
    // Nodo ascendente
    double Omega_calc = atan2(p, q);
    if (Omega_calc < 0) Omega_calc += 2.0 * M_PI;
    std::cout << "✓ Nodo ascendente:\n";
    std::cout << "  Calcolata: Ω = atan2(" << p << ", " << q << ") = " 
              << Omega_calc << " rad = " << (Omega_calc * 180.0 / M_PI) << "°\n";
    std::cout << "  Attesa:    Ω = 94.058°\n";
    std::cout << "  Match: " << (fabs(Omega_calc * 180.0 / M_PI - 94.058) < 1.0 ? "✓" : "✗") << "\n\n";
    
    // Longitudine del perielio (ϖ = ω + Ω)
    double varpi_calc = atan2(h, k);
    if (varpi_calc < 0) varpi_calc += 2.0 * M_PI;
    std::cout << "✓ Longitudine perielio:\n";
    std::cout << "  Calcolata: ϖ = atan2(" << h << ", " << k << ") = " 
              << varpi_calc << " rad = " << (varpi_calc * 180.0 / M_PI) << "°\n";
    std::cout << "  Attesa:    ϖ = 204.34°\n";
    std::cout << "  Match: " << (fabs(varpi_calc * 180.0 / M_PI - 204.34) < 1.0 ? "✓" : "✗") << "\n\n";
    
    // Argomento del perielio
    double omega_calc = varpi_calc - Omega_calc;
    if (omega_calc < 0) omega_calc += 2.0 * M_PI;
    std::cout << "✓ Argomento del perielio:\n";
    std::cout << "  Calcolata: ω = ϖ - Ω = " << (varpi_calc * 180.0 / M_PI) << "° - " 
              << (Omega_calc * 180.0 / M_PI) << "° = " << (omega_calc * 180.0 / M_PI) << "°\n";
    std::cout << "  Attesa:    ω = 110.28°\n";
    std::cout << "  Match: " << (fabs(omega_calc * 180.0 / M_PI - 110.28) < 1.0 ? "✓" : "✗") << "\n\n";
    
    // Anomalia media
    double M_calc = lambda - varpi_calc;
    while (M_calc < 0) M_calc += 2.0 * M_PI;
    while (M_calc >= 2.0 * M_PI) M_calc -= 2.0 * M_PI;
    std::cout << "✓ Anomalia media:\n";
    std::cout << "  Calcolata: M = λ - ϖ = " << lambda << " - " 
              << varpi_calc << " = " << M_calc << " rad = " << (M_calc * 180.0 / M_PI) << "°\n";
    std::cout << "  Attesa:    M = 25.45°\n";
    std::cout << "  Match: " << (fabs(M_calc * 180.0 / M_PI - 25.45) < 1.0 ? "✓" : "✗") << "\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // FASE 3: Risoluzione equazione di Keplero
    // ═══════════════════════════════════════════════════════════════
    std::cout << "FASE 3: Equazione di Keplero: E - e·sin(E) = M\n";
    std::cout << "───────────────────────────────────────────────────────────────\n";
    
    double E = M_calc + e_calc * sin(M_calc);
    for (int iter = 0; iter < 10; iter++) {
        E = M_calc + e_calc * sin(E);
    }
    
    std::cout << "✓ Anomalia eccentrica:\n";
    std::cout << "  Calcolata: E = " << E << " rad = " << (E * 180.0 / M_PI) << "°\n";
    std::cout << "  Attesa:    E = 0.4461 rad\n";
    std::cout << "  Match: " << (fabs(E - 0.4461) < 0.01 ? "✓" : "✗") << "\n\n";
    
    // Anomalia vera
    double sqrt_factor = sqrt(1.0 - e_calc*e_calc);
    double nu = 2.0 * atan2(sqrt_factor * sin(E/2.0), cos(E/2.0) - e_calc * cos(E/2.0) / sqrt_factor);
    // Versione alternativa più stabile:
    nu = 2.0 * atan2(sqrt(1.0 + e_calc) * sin(E/2.0), sqrt(1.0 - e_calc) * cos(E/2.0));
    
    std::cout << "✓ Anomalia vera:\n";
    std::cout << "  Calcolata: ν = " << nu << " rad = " << (nu * 180.0 / M_PI) << "°\n";
    std::cout << "  Attesa:    ν = 0.4593 rad = 26.32°\n";
    std::cout << "  Match: " << (fabs(nu - 0.4593) < 0.01 ? "✓" : "✗") << "\n\n";
    
    // Distanza
    double r = a * (1.0 - e_calc * cos(E));
    std::cout << "✓ Distanza dal Sole:\n";
    std::cout << "  Calcolata: r = a(1 - e·cos(E)) = " << std::setprecision(6) << r << " AU\n";
    std::cout << "  Attesa:    r = 3.0454 AU\n";
    std::cout << "  Match: " << (fabs(r - 3.0454) < 0.01 ? "✓" : "✗") << "\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // FASE 4: Posizione nel piano orbitale
    // ═══════════════════════════════════════════════════════════════
    std::cout << "FASE 4: Coordinate nel piano orbitale\n";
    std::cout << "───────────────────────────────────────────────────────────────\n";
    
    double x_orb = r * cos(nu);
    double y_orb = r * sin(nu);
    double z_orb = 0.0;
    
    std::cout << "✓ Posizione orbitale:\n";
    std::cout << "  Calcolata: (x_orb, y_orb, z_orb) = (" << std::setprecision(6) 
              << x_orb << ", " << y_orb << ", " << z_orb << ") AU\n";
    std::cout << "  Attesa:    (2.7393, 1.3337, 0) AU\n";
    std::cout << "  Match: " << (fabs(x_orb - 2.7393) < 0.01 && fabs(y_orb - 1.3337) < 0.01 ? "✓" : "✗") << "\n\n";
    
    // ═══════════════════════════════════════════════════════════════
    // FASE 5: Rotazioni Euleriane
    // ═══════════════════════════════════════════════════════════════
    std::cout << "FASE 5: Rotazioni Euleriane (piano orbitale → eclittico)\n";
    std::cout << "───────────────────────────────────────────────────────────────\n";
    
    double cos_Omega = cos(Omega_calc);
    double sin_Omega = sin(Omega_calc);
    double cos_omega = cos(omega_calc);
    double sin_omega = sin(omega_calc);
    double cos_i = cos(i_calc);
    double sin_i = sin(i_calc);
    
    std::cout << "Angoli:\n";
    std::cout << "  Ω = " << (Omega_calc * 180.0 / M_PI) << "° → cos(Ω) = " 
              << cos_Omega << ", sin(Ω) = " << sin_Omega << "\n";
    std::cout << "  ω = " << (omega_calc * 180.0 / M_PI) << "° → cos(ω) = " 
              << cos_omega << ", sin(ω) = " << sin_omega << "\n";
    std::cout << "  i = " << (i_calc * 180.0 / M_PI) << "° → cos(i) = " 
              << cos_i << ", sin(i) = " << sin_i << "\n\n";
    
    // Matrice rotazione
    double x_ecl = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * x_orb +
                   (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * y_orb;
    
    double y_ecl = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * x_orb +
                   (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * y_orb;
    
    double z_ecl = sin_omega * sin_i * x_orb + cos_omega * sin_i * y_orb;
    
    std::cout << "✓ Posizione eclittica:\n";
    std::cout << "  Calcolata: (x_ecl, y_ecl, z_ecl) = (" << std::setprecision(6) 
              << x_ecl << ", " << y_ecl << ", " << z_ecl << ") AU\n";
    std::cout << "  Attesa:    (-1.9899, -2.3073, 0.1093) AU\n";
    
    double delta_x = x_ecl - (-1.9899);
    double delta_y = y_ecl - (-2.3073);
    double delta_z = z_ecl - 0.1093;
    
    std::cout << "  Errori:    Δx = " << delta_x << ", Δy = " << delta_y 
              << ", Δz = " << delta_z << "\n";
    std::cout << "  Match: " << (fabs(delta_x) < 0.1 && fabs(delta_y) < 0.1 && fabs(delta_z) < 0.1 ? "✓" : "✗") << "\n\n";
    
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    if (fabs(delta_x) < 0.1 && fabs(delta_y) < 0.1 && fabs(delta_z) < 0.1) {
        std::cout << "✅ TUTTI I VALORI CORRETTI!\n";
    } else {
        std::cout << "⚠️  ERRORE RILEVATO NEI CALCOLI\n";
        std::cout << "Possibili cause:\n";
        std::cout << "1. Ordine delle rotazioni (verificare sequenza R_z, R_x, R_z)\n";
        std::cout << "2. Segni delle componenti\n";
        std::cout << "3. Unità di misura (radianti vs gradi)\n";
    }
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    
    return 0;
}
