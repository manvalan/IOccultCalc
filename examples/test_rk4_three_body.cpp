/**
 * Test RK4 con problema a 3 corpi
 * Verifica propagazione con Sole + Giove + asteroide
 * Confronto energia e conservazione momento angolare
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/types.h"

using namespace ioccultcalc;

void printSeparator(char c = '=') {
    std::cout << std::string(70, c) << "\n";
}

// Calcola energia totale del sistema (asteroide nel campo Sole+Giove)
double totalEnergy(const Vector3D& pos_ast, const Vector3D& vel_ast,
                  const Vector3D& pos_jup, double GM_sun, double GM_jup) {
    // Energia cinetica asteroide
    double v2 = vel_ast.x*vel_ast.x + vel_ast.y*vel_ast.y + vel_ast.z*vel_ast.z;
    double KE = 0.5 * v2;
    
    // Energia potenziale Sole-asteroide
    double r_sun = pos_ast.magnitude();
    double PE_sun = -GM_sun / r_sun;
    
    // Energia potenziale Giove-asteroide
    Vector3D d_jup = pos_ast - pos_jup;
    double r_jup = d_jup.magnitude();
    double PE_jup = -GM_jup / r_jup;
    
    return KE + PE_sun + PE_jup;
}

int main() {
    std::cout << "\n";
    printSeparator();
    std::cout << "TEST RK4: PROBLEMA 3 CORPI (Sole + Giove + Asteroide)\n";
    printSeparator();
    std::cout << "\n";
    
    // Costanti gravitazionali (AU³/day²)
    const double GM_SUN = 2.959122082855911e-04;
    const double GM_JUPITER = 2.825345909524226e-07;  // ~1/1047 della massa solare
    
    // Epoca iniziale
    JulianDate epoch0(2460000.0);  // Epoca arbitraria
    
    // Posizione Giove (orbita circolare a 5.2 AU)
    double a_jup = 5.2;  // AU
    double v_jup = std::sqrt(GM_SUN / a_jup);  // Velocità circolare
    Vector3D pos_jup0(a_jup, 0.0, 0.0);
    Vector3D vel_jup0(0.0, v_jup, 0.0);
    
    std::cout << "CONFIGURAZIONE INIZIALE\n";
    printSeparator('-');
    std::cout << "\nGiove (orbita circolare):\n";
    std::cout << "  a = " << a_jup << " AU\n";
    std::cout << "  r = (" << pos_jup0.x << ", " << pos_jup0.y << ", " << pos_jup0.z << ") AU\n";
    std::cout << "  v = (" << vel_jup0.x << ", " << vel_jup0.y << ", " << vel_jup0.z << ") AU/day\n";
    std::cout << "  |v| = " << vel_jup0.magnitude() << " AU/day (teorica: " << v_jup << ")\n";
    std::cout << "  periodo = " << 2*M_PI*std::sqrt(a_jup*a_jup*a_jup/GM_SUN) << " giorni (~11.86 anni)\n";
    
    // Asteroide tra Marte e Giove (a = 2.8 AU, main belt)
    double a_ast = 2.8;  // AU
    double e_ast = 0.05;  // Bassa eccentricità
    
    // Posizione all'afelio (maggiore distanza)
    double r_ast = a_ast * (1.0 + e_ast);
    Vector3D pos_ast0(r_ast, 0.0, 0.0);
    
    // Velocità all'afelio (perpendicolare al raggio)
    double h = std::sqrt(GM_SUN * a_ast * (1.0 - e_ast*e_ast));  // Momento angolare specifico
    double v_ast = h / r_ast;
    Vector3D vel_ast0(0.0, v_ast, 0.0);
    
    std::cout << "\nAsteroide (orbita ellittica, main belt):\n";
    std::cout << "  a = " << a_ast << " AU\n";
    std::cout << "  e = " << e_ast << "\n";
    std::cout << "  r = (" << pos_ast0.x << ", " << pos_ast0.y << ", " << pos_ast0.z << ") AU (afelio)\n";
    std::cout << "  v = (" << vel_ast0.x << ", " << vel_ast0.y << ", " << vel_ast0.z << ") AU/day\n";
    std::cout << "  |v| = " << vel_ast0.magnitude() << " AU/day\n";
    std::cout << "  periodo = " << 2*M_PI*std::sqrt(a_ast*a_ast*a_ast/GM_SUN) << " giorni (~4.68 anni)\n";
    
    // Energia iniziale
    double E0 = totalEnergy(pos_ast0, vel_ast0, pos_jup0, GM_SUN, GM_JUPITER);
    Vector3D h0 = pos_ast0.cross(vel_ast0);
    
    std::cout << "\nGrandezze conservate:\n";
    std::cout << "  Energia totale: " << std::scientific << std::setprecision(12) << E0 << " AU²/day²\n";
    std::cout << "  |h| asteroide:  " << h0.magnitude() << " AU²/day\n";
    
    // ================================================================
    // TEST 1: Propagazione con RK4 (SOLO GIOVE, no altri pianeti)
    // ================================================================
    std::cout << "\n\n";
    printSeparator();
    std::cout << "TEST 1: Propagazione con RK4 (Sole + Giove)\n";
    printSeparator();
    
    PropagatorOptions options;
    options.integrator = IntegratorType::RK4;
    options.stepSize = 0.01;  // 0.01 giorni = 14.4 minuti
    options.usePlanetaryPerturbations = true;   // Abilita Giove
    options.useRelativisticCorrections = false;
    
    OrbitPropagator propagator(options);
    
    // Propaga per 1/4 periodo asteroide (~1.17 anni)
    double t_prop = 365.25 * 1.17;  // giorni
    JulianDate targetEpoch(epoch0.jd + t_prop);
    
    std::cout << "\nPropagazione per " << std::fixed << std::setprecision(2) 
              << t_prop << " giorni (~1.17 anni)\n";
    std::cout << "Stepsize RK4: " << options.stepSize << " giorni\n";
    std::cout << "Numero step atteso: ~" << (int)(t_prop / options.stepSize) << "\n\n";
    
    OrbitState initialState(epoch0, pos_ast0, vel_ast0);
    
    std::cout << "Propagazione in corso...\n";
    auto start = std::chrono::high_resolution_clock::now();
    
    OrbitState finalState = propagator.propagate(initialState, targetEpoch);
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    std::cout << "✓ Propagazione completata in " << elapsed.count() << " secondi!\n\n";
    
    // Stato finale
    std::cout << "Stato finale asteroide:\n";
    std::cout << "  r = (" << std::setprecision(8) 
              << finalState.position.x << ", " 
              << finalState.position.y << ", " 
              << finalState.position.z << ") AU\n";
    std::cout << "  |r| = " << finalState.position.magnitude() << " AU\n";
    std::cout << "  v = (" << finalState.velocity.x << ", " 
              << finalState.velocity.y << ", " 
              << finalState.velocity.z << ") AU/day\n";
    std::cout << "  |v| = " << finalState.velocity.magnitude() << " AU/day\n";
    
    // Posizione Giove al tempo finale (orbita circolare)
    double omega_jup = v_jup / a_jup;  // Velocità angolare
    double theta_jup = omega_jup * t_prop;
    Vector3D pos_jup_final(
        a_jup * std::cos(theta_jup),
        a_jup * std::sin(theta_jup),
        0.0
    );
    
    // Energia finale
    double Ef = totalEnergy(finalState.position, finalState.velocity, 
                           pos_jup_final, GM_SUN, GM_JUPITER);
    Vector3D hf = finalState.position.cross(finalState.velocity);
    
    std::cout << "\nGrandezze finali:\n";
    std::cout << "  Energia: " << std::scientific << std::setprecision(12) << Ef << " AU²/day²\n";
    std::cout << "  |h|:     " << hf.magnitude() << " AU²/day\n";
    
    // Errori
    double dE = std::abs(Ef - E0);
    double dE_rel = dE / std::abs(E0);
    double dh = (hf - h0).magnitude();
    double dh_rel = dh / h0.magnitude();
    
    std::cout << "\n";
    printSeparator('-');
    std::cout << "ANALISI CONSERVAZIONE\n";
    printSeparator('-');
    std::cout << "ΔE (assoluto): " << std::scientific << dE << " AU²/day²\n";
    std::cout << "ΔE (relativo): " << dE_rel << "\n";
    std::cout << "Δh (assoluto): " << dh << " AU²/day\n";
    std::cout << "Δh (relativo): " << dh_rel << "\n";
    
    // Criteri validazione (più rilassati per 3 corpi)
    // Criteri più rilassati per 3 corpi (perturbazioni reali)
    bool energyOK = dE_rel < 1e-3;    // Energia entro 0.1%
    bool momentumOK = dh_rel < 1e-3;  // Momento angolare entro 0.1%
    
    std::cout << "\n";
    printSeparator('-');
    std::cout << "RISULTATO TEST\n";
    printSeparator('-');
    std::cout << (energyOK ? "✓" : "✗") << " Energia conservata (ΔE/E < 10⁻³)\n";
    std::cout << (momentumOK ? "✓" : "✗") << " Momento angolare conservato (Δh/h < 10⁻³)\n";
    
    if (!energyOK || !momentumOK) {
        std::cout << "\n✗ FAIL: Errore troppo grande (bug nel propagatore?)\n";
    } else {
        std::cout << "\n✓ TEST SUPERATO: RK4 conserva energia nel problema 3 corpi\n";
    }
    
    // ================================================================
    // INFORMAZIONI DIAGNOSTICHE
    // ================================================================
    std::cout << "\n\n";
    printSeparator();
    std::cout << "DIAGNOSTICA\n";
    printSeparator();
    
    // Effetto perturbativo di Giove
    Vector3D d_jup = pos_ast0 - pos_jup0;
    double r_jup_ast = d_jup.magnitude();
    double a_pert = GM_JUPITER / (r_jup_ast * r_jup_ast);
    double a_sole = GM_SUN / (pos_ast0.magnitude() * pos_ast0.magnitude());
    
    std::cout << "\nPerturbazione di Giove:\n";
    std::cout << "  Distanza iniziale Giove-asteroide: " << r_jup_ast << " AU\n";
    std::cout << "  a_Giove / a_Sole: " << std::scientific << (a_pert / a_sole) << "\n";
    std::cout << "  Rapporto masse (M_Jup / M_Sun): " << (GM_JUPITER / GM_SUN) << "\n";
    std::cout << "  → Perturbazione ~" << (a_pert / a_sole * 100) << "% dell'attrazione solare\n";
    
    // Variazione periodo dovuta a Giove
    double delta_T = t_prop * (dE_rel);  // Stima approssimativa
    std::cout << "\nVariazione periodo stimata: " << std::fixed << std::setprecision(3) 
              << delta_T << " giorni\n";
    
    std::cout << "\n";
    printSeparator();
    std::cout << "Test completato!\n";
    printSeparator();
    std::cout << "\n";
    
    return (energyOK && momentumOK) ? 0 : 1;
}
