/**
 * Test di validazione del propagatore orbitale
 * Verifica:
 * 1. Conservazione dell'energia (Hamiltoniana)
 * 2. Conservazione del momento angolare
 * 3. Confronto con soluzione analitica (problema a 2 corpi)
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/types.h"

using namespace ioccultcalc;

void printSeparator(char c = '=') {
    std::cout << std::string(70, c) << "\n";
}

// Calcola energia orbitale specifica: E = v²/2 - GM/r
double orbitalEnergy(const Vector3D& pos, const Vector3D& vel, double GM) {
    double r = pos.magnitude();
    double v2 = vel.x*vel.x + vel.y*vel.y + vel.z*vel.z;
    return v2/2.0 - GM/r;
}

// Calcola momento angolare: h = r × v
Vector3D angularMomentum(const Vector3D& pos, const Vector3D& vel) {
    return pos.cross(vel);
}

int main() {
    std::cout << "\n";
    printSeparator();
    std::cout << "TEST VALIDAZIONE PROPAGATORE ORBITALE\n";
    printSeparator();
    std::cout << "\n";
    
    // Costante gravitazionale del Sole (AU³/day²)
    const double GM_SUN = 2.959122082855911e-04;
    
    // ================================================================
    // TEST 1: Problema a 2 corpi puro (solo Sole)
    // ================================================================
    std::cout << "TEST 1: Propagazione 2-corpi (solo Sole, senza perturbazioni)\n";
    printSeparator('-');
    
    // Stato iniziale: orbita circolare a 1 AU
    JulianDate epoch0(2451545.0);
    Vector3D pos0(1.0, 0.0, 0.0);  // AU
    double v_circular = std::sqrt(GM_SUN / 1.0);  // Velocità circolare
    Vector3D vel0(0.0, v_circular, 0.0);  // AU/day
    
    std::cout << "Stato iniziale (orbita circolare a 1 AU):\n";
    std::cout << "  r = (" << pos0.x << ", " << pos0.y << ", " << pos0.z << ") AU\n";
    std::cout << "  v = (" << vel0.x << ", " << vel0.y << ", " << vel0.z << ") AU/day\n";
    std::cout << "  |v| = " << vel0.magnitude() << " AU/day\n";
    std::cout << "  v_teorica = " << v_circular << " AU/day\n";
    
    // Energia e momento angolare iniziali
    double E0 = orbitalEnergy(pos0, vel0, GM_SUN);
    Vector3D h0 = angularMomentum(pos0, vel0);
    
    std::cout << "\nGrandezze conservate (teoricamente):\n";
    std::cout << "  Energia: " << std::scientific << std::setprecision(12) << E0 << " AU²/day²\n";
    std::cout << "  |h|:     " << h0.magnitude() << " AU²/day\n";
    
    // Configurazione propagatore: SOLO SOLE (nessuna perturbazione)
    PropagatorOptions options;
    options.integrator = IntegratorType::RA15;  // ← Torno a RA15 con fix iterazioni
    options.tolerance = 1e-12;  
    options.stepSize = 1.0;     // ← Step iniziale 1 giorno (RA15 adatta automaticamente)
    options.usePlanetaryPerturbations = false;  // ← NESSUN PIANETA
    options.useRelativisticCorrections = false; // ← NESSUNA RELATIVITÀ
    
    OrbitPropagator propagator(options);
    
    // Test: verifica accelerazione iniziale
    Vector3D acc0 = propagator.computeAcceleration(epoch0, pos0, vel0);
    double acc0_mag = acc0.magnitude();
    double acc0_teorica = GM_SUN / 1.0;  // GM/r² a 1 AU
    
    std::cout << "\nTest accelerazione:\n";
    std::cout << "  a_calcolata = " << std::scientific << std::setprecision(12) 
              << acc0_mag << " AU/day²\n";
    std::cout << "  a_teorica   = " << acc0_teorica << " AU/day²\n";
    std::cout << "  errore rel  = " << std::abs(acc0_mag - acc0_teorica) / acc0_teorica << "\n";
    
    // Propaga per 1 periodo orbitale (~365 giorni)
    std::cout << "\nPropagazione per 1 anno (1 periodo orbitale)...\n";
    
    OrbitState initialState(epoch0, pos0, vel0);
    JulianDate targetEpoch(epoch0.jd + 365.25);
    
    OrbitState finalState = propagator.propagate(initialState, targetEpoch);
    
    std::cout << "✓ Propagazione completata!\n\n";
    
    // Stato finale
    std::cout << "Stato finale:\n";
    std::cout << "  r = (" << finalState.position.x << ", " 
              << finalState.position.y << ", " 
              << finalState.position.z << ") AU\n";
    std::cout << "  |r| = " << finalState.position.magnitude() << " AU\n";
    std::cout << "  v = (" << finalState.velocity.x << ", " 
              << finalState.velocity.y << ", " 
              << finalState.velocity.z << ") AU/day\n";
    std::cout << "  |v| = " << finalState.velocity.magnitude() << " AU/day\n";
    
    // Energia e momento angolare finali
    double Ef = orbitalEnergy(finalState.position, finalState.velocity, GM_SUN);
    Vector3D hf = angularMomentum(finalState.position, finalState.velocity);
    
    std::cout << "\nGrandezze finali:\n";
    std::cout << "  Energia: " << std::scientific << std::setprecision(12) << Ef << " AU²/day²\n";
    std::cout << "  |h|:     " << hf.magnitude() << " AU²/day\n";
    
    // Errori relativi
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
    
    // Verifica ritorno alla posizione iniziale (orbita circolare)
    double dr = (finalState.position - pos0).magnitude();
    std::cout << "\nDistanza dalla posizione iniziale: " << std::fixed << std::setprecision(6) 
              << dr << " AU = " << dr * 149597870.7 << " km\n";
    
    // Criteri di validazione
    bool energyOK = dE_rel < 1e-9;
    bool momentumOK = dh_rel < 1e-9;
    bool positionOK = dr < 1e-6;  // < 150 km
    
    std::cout << "\n";
    printSeparator('-');
    std::cout << "RISULTATO TEST 1:\n";
    printSeparator('-');
    std::cout << (energyOK ? "✓" : "✗") << " Energia conservata (ΔE/E < 10⁻⁹)\n";
    std::cout << (momentumOK ? "✓" : "✗") << " Momento angolare conservato (Δh/h < 10⁻⁹)\n";
    std::cout << (positionOK ? "✓" : "✗") << " Posizione chiusa (dr < 150 km)\n";
    
    if (!energyOK || !momentumOK || !positionOK) {
        std::cout << "\n⚠ FALLIMENTO: Il propagatore NON conserva le grandezze fisiche!\n";
        std::cout << "   Problema nell'integratore o nel force model.\n";
        return 1;
    }
    
    std::cout << "\n✓ TEST 1 SUPERATO: Problema 2-corpi funziona correttamente\n";
    
    // ================================================================
    // TEST 2: Confronto periodo orbitale
    // ================================================================
    std::cout << "\n\n";
    printSeparator();
    std::cout << "TEST 2: Verifica periodo orbitale\n";
    printSeparator();
    
    // Calcola periodo da elementi iniziali
    double a = 1.0;  // semi-asse maggiore
    double T_teorico = 2.0 * M_PI * std::sqrt(a*a*a / GM_SUN);
    
    std::cout << "\nPeriodo orbitale teorico (Keplero): " 
              << std::fixed << std::setprecision(6) << T_teorico << " giorni\n";
    
    // Propaga per 10 periodi
    double propagationTime = 10.0 * T_teorico;
    std::cout << "Propagazione per 10 periodi (" << propagationTime << " giorni)...\n";
    
    JulianDate epoch10(epoch0.jd + propagationTime);
    OrbitState state10 = propagator.propagate(initialState, epoch10);
    
    double dr10 = (state10.position - pos0).magnitude();
    double E10 = orbitalEnergy(state10.position, state10.velocity, GM_SUN);
    double dE10_rel = std::abs(E10 - E0) / std::abs(E0);
    
    std::cout << "✓ Propagazione completata!\n\n";
    std::cout << "Distanza dalla posizione iniziale: " << dr10 << " AU = " 
              << dr10 * 149597870.7 << " km\n";
    std::cout << "Errore energia relativo: " << std::scientific << dE10_rel << "\n";
    
    bool longTermOK = (dr10 < 1e-5) && (dE10_rel < 1e-8);
    
    std::cout << "\n";
    printSeparator('-');
    std::cout << "RISULTATO TEST 2:\n";
    printSeparator('-');
    std::cout << (longTermOK ? "✓" : "✗") << " Stabilità a lungo termine (10 periodi)\n";
    
    if (!longTermOK) {
        std::cout << "\n⚠ FALLIMENTO: Deriva a lungo termine!\n";
        return 1;
    }
    
    std::cout << "\n✓ TEST 2 SUPERATO: Integrazione stabile nel lungo periodo\n";
    
    // ================================================================
    // CONCLUSIONE
    // ================================================================
    std::cout << "\n\n";
    printSeparator();
    std::cout << "TUTTI I TEST SUPERATI ✓\n";
    std::cout << "Il propagatore conserva energia e momento angolare correttamente.\n";
    printSeparator();
    std::cout << "\n";
    
    return 0;
}
