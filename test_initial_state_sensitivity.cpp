/**
 * @file test_initial_state_sensitivity.cpp
 * @brief Test sensibilità a piccole variazioni stato iniziale
 * 
 * Obiettivo: Verificare se l'errore di 2000 km @ 100d
 * può essere spiegato da incertezza ~10-20 km nello stato iniziale
 */

#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/orbit_propagator.h>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

void printSep() { std::cout << std::string(70, '=') << "\n"; }

int main() {
    std::cout << "\n";
    printSep();
    std::cout << "TEST: Sensibilità a Perturbazioni Stato Iniziale\n";
    printSep();
    std::cout << "\n";
    
    // Setup
    JPLHorizonsClient horizons;
    JulianDate epoch0;
    epoch0.jd = 2460615.5;
    
    JulianDate epochF;
    epochF.jd = epoch0.jd + 100.0;
    
    std::cout << "Configurazione:\n";
    std::cout << "  • Asteroide: 433 Eros\n";
    std::cout << "  • Propagazione: 100 giorni\n";
    std::cout << "  • Perturbazioni: Full (pianeti + asteroidi + rel)\n\n";
    
    // Stato nominale
    auto [pos0, vel0] = horizons.getStateVectors("433", epoch0);
    auto [posF_ref, velF_ref] = horizons.getStateVectors("433", epochF);
    
    PropagatorOptions opts;
    opts.integrator = IntegratorType::RK4;
    opts.stepSize = 0.1;
    opts.usePlanetaryPerturbations = true;
    opts.useRelativisticCorrections = true;
    
    OrbitPropagator propagator(opts);
    
    // Test 1: Nominale
    std::cout << "Test 1: Stato nominale\n";
    OrbitState state0;
    state0.position = pos0;
    state0.velocity = vel0;
    state0.epoch = epoch0;
    
    OrbitState stateF = propagator.propagate(state0, epochF);
    Vector3D err_nominal = stateF.position - posF_ref;
    double err_nominal_km = err_nominal.magnitude() * 149597870.7;
    
    std::cout << "  Errore: " << std::fixed << std::setprecision(1) 
              << err_nominal_km << " km\n\n";
    
    // Test 2: +10 km in posizione
    std::cout << "Test 2: Perturbazione +10 km in posizione X\n";
    OrbitState state_pert = state0;
    state_pert.position.x += 10.0 / 149597870.7;  // +10 km in AU
    
    OrbitState stateF_pert = propagator.propagate(state_pert, epochF);
    Vector3D err_pert = stateF_pert.position - posF_ref;
    double err_pert_km = err_pert.magnitude() * 149597870.7;
    
    std::cout << "  Errore: " << err_pert_km << " km\n";
    std::cout << "  Amplificazione: " << err_pert_km / 10.0 << "×\n\n";
    
    // Test 3: +1 km/s in velocità
    std::cout << "Test 3: Perturbazione +1 m/s in velocità X\n";
    state_pert = state0;
    state_pert.velocity.x += 0.001 / 149597870.7 * 86400.0;  // +1 m/s → AU/day
    
    stateF_pert = propagator.propagate(state_pert, epochF);
    err_pert = stateF_pert.position - posF_ref;
    err_pert_km = err_pert.magnitude() * 149597870.7;
    
    std::cout << "  Errore posizione finale: " << err_pert_km << " km\n";
    std::cout << "  Amplificazione: " << err_pert_km / 0.001 << "× (m/s → km @ 100d)\n\n";
    
    // Test 4: Combinazione realistica
    std::cout << "Test 4: Perturbazione realistica (±5 km pos, ±0.5 m/s vel)\n";
    state_pert = state0;
    // Perturbazione random-like
    state_pert.position.x += 3.0 / 149597870.7;
    state_pert.position.y += -4.0 / 149597870.7;
    state_pert.position.z += 2.0 / 149597870.7;
    state_pert.velocity.x += 0.0003 / 149597870.7 * 86400.0;
    state_pert.velocity.y += -0.0004 / 149597870.7 * 86400.0;
    
    stateF_pert = propagator.propagate(state_pert, epochF);
    err_pert = stateF_pert.position - posF_ref;
    err_pert_km = err_pert.magnitude() * 149597870.7;
    
    std::cout << "  Errore finale: " << err_pert_km << " km\n\n";
    
    // Analisi
    printSep();
    std::cout << "ANALISI\n";
    printSep();
    std::cout << "\n";
    
    std::cout << "Errore nominale (Horizons vs IOccultCalc): " 
              << err_nominal_km << " km\n\n";
    
    std::cout << "Se l'errore fosse solo da incertezza elementi:\n";
    std::cout << "  • Incertezza tipica Horizons: ~5-10 km\n";
    std::cout << "  • Amplificazione @ 100d: ~10-20×\n";
    std::cout << "  • Errore atteso: 50-200 km\n\n";
    
    std::cout << "Errore osservato: " << err_nominal_km << " km\n";
    std::cout << "Fattore discrepanza: " << err_nominal_km / 100.0 << "×\n\n";
    
    std::cout << "CONCLUSIONE:\n";
    if (err_nominal_km > 1000) {
        std::cout << "  L'errore di " << err_nominal_km 
                  << " km NON può essere spiegato\n";
        std::cout << "  solo da incertezza elementi iniziali.\n";
        std::cout << "  Altre cause dominanti:\n";
        std::cout << "    - Modello perturbazioni semplificato\n";
        std::cout << "    - Effetti non modellati (J2, J3, tide, ...)\n";
        std::cout << "    - Differenze algoritmo integrazione JPL vs RK4\n";
    } else {
        std::cout << "  L'errore potrebbe essere spiegato da incertezza elementi.\n";
    }
    
    std::cout << "\n";
    printSep();
    std::cout << "\n";
    
    return 0;
}
