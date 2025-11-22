/**
 * @file test_vector_comparison.cpp
 * @brief Verifica: scarichiamo e usiamo correttamente i vettori Horizons?
 */

#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/orbit_propagator.h>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

int main() {
    JPLHorizonsClient horizons;
    
    JulianDate epoch0;
    epoch0.jd = 2460615.5;  // 2024-Jan-01
    
    JulianDate epoch1;
    epoch1.jd = 2460615.6;  // +0.1 giorni = 2.4 ore
    
    std::cout << "\n=================================================================\n";
    std::cout << "TEST: Verifica calcolo velocità da Horizons\n";
    std::cout << "=================================================================\n\n";
    
    // Scarica 2 posizioni ravvicinate
    auto [pos0, vel0_hor] = horizons.getStateVectors("433", epoch0);
    auto [pos1, vel1_hor] = horizons.getStateVectors("433", epoch1);
    
    std::cout << "Epoca 0: JD " << std::fixed << std::setprecision(6) << epoch0.jd << "\n";
    std::cout << "  pos = (" << std::setprecision(8) << pos0.x << ", " << pos0.y << ", " << pos0.z << ") AU\n";
    std::cout << "  vel = (" << vel0_hor.x << ", " << vel0_hor.y << ", " << vel0_hor.z << ") AU/day\n\n";
    
    std::cout << "Epoca 1: JD " << epoch1.jd << " (+0.1 giorni)\n";
    std::cout << "  pos = (" << pos1.x << ", " << pos1.y << ", " << pos1.z << ") AU\n";
    std::cout << "  vel = (" << vel1_hor.x << ", " << vel1_hor.y << ", " << vel1_hor.z << ") AU/day\n\n";
    
    // Calcola velocità numerica
    double dt = epoch1.jd - epoch0.jd;
    Vector3D vel_numeric = (pos1 - pos0) / dt;
    
    std::cout << "Velocità numerica (Δpos/Δt):\n";
    std::cout << "  vel = (" << vel_numeric.x << ", " << vel_numeric.y << ", " << vel_numeric.z << ") AU/day\n\n";
    
    // Confronto
    Vector3D vel_diff = vel0_hor - vel_numeric;
    std::cout << "Differenza (Horizons - numerica):\n";
    std::cout << "  Δvel = (" << vel_diff.x << ", " << vel_diff.y << ", " << vel_diff.z << ") AU/day\n";
    std::cout << "  |Δvel| = " << vel_diff.magnitude() << " AU/day\n";
    std::cout << "        = " << vel_diff.magnitude() * 149597870.7 << " km/day\n";
    std::cout << "        = " << vel_diff.magnitude() * 149597870.7 / 86400.0 << " km/s\n\n";
    
    if (vel_diff.magnitude() < 1e-6) {
        std::cout << "✅ Velocità coerente (<1e-6 AU/day)\n\n";
    } else {
        std::cout << "⚠️  Velocità NON coerente\n\n";
    }
    
    // Test propagazione breve con RK4
    std::cout << "=================================================================\n";
    std::cout << "TEST: Propagazione breve (0.1 giorni)\n";
    std::cout << "=================================================================\n\n";
    
    PropagatorOptions opts;
    opts.integrator = IntegratorType::RK4;
    opts.stepSize = 0.01;  // 0.01 giorni
    opts.usePlanetaryPerturbations = true;
    opts.useRelativisticCorrections = false;
    
    OrbitPropagator propagator(opts);
    
    OrbitState state0;
    state0.position = pos0;
    state0.velocity = vel0_hor;
    state0.epoch = epoch0;
    
    OrbitState state1 = propagator.propagate(state0, epoch1);
    
    std::cout << "Posizione finale (IOccultCalc):\n";
    std::cout << "  pos = (" << state1.position.x << ", " << state1.position.y << ", " << state1.position.z << ") AU\n\n";
    
    std::cout << "Posizione finale (Horizons):\n";
    std::cout << "  pos = (" << pos1.x << ", " << pos1.y << ", " << pos1.z << ") AU\n\n";
    
    Vector3D error = state1.position - pos1;
    std::cout << "ERRORE:\n";
    std::cout << "  Δpos = (" << error.x << ", " << error.y << ", " << error.z << ") AU\n";
    std::cout << "  |Δpos| = " << error.magnitude() << " AU\n";
    std::cout << "         = " << error.magnitude() * 149597870.7 << " km\n\n";
    
    std::cout << "Breakdown:\n";
    std::cout << "  ΔX = " << error.x * 149597870.7 << " km\n";
    std::cout << "  ΔY = " << error.y * 149597870.7 << " km\n";
    std::cout << "  ΔZ = " << error.z * 149597870.7 << " km\n\n";
    
    // Estrapolazione a 100 giorni
    double scale = 100.0 / 0.1;  // 1000×
    std::cout << "Estrapolazione a 100 giorni (×" << scale << "):\n";
    std::cout << "  Errore stimato = " << error.magnitude() * 149597870.7 * scale << " km\n\n";
    
    if (error.magnitude() * 149597870.7 < 1.0) {
        std::cout << "✅ Errore trascurabile (<1 km) → propagazione corretta\n";
    } else {
        std::cout << "⚠️  Errore significativo → problema sistematico\n";
    }
    
    return 0;
}
