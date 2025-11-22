/**
 * @file test_horizons_center.cpp
 * @brief Test diagnostico: verifica center Horizons (@sun vs @ssb)
 * 
 * Ipotesi: errore dovuto a heliocentric vs barycentric mismatch
 */

#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/orbit_propagator.h>
#include <iostream>
#include <iomanip>

using namespace ioccultcalc;

void printVector(const std::string& label, const Vector3D& v) {
    std::cout << label << ": (" 
              << std::fixed << std::setprecision(8)
              << v.x << ", " << v.y << ", " << v.z << ") AU\n";
}

int main() {
    JPLHorizonsClient horizons;
    JulianDate epoch;
    epoch.jd = 2460615.5;  // 2024-Jan-01
    
    std::cout << "\n=================================================================\n";
    std::cout << "TEST: Horizons CENTER (@sun vs @ssb)\n";
    std::cout << "=================================================================\n\n";
    
    std::cout << "Asteroide: 433 Eros\n";
    std::cout << "Epoca: JD " << epoch.jd << " (2024-Jan-01)\n\n";
    
    // Test 1: @sun (heliocentric)
    std::cout << "--- Test 1: CENTER='@sun' (heliocentric) ---\n";
    auto [pos_sun, vel_sun] = horizons.getStateVectors("433", epoch, "@sun");
    printVector("Posizione (@sun)", pos_sun);
    std::cout << "|r| = " << pos_sun.magnitude() << " AU\n\n";
    
    // Test 2: @ssb (barycentric)
    std::cout << "--- Test 2: CENTER='@0' (solar system barycenter) ---\n";
    auto [pos_ssb, vel_ssb] = horizons.getStateVectors("433", epoch, "@0");
    printVector("Posizione (@ssb)", pos_ssb);
    std::cout << "|r| = " << pos_ssb.magnitude() << " AU\n\n";
    
    // Differenza
    std::cout << "--- Differenza (@ssb - @sun) ---\n";
    Vector3D diff = pos_ssb - pos_sun;
    printVector("Δr", diff);
    std::cout << "|Δr| = " << diff.magnitude() << " AU\n";
    std::cout << "     = " << diff.magnitude() * 149597870.7 << " km\n\n";
    
    std::cout << "NOTA: La differenza è la posizione del Sole rispetto al baricentro.\n";
    std::cout << "      Tipicamente ~1000 km (influenza di Giove).\n\n";
    
    // Test propagazione con DE441
    std::cout << "=================================================================\n";
    std::cout << "TEST: DE441 frame (heliocentric vs barycentric)\n";
    std::cout << "=================================================================\n\n";
    
    PropagatorOptions opts;
    opts.integrator = IntegratorType::RK4;
    opts.stepSize = 0.1;
    opts.usePlanetaryPerturbations = true;
    opts.useRelativisticCorrections = false;
    
    OrbitPropagator propagator(opts);
    
    OrbitState state0;
    state0.position = pos_sun;  // Usa heliocentric
    state0.velocity = vel_sun;
    state0.epoch = epoch;
    
    JulianDate epoch_final;
    epoch_final.jd = epoch.jd + 100.0;
    
    std::cout << "Propagazione (100 giorni, heliocentric input)...\n";
    OrbitState state_final = propagator.propagate(state0, epoch_final);
    
    std::cout << "Posizione finale (IOccultCalc):\n";
    printVector("  r", state_final.position);
    std::cout << "  |r| = " << state_final.position.magnitude() << " AU\n\n";
    
    // Confronto con Horizons
    auto [pos_hor, vel_hor] = horizons.getStateVectors("433", epoch_final, "@sun");
    std::cout << "Posizione finale (Horizons @sun):\n";
    printVector("  r", pos_hor);
    std::cout << "  |r| = " << pos_hor.magnitude() << " AU\n\n";
    
    Vector3D error = state_final.position - pos_hor;
    std::cout << "ERRORE:\n";
    printVector("  Δr", error);
    std::cout << "  |Δr| = " << error.magnitude() << " AU\n";
    std::cout << "       = " << error.magnitude() * 149597870.7 << " km\n\n";
    
    std::cout << "Breakdown errore:\n";
    std::cout << "  ΔX = " << error.x * 149597870.7 << " km\n";
    std::cout << "  ΔY = " << error.y * 149597870.7 << " km\n";
    std::cout << "  ΔZ = " << error.z * 149597870.7 << " km (asse perpendicolare)\n\n";
    
    return 0;
}
