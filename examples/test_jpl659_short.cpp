#include <iostream>
#include <iomanip>
#include <cmath>
#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/orbit_propagator.h>
#include <ioccultcalc/orbital_elements.h>
#include <ioccultcalc/types.h>

using namespace ioccultcalc;

int main() {
    std::cout << "=================================================================\n";
    std::cout << "TEST: PROPAGAZIONE BREVE CON JPL#659\n";
    std::cout << "=================================================================\n\n";
    
    std::cout << "Test: Propagare solo 10 giorni da epoch 2004\n\n";
    
    // Elementi JPL#659 @ 2004-Nov-02
    OrbitalElements elem;
    elem.designation = "433";
    elem.name = "Eros";
    elem.epoch.jd = 2453311.5;
    elem.a = 1.458269315549994;
    elem.e = 0.2228078944584026;
    elem.i = 10.8291838260782 * M_PI / 180.0;
    elem.Omega = 304.4010273379536 * M_PI / 180.0;
    elem.omega = 178.665326776373 * M_PI / 180.0;
    elem.M = 326.37047603642 * M_PI / 180.0;
    
    auto eq_elem = elem.toEquinoctial();
    
    PropagatorOptions opts;
    opts.integrator = IntegratorType::RK4;
    opts.stepSize = 0.1;
    
    OrbitPropagator prop(opts);
    auto state_0 = prop.elementsToState(eq_elem);
    
    std::cout << "Stato iniziale @ 2004-Nov-02:\n";
    std::cout << "  IOccultCalc: (" << std::setprecision(6) 
              << state_0.position.x << ", "
              << state_0.position.y << ", "
              << state_0.position.z << ") AU\n";
    std::cout << "  Horizons:    (0.373974, 1.144247, 0.182689) AU\n\n";
    
    double dx0 = (state_0.position.x - 0.373974) * 149597870.7;
    double dy0 = (state_0.position.y - 1.144247) * 149597870.7;
    double dz0 = (state_0.position.z - 0.182689) * 149597870.7;
    double diff0 = std::sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
    
    std::cout << "Differenza @ T=0: " << std::setprecision(1) << diff0 << " km\n\n";
    
    // Propaga 10 giorni
    JulianDate target;
    target.jd = elem.epoch.jd + 10.0;
    
    std::cout << "Propagando 10 giorni...\n\n";
    auto state_10 = prop.propagate(state_0, target);
    
    std::cout << "Stato @ +10 giorni:\n";
    std::cout << "  IOccultCalc: (" << std::setprecision(6)
              << state_10.position.x << ", "
              << state_10.position.y << ", "
              << state_10.position.z << ") AU\n";
    
    // Query Horizons per +10d
    JPLHorizonsClient horizons;
    auto [r_hor, v_hor] = horizons.getStateVectors("433", target, "@sun");
    
    std::cout << "  Horizons:    (" << r_hor.x << ", " << r_hor.y << ", " << r_hor.z << ") AU\n\n";
    
    double dx = (state_10.position.x - r_hor.x) * 149597870.7;
    double dy = (state_10.position.y - r_hor.y) * 149597870.7;
    double dz = (state_10.position.z - r_hor.z) * 149597870.7;
    double diff = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    std::cout << "Differenza @ +10d: " << std::setprecision(1) << diff << " km\n";
    std::cout << "  ΔX = " << dx << " km\n";
    std::cout << "  ΔY = " << dy << " km\n";
    std::cout << "  ΔZ = " << dz << " km\n\n";
    
    if (diff < 100) {
        std::cout << "✅ ECCELLENTE! IOccultCalc accurato su 10 giorni!\n";
    } else if (diff < 500) {
        std::cout << "✓ BUONO! Errore accettabile.\n";
    } else {
        std::cout << "⚠ Ancora " << diff << " km di errore su solo 10 giorni.\n";
    }
    
    return 0;
}
