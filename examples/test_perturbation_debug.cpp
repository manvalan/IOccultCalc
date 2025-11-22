/**
 * Test debug: verifica calcolo perturbazioni planetarie
 * Stampa accelerazioni di ogni pianeta per capire il problema
 */

#include <iostream>
#include <iomanip>
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/types.h"

using namespace ioccultcalc;

int main() {
    std::cout << "\n=== TEST DEBUG PERTURBAZIONI PLANETARIE ===\n\n";
    
    // Costanti
    const double GM_SUN = 2.959122082855911e-04;
    
    // Epoca e posizione asteroide
    JulianDate epoch(2460000.0);
    Vector3D pos_ast(2.8, 0.0, 0.0);  // AU
    Vector3D vel_ast(0.0, 0.01, 0.0);  // AU/day
    
    std::cout << "Asteroide:\n";
    std::cout << "  Posizione: (" << pos_ast.x << ", " << pos_ast.y << ", " << pos_ast.z << ") AU\n";
    std::cout << "  Distanza: " << pos_ast.magnitude() << " AU\n\n";
    
    // Crea propagatore con perturbazioni
    PropagatorOptions opts;
    opts.usePlanetaryPerturbations = true;
    opts.useRelativisticCorrections = false;
    
    OrbitPropagator propagator(opts);
    
    // Calcola accelerazione totale
    Vector3D acc_total = propagator.computeAcceleration(epoch, pos_ast, vel_ast);
    
    std::cout << "Accelerazione TOTALE:\n";
    std::cout << "  a = (" << std::scientific << std::setprecision(12)
              << acc_total.x << ", " << acc_total.y << ", " << acc_total.z << ") AU/day²\n";
    std::cout << "  |a| = " << acc_total.magnitude() << " AU/day²\n\n";
    
    // Accelerazione solo Sole
    double r = pos_ast.magnitude();
    double a_sun = GM_SUN / (r * r);
    std::cout << "Accelerazione SOLE (teorica):\n";
    std::cout << "  |a| = " << a_sun << " AU/day²\n";
    std::cout << "  Direzione: verso (0,0,0)\n\n";
    
    // Differenza = perturbazioni
    Vector3D a_sun_vec(-GM_SUN * pos_ast.x / (r*r*r),
                       -GM_SUN * pos_ast.y / (r*r*r),
                       -GM_SUN * pos_ast.z / (r*r*r));
    
    Vector3D a_pert = acc_total - a_sun_vec;
    
    std::cout << "Accelerazione PERTURBAZIONI (totale - sole):\n";
    std::cout << "  a = (" << a_pert.x << ", " << a_pert.y << ", " << a_pert.z << ") AU/day²\n";
    std::cout << "  |a| = " << a_pert.magnitude() << " AU/day²\n";
    std::cout << "  Rapporto pert/sole: " << (a_pert.magnitude() / a_sun) << "\n\n";
    
    // Test senza perturbazioni
    std::cout << "=== TEST SENZA PERTURBAZIONI ===\n\n";
    
    PropagatorOptions opts_nop;
    opts_nop.usePlanetaryPerturbations = false;
    opts_nop.useRelativisticCorrections = false;
    
    OrbitPropagator propagator_nop(opts_nop);
    
    Vector3D acc_nop = propagator_nop.computeAcceleration(epoch, pos_ast, vel_ast);
    
    std::cout << "Accelerazione (solo Sole):\n";
    std::cout << "  a = (" << acc_nop.x << ", " << acc_nop.y << ", " << acc_nop.z << ") AU/day²\n";
    std::cout << "  |a| = " << acc_nop.magnitude() << " AU/day²\n";
    std::cout << "  Errore vs teorica: " << std::abs(acc_nop.magnitude() - a_sun) / a_sun << "\n\n";
    
    // Propagazione breve con/senza perturbazioni
    std::cout << "=== PROPAGAZIONE 1 GIORNO ===\n\n";
    
    JulianDate epoch1(epoch.jd + 1.0);
    OrbitState state0(epoch, pos_ast, vel_ast);
    
    // Con perturbazioni
    opts.integrator = IntegratorType::RK4;
    opts.stepSize = 0.01;
    OrbitPropagator prop_pert(opts);
    OrbitState state_pert = prop_pert.propagate(state0, epoch1);
    
    // Senza perturbazioni
    opts_nop.integrator = IntegratorType::RK4;
    opts_nop.stepSize = 0.01;
    OrbitPropagator prop_nopert(opts_nop);
    OrbitState state_nop = prop_nopert.propagate(state0, epoch1);
    
    std::cout << "CON perturbazioni:\n";
    std::cout << "  r = (" << std::setprecision(8)
              << state_pert.position.x << ", "
              << state_pert.position.y << ", "
              << state_pert.position.z << ") AU\n";
    std::cout << "  |r| = " << state_pert.position.magnitude() << " AU\n";
    
    std::cout << "\nSENZA perturbazioni:\n";
    std::cout << "  r = (" << state_nop.position.x << ", "
              << state_nop.position.y << ", "
              << state_nop.position.z << ") AU\n";
    std::cout << "  |r| = " << state_nop.position.magnitude() << " AU\n";
    
    Vector3D diff = state_pert.position - state_nop.position;
    std::cout << "\nDifferenza:\n";
    std::cout << "  Δr = (" << diff.x << ", " << diff.y << ", " << diff.z << ") AU\n";
    std::cout << "  |Δr| = " << diff.magnitude() << " AU = " 
              << diff.magnitude() * 149597870.7 << " km\n";
    std::cout << "  Δr/r = " << (diff.magnitude() / state_nop.position.magnitude()) << "\n";
    
    std::cout << "\n=== Test completato ===\n\n";
    
    return 0;
}
