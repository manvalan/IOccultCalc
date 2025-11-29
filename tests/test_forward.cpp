#include <iostream>
#include <iomanip>
#include <cmath>
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/types.h"
#include "ioccultcalc/orbital_elements.h"

using namespace ioccultcalc;

int main() {
    std::cout << "=== TEST PROPAGAZIONE IN AVANTI ===\n\n";
    
    PropagatorOptions opts;
    opts.usePlanetaryPerturbations = true;
    opts.stepSize = 0.1;
    opts.integrator = IntegratorType::RK4;
    OrbitPropagator propagator(opts);
    
    // Elementi Ceres epoca Nov 2026
    OrbitalElements kep;
    kep.epoch = JulianDate(2461000.5);
    kep.a = 2.7656157;
    kep.e = 0.0795763;
    kep.i = 10.58789 * M_PI / 180.0;
    kep.Omega = 80.24963 * M_PI / 180.0;
    kep.omega = 73.29974 * M_PI / 180.0;
    kep.M = 231.53975 * M_PI / 180.0;
    
    EquinoctialElements eq = kep.toEquinoctial();
    
    // Test 1: 30 giorni avanti
    JulianDate target30(kep.epoch.jd + 30);
    OrbitState result30 = propagator.propagate(eq, target30);
    
    constexpr double obliquity = 23.4392911 * M_PI / 180.0;
    double cos_eps = cos(obliquity);
    double sin_eps = sin(obliquity);
    
    double x30 = result30.position.x;
    double y30 = result30.position.y * cos_eps + result30.position.z * sin_eps;
    double z30 = -result30.position.y * sin_eps + result30.position.z * cos_eps;
    
    std::cout << "Dopo +30 giorni (eclittico):\n";
    std::cout << "  X = " << std::setprecision(6) << x30 << " AU\n";
    std::cout << "  Y = " << y30 << " AU\n";
    std::cout << "  Z = " << z30 << " AU\n";
    
    // Test 2: stesso tempo indietro
    JulianDate target_back(kep.epoch.jd);
    OrbitState result_back = propagator.propagate(eq, target_back);
    
    double xb = result_back.position.x;
    double yb = result_back.position.y * cos_eps + result_back.position.z * sin_eps;
    double zb = -result_back.position.y * sin_eps + result_back.position.z * cos_eps;
    
    std::cout << "\nAll'epoca (dt=0, round-trip test):\n";
    std::cout << "  X = " << xb << " AU\n";
    std::cout << "  Y = " << yb << " AU\n";
    std::cout << "  Z = " << zb << " AU\n";
    
    // Test 3: 30 giorni indietro
    JulianDate target_m30(kep.epoch.jd - 30);
    OrbitState result_m30 = propagator.propagate(eq, target_m30);
    
    double xm = result_m30.position.x;
    double ym = result_m30.position.y * cos_eps + result_m30.position.z * sin_eps;
    double zm = -result_m30.position.y * sin_eps + result_m30.position.z * cos_eps;
    
    std::cout << "\nDopo -30 giorni (eclittico):\n";
    std::cout << "  X = " << xm << " AU\n";
    std::cout << "  Y = " << ym << " AU\n";
    std::cout << "  Z = " << zm << " AU\n";
    
    return 0;
}
