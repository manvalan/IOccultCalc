/**
 * Test: OrbitPropagator con perturbazioni per propagazione all'indietro
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/types.h"
#include "ioccultcalc/orbital_elements.h"

using namespace ioccultcalc;

const double RAD2DEG = 180.0 / M_PI;
const double DEG2RAD = M_PI / 180.0;
const double AU_KM = 149597870.7;

int main() {
    std::cout << "===========================================\n";
    std::cout << " TEST ORBIT PROPAGATOR - PROPAGAZIONE INDIETRO\n";
    std::cout << "===========================================\n\n";
    
    // Inizializza OrbitPropagator con opzioni
    PropagatorOptions opts;
    opts.usePlanetaryPerturbations = true;
    opts.useRelativisticCorrections = true;
    opts.stepSize = 0.1;
    opts.integrator = IntegratorType::RK4;
    
    OrbitPropagator propagator(opts);
    
    std::cout << "Configurazione:\n";
    std::cout << "  Perturbazioni: " << (opts.usePlanetaryPerturbations ? "ON" : "OFF") << "\n";
    std::cout << "  Relativita: " << (opts.useRelativisticCorrections ? "ON" : "OFF") << "\n";
    std::cout << "  Step: " << opts.stepSize << " giorni\n";
    std::cout << "  Integratore: RK4\n\n";
    
    // Elementi orbitali Ceres da MPC (epoca 2461000.5 = 18 Nov 2026)
    OrbitalElements ceresKep;
    ceresKep.epoch = JulianDate(2461000.5);
    ceresKep.a = 2.7656157;
    ceresKep.e = 0.0795763;
    ceresKep.i = 10.58789 * DEG2RAD;
    ceresKep.Omega = 80.24963 * DEG2RAD;
    ceresKep.omega = 73.29974 * DEG2RAD;
    ceresKep.M = 231.53975 * DEG2RAD;
    
    std::cout << "Elementi orbitali Ceres (MPC):\n";
    std::cout << "  Epoca: JD " << std::fixed << std::setprecision(1) << ceresKep.epoch.jd << "\n";
    std::cout << "  a = " << std::setprecision(7) << ceresKep.a << " AU\n";
    std::cout << "  e = " << ceresKep.e << "\n";
    std::cout << "  i = " << (ceresKep.i * RAD2DEG) << " deg\n";
    std::cout << "  Omega = " << (ceresKep.Omega * RAD2DEG) << " deg\n";
    std::cout << "  omega = " << (ceresKep.omega * RAD2DEG) << " deg\n";
    std::cout << "  M = " << (ceresKep.M * RAD2DEG) << " deg\n\n";
    
    // Converti in equinoziali
    EquinoctialElements ceresEq = ceresKep.toEquinoctial();
    
    std::cout << "Elementi equinoziali:\n";
    std::cout << "  a = " << std::setprecision(7) << ceresEq.a << " AU\n";
    std::cout << "  h = " << ceresEq.h << "\n";
    std::cout << "  k = " << ceresEq.k << "\n";
    std::cout << "  p = " << ceresEq.p << "\n";
    std::cout << "  q = " << ceresEq.q << "\n";
    std::cout << "  lambda = " << (ceresEq.lambda * RAD2DEG) << " deg\n\n";
    
    // Target: 29 Nov 2025 (JD 2460646.5)
    JulianDate targetDate(2460646.5);
    double dt = targetDate.jd - ceresKep.epoch.jd;
    
    std::cout << "Target: JD " << targetDate.jd << " (29 Nov 2025)\n";
    std::cout << "Delta t = " << dt << " giorni\n\n";
    
    std::cout << "=== PROPAGAZIONE IN CORSO ===\n";
    std::cout << "(questo potrebbe richiedere alcuni secondi...)\n\n";
    
    try {
        // Prima verifichiamo la posizione all'EPOCA
        OrbitState stateAtEpoch = propagator.propagate(ceresEq, ceresKep.epoch);
        
        // Converti in eclittico
        constexpr double obliquity = 23.4392911 * M_PI / 180.0;
        double cos_eps = cos(obliquity);
        double sin_eps = sin(obliquity);
        
        double x0_ecl = stateAtEpoch.position.x;
        double y0_ecl = stateAtEpoch.position.y * cos_eps + stateAtEpoch.position.z * sin_eps;
        double z0_ecl = -stateAtEpoch.position.y * sin_eps + stateAtEpoch.position.z * cos_eps;
        
        std::cout << "=== VERIFICA POSIZIONE ALL'EPOCA (18 Nov 2026) ===\n";
        std::cout << "Posizione ECLITTICA all'epoca degli elementi:\n";
        std::cout << "  X = " << std::setprecision(6) << x0_ecl << " AU\n";
        std::cout << "  Y = " << y0_ecl << " AU\n";
        std::cout << "  Z = " << z0_ecl << " AU\n";
        std::cout << "(NOTA: Questa dovrebbe corrispondere a JPL per 18 Nov 2026)\n\n";
        
        // Ora propaghiamo al target
        OrbitState result = propagator.propagate(ceresEq, targetDate);
        
        // Conversione EQUATORIALE → ECLITTICO (rotazione inversa)
        double x_ecl = result.position.x;
        double y_ecl = result.position.y * cos_eps + result.position.z * sin_eps;
        double z_ecl = -result.position.y * sin_eps + result.position.z * cos_eps;
        
        std::cout << "=== RISULTATO ===\n";
        std::cout << "Posizione EQUATORIALE (output propagatore):\n";
        std::cout << "  X_eq = " << std::setprecision(6) << result.position.x << " AU\n";
        std::cout << "  Y_eq = " << result.position.y << " AU\n";
        std::cout << "  Z_eq = " << result.position.z << " AU\n";
        
        std::cout << "\nPosizione ECLITTICA (convertita):\n";
        std::cout << "  X_ecl = " << x_ecl << " AU\n";
        std::cout << "  Y_ecl = " << y_ecl << " AU\n";
        std::cout << "  Z_ecl = " << z_ecl << " AU\n";
        
        // Riferimento JPL Horizons per Ceres al 29 Nov 2025 (ECLITTICO)
        double x_jpl = 1.8630;
        double y_jpl = -2.1630;
        double z_jpl = -0.4890;
        
        std::cout << "\nRiferimento JPL Horizons (ECLITTICO):\n";
        std::cout << "  X = " << x_jpl << " AU\n";
        std::cout << "  Y = " << y_jpl << " AU\n";
        std::cout << "  Z = " << z_jpl << " AU\n";
        
        double errX = x_ecl - x_jpl;
        double errY = y_ecl - y_jpl;
        double errZ = z_ecl - z_jpl;
        double errTot = sqrt(errX*errX + errY*errY + errZ*errZ);
        
        std::cout << "\nErrore:\n";
        std::cout << "  DeltaX = " << std::setprecision(4) << errX << " AU = " << std::setprecision(0) << (errX * AU_KM) << " km\n";
        std::cout << "  DeltaY = " << errY << " AU = " << (errY * AU_KM) << " km\n";
        std::cout << "  DeltaZ = " << errZ << " AU = " << (errZ * AU_KM) << " km\n";
        std::cout << "  Totale = " << std::setprecision(4) << errTot << " AU = " << std::setprecision(0) << (errTot * AU_KM) << " km\n";
        
        double dist = sqrt(x_jpl*x_jpl + y_jpl*y_jpl + z_jpl*z_jpl);
        double errAngle = atan(errTot / dist) * RAD2DEG * 3600;
        std::cout << "  Errore angolare = " << std::setprecision(1) << errAngle << " arcsec\n";
        
        auto stats = propagator.getLastStats();
        std::cout << "\nStatistiche:\n";
        std::cout << "  Steps: " << stats.nSteps << "\n";
        std::cout << "  Tempo: " << std::setprecision(2) << stats.computeTime << " s\n";
        
        if (errTot * AU_KM < 10000) {
            std::cout << "\n[OK] Errore < 10000 km\n";
        } else if (errTot * AU_KM < 100000) {
            std::cout << "\n[WARN] Errore 10000-100000 km\n";
        } else {
            std::cout << "\n[FAIL] Errore > 100000 km\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Errore: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
