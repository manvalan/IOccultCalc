/**
 * @file test_orbfit_force_model.cpp
 * @brief Test del modello di forze OrbFit-compatible
 * 
 * Verifica che il modello di forze sia correttamente implementato
 * confrontando con valori noti.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "ioccultcalc/orbfit_force_model.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/orbital_elements.h"

using namespace ioccultcalc;

const double RAD2DEG = 180.0 / M_PI;
const double DEG2RAD = M_PI / 180.0;
const double AU_KM = 149597870.7;

// Conversione equatoriale -> eclittico
void equatorialToEcliptic(double x_eq, double y_eq, double z_eq,
                          double& x_ecl, double& y_ecl, double& z_ecl) {
    const double obliquity = 23.4392911 * DEG2RAD;
    double cos_eps = cos(obliquity);
    double sin_eps = sin(obliquity);
    
    x_ecl = x_eq;
    y_ecl = y_eq * cos_eps + z_eq * sin_eps;
    z_ecl = -y_eq * sin_eps + z_eq * cos_eps;
}

int main() {
    std::cout << "=============================================\n";
    std::cout << " TEST MODELLO FORZE ORBFIT-COMPATIBLE\n";
    std::cout << "=============================================\n\n";
    
    // === TEST 1: Verifica inizializzazione ===
    std::cout << "=== TEST 1: Inizializzazione ===\n";
    
    OrbFitForceOptions opts = OrbFitForceOptions::highPrecisionConfig();
    OrbFitForceModel forceModel(opts);
    
    if (!forceModel.initializeEphemeris()) {
        std::cerr << "ERRORE: Impossibile inizializzare effemeridi\n";
        return 1;
    }
    std::cout << "[OK] Modello inizializzato\n\n";
    
    // === TEST 2: Accelerazione su Ceres (posizione nota) ===
    std::cout << "=== TEST 2: Accelerazione su Ceres ===\n";
    
    // Posizione Ceres al JD 2460646.5 (29 Nov 2025) - coordinate eclittiche
    // Valori approssimativi da JPL Horizons
    Vector3D posCeres(1.863, -2.163, -0.489);  // AU, eclittico
    Vector3D velCeres(0.0087, 0.0053, 0.00082);  // AU/day, eclittico
    double jd = 2460646.5;
    
    std::cout << "Posizione Ceres (eclittica):\n";
    std::cout << "  X = " << posCeres.x << " AU\n";
    std::cout << "  Y = " << posCeres.y << " AU\n";
    std::cout << "  Z = " << posCeres.z << " AU\n";
    std::cout << "  r = " << posCeres.magnitude() << " AU\n\n";
    
    // Calcola accelerazione
    Vector3D sunAccel, planetAccel, asteroidAccel, relAccel;
    Vector3D totalAccel = forceModel.computeAccelerationDetailed(
        posCeres, velCeres, jd,
        sunAccel, planetAccel, asteroidAccel, relAccel);
    
    std::cout << "Accelerazioni (AU/day²):\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  Solare:        " << sunAccel.magnitude() << "\n";
    std::cout << "  Planetarie:    " << planetAccel.magnitude() << "\n";
    std::cout << "  Asteroidi:     " << asteroidAccel.magnitude() << "\n";
    std::cout << "  Relativistica: " << relAccel.magnitude() << "\n";
    std::cout << "  TOTALE:        " << totalAccel.magnitude() << "\n\n";
    
    // Verifica ordine di grandezza
    // Per Ceres a ~2.9 AU, accelerazione solare ~GM0/r² ~ 3e-5 AU/day²
    double r = posCeres.magnitude();
    double expectedSunAccel = OrbFitConstants::GM0 / (r * r);
    std::cout << "Verifica ordine grandezza:\n";
    std::cout << "  Atteso (Sole): " << expectedSunAccel << " AU/day²\n";
    std::cout << "  Calcolato:     " << sunAccel.magnitude() << " AU/day²\n";
    
    double relError = std::abs(sunAccel.magnitude() - expectedSunAccel) / expectedSunAccel;
    if (relError < 0.01) {
        std::cout << "  [OK] Errore < 1%\n\n";
    } else {
        std::cout << "  [WARN] Errore = " << (relError * 100) << "%\n\n";
    }
    
    // === TEST 3: Confronto perturbazioni planetarie ===
    std::cout << "=== TEST 3: Perturbazioni Planetarie ===\n";
    
    // Rapporto perturbazioni/solare
    double pertRatio = planetAccel.magnitude() / sunAccel.magnitude();
    std::cout << "Rapporto perturbazioni/solare: " << std::fixed << std::setprecision(6) 
              << (pertRatio * 1e6) << " ppm\n";
    
    // Per Ceres, le perturbazioni dovrebbero essere ~10⁻⁵ - 10⁻⁴ rispetto al Sole
    if (pertRatio > 1e-6 && pertRatio < 1e-3) {
        std::cout << "[OK] Perturbazioni nel range atteso\n\n";
    } else {
        std::cout << "[WARN] Perturbazioni fuori range\n\n";
    }
    
    // === TEST 4: Confronto correzione relativistica ===
    std::cout << "=== TEST 4: Correzione Relativistica ===\n";
    
    // Rapporto relatività/solare
    double relRatio = relAccel.magnitude() / sunAccel.magnitude();
    std::cout << "Rapporto relatività/solare: " << std::scientific << relRatio << "\n";
    
    // Per asteroidi della fascia principale, ~10⁻⁸
    if (relRatio > 1e-10 && relRatio < 1e-6) {
        std::cout << "[OK] Correzione relativistica nel range atteso\n\n";
    } else {
        std::cout << "[WARN] Correzione relativistica fuori range\n\n";
    }
    
    // === TEST 5: Propagazione con nuovo modello ===
    std::cout << "=== TEST 5: Propagazione con OrbFit Force Model ===\n";
    
    // Elementi Ceres epoca Nov 2026 (MPC)
    OrbitalElements ceres;
    ceres.epoch = JulianDate(2461000.5);  // 18 Nov 2026
    ceres.a = 2.7656157;
    ceres.e = 0.0795763;
    ceres.i = 10.58789 * DEG2RAD;
    ceres.Omega = 80.24963 * DEG2RAD;
    ceres.omega = 73.29974 * DEG2RAD;
    ceres.M = 231.53975 * DEG2RAD;
    
    // Configura propagatore con nuovo modello
    PropagatorOptions propOpts;
    propOpts.usePlanetaryPerturbations = true;
    propOpts.useRelativisticCorrections = true;
    propOpts.stepSize = 0.1;
    propOpts.integrator = IntegratorType::RK4;
    
    OrbitPropagator propagator(propOpts);
    
    // Propaga verso Nov 2025 (-354 giorni)
    EquinoctialElements ceresEq = ceres.toEquinoctial();
    JulianDate targetDate(2460646.5);  // 29 Nov 2025
    
    std::cout << "Propagazione da " << ceres.epoch.jd << " a " << targetDate.jd << "\n";
    std::cout << "Delta: " << (targetDate.jd - ceres.epoch.jd) << " giorni\n\n";
    
    OrbitState stateNov2025 = propagator.propagate(ceresEq, targetDate);
    
    double xn, yn, zn;
    equatorialToEcliptic(stateNov2025.position.x, stateNov2025.position.y, 
                         stateNov2025.position.z, xn, yn, zn);
    
    std::cout << "Posizione calcolata (eclittica):\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  X = " << xn << " AU\n";
    std::cout << "  Y = " << yn << " AU\n";
    std::cout << "  Z = " << zn << " AU\n";
    
    // Riferimento JPL Horizons (29 Nov 2025)
    double x_jpl = 1.8630, y_jpl = -2.1630, z_jpl = -0.4890;
    std::cout << "\nRiferimento JPL Horizons:\n";
    std::cout << "  X = " << x_jpl << " AU\n";
    std::cout << "  Y = " << y_jpl << " AU\n";
    std::cout << "  Z = " << z_jpl << " AU\n";
    
    double errJPL = sqrt(pow(xn-x_jpl,2) + pow(yn-y_jpl,2) + pow(zn-z_jpl,2));
    std::cout << "\nErrore vs JPL: " << errJPL << " AU = " << (errJPL * AU_KM) << " km\n";
    
    // === RIEPILOGO ===
    std::cout << "\n=============================================\n";
    std::cout << " RIEPILOGO\n";
    std::cout << "=============================================\n";
    std::cout << "Accelerazione solare: " << std::scientific << sunAccel.magnitude() << " AU/day²\n";
    std::cout << "Perturbazioni planetarie: " << planetAccel.magnitude() << " AU/day²\n";
    std::cout << "Correzione relativistica: " << relAccel.magnitude() << " AU/day²\n";
    std::cout << "Errore propagazione: " << std::fixed << (errJPL * AU_KM) << " km\n";
    
    if (errJPL * AU_KM < 50000) {
        std::cout << "\n[OK] Modello OrbFit funziona correttamente!\n";
    } else if (errJPL * AU_KM < 100000) {
        std::cout << "\n[WARN] Errore ridotto ma ancora significativo\n";
    } else {
        std::cout << "\n[FAIL] Errore ancora troppo grande - verificare modello\n";
    }
    
    return 0;
}
