/**
 * Test propagazione IN AVANTI
 * 
 * Strategia: partiamo da elementi con epoca PASSATA e propaghiamo verso OGGI
 * così possiamo confrontare con posizioni note.
 * 
 * Usiamo elementi AstDyS .eq0 che hanno epoca 1994 e propaghiamo verso 2025.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/types.h"
#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/astdys_client.h"

using namespace ioccultcalc;

const double RAD2DEG = 180.0 / M_PI;
const double DEG2RAD = M_PI / 180.0;
const double AU_KM = 149597870.7;

// Conversione equatoriale -> eclittico
void equatorialToEcliptic(double x_eq, double y_eq, double z_eq,
                          double& x_ecl, double& y_ecl, double& z_ecl) {
    const double obliquity = 23.4392911 * (M_PI / 180.0);
    double cos_eps = cos(obliquity);
    double sin_eps = sin(obliquity);
    
    x_ecl = x_eq;
    y_ecl = y_eq * cos_eps + z_eq * sin_eps;
    z_ecl = -y_eq * sin_eps + z_eq * cos_eps;
}

int main() {
    std::cout << "=============================================\n";
    std::cout << " TEST PROPAGAZIONE IN AVANTI\n";
    std::cout << "=============================================\n\n";
    
    // Elementi Ceres MPC con epoca Nov 2026
    // Li propaghiamo verso OGGI (29 Nov 2025) = -354 giorni
    // POI li propaghiamo verso il FUTURO per vedere se funziona
    
    OrbitalElements ceres;
    ceres.epoch = JulianDate(2461000.5);  // 18 Nov 2026
    ceres.a = 2.7656157;
    ceres.e = 0.0795763;
    ceres.i = 10.58789 * DEG2RAD;
    ceres.Omega = 80.24963 * DEG2RAD;
    ceres.omega = 73.29974 * DEG2RAD;
    ceres.M = 231.53975 * DEG2RAD;
    ceres.name = "Ceres";
    
    std::cout << "Elementi Ceres (MPC, epoca Nov 2026):\n";
    std::cout << "  Epoca: JD " << std::fixed << std::setprecision(1) << ceres.epoch.jd << "\n";
    std::cout << "  a = " << std::setprecision(7) << ceres.a << " AU\n";
    std::cout << "  e = " << ceres.e << "\n";
    std::cout << "  i = " << (ceres.i * RAD2DEG) << " deg\n\n";
    
    EquinoctialElements ceresEq = ceres.toEquinoctial();
    
    // Configura propagatore
    PropagatorOptions opts;
    opts.usePlanetaryPerturbations = true;
    opts.useRelativisticCorrections = true;
    opts.stepSize = 0.1;  // 0.1 giorni
    opts.integrator = IntegratorType::RK4;
    
    OrbitPropagator propagator(opts);
    
    std::cout << "Configurazione propagatore:\n";
    std::cout << "  Perturbazioni planetarie: ON\n";
    std::cout << "  Relativita: ON\n";
    std::cout << "  Step: 0.1 giorni\n";
    std::cout << "  Integratore: RK4\n\n";
    
    // === TEST 1: Posizione all'EPOCA ===
    std::cout << "=== TEST 1: Posizione all'EPOCA (no propagazione) ===\n";
    OrbitState stateEpoch = propagator.propagate(ceresEq, ceres.epoch);
    
    double x0, y0, z0;
    equatorialToEcliptic(stateEpoch.position.x, stateEpoch.position.y, stateEpoch.position.z,
                         x0, y0, z0);
    
    std::cout << "Posizione eclittica al JD " << ceres.epoch.jd << " (18 Nov 2026):\n";
    std::cout << "  X = " << std::setprecision(6) << x0 << " AU\n";
    std::cout << "  Y = " << y0 << " AU\n";
    std::cout << "  Z = " << z0 << " AU\n\n";
    
    // === TEST 2: Propagazione +30 giorni (in avanti verso Dic 2026) ===
    std::cout << "=== TEST 2: Propagazione +30 giorni ===\n";
    JulianDate target30(ceres.epoch.jd + 30);
    OrbitState state30 = propagator.propagate(ceresEq, target30);
    
    double x30, y30, z30;
    equatorialToEcliptic(state30.position.x, state30.position.y, state30.position.z,
                         x30, y30, z30);
    
    std::cout << "Posizione eclittica al JD " << target30.jd << " (18 Dic 2026):\n";
    std::cout << "  X = " << x30 << " AU\n";
    std::cout << "  Y = " << y30 << " AU\n";
    std::cout << "  Z = " << z30 << " AU\n";
    
    // Verifica spostamento ragionevole
    double dist30 = sqrt(pow(x30-x0,2) + pow(y30-y0,2) + pow(z30-z0,2));
    std::cout << "  Spostamento: " << dist30 << " AU = " << (dist30 * AU_KM) << " km\n\n";
    
    // === TEST 3: Propagazione +100 giorni ===
    std::cout << "=== TEST 3: Propagazione +100 giorni ===\n";
    JulianDate target100(ceres.epoch.jd + 100);
    OrbitState state100 = propagator.propagate(ceresEq, target100);
    
    double x100, y100, z100;
    equatorialToEcliptic(state100.position.x, state100.position.y, state100.position.z,
                         x100, y100, z100);
    
    std::cout << "Posizione eclittica al JD " << target100.jd << " (27 Feb 2027):\n";
    std::cout << "  X = " << x100 << " AU\n";
    std::cout << "  Y = " << y100 << " AU\n";
    std::cout << "  Z = " << z100 << " AU\n";
    
    double dist100 = sqrt(pow(x100-x0,2) + pow(y100-y0,2) + pow(z100-z0,2));
    std::cout << "  Spostamento: " << dist100 << " AU = " << (dist100 * AU_KM) << " km\n\n";
    
    // === TEST 4: Round-trip +30 e poi -30 ===
    std::cout << "=== TEST 4: Round-trip (+30 giorni, poi -30 giorni) ===\n";
    
    // Propaga +30
    OrbitState stateForward = propagator.propagate(ceresEq, JulianDate(ceres.epoch.jd + 30));
    
    // Ora crea nuovi elementi da questo stato e propaga -30
    EquinoctialElements eqForward = propagator.stateToElements(stateForward);
    OrbitState stateBack = propagator.propagate(eqForward, ceres.epoch);
    
    double xb, yb, zb;
    equatorialToEcliptic(stateBack.position.x, stateBack.position.y, stateBack.position.z,
                         xb, yb, zb);
    
    std::cout << "Posizione originale (epoca):\n";
    std::cout << "  X = " << x0 << " AU\n";
    std::cout << "  Y = " << y0 << " AU\n";
    std::cout << "  Z = " << z0 << " AU\n";
    
    std::cout << "Posizione dopo round-trip:\n";
    std::cout << "  X = " << xb << " AU\n";
    std::cout << "  Y = " << yb << " AU\n";
    std::cout << "  Z = " << zb << " AU\n";
    
    double errRT = sqrt(pow(xb-x0,2) + pow(yb-y0,2) + pow(zb-z0,2));
    std::cout << "Errore round-trip: " << errRT << " AU = " << (errRT * AU_KM) << " km\n";
    
    if (errRT * AU_KM < 1000) {
        std::cout << "[OK] Round-trip corretto (errore < 1000 km)\n\n";
    } else {
        std::cout << "[WARN] Errore round-trip significativo!\n\n";
    }
    
    // === TEST 5: Propagazione -354 giorni (verso Nov 2025) ===
    std::cout << "=== TEST 5: Propagazione -354 giorni (verso Nov 2025) ===\n";
    JulianDate targetNov2025(2460646.5);  // 29 Nov 2025
    OrbitState stateNov2025 = propagator.propagate(ceresEq, targetNov2025);
    
    double xn, yn, zn;
    equatorialToEcliptic(stateNov2025.position.x, stateNov2025.position.y, stateNov2025.position.z,
                         xn, yn, zn);
    
    std::cout << "Posizione eclittica al JD 2460646.5 (29 Nov 2025):\n";
    std::cout << "  X = " << xn << " AU\n";
    std::cout << "  Y = " << yn << " AU\n";
    std::cout << "  Z = " << zn << " AU\n";
    
    // Riferimento JPL (approssimativo)
    double x_jpl = 1.8630, y_jpl = -2.1630, z_jpl = -0.4890;
    std::cout << "\nRiferimento JPL Horizons:\n";
    std::cout << "  X = " << x_jpl << " AU\n";
    std::cout << "  Y = " << y_jpl << " AU\n";
    std::cout << "  Z = " << z_jpl << " AU\n";
    
    double errJPL = sqrt(pow(xn-x_jpl,2) + pow(yn-y_jpl,2) + pow(zn-z_jpl,2));
    std::cout << "\nErrore vs JPL: " << errJPL << " AU = " << (errJPL * AU_KM) << " km\n";
    
    auto stats = propagator.getLastStats();
    std::cout << "\nStatistiche propagazione:\n";
    std::cout << "  Steps: " << stats.nSteps << "\n";
    std::cout << "  Tempo: " << std::setprecision(2) << stats.computeTime << " s\n";
    
    // Risultato finale
    std::cout << "\n=============================================\n";
    std::cout << " RIEPILOGO\n";
    std::cout << "=============================================\n";
    std::cout << "Round-trip error: " << (errRT * AU_KM) << " km ";
    std::cout << (errRT * AU_KM < 1000 ? "[OK]" : "[FAIL]") << "\n";
    std::cout << "Errore vs JPL (-354d): " << (errJPL * AU_KM) << " km ";
    std::cout << (errJPL * AU_KM < 100000 ? "[OK]" : "[FAIL]") << "\n";
    
    return 0;
}
