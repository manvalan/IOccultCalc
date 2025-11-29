/**
 * Debug: confronto conversione Kepleriano vs Equinoziale
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
const double MU_SUN = 1.32712440018e11;
const double SECONDS_PER_DAY = 86400.0;

// Calcolo Kepleriano diretto (come nel test semplice)
double solveKepler(double M, double e, double tol = 1e-14) {
    M = fmod(M, 2.0 * M_PI);
    if (M < 0) M += 2.0 * M_PI;
    double E = (e < 0.8) ? M : M_PI;
    for (int i = 0; i < 50; ++i) {
        double f = E - e * sin(E) - M;
        double fp = 1.0 - e * cos(E);
        double dE = f / fp;
        E -= dE;
        if (fabs(dE) < tol) break;
    }
    return E;
}

double eccentricToTrue(double E, double e) {
    double cosE = cos(E);
    double sinE = sin(E);
    double cosNu = (cosE - e) / (1.0 - e * cosE);
    double sinNu = sqrt(1.0 - e * e) * sinE / (1.0 - e * cosE);
    return atan2(sinNu, cosNu);
}

void keplerianToPosition(double a, double e, double i, double Omega, double omega, double M,
                         double& x, double& y, double& z) {
    double E = solveKepler(M, e);
    double nu = eccentricToTrue(E, e);
    double r = a * (1.0 - e * cos(E));
    
    double x_orb = r * cos(nu);
    double y_orb = r * sin(nu);
    
    double cosO = cos(Omega), sinO = sin(Omega);
    double coso = cos(omega), sino = sin(omega);
    double cosi = cos(i), sini = sin(i);
    
    double R11 = cosO * coso - sinO * sino * cosi;
    double R12 = -cosO * sino - sinO * coso * cosi;
    double R21 = sinO * coso + cosO * sino * cosi;
    double R22 = -sinO * sino + cosO * coso * cosi;
    double R31 = sino * sini;
    double R32 = coso * sini;
    
    x = R11 * x_orb + R12 * y_orb;
    y = R21 * x_orb + R22 * y_orb;
    z = R31 * x_orb + R32 * y_orb;
}

int main() {
    std::cout << "===========================================\n";
    std::cout << " DEBUG: CONVERSIONE KEPLERIANO VS EQUINOZIALE\n";
    std::cout << "===========================================\n\n";
    
    // Elementi Kepleriani MPC
    OrbitalElements kep;
    kep.epoch = JulianDate(2461000.5);
    kep.a = 2.7656157;
    kep.e = 0.0795763;
    kep.i = 10.58789 * DEG2RAD;
    kep.Omega = 80.24963 * DEG2RAD;
    kep.omega = 73.29974 * DEG2RAD;
    kep.M = 231.53975 * DEG2RAD;
    
    std::cout << "Elementi Kepleriani (MPC):\n";
    std::cout << "  a = " << kep.a << " AU\n";
    std::cout << "  e = " << kep.e << "\n";
    std::cout << "  i = " << (kep.i * RAD2DEG) << " deg\n";
    std::cout << "  Omega = " << (kep.Omega * RAD2DEG) << " deg\n";
    std::cout << "  omega = " << (kep.omega * RAD2DEG) << " deg\n";
    std::cout << "  M = " << (kep.M * RAD2DEG) << " deg\n\n";
    
    // Calcolo posizione DIRETTA da Kepleriani (eclittico)
    double x_kep, y_kep, z_kep;
    keplerianToPosition(kep.a, kep.e, kep.i, kep.Omega, kep.omega, kep.M,
                        x_kep, y_kep, z_kep);
    
    std::cout << "=== METODO 1: Kepleriano diretto ===\n";
    std::cout << "Posizione eclittica all'epoca:\n";
    std::cout << "  X = " << std::setprecision(6) << x_kep << " AU\n";
    std::cout << "  Y = " << y_kep << " AU\n";
    std::cout << "  Z = " << z_kep << " AU\n\n";
    
    // Converti in equinoziali
    EquinoctialElements eq = kep.toEquinoctial();
    
    std::cout << "Elementi equinoziali (convertiti):\n";
    std::cout << "  a = " << eq.a << " AU\n";
    std::cout << "  h = " << eq.h << "\n";
    std::cout << "  k = " << eq.k << "\n";
    std::cout << "  p = " << eq.p << "\n";
    std::cout << "  q = " << eq.q << "\n";
    std::cout << "  lambda = " << (eq.lambda * RAD2DEG) << " deg\n\n";
    
    // Riconverti in Kepleriani per verifica
    OrbitalElements kep_back = eq.toKeplerian();
    
    std::cout << "Elementi Kepleriani (riconvertiti):\n";
    std::cout << "  a = " << kep_back.a << " AU\n";
    std::cout << "  e = " << kep_back.e << "\n";
    std::cout << "  i = " << (kep_back.i * RAD2DEG) << " deg\n";
    std::cout << "  Omega = " << (kep_back.Omega * RAD2DEG) << " deg\n";
    std::cout << "  omega = " << (kep_back.omega * RAD2DEG) << " deg\n";
    std::cout << "  M = " << (kep_back.M * RAD2DEG) << " deg\n\n";
    
    // Usa OrbitPropagator per ottenere posizione
    PropagatorOptions opts;
    opts.usePlanetaryPerturbations = false;  // Solo 2-body per questo test
    OrbitPropagator propagator(opts);
    
    OrbitState state = propagator.propagate(eq, kep.epoch);
    
    // Converti da equatoriale a eclittico
    constexpr double obliquity = 23.4392911 * M_PI / 180.0;
    double cos_eps = cos(obliquity);
    double sin_eps = sin(obliquity);
    
    double x_eq = state.position.x;
    double y_eq_ecl = state.position.y * cos_eps + state.position.z * sin_eps;
    double z_eq_ecl = -state.position.y * sin_eps + state.position.z * cos_eps;
    
    std::cout << "=== METODO 2: OrbitPropagator (elementsToState) ===\n";
    std::cout << "Posizione equatoriale:\n";
    std::cout << "  X = " << state.position.x << " AU\n";
    std::cout << "  Y = " << state.position.y << " AU\n";
    std::cout << "  Z = " << state.position.z << " AU\n";
    std::cout << "Posizione eclittica (convertita):\n";
    std::cout << "  X = " << x_eq << " AU\n";
    std::cout << "  Y = " << y_eq_ecl << " AU\n";
    std::cout << "  Z = " << z_eq_ecl << " AU\n\n";
    
    // Confronto
    std::cout << "=== CONFRONTO ===\n";
    double dx = x_eq - x_kep;
    double dy = y_eq_ecl - y_kep;
    double dz = z_eq_ecl - z_kep;
    double d = sqrt(dx*dx + dy*dy + dz*dz);
    
    std::cout << "Differenza tra i due metodi:\n";
    std::cout << "  dX = " << dx << " AU = " << (dx * AU_KM) << " km\n";
    std::cout << "  dY = " << dy << " AU = " << (dy * AU_KM) << " km\n";
    std::cout << "  dZ = " << dz << " AU = " << (dz * AU_KM) << " km\n";
    std::cout << "  Tot = " << d << " AU = " << (d * AU_KM) << " km\n";
    
    if (d * AU_KM < 100) {
        std::cout << "\n[OK] I due metodi coincidono!\n";
    } else {
        std::cout << "\n[ERRORE] I due metodi NON coincidono!\n";
    }
    
    return 0;
}
