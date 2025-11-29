/**
 * Test: verifichiamo se gli elementi MPC all'EPOCA coincidono con JPL
 * 
 * Se gli elementi MPC descrivono correttamente l'orbita all'epoca,
 * la posizione calcolata all'epoca dovrebbe coincidere con JPL all'epoca.
 */

#include <iostream>
#include <iomanip>
#include <cmath>

const double AU = 149597870.7;
const double MU_SUN = 1.32712440018e11;
const double DEG2RAD = M_PI / 180.0;
const double RAD2DEG = 180.0 / M_PI;
const double SECONDS_PER_DAY = 86400.0;

struct OrbitalElements {
    double a, e, i, Omega, omega, M, epoch;
};

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

void elementsToPosition(const OrbitalElements& elem, double jd,
                        double& x, double& y, double& z) {
    double a_km = elem.a * AU;
    double n = sqrt(MU_SUN / (a_km * a_km * a_km));
    double dt_seconds = (jd - elem.epoch) * SECONDS_PER_DAY;
    double M = elem.M + n * dt_seconds;
    
    double E = solveKepler(M, elem.e);
    double nu = eccentricToTrue(E, elem.e);
    double r = elem.a * (1.0 - elem.e * cos(E));
    
    double x_orb = r * cos(nu);
    double y_orb = r * sin(nu);
    
    double cosO = cos(elem.Omega), sinO = sin(elem.Omega);
    double coso = cos(elem.omega), sino = sin(elem.omega);
    double cosi = cos(elem.i), sini = sin(elem.i);
    
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
    std::cout << " VERIFICA: ELEMENTI MPC COINCIDONO CON JPL?\n";
    std::cout << "===========================================\n\n";
    
    // Elementi MPC per Ceres - epoca 2461000.5 (18 Nov 2026)
    OrbitalElements ceres_mpc;
    ceres_mpc.a = 2.7660992;
    ceres_mpc.e = 0.0785126;
    ceres_mpc.i = 10.58655 * DEG2RAD;
    ceres_mpc.Omega = 80.26760 * DEG2RAD;
    ceres_mpc.omega = 73.80429 * DEG2RAD;
    ceres_mpc.M = 229.39561 * DEG2RAD;
    ceres_mpc.epoch = 2461000.5;
    
    std::cout << "=== TEST 1: Posizione ALL'EPOCA ===\n";
    std::cout << "Se elementi MPC sono corretti, la posizione all'epoca\n";
    std::cout << "dovrebbe coincidere con JPL Horizons per quella data.\n\n";
    
    double x, y, z;
    elementsToPosition(ceres_mpc, ceres_mpc.epoch, x, y, z);
    
    std::cout << "Data: JD " << std::fixed << std::setprecision(1) << ceres_mpc.epoch << " (18 Nov 2026)\n";
    std::cout << "\nPosizione calcolata da elementi MPC:\n";
    std::cout << "  X = " << std::setprecision(6) << x << " AU\n";
    std::cout << "  Y = " << y << " AU\n";
    std::cout << "  Z = " << z << " AU\n";
    
    // DEVO OTTENERE DA JPL HORIZONS LA POSIZIONE AL 18 NOV 2026!
    // Per ora lascio un placeholder
    std::cout << "\nPosizione JPL Horizons al 18 Nov 2026:\n";
    std::cout << "  (Da verificare con query a Horizons)\n";
    
    std::cout << "\n=== TEST 2: Posizione a diverse epoche ===\n";
    std::cout << "Confronto tra posizione calcolata e posizione attesa\n\n";
    
    // Per Nov 2025, abbiamo già il riferimento JPL
    double jd_nov2025 = 2460646.5;
    elementsToPosition(ceres_mpc, jd_nov2025, x, y, z);
    
    std::cout << "Data: JD 2460646.5 (29 Nov 2025)\n";
    std::cout << "Posizione calcolata: X=" << x << ", Y=" << y << ", Z=" << z << " AU\n";
    std::cout << "Posizione JPL:       X=1.8630, Y=-2.1630, Z=-0.4890 AU\n";
    
    double errX = x - 1.8630;
    double errY = y - (-2.1630);
    double errZ = z - (-0.4890);
    double errTot = sqrt(errX*errX + errY*errY + errZ*errZ);
    std::cout << "Errore totale: " << errTot << " AU = " << (errTot * AU) << " km\n";
    
    std::cout << "\n=== ANALISI ===\n";
    std::cout << "Se l'errore è dovuto alle perturbazioni planetarie:\n";
    std::cout << "- Ceres ha perturbazioni significative da Giove\n";
    std::cout << "- In 354 giorni, le perturbazioni possono accumulare ~0.1 AU di errore\n";
    std::cout << "- Gli elementi MPC sono 'medi', non 'osculanti' all'epoca target\n";
    
    std::cout << "\n=== SOLUZIONE PROPOSTA ===\n";
    std::cout << "1. Usare OrbitPropagator con perturbazioni per propagare dall'epoca\n";
    std::cout << "2. Oppure: scaricare elementi MPC con epoca PIU' VICINA alla data target\n";
    std::cout << "3. Oppure: usare effemeridi pre-calcolate (SPK) per asteroidi principali\n";
    
    // Vediamo quanto sarebbe l'errore con epoca più vicina
    std::cout << "\n=== STIMA ERRORE CON EPOCHE DIVERSE ===\n";
    std::cout << "| Δt (giorni) | Errore stimato (km) |\n";
    std::cout << "|-------------|---------------------|\n";
    
    // Approssimazione: errore ~ k * Δt² (per perturbazioni seconde ordine)
    // Con Δt=354, err=15M km → k = 15M/354² ≈ 120 km/day²
    double k = 15000000.0 / (354.0 * 354.0);
    for (int dt : {2, 10, 30, 50, 100, 200, 354}) {
        double err_est = k * dt * dt;
        std::cout << "| " << std::setw(11) << dt 
                  << " | " << std::setw(19) << std::setprecision(0) << err_est << " |\n";
    }
    
    return 0;
}
