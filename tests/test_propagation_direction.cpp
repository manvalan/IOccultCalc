/**
 * Test diagnostico per capire l'errore nella propagazione Kepleriana
 * 
 * Strategia:
 * 1. Prendiamo elementi con epoca T0
 * 2. Propaghiamo AVANTI di N giorni → posizione P1
 * 3. Prendiamo la posizione P1 e ricostruiamo elementi
 * 4. Propaghiamo INDIETRO di N giorni → posizione P0'
 * 5. Confrontiamo P0 originale con P0' → errore round-trip
 * 
 * Questo ci dice quanto errore accumula la propagazione Kepleriana.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

// Costanti
const double AU = 149597870.7;  // km
const double MU_SUN = 1.32712440018e11;  // km^3/s^2
const double DEG2RAD = M_PI / 180.0;
const double RAD2DEG = 180.0 / M_PI;
const double SECONDS_PER_DAY = 86400.0;

// Struttura per elementi orbitali
struct OrbitalElements {
    double a;      // semiasse maggiore (AU)
    double e;      // eccentricità
    double i;      // inclinazione (rad)
    double Omega;  // longitudine nodo ascendente (rad)
    double omega;  // argomento del perielio (rad)
    double M;      // anomalia media (rad)
    double epoch;  // JD epoca
};

// Risolve equazione di Keplero: M = E - e*sin(E)
double solveKepler(double M, double e, double tol = 1e-14) {
    // Normalizza M in [0, 2π]
    M = fmod(M, 2.0 * M_PI);
    if (M < 0) M += 2.0 * M_PI;
    
    // Stima iniziale
    double E = (e < 0.8) ? M : M_PI;
    
    // Newton-Raphson
    for (int i = 0; i < 50; ++i) {
        double f = E - e * sin(E) - M;
        double fp = 1.0 - e * cos(E);
        double dE = f / fp;
        E -= dE;
        if (fabs(dE) < tol) break;
    }
    return E;
}

// Da anomalia eccentrica a vera
double eccentricToTrue(double E, double e) {
    double cosE = cos(E);
    double sinE = sin(E);
    double cosNu = (cosE - e) / (1.0 - e * cosE);
    double sinNu = sqrt(1.0 - e * e) * sinE / (1.0 - e * cosE);
    return atan2(sinNu, cosNu);
}

// Da anomalia vera a eccentrica
double trueToEccentric(double nu, double e) {
    double cosNu = cos(nu);
    double sinNu = sin(nu);
    double cosE = (e + cosNu) / (1.0 + e * cosNu);
    double sinE = sqrt(1.0 - e * e) * sinNu / (1.0 + e * cosNu);
    return atan2(sinE, cosE);
}

// Calcola posizione e velocità da elementi orbitali (in AU, AU/day)
void elementsToStateVector(const OrbitalElements& elem, double jd,
                           double& x, double& y, double& z,
                           double& vx, double& vy, double& vz) {
    // Calcola anomalia media al tempo jd
    double a_km = elem.a * AU;
    double n = sqrt(MU_SUN / (a_km * a_km * a_km));  // rad/s (moto medio)
    double dt_seconds = (jd - elem.epoch) * SECONDS_PER_DAY;
    double M = elem.M + n * dt_seconds;
    
    // Risolvi Keplero
    double E = solveKepler(M, elem.e);
    double nu = eccentricToTrue(E, elem.e);
    
    // Distanza
    double r = elem.a * (1.0 - elem.e * cos(E));
    
    // Posizione nel piano orbitale
    double x_orb = r * cos(nu);
    double y_orb = r * sin(nu);
    
    // Velocità nel piano orbitale
    double p = elem.a * (1.0 - elem.e * elem.e);
    double h = sqrt(MU_SUN * p * AU);  // momento angolare specifico
    double vx_orb = -h / (r * AU) * sin(nu) / SECONDS_PER_DAY * AU;  // AU/day
    double vy_orb = h / (r * AU) * (elem.e + cos(nu)) / SECONDS_PER_DAY * AU;  // AU/day
    
    // Rotazione dal piano orbitale al riferimento eclittico
    double cosO = cos(elem.Omega);
    double sinO = sin(elem.Omega);
    double coso = cos(elem.omega);
    double sino = sin(elem.omega);
    double cosi = cos(elem.i);
    double sini = sin(elem.i);
    
    // Matrice di rotazione
    double R11 = cosO * coso - sinO * sino * cosi;
    double R12 = -cosO * sino - sinO * coso * cosi;
    double R21 = sinO * coso + cosO * sino * cosi;
    double R22 = -sinO * sino + cosO * coso * cosi;
    double R31 = sino * sini;
    double R32 = coso * sini;
    
    // Posizione eclittica
    x = R11 * x_orb + R12 * y_orb;
    y = R21 * x_orb + R22 * y_orb;
    z = R31 * x_orb + R32 * y_orb;
    
    // Velocità eclittica
    vx = R11 * vx_orb + R12 * vy_orb;
    vy = R21 * vx_orb + R22 * vy_orb;
    vz = R31 * vx_orb + R32 * vy_orb;
}

// Calcola solo posizione (wrapper semplice)
void elementsToPosition(const OrbitalElements& elem, double jd,
                        double& x, double& y, double& z) {
    double vx, vy, vz;
    elementsToStateVector(elem, jd, x, y, z, vx, vy, vz);
}

// Propaga anomalia media
double propagateMeanAnomaly(double M0, double a, double dt_days) {
    double a_km = a * AU;
    double n = sqrt(MU_SUN / (a_km * a_km * a_km));  // rad/s
    double dt_seconds = dt_days * SECONDS_PER_DAY;
    return M0 + n * dt_seconds;
}

// Test propagazione Kepleriana
void testPropagation(const OrbitalElements& elem, double dt_days, const std::string& name) {
    std::cout << "\n=== Test propagazione " << name << " ===\n";
    std::cout << "Delta t = " << dt_days << " giorni\n";
    std::cout << "Epoca elementi: JD " << std::fixed << std::setprecision(1) << elem.epoch << "\n";
    
    // Posizione all'epoca degli elementi
    double x0, y0, z0;
    elementsToPosition(elem, elem.epoch, x0, y0, z0);
    std::cout << "\nPosizione all'epoca (T0):\n";
    std::cout << "  X = " << std::setprecision(6) << x0 << " AU\n";
    std::cout << "  Y = " << y0 << " AU\n";
    std::cout << "  Z = " << z0 << " AU\n";
    
    // Calcola anomalia vera all'epoca
    double E0 = solveKepler(elem.M, elem.e);
    double nu0 = eccentricToTrue(E0, elem.e);
    std::cout << "  Anomalia media M0 = " << (elem.M * RAD2DEG) << "°\n";
    std::cout << "  Anomalia eccentrica E0 = " << (E0 * RAD2DEG) << "°\n";
    std::cout << "  Anomalia vera nu0 = " << (nu0 * RAD2DEG) << "°\n";
    
    // Posizione al tempo target
    double jd_target = elem.epoch + dt_days;
    double x1, y1, z1;
    elementsToPosition(elem, jd_target, x1, y1, z1);
    std::cout << "\nPosizione al target (T0 + " << dt_days << "d):\n";
    std::cout << "  X = " << x1 << " AU\n";
    std::cout << "  Y = " << y1 << " AU\n";
    std::cout << "  Z = " << z1 << " AU\n";
    
    // Calcola anomalia al target
    double M1 = propagateMeanAnomaly(elem.M, elem.a, dt_days);
    double E1 = solveKepler(M1, elem.e);
    double nu1 = eccentricToTrue(E1, elem.e);
    std::cout << "  Anomalia media M1 = " << (fmod(M1, 2*M_PI) * RAD2DEG) << "°\n";
    std::cout << "  Anomalia eccentrica E1 = " << (E1 * RAD2DEG) << "°\n";
    std::cout << "  Anomalia vera nu1 = " << (nu1 * RAD2DEG) << "°\n";
    
    // TEST ROUND-TRIP: propaga avanti e poi indietro
    // Propaga indietro dalla posizione target
    double x_back, y_back, z_back;
    elementsToPosition(elem, elem.epoch, x_back, y_back, z_back);  // Questo dovrebbe essere identico a x0,y0,z0
    
    // Ora il vero test: crea elementi "nuovi" con epoca al target e propaga indietro
    OrbitalElements elem_at_target = elem;
    elem_at_target.epoch = jd_target;
    elem_at_target.M = fmod(M1, 2*M_PI);
    if (elem_at_target.M < 0) elem_at_target.M += 2*M_PI;
    
    double x0_calc, y0_calc, z0_calc;
    elementsToPosition(elem_at_target, elem.epoch, x0_calc, y0_calc, z0_calc);  // Propaga INDIETRO di dt_days
    
    std::cout << "\n=== TEST ROUND-TRIP ===\n";
    std::cout << "Elementi originali → Propaga avanti " << dt_days << "d → nuovi elementi → Propaga indietro " << dt_days << "d\n";
    std::cout << "\nPosizione originale (T0):\n";
    std::cout << "  X = " << x0 << " AU\n";
    std::cout << "  Y = " << y0 << " AU\n";
    std::cout << "  Z = " << z0 << " AU\n";
    std::cout << "\nPosizione ricalcolata (round-trip):\n";
    std::cout << "  X = " << x0_calc << " AU\n";
    std::cout << "  Y = " << y0_calc << " AU\n";
    std::cout << "  Z = " << z0_calc << " AU\n";
    
    double errX = fabs(x0_calc - x0);
    double errY = fabs(y0_calc - y0);
    double errZ = fabs(z0_calc - z0);
    double errTot = sqrt(errX*errX + errY*errY + errZ*errZ);
    
    std::cout << "\nErrore round-trip:\n";
    std::cout << "  ΔX = " << (errX * AU) << " km (" << errX << " AU)\n";
    std::cout << "  ΔY = " << (errY * AU) << " km (" << errY << " AU)\n";
    std::cout << "  ΔZ = " << (errZ * AU) << " km (" << errZ << " AU)\n";
    std::cout << "  Totale = " << (errTot * AU) << " km (" << errTot << " AU)\n";
    
    if (errTot * AU > 1000) {
        std::cout << "\n⚠️  ERRORE ROUND-TRIP SIGNIFICATIVO!\n";
    } else {
        std::cout << "\n✓ Round-trip OK (errore < 1000 km)\n";
    }
}

// Test con dati MPC
void testWithMPCData() {
    // Ceres - elementi da MPC (epoca 2461000.5 = 18 Nov 2026)
    OrbitalElements ceres;
    ceres.a = 2.7660992;
    ceres.e = 0.0785126;
    ceres.i = 10.58655 * DEG2RAD;
    ceres.Omega = 80.26760 * DEG2RAD;
    ceres.omega = 73.80429 * DEG2RAD;
    ceres.M = 229.39561 * DEG2RAD;  // dalla colonna 'M' MPC
    ceres.epoch = 2461000.5;
    
    // Test con dt = -354 giorni (verso Nov 2025)
    testPropagation(ceres, -354.0, "Ceres (indietro 354 giorni)");
    
    // Test con dt = +354 giorni (in avanti)
    testPropagation(ceres, +354.0, "Ceres (avanti 354 giorni)");
    
    // Test con intervalli più brevi
    testPropagation(ceres, -30.0, "Ceres (indietro 30 giorni)");
    testPropagation(ceres, -100.0, "Ceres (indietro 100 giorni)");
}

// Test confronto con posizione JPL nota
void testVsJPL() {
    std::cout << "\n\n========================================\n";
    std::cout << "=== CONFRONTO CON JPL HORIZONS ===\n";
    std::cout << "========================================\n";
    
    // Ceres al 29 Nov 2025 (JD 2460646.5)
    // Posizione eclittica eliocentrica da JPL Horizons:
    // X ≈ 1.863 AU, Y ≈ -2.163 AU, Z ≈ -0.489 AU
    
    double jd_target = 2460646.5;  // 29 Nov 2025
    
    // Elementi MPC epoca 2461000.5
    OrbitalElements ceres_mpc;
    ceres_mpc.a = 2.7660992;
    ceres_mpc.e = 0.0785126;
    ceres_mpc.i = 10.58655 * DEG2RAD;
    ceres_mpc.Omega = 80.26760 * DEG2RAD;
    ceres_mpc.omega = 73.80429 * DEG2RAD;
    ceres_mpc.M = 229.39561 * DEG2RAD;
    ceres_mpc.epoch = 2461000.5;
    
    double dt = jd_target - ceres_mpc.epoch;
    std::cout << "Data target: JD " << jd_target << " (29 Nov 2025)\n";
    std::cout << "Epoca elementi MPC: JD " << ceres_mpc.epoch << " (18 Nov 2026)\n";
    std::cout << "Delta t = " << dt << " giorni\n";
    
    double x, y, z;
    elementsToPosition(ceres_mpc, jd_target, x, y, z);
    
    std::cout << "\nPosizione calcolata (eclittica eliocentrica):\n";
    std::cout << "  X = " << std::setprecision(4) << x << " AU\n";
    std::cout << "  Y = " << y << " AU\n";
    std::cout << "  Z = " << z << " AU\n";
    
    // Riferimento JPL
    double x_jpl = 1.863;
    double y_jpl = -2.163;
    double z_jpl = -0.489;
    
    std::cout << "\nPosizione JPL Horizons (riferimento):\n";
    std::cout << "  X = " << x_jpl << " AU\n";
    std::cout << "  Y = " << y_jpl << " AU\n";
    std::cout << "  Z = " << z_jpl << " AU\n";
    
    double errX = x - x_jpl;
    double errY = y - y_jpl;
    double errZ = z - z_jpl;
    double errTot = sqrt(errX*errX + errY*errY + errZ*errZ);
    
    std::cout << "\nErrore vs JPL:\n";
    std::cout << "  ΔX = " << std::setprecision(4) << errX << " AU = " << (errX * AU) << " km\n";
    std::cout << "  ΔY = " << errY << " AU = " << (errY * AU) << " km\n";
    std::cout << "  ΔZ = " << errZ << " AU = " << (errZ * AU) << " km\n";
    std::cout << "  Totale = " << errTot << " AU = " << (errTot * AU) << " km\n";
    
    // Errore angolare approssimativo
    double dist_jpl = sqrt(x_jpl*x_jpl + y_jpl*y_jpl + z_jpl*z_jpl);
    double err_angle = atan(errTot / dist_jpl) * RAD2DEG * 3600;  // arcsec
    std::cout << "  Errore angolare ≈ " << err_angle << " arcsec\n";
    
    if (errTot * AU > 1e6) {
        std::cout << "\n⚠️  ERRORE > 1 milione km - Problema nella propagazione!\n";
    }
}

// Analisi dell'errore per diversi dt
void analyzeErrorVsDt() {
    std::cout << "\n\n========================================\n";
    std::cout << "=== ANALISI ERRORE VS TEMPO ===\n";
    std::cout << "========================================\n";
    
    OrbitalElements ceres;
    ceres.a = 2.7660992;
    ceres.e = 0.0785126;
    ceres.i = 10.58655 * DEG2RAD;
    ceres.Omega = 80.26760 * DEG2RAD;
    ceres.omega = 73.80429 * DEG2RAD;
    ceres.M = 229.39561 * DEG2RAD;
    ceres.epoch = 2461000.5;
    
    std::cout << "\n| Δt (giorni) | Direzione | X calc | X JPL | Errore (AU) | Errore (km) |\n";
    std::cout << "|-------------|-----------|--------|-------|-------------|-------------|\n";
    
    // Stima JPL interpolata (approssimativa)
    // Usiamo X_JPL ≈ 1.86 AU a -354 giorni come riferimento
    double x_jpl_ref = 1.863;
    
    for (int dt : {-354, -300, -200, -100, -50, -10, 0, 10, 50, 100}) {
        double jd = ceres.epoch + dt;
        double x, y, z;
        elementsToPosition(ceres, jd, x, y, z);
        
        // Per dt=0, la posizione dovrebbe essere quella degli elementi originali
        // Per dt=-354, confrontiamo con JPL
        double x_expected = (dt == -354) ? 1.863 : x;  // Solo per dt=-354 abbiamo il riferimento JPL
        double err = (dt == -354) ? fabs(x - x_jpl_ref) : 0.0;
        
        std::cout << "| " << std::setw(11) << dt 
                  << " | " << std::setw(9) << (dt < 0 ? "indietro" : (dt > 0 ? "avanti" : "epoca"))
                  << " | " << std::setw(6) << std::setprecision(3) << x 
                  << " | " << std::setw(5) << (dt == -354 ? "1.863" : "N/A")
                  << " | " << std::setw(11) << std::setprecision(4) << err
                  << " | " << std::setw(11) << std::setprecision(0) << (err * AU) << " |\n";
    }
}

// Test del moto medio
void testMeanMotion() {
    std::cout << "\n\n========================================\n";
    std::cout << "=== VERIFICA MOTO MEDIO ===\n";
    std::cout << "========================================\n";
    
    // Ceres: semiasse maggiore a = 2.766 AU
    double a = 2.7660992;
    double a_km = a * AU;
    
    // Moto medio n = sqrt(GM/a³)
    double n_rad_s = sqrt(MU_SUN / (a_km * a_km * a_km));  // rad/s
    double n_deg_day = n_rad_s * RAD2DEG * SECONDS_PER_DAY;  // °/giorno
    
    std::cout << "Semiasse maggiore a = " << a << " AU = " << a_km << " km\n";
    std::cout << "Moto medio n = " << std::scientific << n_rad_s << " rad/s\n";
    std::cout << "            = " << std::fixed << std::setprecision(4) << n_deg_day << " °/giorno\n";
    
    // Periodo orbitale
    double T_seconds = 2 * M_PI / n_rad_s;
    double T_days = T_seconds / SECONDS_PER_DAY;
    double T_years = T_days / 365.25;
    std::cout << "\nPeriodo orbitale T = " << T_days << " giorni = " << std::setprecision(2) << T_years << " anni\n";
    
    // Propagazione di 354 giorni = quanti gradi?
    double dt = 354;
    double dM = n_deg_day * dt;
    std::cout << "\nIn " << dt << " giorni, l'anomalia media cambia di " << dM << "°\n";
    std::cout << "Numero di orbite: " << (dt / T_days) << "\n";
    
    // Verifica con MPC: M all'epoca = 229.4°
    // Dopo 354 giorni indietro: M = 229.4 - 75.5 = 153.9° circa
    double M0 = 229.39561;
    double M_after = M0 - dM;
    std::cout << "\nAnomalia media:\n";
    std::cout << "  All'epoca (Nov 2026): M = " << M0 << "°\n";
    std::cout << "  A Nov 2025 (calc): M = " << M0 << " - " << dM << " = " << M_after << "°\n";
    std::cout << "  Normalizzata: " << fmod(M_after + 720, 360) << "°\n";
}

int main() {
    std::cout << "===============================================\n";
    std::cout << " TEST DIAGNOSTICO PROPAGAZIONE KEPLERIANA\n";
    std::cout << "===============================================\n";
    
    testMeanMotion();
    testVsJPL();
    analyzeErrorVsDt();
    
    // Test round-trip dettagliato
    testWithMPCData();
    
    return 0;
}
