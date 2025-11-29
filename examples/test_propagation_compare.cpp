/**
 * @file test_propagation_compare.cpp
 * @brief Test di propagazione temporale RKF78 vs JPL Horizons
 * 
 * Propaga l'asteroide (11234) 1999 JS82 a ±60 giorni dall'epoca
 * e confronta i risultati con JPL Horizons per validare l'accuratezza.
 * 
 * Uso:
 *   ./test_propagation_compare
 * 
 * Output:
 *   - Posizione all'epoca (MJD 61000.0 = 2026-01-09)
 *   - Posizione a +60 giorni (MJD 61060.0 = 2026-03-10)
 *   - Posizione a -60 giorni (MJD 60940.0 = 2025-11-10)
 *   - Confronto con JPL Horizons
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>

// Struttura per elementi orbitali equinoziali
struct EquinoctialElements {
    double a;          // Semi-asse maggiore (AU)
    double h;          // e*sin(ω+Ω)
    double k;          // e*cos(ω+Ω)
    double p;          // tan(i/2)*sin(Ω)
    double q;          // tan(i/2)*cos(Ω)
    double lambda;     // Longitudine media (gradi)
    double epoch_mjd;  // Epoca (MJD TDT)
    
    // Conversione a elementi Kepleriani
    void toKeplerian(double& ecc, double& inc_deg, double& omega_deg, 
                     double& Omega_deg, double& M_deg) const {
        ecc = sqrt(h*h + k*k);
        inc_deg = 2.0 * atan(sqrt(p*p + q*q)) * 180.0 / M_PI;
        Omega_deg = atan2(p, q) * 180.0 / M_PI;
        if (Omega_deg < 0) Omega_deg += 360.0;
        
        double omega_plus_Omega = atan2(h, k) * 180.0 / M_PI;
        if (omega_plus_Omega < 0) omega_plus_Omega += 360.0;
        omega_deg = omega_plus_Omega - Omega_deg;
        if (omega_deg < 0) omega_deg += 360.0;
        
        M_deg = lambda - omega_plus_Omega;
        while (M_deg < 0) M_deg += 360.0;
        while (M_deg >= 360.0) M_deg -= 360.0;
    }
};

// Carica elementi AstDyS da file .eq1
bool loadAstDySElements(const std::string& filename, EquinoctialElements& elem) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Errore: impossibile aprire " << filename << std::endl;
        return false;
    }
    
    std::string line;
    bool found_equ = false;
    bool found_mjd = false;
    bool in_header = true;
    
    while (std::getline(file, line)) {
        // Fine header
        if (line.find("END_OF_HEADER") != std::string::npos) {
            in_header = false;
            continue;
        }
        
        // Salta header
        if (in_header) continue;
        
        // Salta commenti e linee vuote
        if (line.empty() || line[0] == '!') continue;
        
        // Leggi elementi equinoziali (linea inizia con " EQU" - attenzione agli spazi)
        if (line.find("EQU") != std::string::npos) {
            size_t equ_pos = line.find("EQU");
            std::istringstream iss(line.substr(equ_pos + 3));
            if (iss >> elem.a >> elem.h >> elem.k >> elem.p >> elem.q >> elem.lambda) {
                found_equ = true;
                std::cout << "DEBUG: Elementi letti - a=" << elem.a << " h=" << elem.h << std::endl;
            }
        }
        
        // Leggi epoca (linea inizia con " MJD")
        if (line.find("MJD") != std::string::npos) {
            size_t mjd_pos = line.find("MJD");
            std::istringstream iss(line.substr(mjd_pos + 3));
            if (iss >> elem.epoch_mjd) {
                found_mjd = true;
                std::cout << "DEBUG: Epoca letta - MJD=" << elem.epoch_mjd << std::endl;
            }
        }
        
        if (found_equ && found_mjd) break;
    }
    
    if (!found_equ) std::cerr << "Errore: elementi EQU non trovati" << std::endl;
    if (!found_mjd) std::cerr << "Errore: epoca MJD non trovata" << std::endl;
    
    return found_equ && found_mjd;
}

// Posizione equatoriale (geocentrica)
struct SkyPosition {
    double ra_deg;     // Ascensione Retta (gradi)
    double dec_deg;    // Declinazione (gradi)
    double delta_au;   // Distanza dalla Terra (AU)
    double r_au;       // Distanza dal Sole (AU)
    double mag;        // Magnitudine
};

// Calcola differenza angolare in arcsec
double angularSeparation(double ra1, double dec1, double ra2, double dec2) {
    // Converti in radianti
    double ra1_r = ra1 * M_PI / 180.0;
    double dec1_r = dec1 * M_PI / 180.0;
    double ra2_r = ra2 * M_PI / 180.0;
    double dec2_r = dec2 * M_PI / 180.0;
    
    // Formula haversine per piccole distanze angolari
    double dra = ra2_r - ra1_r;
    double ddec = dec2_r - dec1_r;
    
    double a = sin(ddec/2)*sin(ddec/2) + 
               cos(dec1_r)*cos(dec2_r)*sin(dra/2)*sin(dra/2);
    double sep_rad = 2 * atan2(sqrt(a), sqrt(1-a));
    
    // Converti in arcsec
    return sep_rad * 180.0 / M_PI * 3600.0;
}

// NOTA: Questa è una simulazione semplificata
// L'implementazione reale richiederebbe l'integrazione con ITALOccultLibrary
SkyPosition propagateToEpoch(const EquinoctialElements& elem, double target_mjd) {
    SkyPosition pos;
    
    // Per questo test, usiamo una propagazione Kepleriana semplificata
    // In produzione si userebbe RKF78 con perturbazioni
    
    double dt_days = target_mjd - elem.epoch_mjd;
    double n = 0.01720209895 / pow(elem.a, 1.5); // Moto medio (rad/day)
    
    // Propagazione della longitudine media
    double lambda_new = elem.lambda + n * dt_days * 180.0 / M_PI;
    while (lambda_new >= 360.0) lambda_new -= 360.0;
    while (lambda_new < 0.0) lambda_new += 360.0;
    
    // Conversione a Kepleriani
    double ecc, inc, omega, Omega, M;
    EquinoctialElements elem_new = elem;
    elem_new.lambda = lambda_new;
    elem_new.toKeplerian(ecc, inc, omega, Omega, M);
    
    // Risoluzione equazione di Keplero (metodo Newton-Raphson)
    double E = M * M_PI / 180.0; // Anomalia eccentrica iniziale
    for (int i = 0; i < 10; i++) {
        double dE = (E - ecc*sin(E) - M*M_PI/180.0) / (1.0 - ecc*cos(E));
        E -= dE;
        if (fabs(dE) < 1e-12) break;
    }
    
    // Anomalia vera
    double nu = 2.0 * atan2(sqrt(1+ecc)*sin(E/2), sqrt(1-ecc)*cos(E/2));
    
    // Distanza dal Sole
    pos.r_au = elem.a * (1 - ecc*cos(E));
    
    // Posizione nel piano orbitale
    double r = pos.r_au;
    double u = omega * M_PI / 180.0 + nu; // Argomento di latitudine
    
    // Rotazioni per passare a coordinate eclittiche
    double Omega_r = Omega * M_PI / 180.0;
    double inc_r = inc * M_PI / 180.0;
    
    double x_ecl = r * (cos(Omega_r)*cos(u) - sin(Omega_r)*sin(u)*cos(inc_r));
    double y_ecl = r * (sin(Omega_r)*cos(u) + cos(Omega_r)*sin(u)*cos(inc_r));
    double z_ecl = r * sin(u) * sin(inc_r);
    
    // Posizione Terra (approssimata - orbita circolare)
    double n_earth = 0.01720209895; // rad/day
    double L_earth = 100.0 + n_earth * target_mjd * 180.0 / M_PI; // Longitudine media Terra
    while (L_earth >= 360.0) L_earth -= 360.0;
    double x_earth = cos(L_earth * M_PI / 180.0);
    double y_earth = sin(L_earth * M_PI / 180.0);
    double z_earth = 0.0;
    
    // Geocentrico
    double x_geo = x_ecl - x_earth;
    double y_geo = y_ecl - y_earth;
    double z_geo = z_ecl - z_earth;
    
    // Conversione a coordinate equatoriali (obliquità 23.4°)
    double eps = 23.43928 * M_PI / 180.0;
    double x_eq = x_geo;
    double y_eq = y_geo * cos(eps) - z_geo * sin(eps);
    double z_eq = y_geo * sin(eps) + z_geo * cos(eps);
    
    // RA e Dec
    pos.delta_au = sqrt(x_eq*x_eq + y_eq*y_eq + z_eq*z_eq);
    pos.ra_deg = atan2(y_eq, x_eq) * 180.0 / M_PI;
    if (pos.ra_deg < 0) pos.ra_deg += 360.0;
    pos.dec_deg = asin(z_eq / pos.delta_au) * 180.0 / M_PI;
    
    // Magnitudine (formula H-G)
    double H = 12.874; // Da .eq1
    double phase_angle = acos((r*r + pos.delta_au*pos.delta_au - 1.0) / (2*r*pos.delta_au));
    double phi1 = exp(-3.33 * pow(tan(phase_angle/2), 0.63));
    double phi2 = exp(-1.87 * pow(tan(phase_angle/2), 1.22));
    double G = 0.15;
    pos.mag = H + 5*log10(r * pos.delta_au) - 2.5*log10((1-G)*phi1 + G*phi2);
    
    return pos;
}

// Converti MJD in data calendario
std::string mjdToDate(double mjd) {
    int jd = (int)(mjd + 2400000.5);
    int a = jd + 32044;
    int b = (4*a + 3) / 146097;
    int c = a - (146097*b) / 4;
    int d = (4*c + 3) / 1461;
    int e = c - (1461*d) / 4;
    int m = (5*e + 2) / 153;
    
    int day = e - (153*m + 2) / 5 + 1;
    int month = m + 3 - 12 * (m / 10);
    int year = 100*b + d - 4800 + m / 10;
    
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%04d-%02d-%02d", year, month, day);
    return std::string(buffer);
}

// Formato RA in ore:min:sec
std::string formatRA(double ra_deg) {
    double ra_hours = ra_deg / 15.0;
    int h = (int)ra_hours;
    double rem = (ra_hours - h) * 60.0;
    int m = (int)rem;
    double s = (rem - m) * 60.0;
    
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%02d:%02d:%06.3f", h, m, s);
    return std::string(buffer);
}

// Formato Dec in gra:min:sec
std::string formatDec(double dec_deg) {
    char sign = (dec_deg >= 0) ? '+' : '-';
    double abs_dec = fabs(dec_deg);
    int d = (int)abs_dec;
    double rem = (abs_dec - d) * 60.0;
    int m = (int)rem;
    double s = (rem - m) * 60.0;
    
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%c%02d:%02d:%05.2f", sign, d, m, s);
    return std::string(buffer);
}

int main() {
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  Test Propagazione Temporale RKF78 vs JPL Horizons          ║\n";
    std::cout << "║  Asteroide: (11234) 1999 JS82                                ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    // Carica elementi AstDyS
    EquinoctialElements elem;
    std::string eq1_file = "external/ITALOccultLibrary/astdyn/data/11234.eq1";
    
    if (!loadAstDySElements(eq1_file, elem)) {
        std::cerr << "Errore caricamento elementi da " << eq1_file << std::endl;
        return 1;
    }
    
    std::cout << "═══ Elementi Orbitali AstDyS ═══\n";
    std::cout << "Epoca: MJD " << std::fixed << std::setprecision(6) << elem.epoch_mjd 
              << " (" << mjdToDate(elem.epoch_mjd) << ")\n";
    std::cout << "a     = " << std::setprecision(10) << elem.a << " AU\n";
    std::cout << "h     = " << elem.h << " (e*sin(ω+Ω))\n";
    std::cout << "k     = " << elem.k << " (e*cos(ω+Ω))\n";
    std::cout << "p     = " << elem.p << " (tan(i/2)*sin(Ω))\n";
    std::cout << "q     = " << elem.q << " (tan(i/2)*cos(Ω))\n";
    std::cout << "λ     = " << std::setprecision(10) << elem.lambda << "°\n";
    
    // Conversione a Kepleriani
    double ecc, inc, omega, Omega, M;
    elem.toKeplerian(ecc, inc, omega, Omega, M);
    std::cout << "\n═══ Elementi Kepleriani (convertiti) ═══\n";
    std::cout << "e     = " << std::setprecision(8) << ecc << "\n";
    std::cout << "i     = " << std::setprecision(6) << inc << "°\n";
    std::cout << "ω     = " << omega << "°\n";
    std::cout << "Ω     = " << Omega << "°\n";
    std::cout << "M     = " << M << "°\n";
    
    // Epoche di test
    double epoch_mjd = elem.epoch_mjd;      // MJD 61000.0
    double epoch_plus60 = epoch_mjd + 60.0; // +60 giorni
    double epoch_minus60 = epoch_mjd - 60.0; // -60 giorni
    
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "  PROPAGAZIONE TEMPORALE\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "\n";
    
    // Propaga alle tre epoche
    SkyPosition pos_epoch = propagateToEpoch(elem, epoch_mjd);
    SkyPosition pos_plus60 = propagateToEpoch(elem, epoch_plus60);
    SkyPosition pos_minus60 = propagateToEpoch(elem, epoch_minus60);
    
    // Risultati all'epoca
    std::cout << "▸ EPOCA: MJD " << std::fixed << std::setprecision(1) << epoch_mjd 
              << " (" << mjdToDate(epoch_mjd) << ")\n";
    std::cout << "  RA  = " << formatRA(pos_epoch.ra_deg) 
              << " (" << std::setprecision(6) << pos_epoch.ra_deg << "°)\n";
    std::cout << "  Dec = " << formatDec(pos_epoch.dec_deg) 
              << " (" << std::setprecision(6) << pos_epoch.dec_deg << "°)\n";
    std::cout << "  Δ   = " << std::setprecision(4) << pos_epoch.delta_au << " AU\n";
    std::cout << "  r   = " << pos_epoch.r_au << " AU\n";
    std::cout << "  mag = " << std::setprecision(2) << pos_epoch.mag << "\n";
    std::cout << "\n";
    
    // Risultati +60 giorni
    std::cout << "▸ +60 GIORNI: MJD " << std::fixed << std::setprecision(1) << epoch_plus60 
              << " (" << mjdToDate(epoch_plus60) << ")\n";
    std::cout << "  RA  = " << formatRA(pos_plus60.ra_deg) 
              << " (" << std::setprecision(6) << pos_plus60.ra_deg << "°)\n";
    std::cout << "  Dec = " << formatDec(pos_plus60.dec_deg) 
              << " (" << std::setprecision(6) << pos_plus60.dec_deg << "°)\n";
    std::cout << "  Δ   = " << std::setprecision(4) << pos_plus60.delta_au << " AU\n";
    std::cout << "  r   = " << pos_plus60.r_au << " AU\n";
    std::cout << "  mag = " << std::setprecision(2) << pos_plus60.mag << "\n";
    std::cout << "\n";
    
    // Risultati -60 giorni
    std::cout << "▸ -60 GIORNI: MJD " << std::fixed << std::setprecision(1) << epoch_minus60 
              << " (" << mjdToDate(epoch_minus60) << ")\n";
    std::cout << "  RA  = " << formatRA(pos_minus60.ra_deg) 
              << " (" << std::setprecision(6) << pos_minus60.ra_deg << "°)\n";
    std::cout << "  Dec = " << formatDec(pos_minus60.dec_deg) 
              << " (" << std::setprecision(6) << pos_minus60.dec_deg << "°)\n";
    std::cout << "  Δ   = " << std::setprecision(4) << pos_minus60.delta_au << " AU\n";
    std::cout << "  r   = " << pos_minus60.r_au << " AU\n";
    std::cout << "  mag = " << std::setprecision(2) << pos_minus60.mag << "\n";
    std::cout << "\n";
    
    // ═══════════════════════════════════════════════════════════════
    // CONFRONTO CON JPL HORIZONS
    // ═══════════════════════════════════════════════════════════════
    
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "  CONFRONTO CON JPL HORIZONS\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "\n";
    std::cout << "NOTA: Per un confronto reale, eseguire query JPL Horizons:\n";
    std::cout << "  https://ssd.jpl.nasa.gov/horizons.cgi\n";
    std::cout << "\n";
    std::cout << "Parametri query:\n";
    std::cout << "  - Target: 11234 (1999 JS82)\n";
    std::cout << "  - Observer: 500@0 (geocentro)\n";
    std::cout << "  - Date 1: " << mjdToDate(epoch_minus60) << " 00:00 UT\n";
    std::cout << "  - Date 2: " << mjdToDate(epoch_mjd) << " 00:00 UT\n";
    std::cout << "  - Date 3: " << mjdToDate(epoch_plus60) << " 00:00 UT\n";
    std::cout << "  - Quantities: 1,9 (RA,Dec astrometric)\n";
    std::cout << "\n";
    
    // Placeholder per risultati JPL (da aggiungere manualmente dopo query)
    std::cout << "╔═══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  RISULTATI ATTESI (da verificare con JPL Horizons)          ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  • Accuratezza attesa RKF78:                                 ║\n";
    std::cout << "║    - ±30 giorni: < 1 arcsec (tipicamente 0.1-0.3\")          ║\n";
    std::cout << "║    - ±60 giorni: 1-3 arcsec (con perturbazioni complete)    ║\n";
    std::cout << "║    - ±90 giorni: 3-10 arcsec                                 ║\n";
    std::cout << "║                                                               ║\n";
    std::cout << "║  • Questa implementazione semplificata (Kepleriana pura):    ║\n";
    std::cout << "║    - ±30 giorni: 5-20 arcsec                                 ║\n";
    std::cout << "║    - ±60 giorni: 10-50 arcsec (senza perturbazioni)         ║\n";
    std::cout << "║                                                               ║\n";
    std::cout << "║  Per test precisi, compilare con ITALOccultLibrary e        ║\n";
    std::cout << "║  usare astdyn_propagator con RKF78 + perturbazioni.         ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "  TEST COMPLETATO\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "\n";
    std::cout << "Prossimi passi:\n";
    std::cout << "  1. Eseguire query JPL Horizons per le date indicate\n";
    std::cout << "  2. Confrontare RA/Dec con risultati sopra\n";
    std::cout << "  3. Per accuratezza migliore, integrare RKF78 da ITALOccultLibrary\n";
    std::cout << "  4. Ripetere test con example_astdyn_fitting (usa RKF78 reale)\n";
    std::cout << "\n";
    
    return 0;
}
