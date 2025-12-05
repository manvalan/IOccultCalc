/**
 * @file test_17030_vs_jpl.cpp
 * @brief Confronto posizioni asteroide 17030: Nostro codice vs JPL Horizons
 * 
 * Periodo: 26 Novembre 2025 - 29 Novembre 2025
 * Intervallo: 6 ore
 * 
 * Valori di riferimento JPL Horizons:
 * https://ssd.jpl.nasa.gov/horizons.cgi
 * 
 * @date 1 Dicembre 2025
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

// Costanti astronomiche
const double AU = 1.495978707e8;           // km
const double GM_SUN = 1.32712440018e11;    // km³/s²
const double DEG_TO_RAD = M_PI / 180.0;
const double RAD_TO_DEG = 180.0 / M_PI;

struct EquinoctialElements {
    double a;          // AU
    double h, k, p, q; 
    double lambda;     // gradi
    double epoch_mjd;  // MJD TDT
};

struct StateVector {
    double x, y, z;       // km
    double vx, vy, vz;    // km/s
};

struct RADec {
    double ra_deg;
    double dec_deg;
};

struct JPLData {
    double mjd;
    double ra_deg;
    double dec_deg;
    double delta_au;  // Distanza Terra-asteroide
    std::string date_str;
};

// ============================================================================
// DATI JPL HORIZONS PER 17030 (26-29 Nov 2025)
// ============================================================================

std::vector<JPLData> getJPLData() {
    std::vector<JPLData> jpl;
    
    // Target: 17030 Sierks
    // Observer: Geocentric (500@399), ICRF J2000
    // Epoch: 2025-Nov-26 00:00 to 2025-Nov-29 00:00, step 6h
    // Source: JPL Horizons API https://ssd.jpl.nasa.gov/api/horizons.api
    // Retrieved: 2025-Dec-01
    
    jpl.push_back({61005.000000, 73.837820, 20.361170, 2.85, "2025-Nov-26 00:00"});
    jpl.push_back({61005.250000, 73.786160, 20.357610, 2.85, "2025-Nov-26 06:00"});
    jpl.push_back({61005.500000, 73.734370, 20.354040, 2.85, "2025-Nov-26 12:00"});
    jpl.push_back({61005.750000, 73.682460, 20.350460, 2.85, "2025-Nov-26 18:00"});
    jpl.push_back({61006.000000, 73.630440, 20.346870, 2.85, "2025-Nov-27 00:00"});
    jpl.push_back({61006.250000, 73.578300, 20.343270, 2.85, "2025-Nov-27 06:00"});
    jpl.push_back({61006.500000, 73.526060, 20.339670, 2.85, "2025-Nov-27 12:00"});
    jpl.push_back({61006.750000, 73.473710, 20.336050, 2.85, "2025-Nov-27 18:00"});
    jpl.push_back({61007.000000, 73.421250, 20.332430, 2.85, "2025-Nov-28 00:00"});
    jpl.push_back({61007.250000, 73.368700, 20.328800, 2.85, "2025-Nov-28 06:00"});
    jpl.push_back({61007.500000, 73.316050, 20.325160, 2.85, "2025-Nov-28 12:00"});
    jpl.push_back({61007.750000, 73.263300, 20.321520, 2.85, "2025-Nov-28 18:00"});
    jpl.push_back({61008.000000, 73.210470, 20.317860, 2.85, "2025-Nov-29 00:00"});
    
    return jpl;
}

// ============================================================================
// CONVERSIONI (dal test standalone funzionante)
// ============================================================================

void equinoctial_to_keplerian(const EquinoctialElements& eq, 
                               double& a, double& e, double& inc, 
                               double& omega, double& Omega, double& M0) {
    a = eq.a * AU;
    e = std::sqrt(eq.h*eq.h + eq.k*eq.k);
    
    double lp = std::atan2(eq.h, eq.k);
    if (lp < 0) lp += 2*M_PI;
    
    double tan_half_i = std::sqrt(eq.p*eq.p + eq.q*eq.q);
    inc = 2.0 * std::atan(tan_half_i);
    
    Omega = std::atan2(eq.p, eq.q);
    if (Omega < 0) Omega += 2*M_PI;
    
    omega = lp - Omega;
    while (omega < 0) omega += 2*M_PI;
    while (omega >= 2*M_PI) omega -= 2*M_PI;
    
    M0 = (eq.lambda * DEG_TO_RAD) - lp;
    while (M0 < 0) M0 += 2*M_PI;
    while (M0 >= 2*M_PI) M0 -= 2*M_PI;
}

double solve_kepler(double M, double e, double tol = 1e-12) {
    double E = M;
    for (int i = 0; i < 30; i++) {
        double dE = (E - e*sin(E) - M) / (1.0 - e*cos(E));
        E -= dE;
        if (std::abs(dE) < tol) break;
    }
    return E;
}

StateVector keplerian_to_cartesian(double a, double e, double inc, 
                                   double omega, double Omega, double M,
                                   double dt_days = 0.0) {
    // Propaga l'anomalia media
    double n = std::sqrt(GM_SUN / (a*a*a)) / 86400.0; // rad/day
    double M_prop = M + n * dt_days;
    while (M_prop < 0) M_prop += 2*M_PI;
    while (M_prop >= 2*M_PI) M_prop -= 2*M_PI;
    
    double E = solve_kepler(M_prop, e);
    
    double cos_E = cos(E);
    double sin_E = sin(E);
    double sqrt_1_e2 = sqrt(1.0 - e*e);
    
    double x_orb = a * (cos_E - e);
    double y_orb = a * sqrt_1_e2 * sin_E;
    
    double n_sec = sqrt(GM_SUN / (a*a*a));
    double vx_orb = -n_sec * a * sin_E / (1.0 - e*cos_E);
    double vy_orb = n_sec * a * sqrt_1_e2 * cos_E / (1.0 - e*cos_E);
    
    double cos_omega = cos(omega);
    double sin_omega = sin(omega);
    double cos_Omega = cos(Omega);
    double sin_Omega = sin(Omega);
    double cos_inc = cos(inc);
    double sin_inc = sin(inc);
    
    StateVector state;
    state.x = (cos_Omega*cos_omega - sin_Omega*sin_omega*cos_inc) * x_orb +
              (-cos_Omega*sin_omega - sin_Omega*cos_omega*cos_inc) * y_orb;
    state.y = (sin_Omega*cos_omega + cos_Omega*sin_omega*cos_inc) * x_orb +
              (-sin_Omega*sin_omega + cos_Omega*cos_omega*cos_inc) * y_orb;
    state.z = sin_omega*sin_inc * x_orb + cos_omega*sin_inc * y_orb;
    
    state.vx = (cos_Omega*cos_omega - sin_Omega*sin_omega*cos_inc) * vx_orb +
               (-cos_Omega*sin_omega - sin_Omega*cos_omega*cos_inc) * vy_orb;
    state.vy = (sin_Omega*cos_omega + cos_Omega*sin_omega*cos_inc) * vx_orb +
               (-sin_Omega*sin_omega + cos_Omega*cos_omega*cos_inc) * vy_orb;
    state.vz = sin_omega*sin_inc * vx_orb + cos_omega*sin_inc * vy_orb;
    
    return state;
}

void ecliptic_to_equatorial(double x_ecl, double y_ecl, double z_ecl,
                            double& x_eq, double& y_eq, double& z_eq) {
    const double eps = 23.43929111 * DEG_TO_RAD;
    double cos_eps = cos(eps);
    double sin_eps = sin(eps);
    
    x_eq = x_ecl;
    y_eq = y_ecl * cos_eps - z_ecl * sin_eps;
    z_eq = y_ecl * sin_eps + z_ecl * cos_eps;
}

// Posizione Terra (approssimata - per geocentrico)
StateVector earth_position(double mjd) {
    // Semplificazione: orbita circolare Terra
    double days_from_j2000 = mjd - 51544.5;
    double L = 280.460 + 0.9856474 * days_from_j2000; // Longitudine media
    while (L < 0) L += 360.0;
    while (L >= 360.0) L -= 360.0;
    double L_rad = L * DEG_TO_RAD;
    
    const double a_earth = 1.0 * AU; // 1 AU
    
    StateVector earth;
    earth.x = a_earth * cos(L_rad);
    earth.y = a_earth * sin(L_rad);
    earth.z = 0.0; // Trascuriamo inclinazione eclittica
    earth.vx = -29.78 * sin(L_rad); // km/s
    earth.vy = 29.78 * cos(L_rad);
    earth.vz = 0.0;
    
    return earth;
}

RADec state_to_radec_geocentric(const StateVector& state_helio, double mjd) {
    // Posizione Terra
    StateVector earth = earth_position(mjd);
    
    // Posizione geocentrica = asteroide - Terra
    StateVector state_geo;
    state_geo.x = state_helio.x - earth.x;
    state_geo.y = state_helio.y - earth.y;
    state_geo.z = state_helio.z - earth.z;
    
    // Converti da eclittico a equatoriale
    StateVector state_eq;
    ecliptic_to_equatorial(state_geo.x, state_geo.y, state_geo.z,
                          state_eq.x, state_eq.y, state_eq.z);
    
    // Calcola RA/Dec
    double r = sqrt(state_eq.x*state_eq.x + state_eq.y*state_eq.y + state_eq.z*state_eq.z);
    
    RADec result;
    result.dec_deg = asin(state_eq.z / r) * RAD_TO_DEG;
    result.ra_deg = atan2(state_eq.y, state_eq.x) * RAD_TO_DEG;
    
    if (result.ra_deg < 0) result.ra_deg += 360.0;
    
    return result;
}

// ============================================================================
// MAIN
// ============================================================================

int main() {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  CONFRONTO 17030: NOSTRO CODICE vs JPL HORIZONS              ║\n";
    std::cout << "║  Periodo: 26-29 Novembre 2025                                 ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    // Elementi orbitali 17030 da file 17030_astdys.eq1
    // Formato OEF2.0, sistema ECLM J2000, epoca MJD 61000.0 TDT (21 Nov 2025)
    EquinoctialElements eq;
    eq.a = 3.1754732060579491;  // AU
    eq.h = -0.018962873482153;   // e*sin(LP)
    eq.k = -0.041272817500319;   // e*cos(LP)
    eq.p = 0.024582276916386;    // tan(i/2)*sin(LN)
    eq.q = -0.006203125871476;   // tan(i/2)*cos(LN)
    eq.lambda = 74.4674157271250;  // gradi, longitudine media
    eq.epoch_mjd = 61000.0;      // MJD TDT
    
    std::cout << "📋 Elementi orbitali 17030 da 17030_astdys.eq1 (epoca MJD " << eq.epoch_mjd << " = 21 Nov 2025):\n";
    std::cout << "  a = " << eq.a << " AU\n";
    std::cout << "  e = " << std::sqrt(eq.h*eq.h + eq.k*eq.k) << "\n";
    std::cout << "\n";
    
    // Converti in kepleriani
    double a, e, inc, omega, Omega, M;
    equinoctial_to_keplerian(eq, a, e, inc, omega, Omega, M);
    
    std::cout << "🔄 Elementi kepleriani:\n";
    std::cout << "  a = " << (a/AU) << " AU\n";
    std::cout << "  e = " << e << "\n";
    std::cout << "  i = " << (inc * RAD_TO_DEG) << "°\n";
    std::cout << "\n";
    
    // Dati JPL
    auto jpl_data = getJPLData();
    
    std::cout << "╔════════════════════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║ Data/Ora           │ Nostro RA  │ Nostro Dec │ JPL RA     │ JPL Dec    │ ΔRA      │ ΔDec     ║\n";
    std::cout << "╠════════════════════╪════════════╪════════════╪════════════╪════════════╪══════════╪══════════╣\n";
    
    double total_ra_error = 0.0;
    double total_dec_error = 0.0;
    double max_ra_error = 0.0;
    double max_dec_error = 0.0;
    int count = 0;
    
    for (const auto& jpl : jpl_data) {
        // Calcola posizione con il nostro codice
        double dt_days = jpl.mjd - eq.epoch_mjd;
        StateVector state_helio = keplerian_to_cartesian(a, e, inc, omega, Omega, M, dt_days);
        RADec our_radec = state_to_radec_geocentric(state_helio, jpl.mjd);
        
        // Calcola errori
        double delta_ra = our_radec.ra_deg - jpl.ra_deg;
        double delta_dec = our_radec.dec_deg - jpl.dec_deg;
        
        // Converti in arcsec
        double delta_ra_arcsec = delta_ra * 3600.0;
        double delta_dec_arcsec = delta_dec * 3600.0;
        
        // Statistiche
        total_ra_error += std::abs(delta_ra_arcsec);
        total_dec_error += std::abs(delta_dec_arcsec);
        max_ra_error = std::max(max_ra_error, std::abs(delta_ra_arcsec));
        max_dec_error = std::max(max_dec_error, std::abs(delta_dec_arcsec));
        count++;
        
        // Output riga tabella
        std::cout << "║ " << std::left << std::setw(18) << jpl.date_str << " │ ";
        std::cout << std::right << std::fixed << std::setprecision(4);
        std::cout << std::setw(10) << our_radec.ra_deg << "° │ ";
        std::cout << std::setw(10) << our_radec.dec_deg << "° │ ";
        std::cout << std::setw(10) << jpl.ra_deg << "° │ ";
        std::cout << std::setw(10) << jpl.dec_deg << "° │ ";
        std::cout << std::setprecision(2);
        std::cout << std::setw(8) << delta_ra_arcsec << "\" │ ";
        std::cout << std::setw(8) << delta_dec_arcsec << "\" ║\n";
    }
    
    std::cout << "╚════════════════════════════════════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    // Statistiche finali
    double mean_ra_error = total_ra_error / count;
    double mean_dec_error = total_dec_error / count;
    double rms_error = std::sqrt(mean_ra_error*mean_ra_error + mean_dec_error*mean_dec_error);
    
    std::cout << "📊 STATISTICHE ERRORI:\n";
    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
    std::cout << "  Media errore RA:  " << std::setw(8) << std::setprecision(2) << mean_ra_error << " arcsec\n";
    std::cout << "  Media errore Dec: " << std::setw(8) << mean_dec_error << " arcsec\n";
    std::cout << "  RMS combinato:    " << std::setw(8) << rms_error << " arcsec\n";
    std::cout << "  Max errore RA:    " << std::setw(8) << max_ra_error << " arcsec\n";
    std::cout << "  Max errore Dec:   " << std::setw(8) << max_dec_error << " arcsec\n";
    std::cout << "\n";
    
    // Valutazione
    std::cout << "🎯 VALUTAZIONE:\n";
    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
    
    if (rms_error < 0.1) {
        std::cout << "  ✅ ECCELLENTE: Errore < 0.1\" (sub-arcsec precision)\n";
    } else if (rms_error < 1.0) {
        std::cout << "  ✅ OTTIMO: Errore < 1\" (arcsec precision)\n";
    } else if (rms_error < 10.0) {
        std::cout << "  ⚠️  BUONO: Errore < 10\" (accettabile per molte applicazioni)\n";
    } else {
        std::cout << "  ❌ PROBLEMATICO: Errore > 10\" (richiede investigazione)\n";
    }
    
    std::cout << "\n";
    std::cout << "📝 NOTE:\n";
    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
    std::cout << "  • Propagazione con solo gravità solare (no perturbazioni planetarie)\n";
    std::cout << "  • Posizione Terra approssimata (orbita circolare)\n";
    std::cout << "  • Per risultati migliori: integrare RKF78 con perturbazioni\n";
    std::cout << "  • JPL usa: DE441 efemeride + perturbazioni complete\n";
    std::cout << "\n";
    std::cout << "⚡ Per migliorare la precisione:\n";
    std::cout << "  1. Usare AstDyn RKF78 con perturbazioni planetarie\n";
    std::cout << "  2. Usare efemeride DE441 per posizione Terra\n";
    std::cout << "  3. Includere perturbazioni J2 Terra\n";
    std::cout << "  4. Applicare correzioni relativistiche\n";
    std::cout << "\n";
    
    return 0;
}
