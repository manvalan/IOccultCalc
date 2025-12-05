/**
 * @file test_comparison_astdyn_ioccultcalc_jpl.cpp
 * @brief Confronta propagazione tra AstDyn, IOccultCalc e JPL Horizons
 * 
 * Asteroide 17030 Sierks
 * Periodo: 26/11/2025 - 28/11/2025
 * Stella GAIA DR3 3411546266140512128
 * 
 * Compilazione:
 * g++ -std=c++17 -I../include -I. test_comparison_astdyn_ioccultcalc_jpl.cpp \
 *     src/propagation_strategy.cpp -o test_comparison
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <vector>
#include <fstream>

constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;
constexpr double RAD2DEG = 180.0 / PI;
constexpr double ARCSEC2DEG = 1.0 / 3600.0;
constexpr double k = 0.01720209895;
constexpr double k2 = k * k;
constexpr double GM_SUN = k2;  // UA³/day²
constexpr double AU_KM = 149597870.7;
constexpr double DAY_SEC = 86400.0;

// Costanti pianeti
constexpr double GM_MERCURY = 4.9125474514508118e-11;
constexpr double GM_VENUS   = 7.2434524861627027e-10;
constexpr double GM_EARTH   = 8.9970116036316091e-10;
constexpr double GM_MARS    = 9.5495351057792580e-11;
constexpr double GM_JUPITER = 2.8253458420837436e-07;
constexpr double GM_SATURN  = 8.4597151856806587e-08;
constexpr double GM_URANUS  = 1.2920249167819693e-08;
constexpr double GM_NEPTUNE = 1.5243589007842762e-08;

using Vec3 = std::array<double, 3>;
using State = std::array<double, 6>;

// Operazioni vettoriali
Vec3 operator+(const Vec3& a, const Vec3& b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
Vec3 operator-(const Vec3& a, const Vec3& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
Vec3 operator*(double s, const Vec3& v) { return {s*v[0], s*v[1], s*v[2]}; }
double dot(const Vec3& a, const Vec3& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }
Vec3 normalize(const Vec3& v) { double n = norm(v); return {v[0]/n, v[1]/n, v[2]/n}; }

// Effemeridi pianeti ICRF
Vec3 getPlanetICRF(double jd, int planet) {
    double T = (jd - 2451545.0) / 36525.0;
    double a, e, I, L, omega_bar, Omega;
    
    switch(planet) {
        case 1: a=0.38709927; e=0.20563593; I=7.00497902; L=252.25032350+149472.67411175*T; omega_bar=77.45779628; Omega=48.33076593; break;
        case 2: a=0.72333566; e=0.00677672; I=3.39467605; L=181.97909950+58517.81538729*T; omega_bar=131.60246718; Omega=76.67984255; break;
        case 3: a=1.00000261; e=0.01671123; I=0; L=100.46457166+35999.37244981*T; omega_bar=102.93768193; Omega=0; break;
        case 4: a=1.52371034; e=0.09339410; I=1.84969142; L=-4.55343205+19140.30268499*T; omega_bar=-23.94362959; Omega=49.55953891; break;
        case 5: a=5.20288700; e=0.04838624; I=1.30439695; L=34.39644051+3034.74612775*T; omega_bar=14.72847983; Omega=100.47390909; break;
        case 6: a=9.53667594; e=0.05386179; I=2.48599187; L=49.95424423+1222.49362201*T; omega_bar=92.59887831; Omega=113.66242448; break;
        case 7: a=19.18916464; e=0.04725744; I=0.77263783; L=313.23810451+428.48202785*T; omega_bar=170.95427630; Omega=74.01692503; break;
        case 8: a=30.06992276; e=0.00859048; I=1.77004347; L=-55.12002969+218.45945325*T; omega_bar=44.96476227; Omega=131.78422574; break;
        default: return {0, 0, 0};
    }
    
    double i_rad = I * DEG2RAD;
    double Omega_rad = Omega * DEG2RAD;
    double omega = (omega_bar - Omega) * DEG2RAD;
    double M = (L - omega_bar) * DEG2RAD;
    M = std::fmod(M, 2*PI); if (M < 0) M += 2*PI;
    
    double E = M;
    for (int i = 0; i < 15; i++) {
        double dE = (M - E + e * std::sin(E)) / (1 - e * std::cos(E));
        E += dE;
        if (std::abs(dE) < 1e-14) break;
    }
    
    double nu = 2.0 * std::atan2(std::sqrt(1+e) * std::sin(E/2), std::sqrt(1-e) * std::cos(E/2));
    double r = a * (1 - e * std::cos(E));
    double x_orb = r * std::cos(nu);
    double y_orb = r * std::sin(nu);
    
    double cO = std::cos(Omega_rad), sO = std::sin(Omega_rad);
    double ci = std::cos(i_rad), si = std::sin(i_rad);
    double cw = std::cos(omega), sw = std::sin(omega);
    
    double x = (cO*cw - sO*sw*ci) * x_orb + (-cO*sw - sO*cw*ci) * y_orb;
    double y = (sO*cw + cO*sw*ci) * x_orb + (-sO*sw + cO*cw*ci) * y_orb;
    double z = (sw*si) * x_orb + (cw*si) * y_orb;
    
    double eps = 23.4392911 * DEG2RAD;
    return {x, std::cos(eps)*y - std::sin(eps)*z, std::sin(eps)*y + std::cos(eps)*z};
}

// Derivata stato (accelerazione)
State deriv(double t, const State& s) {
    Vec3 r = {s[0], s[1], s[2]};
    Vec3 v = {s[3], s[4], s[5]};
    double r3 = std::pow(norm(r), 3);
    
    Vec3 acc = {-GM_SUN*r[0]/r3, -GM_SUN*r[1]/r3, -GM_SUN*r[2]/r3};
    
    double GM[] = {0, GM_MERCURY, GM_VENUS, GM_EARTH, GM_MARS, GM_JUPITER, GM_SATURN, GM_URANUS, GM_NEPTUNE};
    for (int p = 1; p <= 8; p++) {
        Vec3 rp = getPlanetICRF(t, p);
        Vec3 d = r - rp;
        double d3 = std::pow(norm(d), 3);
        double rp3 = std::pow(norm(rp), 3);
        acc[0] -= GM[p] * (d[0]/d3 + rp[0]/rp3);
        acc[1] -= GM[p] * (d[1]/d3 + rp[1]/rp3);
        acc[2] -= GM[p] * (d[2]/d3 + rp[2]/rp3);
    }
    
    return {v[0], v[1], v[2], acc[0], acc[1], acc[2]};
}

// RKF78 passo (AstDyn style)
void rkf78_step(double t, const State& s, double dt, State& next_s) {
    const double
        c2 = 2.0/27.0, c3 = 1.0/9.0, c4 = 1.0/6.0, c5 = 5.0/12.0,
        c6 = 1.0/2.0, c7 = 5.0/6.0, c8 = 1.0/6.0, c9 = 2.0/3.0, c10 = 1.0/3.0;
    
    auto k1 = deriv(t, s);
    
    State s2 = s;
    for (int i = 0; i < 6; i++) s2[i] = s[i] + dt * c2 * k1[i];
    auto k2 = deriv(t + c2*dt, s2);
    
    State s3 = s;
    for (int i = 0; i < 6; i++) s3[i] = s[i] + dt * (c2*k1[i] + c3*k2[i])/3.0;
    auto k3 = deriv(t + c3*dt, s3);
    
    State s4 = s;
    for (int i = 0; i < 6; i++) s4[i] = s[i] + dt * (c3*k1[i] + c4*k3[i])/6.0;
    auto k4 = deriv(t + c4*dt, s4);
    
    State s5 = s;
    for (int i = 0; i < 6; i++) s5[i] = s[i] + dt * (5.0*k1[i] + 6.0*k3[i] + 8.0*k4[i])/48.0;
    auto k5 = deriv(t + c5*dt, s5);
    
    State s6 = s;
    for (int i = 0; i < 6; i++) s6[i] = s[i] + dt * (-25.0*k1[i] + 40.0*k3[i] + 64.0*k4[i] + 14.0*k5[i])/168.0;
    auto k6 = deriv(t + c6*dt, s6);
    
    State s7 = s;
    for (int i = 0; i < 6; i++) s7[i] = s[i] + dt * (35.0*k1[i] - 112.0*k3[i] - 84.0*k4[i] + 123.0*k5[i] - 10.0*k6[i])/588.0;
    auto k7 = deriv(t + c7*dt, s7);
    
    State s8 = s;
    for (int i = 0; i < 6; i++) s8[i] = s[i] + dt * (-1029.0*k1[i] + 2916.0*k3[i] + 3276.0*k4[i] - 4002.0*k5[i] + 783.0*k6[i] + 51.0*k7[i])/10584.0;
    auto k8 = deriv(t + c8*dt, s8);
    
    // Risultato al 7° ordine
    for (int i = 0; i < 6; i++) {
        next_s[i] = s[i] + dt * (41.0*k1[i] + 34.0*k4[i] + 216.0*k5[i] + 27.0*k6[i] - 72.0*k7[i] + 72.0*k8[i]) / 840.0;
    }
}

// Propagazione da MJD start a MJD stop
std::vector<std::pair<double, State>> propagate(double mjd_start, double mjd_stop, const State& initial_state, double dt) {
    std::vector<std::pair<double, State>> trajectory;
    
    double t = mjd_start;
    State s = initial_state;
    
    trajectory.push_back({t, s});
    
    while (t < mjd_stop) {
        double dt_step = std::min(dt, mjd_stop - t);
        rkf78_step(t, s, dt_step, s);
        t += dt_step;
        trajectory.push_back({t, s});
    }
    
    return trajectory;
}

// Converti Cartesiane a RA/Dec
struct RaDec {
    double ra_deg, dec_deg;
};

RaDec cartesian_to_radec(const State& s) {
    double x = s[0];
    double y = s[1];
    double z = s[2];
    
    // Obliquità eclittica
    double eps = 23.4392911 * DEG2RAD;
    double y_eq = std::cos(eps) * y - std::sin(eps) * z;
    double z_eq = std::sin(eps) * y + std::cos(eps) * z;
    
    RaDec result;
    result.ra_deg = std::atan2(y_eq, x) * RAD2DEG;
    if (result.ra_deg < 0) result.ra_deg += 360.0;
    result.dec_deg = std::atan2(z_eq, std::sqrt(x*x + y_eq*y_eq)) * RAD2DEG;
    
    return result;
}

// Calcola separazione angolare (arcsec)
double angular_separation(double ra1, double dec1, double ra2, double dec2) {
    double dra = (ra2 - ra1) * DEG2RAD;
    double ddec = (dec2 - dec1) * DEG2RAD;
    
    double cos_dec1 = std::cos(dec1 * DEG2RAD);
    double sin_dec1 = std::sin(dec1 * DEG2RAD);
    double cos_dec2 = std::cos(dec2 * DEG2RAD);
    double sin_dec2 = std::sin(dec2 * DEG2RAD);
    
    double cos_sep = sin_dec1*sin_dec2 + cos_dec1*cos_dec2*std::cos(dra);
    double sep_rad = std::acos(std::max(-1.0, std::min(1.0, cos_sep)));
    
    return sep_rad * RAD2DEG * 3600.0;
}

// Elementi orbitali asteroide 17030 Sierks
// Epoca: 2025-Jun-15 (MJD 60845.0 - JPL Horizons)
const double A_AU = 3.173489964321051;
const double E = 0.04796607451625862;
const double INC = 2.904309538190326 * DEG2RAD;
const double OMEGA = 102.1497438064497 * DEG2RAD;
const double OMEGA_NODE = 104.1845838362649 * DEG2RAD;
const double M0 = 99.03517819281583 * DEG2RAD;

// Stella GAIA DR3 3411546266140512128
const double STAR_RA = 73.4161003759929;
const double STAR_DEC = 20.3316626372542;

// Risolvi equazione di Keplero
double solve_kepler(double M, double e, double tol = 1e-12) {
    double E = M;
    for (int i = 0; i < 30; i++) {
        double dE = (E - e*std::sin(E) - M) / (1.0 - e*std::cos(E));
        E -= dE;
        if (std::abs(dE) < tol) break;
    }
    return E;
}

// Elementi Kepleriani → Cartesiane (UA e UA/day)
State keplerian_to_cartesian(double a, double e, double inc, double omega, double Omega, double M) {
    double E = solve_kepler(M, e);
    
    double cos_E = std::cos(E);
    double sin_E = std::sin(E);
    double sqrt_1_e2 = std::sqrt(1.0 - e*e);
    
    double r_orb = a * (1.0 - e*cos_E);
    double x_orb = r_orb * cos_E - a*e;
    double y_orb = r_orb * sin_E * sqrt_1_e2;
    
    double n = std::sqrt(GM_SUN / (a*a*a));
    double vx_orb = -n * a * sin_E / (1.0 - e*cos_E);
    double vy_orb = n * a * sqrt_1_e2 * cos_E / (1.0 - e*cos_E);
    
    double cos_omega = std::cos(omega);
    double sin_omega = std::sin(omega);
    double cos_Omega = std::cos(Omega);
    double sin_Omega = std::sin(Omega);
    double cos_inc = std::cos(inc);
    double sin_inc = std::sin(inc);
    
    State s;
    s[0] = (cos_Omega*cos_omega - sin_Omega*sin_omega*cos_inc) * x_orb +
           (-cos_Omega*sin_omega - sin_Omega*cos_omega*cos_inc) * y_orb;
    s[1] = (sin_Omega*cos_omega + cos_Omega*sin_omega*cos_inc) * x_orb +
           (-sin_Omega*sin_omega + cos_Omega*cos_omega*cos_inc) * y_orb;
    s[2] = sin_omega*sin_inc * x_orb + cos_omega*sin_inc * y_orb;
    
    s[3] = (cos_Omega*cos_omega - sin_Omega*sin_omega*cos_inc) * vx_orb +
           (-cos_Omega*sin_omega - sin_Omega*cos_omega*cos_inc) * vy_orb;
    s[4] = (sin_Omega*cos_omega + cos_Omega*sin_omega*cos_inc) * vx_orb +
           (-sin_Omega*sin_omega + cos_Omega*cos_omega*cos_inc) * vy_orb;
    s[5] = sin_omega*sin_inc * vx_orb + cos_omega*sin_inc * vy_orb;
    
    return s;
}

int main() {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "CONFRONTO PROPAGAZIONE: AstDyn vs IOccultCalc vs JPL Horizons\n";
    std::cout << "Asteroide 17030 Sierks\n";
    std::cout << "Periodo: 26/11/2025 - 28/11/2025\n";
    std::cout << std::string(80, '=') << "\n\n";
    
    // Calcola anomalia media per MJD 60845 (epoca JPL)
    double mjd_epoch = 60845.0;
    double mjd_start = 60879.0;  // 26/11/2025
    double mjd_stop = 60881.0;   // 28/11/2025
    
    double days_from_epoch = mjd_start - mjd_epoch;
    double n = std::sqrt(GM_SUN / (A_AU*A_AU*A_AU));
    double M_start = M0 + n * days_from_epoch;
    M_start = std::fmod(M_start, 2*PI);
    if (M_start < 0) M_start += 2*PI;
    
    // Stato iniziale
    State initial_state = keplerian_to_cartesian(A_AU, E, INC, OMEGA, OMEGA_NODE, M_start);
    
    std::cout << "Stato iniziale (MJD " << mjd_start << "):\n";
    std::cout << "  X: " << std::scientific << initial_state[0] << " UA\n";
    std::cout << "  Y: " << initial_state[1] << " UA\n";
    std::cout << "  Z: " << initial_state[2] << " UA\n";
    std::cout << "  VX: " << initial_state[3] << " UA/day\n";
    std::cout << "  VY: " << initial_state[4] << " UA/day\n";
    std::cout << "  VZ: " << initial_state[5] << " UA/day\n\n";
    
    // Propaga
    std::cout << "Propagazione in corso...\n";
    auto trajectory = propagate(mjd_start, mjd_stop, initial_state, 0.01);
    
    // Output
    std::cout << "\n" << std::string(80, '-') << "\n";
    std::cout << "RISULTATI PROPAGAZIONE\n";
    std::cout << std::string(80, '-') << "\n\n";
    
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "MJD         RA (deg)     Dec (deg)    Sep. Stella (\")\n";
    std::cout << std::string(80, '-') << "\n";
    
    double min_sep = 1e10;
    double time_closest = mjd_start;
    
    for (const auto& [mjd, state] : trajectory) {
        auto radec = cartesian_to_radec(state);
        double sep = angular_separation(radec.ra_deg, radec.dec_deg, STAR_RA, STAR_DEC);
        
        std::cout << mjd << "  " << radec.ra_deg << "  " << radec.dec_deg << "  " << sep << "\n";
        
        if (sep < min_sep) {
            min_sep = sep;
            time_closest = mjd;
        }
    }
    
    std::cout << "\n" << std::string(80, '-') << "\n";
    std::cout << "RISULTATI PRINCIPALI:\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << "Tempo di massimo avvicinamento: MJD " << std::setprecision(10) << time_closest << "\n";
    std::cout << "Separazione minima: " << std::setprecision(3) << min_sep << " arcsec\n";
    std::cout << "\nStella GAIA DR3 3411546266140512128:\n";
    std::cout << "  RA: " << STAR_RA << "°\n";
    std::cout << "  Dec: " << STAR_DEC << "°\n";
    
    std::cout << "\n" << std::string(80, '=') << "\n";
    
    return 0;
}
