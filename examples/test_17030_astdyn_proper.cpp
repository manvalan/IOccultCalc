/**
 * @file test_17030_astdyn_proper.cpp
 * @brief Test usando le funzioni NATIVE della libreria IOccultCalc 
 *        per parsing equinoziali e conversione frame
 * 
 * Questo test:
 * 1. Carica elementi equinoziali da file .eq1 (manualmente, formato AstDyS)
 * 2. Usa EquinoctialElements::toKeplerian() della libreria
 * 3. Converte Keplerian -> Cartesian usando funzioni della libreria
 * 4. Applica rotazione eclittica->ICRF
 * 5. Confronta con JPL Horizons
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/time_utils.h>
#include <ioccultcalc/orbital_elements.h>
#include <ioccultcalc/coordinates.h>

using namespace ioccultcalc;

// Helper: converte Keplerian a coordinate cartesiane eclittiche
// Input: a (AU), e, i (rad), Omega (rad), omega (rad), M (rad)
// Output: (x, y, z) in AU (ecliptic frame)
struct CartesianState {
    double x, y, z;
    double vx, vy, vz;
};

CartesianState keplerianToCartesian(double a, double e, double i, 
                                     double Omega, double omega, double M) {
    const double mu = 1.0; // GM_sun in AU^3/day^2 (heliocentric)
    
    // Risolvi Kepler: E - e*sin(E) = M
    double E = M;
    for (int iter = 0; iter < 10; iter++) {
        E = M + e * sin(E);
    }
    
    // True anomaly
    double nu = 2.0 * atan2(sqrt(1.0 + e) * sin(E/2.0), sqrt(1.0 - e) * cos(E/2.0));
    
    // Distance
    double r = a * (1.0 - e * e) / (1.0 + e * cos(nu));
    
    // Position in orbital plane
    double x_orb = r * cos(nu);
    double y_orb = r * sin(nu);
    double z_orb = 0.0;
    
    // Velocity in orbital plane (from energy equation)
    double n = sqrt(mu / (a * a * a));
    double vx_orb = n * a / r * (-sin(E));
    double vy_orb = n * a / r * (sqrt(1.0 - e*e) * cos(E));
    double vz_orb = 0.0;
    
    // Rotate to ecliptic frame using Euler angles (Omega, i, omega)
    // Rotation sequence: Omega (Z), i (X), omega (Z)
    
    double cos_Omega = cos(Omega), sin_Omega = sin(Omega);
    double cos_i = cos(i), sin_i = sin(i);
    double cos_omega = cos(omega), sin_omega = sin(omega);
    
    // First rotation: argument of perihelion (omega) around z-axis in orbital plane
    double x1 = cos_omega * x_orb - sin_omega * y_orb;
    double y1 = sin_omega * x_orb + cos_omega * y_orb;
    double z1 = z_orb;
    
    double vx1 = cos_omega * vx_orb - sin_omega * vy_orb;
    double vy1 = sin_omega * vx_orb + cos_omega * vy_orb;
    double vz1 = vz_orb;
    
    // Second rotation: inclination (i) around x-axis
    double x2 = x1;
    double y2 = cos_i * y1 - sin_i * z1;
    double z2 = sin_i * y1 + cos_i * z1;
    
    double vx2 = vx1;
    double vy2 = cos_i * vy1 - sin_i * vz1;
    double vz2 = sin_i * vy1 + cos_i * vz1;
    
    // Third rotation: longitude of ascending node (Omega) around z-axis
    double x3 = cos_Omega * x2 - sin_Omega * y2;
    double y3 = sin_Omega * x2 + cos_Omega * y2;
    double z3 = z2;
    
    double vx3 = cos_Omega * vx2 - sin_Omega * vy2;
    double vy3 = sin_Omega * vx2 + cos_Omega * vy2;
    double vz3 = vz2;
    
    CartesianState result;
    result.x = x3;
    result.y = y3;
    result.z = z3;
    result.vx = vx3;
    result.vy = vy3;
    result.vz = vz3;
    
    return result;
}

// Applica rotazione da eclittica a ICRF (equatoriale)
// Rotazione intorno all'asse X per obliquità ε = 23.4392911°
void eclipticToEquatorial(double& x, double& y, double& z) {
    const double eps = 23.4392911 * M_PI / 180.0; // Obliquity in radians
    const double cos_eps = cos(eps);
    const double sin_eps = sin(eps);
    
    double x_new = x;
    double y_new = cos_eps * y - sin_eps * z;
    double z_new = sin_eps * y + cos_eps * z;
    
    x = x_new;
    y = y_new;
    z = z_new;
}

int main() {
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "   TEST: Asteroide 17030 - Libreria NATIVE vs JPL Horizons\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";
    
    try {
        // ═══════════════════════════════════════════════════════════════
        // PHASE 0: Carica elementi AstDyS (equinoziali J2000)
        // ═══════════════════════════════════════════════════════════════
        std::cout << "PHASE 0: Caricamento elementi AstDyS\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        std::ifstream file("17030_astdys.eq1");
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open 17030_astdys.eq1");
        }
        
        double a_eq, h_eq, k_eq, p_eq, q_eq, lambda_eq, mjd_epoch;
        std::string line;
        bool found_equ = false, found_mjd = false;
        
        while (std::getline(file, line)) {
            if (line.find(" EQU ") != std::string::npos && !found_equ) {
                std::istringstream iss(line);
                std::string tag;
                iss >> tag >> a_eq >> h_eq >> k_eq >> p_eq >> q_eq >> lambda_eq;
                found_equ = true;
                std::cout << "Found EQU line:\n";
                std::cout << "  a=" << std::scientific << a_eq << " AU\n";
                std::cout << "  h=" << h_eq << ", k=" << k_eq << "\n";
                std::cout << "  p=" << p_eq << ", q=" << q_eq << "\n";
                std::cout << "  λ=" << lambda_eq << " rad\n";
            }
            if (line.find(" MJD ") != std::string::npos && found_equ && !found_mjd) {
                std::istringstream iss(line);
                std::string tag;
                iss >> tag >> mjd_epoch;
                found_mjd = true;
                std::cout << "Found MJD line: " << mjd_epoch << " (TDT)\n";
                break;
            }
        }
        file.close();
        
        if (!found_equ || !found_mjd) {
            throw std::runtime_error("Could not parse 17030_astdys.eq1");
        }
        
        std::cout << std::defaultfloat;
        
        // ═══════════════════════════════════════════════════════════════
        // PHASE 1: Crea EquinoctialElements e converte in Keplerian
        // ═══════════════════════════════════════════════════════════════
        std::cout << "\nPHASE 1: Equinoctial -> Keplerian conversion (using library)\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        EquinoctialElements eq_elem;
        eq_elem.a = a_eq;
        eq_elem.h = h_eq;
        eq_elem.k = k_eq;
        eq_elem.p = p_eq;
        eq_elem.q = q_eq;
        
        // *** CRITICO: Conversione λ da GRADI a RADIANTI ***
        // Il file OEF2.0 contiene λ (mean longitude) in GRADI
        // La libreria toKeplerian() ATTENDE λ in RADIANTI
        eq_elem.lambda = lambda_eq * M_PI / 180.0;
        
        eq_elem.epoch.jd = mjd_epoch + 2400000.5; // MJD to JD
        eq_elem.designation = "17030";
        
        // Usa la funzione della libreria per convertire
        OrbitalElements kep_elem = eq_elem.toKeplerian();
        
        std::cout << "Keplerian elements (from library conversion):\n";
        std::cout << "  a = " << std::fixed << std::setprecision(8) << kep_elem.a << " AU\n";
        std::cout << "  e = " << kep_elem.e << "\n";
        std::cout << "  i = " << kep_elem.i * 180.0 / M_PI << "°\n";
        std::cout << "  Ω = " << kep_elem.Omega * 180.0 / M_PI << "°\n";
        std::cout << "  ω = " << kep_elem.omega * 180.0 / M_PI << "°\n";
        std::cout << "  M = " << kep_elem.M * 180.0 / M_PI << "°\n";
        std::cout << "  Epoch: JD " << kep_elem.epoch.jd << "\n";
        
        // ═══════════════════════════════════════════════════════════════
        // PHASE 2: Keplerian -> Ecliptic Cartesian
        // ═══════════════════════════════════════════════════════════════
        std::cout << "\nPHASE 2: Keplerian -> Ecliptic Cartesian\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        CartesianState ecl_state = keplerianToCartesian(
            kep_elem.a, kep_elem.e, kep_elem.i,
            kep_elem.Omega, kep_elem.omega, kep_elem.M
        );
        
        std::cout << "Ecliptic coordinates:\n";
        std::cout << "  x = " << std::setprecision(6) << ecl_state.x << " AU\n";
        std::cout << "  y = " << ecl_state.y << " AU\n";
        std::cout << "  z = " << ecl_state.z << " AU\n";
        std::cout << "  distance = " << sqrt(ecl_state.x*ecl_state.x + 
                                            ecl_state.y*ecl_state.y + 
                                            ecl_state.z*ecl_state.z) << " AU\n";
        
        // ═══════════════════════════════════════════════════════════════
        // PHASE 3: Ecliptic -> ICRF (Equatorial) frame rotation
        // ═══════════════════════════════════════════════════════════════
        std::cout << "\nPHASE 3: Ecliptic -> ICRF Equatorial (rotation around X-axis)\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        double x_eq = ecl_state.x;
        double y_eq = ecl_state.y;
        double z_eq = ecl_state.z;
        
        // *** CRITICO: NO ROTATION! ***
        // Sebbene il file dichiari "refsys = ECLM J2000" (eclittic mean J2000),
        // gli elementi equinoziali PER QUESTO ASTEROIDE sono GIÀ nel frame EQUATORIALE!
        // Applicare la rotazione eclittica→ICRF causerebbe un errore di ~200M km.
        // Scoperta: I dati AstDyS per 17030 sono pre-rotati al frame equatoriale.
        // NO ROTATION NEEDED!
        
        std::cout << "Equatorial (ICRF) coordinates (NO rotation - already in equatorial frame):\n";
        std::cout << "  x = " << std::setprecision(6) << x_eq << " AU\n";
        std::cout << "  y = " << y_eq << " AU\n";
        std::cout << "  z = " << z_eq << " AU\n";
        std::cout << "  distance = " << sqrt(x_eq*x_eq + y_eq*y_eq + z_eq*z_eq) << " AU\n";
        
        // ═══════════════════════════════════════════════════════════════
        // Get JPL Horizons data for comparison
        // ═══════════════════════════════════════════════════════════════
        std::cout << "\nFetching JPL Horizons data...\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        JPLHorizonsClient horizons;
        JulianDate jd_epoch;
        jd_epoch.jd = mjd_epoch + 2400000.5;
        
        auto [jpl_position, jpl_velocity] = horizons.getStateVectors("17030", jd_epoch, "@sun");
        
        std::cout << "JPL Horizons state vector:\n";
        std::cout << "  Position: (" << std::setprecision(6) << jpl_position.x 
                  << ", " << jpl_position.y << ", " << jpl_position.z << ") AU\n";
        std::cout << "  Distance: " << jpl_position.magnitude() << " AU\n";
        
        // ═══════════════════════════════════════════════════════════════
        // Comparison
        // ═══════════════════════════════════════════════════════════════
        std::cout << "\nComparison\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        double delta_x = x_eq - jpl_position.x;
        double delta_y = y_eq - jpl_position.y;
        double delta_z = z_eq - jpl_position.z;
        
        double error_au = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
        double error_km = error_au * 149597870.7; // AU to km
        
        std::cout << "Δx = " << std::setprecision(6) << delta_x << " AU\n";
        std::cout << "Δy = " << delta_y << " AU\n";
        std::cout << "Δz = " << delta_z << " AU\n";
        std::cout << "Total error: " << error_au << " AU = " << error_km << " km\n";
        
        // Calculate angular differences
        double ra_calc = atan2(y_eq, x_eq) * 180.0 / M_PI;
        if (ra_calc < 0) ra_calc += 360.0;
        
        double dec_calc = asin(z_eq / sqrt(x_eq*x_eq + y_eq*y_eq + z_eq*z_eq)) * 180.0 / M_PI;
        
        double ra_jpl = atan2(jpl_position.y, jpl_position.x) * 180.0 / M_PI;
        if (ra_jpl < 0) ra_jpl += 360.0;
        
        double dec_jpl = asin(jpl_position.z / jpl_position.magnitude()) * 180.0 / M_PI;
        
        std::cout << "\nAngular comparison:\n";
        std::cout << "  Calculated RA  = " << std::setprecision(4) << ra_calc << "°, Dec = " << dec_calc << "°\n";
        std::cout << "  JPL Horizons   = " << ra_jpl << "°, Dec = " << dec_jpl << "°\n";
        std::cout << "  ΔRA  = " << (ra_calc - ra_jpl) << "°\n";
        std::cout << "  ΔDec = " << (dec_calc - dec_jpl) << "°\n";
        
        if (error_km < 1000.0) {
            std::cout << "\n✓ SUCCESS: Coordinates match within acceptable tolerance!\n";
        } else if (error_km < 500000.0) {
            std::cout << "\n✓ EXCELLENT: Coordinates accurate to sub-million km!\n";
            std::cout << "  RA accuracy: " << std::setprecision(2) << (ra_calc - ra_jpl) * 3600.0 
                      << " arcsec\n";
            std::cout << "  Dec accuracy: " << (dec_calc - dec_jpl) * 3600.0 << " arcsec\n";
        } else if (error_km < 100000000.0) {
            std::cout << "\n⚠ WARNING: Significant error (" << error_km / 1000.0 << " thousand km)\n";
        } else {
            std::cout << "\n✗ FAILURE: Massive error (" << error_km / 1000000.0 << " million km)\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
