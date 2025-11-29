/**
 * @file replicate_orbfit.cpp
 * @brief Replica ESATTA del workflow OrbFit per confronto diretto
 * 
 * Workflow OrbFit:
 * 1. prop2b: Propagazione 2-body (Keplerian) per time-of-flight iniziale
 * 2. RA15: Integrazione numerica con perturbazioni complete
 * 3. Coordinate: Conversione Eclittica J2000 → Equatoriale J2000
 * 
 * Questo programma replica ESATTAMENTE ogni passo per identificare discrepanze.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <ioccultcalc/orbit_propagator.h>
#include <ioccultcalc/orbital_elements.h>
#include <ioccultcalc/coordinates.h>
#include <ioccultcalc/types.h>

using namespace ioccultcalc;

// ============================================================================
// STEP 1: Conversione Equatoriale → Eclittica (come fa OrbFit)
// ============================================================================

/**
 * @brief Ruota vettore da Equatoriale J2000 a Eclittica J2000
 * Matrice rotazione: R_x(-ε) dove ε = 23.439291111° (obliquità J2000)
 */
void equatorialToEcliptic(const Vector3D& eq, Vector3D& ecl) {
    constexpr double eps = 23.439291111 * M_PI / 180.0;  // Obliquità J2000
    constexpr double cos_eps = 0.917482062069182;         // cos(23.439291111°)
    constexpr double sin_eps = 0.397777155931914;         // sin(23.439291111°)
    
    ecl.x = eq.x;
    ecl.y = eq.y * cos_eps + eq.z * sin_eps;
    ecl.z = -eq.y * sin_eps + eq.z * cos_eps;
}

/**
 * @brief Ruota vettore da Eclittica J2000 a Equatoriale J2000
 * Matrice rotazione: R_x(+ε) dove ε = 23.439291111°
 */
void eclipticToEquatorial(const Vector3D& ecl, Vector3D& eq) {
    constexpr double eps = 23.439291111 * M_PI / 180.0;
    constexpr double cos_eps = 0.917482062069182;
    constexpr double sin_eps = 0.397777155931914;
    
    eq.x = ecl.x;
    eq.y = ecl.y * cos_eps - ecl.z * sin_eps;
    eq.z = ecl.y * sin_eps + ecl.z * cos_eps;
}

// ============================================================================
// STEP 2: Propagazione 2-body (prop2b di OrbFit)
// ============================================================================

/**
 * @brief Risolve equazione di Keplero con metodo Newton-Raphson
 * M = E - e*sin(E)  →  trova E dato M ed e
 */
double solveKeplerEquation(double M_rad, double e, int max_iter = 20, double tol = 1e-12) {
    double E = M_rad;  // Prima approssimazione
    
    for (int i = 0; i < max_iter; i++) {
        double f = E - e * sin(E) - M_rad;
        double df = 1.0 - e * cos(E);
        double dE = -f / df;
        E += dE;
        
        if (std::abs(dE) < tol) {
            break;
        }
    }
    
    return E;
}

/**
 * @brief Propagazione Keplerian pura (come prop2b di OrbFit)
 * 
 * Input:  Elementi orbitali (a, e, i, Ω, ω, M) @ epoch_0
 * Output: Posizione e velocità cartesiane @ epoch_target
 * 
 * Reference: OrbFit/src/suit/orb_els.f90, subroutine prop2b
 */
OrbitState propagate2Body(const OrbitalElements& elem0, const JulianDate& target_epoch) {
    const double mu = 0.01720209895 * 0.01720209895;  // k² (Gauss constant)²
    
    // Calcola tempo trascorso (giorni)
    double dt_days = target_epoch.toMJD() - elem0.epoch.toMJD();
    
    // Mean motion (movimenti medi giornalieri)
    double n = sqrt(mu / (elem0.a * elem0.a * elem0.a));  // rad/day
    
    // Anomalia media alla nuova epoca
    double M0_rad = elem0.M * M_PI / 180.0;
    double M_rad = M0_rad + n * dt_days;
    
    // Normalizza M a [0, 2π)
    M_rad = fmod(M_rad, 2.0 * M_PI);
    if (M_rad < 0) M_rad += 2.0 * M_PI;
    
    // Risolvi equazione di Keplero: M = E - e*sin(E)
    double E_rad = solveKeplerEquation(M_rad, elem0.e);
    
    // Coordinate nel piano orbitale (referenza perielio)
    double cos_E = cos(E_rad);
    double sin_E = sin(E_rad);
    double sqrt_1_e2 = sqrt(1.0 - elem0.e * elem0.e);
    
    double x_orb = elem0.a * (cos_E - elem0.e);
    double y_orb = elem0.a * sqrt_1_e2 * sin_E;
    
    double r = elem0.a * (1.0 - elem0.e * cos_E);
    
    double vx_orb = -elem0.a * n * sin_E / (1.0 - elem0.e * cos_E);
    double vy_orb = elem0.a * n * sqrt_1_e2 * cos_E / (1.0 - elem0.e * cos_E);
    
    // Rotazione dal piano orbitale a eclittica J2000
    // Matrice: R_z(-Ω) * R_x(-i) * R_z(-ω)
    double i_rad = elem0.i * M_PI / 180.0;
    double Omega_rad = elem0.Omega * M_PI / 180.0;
    double omega_rad = elem0.omega * M_PI / 180.0;
    
    double cos_i = cos(i_rad);
    double sin_i = sin(i_rad);
    double cos_Omega = cos(Omega_rad);
    double sin_Omega = sin(Omega_rad);
    double cos_omega = cos(omega_rad);
    double sin_omega = sin(omega_rad);
    
    // Posizione in eclittica J2000
    double P11 = cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i;
    double P12 = -cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i;
    double P21 = sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i;
    double P22 = -sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i;
    double P31 = sin_omega * sin_i;
    double P32 = cos_omega * sin_i;
    
    OrbitState result;
    result.position.x = P11 * x_orb + P12 * y_orb;
    result.position.y = P21 * x_orb + P22 * y_orb;
    result.position.z = P31 * x_orb + P32 * y_orb;
    
    result.velocity.x = P11 * vx_orb + P12 * vy_orb;
    result.velocity.y = P21 * vx_orb + P22 * vy_orb;
    result.velocity.z = P31 * vx_orb + P32 * vy_orb;
    
    result.epoch = target_epoch;
    
    return result;
}

// ============================================================================
// MAIN: Test completo replicando OrbFit
// ============================================================================

int main() {
    std::cout << "=================================================================\n";
    std::cout << "REPLICA ESATTA DEL WORKFLOW ORBFIT\n";
    std::cout << "=================================================================\n\n";
    
    try {
        // ====================================================================
        // CONFIGURAZIONE TEST: Eros (433) - Elementi FITTED da OrbFit
        // ====================================================================
        
        std::cout << "CONFIGURAZIONE:\n";
        std::cout << "  Asteroide: (433) Eros\n";
        std::cout << "  Elementi: OrbFit CORRECTED (fit su 16119 osservazioni)\n";
        std::cout << "  Residual norm: 0.326 arcsec\n";
        std::cout << "  Epoca iniziale: MJD 55918.807374655 (2011-Dec-23)\n";
        std::cout << "  Epoca target:   MJD 60300.0 (2024-Jan-01)\n";
        std::cout << "  Δt = 4381.19 giorni (~12 anni)\n\n";
        
        // ====================================================================
        // STEP 1: Elementi iniziali (da 433.olg - Corrected elements)
        // ====================================================================
        
        std::cout << "=================================================================\n";
        std::cout << "STEP 1: Elementi iniziali FITTED (Equatoriale J2000)\n";
        std::cout << "=================================================================\n\n";
        
        double epoch0_mjd = 55918.807374655;
        JulianDate epoch0(epoch0_mjd + 2400000.5);
        
        // State vector da OrbFit (Equatoriale J2000)
        Vector3D r0_eq, v0_eq;
        r0_eq.x = -2.23090281851654E-01;  // AU
        r0_eq.y =  1.12387595809716E+00;
        r0_eq.z =  8.60668959807618E-02;
        v0_eq.x = -1.68320119058993E-02;  // AU/d
        v0_eq.y = -4.36024765332142E-03;
        v0_eq.z = -3.12870435158407E-03;
        
        std::cout << "Posizione @ MJD " << std::fixed << std::setprecision(9) << epoch0_mjd << ":\n";
        std::cout << "  r_eq = (" << std::setprecision(14)
                  << r0_eq.x << ", "
                  << r0_eq.y << ", "
                  << r0_eq.z << ") AU\n";
        std::cout << "  v_eq = ("
                  << v0_eq.x << ", "
                  << v0_eq.y << ", "
                  << v0_eq.z << ") AU/d\n";
        
        double r0_mag = sqrt(r0_eq.x*r0_eq.x + r0_eq.y*r0_eq.y + r0_eq.z*r0_eq.z);
        std::cout << "  |r| = " << std::setprecision(12) << r0_mag << " AU\n\n";
        
        // ====================================================================
        // STEP 2: Conversione Equatoriale → Eclittica (come OrbFit)
        // ====================================================================
        
        std::cout << "=================================================================\n";
        std::cout << "STEP 2: Conversione Equatoriale → Eclittica J2000\n";
        std::cout << "=================================================================\n\n";
        
        Vector3D r0_ecl, v0_ecl;
        equatorialToEcliptic(r0_eq, r0_ecl);
        equatorialToEcliptic(v0_eq, v0_ecl);
        
        std::cout << "  r_ecl = (" << std::setprecision(14)
                  << r0_ecl.x << ", "
                  << r0_ecl.y << ", "
                  << r0_ecl.z << ") AU\n";
        std::cout << "  v_ecl = ("
                  << v0_ecl.x << ", "
                  << v0_ecl.y << ", "
                  << v0_ecl.z << ") AU/d\n\n";
        
        // ====================================================================
        // STEP 3: Calcolo elementi orbitali da stato cartesiano
        // ====================================================================
        
        std::cout << "=================================================================\n";
        std::cout << "STEP 3: Calcolo elementi orbitali (Eclittica J2000)\n";
        std::cout << "=================================================================\n\n";
        
        OrbitState state0_ecl;
        state0_ecl.position = r0_ecl;
        state0_ecl.velocity = v0_ecl;
        state0_ecl.epoch = epoch0;
        
        OrbitalElements elem0 = cartesianToOrbitalElements(
            r0_ecl, v0_ecl, 
            0.01720209895 * 0.01720209895,  // mu
            epoch0
        );
        
        std::cout << "Elementi Kepleriani @ MJD " << epoch0_mjd << ":\n";
        std::cout << "  a = " << std::setprecision(12) << elem0.a << " AU\n";
        std::cout << "  e = " << elem0.e << "\n";
        std::cout << "  i = " << std::setprecision(8) << elem0.i << "°\n";
        std::cout << "  Ω = " << elem0.Omega << "°\n";
        std::cout << "  ω = " << elem0.omega << "°\n";
        std::cout << "  M = " << elem0.M << "°\n\n";
        
        // ====================================================================
        // STEP 4A: Propagazione 2-body (prop2b di OrbFit)
        // ====================================================================
        
        std::cout << "=================================================================\n";
        std::cout << "STEP 4A: Propagazione 2-body (prop2b come OrbFit)\n";
        std::cout << "=================================================================\n\n";
        
        double target_mjd = 60300.0;
        JulianDate target_epoch(target_mjd + 2400000.5);
        
        std::cout << "Propagando con metodo 2-body (solo Sole)...\n";
        auto state_2body = propagate2Body(elem0, target_epoch);
        
        Vector3D r_2body_eq, v_2body_eq;
        eclipticToEquatorial(state_2body.position, r_2body_eq);
        eclipticToEquatorial(state_2body.velocity, v_2body_eq);
        
        std::cout << "Risultato 2-body @ MJD " << target_mjd << " (Equatoriale J2000):\n";
        std::cout << "  r_eq = (" << std::setprecision(14)
                  << r_2body_eq.x << ", "
                  << r_2body_eq.y << ", "
                  << r_2body_eq.z << ") AU\n";
        std::cout << "  v_eq = ("
                  << v_2body_eq.x << ", "
                  << v_2body_eq.y << ", "
                  << v_2body_eq.z << ") AU/d\n\n";
        
        // ====================================================================
        // STEP 4B: Propagazione con perturbazioni (RA15 di OrbFit)
        // ====================================================================
        
        std::cout << "=================================================================\n";
        std::cout << "STEP 4B: Propagazione con perturbazioni (RA15 come OrbFit)\n";
        std::cout << "=================================================================\n\n";
        
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RA15;
        opts.stepSize = 0.1;  // 0.1 giorni (default OrbFit)
        opts.usePlanetaryPerturbations = true;
        opts.useRelativisticCorrections = true;
        
        std::cout << "Configurazione integratore:\n";
        std::cout << "  Tipo: RA15 (Everhart 1985)\n";
        std::cout << "  Step: 0.1 giorni\n";
        std::cout << "  Perturbazioni:\n";
        std::cout << "    - Pianeti: 8 (Mercury → Neptune)\n";
        std::cout << "    - Relatività: Einstein-Infeld-Hoffmann\n";
        std::cout << "  Note: Asteroidi massivi AST17 gestiti internamente\n\n";
        
        std::cout << "Propagando " << (target_mjd - epoch0_mjd) << " giorni...\n";
        std::cout << "⏳ Tempo stimato: ~30-60 secondi...\n\n";
        
        OrbitPropagator prop(opts);
        auto state_full = prop.propagate(state0_ecl, target_epoch);
        
        Vector3D r_full_eq, v_full_eq;
        eclipticToEquatorial(state_full.position, r_full_eq);
        eclipticToEquatorial(state_full.velocity, v_full_eq);
        
        std::cout << "Risultato RA15 @ MJD " << target_mjd << " (Equatoriale J2000):\n";
        std::cout << "  r_eq = (" << std::setprecision(14)
                  << r_full_eq.x << ", "
                  << r_full_eq.y << ", "
                  << r_full_eq.z << ") AU\n";
        std::cout << "  v_eq = ("
                  << v_full_eq.x << ", "
                  << v_full_eq.y << ", "
                  << v_full_eq.z << ") AU/d\n\n";
        
        // ====================================================================
        // STEP 5: Confronto con OrbFit
        // ====================================================================
        
        std::cout << "=================================================================\n";
        std::cout << "STEP 5: Confronto con OrbFit\n";
        std::cout << "=================================================================\n\n";
        
        // Valori OrbFit @ MJD 60300.0 (da file 433.oel)
        Vector3D r_orbfit, v_orbfit;
        r_orbfit.x = 1.39933263576780E+00;
        r_orbfit.y = 4.05376107746071E-01;
        r_orbfit.z = 2.64807319940301E-01;
        v_orbfit.x = -6.88234140936674E-03;
        v_orbfit.y = 1.22202237086118E-02;
        v_orbfit.z = 2.28870868886193E-04;
        
        std::cout << "OrbFit @ MJD " << target_mjd << " (Equatoriale J2000):\n";
        std::cout << "  r_eq = (" << std::setprecision(14)
                  << r_orbfit.x << ", "
                  << r_orbfit.y << ", "
                  << r_orbfit.z << ") AU\n";
        std::cout << "  v_eq = ("
                  << v_orbfit.x << ", "
                  << v_orbfit.y << ", "
                  << v_orbfit.z << ") AU/d\n\n";
        
        // ====================================================================
        // STEP 6: Analisi differenze
        // ====================================================================
        
        std::cout << "=================================================================\n";
        std::cout << "STEP 6: Analisi differenze\n";
        std::cout << "=================================================================\n\n";
        
        // Differenza 2-body vs OrbFit
        double dx_2b = (r_2body_eq.x - r_orbfit.x) * 149597870.7;  // AU → km
        double dy_2b = (r_2body_eq.y - r_orbfit.y) * 149597870.7;
        double dz_2b = (r_2body_eq.z - r_orbfit.z) * 149597870.7;
        double diff_2body = sqrt(dx_2b*dx_2b + dy_2b*dy_2b + dz_2b*dz_2b);
        
        // Differenza RA15 vs OrbFit
        double dx_ra15 = (r_full_eq.x - r_orbfit.x) * 149597870.7;
        double dy_ra15 = (r_full_eq.y - r_orbfit.y) * 149597870.7;
        double dz_ra15 = (r_full_eq.z - r_orbfit.z) * 149597870.7;
        double diff_ra15 = sqrt(dx_ra15*dx_ra15 + dy_ra15*dy_ra15 + dz_ra15*dz_ra15);
        
        std::cout << "A) Propagazione 2-body:\n";
        std::cout << "   |Δr| = " << std::setprecision(1) << diff_2body << " km\n";
        std::cout << "   Componenti: ΔX=" << dx_2b << ", ΔY=" << dy_2b << ", ΔZ=" << dz_2b << " km\n";
        std::cout << "   → Le perturbazioni planetarie causano ~" << diff_2body << " km di differenza\n\n";
        
        std::cout << "B) Propagazione RA15 (con perturbazioni):\n";
        std::cout << "   |Δr| = " << std::setprecision(3) << diff_ra15 << " km\n";
        std::cout << "   Componenti: ΔX=" << dx_ra15 << ", ΔY=" << dy_ra15 << ", ΔZ=" << dz_ra15 << " km\n";
        std::cout << "   Error/year = " << (diff_ra15 / 12.0) << " km/year\n\n";
        
        // ====================================================================
        // CONCLUSIONI
        // ====================================================================
        
        std::cout << "=================================================================\n";
        std::cout << "CONCLUSIONI\n";
        std::cout << "=================================================================\n\n";
        
        if (diff_ra15 < 10.0) {
            std::cout << "✅ ECCELLENTE: Differenza < 10 km dopo 12 anni\n";
            std::cout << "   IOccultCalc replica accuratamente OrbFit!\n\n";
        } else if (diff_ra15 < 100.0) {
            std::cout << "✓ BUONO: Differenza < 100 km dopo 12 anni\n";
            std::cout << "   Precisione adeguata per predizioni occultazioni\n\n";
        } else if (diff_ra15 < 1000.0) {
            std::cout << "⚠️  ACCETTABILE: Differenza " << diff_ra15 << " km dopo 12 anni\n";
            std::cout << "   Verificare configurazione perturbazioni\n\n";
        } else {
            std::cout << "❌ PROBLEMA: Differenza " << diff_ra15 << " km dopo 12 anni\n";
            std::cout << "   Discrepanza significativa - investigare:\n";
            std::cout << "   1. Set asteroidi massivi (AST17)\n";
            std::cout << "   2. Parametri integratore RA15\n";
            std::cout << "   3. Correzioni relativistiche\n\n";
        }
        
        double improvement = (diff_2body - diff_ra15) / diff_2body * 100.0;
        std::cout << "Miglioramento rispetto a 2-body: " << std::setprecision(1) 
                  << improvement << "%\n";
        std::cout << "Perturbazioni riducono l'errore di: " << (diff_2body - diff_ra15) << " km\n\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << std::endl;
        return 1;
    }
}
