/**
 * Test: Confronto Propagazione (203) Pompeja
 * IOccultCalc vs OrbFit vs JPL Horizons
 * 
 * Epoca iniziale: MJD 61000.0 (2025-11-21)
 * Epoca finale: MJD 61192.0 (2026-06-01) - 192 giorni dopo
 * 
 * Elementi da OrbFit (eclittica media J2000 - ECLM J2000):
 * - a = 2.7385249933616391 AU
 * - h = e*sin(ω+Ω) = 0.045087089252389
 * - k = e*cos(ω+Ω) = 0.041231297793564
 * - p = tan(i/2)*sin(Ω) = -0.005947645824719
 * - q = tan(i/2)*cos(Ω) = 0.027042352297741
 * - λ = mean longitude = 112.3228065415555°
 */

#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/orbit_propagator.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

// Converte elementi equinoziali (h,k,p,q,λ) a elementi Kepleriani
OrbitalElements equinoctialToKeplerian(double a, double h, double k, double p, double q, double lambda_deg, double mjd_epoch) {
    OrbitalElements elem;
    elem.epoch.jd = mjd_epoch + 2400000.5;  // MJD → JD
    elem.a = a;
    
    // Eccentricità
    elem.e = sqrt(h*h + k*k);
    
    // Inclinazione
    double tan_half_i = sqrt(p*p + q*q);
    elem.i = 2.0 * atan(tan_half_i);
    
    // Longitudine nodo ascendente
    elem.Omega = atan2(p, q);
    if (elem.Omega < 0) elem.Omega += 2.0 * M_PI;
    
    // Argomento pericentro + longitudine nodo (ω + Ω)
    double omega_plus_Omega = atan2(h, k);
    if (omega_plus_Omega < 0) omega_plus_Omega += 2.0 * M_PI;
    
    // Argomento pericentro
    elem.omega = omega_plus_Omega - elem.Omega;
    if (elem.omega < 0) elem.omega += 2.0 * M_PI;
    
    // Anomalia media
    double lambda_rad = lambda_deg * M_PI / 180.0;
    elem.M = lambda_rad - omega_plus_Omega;
    if (elem.M < 0) elem.M += 2.0 * M_PI;
    
    elem.name = "(203) Pompeja";
    elem.designation = "00203";
    
    return elem;
}

int main() {
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "TEST: PROPAGAZIONE (203) POMPEJA\n";
    std::cout << "IOccultCalc vs OrbFit vs JPL Horizons\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";
    
    // Elementi equinoziali da OrbFit (epoch/203.eq1)
    std::cout << "--- ELEMENTI ORBITALI DA ORBFIT ---\n";
    std::cout << "Epoca: MJD 61000.0 TDT (2025-11-21)\n";
    std::cout << "Sistema: Eclittica Media J2000 (ECLM J2000)\n\n";
    
    double a = 2.7385249933616391;
    double h = 0.045087089252389;   // e*sin(ω+Ω)
    double k = 0.041231297793564;   // e*cos(ω+Ω)
    double p = -0.005947645824719;  // tan(i/2)*sin(Ω)
    double q = 0.027042352297741;   // tan(i/2)*cos(Ω)
    double lambda = 112.3228065415555;  // mean longitude (deg)
    double mjd_epoch = 61000.0;
    
    std::cout << "Equinoctial elements:\n";
    std::cout << "  a = " << a << " AU\n";
    std::cout << "  h = " << h << "\n";
    std::cout << "  k = " << k << "\n";
    std::cout << "  p = " << p << "\n";
    std::cout << "  q = " << q << "\n";
    std::cout << "  λ = " << lambda << "°\n\n";
    
    // Crea elementi equinoziali per il propagatore
    EquinoctialElements eq_elem;
    eq_elem.a = a;
    eq_elem.h = h;
    eq_elem.k = k;
    eq_elem.p = p;
    eq_elem.q = q;
    eq_elem.lambda = lambda * M_PI / 180.0;  // Converti in radianti
    eq_elem.epoch.jd = mjd_epoch + 2400000.5;
    
    // Converte a elementi Kepleriani per display
    OrbitalElements elem = equinoctialToKeplerian(a, h, k, p, q, lambda, mjd_epoch);
    
    std::cout << "Keplerian elements (convertiti per display):\n";
    std::cout << "  a = " << elem.a << " AU\n";
    std::cout << "  e = " << elem.e << "\n";
    std::cout << "  i = " << elem.i * 180.0/M_PI << "°\n";
    std::cout << "  Ω = " << elem.Omega * 180.0/M_PI << "°\n";
    std::cout << "  ω = " << elem.omega * 180.0/M_PI << "°\n";
    std::cout << "  M = " << elem.M * 180.0/M_PI << "°\n\n";
    
    // Epoca target: MJD 61192.0 (da 203.oel)
    double mjd_target = 61192.0;
    JulianDate target_epoch;
    target_epoch.jd = mjd_target + 2400000.5;
    
    double dt_days = mjd_target - mjd_epoch;
    
    std::cout << "--- PROPAGAZIONE ---\n";
    std::cout << "Da: MJD " << mjd_epoch << " (JD " << std::fixed << std::setprecision(6) << elem.epoch.jd << ")\n";
    std::cout << "A:  MJD " << mjd_target << " (JD " << target_epoch.jd << ")\n";
    std::cout << "Δt: " << dt_days << " giorni (" << dt_days/365.25 << " anni)\n\n";
    
    // Propaga con RA15 + AST17 (come OrbFit)
    std::cout << "Configurazione propagatore:\n";
    std::cout << "  Integratore: RA15 (Everhart)\n";
    std::cout << "  Step size: 0.05 giorni\n";
    std::cout << "  Perturbazioni: 8 pianeti + AST17\n";
    std::cout << "  Relatività: OFF (match OrbFit default)\n\n";
    
    PropagatorOptions opts;
    opts.integrator = IntegratorType::RA15;
    opts.stepSize = 0.05;  // Come nel .oop
    opts.usePlanetaryPerturbations = true;
    opts.useRelativisticCorrections = false;
    
    OrbitPropagator propagator(opts);
    
    std::cout << "Propagazione in corso...\n";
    OrbitState state_target = propagator.propagate(eq_elem, target_epoch);
    
    std::cout << "\n--- RISULTATI IOCCULTCALC ---\n";
    std::cout << "Posizione @ MJD " << mjd_target << " (J2000 EQUATORIALE):\n";
    std::cout << std::setprecision(12);
    std::cout << "  X = " << state_target.position.x << " AU\n";
    std::cout << "  Y = " << state_target.position.y << " AU\n";
    std::cout << "  Z = " << state_target.position.z << " AU\n";
    std::cout << "Velocità:\n";
    std::cout << "  VX = " << state_target.velocity.x << " AU/day\n";
    std::cout << "  VY = " << state_target.velocity.y << " AU/day\n";
    std::cout << "  VZ = " << state_target.velocity.z << " AU/day\n\n";
    
    // Converti da EQUATORIALE a ECLITTICO per confronto con OrbFit
    constexpr double obliquity = 23.4392911 * M_PI / 180.0;
    double cos_eps = std::cos(obliquity);
    double sin_eps = std::sin(obliquity);
    
    // Rotazione inversa: Equatoriale → Eclittico
    Vector3D pos_ecl, vel_ecl;
    pos_ecl.x = state_target.position.x;
    pos_ecl.y = state_target.position.y * cos_eps + state_target.position.z * sin_eps;
    pos_ecl.z = -state_target.position.y * sin_eps + state_target.position.z * cos_eps;
    
    vel_ecl.x = state_target.velocity.x;
    vel_ecl.y = state_target.velocity.y * cos_eps + state_target.velocity.z * sin_eps;
    vel_ecl.z = -state_target.velocity.y * sin_eps + state_target.velocity.z * cos_eps;
    
    std::cout << "Posizione @ MJD " << mjd_target << " (ECLM J2000 - per confronto OrbFit):\n";
    std::cout << "  X = " << pos_ecl.x << " AU\n";
    std::cout << "  Y = " << pos_ecl.y << " AU\n";
    std::cout << "  Z = " << pos_ecl.z << " AU\n";
    std::cout << "Velocità (ECLM J2000):\n";
    std::cout << "  VX = " << vel_ecl.x << " AU/day\n";
    std::cout << "  VY = " << vel_ecl.y << " AU/day\n";
    std::cout << "  VZ = " << vel_ecl.z << " AU/day\n\n";
    
    // Leggi risultati OrbFit dal file .oel
    std::cout << "--- RISULTATI ORBFIT (da 203.oel) ---\n";
    std::cout << "Elementi equinoziali @ MJD 61192.0:\n";
    std::cout << "  a = 2.7368706317538978 AU\n";
    std::cout << "  h = 0.044799304244679\n";
    std::cout << "  k = 0.041830118242835\n";
    std::cout << "  p = -0.005958715738449\n";
    std::cout << "  q = 0.027053268901305\n";
    std::cout << "  λ = 154.0760142434613°\n\n";
    
    // Converte risultati OrbFit a Kepleriani per confronto
    OrbitalElements elem_orbfit = equinoctialToKeplerian(
        2.7368706317538978,
        0.044799304244679,
        0.041830118242835,
        -0.005958715738449,
        0.027053268901305,
        154.0760142434613,
        mjd_target
    );
    
    std::cout << "Keplerian (convertiti da equinoziali):\n";
    std::cout << "  a = " << elem_orbfit.a << " AU\n";
    std::cout << "  e = " << elem_orbfit.e << "\n";
    std::cout << "  i = " << elem_orbfit.i * 180.0/M_PI << "°\n";
    std::cout << "  Ω = " << elem_orbfit.Omega * 180.0/M_PI << "°\n";
    std::cout << "  ω = " << elem_orbfit.omega * 180.0/M_PI << "°\n";
    std::cout << "  M = " << elem_orbfit.M * 180.0/M_PI << "°\n\n";
    
    // Confronto
    std::cout << "--- CONFRONTO ---\n";
    std::cout << "Differenze elementi orbitali:\n";
    std::cout << "  Δa = " << (elem_orbfit.a - elem.a) * 149597870.7 << " km\n";
    std::cout << "  Δe = " << (elem_orbfit.e - elem.e) << "\n";
    std::cout << "  Δi = " << (elem_orbfit.i - elem.i) * 180.0/M_PI * 3600.0 << " arcsec\n\n";
    
    std::cout << "NOTA: Per confronto posizioni cartesiane, serve propagare\n";
    std::cout << "      gli elementi OrbFit target o query Horizons\n\n";
    
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "Per query JPL Horizons:\n";
    std::cout << "  URL: https://ssd.jpl.nasa.gov/horizons.cgi\n";
    std::cout << "  Target: 203 Pompeja\n";
    std::cout << "  Epoch: 2026-Jun-01 00:00 (MJD 61192.0)\n";
    std::cout << "  Output: State vectors (position & velocity)\n";
    std::cout << "  Center: Solar System Barycenter (@0)\n";
    std::cout << "  Reference: J2000 ecliptic\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    
    return 0;
}
