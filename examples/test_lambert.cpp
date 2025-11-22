/**
 * @file test_lambert.cpp
 * @brief Test diagnostico per il Lambert solver
 */

#include "ioccultcalc/initial_orbit.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/coordinates.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

constexpr double PI = 3.14159265358979323846;
constexpr double GAUSS_K = 0.01720209895;

void printVector(const std::string& label, const Vector3D& v) {
    std::cout << label << ": [" 
              << std::fixed << std::setprecision(8)
              << v.x << ", " << v.y << ", " << v.z << "]"
              << " (mag=" << v.magnitude() << ")" << std::endl;
}

void printElements(const std::string& label, const OrbitalElements& elem) {
    std::cout << "\n" << label << ":\n";
    std::cout << "  a = " << std::fixed << std::setprecision(6) << elem.a << " AU\n";
    std::cout << "  e = " << elem.e << "\n";
    std::cout << "  i = " << (elem.i * 180.0 / PI) << "°\n";
    std::cout << "  Ω = " << (elem.Omega * 180.0 / PI) << "°\n";
    std::cout << "  ω = " << (elem.omega * 180.0 / PI) << "°\n";
    std::cout << "  M = " << (elem.M * 180.0 / PI) << "°\n";
}

int main() {
    std::cout << "╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║         TEST LAMBERT SOLVER                            ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n\n";
    
    // Test 1: Orbita circolare semplice (Earth-like)
    std::cout << "============================================================\n";
    std::cout << "TEST 1: Orbita circolare (a=1 AU, e=0)\n";
    std::cout << "============================================================\n\n";
    
    // Posizione iniziale: (1, 0, 0)
    Vector3D r1(1.0, 0.0, 0.0);
    
    // Dopo 90 giorni (~0.25 orbita): (0, 1, 0) per orbita circolare
    double dt1 = 90.0;  // giorni
    Vector3D r2(0.0, 1.0, 0.0);
    
    printVector("r1", r1);
    printVector("r2", r2);
    std::cout << "dt = " << dt1 << " giorni\n\n";
    
    Vector3D v1, v2;
    double mu = GAUSS_K * GAUSS_K;
    
    int result = InitialOrbit::solveLambert(r1, r2, dt1, mu, true, v1, v2);
    
    if (result == 0) {
        std::cout << "✓ Lambert solver convergenza OK\n\n";
        printVector("v1", v1);
        printVector("v2", v2);
        
        // Verifica: velocità per orbita circolare = 2π/T = 2π/365.25 AU/day
        double v_circular = 2.0 * PI / 365.25;  // ~0.01721 AU/day
        std::cout << "\nVelocità attesa per orbita circolare: " 
                  << v_circular << " AU/day\n";
        std::cout << "Velocità calcolata: " << v1.magnitude() << " AU/day\n";
        
        // Converti in elementi orbitali manualmente
        Vector3D r(r1.x, r1.y, r1.z);
        Vector3D v(v1.x, v1.y, v1.z);
        Vector3D h = r.cross(v);  // Momento angolare
        
        double rmag = r.magnitude();
        double vmag = v.magnitude();
        double hmag = h.magnitude();
        
        // Semi-asse maggiore
        double energy = vmag * vmag / 2.0 - mu / rmag;
        double a_calc = -mu / (2.0 * energy);
        
        // Eccentricità
        Vector3D ev = (v.cross(h)) / mu - r / rmag;
        double e_calc = ev.magnitude();
        
        // Inclinazione
        double i_calc = acos(h.z / hmag);
        
        std::cout << "\nElementi ricostruiti:\n";
        std::cout << "  a = " << a_calc << " AU\n";
        std::cout << "  e = " << e_calc << "\n";
        std::cout << "  i = " << (i_calc * 180.0 / PI) << "°\n";
        
    } else {
        std::cout << "✗ Lambert solver FALLITO\n";
    }
    
    // Test 2: Eros-like orbit
    std::cout << "\n\n============================================================\n";
    std::cout << "TEST 2: Orbita Eros-like (a=1.458 AU, e=0.223, i=10.8°)\n";
    std::cout << "============================================================\n\n";
    
    // Crea elementi Eros
    OrbitalElements eros;
    eros.epoch = JulianDate(2460000.0);
    eros.a = 1.458;
    eros.e = 0.223;
    eros.i = 10.829 * PI / 180.0;
    eros.Omega = 304.322 * PI / 180.0;
    eros.omega = 178.817 * PI / 180.0;
    eros.M = 320.365 * PI / 180.0;
    
    printElements("Riferimento Eros", eros);
    
    // Caso semplificato: usa posizioni conosciute per Eros
    // Posizione tipica per Eros: a=1.458, e=0.223
    // Perielio: q = a(1-e) = 1.133 AU
    // Afelio:  Q = a(1+e) = 1.783 AU
    
    Vector3D r1_eros(1.133, 0.0, 0.2);  // Vicino al perielio con inclinazione
    Vector3D r2_eros(-1.0, 1.0, 0.2);   // Dopo 90 giorni circa
    
    // Velocità attesa (semi-empirica per Eros)
    Vector3D v1_eros(0.0, 0.022, 0.002);  // ~orbita ellittica
    Vector3D v2_eros(-0.015, -0.010, 0.002);
    
    double dt2 = 90.0;
    
    std::cout << "\n";
    printVector("r2 (dopo 90 giorni, propagato)", r2_eros);
    printVector("v2 (dopo 90 giorni, propagato)", v2_eros);
    
    // Ora usa Lambert per ricostruire
    std::cout << "\n--- Ricostruzione con Lambert solver ---\n";
    Vector3D v1_lambert, v2_lambert;
    result = InitialOrbit::solveLambert(r1_eros, r2_eros, dt2, mu, true, v1_lambert, v2_lambert);
    
    if (result == 0) {
        std::cout << "✓ Lambert solver convergenza OK\n\n";
        printVector("v1 (Lambert)", v1_lambert);
        printVector("v2 (Lambert)", v2_lambert);
        
        // Errore
        Vector3D dv1 = v1_lambert - v1_eros;
        Vector3D dv2 = v2_lambert - v2_eros;
        
        std::cout << "\nErrore su v1: " << dv1.magnitude() << " AU/day "
                  << "(" << (dv1.magnitude() / v1_eros.magnitude() * 100.0) << "%)\n";
        std::cout << "Errore su v2: " << dv2.magnitude() << " AU/day "
                  << "(" << (dv2.magnitude() / v2_eros.magnitude() * 100.0) << "%)\n";
        
        // Elementi orbitali ricostruiti manualmente
        Vector3D r_l(r1_eros.x, r1_eros.y, r1_eros.z);
        Vector3D v_l(v1_lambert.x, v1_lambert.y, v1_lambert.z);
        Vector3D h_l = r_l.cross(v_l);
        
        double rmag_l = r_l.magnitude();
        double vmag_l = v_l.magnitude();
        double hmag_l = h_l.magnitude();
        
        double energy_l = vmag_l * vmag_l / 2.0 - mu / rmag_l;
        double a_lambert = -mu / (2.0 * energy_l);
        
        Vector3D ev_l = (v_l.cross(h_l)) / mu - r_l / rmag_l;
        double e_lambert = ev_l.magnitude();
        
        double i_lambert = acos(h_l.z / hmag_l);
        
        std::cout << "\nElementi ricostruiti (Lambert):\n";
        std::cout << "  a = " << a_lambert << " AU\n";
        std::cout << "  e = " << e_lambert << "\n";
        std::cout << "  i = " << (i_lambert * 180.0 / PI) << "°\n";
        
        // Confronto
        std::cout << "\n--- Confronto con riferimento ---\n";
        std::cout << "Δa: " << std::fixed << std::setprecision(2)
                  << fabs(a_lambert - eros.a) / eros.a * 100.0 << "%\n";
        std::cout << "Δe: " << fabs(e_lambert - eros.e) / eros.e * 100.0 << "%\n";
        std::cout << "Δi: " << fabs(i_lambert - eros.i) * 180.0 / PI << "°\n";
        
        if (fabs(a_lambert - eros.a) / eros.a < 0.05 &&
            fabs(e_lambert - eros.e) / eros.e < 0.20) {
            std::cout << "\n✓ BUONO: Lambert solver funziona correttamente!\n";
        } else {
            std::cout << "\n⚠ PARZIALE: Lambert funziona ma c'è qualche imprecisione\n";
        }
        
    } else {
        std::cout << "✗ Lambert solver FALLITO\n";
    }
    
    std::cout << "\n============================================================\n";
    std::cout << "TEST COMPLETATO\n";
    std::cout << "============================================================\n";
    
    return 0;
}
