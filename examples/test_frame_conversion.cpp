/**
 * @file test_frame_conversion.cpp
 * @brief Test conversione frame equatoriale ↔ eclittico
 * 
 * Verifica che le conversioni siano corrette confrontando con valori noti
 */

#include "ioccultcalc/coordinates.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

constexpr double GAUSS_K = 0.01720209895;  // Costante di Gauss

void testConversion() {
    std::cout << "╔══════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST CONVERSIONE EQUATORIALE ↔ ECLITTICO        ║\n";
    std::cout << "╚══════════════════════════════════════════════════╝\n\n";
    
    // Test 1: Punto sull'equatore celeste dovrebbe avere z_ecl ≈ 0
    std::cout << "Test 1: Punto sull'equatore celeste (Dec=0°)\n";
    std::cout << "──────────────────────────────────────────────────\n";
    double ra1 = 0.0;  // RA = 0h
    double dec1 = 0.0; // Dec = 0°
    
    Vector3D eq1(cos(dec1) * cos(ra1), cos(dec1) * sin(ra1), sin(dec1));
    Vector3D ecl1 = Coordinates::equatorialToEcliptic(eq1);
    
    std::cout << "  Equatoriale: [" << eq1.x << ", " << eq1.y << ", " << eq1.z << "]\n";
    std::cout << "  Eclittico:   [" << ecl1.x << ", " << ecl1.y << ", " << ecl1.z << "]\n";
    std::cout << "  z_eclittico dovrebbe essere piccolo: |z| = " << std::abs(ecl1.z) << "\n";
    
    if (std::abs(ecl1.z) < 0.01) {
        std::cout << "  ✅ CORRETTO\n";
    } else {
        std::cout << "  ❌ ERRORE: z troppo grande!\n";
    }
    
    // Test 2: Polo nord celeste → eclittico
    std::cout << "\nTest 2: Polo nord celeste (Dec=+90°)\n";
    std::cout << "──────────────────────────────────────────────────\n";
    double ra2 = 0.0;
    double dec2 = M_PI / 2.0; // +90°
    
    Vector3D eq2(cos(dec2) * cos(ra2), cos(dec2) * sin(ra2), sin(dec2));
    Vector3D ecl2 = Coordinates::equatorialToEcliptic(eq2);
    
    std::cout << "  Equatoriale: [" << eq2.x << ", " << eq2.y << ", " << eq2.z << "]\n";
    std::cout << "  Eclittico:   [" << ecl2.x << ", " << ecl2.y << ", " << ecl2.z << "]\n";
    
    // Polo nord celeste in eclittico dovrebbe avere:
    // z_ecl = cos(ε) ≈ cos(23.44°) ≈ 0.917
    // y_ecl = sin(ε) ≈ sin(23.44°) ≈ 0.398
    double expected_y = sin(23.4392794444 * M_PI / 180.0);
    double expected_z = cos(23.4392794444 * M_PI / 180.0);
    
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  Atteso:      y ≈ " << expected_y << ", z ≈ " << expected_z << "\n";
    std::cout << "  Errore y:    " << std::abs(ecl2.y - expected_y) << "\n";
    std::cout << "  Errore z:    " << std::abs(ecl2.z - expected_z) << "\n";
    
    if (std::abs(ecl2.y - expected_y) < 0.001 && std::abs(ecl2.z - expected_z) < 0.001) {
        std::cout << "  ✅ CORRETTO\n";
    } else {
        std::cout << "  ❌ ERRORE nella conversione!\n";
    }
    
    // Test 3: Roundtrip equatoriale → eclittico → equatoriale
    std::cout << "\nTest 3: Roundtrip conversion\n";
    std::cout << "──────────────────────────────────────────────────\n";
    double ra3 = 1.234;  // radianti
    double dec3 = 0.567; // radianti
    
    Vector3D eq3(cos(dec3) * cos(ra3), cos(dec3) * sin(ra3), sin(dec3));
    Vector3D ecl3 = Coordinates::equatorialToEcliptic(eq3);
    Vector3D eq3_back = Coordinates::eclipticToEquatorial(ecl3);
    
    std::cout << "  Originale:   [" << eq3.x << ", " << eq3.y << ", " << eq3.z << "]\n";
    std::cout << "  → Eclittico: [" << ecl3.x << ", " << ecl3.y << ", " << ecl3.z << "]\n";
    std::cout << "  → Eq (back): [" << eq3_back.x << ", " << eq3_back.y << ", " << eq3_back.z << "]\n";
    
    double diff = (eq3 - eq3_back).magnitude();
    std::cout << "  Differenza:  " << diff << "\n";
    
    if (diff < 1e-10) {
        std::cout << "  ✅ CORRETTO (roundtrip perfetto)\n";
    } else {
        std::cout << "  ❌ ERRORE: roundtrip non conserva il vettore!\n";
    }
    
    // Test 4: raDecToEclipticUnitVector vs conversione manuale
    std::cout << "\nTest 4: raDecToEclipticUnitVector\n";
    std::cout << "──────────────────────────────────────────────────\n";
    double ra4 = 2.345;
    double dec4 = 0.789;
    
    // Metodo 1: funzione diretta
    Vector3D ecl4_direct = Coordinates::raDecToEclipticUnitVector(ra4, dec4);
    
    // Metodo 2: manuale
    Vector3D eq4(cos(dec4) * cos(ra4), cos(dec4) * sin(ra4), sin(dec4));
    Vector3D ecl4_manual = Coordinates::equatorialToEcliptic(eq4);
    
    std::cout << "  Diretto:  [" << ecl4_direct.x << ", " << ecl4_direct.y << ", " << ecl4_direct.z << "]\n";
    std::cout << "  Manuale:  [" << ecl4_manual.x << ", " << ecl4_manual.y << ", " << ecl4_manual.z << "]\n";
    
    double diff4 = (ecl4_direct - ecl4_manual).magnitude();
    std::cout << "  Differenza: " << diff4 << "\n";
    
    if (diff4 < 1e-10) {
        std::cout << "  ✅ CORRETTO\n";
    } else {
        std::cout << "  ❌ ERRORE: raDecToEclipticUnitVector dà risultato diverso!\n";
    }
    
    // Test 5: Inclinazione da momento angolare
    std::cout << "\nTest 5: Calcolo inclinazione da (r, v)\n";
    std::cout << "──────────────────────────────────────────────────\n";
    
    // Orbita test: i=10° nel piano eclittico
    double i_test = 10.0 * M_PI / 180.0;
    double a_test = 1.5;
    double e_test = 0.2;
    
    // Posizione sul nodo ascendente (Omega=0, omega=0, nu=0)
    Vector3D r_test(a_test * (1 - e_test), 0, 0);
    
    // Velocità: perpendicular to r, inclined by i
    double v_mag = sqrt(GAUSS_K * GAUSS_K / a_test * (1 + e_test) / (1 - e_test));
    Vector3D v_test(0, v_mag * cos(i_test), v_mag * sin(i_test));
    
    // Momento angolare
    Vector3D h_test = r_test.cross(v_test);
    double i_computed = acos(h_test.z / h_test.magnitude()) * 180.0 / M_PI;
    
    std::cout << "  Inclinazione impostata: " << (i_test * 180.0 / M_PI) << "°\n";
    std::cout << "  Inclinazione calcolata: " << i_computed << "°\n";
    std::cout << "  Errore: " << std::abs(i_computed - i_test * 180.0 / M_PI) << "°\n";
    
    if (std::abs(i_computed - i_test * 180.0 / M_PI) < 0.01) {
        std::cout << "  ✅ CORRETTO\n";
    } else {
        std::cout << "  ❌ ERRORE nel calcolo inclinazione!\n";
    }
}

int main() {
    try {
        testConversion();
        
        std::cout << "\n╔══════════════════════════════════════════════════╗\n";
        std::cout << "║  RISULTATO TEST                                  ║\n";
        std::cout << "╚══════════════════════════════════════════════════╝\n";
        std::cout << "\nSe tutti i test sono ✅, la conversione frame è corretta.\n";
        std::cout << "Il bug di inclinazione deve essere altrove (Lambert solver?).\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Eccezione: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
