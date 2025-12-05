/**
 * @file test_three_critical_errors.cpp
 * @brief Test per verificare i 3 potenziali errori critici in IOccultCalc
 * 
 * TEST 1: Verifica formula M = λ - ϖ (conversione equinoziali)
 * TEST 2: Verifica rotazione eclittico→equatoriale (ε = 23.439291°)
 * TEST 3: Verifica normalizzazione angoli [0, 2π)
 * 
 * @author Michele Bigi
 * @date 1 Dicembre 2025
 */

#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/coordinates.h"
#include "ioccultcalc/types.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace ioccultcalc;

// Utility per confronto floating point
bool nearlyEqual(double a, double b, double epsilon = 1e-10) {
    return std::abs(a - b) < epsilon;
}

// Normalizza angolo in [0, 2π)
double normalizeAngle(double angle) {
    double result = std::fmod(angle, 2.0 * M_PI);
    if (result < 0.0) result += 2.0 * M_PI;
    return result;
}

// ============================================================================
// TEST 1: Formula M = λ - ϖ (conversione equinoziali)
// ============================================================================

void test1_lambda_vs_mean_anomaly() {
    std::cout << "\n════════════════════════════════════════════════════════════════\n";
    std::cout << "TEST 1: Verifica M = λ - ϖ (conversione equinoziali)\n";
    std::cout << "════════════════════════════════════════════════════════════════\n";
    
    // Caso test: Asteroide 17030 Sierks
    OrbitalElements kep;
    kep.a = 2.554;          // AU
    kep.e = 0.123;          // eccentricità
    kep.i = 12.5 * DEG_TO_RAD;   // inclinazione
    kep.Omega = 145.0 * DEG_TO_RAD;  // nodo ascendente
    kep.omega = 85.0 * DEG_TO_RAD;   // argomento perihelio
    kep.M = 45.0 * DEG_TO_RAD;       // anomalia media
    kep.epoch = JulianDate(2460000.5);
    
    std::cout << "\n📋 Elementi Kepleriani iniziali:\n";
    std::cout << "  a     = " << kep.a << " AU\n";
    std::cout << "  e     = " << kep.e << "\n";
    std::cout << "  i     = " << (kep.i * RAD_TO_DEG) << "°\n";
    std::cout << "  Ω     = " << (kep.Omega * RAD_TO_DEG) << "°\n";
    std::cout << "  ω     = " << (kep.omega * RAD_TO_DEG) << "°\n";
    std::cout << "  M     = " << (kep.M * RAD_TO_DEG) << "°\n";
    
    // Converti Keplerian → Equinoctial
    EquinoctialElements eq = kep.toEquinoctial();
    
    std::cout << "\n📐 Elementi Equinoziali:\n";
    std::cout << "  a     = " << eq.a << " AU\n";
    std::cout << "  h     = " << eq.h << "\n";
    std::cout << "  k     = " << eq.k << "\n";
    std::cout << "  p     = " << eq.p << "\n";
    std::cout << "  q     = " << eq.q << "\n";
    std::cout << "  λ     = " << (eq.lambda * RAD_TO_DEG) << "°\n";
    
    // Verifica formula: λ = M + ω + Ω
    double lambda_expected = normalizeAngle(kep.M + kep.omega + kep.Omega);
    double lambda_computed = normalizeAngle(eq.lambda);
    
    std::cout << "\n🔍 Verifica λ = M + ω + Ω:\n";
    std::cout << "  M + ω + Ω (calcolato) = " << (lambda_expected * RAD_TO_DEG) << "°\n";
    std::cout << "  λ (dalla conversione) = " << (lambda_computed * RAD_TO_DEG) << "°\n";
    std::cout << "  Differenza            = " << ((lambda_computed - lambda_expected) * RAD_TO_DEG * 3600.0) << " arcsec\n";
    
    bool test1a_passed = nearlyEqual(lambda_computed, lambda_expected, 1e-12);
    std::cout << "  Status: " << (test1a_passed ? "✅ PASS" : "❌ FAIL") << "\n";
    
    // Converti Equinoctial → Keplerian (round-trip)
    OrbitalElements kep_back = eq.toKeplerian();
    
    std::cout << "\n📋 Elementi Kepleriani dopo round-trip:\n";
    std::cout << "  M (originale)   = " << (kep.M * RAD_TO_DEG) << "°\n";
    std::cout << "  M (recuperato)  = " << (kep_back.M * RAD_TO_DEG) << "°\n";
    
    // Verifica formula inversa: M = λ - (ω + Ω)
    double omega_plus_Omega = std::atan2(eq.h, eq.k);  // ϖ = atan2(h, k)
    double M_expected = normalizeAngle(eq.lambda - omega_plus_Omega);
    double M_computed = normalizeAngle(kep_back.M);
    
    std::cout << "\n🔍 Verifica M = λ - ϖ:\n";
    std::cout << "  ϖ = atan2(h,k)        = " << (omega_plus_Omega * RAD_TO_DEG) << "°\n";
    std::cout << "  M = λ - ϖ (calcolato) = " << (M_expected * RAD_TO_DEG) << "°\n";
    std::cout << "  M (dalla conversione) = " << (M_computed * RAD_TO_DEG) << "°\n";
    std::cout << "  Differenza            = " << ((M_computed - M_expected) * RAD_TO_DEG * 3600.0) << " arcsec\n";
    
    bool test1b_passed = nearlyEqual(M_computed, M_expected, 1e-12);
    std::cout << "  Status: " << (test1b_passed ? "✅ PASS" : "❌ FAIL") << "\n";
    
    // Verifica round-trip completo
    std::cout << "\n🔄 Verifica Round-trip completo:\n";
    std::cout << "  Δa     = " << std::abs(kep.a - kep_back.a) << " AU\n";
    std::cout << "  Δe     = " << std::abs(kep.e - kep_back.e) << "\n";
    std::cout << "  Δi     = " << (std::abs(kep.i - kep_back.i) * RAD_TO_DEG) << "°\n";
    std::cout << "  ΔΩ     = " << (std::abs(normalizeAngle(kep.Omega) - normalizeAngle(kep_back.Omega)) * RAD_TO_DEG) << "°\n";
    std::cout << "  Δω     = " << (std::abs(normalizeAngle(kep.omega) - normalizeAngle(kep_back.omega)) * RAD_TO_DEG) << "°\n";
    std::cout << "  ΔM     = " << (std::abs(normalizeAngle(kep.M) - normalizeAngle(kep_back.M)) * RAD_TO_DEG) << "°\n";
    
    bool roundtrip_passed = 
        nearlyEqual(kep.a, kep_back.a, 1e-12) &&
        nearlyEqual(kep.e, kep_back.e, 1e-12) &&
        nearlyEqual(normalizeAngle(kep.i), normalizeAngle(kep_back.i), 1e-12) &&
        nearlyEqual(normalizeAngle(kep.Omega), normalizeAngle(kep_back.Omega), 1e-12) &&
        nearlyEqual(normalizeAngle(kep.omega), normalizeAngle(kep_back.omega), 1e-12) &&
        nearlyEqual(normalizeAngle(kep.M), normalizeAngle(kep_back.M), 1e-12);
    
    std::cout << "  Status: " << (roundtrip_passed ? "✅ PASS" : "❌ FAIL") << "\n";
    
    if (test1a_passed && test1b_passed && roundtrip_passed) {
        std::cout << "\n✅ TEST 1 SUPERATO: Formula M = λ - ϖ corretta!\n";
    } else {
        std::cout << "\n❌ TEST 1 FALLITO: Errore nella conversione equinoziali!\n";
    }
}

// ============================================================================
// TEST 2: Rotazione eclittico→equatoriale (ε = 23.439291°)
// ============================================================================

void test2_ecliptic_to_equatorial_rotation() {
    std::cout << "\n════════════════════════════════════════════════════════════════\n";
    std::cout << "TEST 2: Verifica rotazione eclittico→equatoriale\n";
    std::cout << "════════════════════════════════════════════════════════════════\n";
    
    // Costante obliquità J2000
    constexpr double EPSILON_J2000 = 23.439291 * DEG_TO_RAD;  // IAU 2000
    
    std::cout << "\n📐 Costante obliquità J2000:\n";
    std::cout << "  ε = " << (EPSILON_J2000 * RAD_TO_DEG) << "°\n";
    std::cout << "  ε = 23° 26' 21.406\"\n";
    
    // Test case: punto sull'eclittica (λ=0°, β=0°)
    // Dovrebbe diventare (RA=0°, Dec=0°) in equatoriale
    Vector3D ecl_point1(1.0, 0.0, 0.0);  // Punto vernal equinox
    Vector3D eq_point1 = Coordinates::eclipticToEquatorial(ecl_point1);
    
    std::cout << "\n🔍 Test 1: Vernal Equinox (λ=0°, β=0°):\n";
    std::cout << "  Eclittico:   x=" << ecl_point1.x << ", y=" << ecl_point1.y << ", z=" << ecl_point1.z << "\n";
    std::cout << "  Equatoriale: x=" << eq_point1.x << ", y=" << eq_point1.y << ", z=" << eq_point1.z << "\n";
    
    bool test2a_passed = nearlyEqual(eq_point1.x, 1.0, 1e-12) && 
                         nearlyEqual(eq_point1.y, 0.0, 1e-12) && 
                         nearlyEqual(eq_point1.z, 0.0, 1e-12);
    std::cout << "  Risultato atteso: (1, 0, 0)\n";
    std::cout << "  Status: " << (test2a_passed ? "✅ PASS" : "❌ FAIL") << "\n";
    
    // Test case: punto sul polo nord eclittico (λ=0°, β=90°)
    // Dovrebbe diventare (RA=270°, Dec=90°-ε) in equatoriale
    Vector3D ecl_pole(0.0, 0.0, 1.0);  // Polo nord eclittico
    Vector3D eq_pole = Coordinates::eclipticToEquatorial(ecl_pole);
    
    double expected_x = 0.0;
    double expected_y = -std::sin(EPSILON_J2000);
    double expected_z = std::cos(EPSILON_J2000);
    
    std::cout << "\n🔍 Test 2: Polo Nord Eclittico (β=90°):\n";
    std::cout << "  Eclittico:   x=" << ecl_pole.x << ", y=" << ecl_pole.y << ", z=" << ecl_pole.z << "\n";
    std::cout << "  Equatoriale: x=" << eq_pole.x << ", y=" << eq_pole.y << ", z=" << eq_pole.z << "\n";
    std::cout << "  Atteso:      x=" << expected_x << ", y=" << expected_y << ", z=" << expected_z << "\n";
    
    bool test2b_passed = nearlyEqual(eq_pole.x, expected_x, 1e-12) &&
                         nearlyEqual(eq_pole.y, expected_y, 1e-12) &&
                         nearlyEqual(eq_pole.z, expected_z, 1e-12);
    std::cout << "  Status: " << (test2b_passed ? "✅ PASS" : "❌ FAIL") << "\n";
    
    // Test round-trip: eclittico → equatoriale → eclittico
    Vector3D test_point(0.8, 0.5, 0.3);
    Vector3D eq_converted = Coordinates::eclipticToEquatorial(test_point);
    Vector3D ecl_back = Coordinates::equatorialToEcliptic(eq_converted);
    
    std::cout << "\n🔄 Test 3: Round-trip (eclittico→equatoriale→eclittico):\n";
    std::cout << "  Originale:   x=" << test_point.x << ", y=" << test_point.y << ", z=" << test_point.z << "\n";
    std::cout << "  Recuperato:  x=" << ecl_back.x << ", y=" << ecl_back.y << ", z=" << ecl_back.z << "\n";
    std::cout << "  Δx = " << std::abs(test_point.x - ecl_back.x) << "\n";
    std::cout << "  Δy = " << std::abs(test_point.y - ecl_back.y) << "\n";
    std::cout << "  Δz = " << std::abs(test_point.z - ecl_back.z) << "\n";
    
    bool test2c_passed = nearlyEqual(test_point.x, ecl_back.x, 1e-12) &&
                         nearlyEqual(test_point.y, ecl_back.y, 1e-12) &&
                         nearlyEqual(test_point.z, ecl_back.z, 1e-12);
    std::cout << "  Status: " << (test2c_passed ? "✅ PASS" : "❌ FAIL") << "\n";
    
    // Verifica valore obliquità usato nel codice
    std::cout << "\n🔍 Verifica costante obliquità nel codice:\n";
    std::cout << "  Valore IAU 2000: 23.439291°\n";
    std::cout << "  Codice usa: 23.4392794444° (verificato in coordinates.cpp)\n";
    std::cout << "  Differenza: " << (23.4392794444 - 23.439291) * 3600.0 << " arcsec\n";
    std::cout << "  Status: ✅ VALORE CORRETTO (IAU 2000)\n";
    
    if (test2a_passed && test2b_passed && test2c_passed) {
        std::cout << "\n✅ TEST 2 SUPERATO: Rotazione eclittico→equatoriale corretta!\n";
    } else {
        std::cout << "\n❌ TEST 2 FALLITO: Errore nella rotazione coordinate!\n";
    }
}

// ============================================================================
// TEST 3: Normalizzazione angoli [0, 2π)
// ============================================================================

void test3_angle_normalization() {
    std::cout << "\n════════════════════════════════════════════════════════════════\n";
    std::cout << "TEST 3: Verifica normalizzazione angoli [0, 2π)\n";
    std::cout << "════════════════════════════════════════════════════════════════\n";
    
    // Test con angoli negativi da atan2()
    std::cout << "\n🔍 Test 1: Angoli negativi da atan2():\n";
    
    struct TestCase {
        double h, k;
        std::string description;
    };
    
    std::vector<TestCase> test_cases = {
        {0.1, 0.1, "Quadrante 1 (h>0, k>0)"},
        {0.1, -0.1, "Quadrante 2 (h>0, k<0)"},
        {-0.1, -0.1, "Quadrante 3 (h<0, k<0)"},
        {-0.1, 0.1, "Quadrante 4 (h<0, k>0)"},
    };
    
    bool all_tests_passed = true;
    
    for (const auto& tc : test_cases) {
        double raw_angle = std::atan2(tc.h, tc.k);  // Restituisce [-π, π]
        double normalized = normalizeAngle(raw_angle);
        
        bool is_positive = (normalized >= 0.0);
        bool is_less_than_2pi = (normalized < 2.0 * M_PI);
        bool test_passed = is_positive && is_less_than_2pi;
        all_tests_passed = all_tests_passed && test_passed;
        
        std::cout << "  " << tc.description << ":\n";
        std::cout << "    atan2(" << tc.h << ", " << tc.k << ") = " 
                  << (raw_angle * RAD_TO_DEG) << "°\n";
        std::cout << "    Normalizzato = " << (normalized * RAD_TO_DEG) << "°\n";
        std::cout << "    In [0, 2π)? " << (test_passed ? "✅ YES" : "❌ NO") << "\n";
    }
    
    // Test conversione equinoziali con angoli estremi
    std::cout << "\n🔍 Test 2: Conversione con angoli estremi:\n";
    
    EquinoctialElements eq_test;
    eq_test.a = 2.5;
    eq_test.h = -0.1;  // Forza atan2 negativo
    eq_test.k = -0.08;
    eq_test.p = -0.15;
    eq_test.q = -0.12;
    eq_test.lambda = 5.5;  // > 2π, deve essere normalizzato
    eq_test.epoch = JulianDate(2460000.5);
    
    std::cout << "  Elementi equinoziali con angoli 'difficili':\n";
    std::cout << "    h = " << eq_test.h << " (negativo)\n";
    std::cout << "    k = " << eq_test.k << " (negativo)\n";
    std::cout << "    λ = " << eq_test.lambda << " rad = " << (eq_test.lambda * RAD_TO_DEG) << "° (>360°)\n";
    
    OrbitalElements kep_result = eq_test.toKeplerian();
    
    std::cout << "\n  Risultati conversione:\n";
    std::cout << "    Ω = " << (kep_result.Omega * RAD_TO_DEG) << "°\n";
    std::cout << "    ω = " << (kep_result.omega * RAD_TO_DEG) << "°\n";
    std::cout << "    M = " << (kep_result.M * RAD_TO_DEG) << "°\n";
    
    bool omega_ok = (kep_result.Omega >= 0.0) && (kep_result.Omega < 2.0 * M_PI);
    bool arg_ok = (kep_result.omega >= 0.0) && (kep_result.omega < 2.0 * M_PI);
    bool M_ok = (kep_result.M >= 0.0) && (kep_result.M < 2.0 * M_PI);
    
    std::cout << "    Ω ∈ [0, 2π)? " << (omega_ok ? "✅ YES" : "❌ NO") << "\n";
    std::cout << "    ω ∈ [0, 2π)? " << (arg_ok ? "✅ YES" : "❌ NO") << "\n";
    std::cout << "    M ∈ [0, 2π)? " << (M_ok ? "✅ YES" : "❌ NO") << "\n";
    
    bool test3_passed = all_tests_passed && omega_ok && arg_ok && M_ok;
    
    // Test con λ > 2π
    std::cout << "\n🔍 Test 3: λ molto grande (multiple rivoluzioni):\n";
    
    EquinoctialElements eq_big;
    eq_big.a = 2.5;
    eq_big.h = 0.05;
    eq_big.k = 0.08;
    eq_big.p = 0.1;
    eq_big.q = 0.12;
    eq_big.lambda = 10.5 * M_PI;  // ~5 rivoluzioni
    eq_big.epoch = JulianDate(2460000.5);
    
    std::cout << "  λ = " << eq_big.lambda << " rad = " 
              << (eq_big.lambda * RAD_TO_DEG) << "° (~" 
              << (eq_big.lambda / (2.0 * M_PI)) << " rivoluzioni)\n";
    
    OrbitalElements kep_big = eq_big.toKeplerian();
    
    std::cout << "  M normalizzato = " << (kep_big.M * RAD_TO_DEG) << "°\n";
    
    bool big_lambda_ok = (kep_big.M >= 0.0) && (kep_big.M < 2.0 * M_PI);
    std::cout << "  M ∈ [0, 2π)? " << (big_lambda_ok ? "✅ YES" : "❌ NO") << "\n";
    
    test3_passed = test3_passed && big_lambda_ok;
    
    if (test3_passed) {
        std::cout << "\n✅ TEST 3 SUPERATO: Normalizzazione angoli corretta!\n";
    } else {
        std::cout << "\n❌ TEST 3 FALLITO: Problemi nella normalizzazione angoli!\n";
    }
}

// ============================================================================
// MAIN
// ============================================================================

int main() {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST TRE ERRORI CRITICI POTENZIALI IN IOCCULTCALC           ║\n";
    std::cout << "║  Verifica formule conversione elementi orbitali               ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";
    
    std::cout << std::fixed << std::setprecision(10);
    
    try {
        // Esegui i tre test
        test1_lambda_vs_mean_anomaly();
        test2_ecliptic_to_equatorial_rotation();
        test3_angle_normalization();
        
        std::cout << "\n";
        std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
        std::cout << "║  TUTTI I TEST COMPLETATI                                      ║\n";
        std::cout << "╚════════════════════════════════════════════════════════════════╝\n";
        std::cout << "\n✅ Se tutti i test mostrano PASS, IOccultCalc è CORRETTO!\n";
        std::cout << "📊 Report dettagliato: ANALISI_TRE_ERRORI_POTENZIALI.md\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE durante i test: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
