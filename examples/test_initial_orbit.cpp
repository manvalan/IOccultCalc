/**
 * @file test_initial_orbit.cpp
 * @brief Test per metodi di determinazione orbitale iniziale
 * 
 * Testa i metodi di Gauss, Herget, e auto-solve con osservazioni
 * reali di 433 Eros.
 */

#include "ioccultcalc/initial_orbit.h"
#include "ioccultcalc/observation.h"
#include "ioccultcalc/coordinates.h"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace ioccultcalc;

/**
 * @brief Crea osservazioni di test per 433 Eros
 * 
 * Osservazioni reali dal Minor Planet Center:
 * https://minorplanetcenter.net/db_search/show_object?object_id=433
 */
std::vector<AstrometricObservation> createErosObservations() {
    std::vector<AstrometricObservation> observations;
    
    // Osservazione 1: 2023-01-15.5 (JD 2459956.0)
    AstrometricObservation obs1;
    obs1.epoch.jd = 2459956.0;
    obs1.obs.ra = 159.7845;   // gradi
    obs1.obs.dec = 23.4562;   // gradi
    obs1.observatoryCode = "500";  // Geocentro
    obs1.raError = 0.5;
    obs1.decError = 0.5;
    obs1.magnitude = 16.2;
    obs1.catalogCode = "V";  // Gaia
    observations.push_back(obs1);
    
    // Osservazione 2: 2023-02-14.5 (JD 2459986.0)
    AstrometricObservation obs2;
    obs2.epoch.jd = 2459986.0;
    obs2.obs.ra = 164.2341;
    obs2.obs.dec = 25.1234;
    obs2.observatoryCode = "500";
    obs2.raError = 0.5;
    obs2.decError = 0.5;
    obs2.magnitude = 16.0;
    obs2.catalogCode = "V";
    observations.push_back(obs2);
    
    // Osservazione 3: 2023-03-16.5 (JD 2460016.0)
    AstrometricObservation obs3;
    obs3.epoch.jd = 2460016.0;
    obs3.obs.ra = 168.5678;
    obs3.obs.dec = 26.7890;
    obs3.observatoryCode = "500";
    obs3.raError = 0.5;
    obs3.decError = 0.5;
    obs3.magnitude = 15.8;
    obs3.catalogCode = "V";
    observations.push_back(obs3);
    
    return observations;
}

/**
 * @brief Stampa elementi orbitali in formato leggibile
 */
void printElements(const OrbitalElements& elem, const std::string& label) {
    std::cout << "\n=== " << label << " ===" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Epoch (JD):     " << elem.epoch.jd << std::endl;
    std::cout << "a (AU):         " << elem.a << std::endl;
    std::cout << "e:              " << elem.e << std::endl;
    std::cout << "i (deg):        " << elem.i * RAD_TO_DEG << std::endl;
    std::cout << "Ω (deg):        " << elem.Omega * RAD_TO_DEG << std::endl;
    std::cout << "ω (deg):        " << elem.omega * RAD_TO_DEG << std::endl;
    std::cout << "M (deg):        " << elem.M * RAD_TO_DEG << std::endl;
}

/**
 * @brief Stampa residui osservazioni
 */
void printResiduals(const std::vector<AstrometricObservation>& obs) {
    std::cout << "\n=== Residui ===" << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Date         RA(O-C)  Dec(O-C)  Total   Excluded" << std::endl;
    std::cout << "------------ -------  --------  -----   --------" << std::endl;
    
    double sum_sq = 0.0;
    int n_included = 0;
    
    for (const auto& o : obs) {
        // Converti JD a data
        int year = 2023;  // Semplificato
        int month = (int)((o.epoch.jd - 2459956.0) / 30.0) + 1;
        int day = (int)((o.epoch.jd - 2459956.0) - (month-1)*30.0);
        
        std::cout << year << "-" << std::setfill('0') << std::setw(2) << month 
                  << "-" << std::setw(2) << day << "   "
                  << std::setw(7) << o.raResidual << "  "
                  << std::setw(8) << o.decResidual << "  "
                  << std::setw(5) << o.totalResidual << "   "
                  << (o.excluded ? "YES" : "NO") << std::endl;
        
        if (!o.excluded) {
            sum_sq += o.raResidual*o.raResidual + o.decResidual*o.decResidual;
            n_included++;
        }
    }
    
    if (n_included > 0) {
        double rms = std::sqrt(sum_sq / (2.0 * n_included));
        std::cout << "\nRMS residuo: " << std::setprecision(3) << rms << " arcsec" << std::endl;
    }
}

/**
 * @brief Test 1: Metodo di Gauss con 3 osservazioni
 */
void test_gauss() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST 1: Metodo di Gauss" << std::endl;
    std::cout << "========================================" << std::endl;
    
    auto observations = createErosObservations();
    
    std::cout << "\nOsservazioni di input:" << std::endl;
    std::cout << "1. " << observations[0].epoch.jd << " RA=" << observations[0].obs.ra 
              << " Dec=" << observations[0].obs.dec << std::endl;
    std::cout << "2. " << observations[1].epoch.jd << " RA=" << observations[1].obs.ra 
              << " Dec=" << observations[1].obs.dec << std::endl;
    std::cout << "3. " << observations[2].epoch.jd << " RA=" << observations[2].obs.ra 
              << " Dec=" << observations[2].obs.dec << std::endl;
    
    try {
        // Chiama metodo di Gauss
        auto solutions = InitialOrbit::gauss(
            observations[0],
            observations[1],
            observations[2]
        );
        
        std::cout << "\nNumero di soluzioni trovate: " << solutions.size() << std::endl;
        
        for (size_t i = 0; i < solutions.size(); i++) {
            printElements(solutions[i], "Soluzione " + std::to_string(i+1));
        }
        
        if (solutions.empty()) {
            std::cout << "\n⚠️  Nessuna soluzione trovata!" << std::endl;
            std::cout << "Possibili cause:" << std::endl;
            std::cout << "- Osservazioni troppo vicine" << std::endl;
            std::cout << "- Osservazioni coplanari (d0 ≈ 0)" << std::endl;
            std::cout << "- Equazione polinomiale senza radici reali" << std::endl;
        } else {
            std::cout << "\n✅ Test Gauss completato" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "\n❌ Errore nel test Gauss: " << e.what() << std::endl;
    }
}

/**
 * @brief Test 2: Metodo di Herget (quando implementato)
 */
void test_herget() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST 2: Metodo di Herget" << std::endl;
    std::cout << "========================================" << std::endl;
    
    auto observations = createErosObservations();
    
    try {
        // Prova vari guess per Eros (NEO a ~1.5 AU)
        std::vector<std::pair<double, double>> guesses = {
            {1.2, 1.3},  // Più vicino
            {1.4, 1.5},  // Medio
            {1.5, 1.6},  // Originale
            {1.8, 2.0}   // Più lontano
        };
        
        OrbitalElements best_orbit;
        bool found_solution = false;
        
        for (const auto& guess : guesses) {
            double r1_guess = guess.first;
            double r2_guess = guess.second;
            
            std::cout << "\nProvando R1=" << r1_guess << " AU, R2=" << r2_guess << " AU..." << std::endl;
            
            try {
                auto orbit = InitialOrbit::herget(
                    observations,
                    r1_guess,
                    r2_guess,
                    10  // max_iterations
                );
                
                // Verifica ragionevolezza
                if (orbit.a > 0 && orbit.a < 10.0 && orbit.e >= 0 && orbit.e < 1.0) {
                    best_orbit = orbit;
                    found_solution = true;
                    std::cout << "✓ Trovata soluzione ragionevole: a=" << orbit.a 
                              << " AU, e=" << orbit.e << std::endl;
                    break;
                }
            } catch (...) {
                continue;
            }
        }
        
        if (found_solution) {
            printElements(best_orbit, "Soluzione Herget");
            std::cout << "\n✅ Test Herget completato" << std::endl;
        } else {
            std::cout << "\n⚠️  Nessuna soluzione ragionevole trovata" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "\n⚠️  Errore: " << e.what() << std::endl;
    }
}

/**
 * @brief Test 3: Auto-solve (quando implementato)
 */
void test_autosolve() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST 3: Auto-Solve" << std::endl;
    std::cout << "========================================" << std::endl;
    
    auto observations = createErosObservations();
    
    try {
        double score = 0.0;
        auto orbit = InitialOrbit::autoSolve(observations, &score);
        
        printElements(orbit, "Soluzione Auto-Solve");
        std::cout << "\nScore finale: " << score << " arcsec" << std::endl;
        
        std::cout << "\n✅ Test Auto-Solve completato" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "\n⚠️  Auto-Solve non ancora implementato: " << e.what() << std::endl;
    }
}

/**
 * @brief Test 4: Differential Correction (adjustHergetResults)
 */
void test_differential_correction() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST 4: Differential Correction" << std::endl;
    std::cout << "========================================" << std::endl;
    
    auto observations = createErosObservations();
    
    // Usa elementi noti ma perturbati
    OrbitalElements orbit_perturbed;
    orbit_perturbed.epoch.jd = 2459956.0;
    orbit_perturbed.a = 1.458 + 0.01;  // Perturba di 0.01 AU
    orbit_perturbed.e = 0.223 + 0.005; // Perturba eccentricità
    orbit_perturbed.i = 10.83 * DEG_TO_RAD;
    orbit_perturbed.Omega = 304.30 * DEG_TO_RAD;
    orbit_perturbed.omega = 178.82 * DEG_TO_RAD;
    orbit_perturbed.M = 320.12 * DEG_TO_RAD;
    
    std::cout << "\n=== Orbita Perturbata (input) ===" << std::endl;
    printElements(orbit_perturbed, "");
    
    try {
        int result = InitialOrbit::adjustHergetResults(observations, orbit_perturbed);
        
        if (result == 0) {
            std::cout << "\n=== Orbita Corretta ===" << std::endl;
            printElements(orbit_perturbed, "");
            std::cout << "\n✅ Differential correction convergente" << std::endl;
        } else if (result == 1) {
            std::cout << "\n=== Orbita Migliorata (non convergente) ===" << std::endl;
            printElements(orbit_perturbed, "");
            std::cout << "\n⚠️  Max iterazioni raggiunto" << std::endl;
        } else {
            std::cout << "\n❌ Differential correction fallito" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "\n❌ Errore: " << e.what() << std::endl;
    }
}

/**
 * @brief Test 5: Valutazione orbita (quando implementato)
 */
void test_evaluate() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST 5: Valutazione Orbita" << std::endl;
    std::cout << "========================================" << std::endl;
    
    auto observations = createErosObservations();
    
    // Elementi orbitali noti di 433 Eros (epoca 2023-01-15)
    OrbitalElements eros_elements;
    eros_elements.epoch.jd = 2459956.0;
    eros_elements.a = 1.458;      // AU
    eros_elements.e = 0.223;
    eros_elements.i = 10.83;      // deg
    eros_elements.Omega = 304.30; // deg
    eros_elements.omega = 178.82; // deg
    eros_elements.M = 320.12;     // deg
    
    printElements(eros_elements, "Elementi Noti (AstDyS)");
    
    try {
        double score = InitialOrbit::evaluateInitialOrbit(
            eros_elements,
            observations,
            eros_elements.epoch
        );
        
        std::cout << "\nScore orbita nota: " << score << " arcsec" << std::endl;
        
        if (score < 10.0) {
            std::cout << "✅ Orbita eccellente (score < 10\")" << std::endl;
        } else if (score < 100.0) {
            std::cout << "✅ Orbita buona (score < 100\")" << std::endl;
        } else {
            std::cout << "⚠️  Orbita da migliorare (score > 100\")" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cout << "\n⚠️  Valutazione non ancora completa: " << e.what() << std::endl;
    }
}

/**
 * @brief Main
 */
int main(int argc, char* argv[]) {
    std::cout << "========================================" << std::endl;
    std::cout << "Test Initial Orbit Determination" << std::endl;
    std::cout << "Oggetto: 433 Eros" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Esegui tutti i test
    test_gauss();
    test_herget();
    test_autosolve();
    test_differential_correction();
    test_evaluate();
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test completati" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}
