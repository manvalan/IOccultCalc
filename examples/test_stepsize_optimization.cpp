/**
 * @file test_stepsize_optimization.cpp
 * @brief Ottimizzazione step size per integratore RK4
 * 
 * Test sistematico per trovare il miglior trade-off tra:
 * - Precisione (errore vs JPL Horizons)
 * - Velocità (tempo CPU)
 * - Stabilità (accumulo errori)
 * 
 * Testa step size: 0.1, 0.05, 0.01, 0.005, 0.001 giorni
 */

#include <ioccultcalc/orbit_propagator.h>
#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/time_utils.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>

using namespace ioccultcalc;

struct StepSizeResult {
    double stepDays;
    int numSteps;
    double errorKm;
    double errorArcsec;
    double cpuSeconds;
    double errorPerStep;  // Errore medio per step (km)
};

void printSeparator(char c = '=', int len = 70) {
    std::cout << std::string(len, c) << "\n";
}

/**
 * @brief Propaga e confronta con Horizons
 */
StepSizeResult testStepSize(
    const OrbitState& state0,
    double propagationDays,
    double stepDays,
    const Vector3D& horizonsPos
) {
    StepSizeResult result;
    result.stepDays = stepDays;
    result.numSteps = static_cast<int>(propagationDays / stepDays);
    
    // Configura propagatore
    PropagatorOptions opts;
    opts.integrator = IntegratorType::RK4;
    opts.stepSize = stepDays;
    opts.usePlanetaryPerturbations = true;
    opts.useRelativisticCorrections = true;
    
    OrbitPropagator prop(opts);
    
    // Misura tempo esecuzione
    auto start = std::chrono::high_resolution_clock::now();
    
    JulianDate targetEpoch;
    targetEpoch.jd = state0.epoch.jd + propagationDays;
    OrbitState stateFinal = prop.propagate(state0, targetEpoch);
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    result.cpuSeconds = elapsed.count();
    
    // Calcola errore
    Vector3D diff = stateFinal.position - horizonsPos;
    result.errorKm = diff.magnitude() * 149597870.7;  // AU → km
    
    // Errore angolare @ 1.5 AU (tipico main belt)
    double distance = horizonsPos.magnitude();
    result.errorArcsec = (result.errorKm / (distance * 149597870.7)) * 206265.0;
    
    // Errore medio per step
    result.errorPerStep = result.errorKm / result.numSteps;
    
    return result;
}

int main(int argc, char* argv[]) {
    std::cout << "\n";
    printSeparator();
    std::cout << "TEST: OTTIMIZZAZIONE STEP SIZE RK4\n";
    printSeparator();
    
    // Parametri test
    int asteroidId = (argc > 1) ? std::stoi(argv[1]) : 433;  // Default: Eros
    double propagationDays = (argc > 2) ? std::stod(argv[2]) : 100.0;
    
    std::cout << "\nConfigurazione:\n";
    std::cout << "  • Asteroide: (" << asteroidId << ")\n";
    std::cout << "  • Propagazione: " << propagationDays << " giorni\n";
    std::cout << "  • Epoca iniziale: 2024-Jan-01 (JD 2460615.5)\n";
    std::cout << "  • Perturbazioni: Pianeti + Asteroidi + Relatività\n\n";
    
    // Step sizes da testare
    std::vector<double> stepSizes = {0.1, 0.05, 0.01, 0.005, 0.001};
    
    // ================================================================
    // STEP 1: Scarica stato iniziale e finale da Horizons
    // ================================================================
    printSeparator('-');
    std::cout << "STEP 1: Download dati JPL Horizons\n";
    printSeparator('-');
    
    JPLHorizonsClient horizons;
    
    JulianDate epoch0_utc;
    epoch0_utc.jd = 2460615.5;  // 2024-Jan-01 00:00 UTC
    
    // NOTA: Horizons interpreta JD input come TDB automaticamente
    // Quindi quando gli passiamo 2460615.5, lui lo legge come TDB
    // Per avere coerenza, convertiamo epoch a TDB anche noi
    JulianDate epoch0_tdb = TimeUtils::utcToTDB(epoch0_utc);
    
    std::cout << "Epoch:\n";
    std::cout << "  • Input UTC: JD " << std::fixed << std::setprecision(6) << epoch0_utc.jd << "\n";
    std::cout << "  • TDB (per DE441): JD " << epoch0_tdb.jd << "\n";
    std::cout << "  • Differenza: " << (epoch0_tdb.jd - epoch0_utc.jd) * 86400.0 << " sec\n\n";
    
    std::cout << "Download stato iniziale (Horizons con TDB)...\n";
    auto [posInit, velInit] = horizons.getStateVectors(std::to_string(asteroidId), epoch0_tdb);
    
    OrbitState stateInitial;
    stateInitial.position = posInit;
    stateInitial.velocity = velInit;
    stateInitial.epoch = epoch0_tdb;  // Usa TDB per propagazione!
    
    std::cout << "✓ Stato iniziale scaricato\n";
    std::cout << "  r = (" << stateInitial.position.x << ", " 
              << stateInitial.position.y << ", " 
              << stateInitial.position.z << ") AU\n";
    std::cout << "  |r| = " << stateInitial.position.magnitude() << " AU\n\n";
    
    JulianDate epochFinal_tdb;
    epochFinal_tdb.jd = epoch0_tdb.jd + propagationDays;
    
    std::cout << "Download stato finale (TDB)...\n";
    auto [posFinal, velFinal] = horizons.getStateVectors(std::to_string(asteroidId), epochFinal_tdb);
    
    std::cout << "✓ Stato finale scaricato (riferimento JPL)\n";
    std::cout << "  r = (" << posFinal.x << ", " 
              << posFinal.y << ", " 
              << posFinal.z << ") AU\n";
    std::cout << "  |r| = " << posFinal.magnitude() << " AU\n\n";
    
    // ================================================================
    // STEP 2: Test tutti gli step sizes
    // ================================================================
    printSeparator('-');
    std::cout << "STEP 2: Test diversi step sizes\n";
    printSeparator('-');
    std::cout << "\n";
    
    std::vector<StepSizeResult> results;
    
    for (double step : stepSizes) {
        std::cout << "Test step = " << step << " giorni";
        std::cout << " (" << static_cast<int>(propagationDays / step) << " steps)... ";
        std::cout.flush();
        
        StepSizeResult res = testStepSize(
            stateInitial, 
            propagationDays, 
            step, 
            posFinal
        );
        
        results.push_back(res);
        
        std::cout << "✓ " << std::fixed << std::setprecision(2) 
                  << res.cpuSeconds << "s, " 
                  << res.errorKm << " km\n";
    }
    
    // ================================================================
    // STEP 3: Analisi risultati
    // ================================================================
    std::cout << "\n";
    printSeparator();
    std::cout << "RISULTATI\n";
    printSeparator();
    
    // Tabella comparativa
    std::cout << "\n";
    std::cout << std::left;
    std::cout << std::setw(12) << "Step [day]"
              << std::setw(12) << "N Steps"
              << std::setw(14) << "CPU [s]"
              << std::setw(16) << "Error [km]"
              << std::setw(16) << "Error [\"]"
              << std::setw(18) << "Error/step [m]\n";
    printSeparator('-');
    
    for (const auto& res : results) {
        std::cout << std::fixed
                  << std::setw(12) << std::setprecision(3) << res.stepDays
                  << std::setw(12) << res.numSteps
                  << std::setw(14) << std::setprecision(3) << res.cpuSeconds
                  << std::setw(16) << std::setprecision(1) << res.errorKm
                  << std::setw(16) << std::setprecision(4) << res.errorArcsec
                  << std::setw(18) << std::setprecision(1) << (res.errorPerStep * 1000.0)
                  << "\n";
    }
    
    // ================================================================
    // STEP 4: Raccomandazioni
    // ================================================================
    std::cout << "\n";
    printSeparator('-');
    std::cout << "ANALISI & RACCOMANDAZIONI\n";
    printSeparator('-');
    
    // Trova miglior trade-off (< 100 km errore, tempo ragionevole)
    StepSizeResult* best = nullptr;
    for (auto& res : results) {
        if (res.errorKm < 100.0) {
            if (best == nullptr || res.cpuSeconds < best->cpuSeconds) {
                best = &res;
            }
        }
    }
    
    // Trova step più veloce
    auto fastest = std::min_element(results.begin(), results.end(),
        [](const StepSizeResult& a, const StepSizeResult& b) {
            return a.cpuSeconds < b.cpuSeconds;
        });
    
    // Trova step più preciso
    auto mostPrecise = std::min_element(results.begin(), results.end(),
        [](const StepSizeResult& a, const StepSizeResult& b) {
            return a.errorKm < b.errorKm;
        });
    
    std::cout << "\nPiù veloce:\n";
    std::cout << "  • Step: " << fastest->stepDays << " giorni\n";
    std::cout << "  • Tempo: " << fastest->cpuSeconds << " s\n";
    std::cout << "  • Errore: " << fastest->errorKm << " km (" 
              << fastest->errorArcsec << "\")\n\n";
    
    std::cout << "Più preciso:\n";
    std::cout << "  • Step: " << mostPrecise->stepDays << " giorni\n";
    std::cout << "  • Tempo: " << mostPrecise->cpuSeconds << " s\n";
    std::cout << "  • Errore: " << mostPrecise->errorKm << " km (" 
              << mostPrecise->errorArcsec << "\")\n\n";
    
    if (best != nullptr) {
        std::cout << "✓ RACCOMANDATO (miglior trade-off):\n";
        std::cout << "  • Step: " << best->stepDays << " giorni\n";
        std::cout << "  • Tempo: " << best->cpuSeconds << " s\n";
        std::cout << "  • Errore: " << best->errorKm << " km (" 
                  << best->errorArcsec << "\")\n";
        std::cout << "  • Precisione: <100 km con tempo ottimale\n\n";
    } else {
        std::cout << "⚠ ATTENZIONE: Nessuno step raggiunge <100 km errore\n";
        std::cout << "  → Considerare step più piccoli o propagazione più breve\n\n";
    }
    
    // Analisi scaling
    std::cout << "Scaling error vs step size:\n";
    if (results.size() >= 2) {
        double ratio = results[0].errorKm / results.back().errorKm;
        double stepRatio = results[0].stepDays / results.back().stepDays;
        double power = std::log(ratio) / std::log(stepRatio);
        std::cout << "  • Errore ~ step^" << std::fixed << std::setprecision(2) 
                  << power << "\n";
        if (power >= 3.5 && power <= 4.5) {
            std::cout << "  • ✓ Comportamento RK4 atteso (O(h⁴))\n";
        } else if (power < 3.5) {
            std::cout << "  • ⚠ Errore cresce più lentamente del previsto\n";
            std::cout << "    → Dominano perturbazioni o accumulo numerico\n";
        } else {
            std::cout << "  • ⚠ Errore cresce più rapidamente del previsto\n";
            std::cout << "    → Possibile instabilità numerica\n";
        }
    }
    
    std::cout << "\n";
    printSeparator();
    std::cout << "Test completato!\n";
    printSeparator();
    std::cout << "\n";
    
    return 0;
}
