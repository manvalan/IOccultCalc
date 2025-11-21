/**
 * @file advanced_propagation_example.cpp
 * @brief Esempio di propagazione orbitale con modello completo di forze
 * 
 * Dimostra:
 * - Configurazione ForceModel con vari livelli di precisione
 * - Propagazione con perturbazioni N-body complete
 * - Analisi contributi delle varie perturbazioni
 * - Confronto tra diversi modelli di forze
 */

#include <ioccultcalc/orbital_elements.h>
#include <ioccultcalc/force_model.h>
#include <ioccultcalc/numerical_integrator.h>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace ioccultcalc;

void printHeader(const std::string& title) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(70, '=') << "\n";
}

void example1_basicPropagation() {
    printHeader("Example 1: Basic N-body Propagation");
    
    // (1) Crea elementi orbitali di un asteroide (esempio: 472 Roma)
    EquinoctialElements elements;
    elements.a = 2.572;           // AU - semiasse maggiore
    elements.h = 0.0821;          // e*sin(ω+Ω)
    elements.k = 0.1234;          // e*cos(ω+Ω)  
    elements.p = 0.0512;          // tan(i/2)*sin(Ω)
    elements.q = 0.0378;          // tan(i/2)*cos(Ω)
    elements.lambda = 2.145;      // rad - longitudine media
    elements.epoch = JulianDate(2460000.0); // J2000
    
    std::cout << "\nInitial orbit (472 Roma-like):\n";
    std::cout << "  a = " << elements.a << " AU\n";
    std::cout << "  e = " << elements.getEccentricity() << "\n";
    std::cout << "  i = " << elements.getInclination() * 180/M_PI << " deg\n";
    
    // (2) Crea stato orbitale
    OrbitalState initialState = OrbitalState::fromElements(elements);
    
    // (3) Configura modello di forze (standard = tutti i pianeti)
    ForceModelConfig forceConfig = ForceModelConfig::standardConfig();
    
    std::cout << "\nForce model: STANDARD (all 8 planets + Moon)\n";
    
    // (4) Configura integratore
    IntegratorOptions integOptions;
    integOptions.relTolerance = 1e-12;
    integOptions.maxStep = 10.0; // giorni
    
    // (5) Propaga 1 anno
    JulianDate finalTime(elements.epoch.jd + 365.25);
    
    std::cout << "\nPropagating for 365 days...\n";
    auto start = std::chrono::high_resolution_clock::now();
    
    IntegrationResult result = IntegratorUtils::propagateWithForceModel(
        initialState, finalTime, forceConfig, integOptions);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    // (6) Risultati
    if (result.success) {
        std::cout << "\n✓ Propagation successful!\n";
        std::cout << "  Number of steps: " << result.numberOfSteps << "\n";
        std::cout << "  Rejected steps: " << result.numberOfRejectedSteps << "\n";
        std::cout << "  Energy error: " << std::scientific << result.energyError << "\n";
        std::cout << "  Computation time: " << duration.count() << " ms\n";
        
        // Posizione finale
        Vector3D finalPos = result.finalState.position;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "\nFinal position (AU):\n";
        std::cout << "  x = " << finalPos.x << "\n";
        std::cout << "  y = " << finalPos.y << "\n";
        std::cout << "  z = " << finalPos.z << "\n";
    } else {
        std::cout << "\n✗ Propagation failed: " << result.errorMessage << "\n";
    }
}

void example2_perturbationAnalysis() {
    printHeader("Example 2: Perturbation Analysis");
    
    // Stato di un asteroide main-belt a 2.5 AU
    Vector3D pos(2.5, 0.0, 0.1);  // AU
    Vector3D vel(0.0, 0.012, 0.0); // AU/day (~17 km/s tipico main-belt)
    
    double jd = 2460000.0; // Epoca
    
    std::cout << "\nAsteroid state:\n";
    std::cout << "  Position: (" << pos.x << ", " << pos.y << ", " << pos.z << ") AU\n";
    std::cout << "  Distance: " << pos.magnitude() << " AU\n";
    
    // Crea force model completo
    ForceModelConfig config = ForceModelConfig::highPrecisionConfig();
    ForceModel model(config);
    
    // Analizza contributi
    std::vector<ForceModelAnalyzer::PerturbationAnalysis> analysis =
        ForceModelAnalyzer::analyze(model, jd, pos, vel);
    
    std::cout << "\n" << ForceModelAnalyzer::generateReport(analysis);
    
    // Stima errori se omettiamo alcuni corpi
    std::cout << "\nEstimated position error if body is omitted (1 year propagation):\n";
    
    std::vector<PerturbingBody> bodiesToTest = {
        PerturbingBody::JUPITER,
        PerturbingBody::SATURN,
        PerturbingBody::EARTH,
        PerturbingBody::MARS,
        PerturbingBody::VENUS
    };
    
    for (auto body : bodiesToTest) {
        double error_km = ForceModelAnalyzer::estimateOmissionError(
            body, jd, pos, 365.25);
        
        std::string bodyName;
        switch(body) {
            case PerturbingBody::JUPITER: bodyName = "Jupiter"; break;
            case PerturbingBody::SATURN: bodyName = "Saturn"; break;
            case PerturbingBody::EARTH: bodyName = "Earth"; break;
            case PerturbingBody::MARS: bodyName = "Mars"; break;
            case PerturbingBody::VENUS: bodyName = "Venus"; break;
            default: bodyName = "Unknown"; break;
        }
        
        std::cout << "  " << std::setw(10) << bodyName << ": "
                  << std::setw(12) << std::fixed << std::setprecision(2)
                  << error_km << " km\n";
    }
}

void example3_configurationComparison() {
    printHeader("Example 3: Force Model Configuration Comparison");
    
    // Stato iniziale
    EquinoctialElements elements;
    elements.a = 2.572;
    elements.h = 0.0821;
    elements.k = 0.1234;
    elements.p = 0.0512;
    elements.q = 0.0378;
    elements.lambda = 2.145;
    elements.epoch = JulianDate(2460000.0);
    
    OrbitalState initialState = OrbitalState::fromElements(elements);
    JulianDate finalTime(elements.epoch.jd + 365.25);
    
    IntegratorOptions integOptions;
    integOptions.relTolerance = 1e-12;
    
    // Test diverse configurazioni
    struct TestCase {
        std::string name;
        ForceModelConfig config;
    };
    
    std::vector<TestCase> testCases = {
        {"FAST (giants only)", ForceModelConfig::fastConfig()},
        {"STANDARD (all planets)", ForceModelConfig::standardConfig()},
        {"HIGH (planets + asteroids)", ForceModelConfig::highPrecisionConfig()}
    };
    
    std::cout << "\nPropagating 365 days with different force models:\n\n";
    std::cout << std::setw(30) << "Configuration"
              << std::setw(15) << "Time (ms)"
              << std::setw(15) << "Steps"
              << std::setw(20) << "Final X (AU)\n";
    std::cout << std::string(80, '-') << "\n";
    
    for (const auto& test : testCases) {
        auto start = std::chrono::high_resolution_clock::now();
        
        IntegrationResult result = IntegratorUtils::propagateWithForceModel(
            initialState, finalTime, test.config, integOptions);
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        if (result.success) {
            std::cout << std::setw(30) << test.name
                      << std::setw(15) << duration.count()
                      << std::setw(15) << result.numberOfSteps
                      << std::setw(20) << std::fixed << std::setprecision(8)
                      << result.finalState.position.x << "\n";
        } else {
            std::cout << std::setw(30) << test.name
                      << " FAILED\n";
        }
    }
    
    std::cout << "\nNote: Differences in final position show impact of perturbations\n";
}

void example4_outputTrajectory() {
    printHeader("Example 4: Generate Trajectory with Output Points");
    
    // Stato iniziale
    EquinoctialElements elements;
    elements.a = 2.572;
    elements.h = 0.0821;
    elements.k = 0.1234;
    elements.p = 0.0512;
    elements.q = 0.0378;
    elements.lambda = 2.145;
    elements.epoch = JulianDate(2460000.0);
    
    OrbitalState initialState = OrbitalState::fromElements(elements);
    
    // Genera output ogni 10 giorni per 1 anno
    std::vector<JulianDate> outputTimes;
    for (int i = 0; i <= 36; ++i) {
        outputTimes.push_back(JulianDate(elements.epoch.jd + i * 10.0));
    }
    
    std::cout << "\nGenerating trajectory with " << outputTimes.size()
              << " output points (10-day intervals)\n";
    
    ForceModelConfig config = ForceModelConfig::standardConfig();
    IntegratorOptions options;
    options.relTolerance = 1e-12;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    std::vector<OrbitalState> trajectory =
        IntegratorUtils::propagateWithOutputAndForceModel(
            initialState, outputTimes, config, options);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "\n✓ Trajectory computed in " << duration.count() << " ms\n";
    std::cout << "\nSample points:\n";
    std::cout << std::setw(10) << "Day"
              << std::setw(15) << "X (AU)"
              << std::setw(15) << "Y (AU)"
              << std::setw(15) << "Z (AU)"
              << std::setw(15) << "Distance (AU)\n";
    std::cout << std::string(70, '-') << "\n";
    
    // Mostra solo alcuni punti
    for (size_t i = 0; i < trajectory.size(); i += 6) {
        const auto& state = trajectory[i];
        double day = (state.epoch.jd - initialState.epoch.jd);
        double dist = state.position.magnitude();
        
        std::cout << std::setw(10) << std::fixed << std::setprecision(0) << day
                  << std::setw(15) << std::setprecision(6) << state.position.x
                  << std::setw(15) << state.position.y
                  << std::setw(15) << state.position.z
                  << std::setw(15) << dist << "\n";
    }
}

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║   IOccultCalc - Advanced Orbital Propagation Examples           ║\n";
    std::cout << "║   N-body Force Model with Planetary Perturbations                ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";
    
    try {
        example1_basicPropagation();
        example2_perturbationAnalysis();
        example3_configurationComparison();
        example4_outputTrajectory();
        
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "All examples completed successfully!\n";
        std::cout << std::string(70, '=') << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
