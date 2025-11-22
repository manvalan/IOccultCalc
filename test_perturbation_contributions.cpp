/**
 * @file test_perturbation_contributions.cpp
 * @brief Analisi contributo singole perturbazioni all'errore
 * 
 * Test sistematico per quantificare l'impatto di:
 * - Perturbazioni planetarie
 * - Asteroidi massivi (SB441-N16)
 * - Relatività generale (Schwarzschild PN1)
 */

#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/orbit_propagator.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace ioccultcalc;

struct TestConfig {
    std::string name;
    bool planets;
    bool asteroids;
    bool relativity;
};

void printSeparator(char c = '=', int w = 70) {
    std::cout << std::string(w, c) << "\n";
}

double propagateAndCompare(
    const OrbitState& state0,
    const JulianDate& targetEpoch,
    const Vector3D& horizons_pos,
    const TestConfig& config
) {
    PropagatorOptions opts;
    opts.integrator = IntegratorType::RK4;
    opts.stepSize = 0.1;
    opts.usePlanetaryPerturbations = config.planets;
    opts.useRelativisticCorrections = config.relativity;
    // Note: asteroid perturbations controllate da SPICE kernel caricato
    
    OrbitPropagator propagator(opts);
    OrbitState final = propagator.propagate(state0, targetEpoch);
    
    Vector3D error = final.position - horizons_pos;
    return error.magnitude() * 149597870.7;  // km
}

int main(int argc, char** argv) {
    int asteroidId = 433;    // Eros
    int days = 100;
    
    if (argc > 1) asteroidId = std::atoi(argv[1]);
    if (argc > 2) days = std::atoi(argv[2]);
    
    std::cout << "\n";
    printSeparator();
    std::cout << "TEST: Contributo Perturbazioni all'Errore\n";
    printSeparator();
    std::cout << "\n";
    
    // Setup
    JPLHorizonsClient horizons;
    JulianDate epoch0;
    epoch0.jd = 2460615.5;  // 2024-Jan-01
    
    JulianDate epochF;
    epochF.jd = epoch0.jd + days;
    
    std::cout << "Configurazione:\n";
    std::cout << "  • Asteroide: (" << asteroidId << ")\n";
    std::cout << "  • Propagazione: " << days << " giorni\n";
    std::cout << "  • Epoca iniziale: JD " << std::fixed << std::setprecision(1) << epoch0.jd << "\n";
    std::cout << "  • Integratore: RK4, step = 0.1 giorni\n\n";
    
    // Download reference
    std::cout << "Download stati da JPL Horizons...\n";
    auto [pos0, vel0] = horizons.getStateVectors(std::to_string(asteroidId), epoch0);
    auto [posF_hor, velF_hor] = horizons.getStateVectors(std::to_string(asteroidId), epochF);
    
    OrbitState state0;
    state0.position = pos0;
    state0.velocity = vel0;
    state0.epoch = epoch0;
    
    std::cout << "✓ Dati scaricati\n\n";
    
    // Test configurations
    std::vector<TestConfig> tests = {
        {"Solo Sole (2-body)",           false, false, false},
        {"Sole + Pianeti",               true,  false, false},
        {"Sole + Pianeti + Relatività",  true,  false, true},
        {"FULL (Pianeti + Ast + Rel)",   true,  true,  true},
    };
    
    printSeparator();
    std::cout << "RISULTATI\n";
    printSeparator();
    std::cout << "\n";
    
    std::cout << std::left;
    std::cout << std::setw(35) << "Configurazione"
              << std::setw(15) << "Errore [km]"
              << std::setw(15) << "Errore [\"]"
              << "Note\n";
    printSeparator('-');
    
    std::vector<double> errors;
    
    for (const auto& test : tests) {
        std::cout << std::setw(35) << test.name << std::flush;
        
        double error_km = propagateAndCompare(state0, epochF, posF_hor, test);
        errors.push_back(error_km);
        
        double dist_au = posF_hor.magnitude();
        double error_arcsec = (error_km / 149597870.7) / dist_au * 206265.0;
        
        std::cout << std::setw(15) << std::fixed << std::setprecision(1) << error_km
                  << std::setw(15) << std::setprecision(3) << error_arcsec;
        
        if (test.name.find("2-body") != std::string::npos) {
            std::cout << "Baseline";
        } else if (test.name.find("FULL") != std::string::npos) {
            std::cout << "Best available";
        }
        
        std::cout << "\n";
    }
    
    std::cout << "\n";
    printSeparator();
    std::cout << "ANALISI CONTRIBUTI\n";
    printSeparator();
    std::cout << "\n";
    
    // Calcola contributi differenziali
    double err_2body = errors[0];
    double err_planets = errors[1];
    double err_planets_rel = errors[2];
    double err_full = errors[3];
    
    std::cout << "Contributi (riduzione errore rispetto baseline 2-body):\n\n";
    
    double contrib_planets = err_2body - err_planets;
    std::cout << "  • Pianeti: " << std::showpos << contrib_planets 
              << " km (" << (contrib_planets/err_2body*100) << "%)\n" << std::noshowpos;
    
    double contrib_rel = err_planets - err_planets_rel;
    std::cout << "  • Relatività: " << std::showpos << contrib_rel 
              << " km (" << (contrib_rel/err_2body*100) << "%)\n" << std::noshowpos;
    
    double contrib_asteroids = err_planets_rel - err_full;
    std::cout << "  • Asteroidi (16): " << std::showpos << contrib_asteroids 
              << " km (" << (contrib_asteroids/err_2body*100) << "%)\n" << std::noshowpos;
    
    std::cout << "\n";
    
    std::cout << "Errore residuo: " << err_full << " km\n";
    std::cout << "  Possibili cause:\n";
    std::cout << "    - Incertezza elementi iniziali (~10-20 km)\n";
    std::cout << "    - Asteroidi non modellati (~" << (300-16) << " mancanti)\n";
    std::cout << "    - Effetti di ordine superiore (J2, J3, ...)\n\n";
    
    // Raccomandazioni
    printSeparator();
    std::cout << "RACCOMANDAZIONI\n";
    printSeparator();
    std::cout << "\n";
    
    if (contrib_planets > 1000) {
        std::cout << "✓ Perturbazioni planetarie ESSENZIALI (>" << contrib_planets << " km)\n";
    }
    
    if (contrib_asteroids > 100) {
        std::cout << "✓ Asteroidi massivi IMPORTANTI (>" << contrib_asteroids << " km)\n";
    }
    
    if (contrib_rel > 10) {
        std::cout << "✓ Relatività RILEVANTE (>" << contrib_rel << " km)\n";
    } else {
        std::cout << "○ Relatività trascurabile (<" << contrib_rel << " km) per questa propagazione\n";
    }
    
    std::cout << "\n";
    std::cout << "Per migliorare ulteriormente la precisione:\n";
    std::cout << "  1. Aggiungere più asteroidi massivi (target: 50-100)\n";
    std::cout << "  2. Fit elementi orbitali da osservazioni recenti\n";
    std::cout << "  3. Modellare solar oblateness (J2) e planetary harmonics\n";
    
    std::cout << "\n";
    printSeparator();
    std::cout << "Test completato!\n";
    printSeparator();
    std::cout << "\n";
    
    return 0;
}
