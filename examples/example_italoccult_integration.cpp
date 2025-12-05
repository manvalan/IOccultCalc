/**
 * @file example_italoccult_integration.cpp
 * @brief Esempio di integrazione ITALOccultLibrary in IOccultCalc
 * @author IOccultCalc Integration Team
 * @date 1 Dicembre 2025
 * 
 * Questo esempio mostra come usare ITALOccultLibrary per propagare
 * asteroidi in IOccultCalc con precisione AstDyn/RKF78.
 * 
 * Compilazione:
 *   g++ -std=c++17 -O2 -o example_italoccult \
 *       example_italoccult_integration.cpp \
 *       src/integration/italoccult_integration.cpp \
 *       -I./include -I/usr/local/include -I/usr/local/include/eigen3 \
 *       -L/usr/local/lib -litaloccultlib -lastdyn
 * 
 * Esecuzione:
 *   ./example_italoccult data/17030.eq1 61007.0
 */

#include <ioccultcalc/integration/italoccult_integration.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace ioccultcalc;

// Helper per stampare un AsteroidState
void printAsteroidState(const AsteroidState& state) {
    std::cout << std::fixed << std::setprecision(9);
    
    std::cout << "\n=== Stato Asteroide ===" << std::endl;
    std::cout << "Nome:  " << state.name << std::endl;
    std::cout << "Epoca: MJD " << std::setprecision(6) << state.epoch_mjd_tdb << " TDB" << std::endl;
    
    std::cout << "\nPosizione ICRF (AU):" << std::endl;
    std::cout << "  X = " << state.x_icrf_au << std::endl;
    std::cout << "  Y = " << state.y_icrf_au << std::endl;
    std::cout << "  Z = " << state.z_icrf_au << std::endl;
    
    double r = std::sqrt(
        state.x_icrf_au * state.x_icrf_au +
        state.y_icrf_au * state.y_icrf_au +
        state.z_icrf_au * state.z_icrf_au
    );
    std::cout << "  R = " << r << " AU (" 
              << (r * 149.597870700) << " milioni km)" << std::endl;
    
    std::cout << "\nVelocità ICRF (AU/day):" << std::endl;
    std::cout << "  VX = " << state.vx_icrf_au_day << std::endl;
    std::cout << "  VY = " << state.vy_icrf_au_day << std::endl;
    std::cout << "  VZ = " << state.vz_icrf_au_day << std::endl;
    
    double v = std::sqrt(
        state.vx_icrf_au_day * state.vx_icrf_au_day +
        state.vy_icrf_au_day * state.vy_icrf_au_day +
        state.vz_icrf_au_day * state.vz_icrf_au_day
    );
    std::cout << "  V = " << v << " AU/day" << std::endl;
    
    std::cout << "\nParametri Orbitali (approssimati):" << std::endl;
    std::cout << "  a = " << std::setprecision(6) << state.semi_major_axis_au << " AU" << std::endl;
    std::cout << "  e = " << std::setprecision(9) << state.eccentricity << std::endl;
    std::cout << "  i = " << std::setprecision(6) << state.inclination_deg << " deg" << std::endl;
    
    if (!state.propagation_info.empty()) {
        std::cout << "\nInfo: " << state.propagation_info << std::endl;
    }
}

// Esempio 1: Propagazione singola con helper function
void example1_quick_propagate(const std::string& eq1_file, double target_mjd) {
    std::cout << "\n\n========================================" << std::endl;
    std::cout << "  ESEMPIO 1: Quick Propagate" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Propagazione one-liner con alta precisione
    auto state = quickPropagateFromEQ1(
        eq1_file,
        target_mjd,
        PropagationSettings::highAccuracy()
    );
    
    printAsteroidState(state);
}

// Esempio 2: Propagazione con oggetto integrator
void example2_integrator_object(const std::string& eq1_file, double target_mjd) {
    std::cout << "\n\n========================================" << std::endl;
    std::cout << "  ESEMPIO 2: Uso Integrator Object" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Crea integrator con configurazione personalizzata
    PropagationSettings settings;
    settings.tolerance = 1e-12;           // Alta precisione (1e-12 AU)
    settings.initial_step = 0.01;         // Step piccolo (0.01 giorni)
    settings.include_planets = true;      // Perturbazioni planetarie
    settings.include_relativity = true;   // Effetti relativistici
    
    ITALOccultIntegration integrator(settings);
    
    // Carica asteroide
    std::cout << "\nCaricamento " << eq1_file << "..." << std::endl;
    if (!integrator.loadAsteroidFromEQ1(eq1_file)) {
        std::cerr << "❌ Errore nel caricamento!" << std::endl;
        return;
    }
    std::cout << "✓ Asteroide caricato" << std::endl;
    
    // Propaga
    std::cout << "Propagazione a MJD " << target_mjd << "..." << std::endl;
    auto state = integrator.propagateToEpoch(target_mjd);
    
    printAsteroidState(state);
}

// Esempio 3: Propagazione multiple epoche (traiettoria)
void example3_multiple_epochs(const std::string& eq1_file) {
    std::cout << "\n\n========================================" << std::endl;
    std::cout << "  ESEMPIO 3: Traiettoria Multi-Epoca" << std::endl;
    std::cout << "========================================" << std::endl;
    
    ITALOccultIntegration integrator(PropagationSettings::highAccuracy());
    
    if (!integrator.loadAsteroidFromEQ1(eq1_file)) {
        std::cerr << "❌ Errore nel caricamento!" << std::endl;
        return;
    }
    
    // Genera 30 epoche (un'epoca al giorno per un mese)
    std::vector<double> epochs;
    double start_mjd = 61000.0;
    for (int i = 0; i < 30; ++i) {
        epochs.push_back(start_mjd + i);
    }
    
    std::cout << "\nPropagazione " << epochs.size() << " epoche..." << std::endl;
    auto states = integrator.propagateToEpochs(epochs);
    
    std::cout << "\nTraiettoria calcolata:" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n  MJD      |  R (AU)  | Delta R (km)" << std::endl;
    std::cout << "-----------|----------|---------------" << std::endl;
    
    double prev_r = 0.0;
    for (size_t i = 0; i < states.size(); ++i) {
        const auto& s = states[i];
        double r = std::sqrt(s.x_icrf_au * s.x_icrf_au +
                            s.y_icrf_au * s.y_icrf_au +
                            s.z_icrf_au * s.z_icrf_au);
        
        double delta_km = (i > 0) ? (r - prev_r) * 149597870.7 : 0.0;
        
        std::cout << " " << s.epoch_mjd_tdb << " | " 
                  << r << " | ";
        if (i > 0) {
            std::cout << std::setw(10) << std::setprecision(1) << delta_km;
        } else {
            std::cout << "     -    ";
        }
        std::cout << std::endl;
        
        prev_r = r;
    }
    
    // Statistiche
    double r_first = std::sqrt(
        states.front().x_icrf_au * states.front().x_icrf_au +
        states.front().y_icrf_au * states.front().y_icrf_au +
        states.front().z_icrf_au * states.front().z_icrf_au
    );
    double r_last = std::sqrt(
        states.back().x_icrf_au * states.back().x_icrf_au +
        states.back().y_icrf_au * states.back().y_icrf_au +
        states.back().z_icrf_au * states.back().z_icrf_au
    );
    
    std::cout << "\nStatistiche traiettoria:" << std::endl;
    std::cout << "  Intervallo: " << (epochs.back() - epochs.front()) << " giorni" << std::endl;
    std::cout << "  R iniziale: " << r_first << " AU" << std::endl;
    std::cout << "  R finale:   " << r_last << " AU" << std::endl;
    std::cout << "  Delta R:    " << std::setprecision(9) << (r_last - r_first) << " AU" << std::endl;
    std::cout << "              " << std::setprecision(1) 
              << ((r_last - r_first) * 149597870.7) << " km" << std::endl;
}

// Esempio 4: Confronto modalità Fast vs High Accuracy
void example4_compare_modes(const std::string& eq1_file, double target_mjd) {
    std::cout << "\n\n========================================" << std::endl;
    std::cout << "  ESEMPIO 4: Fast vs High Accuracy" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Fast mode
    auto state_fast = quickPropagateFromEQ1(
        eq1_file, target_mjd, PropagationSettings::fast()
    );
    
    // High accuracy mode
    auto state_accurate = quickPropagateFromEQ1(
        eq1_file, target_mjd, PropagationSettings::highAccuracy()
    );
    
    // Calcola differenze
    double dx = state_accurate.x_icrf_au - state_fast.x_icrf_au;
    double dy = state_accurate.y_icrf_au - state_fast.y_icrf_au;
    double dz = state_accurate.z_icrf_au - state_fast.z_icrf_au;
    double diff_au = std::sqrt(dx*dx + dy*dy + dz*dz);
    double diff_km = diff_au * 149597870.7;
    
    std::cout << "\nModalità FAST:" << std::endl;
    std::cout << "  Posizione: [" << std::setprecision(6)
              << state_fast.x_icrf_au << ", "
              << state_fast.y_icrf_au << ", "
              << state_fast.z_icrf_au << "] AU" << std::endl;
    
    std::cout << "\nModalità HIGH ACCURACY:" << std::endl;
    std::cout << "  Posizione: [" << std::setprecision(9)
              << state_accurate.x_icrf_au << ", "
              << state_accurate.y_icrf_au << ", "
              << state_accurate.z_icrf_au << "] AU" << std::endl;
    
    std::cout << "\nDifferenza:" << std::endl;
    std::cout << "  " << std::scientific << diff_au << " AU" << std::endl;
    std::cout << "  " << std::fixed << std::setprecision(3) 
              << diff_km << " km" << std::endl;
    
    if (diff_km < 1.0) {
        std::cout << "  ✓ Differenza < 1 km (entrambe le modalità adeguate)" << std::endl;
    } else if (diff_km < 100.0) {
        std::cout << "  ⚠ Differenza " << diff_km << " km (considera high accuracy)" << std::endl;
    } else {
        std::cout << "  ❌ Differenza significativa! Usa high accuracy" << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::cout << "========================================" << std::endl;
    std::cout << "  ITALOccultLibrary Integration Examples" << std::endl;
    std::cout << "  IOccultCalc + AstDyn/RKF78" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Parametri di default
    std::string eq1_file = "astdyn/data/17030.eq1";
    double target_mjd = 61007.0;
    
    // Override da command line
    if (argc >= 2) {
        eq1_file = argv[1];
    }
    if (argc >= 3) {
        target_mjd = std::stod(argv[2]);
    }
    
    std::cout << "\nParametri:" << std::endl;
    std::cout << "  File .eq1: " << eq1_file << std::endl;
    std::cout << "  Target MJD: " << target_mjd << std::endl;
    
    try {
        // Esegui tutti gli esempi
        example1_quick_propagate(eq1_file, target_mjd);
        example2_integrator_object(eq1_file, target_mjd);
        example3_multiple_epochs(eq1_file);
        example4_compare_modes(eq1_file, target_mjd);
        
        std::cout << "\n\n========================================" << std::endl;
        std::cout << "  ✓ Tutti gli esempi completati" << std::endl;
        std::cout << "========================================\n" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
