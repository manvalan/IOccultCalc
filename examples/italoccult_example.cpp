/**
 * @file examples/italoccult_example.cpp
 * @brief Example di utilizzo di ITALOccultLibrary in IOccultCalc
 * @date 1 Dicembre 2025
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "ioccultcalc/integration/italoccult_integration.h"

void printAsteroidState(const AsteroidState& state) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Asteroide: " << state.name << std::endl;
    std::cout << "Epoca: MJD " << std::fixed << std::setprecision(5) 
              << state.epoch_mjd_tdb << " TDB" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    // Posizione
    double distance = std::hypot({state.pos_x, state.pos_y, state.pos_z});
    std::cout << "\nPosizione ICRF (AU):" << std::endl;
    std::cout << "  X = " << std::setprecision(9) << state.pos_x << std::endl;
    std::cout << "  Y = " << std::setprecision(9) << state.pos_y << std::endl;
    std::cout << "  Z = " << std::setprecision(9) << state.pos_z << std::endl;
    std::cout << "  R = " << std::setprecision(6) << distance 
              << " AU (" << std::setprecision(1) 
              << (distance * 149597870.7) << " milioni km)" << std::endl;
    
    // Velocità
    double speed = std::hypot({state.vel_x, state.vel_y, state.vel_z});
    std::cout << "\nVelocità ICRF (AU/day):" << std::endl;
    std::cout << "  VX = " << std::setprecision(9) << state.vel_x << std::endl;
    std::cout << "  VY = " << std::setprecision(9) << state.vel_y << std::endl;
    std::cout << "  VZ = " << std::setprecision(9) << state.vel_z << std::endl;
    std::cout << "  V = " << std::setprecision(6) << speed << " AU/day" << std::endl;
    
    // Parametri orbitali
    std::cout << "\nParametri orbitali (approssimati):" << std::endl;
    std::cout << "  Semi-major axis: " << std::setprecision(6) 
              << state.semi_major_axis << " AU" << std::endl;
    std::cout << "  Eccentricity: " << std::setprecision(6) 
              << state.eccentricity << std::endl;
    std::cout << "  Inclination: " << std::setprecision(6) 
              << state.inclination << "°" << std::endl;
    
    if (!state.prop_stats.empty()) {
        std::cout << "\nStatistiche propagazione:" << std::endl;
        std::cout << state.prop_stats << std::endl;
    }
    
    std::cout << std::string(60, '=') << std::endl;
}

int main() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "  ITALOccultLibrary - Esempio di utilizzo" << std::endl;
    std::cout << "  Integrazione con IOccultCalc" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    try {
        // =====================================================
        // ESEMPIO 1: Propagazione singola
        // =====================================================
        std::cout << "\n[ESEMPIO 1] Propagazione singola epoca" << std::endl;
        std::cout << "--------------------------------------------\n";
        
        ITALOccultIntegration integrator(PropagationSettings::highAccuracy());
        
        // Carica asteroide 17030 (Sierks) dal file .eq1
        std::cout << "Caricamento asteroide da file .eq1..." << std::endl;
        if (!integrator.loadAsteroidFromEQ1("data/17030.eq1")) {
            std::cerr << "Errore: Non riesco a caricare il file .eq1" << std::endl;
            return 1;
        }
        std::cout << "✓ Asteroide caricato\n" << std::endl;
        
        // Propaga a MJD 61007 (7 giorni dopo l'epoca iniziale 61000)
        std::cout << "Propagazione a MJD 61007..." << std::endl;
        auto state = integrator.propagateToEpoch(61007.0);
        printAsteroidState(state);
        
        // =====================================================
        // ESEMPIO 2: Propagazione multipla
        // =====================================================
        std::cout << "\n[ESEMPIO 2] Propagazione a multiple epoche" << std::endl;
        std::cout << "--------------------------------------------\n";
        
        std::vector<double> epochs = {
            61000.0,   // Epoca iniziale
            61001.0,   // +1 giorno
            61003.0,   // +3 giorni
            61007.0,   // +7 giorni
            61014.0    // +14 giorni
        };
        
        std::cout << "Propagazione a " << epochs.size() << " epoche..." << std::endl;
        auto states = integrator.propagateToEpochs(epochs);
        
        std::cout << "\nTraiettoria dell'asteroide 17030 (Sierks):\n" << std::endl;
        std::cout << std::left << std::setw(12) << "Epoca (MJD)"
                  << std::right << std::setw(15) << "Distanza (AU)"
                  << std::setw(15) << "Lat eclittica"
                  << std::setw(15) << "Lon eclittica" << std::endl;
        std::cout << std::string(57, '-') << std::endl;
        
        for (const auto& s : states) {
            double distance = std::hypot({s.pos_x, s.pos_y, s.pos_z});
            double eclip_lat = std::atan2(s.pos_z, 
                                std::hypot(s.pos_x, s.pos_y)) * 180.0 / M_PI;
            double eclip_lon = std::atan2(s.pos_y, s.pos_x) * 180.0 / M_PI;
            if (eclip_lon < 0) eclip_lon += 360.0;
            
            std::cout << std::fixed << std::setprecision(2) 
                      << std::left << std::setw(12) << s.epoch_mjd_tdb
                      << std::right << std::setw(15) << distance
                      << std::setw(15) << eclip_lat
                      << std::setw(15) << eclip_lon << std::endl;
        }
        
        // =====================================================
        // ESEMPIO 3: Utilizzo funzione helper
        // =====================================================
        std::cout << "\n[ESEMPIO 3] Utilizzo funzione helper quickPropagateFromEQ1" << std::endl;
        std::cout << "--------------------------------------------\n";
        
        std::cout << "One-liner per propagazione veloce..." << std::endl;
        auto quick_state = quickPropagateFromEQ1(
            "data/17030.eq1",
            61007.0,
            PropagationSettings::highAccuracy()
        );
        
        std::cout << "✓ Propagazione completata" << std::endl;
        std::cout << "Asteroide: " << quick_state.name << std::endl;
        std::cout << "Epoca: MJD " << std::fixed << std::setprecision(2) 
                  << quick_state.epoch_mjd_tdb << std::endl;
        
        double dist = std::hypot({quick_state.pos_x, quick_state.pos_y, quick_state.pos_z});
        std::cout << "Distanza: " << std::setprecision(4) << dist << " AU" << std::endl;
        
        // =====================================================
        // ESEMPIO 4: Modalità Fast per performance
        // =====================================================
        std::cout << "\n[ESEMPIO 4] Propagazione in modalità FAST" << std::endl;
        std::cout << "--------------------------------------------\n";
        
        ITALOccultIntegration fast_integrator(PropagationSettings::fast());
        fast_integrator.loadAsteroidFromEQ1("data/17030.eq1");
        
        std::cout << "Propagazione a 100 epoche con modalità fast..." << std::endl;
        std::vector<double> many_epochs;
        for (int i = 0; i < 100; ++i) {
            many_epochs.push_back(61000.0 + i * 0.1);  // Ogni 2.4 ore
        }
        
        auto many_states = fast_integrator.propagateToEpochs(many_epochs);
        
        std::cout << "✓ Propagate " << many_states.size() << " epoche" << std::endl;
        
        // Statistiche
        double min_dist = 1e9, max_dist = 0;
        for (const auto& s : many_states) {
            double d = std::hypot({s.pos_x, s.pos_y, s.pos_z});
            min_dist = std::min(min_dist, d);
            max_dist = std::max(max_dist, d);
        }
        
        std::cout << "Distanza minima: " << std::setprecision(6) 
                  << min_dist << " AU" << std::endl;
        std::cout << "Distanza massima: " << std::setprecision(6) 
                  << max_dist << " AU" << std::endl;
        std::cout << "Escursione: " << std::setprecision(6) 
                  << (max_dist - min_dist) << " AU" << std::endl;
        
        // =====================================================
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "✓ Tutti gli esempi completati con successo!" << std::endl;
        std::cout << std::string(60, '=') << std::endl << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ Errore: " << e.what() << std::endl;
        return 1;
    }
}
