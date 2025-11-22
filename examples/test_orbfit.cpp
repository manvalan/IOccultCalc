/**
 * @file test_orbfit.cpp
 * @brief Test integrazione OrbFit wrapper
 */

#include "ioccultcalc/orbfit_wrapper.h"
#include "ioccultcalc/mpc_client.h"
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace ioccultcalc;

void printElements(const OrbitalElements& elem, const std::string& title) {
    std::cout << "\n" << title << "\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Epoch (JD): " << elem.epoch.jd << "\n";
    std::cout << "  a (AU):     " << elem.a << "\n";
    std::cout << "  e:          " << elem.e << "\n";
    std::cout << "  i (deg):    " << elem.i << "\n";
    std::cout << "  Ω (deg):    " << elem.Omega << "\n";
    std::cout << "  ω (deg):    " << elem.omega << "\n";
    std::cout << "  M (deg):    " << elem.M << "\n";
    std::cout << std::string(60, '=') << "\n";
}

int main() {
    std::cout << "\n╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║     TEST ORBFIT WRAPPER - GAUSS METHOD                 ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n\n";

    // Initialize OrbFit
    std::cout << "Inizializzazione OrbFit...\n";
    if (!OrbFitWrapper::initialize("/Users/michelebigi/Astro/OrbFit/lib")) {
        std::cerr << "✗ Errore inizializzazione OrbFit\n";
        return 1;
    }
    std::cout << "✓ OrbFit inizializzato\n\n";

    // Carica osservazioni Eros da file RWO (già testato)
    std::cout << "Caricamento osservazioni 433 Eros...\n";
    
    // Usa file cache se esiste
    std::string cache_file = "/tmp/eros_obs.rwo";
    std::vector<AstrometricObservation> observations;
    
    if (std::ifstream(cache_file).good()) {
        // Carica da cache
        std::cout << "Caricamento da cache...\n";
        // TODO: implementare parsing RWO
        // Per ora usiamo MPC client
    }
    
    // Fallback: scarica da MPC
    MPCClient mpc_client;
    auto obs_set = mpc_client.getObservations("433");
    observations = obs_set.observations;
    
    if (observations.empty()) {
        std::cerr << "✗ Nessuna osservazione scaricata\n";
        return 1;
    }
    
    std::cout << "✓ Caricate " << observations.size() << " osservazioni\n\n";

    // Filtra osservazioni per span di 90 giorni
    double jd_start = observations.front().epoch.jd;
    double jd_end = jd_start + 90.0;
    
    std::vector<AstrometricObservation> filtered;
    for (const auto& obs : observations) {
        if (obs.epoch.jd >= jd_start && obs.epoch.jd <= jd_end) {
            filtered.push_back(obs);
        }
    }
    
    if (filtered.size() < 3) {
        std::cerr << "✗ Troppo poche osservazioni in 90 giorni\n";
        return 1;
    }
    
    // Seleziona 3 osservazioni: prima, metà, ultima
    AstrometricObservation obs1 = filtered.front();
    AstrometricObservation obs2 = filtered[filtered.size() / 2];
    AstrometricObservation obs3 = filtered.back();
    
    std::cout << "Osservazioni selezionate (span 90 giorni):\n";
    std::cout << "  1. JD " << std::fixed << std::setprecision(2) << obs1.epoch.jd 
              << " (RA=" << std::setprecision(4) << obs1.obs.ra 
              << "°, Dec=" << obs1.obs.dec << "°)\n";
    std::cout << "  2. JD " << obs2.epoch.jd 
              << " (RA=" << obs2.obs.ra << "°, Dec=" << obs2.obs.dec << "°)\n";
    std::cout << "  3. JD " << obs3.epoch.jd 
              << " (RA=" << obs3.obs.ra << "°, Dec=" << obs3.obs.dec << "°)\n\n";

    // Test metodo di Gauss di OrbFit
    std::cout << "Esecuzione metodo di Gauss (OrbFit)...\n";
    
    auto solutions = OrbFitWrapper::gauss(obs1, obs2, obs3, true);
    
    if (solutions.empty()) {
        std::cerr << "✗ Nessuna soluzione trovata\n";
        return 1;
    }
    
    std::cout << "\n✓ Trovate " << solutions.size() << " soluzioni\n";
    
    // Stampa tutte le soluzioni
    for (size_t i = 0; i < solutions.size(); i++) {
        printElements(solutions[i].elements, 
                     "Soluzione OrbFit #" + std::to_string(i+1) + 
                     " (root " + std::to_string(solutions[i].root_index) + ")");
        
        std::cout << "  Distanza topocentric: " << std::setprecision(4) 
                  << solutions[i].topocentric_distance << " AU\n";
        
        // Check fisica
        bool physical = (solutions[i].elements.a > 0.5 && 
                        solutions[i].elements.a < 5.0 &&
                        solutions[i].elements.e >= 0 && 
                        solutions[i].elements.e < 1.0);
        
        if (physical) {
            std::cout << "  ✓ Orbita fisica\n";
        } else {
            std::cout << "  ✗ Orbita non fisica\n";
        }
    }
    
    // Confronto con Eros reale
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "CONFRONTO CON EROS REALE\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "Eros reale:\n";
    std::cout << "  a = 1.458 AU\n";
    std::cout << "  e = 0.223\n";
    std::cout << "  i = 10.83°\n\n";
    
    // Trova soluzione migliore
    size_t best_idx = 0;
    double best_score = 1e9;
    
    for (size_t i = 0; i < solutions.size(); i++) {
        if (!solutions[i].is_valid) continue;
        
        double err_a = std::abs(solutions[i].elements.a - 1.458) / 1.458;
        double err_e = std::abs(solutions[i].elements.e - 0.223) / 0.223;
        double err_i = std::abs(solutions[i].elements.i - 10.83) / 10.83;
        double score = err_a + err_e + err_i;
        
        if (score < best_score) {
            best_score = score;
            best_idx = i;
        }
    }
    
    std::cout << "Soluzione migliore: #" << (best_idx + 1) << "\n";
    std::cout << "  Errore a: " << std::setprecision(1) 
              << std::abs(solutions[best_idx].elements.a - 1.458) / 1.458 * 100 << "%\n";
    std::cout << "  Errore e: " 
              << std::abs(solutions[best_idx].elements.e - 0.223) / 0.223 * 100 << "%\n";
    std::cout << "  Errore i: " 
              << std::abs(solutions[best_idx].elements.i - 10.83) / 10.83 * 100 << "%\n";
    
    // Cleanup
    OrbFitWrapper::cleanup();
    
    std::cout << "\n✓ Test completato\n\n";
    
    return 0;
}
