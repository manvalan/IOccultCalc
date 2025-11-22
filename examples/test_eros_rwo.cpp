/**
 * @file test_eros_rwo.cpp
 * @brief Test Initial Orbit Determination con dati reali di 433 Eros da RWO
 */

#include "ioccultcalc/initial_orbit.h"
#include "ioccultcalc/mpc_client.h"
#include "ioccultcalc/coordinates.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace ioccultcalc;

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
 * @brief Seleziona subset di osservazioni ben distribuite
 */
std::vector<AstrometricObservation> selectObservations(
    const std::vector<AstrometricObservation>& all_obs,
    int n_desired,
    double jd_start,
    double jd_end)
{
    std::vector<AstrometricObservation> filtered;
    
    // Filtra per intervallo temporale
    for (const auto& obs : all_obs) {
        if (obs.epoch.jd >= jd_start && obs.epoch.jd <= jd_end) {
            filtered.push_back(obs);
        }
    }
    
    if (filtered.size() <= n_desired) {
        return filtered;
    }
    
    // Seleziona uniformemente
    std::vector<AstrometricObservation> selected;
    double step = (double)filtered.size() / n_desired;
    
    for (int i = 0; i < n_desired; i++) {
        int idx = (int)(i * step);
        if (idx < filtered.size()) {
            selected.push_back(filtered[idx]);
        }
    }
    
    return selected;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Initial Orbit Determination Test" << std::endl;
    std::cout << "Dati reali: 433 Eros (RWO)" << std::endl;
    std::cout << "========================================" << std::endl;
    
    try {
        // Scarica osservazioni reali
        std::cout << "\n1. Scaricamento osservazioni da AstDyS..." << std::endl;
        MPCClient client;
        client.setTimeout(30);
        
        ObservationSet obsSet = client.getObservations("433");
        std::cout << "   ✓ Scaricate " << obsSet.observations.size() << " osservazioni" << std::endl;
        std::cout << "   Arc: " << std::fixed << std::setprecision(1) 
                  << obsSet.arcLength << " giorni" << std::endl;
        
        if (obsSet.observations.empty()) {
            std::cerr << "Errore: Nessuna osservazione disponibile!" << std::endl;
            return 1;
        }
        
        // Seleziona osservazioni su span più lungo (180 giorni per main belt)
        double jd_last = obsSet.lastObservation.jd;
        double jd_start = jd_last - 180.0;
        
        auto recent_obs = selectObservations(obsSet.observations, 20, jd_start, jd_last);
        
        std::cout << "\n2. Osservazioni selezionate: " << recent_obs.size() << std::endl;
        std::cout << "   Periodo: JD " << std::fixed << std::setprecision(2) 
                  << recent_obs.front().epoch.jd << " - " 
                  << recent_obs.back().epoch.jd << std::endl;
        std::cout << "   Span: " << (recent_obs.back().epoch.jd - recent_obs.front().epoch.jd) 
                  << " giorni" << std::endl;
        
        // Mostra prime osservazioni
        std::cout << "\n   Prima osservazione:" << std::endl;
        std::cout << "     JD: " << recent_obs[0].epoch.jd << std::endl;
        std::cout << "     RA: " << std::setprecision(6) << recent_obs[0].obs.ra << "°" << std::endl;
        std::cout << "     Dec: " << recent_obs[0].obs.dec << "°" << std::endl;
        std::cout << "     Obs: " << recent_obs[0].observatoryCode << std::endl;
        
        // TEST 1: Gauss con prime 3 osservazioni
        std::cout << "\n========================================" << std::endl;
        std::cout << "TEST 1: Metodo di Gauss" << std::endl;
        std::cout << "========================================" << std::endl;
        
        std::vector<AstrometricObservation> gauss_obs = {
            recent_obs[0], 
            recent_obs[recent_obs.size()/2], 
            recent_obs.back()
        };
        
        std::cout << "\nOsservazioni (prima, metà, ultima):" << std::endl;
        for (size_t i = 0; i < gauss_obs.size(); i++) {
            std::cout << "  " << (i+1) << ". JD " << std::fixed << std::setprecision(2) 
                      << gauss_obs[i].epoch.jd 
                      << " RA=" << std::setprecision(4) << gauss_obs[i].obs.ra 
                      << "° Dec=" << gauss_obs[i].obs.dec << "°" << std::endl;
        }
        
        try {
            auto gauss_solutions = InitialOrbit::gauss(gauss_obs[0], gauss_obs[1], gauss_obs[2]);
            
            std::cout << "\n✓ Gauss trovate " << gauss_solutions.size() << " soluzioni" << std::endl;
            
            for (size_t i = 0; i < gauss_solutions.size(); i++) {
                printElements(gauss_solutions[i], "Soluzione Gauss #" + std::to_string(i+1));
                
                // Check ragionevolezza
                bool reasonable = (gauss_solutions[i].a > 0 && 
                                  gauss_solutions[i].a < 10.0 &&
                                  gauss_solutions[i].e >= 0 && 
                                  gauss_solutions[i].e < 1.0);
                
                if (reasonable) {
                    std::cout << "✓ Orbita ragionevole" << std::endl;
                } else {
                    std::cout << "⚠ Orbita non fisica" << std::endl;
                }
            }
            
        } catch (const std::exception& e) {
            std::cout << "⚠ Gauss fallito: " << e.what() << std::endl;
        }
        
        // TEST 2: Herget con vari guess
        std::cout << "\n========================================" << std::endl;
        std::cout << "TEST 2: Metodo di Herget" << std::endl;
        std::cout << "========================================" << std::endl;
        
        std::cout << "\nUsando tutte " << recent_obs.size() << " osservazioni" << std::endl;
        
        // Prova vari guess di distanza
        std::vector<std::pair<double, double>> guesses = {
            {1.2, 1.3},  // Vicino
            {1.4, 1.5},  // Medio
            {1.6, 1.7},  // Lontano
            {1.8, 2.0}   // Molto lontano
        };
        
        OrbitalElements best_herget;
        bool found_herget = false;
        
        for (const auto& guess : guesses) {
            std::cout << "\nProvando R1=" << std::fixed << std::setprecision(2) 
                      << guess.first << " AU, R2=" << guess.second << " AU..." << std::endl;
            
            try {
                auto orbit = InitialOrbit::herget(recent_obs, guess.first, guess.second, 5);
                
                // Verifica ragionevolezza
                if (orbit.a > 0.5 && orbit.a < 3.0 && orbit.e >= 0 && orbit.e < 0.5) {
                    std::cout << "  ✓ Soluzione: a=" << orbit.a << " AU, e=" << orbit.e << std::endl;
                    best_herget = orbit;
                    found_herget = true;
                    break;
                } else {
                    std::cout << "  ✗ Non fisica: a=" << orbit.a << " AU, e=" << orbit.e << std::endl;
                }
            } catch (...) {
                std::cout << "  ✗ Fallito" << std::endl;
            }
        }
        
        if (found_herget) {
            printElements(best_herget, "Migliore Soluzione Herget");
        } else {
            std::cout << "\n⚠ Nessuna soluzione Herget ragionevole trovata" << std::endl;
        }
        
        // Confronto con elementi noti
        std::cout << "\n========================================" << std::endl;
        std::cout << "RIFERIMENTO: Elementi Noti" << std::endl;
        std::cout << "========================================" << std::endl;
        
        std::cout << "\n433 Eros (epoca osculating 2024-01-01):" << std::endl;
        std::cout << "  a = 1.458 AU" << std::endl;
        std::cout << "  e = 0.223" << std::endl;
        std::cout << "  i = 10.83°" << std::endl;
        std::cout << "  Ω = 304.3°" << std::endl;
        std::cout << "  ω = 178.8°" << std::endl;
        
        std::cout << "\n========================================" << std::endl;
        std::cout << "Test completato" << std::endl;
        std::cout << "========================================" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\nErrore fatale: " << e.what() << std::endl;
        return 1;
    }
}
