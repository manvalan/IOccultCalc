/**
 * @file debug_herget_inclination.cpp
 * @brief Debug sistematico del bug di inclinazione nel metodo di Herget
 * 
 * BUG: Herget calcola i=45° invece di 11° per Eros
 * Questo programma testa ogni step dell'algoritmo per trovare dove si introduce l'errore
 */

#include "ioccultcalc/initial_orbit.h"
#include "ioccultcalc/mpc_client.h"
#include "ioccultcalc/coordinates.h"
#include "ioccultcalc/spice_spk_reader.h"
#include "ioccultcalc/orbital_elements.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

// Orbital elements di riferimento per Eros (JPL Horizons)
struct ReferenceOrbit {
    double a = 1.458;      // AU
    double e = 0.223;
    double i = 10.83;      // deg
    double Om = 304.3;     // deg  
    double om = 178.8;     // deg
    double M = 320.3;      // deg (at epoch)
    double epoch = 2460000.5; // JD
};

void printVector(const std::string& name, const std::array<double, 3>& v) {
    std::cout << name << ": [" 
              << std::fixed << std::setprecision(8)
              << v[0] << ", " << v[1] << ", " << v[2] << "]\n";
}

void printElements(const std::string& title, const OrbitalElements& el) {
    std::cout << "\n" << title << ":\n";
    std::cout << "  a  = " << std::setw(12) << std::setprecision(6) << el.a << " AU\n";
    std::cout << "  e  = " << std::setw(12) << std::setprecision(6) << el.e << "\n";
    std::cout << "  i  = " << std::setw(12) << std::setprecision(4) << el.i * 180.0 / M_PI << " deg\n";
    std::cout << "  Om = " << std::setw(12) << std::setprecision(4) << el.Omega * 180.0 / M_PI << " deg\n";
    std::cout << "  om = " << std::setw(12) << std::setprecision(4) << el.omega * 180.0 / M_PI << " deg\n";
    std::cout << "  M  = " << std::setw(12) << std::setprecision(4) << el.M * 180.0 / M_PI << " deg\n";
}

void compareWithReference(const OrbitalElements& computed, const ReferenceOrbit& ref) {
    std::cout << "\n╔════════════════════════════════════════════════╗\n";
    std::cout << "║  CONFRONTO CON ORBITA DI RIFERIMENTO (JPL)    ║\n";
    std::cout << "╚════════════════════════════════════════════════╝\n";
    
    double da = std::abs(computed.a - ref.a);
    double de = std::abs(computed.e - ref.e);
    double di = std::abs(computed.i * 180.0 / M_PI - ref.i);
    
    std::cout << "\nErrori assoluti:\n";
    std::cout << "  Δa  = " << std::setw(12) << std::setprecision(6) << da << " AU     ";
    std::cout << "(" << std::setprecision(2) << (da / ref.a * 100) << "%)\n";
    
    std::cout << "  Δe  = " << std::setw(12) << std::setprecision(6) << de << "        ";
    std::cout << "(" << std::setprecision(2) << (de / ref.e * 100) << "%)\n";
    
    std::cout << "  Δi  = " << std::setw(12) << std::setprecision(4) << di << " deg   ";
    std::cout << "(" << std::setprecision(2) << (di / ref.i * 100) << "%)\n";
    
    // Verifica
    std::cout << "\nStatus:\n";
    if (di < 1.0) {
        std::cout << "  ✅ Inclinazione CORRETTA (errore < 1°)\n";
    } else if (di < 5.0) {
        std::cout << "  ⚠️  Inclinazione ACCETTABILE (errore < 5°)\n";
    } else {
        std::cout << "  ❌ Inclinazione ERRATA (errore > 5°)\n";
    }
}

int main() {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════╗\n";
    std::cout << "║  DEBUG HERGET INCLINATION - STEP BY STEP ANALYSIS   ║\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n";
    
    ReferenceOrbit ref;
    
    try {
        // Initialize SPICE ephemeris
        std::cout << "\n[1] Inizializzazione SPICE SPK...\n";
        auto spk = std::make_shared<SpiceSPKReader>();
        std::cout << "✅ SPICE SPK ready\n";
        
        // Load observations for Eros
        std::cout << "\n[2] Caricamento osservazioni Eros...\n";
        MPCClient mpc;
        auto obs_set = mpc.getObservations("433");
        
        if (obs_set.observations.empty()) {
            std::cerr << "❌ Nessuna osservazione trovata per 433 Eros\n";
            return 1;
        }
        
        std::cout << "✅ Caricate " << obs_set.observations.size() << " osservazioni\n";
        std::cout << "   Prima: " << obs_set.observations.front().jd_utc << " JD\n";
        std::cout << "   Ultima: " << obs_set.observations.back().jd_utc << " JD\n";
        
        // Select 3 observations with good spacing
        std::cout << "\n[3] Selezione tripla osservazioni...\n";
        
        // Usiamo osservazioni con spacing ottimale (~30-60 giorni)
        size_t n = obs_set.observations.size();
        size_t idx1 = 0;
        size_t idx2 = n / 2;
        size_t idx3 = n - 1;
        
        auto& obs1 = obs_set.observations[idx1];
        auto& obs2 = obs_set.observations[idx2];
        auto& obs3 = obs_set.observations[idx3];
        
        double span12 = obs2.jd_utc - obs1.jd_utc;
        double span23 = obs3.jd_utc - obs2.jd_utc;
        double span_total = obs3.jd_utc - obs1.jd_utc;
        
        std::cout << "  Obs 1: JD " << std::fixed << std::setprecision(1) << obs1.jd_utc 
                  << " (RA=" << obs1.ra << ", Dec=" << obs1.dec << ")\n";
        std::cout << "  Obs 2: JD " << obs2.jd_utc << " (+" << span12 << " giorni)\n";
        std::cout << "  Obs 3: JD " << obs3.jd_utc << " (+" << span23 << " giorni)\n";
        std::cout << "  Span totale: " << span_total << " giorni\n";
        
        // Get Earth positions
        std::cout << "\n[4] Calcolo posizioni Terra (SPICE SPK)...\n";
        std::array<double, 3> earth1, earth2, earth3;
        std::array<double, 3> earth_vel; // dummy
        
        if (!spk->getBodyPosition(3, obs1.jd_utc, earth1.data(), earth_vel.data())) {
            std::cerr << "❌ Errore calcolo posizione Terra (obs1)\n";
            return 1;
        }
        if (!spk->getBodyPosition(3, obs2.jd_utc, earth2.data(), earth_vel.data())) {
            std::cerr << "❌ Errore calcolo posizione Terra (obs2)\n";
            return 1;
        }
        if (!spk->getBodyPosition(3, obs3.jd_utc, earth3.data(), earth_vel.data())) {
            std::cerr << "❌ Errore calcolo posizione Terra (obs3)\n";
            return 1;
        }
        
        printVector("  Terra 1", earth1);
        printVector("  Terra 2", earth2);
        printVector("  Terra 3", earth3);
        
        // Convert RA/Dec to unit vectors
        std::cout << "\n[5] Conversione RA/Dec → vettori direzione...\n";
        std::array<double, 3> rho_hat1, rho_hat2, rho_hat3;
        
        rho_hat1[0] = std::cos(obs1.dec) * std::cos(obs1.ra);
        rho_hat1[1] = std::cos(obs1.dec) * std::sin(obs1.ra);
        rho_hat1[2] = std::sin(obs1.dec);
        
        rho_hat2[0] = std::cos(obs2.dec) * std::cos(obs2.ra);
        rho_hat2[1] = std::cos(obs2.dec) * std::sin(obs2.ra);
        rho_hat2[2] = std::sin(obs2.dec);
        
        rho_hat3[0] = std::cos(obs3.dec) * std::cos(obs3.ra);
        rho_hat3[1] = std::cos(obs3.dec) * std::sin(obs3.ra);
        rho_hat3[2] = std::sin(obs3.dec);
        
        printVector("  rho_hat 1", rho_hat1);
        printVector("  rho_hat 2", rho_hat2);
        printVector("  rho_hat 3", rho_hat3);
        
        // Run Herget method
        std::cout << "\n[6] Esecuzione Herget method...\n";
        std::cout << "──────────────────────────────────────────────────────\n";
        
        // Initial guesses: ~1 AU for NEA
        double r1_guess = 1.2;
        double r2_guess = 1.5;
        
        OrbitalElements elements = InitialOrbit::herget(
            obs_set.observations,
            r1_guess,
            r2_guess,
            10  // max iterations
        );
        
        std::cout << "──────────────────────────────────────────────────────\n";
        std::cout << "✅ Herget method completato\n";
        
        // Print results
        printElements("ELEMENTI ORBITALI CALCOLATI", elements);
        
        // Compare with reference
        compareWithReference(elements, ref);
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ Eccezione: " << e.what() << "\n";
        return 1;
    }
    
    std::cout << "\n";
    return 0;
}
