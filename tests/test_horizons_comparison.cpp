/**
 * @file test_horizons_comparison.cpp
 * @brief Confronto propagazione con JPL Horizons per asteroide 17030
 */

#include "ioccultcalc/astdyn_propagation_helper.h"
#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/time_utils.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace ioccultcalc;

// Colori per output
#define RESET   "\033[0m"
#define GREEN   "\033[32m"
#define RED     "\033[31m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define CYAN    "\033[36m"

void printHeader(const std::string& title) {
    std::cout << "\n" << BLUE << "═══════════════════════════════════════════════════════\n";
    std::cout << title << "\n";
    std::cout << "═══════════════════════════════════════════════════════" << RESET << "\n\n";
}

void printValue(const std::string& label, double value, int precision = 6) {
    std::cout << "  " << CYAN << label << ": " << RESET 
              << std::fixed << std::setprecision(precision) << value << "\n";
}

int main() {
    std::cout << GREEN << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  CONFRONTO CON JPL HORIZONS - ASTEROIDE 17030        ║\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n" << RESET;
    
    try {
        // Carica elementi
        printHeader("Caricamento Elementi");
        AstDySElements elements = AstDySClient::downloadElements(17030);
        if (elements.a <= 0) {
            std::cerr << RED << "✗ Errore caricamento elementi" << RESET << "\n";
            return 1;
        }
        std::cout << GREEN << "✓ Elementi caricati\n" << RESET;
        printValue("Epoch (MJD TDB)", elements.epoch_mjd, 5);
        printValue("a (AU)", elements.a, 8);
        printValue("e", elements.e, 8);
        
        // Date di test: 26-30 novembre 2025
        std::vector<double> test_mjds = {
            61005.0,  // 26 nov 00:00
            61006.0,  // 27 nov 00:00
            61007.5,  // 28 nov 12:00
            61008.0,  // 29 nov 00:00
            61009.0   // 30 nov 00:00
        };
        
        printHeader("Propagazione e Confronto con Horizons");
        std::cout << "  " << YELLOW << "NOTA: Confronta manualmente con JPL Horizons\n" << RESET;
        std::cout << "  " << YELLOW << "Horizons query: 17030, 2025-11-26 00:00 to 2025-11-30 00:00\n\n" << RESET;
        
        auto& helper = AstDynPropagationHelper::getInstance();
        
        std::cout << std::left << std::setw(12) << "MJD TDB" 
                  << std::setw(15) << "RA (deg)" 
                  << std::setw(15) << "Dec (deg)"
                  << std::setw(20) << "Distance (AU)" << "\n";
        std::cout << std::string(62, '-') << "\n";
        
        for (double mjd : test_mjds) {
            auto radec = helper.getRADec(elements, mjd);
            
            // Converti MJD in data leggibile
            JulianDate jd(mjd + 2400000.5);
            int year, month, day, hour, minute;
            double second;
            TimeUtils::jdToCalendar(jd, year, month, day, hour, minute, second);
            
            std::cout << std::fixed << std::setprecision(5) << std::setw(12) << mjd
                      << std::setprecision(8) << std::setw(15) << radec.first
                      << std::setprecision(8) << std::setw(15) << radec.second;
            
            // Calcola distanza geocentrica (approssimativa)
            auto propagated = helper.propagate(elements, mjd);
            double r_approx = propagated.a * (1.0 - propagated.e);  // Approssimazione
            std::cout << std::setprecision(8) << std::setw(20) << r_approx << "\n";
        }
        
        std::cout << "\n  " << CYAN << "Per confrontare con Horizons:\n" << RESET;
        std::cout << "  1. Vai su https://ssd.jpl.nasa.gov/horizons/\n";
        std::cout << "  2. Target Body: 17030\n";
        std::cout << "  3. Observer: Geocentric\n";
        std::cout << "  4. Time Span: 2025-11-26 00:00 to 2025-11-30 00:00\n";
        std::cout << "  5. Table Settings: RA/Dec (ICRF)\n\n";
        
        // Verifica quale propagatore viene usato
        printHeader("Configurazione Propagatore");
        std::cout << "  " << CYAN << "AstDynPropagationHelper usa:\n" << RESET;
        std::cout << "    - AstDynPropagator con tolleranza 1e-12\n";
        std::cout << "    - RKF78 integrator\n";
        std::cout << "    - Tutte le perturbazioni planetarie\n";
        std::cout << "    - Correzioni relativistiche\n";
        std::cout << "    - Perturbazioni asteroidi (AST17)\n";
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Exception: " << e.what() << RESET << "\n";
        return 1;
    }
    
    std::cout << "\n" << GREEN << "✓ Test completato\n" << RESET;
    return 0;
}

