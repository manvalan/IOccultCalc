/**
 * Test: Confronto posizione 17030 tra RKF78 e JPL Horizons
 * Data: 2025-11-28 per verificare occultazione
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/coordinates.h"

using namespace ioccultcalc;

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    // Carica elementi orbitali di 17030 dal catalogo
    std::ifstream file("/Users/michelebigi/.ioccultcalc/data/all_numbered_asteroids.json");
    if (!file.is_open()) {
        std::cerr << "Errore: impossibile aprire catalogo asteroidi\n";
        return 1;
    }
    
    std::string line;
    ioccultcalc::OrbitalElements elements;
    bool found = false;
    
    // Parse JSON manuale per trovare 17030
    while (std::getline(file, line)) {
        if (line.find("\"number\": 17030") != std::string::npos) {
            found = true;
            std::cout << "=== ELEMENTI ORBITALI 17030 Sierks ===\n";
            
            // Leggi le prossime righe per gli elementi
            for (int i = 0; i < 30 && std::getline(file, line); i++) {
                if (line.find("\"epoch_jd\"") != std::string::npos) {
                    sscanf(line.c_str(), " \"epoch_jd\": %lf", &elements.epoch.jd);
                    std::cout << "Epoca: JD " << elements.epoch.jd << "\n";
                }
                else if (line.find("\"a\"") != std::string::npos && line.find("\"a\"") == line.find_first_not_of(" \t")) {
                    sscanf(line.c_str(), " \"a\": %lf", &elements.a);
                    std::cout << "a = " << elements.a << " AU\n";
                }
                else if (line.find("\"e\"") != std::string::npos && line.find("\"e\"") == line.find_first_not_of(" \t")) {
                    sscanf(line.c_str(), " \"e\": %lf", &elements.e);
                    std::cout << "e = " << elements.e << "\n";
                }
                else if (line.find("\"i\"") != std::string::npos && line.find("\"i\"") == line.find_first_not_of(" \t")) {
                    sscanf(line.c_str(), " \"i\": %lf", &elements.i);
                    elements.i *= M_PI / 180.0;  // Converti in radianti
                    std::cout << "i = " << (elements.i * 180.0 / M_PI) << "°\n";
                }
                else if (line.find("\"Omega\"") != std::string::npos) {
                    sscanf(line.c_str(), " \"Omega\": %lf", &elements.Omega);
                    elements.Omega *= M_PI / 180.0;
                    std::cout << "Ω = " << (elements.Omega * 180.0 / M_PI) << "°\n";
                }
                else if (line.find("\"omega\"") != std::string::npos && line.find("\"omega\"") == line.find_first_not_of(" \t")) {
                    sscanf(line.c_str(), " \"omega\": %lf", &elements.omega);
                    elements.omega *= M_PI / 180.0;
                    std::cout << "ω = " << (elements.omega * 180.0 / M_PI) << "°\n";
                }
                else if (line.find("\"M\"") != std::string::npos && line.find("\"M\"") == line.find_first_not_of(" \t")) {
                    sscanf(line.c_str(), " \"M\": %lf", &elements.M);
                    elements.M *= M_PI / 180.0;
                    std::cout << "M = " << (elements.M * 180.0 / M_PI) << "°\n";
                }
                else if (line.find("\"H\"") != std::string::npos && line.find("\"H\"") == line.find_first_not_of(" \t")) {
                    sscanf(line.c_str(), " \"H\": %lf", &elements.H);
                    std::cout << "H = " << elements.H << "\n";
                }
                else if (line.find("\"name\"") != std::string::npos) {
                    size_t start = line.find(":") + 1;
                    size_t end = line.find_last_of("\"");
                    elements.name = line.substr(start, end - start);
                    // Rimuovi spazi e virgolette
                    elements.name.erase(0, elements.name.find_first_not_of(" \t\""));
                    elements.name.erase(elements.name.find_last_not_of(" \t\",") + 1);
                    std::cout << "Nome: " << elements.name << "\n";
                }
            }
            break;
        }
    }
    file.close();
    
    if (!found) {
        std::cerr << "Asteroide 17030 non trovato nel catalogo!\n";
        return 1;
    }
    
    elements.designation = "17030";
    
    // Converti a elementi equinoziali per propagazione
    EquinoctialElements eqElements = EquinoctialElements::fromKeplerian(elements);
    
    // Crea Ephemeris (usa RKF78 internamente)
    Ephemeris ephemeris(eqElements);
    
    std::cout << "\n=== PROPAGAZIONE RKF78 via Ephemeris ===\n";
    std::cout << "Configurazione:\n";
    std::cout << "  Integratore: RKF78 (default)\n";
    std::cout << "  Perturbazioni planetarie: SI\n\n";
    
    // Date di test
    struct TestDate {
        std::string label;
        double jd;
    };
    
    std::vector<TestDate> dates = {
        {"2025-11-27 00:00", 2460641.5},
        {"2025-11-28 00:00", 2460642.5},
        {"2025-11-28 12:00", 2460643.0},
        {"2025-11-29 00:00", 2460643.5}
    };
    
    std::cout << "Data            JD          RA (deg)    Dec (deg)   RA (h:m:s)        Distanza (AU)\n";
    std::cout << "================================================================================\n";
    
    for (const auto& date : dates) {
        JulianDate epoch;
        epoch.jd = date.jd;
        
        // Calcola ephemeris
        EphemerisData ephData = ephemeris.compute(epoch);
        
        // Posizione geocentrica già in equatoriale
        double ra_deg = ephData.geocentricPos.ra * 180.0 / M_PI;
        double dec_deg = ephData.geocentricPos.dec * 180.0 / M_PI;
        double ra_hours = ra_deg / 15.0;
        
        int h = (int)ra_hours;
        int m = (int)((ra_hours - h) * 60);
        double s = ((ra_hours - h) * 60 - m) * 60;
        
        double dist = ephData.distance;
        
        std::cout << std::setw(15) << date.label 
                  << std::setw(12) << std::fixed << std::setprecision(2) << date.jd
                  << std::setw(12) << std::setprecision(4) << ra_deg
                  << std::setw(12) << std::setprecision(4) << dec_deg
                  << "   " << std::setw(2) << std::setfill('0') << h << ":"
                  << std::setw(2) << m << ":"
                  << std::setw(5) << std::setfill(' ') << std::setprecision(2) << s
                  << std::setw(12) << std::setprecision(6) << dist << "\n";
    }
    
    std::cout << "\n";
    std::cout << "=== CONFRONTA CON JPL HORIZONS ===\n";
    std::cout << "Esegui questa query su https://ssd.jpl.nasa.gov/horizons.cgi :\n";
    std::cout << "  Target: 17030 Sierks\n";
    std::cout << "  Observer: Geocentric (500@0)\n";
    std::cout << "  Date: 2025-11-28 00:00\n";
    std::cout << "  Quantities: 1 (RA/Dec)\n\n";
    
    std::cout << "Stella target UCAC4 552-011427 (Gaia DR3 3411546266140512128):\n";
    std::cout << "  RA  = 73.4161° = 4h 53m 39.9s\n";
    std::cout << "  Dec = +20.3317° = +20° 19' 54\"\n";
    
    return 0;
}
