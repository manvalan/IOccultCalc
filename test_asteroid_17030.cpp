#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/time_utils.h"
#include <iostream>
#include <iomanip>

using namespace ioccultcalc;

int main() {
    try {
        // Carica elementi asteroide 17030
        AstDySClient client;
        std::string dataDir = std::string(getenv("HOME")) + "/.ioccultcalc/data";
        
        std::cout << "Caricamento elementi asteroide 17030 Sierks..." << std::endl;
        OrbitalElements elem = client.loadFromLocalFile(dataDir, 17030);
        
        std::cout << "\nElementi orbitali:" << std::endl;
        std::cout << "  Numero: " << elem.number << std::endl;
        std::cout << "  Nome: " << elem.name << std::endl;
        std::cout << "  H: " << elem.H << std::endl;
        std::cout << "  Semiasse: " << elem.a << " AU" << std::endl;
        std::cout << "  Eccentricità: " << elem.e << std::endl;
        std::cout << "  Inclinazione: " << elem.i * 180/M_PI << "°" << std::endl;
        
        // Calcola posizione per oggi (28 Nov 2025)
        JulianDate today = TimeUtils::calendarToJulian(2025, 11, 28, 12, 0, 0);
        std::cout << "\nCalcolo posizione per 28 Nov 2025 12:00 UTC (JD " 
                  << std::fixed << std::setprecision(2) << today.jd << ")..." << std::endl;
        
        Ephemeris eph(elem);
        EphemerisData data = eph.compute(today);
        
        std::cout << "\nPosizione geocentrica:" << std::endl;
        std::cout << "  RA: " << std::fixed << std::setprecision(4) 
                  << data.geocentricPos.ra * 180/M_PI << "°" << std::endl;
        std::cout << "  Dec: " << data.geocentricPos.dec * 180/M_PI << "°" << std::endl;
        std::cout << "  Distanza: " << std::setprecision(3) << data.distance << " AU" << std::endl;
        std::cout << "  Magnitudine: " << std::setprecision(1) << data.magnitude << std::endl;
        std::cout << "  Elongazione solare: " << std::setprecision(1) << data.elongation << "°" << std::endl;
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Errore: " << e.what() << std::endl;
        return 1;
    }
}
