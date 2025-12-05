/**
 * @file test_17030_diagnostics.cpp
 * @brief Test diagnostico per capire il problema di accuratezza
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ioccultcalc/orbit_propagator.h>
#include <ioccultcalc/orbital_elements.h>
#include <ioccultcalc/time_utils.h>
#include <ioccultcalc/jpl_horizons_client.h>

using namespace ioccultcalc;

bool loadAstDysElements(const std::string& filename, EquinoctialElements& elements) {
    std::ifstream file(filename);
    if (!file.good()) {
        std::cerr << "Impossibile leggere " << filename << "\n";
        return false;
    }
    
    std::string line;
    bool foundEqu = false, foundMjd = false;
    
    while (std::getline(file, line)) {
        if (line.find(" EQU ") != std::string::npos && !foundEqu) {
            std::istringstream iss(line);
            std::string tag;
            double a, h, k, p, q, lambda;
            if (iss >> tag >> a >> h >> k >> p >> q >> lambda) {
                elements.a = a;
                elements.h = h;
                elements.k = k;
                elements.p = p;
                elements.q = q;
                elements.lambda = lambda;
                foundEqu = true;
            }
        }
        if (line.find(" MJD ") != std::string::npos && foundEqu && !foundMjd) {
            std::istringstream iss(line);
            std::string tag;
            double mjd;
            if (iss >> tag >> mjd) {
                elements.epoch.jd = mjd + 2400000.5;
                foundMjd = true;
                break;
            }
        }
    }
    file.close();
    
    return foundEqu && foundMjd;
}

int main() {
    std::cout << "===== TEST DIAGNOSTICO 17030 =====\n\n";
    
    try {
        // Carica elementi
        EquinoctialElements astElements;
        if (!loadAstDysElements("17030_astdys.eq1", astElements)) {
            std::cerr << "Errore caricamento elementi\n";
            return 1;
        }
        
        auto kepElements = astElements.toKeplerian();
        std::cout << "Elementi AstDyS caricati:\n";
        std::cout << "  Epoca: JD " << std::fixed << std::setprecision(2) << astElements.epoch.jd << "\n";
        std::cout << "  a = " << std::setprecision(8) << kepElements.a << " AU\n";
        std::cout << "  e = " << kepElements.e << "\n";
        std::cout << "  i = " << (kepElements.i * 180.0 / M_PI) << "°\n\n";
        
        // Setup propagatore
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RKF78;
        opts.stepSize = 0.01;
        opts.tolerance = 1e-12;
        opts.usePlanetaryPerturbations = true;
        opts.useRelativisticCorrections = true;
        
        OrbitPropagator propagator(opts);
        
        // Propaga a una data specifica
        JulianDate testDate;
        testDate.jd = 2461000.5;  // 26 Nov 2025
        
        std::cout << "Propagazione da JD " << astElements.epoch.jd 
                  << " a JD " << testDate.jd << "\n";
        std::cout << "Delta = " << (testDate.jd - astElements.epoch.jd) << " giorni\n\n";
        
        auto stateIOC = propagator.propagate(astElements, testDate);
        
        std::cout << "RISULTATO IOccultCalc (RKF78):\n";
        std::cout << "  Posizione (AU):       (" << std::setprecision(10)
                  << stateIOC.position.x << ", "
                  << stateIOC.position.y << ", "
                  << stateIOC.position.z << ")\n";
        std::cout << "  Modulo posizione:     " << stateIOC.position.magnitude() << " AU\n";
        std::cout << "  Velocità (AU/day):    ("
                  << stateIOC.velocity.x << ", "
                  << stateIOC.velocity.y << ", "
                  << stateIOC.velocity.z << ")\n";
        std::cout << "  Modulo velocità:      " << stateIOC.velocity.magnitude() << " AU/day\n\n";
        
        // Scarica da JPL
        JPLHorizonsClient horizons;
        horizons.setTimeout(60);
        
        std::cout << "Scaricamento da JPL Horizons...\n";
        auto [posJPL, velJPL] = horizons.getStateVectors("17030", testDate, "@sun");
        
        std::cout << "RISULTATO JPL Horizons:\n";
        std::cout << "  Posizione (AU):       (" << std::setprecision(10)
                  << posJPL.x << ", "
                  << posJPL.y << ", "
                  << posJPL.z << ")\n";
        std::cout << "  Modulo posizione:     " << posJPL.magnitude() << " AU\n";
        std::cout << "  Velocità (AU/day):    ("
                  << velJPL.x << ", "
                  << velJPL.y << ", "
                  << velJPL.z << ")\n";
        std::cout << "  Modulo velocità:      " << velJPL.magnitude() << " AU/day\n\n";
        
        // Confronto
        std::cout << "CONFRONTO:\n";
        std::cout << "  ΔX = " << std::scientific << std::setprecision(6) 
                  << (stateIOC.position.x - posJPL.x) << " AU = "
                  << std::fixed << (stateIOC.position.x - posJPL.x) * 149597870.7 << " km\n";
        std::cout << "  ΔY = " << std::scientific << (stateIOC.position.y - posJPL.y) 
                  << " AU = " << std::fixed 
                  << (stateIOC.position.y - posJPL.y) * 149597870.7 << " km\n";
        std::cout << "  ΔZ = " << std::scientific << (stateIOC.position.z - posJPL.z) 
                  << " AU = " << std::fixed 
                  << (stateIOC.position.z - posJPL.z) * 149597870.7 << " km\n\n";
        
        Vector3D deltaPos = stateIOC.position - posJPL;
        std::cout << "  |Δpos| = " << std::scientific << deltaPos.magnitude() 
                  << " AU = " << std::fixed << deltaPos.magnitude() * 149597870.7 << " km\n";
        std::cout << "  Percentuale: " << (deltaPos.magnitude() / posJPL.magnitude() * 100.0) 
                  << "%\n\n";
        
        // Converti in coordinate sferiche per capire meglio
        std::cout << "COORDINATE SFERICHE:\n";
        
        double raIOC = atan2(stateIOC.position.y, stateIOC.position.x) * 180.0 / M_PI;
        double decIOC = asin(stateIOC.position.z / stateIOC.position.magnitude()) * 180.0 / M_PI;
        
        double raJPL = atan2(posJPL.y, posJPL.x) * 180.0 / M_PI;
        double decJPL = asin(posJPL.z / posJPL.magnitude()) * 180.0 / M_PI;
        
        std::cout << "  IOC:       RA = " << raIOC << "°,  Dec = " << decIOC << "°\n";
        std::cout << "  JPL:       RA = " << raJPL << "°,  Dec = " << decJPL << "°\n";
        std::cout << "  ΔRA = " << (raIOC - raJPL) * 3600.0 << " arcsec\n";
        std::cout << "  ΔDec = " << (decIOC - decJPL) * 3600.0 << " arcsec\n\n";
        
        // Verifica velocità radiale
        std::cout << "VELOCITA' RADIALE:\n";
        double vr_ioc = (stateIOC.position.x * stateIOC.velocity.x + 
                         stateIOC.position.y * stateIOC.velocity.y +
                         stateIOC.position.z * stateIOC.velocity.z) / stateIOC.position.magnitude();
        double vr_jpl = (posJPL.x * velJPL.x + 
                         posJPL.y * velJPL.y +
                         posJPL.z * velJPL.z) / posJPL.magnitude();
        
        std::cout << "  IOC: " << vr_ioc << " AU/day\n";
        std::cout << "  JPL: " << vr_jpl << " AU/day\n";
        std::cout << "  ΔVr = " << (vr_ioc - vr_jpl) << " AU/day\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
