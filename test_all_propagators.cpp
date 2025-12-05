#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/orbit_fitter.h"
#include "ioccultcalc/time_utils.h"
#include "ioccultcalc/jpl_horizons_client.h"

using namespace ioccultcalc;

// Carica elementi da file
EquinoctialElements loadAstDysElements() {
    EquinoctialElements astElem;
    std::ifstream eq1File("17030_astdys.eq1");
    
    if (!eq1File.good()) {
        throw std::runtime_error("File 17030_astdys.eq1 non trovato!");
    }
    
    std::string line;
    bool found = false;
    while (std::getline(eq1File, line)) {
        if (line.find(" EQU ") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            double a, h, k, p, q, lambda;
            if (iss >> tag >> a >> h >> k >> p >> q >> lambda) {
                astElem.a = a;
                astElem.h = h;
                astElem.k = k;
                astElem.p = p;
                astElem.q = q;
                astElem.lambda = lambda;
                found = true;
            }
        } else if (found && line.find(" MJD ") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            double mjd;
            if (iss >> tag >> mjd) {
                astElem.epoch.jd = mjd + 2400000.5;
                break;
            }
        }
    }
    eq1File.close();
    
    if (!found) {
        throw std::runtime_error("Elementi non trovati nel file!");
    }
    
    return astElem;
}

int main() {
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "TEST CONFRONTO PROPAGATORI vs JPL HORIZONS\n";
    std::cout << "═══════════════════════════════════════════════════════════\n\n";
    
    // Epoca target
    JulianDate targetEpoch;
    targetEpoch.jd = 2461007.5;  // 28 Nov 2025
    
    int year, month, day, hour, minute;
    double second;
    TimeUtils::jdToCalendar(targetEpoch, year, month, day, hour, minute, second);
    
    std::cout << "Data target: " << year << "-" << std::setw(2) << std::setfill('0') 
              << month << "-" << day << " 00:00 UTC\n";
    std::cout << "JD: " << std::fixed << std::setprecision(1) << targetEpoch.jd << "\n\n";
    
    // Carica elementi
    std::cout << "Caricamento elementi AstDyS...\n";
    EquinoctialElements astElem = loadAstDysElements();
    
    int elemYear, elemMonth, elemDay, elemHour, elemMinute;
    double elemSecond;
    TimeUtils::jdToCalendar(astElem.epoch, elemYear, elemMonth, elemDay, elemHour, elemMinute, elemSecond);
    
    std::cout << "Epoca elementi: " << elemYear << "-" << std::setw(2) << std::setfill('0') 
              << elemMonth << "-" << elemDay << " (JD " << astElem.epoch.jd << ")\n";
    std::cout << "Delta tempo: " << std::setprecision(1) 
              << (targetEpoch.jd - astElem.epoch.jd) << " giorni\n\n";
    
    auto kepElem = astElem.toKeplerian();
    std::cout << "Elementi Keplerian:\n";
    std::cout << "  a = " << std::setprecision(8) << kepElem.a << " AU\n";
    std::cout << "  e = " << std::setprecision(6) << kepElem.e << "\n";
    std::cout << "  i = " << std::setprecision(4) << (kepElem.i * RAD_TO_DEG) << "°\n";
    std::cout << "  Ω = " << (kepElem.Omega * RAD_TO_DEG) << "°\n";
    std::cout << "  ω = " << (kepElem.omega * RAD_TO_DEG) << "°\n";
    std::cout << "  M = " << (kepElem.M * RAD_TO_DEG) << "°\n\n";
    
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "PROPAGAZIONE CON DIVERSI METODI\n";
    std::cout << "═══════════════════════════════════════════════════════════\n\n";
    
    // ================================================================
    // 1. OrbitPropagator con RKF78
    // ================================================================
    std::cout << "1. OrbitPropagator + RKF78 (step=0.01, tol=1e-12)\n";
    std::cout << "   " << std::string(55, '-') << "\n";
    
    try {
        Ephemeris eph1(astElem);
        eph1.enableNumericalPropagation(true, IntegratorType::RKF78);
        
        PropagatorOptions opts1;
        opts1.integrator = IntegratorType::RKF78;
        opts1.stepSize = 0.01;
        opts1.tolerance = 1e-12;
        opts1.usePlanetaryPerturbations = true;
        eph1.setPropagatorOptions(opts1);
        
        EphemerisData data1 = eph1.compute(targetEpoch);
        
        std::cout << "   RA  = " << std::fixed << std::setprecision(6) 
                  << (data1.geocentricPos.ra * RAD_TO_DEG) << "°\n";
        std::cout << "   Dec = " << (data1.geocentricPos.dec * RAD_TO_DEG) << "°\n";
        std::cout << "   Dist = " << std::setprecision(8) << data1.distance << " AU\n\n";
    } catch (const std::exception& e) {
        std::cout << "   ❌ ERRORE: " << e.what() << "\n\n";
    }
    
    // ================================================================
    // 2. OrbitPropagator con RK4
    // ================================================================
    std::cout << "2. OrbitPropagator + RK4 (step=0.1)\n";
    std::cout << "   " << std::string(55, '-') << "\n";
    
    try {
        Ephemeris eph2(astElem);
        eph2.enableNumericalPropagation(true, IntegratorType::RK4);
        
        PropagatorOptions opts2;
        opts2.integrator = IntegratorType::RK4;
        opts2.stepSize = 0.1;
        opts2.usePlanetaryPerturbations = true;
        eph2.setPropagatorOptions(opts2);
        
        EphemerisData data2 = eph2.compute(targetEpoch);
        
        std::cout << "   RA  = " << std::fixed << std::setprecision(6) 
                  << (data2.geocentricPos.ra * RAD_TO_DEG) << "°\n";
        std::cout << "   Dec = " << (data2.geocentricPos.dec * RAD_TO_DEG) << "°\n";
        std::cout << "   Dist = " << std::setprecision(8) << data2.distance << " AU\n\n";
    } catch (const std::exception& e) {
        std::cout << "   ❌ ERRORE: " << e.what() << "\n\n";
    }
    
    // ================================================================
    // 3. OrbitPropagator con GAUSS_RADAU
    // ================================================================
    std::cout << "3. OrbitPropagator + GAUSS_RADAU\n";
    std::cout << "   " << std::string(55, '-') << "\n";
    
    try {
        Ephemeris eph3(astElem);
        eph3.enableNumericalPropagation(true, IntegratorType::GAUSS_RADAU);
        
        PropagatorOptions opts3;
        opts3.integrator = IntegratorType::GAUSS_RADAU;
        opts3.stepSize = 0.1;
        opts3.tolerance = 1e-12;
        opts3.usePlanetaryPerturbations = true;
        eph3.setPropagatorOptions(opts3);
        
        EphemerisData data3 = eph3.compute(targetEpoch);
        
        std::cout << "   RA  = " << std::fixed << std::setprecision(6) 
                  << (data3.geocentricPos.ra * RAD_TO_DEG) << "°\n";
        std::cout << "   Dec = " << (data3.geocentricPos.dec * RAD_TO_DEG) << "°\n";
        std::cout << "   Dist = " << std::setprecision(8) << data3.distance << " AU\n\n";
    } catch (const std::exception& e) {
        std::cout << "   ❌ ERRORE: " << e.what() << "\n\n";
    }
    
    // ================================================================
    // 4. Propagazione Kepleriana (nessuna perturbazione)
    // ================================================================
    std::cout << "4. Propagazione Kepleriana (2-body, nessuna perturbazione)\n";
    std::cout << "   " << std::string(55, '-') << "\n";
    
    try {
        Ephemeris eph4(astElem);
        eph4.enableNumericalPropagation(false);  // Usa Keplerian
        
        EphemerisData data4 = eph4.compute(targetEpoch);
        
        std::cout << "   RA  = " << std::fixed << std::setprecision(6) 
                  << (data4.geocentricPos.ra * RAD_TO_DEG) << "°\n";
        std::cout << "   Dec = " << (data4.geocentricPos.dec * RAD_TO_DEG) << "°\n";
        std::cout << "   Dist = " << std::setprecision(8) << data4.distance << " AU\n\n";
    } catch (const std::exception& e) {
        std::cout << "   ❌ ERRORE: " << e.what() << "\n\n";
    }
    
    // ================================================================
    // 5. OrbitFitter::propagate() - SALTATO (causa crash)
    // ================================================================
    std::cout << "5. OrbitFitter::propagate() - [SALTATO]\n";
    std::cout << "   " << std::string(55, '-') << "\n";
    std::cout << "   (Causa crash - skip per ora)\n\n";
    
    // ================================================================
    // 6. JPL Horizons (RIFERIMENTO)
    // ================================================================
    std::cout << "6. JPL Horizons (RIFERIMENTO)\n";
    std::cout << "   " << std::string(55, '-') << "\n";
    
    try {
        JPLHorizonsClient horizons;
        
        // Query per (17030) Sierks - coordinate apparenti geocentriche
        std::cout << "   Query in corso...\n";
        
        auto coords = horizons.getApparentCoordinates("17030", targetEpoch, "500");
        
        double ra_deg = coords.first * RAD_TO_DEG;
        double dec_deg = coords.second * RAD_TO_DEG;
        
        std::cout << "   RA  = " << std::fixed << std::setprecision(6) << ra_deg << "°\n";
        std::cout << "   Dec = " << dec_deg << "°\n";
        
        // Prova anche a ottenere la distanza
        try {
            auto eph = horizons.getEphemeris("17030", targetEpoch, "500@0");
            std::cout << "   Dist = " << std::setprecision(8) << eph.distance << " AU\n\n";
        } catch (...) {
            std::cout << "   Dist = (non disponibile)\n\n";
        }
    } catch (const std::exception& e) {
        std::cout << "   ❌ ERRORE: " << e.what() << "\n\n";
    }
    
    // ================================================================
    // STELLA TARGET
    // ================================================================
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "STELLA TARGET: Gaia DR3 3411546266140512128\n";
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "RA  = 73.416106°\n";
    std::cout << "Dec = 20.331662°\n\n";
    
    std::cout << "Per occultazione, l'asteroide deve essere entro ~1\" dalla stella.\n";
    std::cout << "Attualmente vediamo separazioni di ~230°!\n\n";
    
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "TEST COMPLETATO\n";
    std::cout << "═══════════════════════════════════════════════════════════\n";
    
    return 0;
}
