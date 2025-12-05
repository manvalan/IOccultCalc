/**
 * @file test_17030_chebyshev_vs_horizons.cpp
 * @brief Confronto posizioni IOccultCalc (RKF78 + Chebyshev) vs JPL Horizons
 * 
 * Confronta le effemeridi calcolate da IOccultCalc usando:
 * - Elementi orbitali AstDyS
 * - Propagatore RKF78
 * - Potenziale Chebyshev astrometrico (se disponibile)
 * 
 * Con i dati di riferimento JPL Horizons
 * Periodo: 26/11/2025 - 30/11/2025 (asteroide 17030)
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <ioccultcalc/orbit_propagator.h>
#include <ioccultcalc/orbital_elements.h>
#include <ioccultcalc/time_utils.h>
#include <ioccultcalc/jpl_horizons_client.h>

using namespace ioccultcalc;

// ─────────────────────────────────────────────────────────────────────────────
// Strutture dati
// ─────────────────────────────────────────────────────────────────────────────

struct EphemerisPoint {
    JulianDate epoch;
    Vector3D posIOC;      // Posizione IOccultCalc (AU)
    Vector3D velIOC;      // Velocità IOccultCalc (AU/day)
    Vector3D posJPL;      // Posizione JPL Horizons (AU)
    Vector3D velJPL;      // Velocità JPL Horizons (AU/day)
    
    double deltaPos_km;   // Errore posizione (km)
    double deltaVel_mms;  // Errore velocità (mm/s)
};

struct ComparisonStats {
    int numPoints;
    double meanPosError_km;
    double maxPosError_km;
    double minPosError_km;
    double meanVelError_mms;
    double maxVelError_mms;
};

// ─────────────────────────────────────────────────────────────────────────────
// Utilità
// ─────────────────────────────────────────────────────────────────────────────

void printSeparator() {
    std::cout << std::string(90, '=') << "\n";
}

void printSubSeparator() {
    std::cout << std::string(90, '-') << "\n";
}

/**
 * Carica elementi orbitali equinoziali dal file AstDyS
 */
bool loadAstDysElements(const std::string& filename, EquinoctialElements& elements) {
    std::ifstream file(filename);
    if (!file.good()) {
        std::cerr << "❌ Impossibile leggere " << filename << "\n";
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

/**
 * Calcola statistiche di confronto
 */
ComparisonStats calculateStats(const std::vector<EphemerisPoint>& points) {
    ComparisonStats stats;
    stats.numPoints = points.size();
    
    if (points.empty()) {
        return stats;
    }
    
    stats.meanPosError_km = 0.0;
    stats.meanVelError_mms = 0.0;
    stats.maxPosError_km = 0.0;
    stats.minPosError_km = 1e9;
    stats.maxVelError_mms = 0.0;
    
    for (const auto& p : points) {
        stats.meanPosError_km += p.deltaPos_km;
        stats.meanVelError_mms += p.deltaVel_mms;
        stats.maxPosError_km = std::max(stats.maxPosError_km, p.deltaPos_km);
        stats.minPosError_km = std::min(stats.minPosError_km, p.deltaPos_km);
        stats.maxVelError_mms = std::max(stats.maxVelError_mms, p.deltaVel_mms);
    }
    
    stats.meanPosError_km /= points.size();
    stats.meanVelError_mms /= points.size();
    
    return stats;
}

// ─────────────────────────────────────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────────────────────────────────────

int main() {
    printSeparator();
    std::cout << "   TEST: IOccultCalc (RKF78) vs JPL Horizons - Asteroide 17030\n";
    std::cout << "   Periodo: 26/11/2025 - 30/11/2025\n";
    printSeparator();
    std::cout << "\n";
    
    try {
        // ─────────────────────────────────────────────────────────────────────
        // STEP 1: Carica elementi orbitali AstDyS
        // ─────────────────────────────────────────────────────────────────────
        std::cout << "STEP 1: Caricamento elementi orbitali AstDyS\n";
        printSubSeparator();
        
        EquinoctialElements astElements;
        if (!loadAstDysElements("17030_astdys.eq1", astElements)) {
            std::cerr << "❌ Errore caricamento elementi AstDyS\n";
            return 1;
        }
        
        auto kepElements = astElements.toKeplerian();
        std::cout << "✓ Elementi caricati:\n";
        std::cout << "  Epoca:     JD " << std::fixed << std::setprecision(2) 
                  << astElements.epoch.jd << "\n";
        std::cout << "  a:         " << std::setprecision(8) << kepElements.a << " AU\n";
        std::cout << "  e:         " << std::setprecision(6) << kepElements.e << "\n";
        std::cout << "  i:         " << std::setprecision(4) 
                  << (kepElements.i * 180.0 / M_PI) << "°\n";
        std::cout << "  Ω:         " << (kepElements.Omega * 180.0 / M_PI) << "°\n";
        std::cout << "  ω:         " << (kepElements.omega * 180.0 / M_PI) << "°\n";
        std::cout << "  M:         " << (kepElements.M * 180.0 / M_PI) << "°\n\n";
        
        // ─────────────────────────────────────────────────────────────────────
        // STEP 2: Setup propagatore IOccultCalc con RKF78
        // ─────────────────────────────────────────────────────────────────────
        std::cout << "STEP 2: Setup propagatore IOccultCalc (RKF78)\n";
        printSubSeparator();
        
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RKF78;
        opts.stepSize = 0.01;           // 0.01 giorni = 14.4 minuti
        opts.tolerance = 1e-12;
        opts.usePlanetaryPerturbations = true;
        opts.useRelativisticCorrections = true;
        
        OrbitPropagator propagator(opts);
        std::cout << "✓ Propagatore configurato:\n";
        std::cout << "  Integratore:        RKF78 (adattivo)\n";
        std::cout << "  Step size:          " << opts.stepSize << " giorni\n";
        std::cout << "  Tolleranza:         " << opts.tolerance << "\n";
        std::cout << "  Perturbazioni:      Attive\n";
        std::cout << "  Correzioni relativ: Attive\n\n";
        
        // ─────────────────────────────────────────────────────────────────────
        // STEP 3: Setup JPL Horizons
        // ─────────────────────────────────────────────────────────────────────
        std::cout << "STEP 3: Connessione a JPL Horizons\n";
        printSubSeparator();
        
        JPLHorizonsClient horizons;
        horizons.setTimeout(60);
        
        if (!horizons.isTargetAvailable("17030")) {
            std::cerr << "❌ Asteroide 17030 non disponibile su JPL Horizons\n";
            return 1;
        }
        std::cout << "✓ Asteroide 17030 disponibile\n\n";
        
        // ─────────────────────────────────────────────────────────────────────
        // STEP 4: Scarica dati JPL per l'epoca iniziale
        // ─────────────────────────────────────────────────────────────────────
        std::cout << "STEP 4: Scaricamento stato iniziale JPL\n";
        printSubSeparator();
        
        // Usa l'epoca degli elementi AstDyS come riferimento
        JulianDate refEpoch = astElements.epoch;
        std::cout << "✓ Epoca riferimento: JD " << std::fixed << std::setprecision(2) 
                  << refEpoch.jd << "\n\n";
        
        auto [pos0_jpl, vel0_jpl] = horizons.getStateVectors("17030", refEpoch, "@sun");
        std::cout << "  Pos JPL (AU): (" << std::setprecision(10)
                  << pos0_jpl.x << ", " << pos0_jpl.y << ", " << pos0_jpl.z << ")\n";
        std::cout << "  Vel JPL (AU/d): (" << vel0_jpl.x << ", " << vel0_jpl.y << ", "
                  << vel0_jpl.z << ")\n\n";
        
        // ─────────────────────────────────────────────────────────────────────
        // STEP 5: Propaga e confronta su intervallo 26/11 - 30/11/2025
        // ─────────────────────────────────────────────────────────────────────
        std::cout << "STEP 5: Confronto posizioni (26/11/2025 - 30/11/2025)\n";
        printSubSeparator();
        
        // Converte date in JD
        // 26/11/2025 00:00 UTC
        JulianDate startDate;
        startDate.jd = 2461004.5;  // 26 Nov 2025
        
        // 30/11/2025 23:59 UTC
        JulianDate endDate;
        endDate.jd = 2461008.5;    // 30 Nov 2025
        
        std::vector<EphemerisPoint> comparisonPoints;
        
        // Test ogni 12 ore
        double step = 0.5;  // 12 ore
        int pointCount = 0;
        
        std::cout << "Propagazione e confronto... ";
        std::flush(std::cout);
        
        for (double jd = startDate.jd; jd <= endDate.jd; jd += step) {
            JulianDate epoch;
            epoch.jd = jd;
            
            // Propaga con IOccultCalc
            auto stateIOC = propagator.propagate(astElements, epoch);
            
            // Scarica da JPL
            auto [posJPL, velJPL] = horizons.getStateVectors("17030", epoch, "@sun");
            
            // Calcola differenze
            Vector3D deltaPos = stateIOC.position - posJPL;
            Vector3D deltaVel = stateIOC.velocity - velJPL;
            
            double deltaPos_km = deltaPos.magnitude() * 149597870.7;  // AU -> km
            double deltaVel_mms = deltaVel.magnitude() * 149597870.7 / 86400.0 * 1000.0;  // AU/day -> mm/s
            
            EphemerisPoint point;
            point.epoch = epoch;
            point.posIOC = stateIOC.position;
            point.velIOC = stateIOC.velocity;
            point.posJPL = posJPL;
            point.velJPL = velJPL;
            point.deltaPos_km = deltaPos_km;
            point.deltaVel_mms = deltaVel_mms;
            
            comparisonPoints.push_back(point);
            pointCount++;
        }
        
        std::cout << "✓ " << pointCount << " punti calcolati\n\n";
        
        // ─────────────────────────────────────────────────────────────────────
        // STEP 6: Visualizza risultati dettagliati
        // ─────────────────────────────────────────────────────────────────────
        std::cout << "STEP 6: Risultati dettagliati\n";
        printSubSeparator();
        
        std::cout << std::setw(20) << "Data/Ora"
                  << std::setw(15) << "Err Pos (km)"
                  << std::setw(15) << "Err Pos (m)"
                  << std::setw(15) << "Err Vel (mm/s)"
                  << std::setw(15) << "Err RA (mas)"
                  << std::setw(15) << "Err Dec (mas)\n";
        printSubSeparator();
        
        for (const auto& p : comparisonPoints) {
            // Converti JD a data calendario
            int year, month, day, hour, minute;
            double second;
            TimeUtils::jdToCalendar(p.epoch, year, month, day, hour, minute, second);
            
            // Calcola coordinate sferiche
            double raIOC = atan2(p.posIOC.y, p.posIOC.x);
            double decIOC = asin(p.posIOC.z / p.posIOC.magnitude());
            
            double raJPL = atan2(p.posJPL.y, p.posJPL.x);
            double decJPL = asin(p.posJPL.z / p.posJPL.magnitude());
            
            double deltaRA_mas = (raIOC - raJPL) * 206264.806247;     // radianti -> milliarcsec
            double deltaDec_mas = (decIOC - decJPL) * 206264.806247;
            
            // Output
            char dateStr[32];
            sprintf(dateStr, "%04d-%02d-%02d %02d:%02d", year, month, day, hour, minute);
            
            std::cout << std::setw(20) << dateStr
                      << std::setw(15) << std::fixed << std::setprecision(6) << p.deltaPos_km
                      << std::setw(15) << std::setprecision(2) << p.deltaPos_km * 1000.0
                      << std::setw(15) << std::setprecision(3) << p.deltaVel_mms
                      << std::setw(15) << std::setprecision(2) << deltaRA_mas
                      << std::setw(15) << std::setprecision(2) << deltaDec_mas << "\n";
        }
        
        std::cout << "\n";
        
        // ─────────────────────────────────────────────────────────────────────
        // STEP 7: Statistiche
        // ─────────────────────────────────────────────────────────────────────
        std::cout << "STEP 7: Statistiche di errore\n";
        printSubSeparator();
        
        ComparisonStats stats = calculateStats(comparisonPoints);
        
        std::cout << "Errore di Posizione:\n";
        std::cout << "  Media:  " << std::fixed << std::setprecision(6) 
                  << stats.meanPosError_km << " km = "
                  << stats.meanPosError_km * 1000.0 << " m\n";
        std::cout << "  Max:    " << stats.maxPosError_km << " km = "
                  << stats.maxPosError_km * 1000.0 << " m\n";
        std::cout << "  Min:    " << stats.minPosError_km << " km = "
                  << stats.minPosError_km * 1000.0 << " m\n\n";
        
        std::cout << "Errore di Velocità:\n";
        std::cout << "  Media:  " << std::setprecision(4) 
                  << stats.meanVelError_mms << " mm/s\n";
        std::cout << "  Max:    " << stats.maxVelError_mms << " mm/s\n\n";
        
        // ─────────────────────────────────────────────────────────────────────
        // STEP 8: Valutazione accuratezza
        // ─────────────────────────────────────────────────────────────────────
        std::cout << "STEP 8: Valutazione accuratezza\n";
        printSubSeparator();
        
        if (stats.maxPosError_km < 0.001) {
            std::cout << "✓✓✓ ECCEZIONALE: Errore < 1 m (precisione sub-metrica)\n";
        } else if (stats.maxPosError_km < 0.1) {
            std::cout << "✓✓✓ ECCELLENTE: Errore < 100 m\n";
        } else if (stats.maxPosError_km < 1.0) {
            std::cout << "✓✓  OTTIMO: Errore < 1 km\n";
        } else if (stats.maxPosError_km < 100.0) {
            std::cout << "✓   BUONO: Errore < 100 km\n";
        } else if (stats.maxPosError_km < 1000.0) {
            std::cout << "⚠   ACCETTABILE: Errore < 1,000 km\n";
        } else {
            std::cout << "✗   ALTO: Errore > 1,000 km - verificare parametri\n";
        }
        
        std::cout << "\n";
        
        // ─────────────────────────────────────────────────────────────────────
        // STEP 9: Esporta risultati
        // ─────────────────────────────────────────────────────────────────────
        std::cout << "STEP 9: Esportazione risultati\n";
        printSubSeparator();
        
        std::ofstream outFile("comparison_17030_ioc_vs_horizons.csv");
        outFile << "JD,RA_IOC(deg),Dec_IOC(deg),RA_JPL(deg),Dec_JPL(deg),"
                << "DeltaRA(mas),DeltaDec(mas),DeltaPos(km),DeltaVel(mm/s)\n";
        
        for (const auto& p : comparisonPoints) {
            double raIOC = atan2(p.posIOC.y, p.posIOC.x) * 180.0 / M_PI;
            double decIOC = asin(p.posIOC.z / p.posIOC.magnitude()) * 180.0 / M_PI;
            
            double raJPL = atan2(p.posJPL.y, p.posJPL.x) * 180.0 / M_PI;
            double decJPL = asin(p.posJPL.z / p.posJPL.magnitude()) * 180.0 / M_PI;
            
            double deltaRA_mas = (raIOC - raJPL) * 3600.0 * 1000.0;
            double deltaDec_mas = (decIOC - decJPL) * 3600.0 * 1000.0;
            
            outFile << std::fixed << std::setprecision(6) << p.epoch.jd << ","
                    << raIOC << "," << decIOC << ","
                    << raJPL << "," << decJPL << ","
                    << deltaRA_mas << "," << deltaDec_mas << ","
                    << p.deltaPos_km << "," << p.deltaVel_mms << "\n";
        }
        outFile.close();
        
        std::cout << "✓ Risultati esportati in comparison_17030_ioc_vs_horizons.csv\n\n";
        
        printSeparator();
        std::cout << "✓ TEST COMPLETATO CON SUCCESSO\n";
        printSeparator();
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
