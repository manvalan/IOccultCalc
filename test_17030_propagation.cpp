#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/time_utils.h"

using namespace ioccultcalc;

int main() {
    std::cout << "═══════════════════════════════════════════════════════════\n";
    std::cout << "TEST PROPAGAZIONE (17030) Sierks\n";
    std::cout << "═══════════════════════════════════════════════════════════\n\n";
    
    // ================================================================
    // STEP 1: Leggi elementi AstDyS dal file
    // ================================================================
    std::cout << "STEP 1: Caricamento elementi da 17030_astdys.eq1\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    
    EquinoctialElements astElem;
    std::ifstream eq1File("17030_astdys.eq1");
    
    if (!eq1File.good()) {
        std::cerr << "❌ File 17030_astdys.eq1 non trovato!\n";
        return 1;
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
                
                std::cout << "Elementi Equinoziali letti:\n";
                std::cout << "  a      = " << std::fixed << std::setprecision(8) << a << " AU\n";
                std::cout << "  h      = " << h << "\n";
                std::cout << "  k      = " << k << "\n";
                std::cout << "  p      = " << p << "\n";
                std::cout << "  q      = " << q << "\n";
                std::cout << "  lambda = " << lambda << " rad\n";
            }
        } else if (found && line.find(" MJD ") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            double mjd;
            if (iss >> tag >> mjd) {
                astElem.epoch.jd = mjd + 2400000.5;
                std::cout << "  MJD    = " << std::fixed << std::setprecision(1) << mjd << "\n";
                std::cout << "  JD     = " << astElem.epoch.jd << "\n";
                
                // Converti in data calendario
                int year, month, day, hour, minute;
                double second;
                TimeUtils::jdToCalendar(astElem.epoch, year, month, day, hour, minute, second);
                std::cout << "  Data   = " << year << "-" << std::setw(2) << std::setfill('0') 
                          << month << "-" << day << "\n\n";
                break;
            }
        }
    }
    eq1File.close();
    
    if (!found) {
        std::cerr << "❌ Elementi non trovati nel file!\n";
        return 1;
    }
    
    // Converti in Keplerian per visualizzazione
    auto kepEpochOriginale = astElem.toKeplerian();
    std::cout << "Elementi Keplerian (Epoca originale):\n";
    std::cout << "  Epoca = JD " << std::fixed << std::setprecision(1) << kepEpochOriginale.epoch.jd << "\n";
    std::cout << "  a     = " << std::setprecision(8) << kepEpochOriginale.a << " AU\n";
    std::cout << "  e     = " << std::setprecision(6) << kepEpochOriginale.e << "\n";
    std::cout << "  i     = " << std::setprecision(4) << (kepEpochOriginale.i * RAD_TO_DEG) << "°\n";
    std::cout << "  Ω     = " << (kepEpochOriginale.Omega * RAD_TO_DEG) << "°\n";
    std::cout << "  ω     = " << (kepEpochOriginale.omega * RAD_TO_DEG) << "°\n";
    std::cout << "  M     = " << (kepEpochOriginale.M * RAD_TO_DEG) << "°\n\n";
    
    // ================================================================
    // STEP 2: Crea Ephemeris e propaga a epoca target
    // ================================================================
    std::cout << "STEP 2: Propagazione con RKF78\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    
    // Epoca target: 28 Nov 2025 00:00 UTC
    JulianDate targetEpoch;
    targetEpoch.jd = 2461007.5;  // 28 Nov 2025
    
    int targetYear, targetMonth, targetDay, targetHour, targetMinute;
    double targetSecond;
    TimeUtils::jdToCalendar(targetEpoch, targetYear, targetMonth, targetDay, targetHour, targetMinute, targetSecond);
    
    std::cout << "Epoca target = JD " << targetEpoch.jd << " (" 
              << targetYear << "-" << std::setw(2) << std::setfill('0') << targetMonth 
              << "-" << targetDay << ")\n";
    
    double deltaT = targetEpoch.jd - astElem.epoch.jd;
    std::cout << "Delta tempo  = " << std::fixed << std::setprecision(1) << deltaT << " giorni\n";
    std::cout << "             = " << std::setprecision(2) << (deltaT / 365.25) << " anni\n\n";
    
    // Crea Ephemeris
    Ephemeris ephemeris(astElem);
    ephemeris.enableNumericalPropagation(true, IntegratorType::RKF78);
    
    PropagatorOptions propOpts;
    propOpts.integrator = IntegratorType::RKF78;
    propOpts.stepSize = 0.01;
    propOpts.tolerance = 1e-12;
    propOpts.usePlanetaryPerturbations = true;
    ephemeris.setPropagatorOptions(propOpts);
    
    std::cout << "Configurazione propagatore:\n";
    std::cout << "  Integratore: RKF78\n";
    std::cout << "  Step size  : 0.01 giorni\n";
    std::cout << "  Tolleranza : 1e-12\n";
    std::cout << "  Perturbazioni planetarie: SÌ\n\n";
    
    // Propaga
    std::cout << "Propagazione in corso...\n";
    EphemerisData ephData = ephemeris.compute(targetEpoch);
    
    std::cout << "✓ Propagazione completata!\n\n";
    
    // ================================================================
    // STEP 3: Risultati
    // ================================================================
    std::cout << "STEP 3: Risultati propagazione\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    
    std::cout << "Posizione Geocentrica (Equatoriale J2000):\n";
    std::cout << "  RA  = " << std::fixed << std::setprecision(6) 
              << (ephData.geocentricPos.ra * RAD_TO_DEG) << "°\n";
    std::cout << "  Dec = " << (ephData.geocentricPos.dec * RAD_TO_DEG) << "°\n";
    std::cout << "  Distanza = " << std::setprecision(8) << ephData.distance << " AU\n\n";
    
    std::cout << "Confronto con stella target Gaia DR3 3411546266140512128:\n";
    std::cout << "  Stella RA  = 73.416106°\n";
    std::cout << "  Stella Dec = 20.331662°\n";
    double deltaRA = ephData.geocentricPos.ra * RAD_TO_DEG - 73.416106;
    double deltaDec = ephData.geocentricPos.dec * RAD_TO_DEG - 20.331662;
    std::cout << "  ΔRA  = " << std::setprecision(2) << deltaRA << "°\n";
    std::cout << "  ΔDec = " << deltaDec << "°\n";
    std::cout << "  Separazione ~ " << std::setprecision(1) 
              << sqrt(deltaRA*deltaRA + deltaDec*deltaDec) * 3600.0 << " arcsec\n\n";
    
    // ================================================================
    // STEP 4: Test propagazione in vari momenti
    // ================================================================
    std::cout << "STEP 4: Tabella propagazione 25 Nov - 2 Dic 2025\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    std::cout << "Data/Ora         | RA (deg)  | Dec (deg) | Delta dalla stella (deg)\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    
    for (double jd = 2461004.5; jd <= 2461011.5; jd += 1.0) {
        JulianDate epoch;
        epoch.jd = jd;
        
        EphemerisData eph = ephemeris.compute(epoch);
        
        int y, m, d, h, min;
        double s;
        TimeUtils::jdToCalendar(epoch, y, m, d, h, min, s);
        
        double ra_deg = eph.geocentricPos.ra * RAD_TO_DEG;
        double dec_deg = eph.geocentricPos.dec * RAD_TO_DEG;
        double sep_deg = sqrt(pow(ra_deg - 73.416106, 2) + pow(dec_deg - 20.331662, 2));
        
        std::cout << y << "-" << std::setw(2) << std::setfill('0') << m << "-" << d << " 00:00"
                  << " | " << std::setw(9) << std::fixed << std::setprecision(3) << ra_deg
                  << " | " << std::setw(9) << dec_deg
                  << " | " << std::setw(8) << std::setprecision(2) << sep_deg << "\n";
    }
    
    std::cout << "\n═══════════════════════════════════════════════════════════\n";
    std::cout << "TEST COMPLETATO\n";
    std::cout << "═══════════════════════════════════════════════════════════\n";
    
    return 0;
}
