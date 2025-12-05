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
    std::cout << "TEST PROPAGAZIONE (17030) Sierks - ELEMENTI ORBFIT\n";
    std::cout << "═══════════════════════════════════════════════════════════\n\n";
    
    // ================================================================
    // STEP 1: Leggi elementi OrbFit corretti dal file
    // ================================================================
    std::cout << "STEP 1: Caricamento elementi da 17030_orbfit.oel\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    
    EquinoctialElements orbfitElem;
    std::ifstream oelFile("17030_orbfit.oel");
    
    if (!oelFile.good()) {
        std::cerr << "❌ File 17030_orbfit.oel non trovato!\n";
        return 1;
    }
    
    std::string line;
    bool found = false;
    while (std::getline(oelFile, line)) {
        if (line.find(" EQU ") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            double a, h, k, p, q, lambda;
            if (iss >> tag >> a >> h >> k >> p >> q >> lambda) {
                orbfitElem.a = a;
                orbfitElem.h = h;
                orbfitElem.k = k;
                orbfitElem.p = p;
                orbfitElem.q = q;
                // Nel file OrbFit lambda è in GRADI, convertiamo in radianti
                orbfitElem.lambda = lambda * M_PI / 180.0;
                found = true;
                
                std::cout << "Elementi Equinoziali OrbFit letti:\n";
                std::cout << "  a      = " << std::fixed << std::setprecision(12) << a << " AU\n";
                std::cout << "  h      = " << std::setprecision(12) << h << "\n";
                std::cout << "  k      = " << k << "\n";
                std::cout << "  p      = " << p << "\n";
                std::cout << "  q      = " << q << "\n";
                std::cout << "  lambda = " << lambda << "° = " << orbfitElem.lambda << " rad\n";
            }
        } else if (found && line.find(" MJD ") != std::string::npos) {
            std::istringstream iss(line);
            std::string tag;
            double mjd;
            if (iss >> tag >> mjd) {
                orbfitElem.epoch.jd = mjd + 2400000.5;
                std::cout << "\n  Epoca: MJD " << std::setprecision(6) << mjd 
                          << " = JD " << std::setprecision(6) << orbfitElem.epoch.jd << "\n";
                std::cout << "         20 Novembre 2025, 00:00 UTC\n";
                break;
            }
        }
    }
    
    if (!found) {
        std::cerr << "❌ Elementi non trovati nel file!\n";
        return 1;
    }
    
    // Converti in elementi Kepleriani per visualizzazione
    double e = std::sqrt(orbfitElem.h * orbfitElem.h + orbfitElem.k * orbfitElem.k);
    double i_rad = 2.0 * std::atan(std::sqrt(orbfitElem.p * orbfitElem.p + 
                                              orbfitElem.q * orbfitElem.q));
    double i_deg = i_rad * 180.0 / M_PI;
    
    std::cout << "\nElementi Kepleriani equivalenti:\n";
    std::cout << "  e = " << std::setprecision(8) << e << "\n";
    std::cout << "  i = " << std::setprecision(6) << i_deg << "°\n";
    std::cout << "✓ Caricamento completato\n\n";
    
    // ================================================================
    // STEP 2: Setup Propagatore
    // ================================================================
    std::cout << "STEP 2: Setup Propagatore RKF78\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    
    Ephemeris ephemeris(orbfitElem);
    ephemeris.enableNumericalPropagation(true, IntegratorType::RKF78);
    
    PropagatorOptions opts;
    opts.integrator = IntegratorType::RKF78;
    opts.stepSize = 0.01;
    opts.tolerance = 1.0e-12;
    opts.usePlanetaryPerturbations = true;
    opts.useRelativisticCorrections = false;
    
    ephemeris.setPropagatorOptions(opts);
    
    std::cout << "  Integratore: RKF78\n";
    std::cout << "  Step size: " << opts.stepSize << " giorni\n";
    std::cout << "  Tolleranza: " << opts.tolerance << "\n";
    std::cout << "  Perturbazioni planetarie: " << (opts.usePlanetaryPerturbations ? "SÌ" : "NO") << "\n";
    std::cout << "✓ Propagatore configurato\n\n";
    
    // ================================================================
    // STEP 3: Propagazione a date multiple
    // ================================================================
    std::cout << "STEP 3: Propagazione Nov 20-28, 2025\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    std::cout << "Data              JD            RA (°)      Dec (°)    Delta (AU)\n";
    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
    
    // Date da Nov 20 a Nov 28 (ogni giorno)
    struct TestDate {
        int year, month, day;
        const char* name;
    };
    
    TestDate dates[] = {
        {2025, 11, 20, "Nov 20 (epoca)"},
        {2025, 11, 21, "Nov 21"},
        {2025, 11, 22, "Nov 22"},
        {2025, 11, 23, "Nov 23"},
        {2025, 11, 24, "Nov 24"},
        {2025, 11, 25, "Nov 25"},
        {2025, 11, 26, "Nov 26"},
        {2025, 11, 27, "Nov 27"},
        {2025, 11, 28, "Nov 28 (TARGET!)"}
    };
    
    for (const auto& date : dates) {
        // Calcola JD per mezzanotte UTC
        JulianDate targetEpoch = TimeUtils::calendarToJD(date.year, date.month, date.day, 0, 0, 0.0);
        
        // Propaga
        EphemerisData ephData;
        try {
            ephData = ephemeris.compute(targetEpoch);
        } catch (const std::exception& e) {
            std::cerr << "❌ Errore propagazione: " << e.what() << "\n";
            continue;
        }
        
        // Coordinate già in RA/Dec
        double ra_rad = ephData.geocentricPos.ra;
        double dec_rad = ephData.geocentricPos.dec;
        double distance = ephData.distance;
        
        if (ra_rad < 0) ra_rad += 2.0 * M_PI;
        
        double ra_deg = ra_rad * 180.0 / M_PI;
        double dec_deg = dec_rad * 180.0 / M_PI;
        
        std::cout << std::setw(17) << std::left << date.name
                  << std::fixed << std::setprecision(1) << std::setw(13) << targetEpoch.jd
                  << std::setprecision(4) << std::setw(12) << ra_deg
                  << std::setw(11) << dec_deg
                  << std::setprecision(6) << distance << "\n";
    }
    
    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
    
    // ================================================================
    // STEP 4: Propagazione precisa a Nov 28, 2025 00:00 UTC
    // ================================================================
    std::cout << "\nSTEP 4: Posizione precisa il 28 Novembre 2025, 00:00 UTC\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    
    JulianDate targetEpoch = TimeUtils::calendarToJD(2025, 11, 28, 0, 0, 0.0);
    EphemerisData finalData = ephemeris.compute(targetEpoch);
    
    double ra_rad = finalData.geocentricPos.ra;
    double dec_rad = finalData.geocentricPos.dec;
    double distance = finalData.distance;
    if (ra_rad < 0) ra_rad += 2.0 * M_PI;
    
    double ra_deg = ra_rad * 180.0 / M_PI;
    double dec_deg = dec_rad * 180.0 / M_PI;
    
    // Converti RA in ore
    double ra_hours = ra_deg / 15.0;
    int ra_h = (int)ra_hours;
    int ra_m = (int)((ra_hours - ra_h) * 60.0);
    double ra_s = ((ra_hours - ra_h) * 60.0 - ra_m) * 60.0;
    
    // Converti Dec in gradi/arcmin/arcsec
    int dec_sign = (dec_deg >= 0) ? 1 : -1;
    double dec_abs = std::abs(dec_deg);
    int dec_d = (int)dec_abs;
    int dec_m = (int)((dec_abs - dec_d) * 60.0);
    double dec_s = ((dec_abs - dec_d) * 60.0 - dec_m) * 60.0;
    
    std::cout << "JD       = " << std::fixed << std::setprecision(6) << targetEpoch.jd << "\n";
    std::cout << "RA       = " << std::setprecision(6) << ra_deg << "° = "
              << std::setw(2) << std::setfill('0') << ra_h << "h "
              << std::setw(2) << ra_m << "m "
              << std::setprecision(3) << ra_s << "s\n" << std::setfill(' ');
    std::cout << "Dec      = " << std::setprecision(6) << dec_deg << "° = "
              << (dec_sign > 0 ? "+" : "-")
              << std::setw(2) << std::setfill('0') << dec_d << "° "
              << std::setw(2) << dec_m << "' "
              << std::setprecision(2) << dec_s << "\"\n" << std::setfill(' ');
    std::cout << "Distance = " << std::setprecision(8) << distance << " AU\n";
    
    // ================================================================
    // STEP 5: Confronto con stella target
    // ================================================================
    std::cout << "\nSTEP 5: Confronto con Gaia DR3 3411546266140512128\n";
    std::cout << "─────────────────────────────────────────────────────────\n";
    
    // Coordinate stella da IOTA professional
    double star_ra_deg = 73.177317;   // 04h 52m 42.56s
    double star_dec_deg = 20.325531;  // +20° 19' 31.9"
    
    std::cout << "Stella:    RA = " << std::fixed << std::setprecision(6) << star_ra_deg 
              << "°, Dec = " << star_dec_deg << "°\n";
    std::cout << "Asteroide: RA = " << ra_deg << "°, Dec = " << dec_deg << "°\n";
    
    // Calcola separazione angolare
    double dRA = (ra_deg - star_ra_deg) * std::cos(dec_rad);
    double dDec = dec_deg - star_dec_deg;
    double separation_deg = std::sqrt(dRA*dRA + dDec*dDec);
    double separation_arcsec = separation_deg * 3600.0;
    
    std::cout << "\nSeparazione angolare:\n";
    std::cout << "  ΔRA  = " << std::setprecision(6) << (ra_deg - star_ra_deg) 
              << "° = " << (ra_deg - star_ra_deg)*3600.0 << " arcsec\n";
    std::cout << "  ΔDec = " << (dec_deg - star_dec_deg) 
              << "° = " << (dec_deg - star_dec_deg)*3600.0 << " arcsec\n";
    std::cout << "  Separazione totale = " << separation_deg << "° = " 
              << separation_arcsec << " arcsec\n";
    
    if (separation_arcsec < 2.0) {
        std::cout << "\n✓✓✓ OCCULTAZIONE POSSIBILE! (< 2 arcsec) ✓✓✓\n";
    } else if (separation_arcsec < 10.0) {
        std::cout << "\n⚠ Passaggio vicino (< 10 arcsec)\n";
    } else {
        std::cout << "\n✗ Separazione troppo grande per occultazione\n";
    }
    
    std::cout << "\n═══════════════════════════════════════════════════════════\n";
    std::cout << "TEST COMPLETATO\n";
    std::cout << "═══════════════════════════════════════════════════════════\n";
    
    return 0;
}
