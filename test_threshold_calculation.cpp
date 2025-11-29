/**
 * Test numerico del FIX threshold occultazione
 * Verifica che il nuovo threshold LinOccult sia ~5-10× più grande del vecchio
 */

#include <iostream>
#include <iomanip>
#include <cmath>

const double AU = 149597870.7; // km
const double RAD_TO_DEG = 57.29577951308232;
const double R_EARTH_KM = 6371.0;

// VECCHIO METODO (pre-fix)
double oldThreshold(double asteroidDiameterKm, double distanceAU, double uncertaintyKm) {
    double asteroidDistanceKm = distanceAU * AU;
    double asteroidAngularSize = (asteroidDiameterKm / asteroidDistanceKm) * RAD_TO_DEG * 3600.0;
    double threshold = asteroidAngularSize + 3.0 * (uncertaintyKm / asteroidDistanceKm) * RAD_TO_DEG * 3600.0;
    if (threshold == 0) threshold = 10.0;
    return threshold;
}

// NUOVO METODO (post-fix, stile LinOccult)
double newThreshold(double asteroidDiameterKm, double distanceAU, double uncertaintyKm) {
    double asteroidDistanceKm = distanceAU * AU;
    
    const double MAX_SHADOW_RADIUS_KM = 3.0 * R_EARTH_KM;  // ~19,000 km
    double maxDistKm = uncertaintyKm + R_EARTH_KM;
    
    if (maxDistKm > MAX_SHADOW_RADIUS_KM) {
        maxDistKm = MAX_SHADOW_RADIUS_KM;
    }
    
    double thresholdArcsec = (maxDistKm / asteroidDistanceKm) * RAD_TO_DEG * 3600.0;
    if (thresholdArcsec < 10.0) thresholdArcsec = 10.0;
    
    return thresholdArcsec;
}

void testAsteroid(const char* name, double diamKm, double distAU, double uncertKm) {
    double oldThresh = oldThreshold(diamKm, distAU, uncertKm);
    double newThresh = newThreshold(diamKm, distAU, uncertKm);
    double ratio = newThresh / oldThresh;
    
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\n" << name << ":\n";
    std::cout << "  Diametro: " << diamKm << " km | Distanza: " << distAU << " AU | Incertezza: " << uncertKm << " km\n";
    std::cout << "  Vecchio threshold: " << std::setw(8) << oldThresh << " arcsec\n";
    std::cout << "  Nuovo threshold:   " << std::setw(8) << newThresh << " arcsec\n";
    std::cout << "  Rapporto (new/old): " << std::setw(7) << ratio << "×";
    
    if (ratio >= 3.0) {
        std::cout << "  ✓ MIGLIORAMENTO SIGNIFICATIVO\n";
    } else if (ratio >= 1.5) {
        std::cout << "  ✓ Miglioramento moderato\n";
    } else {
        std::cout << "  - Piccolo miglioramento\n";
    }
}

int main() {
    std::cout << "═════════════════════════════════════════════════════════════\n";
    std::cout << "TEST NUMERICO: FIX THRESHOLD OCCULTAZIONE\n";
    std::cout << "Confronto vecchio metodo vs nuovo metodo LinOccult\n";
    std::cout << "═════════════════════════════════════════════════════════════\n";
    
    // Test asteroidi reali in configurazioni tipiche
    testAsteroid("(1) Ceres", 939, 2.8, 1000);
    testAsteroid("(2) Pallas", 512, 2.7, 800);
    testAsteroid("(4) Vesta", 525, 2.5, 500);
    testAsteroid("(10) Hygiea", 434, 3.2, 1200);
    testAsteroid("(15) Eunomia", 255, 2.7, 1500);
    testAsteroid("(52) Europa", 303, 3.1, 2000);
    
    // Test casi estremi
    std::cout << "\n═══ CASI ESTREMI ═══\n";
    testAsteroid("Asteroide piccolo vicino", 50, 1.5, 3000);
    testAsteroid("Asteroide grande lontano", 900, 4.0, 500);
    testAsteroid("Incertezza molto alta", 200, 2.5, 5000);
    
    std::cout << "\n═════════════════════════════════════════════════════════════\n";
    std::cout << "CONCLUSIONI:\n";
    std::cout << "- Il nuovo threshold è tipicamente 3-10× più grande\n";
    std::cout << "- Questo porterà a molti più eventi candidati\n";
    std::cout << "- Il filtro di probabilità eliminerà i falsi positivi\n";
    std::cout << "═════════════════════════════════════════════════════════════\n";
    
    return 0;
}
