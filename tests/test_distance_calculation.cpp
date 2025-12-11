/**
 * @file test_distance_calculation.cpp
 * @brief Test calcolo distanza asteroide-stella per verificare errore
 */

#include "ioccultcalc/astdyn_propagation_helper.h"
#include "ioccultcalc/astdys_client.h"
#include "ioc_gaialib/unified_gaia_catalog.h"
#include "ioccultcalc/time_utils.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <filesystem>

using namespace ioccultcalc;
using namespace ioc::gaia;

constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;

double angularSeparation(double ra1_deg, double dec1_deg, double ra2_deg, double dec2_deg) {
    double ra1_rad = ra1_deg * DEG_TO_RAD;
    double dec1_rad = dec1_deg * DEG_TO_RAD;
    double ra2_rad = ra2_deg * DEG_TO_RAD;
    double dec2_rad = dec2_deg * DEG_TO_RAD;
    
    double delta_ra = ra2_rad - ra1_rad;
    double delta_dec = dec2_rad - dec1_rad;
    
    double a = std::sin(delta_dec / 2) * std::sin(delta_dec / 2) +
               std::cos(dec1_rad) * std::cos(dec2_rad) * 
               std::sin(delta_ra / 2) * std::sin(delta_ra / 2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
    
    return c * RAD_TO_DEG * 3600.0;  // arcsec
}

int main() {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TEST CALCOLO DISTANZA ASTEROIDE-STELLA                  ║\n";
    std::cout << "║  Asteroide 17030 - Stella GAIA:3411546266140512128        ║\n";
    std::cout << "║  28 Novembre 2025 06:33 UTC                               ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    // Inizializza catalogo
    std::string home = std::getenv("HOME") ? std::getenv("HOME") : "";
    std::string catalog_path = home + "/.catalog/gaia_mag18_v2_multifile";
    
    if (std::filesystem::exists(catalog_path) && std::filesystem::is_symlink(catalog_path)) {
        catalog_path = std::filesystem::canonical(catalog_path);
    }
    
    std::string json_config = R"({
        "catalog_type": "multifile_v2",
        "multifile_directory": ")" + catalog_path + R"(",
        "max_cached_chunks": 50,
        "log_level": "warning"
    })";
    
    if (!UnifiedGaiaCatalog::initialize(json_config)) {
        std::cerr << "❌ Failed to initialize UnifiedGaiaCatalog\n";
        return 1;
    }
    
    auto& catalog = UnifiedGaiaCatalog::getInstance();
    
    // Carica elementi asteroide
    std::cout << "Caricamento asteroide 17030...\n";
    AstDySElements elements = AstDySClient::downloadElements(17030);
    std::cout << "✓ Elementi caricati\n";
    std::cout << "  Epoch: MJD " << std::fixed << std::setprecision(2) << elements.epoch_mjd << "\n\n";
    
    // Data evento: 28 Nov 2025 06:33:10 UTC = MJD 61007.27303852
    double event_mjd = 61007.27303852;
    
    std::cout << "Data evento: MJD " << std::setprecision(8) << event_mjd << "\n";
    std::cout << "  (28 Novembre 2025 06:33:10 UTC)\n\n";
    
    // Carica stella
    uint64_t star_source_id = 3411546266140512128ULL;
    auto star_opt = catalog.queryBySourceId(star_source_id);
    
    if (!star_opt.has_value()) {
        std::cerr << "❌ Stella non trovata nel catalogo\n";
        return 1;
    }
    
    const auto& star = star_opt.value();
    std::cout << "Stella GAIA:" << star_source_id << ":\n";
    std::cout << "  RA: " << std::setprecision(10) << star.ra << " deg\n";
    std::cout << "  Dec: " << star.dec << " deg\n";
    std::cout << "  G mag: " << std::setprecision(2) << star.phot_g_mean_mag << "\n\n";
    
    // Calcola posizione asteroide
    std::cout << "Calcolo posizione asteroide...\n";
    auto& helper = AstDynPropagationHelper::getInstance();
    auto asteroid_radec = helper.getRADec(elements, event_mjd);
    
    std::cout << "Posizione asteroide:\n";
    std::cout << "  RA: " << std::setprecision(10) << asteroid_radec.first << " deg\n";
    std::cout << "  Dec: " << asteroid_radec.second << " deg\n\n";
    
    // Calcola distanza angolare
    double separation_arcsec = angularSeparation(
        star.ra, star.dec,
        asteroid_radec.first, asteroid_radec.second
    );
    
    std::cout << "Distanza angolare:\n";
    std::cout << "  Separation: " << std::setprecision(3) << separation_arcsec << " arcsec\n";
    std::cout << "  Separation: " << separation_arcsec / 60.0 << " arcmin\n";
    std::cout << "  Separation: " << separation_arcsec / 3600.0 << " deg\n\n";
    
    // Confronto con Horizons
    std::cout << "Confronto con Horizons (28 Nov 2025 06:33 UTC):\n";
    std::cout << "  Horizons RA: 73.316067 deg (04h 53m 15.856s)\n";
    std::cout << "  Horizons Dec: 20.325164 deg (+20° 19' 30.59\")\n";
    std::cout << "  Calcolo RA: " << std::setprecision(10) << asteroid_radec.first << " deg\n";
    std::cout << "  Calcolo Dec: " << asteroid_radec.second << " deg\n";
    
    double ra_diff = std::abs(asteroid_radec.first - 73.316067);
    if (ra_diff > 180.0) ra_diff = 360.0 - ra_diff;
    double dec_diff = std::abs(asteroid_radec.second - 20.325164);
    
    std::cout << "\n  Differenza RA: " << ra_diff << " deg = " 
              << ra_diff * 3600.0 << " arcsec\n";
    std::cout << "  Differenza Dec: " << dec_diff << " deg = " 
              << dec_diff * 3600.0 << " arcsec\n\n";
    
    // Verifica se la distanza è corretta
    std::cout << "Verifica distanza:\n";
    std::cout << "  Se la posizione asteroide è corretta, la distanza dovrebbe essere:\n";
    std::cout << "    Molto piccola (< 1 arcsec) se c'è occultazione\n";
    std::cout << "    Grande (> 1000 arcsec) se non c'è occultazione\n";
    std::cout << "\n  Distanza calcolata: " << separation_arcsec << " arcsec\n";
    
    if (separation_arcsec < 1.0) {
        std::cout << "  ✓ Distanza molto piccola - possibile occultazione\n";
    } else if (separation_arcsec < 100.0) {
        std::cout << "  ⚠ Distanza piccola ma non zero\n";
    } else {
        std::cout << "  ⚠ Distanza grande - verifica posizione asteroide\n";
    }
    
    // Calcola anche la posizione della stella rispetto a Horizons
    std::cout << "\nPosizione stella (da Gaia):\n";
    std::cout << "  RA: " << std::setprecision(10) << star.ra << " deg\n";
    std::cout << "  Dec: " << star.dec << " deg\n";
    
    // Se conosciamo la posizione della stella da altre fonti, confrontiamo
    // Per ora assumiamo che la stella sia corretta
    
    return 0;
}

