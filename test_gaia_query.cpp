#include <iostream>
#include <iomanip>
#include "ioccultcalc/gaia_adapter.h"

using namespace ioccultcalc;

int main() {
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "   TEST QUERY GAIA - Stella UCAC4 552-011427\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";
    
    // Coordinate UCAC4 della stella target
    double star_ra = 73.4161092;   // gradi
    double star_dec = 20.3316578;  // gradi
    double radius = 0.1;           // gradi (6 arcmin)
    double max_mag = 14.0;
    
    std::cout << "Query parametri:\n";
    std::cout << "  RA center:  " << star_ra << "°\n";
    std::cout << "  Dec center: " << star_dec << "°\n";
    std::cout << "  Radius:     " << radius << "° = " << radius*60 << " arcmin\n";
    std::cout << "  Max mag:    " << max_mag << "\n\n";
    
    // Inizializza client Gaia
    GaiaClient client;
    
    std::cout << "Esecuzione query...\n";
    auto stars = client.queryCone(star_ra, star_dec, radius, max_mag);
    
    std::cout << "\nRisultato: " << stars.size() << " stelle trovate\n";
    std::cout << "─────────────────────────────────────────────────────────────────\n";
    
    for (size_t i = 0; i < std::min(stars.size(), (size_t)10); i++) {
        const auto& s = stars[i];
        
        // Calcola separazione dalla posizione target
        double delta_ra = (s.pos.ra - star_ra) * cos(star_dec * M_PI / 180.0);
        double delta_dec = s.pos.dec - star_dec;
        double sep_arcsec = sqrt(delta_ra*delta_ra + delta_dec*delta_dec) * 3600;
        
        std::cout << "\n" << (i+1) << ". Gaia DR3 " << s.sourceId << "\n";
        std::cout << "   RA  = " << std::fixed << std::setprecision(7) << s.pos.ra << "°\n";
        std::cout << "   Dec = " << std::showpos << s.pos.dec << std::noshowpos << "°\n";
        std::cout << "   G mag = " << std::setprecision(2) << s.gMag << "\n";
        std::cout << "   Sep = " << sep_arcsec << " arcsec\n";
    }
    
    std::cout << "\n═══════════════════════════════════════════════════════════════\n";
    std::cout << "STELLA TARGET (UCAC4 552-011427):\n";
    std::cout << "  RA J2000:  73.4161092°\n";
    std::cout << "  Dec J2000: +20.3316578°\n";
    std::cout << "  V mag:     13.111\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    
    return 0;
}
