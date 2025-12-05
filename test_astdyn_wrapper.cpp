/**
 * @file test_astdyn_wrapper.cpp
 * @brief Test del wrapper AstDynPropagator
 * 
 * Verifica che il wrapper funzioni correttamente per il caso 17030
 */

#include "ioccultcalc/astdyn_interface.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

// Calcola distanza angolare tra due punti
double angularDistance(double ra1_deg, double dec1_deg, 
                       double ra2_deg, double dec2_deg) {
    double ra1 = ra1_deg * M_PI / 180.0;
    double dec1 = dec1_deg * M_PI / 180.0;
    double ra2 = ra2_deg * M_PI / 180.0;
    double dec2 = dec2_deg * M_PI / 180.0;
    
    double cos_sep = sin(dec1) * sin(dec2) + 
                     cos(dec1) * cos(dec2) * cos(ra1 - ra2);
    cos_sep = std::max(-1.0, std::min(1.0, cos_sep));
    double sep_rad = acos(cos_sep);
    return sep_rad * 180.0 / M_PI * 3600.0;  // arcsec
}

int main() {
    std::cout << "================================================================\n";
    std::cout << " Test AstDynPropagator Wrapper\n";
    std::cout << " Asteroide 17030 - Occultazione 28 Nov 2025\n";
    std::cout << "================================================================\n\n";
    
    try {
        // Carica elementi OrbFit
        std::cout << "Carico elementi da 17030_orbfit.oel...\n";
        AstDySElements elem = AstDySElements::fromFile("17030_orbfit.oel");
        
        std::cout << "Elementi caricati:\n";
        std::cout << "  Nome:  " << elem.name << "\n";
        std::cout << "  Epoca: MJD " << std::fixed << std::setprecision(6) << elem.epoch_mjd << "\n";
        std::cout << "  a = " << elem.a << " AU\n";
        std::cout << "  e = " << elem.e << "\n";
        std::cout << "  i = " << elem.i << "°\n";
        std::cout << "  Ω = " << elem.Omega << "°\n";
        std::cout << "  ω = " << elem.omega << "°\n";
        std::cout << "  M = " << elem.M << "°\n";
        std::cout << "\n";
        
        // Crea propagatore
        std::cout << "Creo AstDynPropagator (tolleranza 1e-12)...\n";
        AstDynPropagator propagator(1e-12);
        
        // Stella target
        double star_ra = 73.416108;
        double star_dec = 20.331662;
        std::cout << "Stella GAIA DR3 3411546266140512128:\n";
        std::cout << "  RA  = " << std::fixed << std::setprecision(8) << star_ra << "°\n";
        std::cout << "  Dec = " << star_dec << "°\n";
        std::cout << "\n";
        
        // Test propagazione 28 Nov 2025, 00:35 UTC
        std::cout << "================================================================\n";
        std::cout << " Test 1: 28 Nov 2025, 00:35 UTC (momento minima distanza)\n";
        std::cout << "================================================================\n";
        
        // MJD per 28 Nov 2025, 00:35 UTC
        double test_mjd = 60607.0 + (35.0 / 1440.0);  // MJD base + 35 minuti
        double jd_base = 2400000.5;
        double test_jd = test_mjd + jd_base;
        
        std::cout << "Target: JD " << std::fixed << std::setprecision(6) << test_jd 
                  << " = MJD " << test_mjd << "\n";
        std::cout << "Δt da epoca = " << (test_mjd - elem.epoch_mjd) << " giorni\n\n";
        
        // Calcola posizione
        std::cout << "Calcolo posizione con AstDynPropagator...\n";
        auto [ra_deg, dec_deg] = propagator.getRADec(elem, test_mjd);
        
        std::cout << "\nRisultato:\n";
        std::cout << "  RA  = " << std::fixed << std::setprecision(8) << ra_deg << "°\n";
        std::cout << "  Dec = " << dec_deg << "°\n";
        
        // Calcola distanza da stella
        double sep_arcsec = angularDistance(ra_deg, dec_deg, star_ra, star_dec);
        std::cout << "  Separazione da stella = " << std::fixed << std::setprecision(2) 
                  << sep_arcsec << " arcsec\n";
        
        // Statistiche propagazione
        std::cout << "\nStatistiche propagazione:\n";
        std::cout << "  Passi accettati: " << propagator.getLastStepsAccepted() << "\n";
        std::cout << "  Passi rifiutati: " << propagator.getLastStepsRejected() << "\n";
        std::cout << "  Passo min: " << std::fixed << std::setprecision(6) 
                  << propagator.getLastMinStep() << " giorni\n";
        std::cout << "  Passo max: " << propagator.getLastMaxStep() << " giorni\n";
        
        // Confronto con AstDyn standalone
        std::cout << "\n================================================================\n";
        std::cout << " CONFRONTO con AstDyn standalone\n";
        std::cout << "================================================================\n";
        std::cout << "AstDyn standalone: RA = 73.416167°, Dec = 20.332083°, Sep = 1.53\"\n";
        std::cout << "Wrapper:           RA = " << std::fixed << std::setprecision(6) 
                  << ra_deg << "°, Dec = " << dec_deg << "°, Sep = " 
                  << std::setprecision(2) << sep_arcsec << "\"\n";
        
        double delta_ra = (ra_deg - 73.416167) * 3600.0;
        double delta_dec = (dec_deg - 20.332083) * 3600.0;
        std::cout << "\nDifferenza:\n";
        std::cout << "  ΔRA  = " << std::showpos << std::fixed << std::setprecision(2) 
                  << delta_ra << " arcsec\n";
        std::cout << "  ΔDec = " << delta_dec << " arcsec\n";
        
        if (std::abs(delta_ra) < 1.0 && std::abs(delta_dec) < 1.0) {
            std::cout << "\n✅ WRAPPER FUNZIONA CORRETTAMENTE!\n";
        } else {
            std::cout << "\n⚠️  DIFFERENZA SIGNIFICATIVA - Verificare!\n";
        }
        
        std::cout << "\n================================================================\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ ERRORE: " << e.what() << "\n";
        return 1;
    }
}
