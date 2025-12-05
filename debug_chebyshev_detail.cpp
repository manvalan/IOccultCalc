#include <iostream>
#include <iomanip>
#include <vector>
#include "include/ioccultcalc/astdyn_interface.h"
#include "include/ioccultcalc/propagation_strategy.h"

using namespace ioccultcalc;

int main() {
    std::cout << "🔍 DEBUG DETTAGLIATO: Analisi Chebyshev vs RKF78\n";
    std::cout << "==================================================\n\n";
    
    // Carica elementi reali di 17030
    AstDySElements elements = AstDySElements::fromFile("17030_astdys.eq1");
    
    std::cout << "📂 Elementi 17030 caricati (epoca MJD " << elements.epoch_mjd << ")\n";
    std::cout << "   a=" << elements.a << " AU, e=" << elements.e 
              << ", i=" << elements.i << "°\n\n";
    
    // Configura strategia
    PropagationConfig config = propagation_presets::createBalanced();
    config.verbose_timing = true;
    TwoPhaseStrategy strategy(config);
    strategy.setElements(elements);
    
    // Target: 28 Nov 2025 (MJD 60006.024306)
    double target_mjd = 60006.024306;
    
    std::cout << "🎯 TARGET: MJD " << target_mjd << " (28 Nov 2025)\n";
    std::cout << "⏰ Propagazione: " << (target_mjd - elements.epoch_mjd) << " giorni\n\n";
    
    // PARTE 1: RKF78 diretto (riferimento)
    std::cout << "🔴 PARTE 1: RKF78 DIRETTO (RIFERIMENTO)\n";
    std::cout << "-----------------------------------------\n";
    
    auto rkf78_pos = strategy.getRKF78Position(target_mjd);
    std::cout << "✅ RKF78: RA=" << std::fixed << std::setprecision(6) 
              << rkf78_pos.ra_deg << "°, Dec=" << rkf78_pos.dec_deg << "°\n\n";
    
    // PARTE 2: Analisi generazione Chebyshev
    std::cout << "🔵 PARTE 2: GENERAZIONE POLINOMI CHEBYSHEV\n";
    std::cout << "--------------------------------------------\n";
    
    // Finestra per Chebyshev: 10 giorni centrati su target
    double window_days = 10.0;
    double start_mjd = target_mjd - window_days/2.0;
    double end_mjd = target_mjd + window_days/2.0;
    
    std::cout << "📊 Finestra Chebyshev: [" << start_mjd << ", " << end_mjd << "] MJD\n";
    std::cout << "🔢 Gradi polinomi: 8\n";
    
    // Genera punti di controllo manualmente per debug
    int degree = 8;
    int n_points = degree + 5; // 13 punti
    std::cout << "📍 Punti controllo: " << n_points << "\n\n";
    
    // Inizializza AstDyn direttamente per controllo
    AstDynPropagator astdyn(1e-12);
    astdyn.usePlanetPerturbations(true);
    astdyn.useAsteroidPerturbations(false);
    astdyn.useRelativisticCorrections(false);
    
    std::cout << "🎯 PUNTI DI CONTROLLO GENERATI:\n";
    std::vector<double> times(n_points);
    std::vector<double> ra_values(n_points);
    std::vector<double> dec_values(n_points);
    
    for (int i = 0; i < n_points; i++) {
        double t = start_mjd + (end_mjd - start_mjd) * i / (n_points - 1);
        times[i] = t;
        
        auto radec = astdyn.getRADec(elements, t);
        ra_values[i] = radec.first;
        dec_values[i] = radec.second;
        
        std::cout << "  " << std::setw(2) << i+1 << ": MJD " << std::fixed << std::setprecision(3) 
                  << t << " → RA=" << std::setprecision(6) << ra_values[i] 
                  << "°, Dec=" << dec_values[i] << "°\n";
        
        // Evidenzia il punto target
        if (std::abs(t - target_mjd) < 0.1) {
            std::cout << "      ⭐ PUNTO TARGET\n";
        }
    }
    
    std::cout << "\n🧮 ANALISI VARIAZIONE RA/Dec:\n";
    double ra_min = *std::min_element(ra_values.begin(), ra_values.end());
    double ra_max = *std::max_element(ra_values.begin(), ra_values.end());
    double dec_min = *std::min_element(dec_values.begin(), dec_values.end());
    double dec_max = *std::max_element(dec_values.begin(), dec_values.end());
    
    std::cout << "  RA range: [" << ra_min << "°, " << ra_max << "°] = " 
              << (ra_max - ra_min) << "° variazione\n";
    std::cout << "  Dec range: [" << dec_min << "°, " << dec_max << "°] = " 
              << (dec_max - dec_min) << "° variazione\n\n";
    
    // PARTE 3: Calcolo coefficienti Chebyshev manuale
    std::cout << "⚙️ PARTE 3: CALCOLO COEFFICIENTI CHEBYSHEV\n";
    std::cout << "-------------------------------------------\n";
    
    // Implementazione semplice least squares (come nel codice)
    std::vector<double> ra_coeffs(degree + 1, 0.0);
    std::vector<double> dec_coeffs(degree + 1, 0.0);
    
    for (int j = 0; j <= degree; j++) {
        double ra_sum = 0.0;
        double dec_sum = 0.0;
        
        for (int i = 0; i < n_points; i++) {
            double t_norm = 2.0 * (times[i] - times[0]) / (times.back() - times[0]) - 1.0;
            
            // Funzione base Chebyshev T_j(x)
            double basis = 1.0; // T_0
            if (j == 1) basis = t_norm; // T_1
            else if (j > 1) {
                double T_prev = 1.0;
                double T_curr = t_norm;
                for (int k = 2; k <= j; k++) {
                    double T_next = 2.0 * t_norm * T_curr - T_prev;
                    T_prev = T_curr;
                    T_curr = T_next;
                }
                basis = T_curr;
            }
            
            ra_sum += ra_values[i] * basis;
            dec_sum += dec_values[i] * basis;
        }
        
        ra_coeffs[j] = ra_sum / n_points;
        dec_coeffs[j] = dec_sum / n_points;
        
        std::cout << "  C_" << j << ": RA=" << std::scientific << std::setprecision(3) 
                  << ra_coeffs[j] << ", Dec=" << dec_coeffs[j] << "\n";
    }
    
    // PARTE 4: Valutazione al target
    std::cout << "\n🎯 PARTE 4: VALUTAZIONE CHEBYSHEV AL TARGET\n";
    std::cout << "--------------------------------------------\n";
    
    // Normalizza tempo target
    double t_norm_target = 2.0 * (target_mjd - start_mjd) / (end_mjd - start_mjd) - 1.0;
    std::cout << "📐 Tempo normalizzato: " << t_norm_target << "\n";
    
    // Valuta polinomi
    double ra_cheb = ra_coeffs[0]; // T_0 = 1
    double dec_cheb = dec_coeffs[0];
    
    if (ra_coeffs.size() > 1) {
        double T_prev = 1.0;
        double T_curr = t_norm_target;
        ra_cheb += ra_coeffs[1] * T_curr;
        dec_cheb += dec_coeffs[1] * T_curr;
        
        std::cout << "  T_0=" << 1.0 << ", T_1=" << T_curr << "\n";
        
        for (size_t n = 2; n < ra_coeffs.size(); n++) {
            double T_next = 2.0 * t_norm_target * T_curr - T_prev;
            ra_cheb += ra_coeffs[n] * T_next;
            dec_cheb += dec_coeffs[n] * T_next;
            
            std::cout << "  T_" << n << "=" << T_next << "\n";
            
            T_prev = T_curr;
            T_curr = T_next;
        }
    }
    
    std::cout << "\n✅ RISULTATO CHEBYSHEV: RA=" << std::fixed << std::setprecision(6) 
              << ra_cheb << "°, Dec=" << dec_cheb << "°\n";
    
    // PARTE 5: Confronto finale
    std::cout << "\n🔍 PARTE 5: CONFRONTO E ANALISI ERRORI\n";
    std::cout << "---------------------------------------\n";
    
    double ra_diff = ra_cheb - rkf78_pos.ra_deg;
    double dec_diff = dec_cheb - rkf78_pos.dec_deg;
    double error_arcsec = std::sqrt(ra_diff*ra_diff + dec_diff*dec_diff) * 3600.0;
    
    std::cout << "📊 CONFRONTO:\n";
    std::cout << "  RKF78:     RA=" << std::setprecision(6) << rkf78_pos.ra_deg 
              << "°, Dec=" << rkf78_pos.dec_deg << "°\n";
    std::cout << "  Chebyshev: RA=" << ra_cheb << "°, Dec=" << dec_cheb << "°\n";
    std::cout << "  Differenze: ΔRA=" << ra_diff << "°, ΔDec=" << dec_diff << "°\n";
    std::cout << "  Errore totale: " << std::setprecision(1) << error_arcsec << " arcsec\n\n";
    
    if (error_arcsec > 1000.0) {
        std::cout << "❌ ERRORE ECCESSIVO! Possibili cause:\n";
        std::cout << "   1. Algoritmo calcolo coefficienti errato\n";
        std::cout << "   2. Normalizzazione tempo sbagliata\n";
        std::cout << "   3. Valutazione polinomi incorretta\n";
        std::cout << "   4. Range RA vicino a 0°/360° (discontinuità)\n";
    } else if (error_arcsec > 10.0) {
        std::cout << "⚠️ ERRORE SIGNIFICATIVO - Controllare algoritmo\n";
    } else {
        std::cout << "✅ ERRORE ACCETTABILE per interpolazione\n";
    }
    
    // Test metodo della strategia
    std::cout << "\n🔄 PARTE 6: TEST METODO STRATEGIA\n";
    std::cout << "-----------------------------------\n";
    
    auto cheb_pos = strategy.getChebyshevPosition(target_mjd);
    std::cout << "🔵 Strategia Chebyshev: RA=" << cheb_pos.ra_deg 
              << "°, Dec=" << cheb_pos.dec_deg << "°\n";
    
    double strategy_error = std::sqrt(std::pow(cheb_pos.ra_deg - rkf78_pos.ra_deg, 2) + 
                                    std::pow(cheb_pos.dec_deg - rkf78_pos.dec_deg, 2)) * 3600.0;
    std::cout << "📏 Errore strategia: " << strategy_error << " arcsec\n";
    
    return 0;
}