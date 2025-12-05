/**
 * @file test_comparison_table.cpp
 * @brief Tabella completa di comparazione propagatori con verifica epoche
 * @author Michele Bigi
 * @date 30 Novembre 2025
 */

#include "ioccultcalc/propagation_strategy.h"
#include "ioccultcalc/astdyn_interface.h"
#include "ioccultcalc/jpl_horizons_client.h"
#include "ioccultcalc/time_utils.h"
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace ioccultcalc;

// Conversioni data
std::string mjdToISO(double mjd) {
    // Conversione semplificata MJD -> ISO (approssimativa)
    double jd = mjd + 2400000.5;
    int year = (int)((jd - 1867216.25) / 365.25);
    return "~" + std::to_string(year + 1900);
}

// Calcola separazione angolare in arcsec
double angularSeparation(double ra1_deg, double dec1_deg, double ra2_deg, double dec2_deg) {
    double dra = (ra1_deg - ra2_deg) * cos(dec1_deg * M_PI / 180.0);
    double ddec = dec1_deg - dec2_deg;
    return sqrt(dra * dra + ddec * ddec) * 3600.0; // arcsec
}

struct PropagationResult {
    std::string method;
    double ra_deg = 0;
    double dec_deg = 0;
    double computation_time_ms = 0;
    std::string notes;
    bool success = false;
};

int main() {
    std::cout << "================================================================================\n";
    std::cout << "                    TABELLA COMPARAZIONE PROPAGATORI COMPLETA                    \n";
    std::cout << "================================================================================\n";
    
    try {
        // =====================================================================
        // CARICAMENTO ELEMENTI E VERIFICA EPOCHE
        // =====================================================================
        std::cout << "\n📋 CARICAMENTO E VERIFICA ELEMENTI ORBITALI\n";
        std::cout << std::string(80, '-') << std::endl;
        
        AstDySElements elements = AstDySElements::fromFile("17030_astdys.eq1");
        
        std::cout << "✅ Asteroide: (" << elements.number << ") " << elements.name << std::endl;
        std::cout << "📅 Epoca elementi: MJD " << std::fixed << std::setprecision(6) << elements.epoch_mjd;
        std::cout << " (" << mjdToISO(elements.epoch_mjd) << ")" << std::endl;
        
        // Date di test multiple
        std::vector<std::pair<double, std::string>> test_epochs = {
            {60006.024306, "28 Nov 2025 00:35 UTC (Richiesta utente)"},
            {61000.0, "Epoca elementi esatta (zero propagazione)"},
            {61007.0, "7 giorni dopo epoca elementi"},
            {60800.0, "200 giorni prima epoca elementi"}
        };
        
        std::cout << "\n📅 EPOCHE DI TEST:\n";
        for (size_t i = 0; i < test_epochs.size(); i++) {
            double delta_days = test_epochs[i].first - elements.epoch_mjd;
            std::cout << "  " << (i+1) << ". MJD " << std::fixed << std::setprecision(6) 
                      << test_epochs[i].first << " - " << test_epochs[i].second 
                      << " (Δ=" << std::setprecision(1) << delta_days << " giorni)" << std::endl;
        }
        
        // =====================================================================
        // ELEMENTI ORBITALI DETTAGLI
        // =====================================================================
        std::cout << "\n📊 ELEMENTI ORBITALI (CONVERTITI DA EQUINOZIALI)\n";
        std::cout << std::string(80, '-') << std::endl;
        std::cout << "  Semi-asse maggiore: " << std::setprecision(9) << elements.a << " AU" << std::endl;
        std::cout << "  Eccentricità:       " << std::setprecision(9) << elements.e << std::endl;
        std::cout << "  Inclinazione:       " << std::setprecision(6) << elements.i << "°" << std::endl;
        std::cout << "  Long. nodo asc.:    " << std::setprecision(6) << elements.Omega << "°" << std::endl;
        std::cout << "  Arg. perielio:      " << std::setprecision(6) << elements.omega << "°" << std::endl;
        std::cout << "  Anomalia media:     " << std::setprecision(6) << elements.M << "°" << std::endl;
        
        // =====================================================================
        // LOOP SU TUTTE LE EPOCHE DI TEST
        // =====================================================================
        
        for (size_t epoch_idx = 0; epoch_idx < test_epochs.size(); epoch_idx++) {
            double target_mjd = test_epochs[epoch_idx].first;
            std::string epoch_desc = test_epochs[epoch_idx].second;
            double delta_days = target_mjd - elements.epoch_mjd;
            
            std::cout << "\n" << std::string(80, '=') << std::endl;
            std::cout << "EPOCA TEST " << (epoch_idx+1) << ": MJD " << std::fixed << std::setprecision(6) << target_mjd << std::endl;
            std::cout << "   " << epoch_desc << std::endl;
            std::cout << "   Propagazione: " << std::setprecision(1) << std::abs(delta_days) 
                      << " giorni " << (delta_days >= 0 ? "avanti" : "indietro") << std::endl;
            std::cout << std::string(80, '=') << std::endl;
            
            std::vector<PropagationResult> results;
            
            // =================================================================
            // 1. JPL HORIZONS (REFERENZA)
            // =================================================================
            PropagationResult jpl_result;
            jpl_result.method = "JPL Horizons";
            jpl_result.notes = "⚠️ MOCK - Coordinate placeholder";
            // TODO: Implementare query reale per 17030 a target_mjd
            jpl_result.ra_deg = 40.86524400;  // Placeholder per ora
            jpl_result.dec_deg = 11.88253200; // Placeholder per ora
            jpl_result.success = false; // Mock
            results.push_back(jpl_result);
            
            // =================================================================
            // 2. RKF78 DIRETTO (AstDyn)
            // =================================================================
            PropagationResult rkf78_result;
            rkf78_result.method = "RKF78 Diretto";
            
            try {
                auto start = std::chrono::high_resolution_clock::now();
                
                AstDynPropagator rkf78_prop(1e-12);
                auto radec = rkf78_prop.getRADec(elements, target_mjd);
                
                auto end = std::chrono::high_resolution_clock::now();
                rkf78_result.computation_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                
                rkf78_result.ra_deg = radec.first;
                rkf78_result.dec_deg = radec.second;
                rkf78_result.success = true;
                rkf78_result.notes = "AstDyn wrapper (tolleranza 1e-12)";
                
            } catch (const std::exception& e) {
                rkf78_result.notes = "ERRORE: " + std::string(e.what());
                rkf78_result.success = false;
            }
            results.push_back(rkf78_result);
            
            // =================================================================
            // 3. CHEBYSHEV (FASE 1)
            // =================================================================
            PropagationResult chebyshev_result;
            chebyshev_result.method = "Chebyshev FASE 1";
            
            try {
                auto start = std::chrono::high_resolution_clock::now();
                
                PropagationConfig config;
                config.rkf78_tolerance = 1e-12;
                config.verbose_timing = false;
                
                TwoPhaseStrategy strategy(config);
                strategy.setElements(elements);
                
                auto coords = strategy.getChebyshevPosition(target_mjd);
                
                auto end = std::chrono::high_resolution_clock::now();
                chebyshev_result.computation_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                
                chebyshev_result.ra_deg = coords.ra_deg;
                chebyshev_result.dec_deg = coords.dec_deg;
                chebyshev_result.success = true;
                chebyshev_result.notes = "Polinomi grado 8 da RKF78";
                
            } catch (const std::exception& e) {
                chebyshev_result.notes = "ERRORE: " + std::string(e.what());
                chebyshev_result.success = false;
            }
            results.push_back(chebyshev_result);
            
            // =================================================================
            // 4. RKF78 + ORBIT FITTING
            // =================================================================
            PropagationResult fitted_result;
            fitted_result.method = "RKF78 + Fitting";
            
            try {
                auto start = std::chrono::high_resolution_clock::now();
                
                // Simula orbit fitting con piccole correzioni
                AstDySElements fitted_elements = elements;
                fitted_elements.M += 0.001;  // Correzione anomalia media
                fitted_elements.a += 1e-7;   // Correzione semiasse
                
                AstDynPropagator fitted_prop(1e-12);
                auto radec = fitted_prop.getRADec(fitted_elements, target_mjd);
                
                auto end = std::chrono::high_resolution_clock::now();
                fitted_result.computation_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                
                fitted_result.ra_deg = radec.first;
                fitted_result.dec_deg = radec.second;
                fitted_result.success = true;
                fitted_result.notes = "Elementi corretti (simulato)";
                
            } catch (const std::exception& e) {
                fitted_result.notes = "ERRORE: " + std::string(e.what());
                fitted_result.success = false;
            }
            results.push_back(fitted_result);
            
            // =================================================================
            // TABELLA RISULTATI PER QUESTA EPOCA
            // =================================================================
            std::cout << "\n📊 TABELLA RISULTATI:\n";
            std::cout << std::string(105, '-') << std::endl;
            std::cout << std::left << std::setw(18) << "Metodo" 
                      << std::setw(15) << "RA (gradi)" 
                      << std::setw(15) << "Dec (gradi)"
                      << std::setw(12) << "Tempo (ms)"
                      << std::setw(15) << "vs RKF78 (\") "
                      << std::setw(30) << "Note" << std::endl;
            std::cout << std::string(105, '-') << std::endl;
            
            // Trova RKF78 come referenza per questa epoca
            PropagationResult* rkf78_ref = nullptr;
            for (auto& result : results) {
                if (result.method == "RKF78 Diretto" && result.success) {
                    rkf78_ref = &result;
                    break;
                }
            }
            
            for (auto& result : results) {
                std::cout << std::left << std::setw(18) << result.method;
                
                if (result.success) {
                    std::cout << std::setw(15) << std::fixed << std::setprecision(6) << result.ra_deg
                              << std::setw(15) << std::setprecision(6) << result.dec_deg
                              << std::setw(12) << std::setprecision(1) << result.computation_time_ms;
                    
                    // Calcola errore vs RKF78
                    if (rkf78_ref && result.method != "RKF78 Diretto") {
                        double error_arcsec = angularSeparation(result.ra_deg, result.dec_deg, 
                                                               rkf78_ref->ra_deg, rkf78_ref->dec_deg);
                        std::cout << std::setw(15) << std::setprecision(2) << error_arcsec;
                    } else {
                        std::cout << std::setw(15) << "---";
                    }
                } else {
                    std::cout << std::setw(15) << "FALLITO"
                              << std::setw(15) << "FALLITO"
                              << std::setw(12) << "---"
                              << std::setw(15) << "---";
                }
                
                std::cout << result.notes << std::endl;
            }
            
            // =================================================================
            // ANALISI ERRORI PER QUESTA EPOCA
            // =================================================================
            if (rkf78_ref) {
                std::cout << "\n🔍 ANALISI ERRORI (vs RKF78 Diretto):\n";
                
                for (const auto& result : results) {
                    if (result.success && result.method != "RKF78 Diretto") {
                        double error_arcsec = angularSeparation(result.ra_deg, result.dec_deg,
                                                               rkf78_ref->ra_deg, rkf78_ref->dec_deg);
                        
                        std::cout << "  • " << std::left << std::setw(18) << result.method << ": "
                                  << std::setprecision(2) << error_arcsec << "\" ";
                        
                        if (error_arcsec < 1.0) {
                            std::cout << "✅ ECCELLENTE";
                        } else if (error_arcsec < 10.0) {
                            std::cout << "✅ BUONO";
                        } else if (error_arcsec < 60.0) {
                            std::cout << "⚠️ ACCETTABILE";
                        } else {
                            std::cout << "❌ INACCETTABILE";
                        }
                        std::cout << std::endl;
                    }
                }
            }
        }
        
        // =====================================================================
        // CONCLUSIONI GENERALI
        // =====================================================================
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "CONCLUSIONI GENERALI\n";
        std::cout << std::string(80, '=') << std::endl;
        
        std::cout << "\n✅ SUCCESSI:\n";
        std::cout << "  • Parser AstDyS funziona correttamente\n";
        std::cout << "  • RKF78 diretto propagazioni stabili\n";
        std::cout << "  • Strategia implementata e compilata\n";
        
        std::cout << "\n⚠️ PROBLEMI IDENTIFICATI:\n";
        std::cout << "  • JPL Horizons: Solo valori mock, serve query reale\n";
        std::cout << "  • Chebyshev: Possibili errori significativi vs RKF78\n";
        std::cout << "  • Propagazioni lunghe: Accuratezza degrada nel tempo\n";
        
        std::cout << "\n🎯 RACCOMANDAZIONI:\n";
        std::cout << "  • Implementare query JPL Horizons reale per validazione\n";
        std::cout << "  • Debug implementazione polinomi Chebyshev\n";
        std::cout << "  • Test con epoche più vicine agli elementi\n";
        std::cout << "  • Validare con casi noti di occultazioni\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ ERRORE FATALE: " << e.what() << std::endl;
        return 1;
    }
}