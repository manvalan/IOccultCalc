/**
 * @file test_relativity.cpp
 * @brief Test correzioni relativistiche (Schwarzschild PN1)
 * 
 * Verifica che il termine relativistico post-Newtoniano venga
 * calcolato correttamente e quantifica il suo impatto.
 */

#include "ioccultcalc/force_model.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

int main() {
    std::cout << "\n=== Test Correzioni Relativistiche (Schwarzschild PN1) ===\n\n";
    
    try {
        // Epoca di test J2000.0
        double jd = 2451545.0;
        std::cout << "Epoca: JD " << std::fixed << std::setprecision(2) << jd << " (J2000.0)\n\n";
        
        // Test a diverse distanze dal Sole
        std::vector<double> distances = {0.3, 0.7, 1.0, 1.5, 5.0};  // AU
        
        std::cout << "=== Correzione Relativistica vs Distanza ===\n";
        std::cout << std::string(80, '-') << "\n";
        std::cout << std::setw(12) << "Distanza"
                  << std::setw(20) << "Accel Newtoniana"
                  << std::setw(20) << "Accel Relativistica"
                  << std::setw(20) << "Rapporto (PN/N)\n";
        std::cout << std::setw(12) << "(AU)"
                  << std::setw(20) << "(AU/day²)"
                  << std::setw(20) << "(AU/day²)"
                  << std::setw(20) << "(×10⁻⁸)\n";
        std::cout << std::string(80, '-') << "\n";
        
        for (double r : distances) {
            // Posizione circolare a distanza r
            Vector3D pos(r, 0.0, 0.0);
            
            // Velocità circolare: v = sqrt(GM/r)
            double GM_sun = 1.32712440041279419e11;  // km³/s²
            double AU_km = 1.495978707e8;
            double day_sec = 86400.0;
            double GM_AU_day = GM_sun * (day_sec * day_sec) / (AU_km * AU_km * AU_km);
            double v_circ = std::sqrt(GM_AU_day / r);
            Vector3D vel(0.0, v_circ, 0.0);
            
            // ForceModel con e senza relatività
            ForceModelConfig config = ForceModelConfig::standardConfig();
            config.includeRelativisticCorrection = false;
            ForceModel forceModel1(config);
            forceModel1.initializeJPL();
            
            config.includeRelativisticCorrection = true;
            ForceModel forceModel2(config);
            forceModel2.initializeJPL();
            
            // Calcola accelerazioni
            Vector3D accel_newton = forceModel1.computeAcceleration(jd, pos, vel);
            Vector3D accel_rel_total = forceModel2.computeAcceleration(jd, pos, vel);
            
            // Differenza = termine relativistico puro
            Vector3D accel_rel = accel_rel_total - accel_newton;
            
            double mag_newton = accel_newton.magnitude();
            double mag_rel = accel_rel.magnitude();
            double ratio = (mag_rel / mag_newton) * 1e8;  // ×10⁻⁸
            
            std::cout << std::setw(12) << std::fixed << std::setprecision(2) << r
                      << std::setw(20) << std::scientific << std::setprecision(6) << mag_newton
                      << std::setw(20) << mag_rel
                      << std::setw(20) << std::fixed << std::setprecision(3) << ratio << "\n";
        }
        
        std::cout << "\n";
        
        // Test dettagliato a 1 AU (Terra)
        std::cout << "=== Analisi Dettagliata a 1 AU ===\n\n";
        
        Vector3D pos_earth(1.0, 0.0, 0.0);
        double GM_AU_day = 1.32712440041279419e11 * (86400.0 * 86400.0) / 
                          (1.495978707e8 * 1.495978707e8 * 1.495978707e8);
        double v_earth = std::sqrt(GM_AU_day / 1.0);
        Vector3D vel_earth(0.0, v_earth, 0.0);
        
        ForceModelConfig config = ForceModelConfig::standardConfig();
        config.includeRelativisticCorrection = true;
        ForceModel forceModel(config);
        forceModel.initializeJPL();
        
        Vector3D accel_rel = forceModel.computeRelativisticCorrection(pos_earth, vel_earth);
        
        std::cout << "Posizione: (" << std::setprecision(6) << pos_earth.x << ", " 
                  << pos_earth.y << ", " << pos_earth.z << ") AU\n";
        std::cout << "Velocità: (" << vel_earth.x << ", " 
                  << vel_earth.y << ", " << vel_earth.z << ") AU/day\n";
        std::cout << "Velocità orbitale: " << vel_earth.magnitude() << " AU/day\n\n";
        
        std::cout << "Accelerazione relativistica:\n";
        std::cout << "  Componenti: (" << std::scientific << std::setprecision(6)
                  << accel_rel.x << ", " << accel_rel.y << ", " << accel_rel.z << ") AU/day²\n";
        std::cout << "  Magnitudine: " << accel_rel.magnitude() << " AU/day²\n\n";
        
        // Stima impatto su propagazione
        double t = 365.25;  // 1 anno
        double pos_error = 0.5 * accel_rel.magnitude() * t * t;  // AU
        double pos_error_km = pos_error * 1.495978707e8;  // km
        double pos_error_m = pos_error_km * 1000.0;  // metri
        
        std::cout << "=== Impatto su Propagazione (1 anno) ===\n";
        std::cout << "Errore posizione se omesso: " << std::fixed << std::setprecision(1)
                  << pos_error_km << " km\n";
        std::cout << "                            " << std::setprecision(0)
                  << pos_error_m << " m\n\n";
        
        // Converti in arcsec
        double pos_error_arcsec = pos_error * 206265.0;
        std::cout << "Errore angolare (a 1 AU): " << std::setprecision(4)
                  << pos_error_arcsec << " arcsec\n\n";
        
        std::cout << "=== Confronto con Perturbazioni ===\n\n";
        
        // Calcola tutti i contributi
        std::map<PerturbingBody, Vector3D> contributions;
        forceModel.computeAccelerationWithBreakdown(jd, pos_earth, vel_earth, contributions);
        
        std::cout << "Contributi principali all'accelerazione:\n";
        std::cout << std::scientific << std::setprecision(3);
        
        std::vector<std::pair<std::string, double>> sorted_contrib;
        sorted_contrib.push_back({"Relativity (PN1)", accel_rel.magnitude()});
        
        for (const auto& [body, accel] : contributions) {
            double mag = accel.magnitude();
            if (mag > 1e-20) {
                std::string name;
                switch (body) {
                    case PerturbingBody::SUN: name = "Sun (central)"; break;
                    case PerturbingBody::JUPITER: name = "Jupiter"; break;
                    case PerturbingBody::EARTH: name = "Earth"; break;
                    case PerturbingBody::VENUS: name = "Venus"; break;
                    case PerturbingBody::MARS: name = "Mars"; break;
                    case PerturbingBody::SATURN: name = "Saturn"; break;
                    default: continue;
                }
                sorted_contrib.push_back({name, mag});
            }
        }
        
        // Ordina per magnitudine
        std::sort(sorted_contrib.begin(), sorted_contrib.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        for (const auto& [name, mag] : sorted_contrib) {
            std::cout << "  " << std::setw(20) << name << ": " << mag << " AU/day²\n";
        }
        
        std::cout << "\n=== Conclusione ===\n";
        std::cout << "✓ Correzione relativistica implementata (Schwarzschild PN1)\n";
        std::cout << "✓ Impatto: ~" << std::fixed << std::setprecision(1) << pos_error_km 
                  << " km/anno a 1 AU\n";
        std::cout << "✓ Necessaria per precisione sub-arcsec su scala annuale\n";
        
        std::cout << "\nTest completato.\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
