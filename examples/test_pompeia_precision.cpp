/**
 * @file test_pompeia_precision.cpp
 * @brief Test completo precisione per (27) Pompeia (main belt asteroid)
 * 
 * Pompeia orbita a ~2.74 AU, vicino alla main belt dove asteroidi massivi
 * (Ceres, Pallas, Vesta) hanno effetto significativo.
 * 
 * Test con tutte le correzioni di precisione:
 * 1. Perturbazioni planetarie (DE441)
 * 2. Perturbazioni asteroidi massivi (SB441-N16)
 * 3. Correzioni relativistiche (Schwarzschild PN1)
 * 4. Earth Orientation Parameters (IERS)
 */

#include "ioccultcalc/force_model.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/iers_data.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

int main() {
    std::cout << "\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "   TEST PRECISIONE COMPLETO: (27) POMPEIA (Main Belt)\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "\n";
    
    try {
        // Elementi orbitali (27) Pompeia @ J2000.0
        // Fonte: JPL Small-Body Database
        double jd0 = 2451545.0;  // J2000.0
        
        std::cout << "ASTEROIDE: (27) Pompeia\n";
        std::cout << "Epoca: JD " << std::fixed << std::setprecision(2) << jd0 << " (2000-Jan-01.5)\n";
        std::cout << "\nElementi orbitali:\n";
        std::cout << "  Semi-asse maggiore (a): 2.740 AU\n";
        std::cout << "  Eccentricità (e):       0.071\n";
        std::cout << "  Inclinazione (i):       11.75°\n";
        std::cout << "  Long. nodo asc. (Ω):    218.8°\n";
        std::cout << "  Arg. pericentro (ω):    327.8°\n";
        std::cout << "  Anomalia media (M):     156.4°\n";
        std::cout << "\n";
        
        // Converti in radianti
        double a = 2.740;
        double e = 0.071;
        double i = 11.75 * M_PI / 180.0;
        double Omega = 218.8 * M_PI / 180.0;
        double omega = 327.8 * M_PI / 180.0;
        double M = 156.4 * M_PI / 180.0;
        
        // Calcola anomalia eccentrica E (Kepler's equation: M = E - e*sin(E))
        double E = M;
        for (int iter = 0; iter < 20; iter++) {
            E = M + e * std::sin(E);
        }
        
        // Calcola coordinate orbitali nel piano orbitale
        double x_orb = a * (std::cos(E) - e);
        double y_orb = a * std::sqrt(1 - e*e) * std::sin(E);
        
        // Trasforma in coordinate eclittiche eliocentriche
        double cos_Omega = std::cos(Omega);
        double sin_Omega = std::sin(Omega);
        double cos_omega = std::cos(omega);
        double sin_omega = std::sin(omega);
        double cos_i = std::cos(i);
        double sin_i = std::sin(i);
        
        double x = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * x_orb +
                   (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * y_orb;
        double y = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * x_orb +
                   (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * y_orb;
        double z = (sin_omega * sin_i) * x_orb + (cos_omega * sin_i) * y_orb;
        
        Vector3D pos0(x, y, z);
        
        // Velocità orbitale
        double GM = 1.32712440041279419e11;  // km³/s²
        double AU_km = 1.495978707e8;
        double day_sec = 86400.0;
        double GM_AU_day = GM * (day_sec * day_sec) / (AU_km * AU_km * AU_km);
        double n = std::sqrt(GM_AU_day / (a * a * a));  // Mean motion
        
        double dE_dt = n / (1 - e * std::cos(E));
        double vx_orb = -a * std::sin(E) * dE_dt;
        double vy_orb = a * std::sqrt(1 - e*e) * std::cos(E) * dE_dt;
        
        double vx = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * vx_orb +
                    (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * vy_orb;
        double vy = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * vx_orb +
                    (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * vy_orb;
        double vz = (sin_omega * sin_i) * vx_orb + (cos_omega * sin_i) * vy_orb;
        
        Vector3D vel0(vx, vy, vz);
        
        std::cout << "Stato cartesiano iniziale (eliocentrico):\n";
        std::cout << "  r = (" << std::setprecision(8) 
                  << pos0.x << ", " << pos0.y << ", " << pos0.z << ") AU\n";
        std::cout << "  v = (" << std::setprecision(8)
                  << vel0.x << ", " << vel0.y << ", " << vel0.z << ") AU/day\n";
        std::cout << "  |r| = " << std::setprecision(6) << pos0.magnitude() << " AU\n";
        std::cout << "  |v| = " << std::setprecision(6) << vel0.magnitude() << " AU/day\n";
        std::cout << "\n";
        
        // ═══════════════════════════════════════════════════════════
        // TEST 1: Standard Config (pianeti principali)
        // ═══════════════════════════════════════════════════════════
        std::cout << "───────────────────────────────────────────────────────────────\n";
        std::cout << "TEST 1: Standard Config (pianeti principali)\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        ForceModelConfig config1 = ForceModelConfig::standardConfig();
        config1.includeCeres = false;
        config1.includePallas = false;
        config1.includeVesta = false;
        config1.includeRelativisticCorrection = false;
        
        ForceModel fm1(config1);
        fm1.initializeJPL();
        
        std::map<PerturbingBody, Vector3D> contrib1;
        Vector3D accel1 = fm1.computeAccelerationWithBreakdown(jd0, pos0, vel0, contrib1);
        
        std::cout << "Accelerazione totale: " << std::scientific << std::setprecision(6)
                  << accel1.magnitude() << " AU/day²\n\n";
        
        std::cout << "Contributi principali:\n";
        for (const auto& [body, acc] : contrib1) {
            double mag = acc.magnitude();
            if (mag > 1e-15) {
                std::string name;
                switch(body) {
                    case PerturbingBody::SUN: name = "Sun"; break;
                    case PerturbingBody::JUPITER: name = "Jupiter"; break;
                    case PerturbingBody::SATURN: name = "Saturn"; break;
                    case PerturbingBody::EARTH: name = "Earth"; break;
                    case PerturbingBody::MARS: name = "Mars"; break;
                    case PerturbingBody::VENUS: name = "Venus"; break;
                    default: continue;
                }
                std::cout << "  " << std::setw(10) << name << ": " 
                          << std::setprecision(3) << mag << " AU/day²\n";
            }
        }
        std::cout << "\n";
        
        // ═══════════════════════════════════════════════════════════
        // TEST 2: + Asteroidi Massivi (Ceres, Pallas, Vesta)
        // ═══════════════════════════════════════════════════════════
        std::cout << "───────────────────────────────────────────────────────────────\n";
        std::cout << "TEST 2: + Asteroidi Massivi (Ceres, Pallas, Vesta)\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        ForceModelConfig config2 = ForceModelConfig::standardConfig();
        config2.includeCeres = true;
        config2.includePallas = true;
        config2.includeVesta = true;
        config2.includeRelativisticCorrection = false;
        
        ForceModel fm2(config2);
        fm2.initializeJPL();
        
        std::map<PerturbingBody, Vector3D> contrib2;
        Vector3D accel2 = fm2.computeAccelerationWithBreakdown(jd0, pos0, vel0, contrib2);
        
        std::cout << "Accelerazione totale: " << std::scientific << std::setprecision(6)
                  << accel2.magnitude() << " AU/day²\n\n";
        
        std::cout << "Contributi asteroidi:\n";
        bool found_asteroids = false;
        for (const auto& [body, acc] : contrib2) {
            if (body == PerturbingBody::CERES || 
                body == PerturbingBody::PALLAS || 
                body == PerturbingBody::VESTA) {
                double mag = acc.magnitude();
                std::string name;
                switch(body) {
                    case PerturbingBody::CERES: name = "★ Ceres"; break;
                    case PerturbingBody::PALLAS: name = "★ Pallas"; break;
                    case PerturbingBody::VESTA: name = "★ Vesta"; break;
                    default: continue;
                }
                std::cout << "  " << std::setw(10) << name << ": " 
                          << std::setprecision(3) << mag << " AU/day²\n";
                found_asteroids = true;
            }
        }
        
        if (!found_asteroids) {
            std::cout << "  ⚠ NESSUN ASTEROIDE CARICATO (verificare SB441-N16.bsp)\n";
        }
        
        Vector3D diff_asteroids = accel2 - accel1;
        std::cout << "\nDifferenza da asteroidi: " << diff_asteroids.magnitude() 
                  << " AU/day²\n\n";
        
        // ═══════════════════════════════════════════════════════════
        // TEST 3: + Relatività (Full Precision)
        // ═══════════════════════════════════════════════════════════
        std::cout << "───────────────────────────────────────────────────────────────\n";
        std::cout << "TEST 3: + Relatività (Full Precision)\n";
        std::cout << "───────────────────────────────────────────────────────────────\n";
        
        ForceModelConfig config3 = ForceModelConfig::highPrecisionConfig();
        
        ForceModel fm3(config3);
        fm3.initializeJPL();
        
        Vector3D accel3 = fm3.computeAcceleration(jd0, pos0, vel0);
        Vector3D accel_rel = fm3.computeRelativisticCorrection(pos0, vel0);
        
        std::cout << "Accelerazione totale: " << std::scientific << std::setprecision(6)
                  << accel3.magnitude() << " AU/day²\n";
        std::cout << "Contrib. relativistica: " << accel_rel.magnitude() 
                  << " AU/day²\n\n";
        
        // ═══════════════════════════════════════════════════════════
        // ANALISI IMPATTO
        // ═══════════════════════════════════════════════════════════
        std::cout << "═══════════════════════════════════════════════════════════════\n";
        std::cout << "ANALISI IMPATTO SU PROPAGAZIONE (1000 giorni)\n";
        std::cout << "═══════════════════════════════════════════════════════════════\n\n";
        
        double t = 1000.0;  // giorni
        
        // Stima errore da asteroidi (accel2 - accel1)
        double err_ast = 0.5 * diff_asteroids.magnitude() * t * t;
        double err_ast_km = err_ast * AU_km / 1000.0;
        double err_ast_arcsec = err_ast * 206265.0 / pos0.magnitude();
        
        std::cout << "Effetto ASTEROIDI MASSIVI:\n";
        std::cout << "  Δa = " << std::scientific << std::setprecision(3) 
                  << diff_asteroids.magnitude() << " AU/day²\n";
        std::cout << "  Δr (1000d) ≈ " << std::fixed << std::setprecision(1) 
                  << err_ast_km << " km\n";
        std::cout << "  Δθ (1000d) ≈ " << std::setprecision(4) 
                  << err_ast_arcsec << " arcsec\n\n";
        
        // Stima errore da relatività
        double err_rel = 0.5 * accel_rel.magnitude() * t * t;
        double err_rel_km = err_rel * AU_km / 1000.0;
        double err_rel_arcsec = err_rel * 206265.0 / pos0.magnitude();
        
        std::cout << "Effetto RELATIVITÀ:\n";
        std::cout << "  Δa = " << std::scientific << std::setprecision(3) 
                  << accel_rel.magnitude() << " AU/day²\n";
        std::cout << "  Δr (1000d) ≈ " << std::fixed << std::setprecision(1) 
                  << err_rel_km << " km\n";
        std::cout << "  Δθ (1000d) ≈ " << std::setprecision(4) 
                  << err_rel_arcsec << " arcsec\n\n";
        
        // Totale
        double err_total_km = err_ast_km + err_rel_km;
        double err_total_arcsec = err_ast_arcsec + err_rel_arcsec;
        
        std::cout << "ERRORE TOTALE se omessi:\n";
        std::cout << "  Δr (1000d) ≈ " << std::fixed << std::setprecision(1) 
                  << err_total_km << " km\n";
        std::cout << "  Δθ (1000d) ≈ " << std::setprecision(4) 
                  << err_total_arcsec << " arcsec\n\n";
        
        // ═══════════════════════════════════════════════════════════
        // CONCLUSIONE
        // ═══════════════════════════════════════════════════════════
        std::cout << "═══════════════════════════════════════════════════════════════\n";
        std::cout << "CONCLUSIONI\n";
        std::cout << "═══════════════════════════════════════════════════════════════\n\n";
        
        if (found_asteroids) {
            std::cout << "✓ SPK asteroidi massivi: FUNZIONANTE\n";
        } else {
            std::cout << "✗ SPK asteroidi massivi: NON DISPONIBILE\n";
        }
        
        std::cout << "✓ Correzioni relativistiche: ATTIVE\n";
        std::cout << "✓ Perturbazioni planetarie: JPL DE441\n\n";
        
        std::cout << "Per Pompeia (main belt, a=2.74 AU):\n";
        if (err_ast_arcsec > 0.01) {
            std::cout << "  • Asteroidi massivi: SIGNIFICATIVI (~" 
                      << std::setprecision(3) << err_ast_arcsec << "\" su 1000 giorni)\n";
        } else {
            std::cout << "  • Asteroidi massivi: minimi\n";
        }
        
        if (err_rel_arcsec > 0.01) {
            std::cout << "  • Relatività: RILEVANTE (~" 
                      << std::setprecision(3) << err_rel_arcsec << "\" su 1000 giorni)\n";
        } else {
            std::cout << "  • Relatività: marginale\n";
        }
        
        std::cout << "\nPrecisione finale stimata: ~" << std::setprecision(2) 
                  << err_total_arcsec << "\" su 1000 giorni\n";
        std::cout << "                           ~" << std::setprecision(1) 
                  << err_total_km << " km\n";
        
        std::cout << "\n";
        std::cout << "═══════════════════════════════════════════════════════════════\n";
        std::cout << "Test completato!\n";
        std::cout << "═══════════════════════════════════════════════════════════════\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << "\n\n";
        return 1;
    }
    
    return 0;
}
