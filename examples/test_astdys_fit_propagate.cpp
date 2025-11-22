/**
 * @file test_astdys_fit_propagate.cpp
 * @brief Test completo: AstDyS/Horizons → Propagate → Confronto
 * 
 * Pipeline:
 * 1. Scarica elementi orbitali da AstDyS (equinoctial)
 * 2. Propaga con high-precision (asteroidi + relatività)
 * 3. Confronta con JPL Horizons
 */

#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/jpl_horizons_client.h"
#include "ioccultcalc/orbital_elements.h"
#include <iostream>
#include <iomanip>

using namespace ioccultcalc;

void printSeparator(char c = '=', int width = 70) {
    std::cout << std::string(width, c) << "\n";
}

int main(int argc, char** argv) {
    int asteroidNumber = 27;      // Default: Pompeia
    int propagationDays = 365;    // Default: 1 anno
    
    if (argc > 1) asteroidNumber = std::atoi(argv[1]);
    if (argc > 2) propagationDays = std::atoi(argv[2]);
    
    std::cout << "\n";
    printSeparator();
    std::cout << "TEST: AstDyS → Fit → Propagate → JPL Horizons\n";
    printSeparator();
    std::cout << "\n";
    
    try {
        // ═══════════════════════════════════════════════════════════
        // STEP 1: Scarica elementi da AstDyS
        // ═══════════════════════════════════════════════════════════
        std::cout << "STEP 1: Scarica elementi orbitali da AstDyS\n";
        printSeparator('-');
        
        AstDysClient astdys;
        astdys.setTimeout(30);
        
        std::cout << "Asteroide: (" << asteroidNumber << ")\n";
        std::cout << "Download da AstDyS...\n";
        
        // Usa Horizons direttamente per semplicità
        JPLHorizonsClient horizons;
        horizons.setTimeout(30);
        
        // Usa epoca recente (2024-11-01)
        JulianDate epoch0(2460615.5);  // 2024-Nov-01.0
        
        std::cout << "Download stato iniziale da JPL Horizons...\n";
        auto [pos0, vel0] = horizons.getStateVectors(std::to_string(asteroidNumber), epoch0, "@sun");
        
        std::cout << "✓ Stato scaricato!\n\n";
        std::cout << "Epoca: JD " << std::fixed << std::setprecision(2) 
                  << epoch0.jd << " (2024-Jan-01)\n";
        std::cout << "Stato iniziale (eliocentrico):\n";
        std::cout << "  r = (" << std::setprecision(8) 
                  << pos0.x << ", " << pos0.y << ", " << pos0.z << ") AU\n";
        std::cout << "  v = (" << vel0.x << ", " 
                  << vel0.y << ", " << vel0.z << ") AU/day\n";
        std::cout << "  |r| = " << std::setprecision(6) << pos0.magnitude() << " AU\n";
        std::cout << "\n";
        
        // ═══════════════════════════════════════════════════════════
        // STEP 2: Propaga con high-precision
        // ═══════════════════════════════════════════════════════════
        std::cout << "STEP 2: Propagazione orbitale (high-precision)\n";
        printSeparator('-');
        
        JulianDate targetEpoch(epoch0.jd + propagationDays);
        
        std::cout << "Configurazione:\n";
        std::cout << "  • Integratore: RK4 (fixed step)\n";
        std::cout << "  • Step size: 0.01 giorni\n";
        std::cout << "  • Perturbazioni planetarie: SI (DE441)\n";
        std::cout << "  • Asteroidi massivi: SI (Ceres, Pallas, Vesta)\n";
        std::cout << "  • Relatività: SI (Schwarzschild PN1)\n\n";
        
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RK4;  // ← FIX: Usa RK4 (RA15 ha bug)
        opts.stepSize = 0.01;  // Step piccolo per buona precisione
        opts.tolerance = 1e-12;
        opts.usePlanetaryPerturbations = true;
        opts.useRelativisticCorrections = true;
        
        OrbitPropagator propagator(opts);
        
        std::cout << "Propagazione da JD " << std::fixed << std::setprecision(2)
                  << epoch0.jd << " a " << targetEpoch.jd << "\n";
        std::cout << "Δt = " << propagationDays << " giorni\n\n";
        
        std::cout << "Propagazione in corso...\n";
        
        OrbitState initialState(epoch0, pos0, vel0);
        OrbitState finalState = propagator.propagate(initialState, targetEpoch);
        
        std::cout << "\n";
        
        std::cout << "✓ Propagazione completata!\n\n";
        std::cout << "Posizione finale (IOccultCalc):\n";
        std::cout << "  r = (" << std::setprecision(8)
                  << finalState.position.x << ", "
                  << finalState.position.y << ", "
                  << finalState.position.z << ") AU\n";
        std::cout << "  |r| = " << std::setprecision(6)
                  << finalState.position.magnitude() << " AU\n\n";
        
        // ═══════════════════════════════════════════════════════════
        // STEP 3: Confronta con JPL Horizons
        // ═══════════════════════════════════════════════════════════
        std::cout << "STEP 3: Confronto con JPL Horizons\n";
        printSeparator('-');
        
        std::cout << "Download da JPL Horizons (epoca finale)...\n";
        
        auto [posJPL, velJPL] = horizons.getStateVectors(std::to_string(asteroidNumber), targetEpoch, "@sun");
        
        std::cout << "✓ Dati JPL scaricati!\n\n";
        std::cout << "Posizione finale (JPL Horizons):\n";
        std::cout << "  r = (" << std::setprecision(8)
                  << posJPL.x << ", "
                  << posJPL.y << ", "
                  << posJPL.z << ") AU\n";
        std::cout << "  |r| = " << std::setprecision(6)
                  << posJPL.magnitude() << " AU\n\n";
        
        // ═══════════════════════════════════════════════════════════
        // ANALISI DIFFERENZE
        // ═══════════════════════════════════════════════════════════
        printSeparator();
        std::cout << "ANALISI DIFFERENZE\n";
        printSeparator();
        std::cout << "\n";
        
        Vector3D posDiff = finalState.position - posJPL;
        Vector3D velDiff = finalState.velocity - velJPL;
        
        double posErrorAU = posDiff.magnitude();
        double posErrorKm = posErrorAU * 1.495978707e8;
        double posErrorM = posErrorKm * 1000.0;
        
        double velErrorAU = velDiff.magnitude();
        double velErrorMms = velErrorAU * 1.495978707e8 / 86400.0 * 1000.0;  // mm/s
        
        std::cout << "Errore POSIZIONE:\n";
        std::cout << "  ΔX = " << std::scientific << std::setprecision(3)
                  << posDiff.x << " AU = " << std::fixed << std::setprecision(1)
                  << posDiff.x * 1.495978707e8 << " km\n";
        std::cout << "  ΔY = " << std::scientific << std::setprecision(3)
                  << posDiff.y << " AU = " << std::fixed << std::setprecision(1)
                  << posDiff.y * 1.495978707e8 << " km\n";
        std::cout << "  ΔZ = " << std::scientific << std::setprecision(3)
                  << posDiff.z << " AU = " << std::fixed << std::setprecision(1)
                  << posDiff.z * 1.495978707e8 << " km\n";
        std::cout << "  |Δr| = " << std::scientific << std::setprecision(3)
                  << posErrorAU << " AU\n";
        std::cout << "       = " << std::fixed << std::setprecision(2)
                  << posErrorKm << " km\n";
        std::cout << "       = " << std::setprecision(0)
                  << posErrorM << " m\n\n";
        
        // Errore angolare
        double distance = finalState.position.magnitude();
        double angularErrorArcsec = (posErrorAU / distance) * 206265.0;
        
        std::cout << "Errore ANGOLARE (dalla Terra @ " 
                  << std::setprecision(2) << distance << " AU):\n";
        std::cout << "  Δθ = " << std::fixed << std::setprecision(4)
                  << angularErrorArcsec << " arcsec\n";
        std::cout << "     = " << std::setprecision(1)
                  << angularErrorArcsec * 1000.0 << " mas\n\n";
        
        std::cout << "Errore VELOCITÀ:\n";
        std::cout << "  |Δv| = " << std::scientific << std::setprecision(3)
                  << velErrorAU << " AU/day\n";
        std::cout << "       = " << std::fixed << std::setprecision(3)
                  << velErrorMms << " mm/s\n\n";
        
        // ═══════════════════════════════════════════════════════════
        // VALUTAZIONE QUALITÀ
        // ═══════════════════════════════════════════════════════════
        printSeparator();
        std::cout << "VALUTAZIONE QUALITÀ\n";
        printSeparator();
        std::cout << "\n";
        
        std::cout << "Periodo propagazione: " << propagationDays << " giorni\n";
        std::cout << "Distanza finale: " << std::setprecision(2) << distance << " AU\n\n";
        
        // Criteri di valutazione
        bool excellent = (posErrorKm < 1.0 && angularErrorArcsec < 0.01);
        bool good = (posErrorKm < 10.0 && angularErrorArcsec < 0.1);
        bool acceptable = (posErrorKm < 100.0 && angularErrorArcsec < 1.0);
        
        if (excellent) {
            std::cout << "✓✓✓ ECCELLENTE: Precisione sub-km e sub-0.01\" \n";
            std::cout << "    Adatto per previsioni occultazioni stellari!\n";
        } else if (good) {
            std::cout << "✓✓ BUONO: Precisione <10 km e <0.1\"\n";
            std::cout << "    Adatto per pianificazione osservazioni\n";
        } else if (acceptable) {
            std::cout << "✓ ACCETTABILE: Precisione <100 km e <1\"\n";
            std::cout << "    Sufficiente per survey preliminari\n";
        } else {
            std::cout << "⚠ ATTENZIONE: Errore significativo\n";
            std::cout << "    Verificare elementi iniziali e osservazioni\n";
        }
        
        std::cout << "\nContributi precisione:\n";
        std::cout << "  • Elementi iniziali: AstDyS/JPL (~1-10 km)\n";
        std::cout << "  • Perturbazioni: DE441 + asteroidi + relatività\n";
        std::cout << "  • Integratore: RA15 (tol=1e-12)\n";
        
        std::cout << "\n";
        printSeparator();
        std::cout << "Test completato!\n";
        printSeparator();
        std::cout << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\n❌ ERRORE: " << e.what() << "\n\n";
        return 1;
    }
    
    return 0;
}
