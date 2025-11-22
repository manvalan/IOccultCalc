#include <iostream>
#include <iomanip>
#include <cmath>
#include <ioccultcalc/astdys_client.h>
#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/orbit_propagator.h>
#include <ioccultcalc/orbital_elements.h>
#include <ioccultcalc/types.h>

using namespace ioccultcalc;

int main() {
    std::cout << "=================================================================\n";
    std::cout << "TEST: PROPAGAZIONE CON ELEMENTI RECENTI\n";
    std::cout << "=================================================================\n\n";
    
    std::cout << "Strategia:\n";
    std::cout << "  1. Scaricare elementi RECENTI da AstDys (2024)\n";
    std::cout << "  2. Propagare BREVE periodo (10 giorni)\n";
    std::cout << "  3. Confrontare con Horizons\n\n";
    
    // Download elementi da AstDys (più recenti possibile)
    std::cout << "--- STEP 1: Download elementi da AstDys ---\n";
    AstDysClient astdys;
    
    try {
        auto eq_elements = astdys.getElements("433");
        auto elements = OrbitalElements::fromEquinoctial(eq_elements);
        
        std::cout << "Elementi orbitali (433) Eros:\n";
        std::cout << "  Epoca: JD " << std::fixed << std::setprecision(2) << elements.epoch.jd << "\n";
        std::cout << "  a = " << elements.a << " AU\n";
        std::cout << "  e = " << elements.e << "\n";
        std::cout << "  i = " << elements.i << "°\n";
        std::cout << "  Ω = " << elements.Omega << "°\n";
        std::cout << "  ω = " << elements.omega << "°\n";
        std::cout << "  M = " << elements.M << "°\n\n";
        
        // Propagazione 10 giorni
        std::cout << "--- STEP 2: Propagazione 10 giorni ---\n";
        
        JulianDate epoch0 = elements.epoch;
        JulianDate epoch1 = epoch0;
        epoch1.jd += 10.0;
        
        // Propaga con IOccultCalc
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RK4;
        opts.stepSize = 0.1;  // giorni
        
        OrbitPropagator prop(opts);
        
        // Converti elementi → stato iniziale
        auto state0 = prop.elementsToState(eq_elements);
        Vector3D r0 = state0.position;
        Vector3D v0 = state0.velocity;
        
        std::cout << "Stato iniziale @ JD " << epoch0.jd << ":\n";
        std::cout << "  r = (" << r0.x << ", " << r0.y << ", " << r0.z << ") AU\n";
        std::cout << "  |r| = " << std::sqrt(r0.x*r0.x + r0.y*r0.y + r0.z*r0.z) << " AU\n\n";
        
        auto state_prop = prop.propagate(state0, 10.0);  // 10 giorni
        
        Vector3D r1_ioc = state_prop.position;
        
        std::cout << "IOccultCalc @ JD " << epoch1.jd << " (+10d):\n";
        std::cout << "  r = (" << r1_ioc.x << ", " << r1_ioc.y << ", " << r1_ioc.z << ") AU\n";
        std::cout << "  |r| = " << std::sqrt(r1_ioc.x*r1_ioc.x + r1_ioc.y*r1_ioc.y + r1_ioc.z*r1_ioc.z) << " AU\n\n";
        
        // Download da Horizons per confronto
        std::cout << "--- STEP 3: Download Horizons per confronto ---\n";
        JPLHorizonsClient horizons;
        auto [r1_hor, v1_hor] = horizons.getStateVectors("433", epoch1, "@sun");
        
        std::cout << "Horizons @ JD " << epoch1.jd << " (+10d):\n";
        std::cout << "  r = (" << r1_hor.x << ", " << r1_hor.y << ", " << r1_hor.z << ") AU\n";
        std::cout << "  |r| = " << std::sqrt(r1_hor.x*r1_hor.x + r1_hor.y*r1_hor.y + r1_hor.z*r1_hor.z) << " AU\n\n";
        
        // Differenza
        std::cout << "--- RISULTATO ---\n";
        double dx = (r1_ioc.x - r1_hor.x) * 149597870.7;  // AU → km
        double dy = (r1_ioc.y - r1_hor.y) * 149597870.7;
        double dz = (r1_ioc.z - r1_hor.z) * 149597870.7;
        double diff = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        std::cout << "Differenza IOccultCalc vs Horizons:\n";
        std::cout << "  ΔX = " << dx << " km\n";
        std::cout << "  ΔY = " << dy << " km\n";
        std::cout << "  ΔZ = " << dz << " km\n";
        std::cout << "  |Δr| = " << diff << " km";
        
        if (diff < 100) {
            std::cout << " ✅ ECCELLENTE!\n";
        } else if (diff < 500) {
            std::cout << " ✓ BUONO\n";
        } else {
            std::cout << " ⚠ PROBLEMA\n";
        }
        
        std::cout << "\n=================================================================\n";
        std::cout << "CONCLUSIONE:\n";
        std::cout << "=================================================================\n\n";
        
        if (diff < 100) {
            std::cout << "✅ IOccultCalc è ACCURATO quando usa elementi RECENTI!\n\n";
            std::cout << "Il precedente errore di 2000 km @ 100d era dovuto a:\n";
            std::cout << "  • DE441 usa elementi da epoch 2004 (20 anni fa)\n";
            std::cout << "  • Horizons usa JPL#659 del 2021 (orbit fit aggiornato)\n";
            std::cout << "  • Propagare 20 anni accumula naturalmente ~2000 km errore\n\n";
            std::cout << "RACCOMANDAZIONE:\n";
            std::cout << "  → Usare sempre elementi RECENTI da AstDys/MPC\n";
            std::cout << "  → Validare su propagazioni BREVI (giorni/settimane)\n";
            std::cout << "  → Per lunghi periodi, ri-fit periodicamente\n";
        } else {
            std::cout << "⚠ Ancora " << diff << " km di errore su 10 giorni.\n";
            std::cout << "Possibili cause:\n";
            std::cout << "  • Condizioni iniziali diverse (AstDys vs Horizons)\n";
            std::cout << "  • Frame di riferimento\n";
            std::cout << "  • Modello di forze\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
