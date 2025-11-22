#include <iostream>
#include <iomanip>
#include <cmath>
#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/orbit_propagator.h>
#include <ioccultcalc/orbital_elements.h>
#include <ioccultcalc/types.h>

using namespace ioccultcalc;

int main() {
    std::cout << "=================================================================\n";
    std::cout << "TEST: PROPAGAZIONE CON ELEMENTI JPL#659 (EROS)\n";
    std::cout << "=================================================================\n\n";
    
    std::cout << "Obiettivo:\n";
    std::cout << "  • Usare elementi orbitali JPL#659 (soluzione ufficiale Horizons)\n";
    std::cout << "  • Propagare da epoch 2004-Nov-02 a 2024-Nov-01 (20 anni)\n";
    std::cout << "  • Confrontare con Horizons @ 2024-Nov-01\n\n";
    
    // Elementi orbitali JPL#659 per (433) Eros
    // Source: JPL Horizons, soln ref.= JPL#659 (2021-May-24)
    std::cout << "--- ELEMENTI ORBITALI JPL#659 ---\n";
    
    OrbitalElements elem;
    elem.designation = "433";
    elem.name = "Eros";
    elem.epoch.jd = 2453311.5;  // 2004-Nov-02.00 TDB
    elem.a = 1.458269315549994;  // AU
    elem.e = 0.2228078944584026;
    elem.i = 10.8291838260782 * M_PI / 180.0;  // deg → rad
    elem.Omega = 304.4010273379536 * M_PI / 180.0;  // deg → rad  
    elem.omega = 178.665326776373 * M_PI / 180.0;  // deg → rad
    elem.M = 326.37047603642 * M_PI / 180.0;  // deg → rad
    elem.H = 10.38;
    elem.G = 0.460;
    
    std::cout << "  Epoch: JD " << std::fixed << std::setprecision(1) << elem.epoch.jd 
              << " (2004-Nov-02)\n";
    std::cout << "  a = " << std::setprecision(6) << elem.a << " AU\n";
    std::cout << "  e = " << elem.e << "\n";
    std::cout << "  i = " << std::setprecision(2) << elem.i * 180.0 / M_PI << "°\n";
    std::cout << "  Ω = " << elem.Omega * 180.0 / M_PI << "°\n";
    std::cout << "  ω = " << elem.omega * 180.0 / M_PI << "°\n";
    std::cout << "  M = " << elem.M * 180.0 / M_PI << "°\n\n";
    
    // Converti in equinoctial
    auto eq_elem = elem.toEquinoctial();
    
    // Setup propagatore
    std::cout << "--- PROPAGAZIONE ---\n";
    PropagatorOptions opts;
    opts.integrator = IntegratorType::RK4;
    opts.stepSize = 0.1;  // giorni
    
    OrbitPropagator prop(opts);
    
    // Converti elementi → stato @ epoch 2004
    auto state_2004 = prop.elementsToState(eq_elem);
    
    std::cout << "Stato iniziale @ 2004-Nov-02:\n";
    std::cout << "  r = (" << std::setprecision(6) 
              << state_2004.position.x << ", " 
              << state_2004.position.y << ", " 
              << state_2004.position.z << ") AU\n";
    double r0 = std::sqrt(state_2004.position.x * state_2004.position.x + 
                         state_2004.position.y * state_2004.position.y + 
                         state_2004.position.z * state_2004.position.z);
    std::cout << "  |r| = " << r0 << " AU\n\n";
    
    // Propaga 20 anni (2004 → 2024)
    JulianDate target_epoch;
    target_epoch.jd = 2460615.5;  // 2024-Nov-01
    double dt_days = target_epoch.jd - elem.epoch.jd;  // ~7304 giorni = ~20 anni
    
    std::cout << "Propagando " << std::setprecision(0) << dt_days << " giorni (~20 anni)...\n";
    std::cout << "(Questo potrebbe richiedere alcuni minuti con RK4 step=0.1d)\n\n";
    
    auto state_2024 = prop.propagate(state_2004, target_epoch);
    
    Vector3D r_ioc = state_2024.position;
    
    std::cout << "IOccultCalc @ 2024-Nov-01 (propagato da 2004):\n";
    std::cout << "  r = (" << std::setprecision(6)
              << r_ioc.x << ", " 
              << r_ioc.y << ", " 
              << r_ioc.z << ") AU\n";
    double r1 = std::sqrt(r_ioc.x*r_ioc.x + r_ioc.y*r_ioc.y + r_ioc.z*r_ioc.z);
    std::cout << "  |r| = " << r1 << " AU\n\n";
    
    // Confronta con Horizons @ 2024-Nov-01
    std::cout << "--- HORIZONS REFERENCE ---\n";
    
    // Valore noto da query precedente
    std::cout << "Horizons @ 2024-Nov-01 (JPL#659):\n";
    std::cout << "  r = (-0.701638, -1.362965, -0.257699) AU\n";
    std::cout << "  |r| = 1.554471 AU\n\n";
    
    Vector3D r_hor;
    r_hor.x = -0.7016379003;
    r_hor.y = -1.3629652377;
    r_hor.z = -0.2576985130;
    
    // Differenza
    std::cout << "--- RISULTATO ---\n";
    double dx = (r_ioc.x - r_hor.x) * 149597870.7;  // AU → km
    double dy = (r_ioc.y - r_hor.y) * 149597870.7;
    double dz = (r_ioc.z - r_hor.z) * 149597870.7;
    double diff = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    std::cout << "Differenza IOccultCalc vs Horizons (dopo 20 anni):\n";
    std::cout << "  ΔX = " << std::setprecision(1) << dx << " km\n";
    std::cout << "  ΔY = " << dy << " km\n";
    std::cout << "  ΔZ = " << dz << " km\n";
    std::cout << "  |Δr| = " << diff << " km\n\n";
    
    double pct_x = std::abs(dx) / diff * 100.0;
    double pct_y = std::abs(dy) / diff * 100.0;
    double pct_z = std::abs(dz) / diff * 100.0;
    
    std::cout << "Distribuzione errore:\n";
    std::cout << "  X: " << std::setprecision(1) << pct_x << "%\n";
    std::cout << "  Y: " << pct_y << "%\n";
    std::cout << "  Z: " << pct_z << "%\n\n";
    
    if (pct_z > 80.0) {
        std::cout << "⚠ ATTENZIONE: >" << pct_z << "% dell'errore è su Z!\n";
        std::cout << "   Questo suggerisce un problema di frame eclittico.\n\n";
    }
    
    std::cout << "=================================================================\n";
    std::cout << "ANALISI:\n";
    std::cout << "=================================================================\n\n";
    
    if (diff < 500) {
        std::cout << "✅ ECCELLENTE! Errore <500 km dopo 20 anni di propagazione!\n\n";
        std::cout << "IOccultCalc può replicare accuratamente le propagazioni JPL\n";
        std::cout << "quando si usano gli STESSI elementi iniziali.\n";
    } else if (diff < 2000) {
        std::cout << "✓ BUONO. Errore ~" << (int)diff << " km dopo 20 anni.\n\n";
        std::cout << "Possibili miglioramenti:\n";
        std::cout << "  • Frame di riferimento (ECLIPJ2000_DE405 vs ICRF)\n";
        std::cout << "  • Modello di forze (perturbazioni mancanti)\n";
        std::cout << "  • Step size integrazione\n";
    } else {
        std::cout << "⚠ Errore significativo: " << (int)diff << " km dopo 20 anni.\n\n";
        std::cout << "Cause probabili:\n";
        std::cout << "  1. Frame mismatch (Z-axis " << pct_z << "%)\n";
        std::cout << "  2. Perturbazioni non modellate\n";
        std::cout << "  3. Differenze numeriche accumulate\n";
    }
    
    std::cout << "\n=================================================================\n";
    
    return 0;
}
