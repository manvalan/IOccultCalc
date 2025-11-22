#include <iostream>
#include <fstream>
#include <cmath>
#include <ioccultcalc/jpl_horizons_client.h>
#include <ioccultcalc/spice_spk_reader.h>
#include <ioccultcalc/types.h>

using namespace ioccultcalc;

int main() {
    // Epoca: 2024-Nov-01 = JD 2460615.5
    JulianDate epoch;
    epoch.jd = 2460615.5;
    
    std::cout << "=== CONFRONTO HORIZONS vs SPK ===\n\n";
    std::cout << "Epoca: JD " << epoch.jd << " (2024-Nov-01)\n";
    std::cout << "Body: 433 Eros\n";
    std::cout << "Center: @sun (10)\n\n";
    
    // 1. Horizons (ECLIPTIC+ICRF)
    std::cout << "--- HORIZONS (ECLIPTIC+ICRF) ---\n";
    JPLHorizonsClient horizons;
    auto [r_hor, v_hor] = horizons.getStateVectors("433", epoch, "@sun");
    
    std::cout << "Position: (" << r_hor.x << ", " << r_hor.y << ", " << r_hor.z << ") AU\n";
    double dist_hor = std::sqrt(r_hor.x*r_hor.x + r_hor.y*r_hor.y + r_hor.z*r_hor.z);
    std::cout << "|r| = " << dist_hor << " AU\n\n";
    
    // 2. SPK (ECLIPJ2000_DE441)
    std::cout << "--- SPK DE441 (ECLIPJ2000_DE441) ---\n";
    SPICESPKReader spk;
    const char* home = getenv("HOME");
    std::string spkPath = std::string(home) + "/.ioccultcalc/ephemerides/linux_m13000p17000.441";
    spk.loadFile(spkPath);
    auto state_spk = spk.getState(433, epoch.jd, 10);  // 433=Eros, 10=Sun
    
    std::cout << "Position: (" << state_spk.first.x << ", " << state_spk.first.y << ", " << state_spk.first.z << ") AU\n";
    double dist_spk = std::sqrt(state_spk.first.x*state_spk.first.x + 
                                state_spk.first.y*state_spk.first.y + 
                                state_spk.first.z*state_spk.first.z);
    std::cout << "|r| = " << dist_spk << " AU\n\n";
    
    // 3. Differenza
    std::cout << "--- DIFFERENZA ---\n";
    double dx = (state_spk.first.x - r_hor.x) * 149597870.7;  // AU → km
    double dy = (state_spk.first.y - r_hor.y) * 149597870.7;
    double dz = (state_spk.first.z - r_hor.z) * 149597870.7;
    double diff = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    std::cout << "ΔX = " << dx << " km (" << (dx/diff*100) << "%)\n";
    std::cout << "ΔY = " << dy << " km (" << (dy/diff*100) << "%)\n";
    std::cout << "ΔZ = " << dz << " km (" << (dz/diff*100) << "%)\n";
    std::cout << "|Δr| = " << diff << " km\n\n";
    
    if (std::abs(dz) > 0.5 * diff) {
        std::cout << "⚠ ATTENZIONE: >" << (std::abs(dz)/diff*100) << "% dell'errore è su Z-axis!\n";
        std::cout << "Questo suggerisce un problema di frame eclittico.\n";
    }
    
    return 0;
}
