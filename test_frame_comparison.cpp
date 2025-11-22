#include <iostream>
#include <cmath>
#include <SpiceUsr.h>

int main() {
    // Inizializza SPICE
    const char* home = getenv("HOME");
    if (!home) {
        std::cerr << "HOME not set" << std::endl;
        return 1;
    }
    
    std::string spkPath = std::string(home) + "/.ioccultcalc/ephemerides/linux_m13000p17000.441";
    std::string framePath1 = std::string(home) + "/.ioccultcalc/ephemerides/eclipj2000_de405.tf";
    std::string framePath2 = std::string(home) + "/.ioccultcalc/ephemerides/eclipj2000_de441.tf";
    
    // Carica kernels
    furnsh_c(spkPath.c_str());
    furnsh_c(framePath1.c_str());  // Vecchio DE405
    furnsh_c(framePath2.c_str());  // Nuovo DE441
    
    if (failed_c()) {
        char msg[1841];
        getmsg_c("LONG", 1840, msg);
        std::cerr << "SPICE error: " << msg << std::endl;
        reset_c();
        return 1;
    }
    
    // Epoca: 2024-Nov-01 = JD 2460615.5
    double jd = 2460615.5;
    double et = (jd - 2451545.0) * 86400.0;
    
    // Test Eros (433) rispetto a Sole (10)
    double state_de405[6], state_de441[6], state_j2000[6];
    double lt;
    
    // 1. Frame vecchio (ECLIPJ2000_DE405)
    spkezr_c("433", et, "ECLIPJ2000_DE405", "NONE", "10", state_de405, &lt);
    
    // 2. Frame nuovo (ECLIPJ2000_DE441)
    spkezr_c("433", et, "ECLIPJ2000_DE441", "NONE", "10", state_de441, &lt);
    
    // 3. Frame J2000 equatoriale
    spkezr_c("433", et, "J2000", "NONE", "10", state_j2000, &lt);
    
    if (failed_c()) {
        char msg[1841];
        getmsg_c("LONG", 1840, msg);
        std::cerr << "SPICE error: " << msg << std::endl;
        reset_c();
        return 1;
    }
    
    // Conversione km → AU
    const double KM_TO_AU = 1.0 / 149597870.7;
    
    std::cout << "433 Eros @ JD 2460615.5 (2024-Nov-01):\n\n";
    
    std::cout << "ECLIPJ2000_DE405 (vecchio):\n";
    std::cout << "  X = " << state_de405[0] * KM_TO_AU << " AU\n";
    std::cout << "  Y = " << state_de405[1] * KM_TO_AU << " AU\n";
    std::cout << "  Z = " << state_de405[2] * KM_TO_AU << " AU\n\n";
    
    std::cout << "ECLIPJ2000_DE441 (nuovo IAU 2006):\n";
    std::cout << "  X = " << state_de441[0] * KM_TO_AU << " AU\n";
    std::cout << "  Y = " << state_de441[1] * KM_TO_AU << " AU\n";
    std::cout << "  Z = " << state_de441[2] * KM_TO_AU << " AU\n\n";
    
    std::cout << "J2000 (equatoriale):\n";
    std::cout << "  X = " << state_j2000[0] * KM_TO_AU << " AU\n";
    std::cout << "  Y = " << state_j2000[1] * KM_TO_AU << " AU\n";
    std::cout << "  Z = " << state_j2000[2] * KM_TO_AU << " AU\n\n";
    
    // Differenza DE405 vs DE441
    double dx = (state_de441[0] - state_de405[0]);
    double dy = (state_de441[1] - state_de405[1]);
    double dz = (state_de441[2] - state_de405[2]);
    double diff_km = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    std::cout << "Differenza DE405 → DE441:\n";
    std::cout << "  ΔX = " << dx << " km\n";
    std::cout << "  ΔY = " << dy << " km\n";
    std::cout << "  ΔZ = " << dz << " km\n";
    std::cout << "  |Δr| = " << diff_km << " km\n\n";
    
    // Confronto con Horizons ECLIPTIC+ICRF
    std::cout << "Horizons ECLIPTIC+ICRF (reference):\n";
    std::cout << "  X = -0.7016379003 AU\n";
    std::cout << "  Y = -1.3629652377 AU\n";
    std::cout << "  Z = -0.2576985130 AU\n\n";
    
    // Differenza vs Horizons
    double horizons[3] = {-0.7016379003, -1.3629652377, -0.2576985130};
    double spice_au[3] = {state_de441[0] * KM_TO_AU, state_de441[1] * KM_TO_AU, state_de441[2] * KM_TO_AU};
    
    double hdx = (spice_au[0] - horizons[0]) * 149597870.7;
    double hdy = (spice_au[1] - horizons[1]) * 149597870.7;
    double hdz = (spice_au[2] - horizons[2]) * 149597870.7;
    double hdiff_km = std::sqrt(hdx*hdx + hdy*hdy + hdz*hdz);
    
    std::cout << "Differenza DE441 vs Horizons:\n";
    std::cout << "  ΔX = " << hdx << " km\n";
    std::cout << "  ΔY = " << hdy << " km\n";
    std::cout << "  ΔZ = " << hdz << " km\n";
    std::cout << "  |Δr| = " << hdiff_km << " km\n";
    
    // Cleanup
    unload_c(spkPath.c_str());
    unload_c(framePath1.c_str());
    unload_c(framePath2.c_str());
    
    return 0;
}
