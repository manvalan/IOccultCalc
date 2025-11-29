/**
 * Test per verificare quali frame sono disponibili nei file SPK
 */

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cstdlib>

extern "C" {
    #include "SpiceUsr.h"
}

int main() {
    std::cout << "=== TEST SPK FRAMES ===\n\n";
    
    // Configura SPICE error handling
    erract_c("SET", 0, const_cast<char*>("RETURN"));
    
    // Carica de441 per i pianeti
    const char* home = getenv("HOME");
    if (!home) {
        std::cerr << "HOME not set\n";
        return 1;
    }
    
    std::string de441Path = std::string(home) + "/.ioccultcalc/ephemerides/de430.bsp";
    furnsh_c(de441Path.c_str());
    if (failed_c()) {
        char msg[1841];
        getmsg_c("SHORT", 1840, msg);
        std::cerr << "Error loading de441: " << msg << std::endl;
        reset_c();
        return 1;
    }
    std::cout << "Loaded: " << de441Path << "\n";
    
    // Carica leapseconds kernel (necessario per SPICE)
    std::string lskPath = std::string(home) + "/.ioccultcalc/ephemerides/naif0012.tls";
    FILE* f = fopen(lskPath.c_str(), "r");
    if (f) {
        fclose(f);
        furnsh_c(lskPath.c_str());
        if (failed_c()) {
            char msg[1841];
            getmsg_c("SHORT", 1840, msg);
            std::cerr << "Warning loading leapseconds: " << msg << std::endl;
            reset_c();
        } else {
            std::cout << "Loaded: " << lskPath << "\n";
        }
    } else {
        std::cout << "Leapseconds file not found, will use raw ET\n";
    }
    
    // Carica codes_300ast
    std::string astPath = std::string(home) + "/.ioccultcalc/ephemerides/codes_300ast_20100725.bsp";
    furnsh_c(astPath.c_str());
    if (failed_c()) {
        char msg[1841];
        getmsg_c("SHORT", 1840, msg);
        std::cerr << "Error loading asteroid SPK: " << msg << std::endl;
        reset_c();
        return 1;
    }
    std::cout << "Loaded: " << astPath << "\n\n";
    
    // Data test: 29 Nov 2025
    double jd = 2460646.5;
    double et = (jd - 2451545.0) * 86400.0;  // ET in secondi da J2000
    
    std::cout << "Test date: JD " << std::fixed << std::setprecision(1) << jd << "\n";
    std::cout << "ET = " << std::setprecision(0) << et << " s from J2000\n\n";
    
    // Frame da testare
    const char* frames[] = {"J2000", "ECLIPJ2000", "DE-440", "DE-441", "MEAN_ECLIP_OF_DATE"};
    int nframes = sizeof(frames) / sizeof(frames[0]);
    
    // Target: Ceres (2000001) vs Sun (10)
    int target = 2000001;
    int observer = 10;
    
    std::cout << "Testing Ceres (2000001) vs Sun (10) with different frames:\n\n";
    
    for (int i = 0; i < nframes; i++) {
        reset_c();  // Reset errori precedenti
        
        double state[6];
        double lt;
        
        spkezr_c(std::to_string(target).c_str(), et, frames[i], "NONE", 
                 std::to_string(observer).c_str(), state, &lt);
        
        if (failed_c()) {
            char msg[1841];
            getmsg_c("SHORT", 1840, msg);
            std::cout << "Frame " << frames[i] << ": FAILED - " << msg << "\n";
        } else {
            constexpr double KM_TO_AU = 1.0 / 149597870.7;
            std::cout << "Frame " << frames[i] << ": OK\n";
            std::cout << "  X = " << std::setprecision(6) << state[0] * KM_TO_AU << " AU\n";
            std::cout << "  Y = " << state[1] * KM_TO_AU << " AU\n";
            std::cout << "  Z = " << state[2] * KM_TO_AU << " AU\n\n";
        }
    }
    
    // Test anche con sb441-n16 se esiste
    std::string sb441Path = std::string(home) + "/.ioccultcalc/ephemerides/sb441-n16.bsp";
    f = fopen(sb441Path.c_str(), "r");
    if (f) {
        fclose(f);
        
        std::cout << "\n=== TEST SB441-N16 ===\n";
        reset_c();
        furnsh_c(sb441Path.c_str());
        if (failed_c()) {
            char msg[1841];
            getmsg_c("SHORT", 1840, msg);
            std::cerr << "Error loading sb441: " << msg << std::endl;
        } else {
            std::cout << "Loaded: " << sb441Path << "\n\n";
            
            for (int i = 0; i < 2; i++) {  // Solo J2000 e ECLIPJ2000
                reset_c();
                
                double state[6];
                double lt;
                
                spkezr_c(std::to_string(target).c_str(), et, frames[i], "NONE", 
                         std::to_string(observer).c_str(), state, &lt);
                
                if (failed_c()) {
                    char msg[1841];
                    getmsg_c("SHORT", 1840, msg);
                    std::cout << "Frame " << frames[i] << ": FAILED - " << msg << "\n";
                } else {
                    constexpr double KM_TO_AU = 1.0 / 149597870.7;
                    std::cout << "Frame " << frames[i] << ": OK\n";
                    std::cout << "  X = " << std::setprecision(6) << state[0] * KM_TO_AU << " AU\n";
                    std::cout << "  Y = " << state[1] * KM_TO_AU << " AU\n";
                    std::cout << "  Z = " << state[2] * KM_TO_AU << " AU\n\n";
                }
            }
        }
    }
    
    return 0;
}
