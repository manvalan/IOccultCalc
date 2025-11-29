#include "ioccultcalc/ioc_spice_wrapper.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace ioccultcalc;

void printState(const std::string& name, const StateVector& state) {
    std::cout << name << ":\n";
    std::cout << "  Position: [" << state.x << ", " << state.y << ", " << state.z << "] km\n";
    std::cout << "  Velocity: [" << state.vx << ", " << state.vy << ", " << state.vz << "] km/s\n";
    double r = std::sqrt(state.x*state.x + state.y*state.y + state.z*state.z);
    std::cout << "  Distance: " << r << " km (" << (r/149597870.7) << " AU)\n\n";
}

int main(int argc, char** argv) {
    try {
        std::cout << "========================================\n";
        std::cout << "IOC_SPICE Wrapper Test\n";
        std::cout << "========================================\n\n";
        
        // Inizializza wrapper con parallelismo
        std::cout << "Initializing IOC_SPICE wrapper...\n";
        IOCSpiceWrapper spice(0);  // 0 = auto-detect threads
        std::cout << "Using " << spice.getNumThreads() << " threads\n\n";
        
        // Carica kernel necessari
        std::string home = std::getenv("HOME");
        std::vector<std::string> kernels = {
            home + "/.ioccultcalc/ephemerides/de441.bsp",
            home + "/cspice/lsk/naif0012.tls",
            home + "/cspice/pck/pck00010.tpc"
        };
        
        std::cout << "Loading SPICE kernels...\n";
        for (const auto& kernel : kernels) {
            std::cout << "  " << kernel << "\n";
            try {
                spice.loadKernel(kernel);
            } catch (const SpiceException& e) {
                std::cerr << "Warning: " << e.what() << "\n";
            }
        }
        std::cout << "\n";
        
        // Test 1: Conversione UTC -> ET
        std::cout << "Test 1: UTC to ET conversion\n";
        std::cout << "--------------------------------\n";
        std::string utc = "2025-12-01T00:00:00";
        double et = spice.utcToET(utc);
        std::cout << "UTC: " << utc << "\n";
        std::cout << "ET:  " << et << " seconds from J2000\n\n";
        
        // Test 2: Singolo stato (Earth)
        std::cout << "Test 2: Single state query (Earth at J2000)\n";
        std::cout << "--------------------------------\n";
        StateVector earth = spice.getState("EARTH", 0.0, "SUN", "J2000", "NONE");
        printState("Earth", earth);
        
        // Test 3: Stati multipli tempi (parallelo)
        std::cout << "Test 3: Multiple times for single body (parallel)\n";
        std::cout << "--------------------------------\n";
        std::vector<double> times;
        for (int i = 0; i < 10; ++i) {
            times.push_back(i * 86400.0);  // 10 giorni
        }
        
        auto mars_states = spice.getMultiState("MARS BARYCENTER", times, "SUN", "J2000", "NONE");
        std::cout << "Computed " << mars_states.size() << " states for Mars\n";
        std::cout << "Mars at t=0:   r = " 
                  << std::sqrt(mars_states[0].x*mars_states[0].x + 
                              mars_states[0].y*mars_states[0].y + 
                              mars_states[0].z*mars_states[0].z) / 149597870.7 
                  << " AU\n\n";
        
        // Test 4: Multipli corpi singolo tempo (parallelo)
        std::cout << "Test 4: Multiple bodies at single time (parallel)\n";
        std::cout << "--------------------------------\n";
        std::vector<std::string> planets = {
            "MERCURY BARYCENTER",
            "VENUS BARYCENTER", 
            "EARTH",
            "MARS BARYCENTER"
        };
        
        auto planet_states = spice.getMultiBodyState(planets, et, "SUN", "J2000", "NONE");
        for (size_t i = 0; i < planets.size(); ++i) {
            printState(planets[i], planet_states[i]);
        }
        
        std::cout << "========================================\n";
        std::cout << "All tests completed successfully!\n";
        std::cout << "========================================\n";
        
        return 0;
        
    } catch (const SpiceException& e) {
        std::cerr << "SPICE Error: " << e.what() << "\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
