#include "chebyshev_occultation_manager.h"
#include <iostream>

int main() {
    try {
        ioccultcalc::OccultationSearchConfig config;
        config.start_mjd = 60642.0;
        config.end_mjd = 60643.0;
        config.num_propagation_points = 50;
        config.num_chebyshev_coeffs = 8;
        
        ioccultcalc::ChebyshevOccultationManager manager(config);
        
        std::cout << "1. Load...\n";
        if (!manager.loadAsteroidFromEQ1("17030_astdys.eq1")) {
            return 1;
        }
        
        std::cout << "2. Propagate...\n";
        if (!manager.propagateAndFit()) {
            return 1;
        }
        
        std::cout << "3. Get position...\n";
        double mjd = 60642.5;
        std::cout << "Requesting position at MJD " << mjd << std::endl;
        
        Eigen::Vector3d pos = manager.getPositionAtEpoch(mjd);
        
        std::cout << "Position: " << pos.transpose() << std::endl;
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}
