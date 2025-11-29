#include <iostream>
#include "include/ioccultcalc/asteroid_database.h"

int main() {
    using namespace ioccultcalc;
    
    std::cout << "=== MPC Data Parser Test ===" << std::endl;
    
    AsteroidDatabase db;
    
    std::string mpcFile = std::string(getenv("HOME")) + "/.ioccultcalc/data/all_numbered_asteroids.json";
    
    std::cout << "Loading from: " << mpcFile << std::endl;
    
    int count = db.importFromJson(mpcFile, false);
    
    std::cout << "\n✓ Loaded " << count << " asteroids from MPC" << std::endl;
    
    // Test alcuni asteroidi famosi
    std::vector<int> testAsteroids = {1, 2, 4, 10, 704};
    std::vector<std::string> names = {"Ceres", "Pallas", "Vesta", "Hygiea", "Interamnia"};
    
    std::cout << "\n=== Sample Asteroids ===" << std::endl;
    for (size_t i = 0; i < testAsteroids.size(); i++) {
        auto props = db.getProperties(testAsteroids[i]);
        std::cout << "(" << props.number << ") " << props.name 
                  << " - H=" << props.H 
                  << " D≈" << (int)props.diameter << " km"
                  << " a=" << props.a << " AU"
                  << std::endl;
    }
    
    // Stats
    AsteroidRange range(1, 1000000);
    auto all = db.query(range);
    
    std::cout << "\n=== Statistics ===" << std::endl;
    std::cout << "Total asteroids: " << all.size() << std::endl;
    
    // H magnitude distribution
    int bright = 0, medium = 0, faint = 0;
    for (const auto& a : all) {
        if (a.H < 10.0) bright++;
        else if (a.H < 14.0) medium++;
        else faint++;
    }
    
    std::cout << "H < 10 (very bright): " << bright << std::endl;
    std::cout << "H 10-14 (bright): " << medium << std::endl;
    std::cout << "H > 14 (faint): " << faint << std::endl;
    
    return 0;
}
