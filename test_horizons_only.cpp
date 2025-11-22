#include <ioccultcalc/jpl_horizons_client.h>
#include <iostream>
#include <iomanip>
using namespace ioccultcalc;

int main() {
    JPLHorizonsClient horizons;
    
    JulianDate t0; t0.jd = 2460615.5;  // 2024-Jan-01 UTC
    JulianDate t1; t1.jd = t0.jd + 100;  // +100 giorni
    
    std::cout << "Test: Download 2 epoche da Horizons (senza propagazione)\n";
    std::cout << "Epoca 0: JD " << t0.jd << "\n";
    std::cout << "Epoca 1: JD " << t1.jd << " (+100 giorni)\n\n";
    
    auto [pos0, vel0] = horizons.getStateVectors("433", t0);
    auto [pos1, vel1] = horizons.getStateVectors("433", t1);
    
    std::cout << "Posizione @ t0: (" << pos0.x << ", " << pos0.y << ", " << pos0.z << ") AU\n";
    std::cout << "Posizione @ t1: (" << pos1.x << ", " << pos1.y << ", " << pos1.z << ") AU\n";
    
    double dist_au = (pos1 - pos0).magnitude();
    std::cout << "\nDistanza percorsa: " << dist_au << " AU = " 
              << dist_au * 149597870.7 << " km\n";
    
    return 0;
}
