/**
 * @brief Test minimale PlanetaryEphemeris
 */

#include <astdyn/ephemeris/PlanetaryEphemeris.hpp>
#include <astdyn/ephemeris/PlanetaryData.hpp>
#include <iostream>
#include <Eigen/Dense>

int main() {
    std::cout << "Test PlanetaryEphemeris::getPosition()\n";
    
    // Test con J2000.0 (epoca standard)
    double jd_tdb = 2451545.0;  // J2000.0 = 1 Jan 2000 12:00 TT
    
    std::cout << "JD: " << jd_tdb << " (J2000.0)\n";
    
    // Prima provo a ottenere i dati del corpo
    auto earth_data = astdyn::ephemeris::PlanetaryData::getBodyData(
        astdyn::ephemeris::CelestialBody::EARTH);
    std::cout << "Corpo: " << earth_data.name << "\n";
    std::cout << "Chiamata getPosition(EARTH, " << jd_tdb << ")...\n";
    
    try {
        Eigen::Vector3d earth = astdyn::ephemeris::PlanetaryEphemeris::getPosition(
            astdyn::ephemeris::CelestialBody::EARTH, jd_tdb);
        
        std::cout << "✓ Successo!\n";
        std::cout << "Earth position: " << earth.transpose() << " AU\n";
    } catch (const std::exception& e) {
        std::cout << "✗ Eccezione: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
